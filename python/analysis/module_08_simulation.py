# Phys 39 Module 8: numerical TEC-rod model
#
# Used for the S18.3 "Check The Derivation" and simulation tables, the S18.4 simulation
# table, the focused 100 s P/PI simulation (S19), and the Stage 6 integrated
# theory-experiment-simulation comparison. It stands in for the course's Module 6
# one-lump app and six-sensor Angstrom simulator, with every parameter taken from our
# own measurements.
#
# Model (theta = T - T_room)
#   TEC, one lump (Module 6):
#       tau dthetaT/dt = chi(u) u - thetaT - g (thetaT - theta(X=0))
#       chi(u) = chi_h for u > 0 (heating), chi_c for u < 0 (cooling); g = G/H is the
#       TEC-to-rod contact conductance over the TEC-to-room conductance. g = 0 and no
#       rod is the plain Module 6 one-lump model.
#   Rod, one-dimensional fin (Stages 5-6), X from the TEC face to the tip (28.66 cm):
#       dtheta/dt = kappa d2theta/dX2 - nu theta
#       base: -dtheta/dX = (thetaT - theta)/l_c (contact, l_c = kA/G); tip: dtheta/dX = 0
#       T1 sits at X = X1 and the other sensors at X1 + x_i (x_i from the sensor map).
#   Finite differences in X (2 mm), integrated exactly over each record interval with
#   the matrix exponential: the input is held between records (as the Arduino holds
#   its PWM), and each simulated record is the average over its interval (as the
#   Arduino averages 500 ADC readings over about 1 s).
#   Controller: the sketch's updateControl() (PWM_SINE, FF, P, PI with conditional
#   integration, cooling commands multiplied by r, clamp at +-u_max, integer PWM),
#   fed back from the TEC node like the apparatus (the course simulator feeds back
#   from the rod base T0; ours does not, so no gain conversion is needed).
#
# Parameters, all measured:
#   kappa, nu      weighted mean of the Angstrom results (analysis/angstrom_<name>.csv)
#   l_c, X1        fitted to the measured TEC -> T1 transfer (amplitude ratio and lag)
#                  of the PI runs at 100, 200, 300 s, using the finite-rod transfer
#                  function cosh[s(Lr - X1)] / (cosh(s Lr) + l_c s sinh(s Lr))
#   chi_h, chi_c,  fitted by replaying the recorded PWM of the raw-PWM run
#   tau, g         (S18.3-S18.4, then the cool-down) and matching the TEC record
#   one-lump only  chi_h 0.41 (S18.4 static), chi_c 0.17 (S18.3 mean shift),
#                  tau = mean S18.3 tau_eff (control_s18_3_7.csv)
#
# Subcommands (default: all):
#   calibrate   l_c, X1 and the coupled-model chi_h, chi_c, tau, g -> sim_parameters.csv
#   s18-3       course six-run one-lump table (tau 80 s, chi_h 0.46, chi_c 0.23/0.46,
#               u_b 50, room 22 C); 300 s numerical vs analytical overlay; 100 s
#               model-experiment comparison with the measured parameters
#   s18-4       zero-centered and positive-only 100 s drives: settled fundamental and
#               residual amplitudes (one-lump and coupled) next to the measured ones
#   control     focused 100 s P vs PI at the same Kp, the targeted adjustment (Kp 30,
#               Ki 0.3), and the selected gains at 200 and 300 s, against the measured runs
#   rod LABEL   integrated comparison for one production interval: simulated PI run at
#               that period vs the measured TEC and T1-T5, aligned on the setpoint phase,
#               with the analytical rod prediction from the measured T1; --independent
#               takes kappa and nu from the other periods only (an independent test)
#
# Outputs: data/module_08/analysis/sim_*.csv, docs/figures/module_08/sim_*.png
#
# Usage (from the repository root):
#   python python/analysis/module_08_simulation.py            # everything
#   python python/analysis/module_08_simulation.py rod p300   # one integrated comparison

import argparse
import csv
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
import angstrom_lib as al  # noqa: E402
from module_08_common import (ROD_ANALYSIS_LENGTH_M, ROD_COLUMNS, ROD_OVERALL_LENGTH_M,  # noqa: E402
                              SENSOR_LABELS, SENSOR_POSITIONS_CM, TEMPERATURE_COLUMNS)

DATA = al.MODULE8_DATA_DIR
ANALYSIS = al.ANALYSIS_DIR
FIGS = al.REPO_ROOT / "docs" / "figures" / "module_08"

RAW_RUN = "data/module_08/angstrom_20261116-211726.csv"   # S18.3-S18.4 raw PWM
RAW_BASELINE = (1390.0, 1691.0)                            # step 28, TEC off
RAW_WINDOW = (1691.0, 9310.0)                              # steps 32-36 (36: cool-down)
RAW_STEPS = {300: 32, 200: 33, 100: 34, "pos100": 35}

# Control intervals whose TEC -> T1 transfer calibrates the contact (selected PI gains)
CONTACT_LABELS = ["s18_7_pi100_t2", "s18_7_pi200_sel", "s18_7_pi300_sel"]


# ------------------------------------------------------------------
# MODEL
# ------------------------------------------------------------------
@dataclass
class Lump:
    chi_h: float
    chi_c: float
    tau: float
    g: float = 0.0


@dataclass
class Rod:
    kappa: float
    nu: float
    contact_m: float
    x1_m: float
    length_m: float = ROD_OVERALL_LENGTH_M
    dx_m: float = 0.002
    positions_cm: dict = field(default_factory=lambda: dict(SENSOR_POSITIONS_CM))


class Model:
    """Linear state model y' = A y + b v with v = chi(u) u, discretized for a record
    interval dt; outputs TEC and T1-T5 (theta, averaged over the interval)."""

    def __init__(self, lump: Lump, rod: Rod = None, dt=1.0):
        self.lump, self.rod, self.dt = lump, rod, dt
        tau = lump.tau
        if rod is None:
            A = np.array([[-1.0 / tau]])
            C = np.zeros((len(TEMPERATURE_COLUMNS), 1))
            C[0, 0] = 1.0
            self.x_nodes = None
        else:
            N = max(10, int(round(rod.length_m / rod.dx_m)))
            dx = rod.length_m / N
            n = N + 2                        # TEC + nodes 0..N
            A = np.zeros((n, n))
            A[0, 0] = -(1 + lump.g) / tau
            A[0, 1] = lump.g / tau
            k = rod.kappa / dx**2
            beta = 2 * dx / rod.contact_m
            A[1, 0] = k * beta
            A[1, 1] = -2 * k - k * beta - rod.nu
            A[1, 2] = 2 * k
            for j in range(1, N):
                A[j + 1, j] = k
                A[j + 1, j + 1] = -2 * k - rod.nu
                A[j + 1, j + 2] = k
            A[N + 1, N] = 2 * k
            A[N + 1, N + 1] = -2 * k - rod.nu
            self.x_nodes = np.linspace(0, rod.length_m, N + 1)
            C = np.zeros((len(TEMPERATURE_COLUMNS), n))
            C[0, 0] = 1.0
            for i, c in enumerate(ROD_COLUMNS, start=1):
                X = rod.x1_m + rod.positions_cm[c] / 100
                j = min(int(X / dx), N - 1)
                f = X / dx - j
                C[i, j + 1] += 1 - f
                C[i, j + 2] += f
        b = np.zeros(A.shape[0])
        b[0] = 1.0 / tau
        n = A.shape[0]
        K = np.zeros((3 * n, 3 * n))
        K[:n, :n] = A
        K[:n, n:2 * n] = np.eye(n)
        K[n:2 * n, 2 * n:] = np.eye(n)
        E = expm(K * dt)
        S1, S2 = E[:n, n:2 * n], E[:n, 2 * n:]
        self.A, self.b, self.C = A, b, C
        self.Phi, self.Gam = E[:n, :n], S1 @ b
        self.PsiY, self.PsiV = S1 / dt, S2 @ b / dt
        self.CPhi, self.CGam = C @ self.Phi, C @ self.Gam
        self.CPsiY, self.CPsiV = C @ self.PsiY, C @ self.PsiV

    def drive(self, u):
        return (self.lump.chi_h if u > 0 else self.lump.chi_c) * u

    def initial_state(self, theta_tec, theta_rod=None):
        """State from the TEC theta and the five rod-sensor thetas (interpolated along
        the rod, constant beyond the end sensors)."""
        y = np.zeros(self.A.shape[0])
        y[0] = theta_tec
        if self.x_nodes is not None:
            if theta_rod is None:
                y[1:] = theta_tec
            else:
                xs = np.array([self.rod.x1_m + self.rod.positions_cm[c] / 100 for c in ROD_COLUMNS])
                y[1:] = np.interp(self.x_nodes, xs, theta_rod)
        return y

    def steady_gain(self):
        """Static outputs per unit drive v (theta = -C A^-1 b v)."""
        return -self.C @ np.linalg.solve(self.A, self.b)


class Controller:
    """Port of updateControl() in module_08_angstrom_controller.ino."""

    def __init__(self, mode, period, Ta=40.0, Tb=16.0, kp=0.0, ki=0.0, u_a=0.0, u_b=0.0,
                 r=2.82, chi_H=0.485, T_room_set=22.0, u_max=255.0, ff_add=False,
                 integer=True, dt=1.0):
        self.__dict__.update(locals())
        del self.__dict__["self"]
        self.phase, self.q = 0.0, 0.0

    def scale(self, u):
        return u * self.r if u < 0 else u

    def __call__(self, t_tec):
        self.phase = (self.phase + 2 * math.pi * self.dt / self.period) % (2 * math.pi)
        s = math.sin(self.phase)
        tset = float("nan")
        if self.mode == "PWM_SINE":
            u = self.u_a + self.u_b * s
        else:
            tset = self.Ta + self.Tb * s
            uff = (tset - self.T_room_set) / self.chi_H
            if self.mode == "FF":
                u = self.scale(uff)
            else:
                e = tset - t_tec
                ffp = uff if self.ff_add else 0.0
                if self.mode == "PI":
                    q_next = self.q + e * self.dt
                    u_next = self.scale(ffp + self.kp * e + self.ki * q_next)
                    hold = (u_next > self.u_max and e > 0) or (u_next < -self.u_max and e < 0)
                    if not hold:
                        self.q = q_next
                else:
                    self.q = 0.0
                u = self.scale(ffp + self.kp * e + self.ki * self.q)
        u = max(-self.u_max, min(self.u_max, u))
        if self.integer:
            u = float(int(u + 0.5) if u >= 0 else int(u - 0.5))
        return u, tset


def simulate(model: Model, room, n_records, controller=None, pwm=None, y0=None):
    """Records 1..n: time, averaged temperatures (C), the PWM computed at that record
    (applied over the next interval, as in the sketch), and the setpoint.
    Give either a controller or a PWM sequence pwm[k] applied over interval k."""
    y = np.zeros(model.A.shape[0]) if y0 is None else y0.copy()
    out = {c: np.empty(n_records) for c in TEMPERATURE_COLUMNS}
    out["signed_pwm"] = np.empty(n_records)
    out["setpoint_C"] = np.full(n_records, np.nan)
    u = 0.0
    if controller is not None:
        u, _ = controller(room + y[0])
    for k in range(n_records):
        if pwm is not None:
            u = pwm[k]
        v = model.drive(u)
        temps = model.CPsiY @ y + model.CPsiV * v + room
        y = model.Phi @ y + model.Gam * v
        for i, c in enumerate(TEMPERATURE_COLUMNS):
            out[c][k] = temps[i]
        if controller is not None:
            u, tset = controller(temps[0])
            out["setpoint_C"][k] = tset
        out["signed_pwm"][k] = u
    out["time_s"] = model.dt * np.arange(1, n_records + 1)
    out["state"] = y
    return out


def fit_cycles(rec, period, n_cycles, columns=TEMPERATURE_COLUMNS, dt=None):
    """Common time fit over the last n_cycles of a simulated record; temperatures at
    t - dt/2 and the PWM at t + dt/2 (the alignment used for the measured data)."""
    t = rec["time_s"]
    dt = dt or float(t[1] - t[0])
    n = int(round(n_cycles * period / dt))
    sl = slice(len(t) - n, len(t))
    w = 2 * math.pi / period
    fits = {c: al.fit_sine(t[sl] - dt / 2, rec[c][sl], w) for c in columns}
    fits["signed_pwm"] = al.fit_sine(t[sl] + dt / 2, rec["signed_pwm"][sl], w)
    if np.isfinite(rec["setpoint_C"][sl]).all():
        fits["setpoint_C"] = al.fit_sine(t[sl], rec["setpoint_C"][sl], w)
    return fits, sl


# ------------------------------------------------------------------
# PARAMETERS
# ------------------------------------------------------------------
def kappa_nu_from(name="production", exclude=()):
    path = ANALYSIS / f"angstrom_{name}.csv"
    rows = [r for r in csv.DictReader(path.open()) if r.get("kappa_used")
            and r["label"] not in exclude]
    k, sk, *_ = al.weighted_mean([float(r["kappa_used"]) for r in rows],
                                 [float(r["kappa_used_sigma"]) for r in rows])
    n, sn, *_ = al.weighted_mean([float(r["nu_used"]) for r in rows],
                                 [float(r["nu_used_sigma"]) for r in rows])
    return k, sk, n, sn, [r["label"] for r in rows]


def contact_transfer(omega, kappa, nu, contact_m, x1_m, Lr=ROD_OVERALL_LENGTH_M):
    """Complex T1/TEC ratio of the fin with a contact resistance at its base."""
    s = np.sqrt((nu + 1j * omega) / kappa)
    return np.cosh(s * (Lr - x1_m)) / (np.cosh(s * Lr) + contact_m * s * np.sinh(s * Lr))


def measured_transfer(labels, path=DATA / "intervals_control.csv"):
    ivs = {iv.label: iv for iv in al.read_intervals(path)}
    rows = []
    for lab in labels:
        seg = al.select_segment(ivs[lab])
        t = seg.data["time_s"]
        sel, n = al.complete_cycles(t, seg.period)
        tec = al.fit_sine(t[sel], seg.data["TEC_C"][sel], seg.omega)
        t1 = al.fit_sine(t[sel], seg.data["T1_C"][sel], seg.omega)
        rows.append({"label": lab, "period_s": seg.period, "omega": seg.omega,
                     "ratio": t1.amp / tec.amp, "lag_rad": al.wrap(tec.phase - t1.phase),
                     "tec_mean": tec.mean, "t1_mean": t1.mean, "n_cycles": n})
    return rows


def calibrate_contact(kappa, nu, rows):
    def res(p):
        out = []
        for r in rows:
            h = contact_transfer(r["omega"], kappa, nu, p[0], p[1])
            out += [math.log(abs(h) / r["ratio"]), -np.angle(h) - r["lag_rad"]]
        return out
    f = least_squares(res, [0.04, ROD_OVERALL_LENGTH_M - ROD_ANALYSIS_LENGTH_M],
                      bounds=([0.0, 0.0], [0.3, 0.06]))
    J = f.jac
    s2 = float(f.fun @ f.fun) / max(1, len(f.fun) - 2)
    cov = s2 * np.linalg.inv(J.T @ J)
    return f.x[0], f.x[1], np.sqrt(np.diag(cov)), np.array(f.fun)


def load_raw_window():
    parts = al.load_parts(al.resolve_path(RAW_RUN))
    p = parts[0]
    base = (p["time_s"] >= RAW_BASELINE[0]) & (p["time_s"] <= RAW_BASELINE[1])
    t_inf = {c: float(np.nanmean(p[c][base])) for c in TEMPERATURE_COLUMNS}
    m = (p["time_s"] >= RAW_WINDOW[0]) & (p["time_s"] <= RAW_WINDOW[1])
    d = {k: v[m] for k, v in p.items()}
    return d, t_inf


def replay(model, d, t_inf, room):
    """Drive the model with the recorded PWM (each record's PWM held until the next)
    from the measured initial temperatures; returns simulated records aligned with
    d[1:]."""
    th0 = [d[c][0] - t_inf[c] for c in TEMPERATURE_COLUMNS]
    y0 = model.initial_state(th0[0], np.array(th0[1:]))
    pwm = d["signed_pwm"].astype(float)
    rec = simulate(model, room, len(pwm) - 1, pwm=pwm[:-1], y0=y0)
    # the model is in theta; express it per channel against that channel's baseline
    for c in TEMPERATURE_COLUMNS:
        rec[c] = rec[c] - room + t_inf[c]
    rec["time_s"] = d["time_s"][1:]
    rec["signed_pwm"] = pwm[1:]
    return rec


def fit_coupled_lump(rod, d, t_inf, dt):
    """chi_h, chi_c, tau, g of the TEC lump with the rod attached, from the raw-PWM
    record (TEC channel)."""
    room = t_inf["TEC_C"]
    meas = d["TEC_C"][1:]

    def res(p):
        lump = Lump(*np.exp(p))
        rec = replay(Model(lump, rod, dt), d, t_inf, room)
        return rec["TEC_C"] - meas
    p0 = np.log([0.45, 0.17, 60.0, 0.5])
    f = least_squares(res, p0, x_scale=1.0, diff_step=1e-3)
    J = f.jac
    s2 = float(f.fun @ f.fun) / (len(f.fun) - 4)
    # residuals are strongly correlated (smooth model error), so the formal
    # uncertainties are lower bounds; scale by the lag-1 correlation as for the fits
    r1 = float(np.corrcoef(f.fun[:-1], f.fun[1:])[0, 1])
    infl = (1 + r1) / max(1e-3, 1 - r1)
    cov = s2 * infl * np.linalg.inv(J.T @ J)
    vals = np.exp(f.x)
    sig = vals * np.sqrt(np.diag(cov))
    return Lump(*vals), sig, math.sqrt(float(np.mean(f.fun**2)))


def control_rows():
    path = ANALYSIS / "control_s18_3_7.csv"
    return {r["label"]: r for r in csv.DictReader(path.open())}


def one_lump_measured():
    """Measured one-lump parameters with the rod mounted: chi_h from the S18.4 static
    rise, chi_c from the S18.3 mean shift (chi_h - chi_c) u_b / pi, tau = mean tau_eff."""
    rows = control_rows()
    room = float(rows["s18_4_pos100"]["room_C"])
    chi_h = (float(rows["s18_4_pos100"]["TEC_mean_C"]) - room) / float(rows["s18_4_pos100"]["u_a"])
    shift = np.mean([float(rows[f"s18_3_pwm{P}"]["TEC_mean_minus_room_C"]) for P in (100, 200, 300)])
    u_b = float(rows["s18_3_pwm100"]["u_b"])
    chi_c = chi_h - math.pi * shift / u_b
    tau = float(np.mean([float(rows[f"s18_3_pwm{P}"]["tau_eff_s"]) for P in (100, 200, 300)]))
    return Lump(chi_h, chi_c, tau)


PARAM_FILE = ANALYSIS / "sim_parameters.csv"


def save_parameters(par):
    with PARAM_FILE.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "value", "sigma", "unit", "source"])
        for row in par:
            w.writerow(row)


def load_parameters():
    if not PARAM_FILE.exists():
        raise SystemExit(f"{PARAM_FILE} missing: run the 'calibrate' subcommand first")
    p = {r["name"]: float(r["value"]) for r in csv.DictReader(PARAM_FILE.open())}
    rod = Rod(p["kappa"], p["nu"], p["contact_m"], p["x1_m"])
    coupled = Lump(p["chi_h"], p["chi_c"], p["tau_s"], p["g"])
    one = Lump(p["one_lump_chi_h"], p["one_lump_chi_c"], p["one_lump_tau_s"])
    return rod, coupled, one, p


# ------------------------------------------------------------------
# PLOT HELPERS
# ------------------------------------------------------------------
def save(fig, name):
    FIGS.mkdir(parents=True, exist_ok=True)
    path = FIGS / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path.relative_to(al.REPO_ROOT)


def write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(rows[0].keys())
    for r in rows[1:]:
        keys += [k for k in r if k not in keys]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    return path.relative_to(al.REPO_ROOT)


def lag_deg(ref, resp):
    return math.degrees(al.wrap(ref.phase - resp.phase))


# ------------------------------------------------------------------
# CALIBRATE
# ------------------------------------------------------------------
def cmd_calibrate(args):
    kappa, sk, nu, sn, labels = kappa_nu_from(args.name)
    print(f"kappa = {kappa:.4e} ± {sk:.1e} m^2/s, nu = {nu:.3e} ± {sn:.1e} 1/s "
          f"(weighted mean of {', '.join(labels)} in angstrom_{args.name}.csv)")
    rows = measured_transfer(CONTACT_LABELS)
    lc, x1, sig, res = calibrate_contact(kappa, nu, rows)
    print(f"contact: l_c = {100 * lc:.2f} ± {100 * sig[0]:.2f} cm, X1 = {100 * x1:.2f} ± "
          f"{100 * sig[1]:.2f} cm (geometric: {100 * (ROD_OVERALL_LENGTH_M - ROD_ANALYSIS_LENGTH_M):.2f} cm)")
    for r, (ra, rp) in zip(rows, res.reshape(-1, 2)):
        h = contact_transfer(r["omega"], kappa, nu, lc, x1)
        print(f"  {r['label']:16s} {r['period_s']:5.0f} s  T1/TEC {r['ratio']:.3f} (model {abs(h):.3f}), "
              f"lag {r['lag_rad']:.3f} rad (model {-np.angle(h):.3f})")
    m = math.sqrt(nu / kappa)
    Lr = ROD_OVERALL_LENGTH_M
    mean_ratio_model = math.cosh(m * (Lr - x1)) / (math.cosh(m * Lr) + lc * m * math.sinh(m * Lr))
    rod = Rod(kappa, nu, lc, x1)

    d, t_inf = load_raw_window()
    dt = float(np.median(np.diff(d["time_s"])))
    print(f"fitting the coupled TEC lump to the raw-PWM record ({len(d['time_s'])} records, "
          f"dt {dt:.3f} s); this takes about a minute")
    lump, lsig, rms = fit_coupled_lump(rod, d, t_inf, dt)
    print(f"coupled lump: chi_h = {lump.chi_h:.3f} ± {lsig[0]:.3f}, chi_c = {lump.chi_c:.3f} ± "
          f"{lsig[1]:.3f} C/PWM, tau = {lump.tau:.1f} ± {lsig[2]:.1f} s, g = {lump.g:.3f} ± "
          f"{lsig[3]:.3f}; TEC rms residual {rms:.2f} C")
    one = one_lump_measured()
    print(f"one-lump (measured with the rod): chi_h {one.chi_h:.3f}, chi_c {one.chi_c:.3f} C/PWM, "
          f"tau {one.tau:.1f} s")
    par = [["kappa", kappa, sk, "m^2/s", f"angstrom_{args.name}.csv weighted mean"],
           ["nu", nu, sn, "1/s", f"angstrom_{args.name}.csv weighted mean"],
           ["contact_m", lc, sig[0], "m", "TEC->T1 transfer of " + " ".join(CONTACT_LABELS)],
           ["x1_m", x1, sig[1], "m", "same fit; geometric TEC face to hole 1 = 0.0176 m"],
           ["chi_h", lump.chi_h, lsig[0], "C/PWM", "replay of the raw-PWM run, TEC channel"],
           ["chi_c", lump.chi_c, lsig[1], "C/PWM", "same"],
           ["tau_s", lump.tau, lsig[2], "s", "same"],
           ["g", lump.g, lsig[3], "", "same (G/H)"],
           ["replay_rms_C", rms, "", "C", "TEC rms residual of that fit"],
           ["one_lump_chi_h", one.chi_h, "", "C/PWM", "S18.4 static rise"],
           ["one_lump_chi_c", one.chi_c, "", "C/PWM", "S18.3 mean shift"],
           ["one_lump_tau_s", one.tau, "", "s", "mean S18.3 tau_eff"],
           ["mean_ratio_T1_TEC_model", mean_ratio_model, "", "", "static (T1-room)/(TEC-room)"]]
    save_parameters(par)
    print(f"wrote {PARAM_FILE.relative_to(al.REPO_ROOT)}")

    # validation figure: the replay with all six channels
    model = Model(lump, rod, dt)
    rec = replay(model, d, t_inf, t_inf["TEC_C"])
    t = rec["time_s"]
    fig, axes = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True,
                             gridspec_kw={"height_ratios": [1, 2.2, 2.2]})
    axes[0].plot(t, rec["signed_pwm"], color=al.REFERENCE_COLOR, lw=0.7)
    axes[0].set_ylabel("signed PWM")
    for c in ["TEC_C", "T1_C"]:
        axes[1].plot(t, d[c][1:], color=al.SENSOR_COLORS[c], lw=0.8, label=f"{SENSOR_LABELS[c]} measured")
        axes[1].plot(t, rec[c], color=al.SENSOR_COLORS[c], lw=1.0, ls="--", label=f"{SENSOR_LABELS[c]} model")
    for c in ROD_COLUMNS[1:]:
        axes[2].plot(t, d[c][1:], color=al.SENSOR_COLORS[c], lw=0.8, label=f"{SENSOR_LABELS[c]} measured")
        axes[2].plot(t, rec[c], color=al.SENSOR_COLORS[c], lw=1.0, ls="--", label=f"{SENSOR_LABELS[c]} model")
    for ax in axes[1:]:
        ax.set_ylabel("Temperature (C)")
        ax.legend(ncol=4, fontsize=7, frameon=False)
    for ax in axes:
        al.style_axes(ax)
    for P, s in RAW_STEPS.items():
        t0 = d["time_s"][d["step"] == s][0]
        for ax in axes:
            ax.axvline(t0, color=al.MUTED, lw=0.6, ls=":")
        axes[0].text(t0 + 20, 45, f"{P} s" if P != "pos100" else "S18.4", fontsize=8)
    axes[2].set_xlabel("Arduino time (s), run angstrom_20261116-211726")
    axes[0].set_title("Coupled TEC-rod model driven by the recorded PWM (S18.3, S18.4, cool-down). "
                      "Fitted to TEC only; T1-T5 are predictions", fontsize=9)
    print(f"figure: {save(fig, 'sim_raw_pwm_replay.png')}")
    cols = {}
    for c in TEMPERATURE_COLUMNS:
        cols[c] = math.sqrt(float(np.mean((rec[c] - d[c][1:]) ** 2)))
    print("replay rms residual by channel (C): " +
          ", ".join(f"{SENSOR_LABELS[c]} {v:.2f}" for c, v in cols.items()))


# ------------------------------------------------------------------
# S18.3
# ------------------------------------------------------------------
def periodic_run(model, room, ctrl, period, n_cycles_total, n_fit=3):
    n = int(round(n_cycles_total * period / model.dt))
    rec = simulate(model, room, n, controller=ctrl)
    fits, sl = fit_cycles(rec, period, n_fit, dt=model.dt)
    return rec, fits, sl


def cmd_s18_3(args):
    rod, coupled, one, p = load_parameters()
    # (a) course table: tau 80 s, chi_h 0.46, chi_c 0.23 / 0.46, u_b 50, room 22 C
    rows = []
    store = {}
    for P in (100, 200, 300):
        for chi_c in (0.23, 0.46):
            model = Model(Lump(0.46, chi_c, 80.0), None, dt=0.1)
            ctrl = Controller("PWM_SINE", P, u_a=0, u_b=50, integer=False, dt=0.1)
            rec, fits, sl = periodic_run(model, 22.0, ctrl, P, max(12, int(1200 / P)))
            T, U = fits["TEC_C"], fits["signed_pwm"]
            rows.append({"period_s": P, "chi_c": chi_c, "mean_minus_room_C": T.mean - 22.0,
                         "amp_C": T.amp, "lag_deg": lag_deg(U, T), "resid_amp_C": T.half_pp_resid})
            store[(P, chi_c)] = (rec, fits, sl)
            w = 2 * math.pi / P
            pred_amp = (0.46 + chi_c) / 2 * 50 / math.sqrt(1 + (w * 80) ** 2)
            pred_mean = (0.46 - chi_c) * 50 / math.pi
            rows[-1].update({"analytic_amp_C": pred_amp, "analytic_mean_C": pred_mean,
                             "analytic_lag_deg": math.degrees(math.atan(w * 80))})
    print("S18.3 course table (one-lump, tau 80 s, chi_h 0.46 C/PWM, u_b 50, room 22 C)")
    print("| Period (s) | chi_c (C/PWM) | Mean - room (C) | Amplitude (C) | Lag (deg) | Residual (C) |")
    print("| ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in rows:
        print(f"| {r['period_s']} | {r['chi_c']:.2f} | {r['mean_minus_room_C']:.2f} | {r['amp_C']:.2f} | "
              f"{r['lag_deg']:.1f} | {r['resid_amp_C']:.3f} |")
    print(f"table: {write_rows(ANALYSIS / 'sim_s18_3_course.csv', rows)}")

    # aligned plots for both chi_c at 100 s: PWM, temperature, residual
    fig, axes = plt.subplots(3, 2, figsize=(11, 7), sharex="col")
    for j, chi_c in enumerate((0.23, 0.46)):
        rec, fits, sl = store[(100, chi_c)]
        t = rec["time_s"][sl] - rec["time_s"][sl][0]
        T = fits["TEC_C"]
        axes[0, j].plot(t, rec["signed_pwm"][sl], color=al.REFERENCE_COLOR)
        axes[1, j].plot(t, rec["TEC_C"][sl], color=al.INK, label="simulated")
        axes[1, j].plot(t, T.mean + T.amp * np.cos(2 * math.pi / 100 * (rec["time_s"][sl] - 0.05) + T.phase),
                        color=al.SENSOR_COLORS["T1_C"], ls="--", label="mean + fundamental")
        axes[2, j].plot(T.t - T.t[0], T.resid, color=al.SENSOR_COLORS["T2_C"])
        axes[0, j].set_title(f"100 s, chi_h 0.46, chi_c {chi_c:.2f} C/PWM: residual "
                             f"{T.half_pp_resid:.3f} C", fontsize=9)
        axes[2, j].set_xlabel("time in the settled window (s)")
        axes[1, j].legend(fontsize=7, frameon=False)
    axes[0, 0].set_ylabel("signed PWM")
    axes[1, 0].set_ylabel("T (C)")
    axes[2, 0].set_ylabel("residual (C)")
    for ax in axes.flat:
        al.style_axes(ax)
    print(f"figure: {save(fig, 'sim_s18_3_course_100s.png')}")

    # (b) numerical vs analytical, equal chi = 0.46, 300 s, T(0) = T_room
    P, chi, tau, ub = 300.0, 0.46, 80.0, 50.0
    model = Model(Lump(chi, chi, tau), None, dt=0.1)
    ctrl = Controller("PWM_SINE", P, u_a=0, u_b=ub, integer=False, dt=0.1)
    n = int(5 * P / 0.1)
    # numerical: hold u at the mid-interval value (second-order accurate)
    t_mid = (np.arange(n) + 0.5) * 0.1
    pwm = ub * np.sin(2 * math.pi * t_mid / P)
    y = np.zeros(1)
    num = np.empty(n)
    for k in range(n):
        y = model.Phi @ y + model.Gam * model.drive(pwm[k])
        num[k] = y[0] + 22.0
    t = (np.arange(n) + 1) * 0.1
    w = 2 * math.pi / P
    amp = chi * ub / math.sqrt(1 + (w * tau) ** 2)
    delta = math.atan(w * tau)
    Cc = chi * ub * w * tau / (1 + (w * tau) ** 2)
    ana = 22.0 + amp * np.sin(w * t - delta) + Cc * np.exp(-t / tau)
    err = float(np.max(np.abs(num - ana)))
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, (a, b), title in zip(axes, [(0, 2 * P), (4 * P, 5 * P)], ["startup (two cycles)", "settled cycle"]):
        m = (t >= a) & (t <= b)
        ax.plot(t[m], num[m], color=al.INK, lw=2, label="numerical (one-lump, dt 0.1 s)")
        ax.plot(t[m], ana[m], color=al.SENSOR_COLORS["T2_C"], ls="--", lw=1.4, label="analytical (boxed solution)")
        ax.plot(t[m], 22.0 + Cc * np.exp(-t[m] / tau), color=al.MUTED, ls=":", lw=1, label="room + C exp(-t/tau)")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("time (s)")
        al.style_axes(ax)
    axes[0].set_ylabel("T (C)")
    axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle(f"S18.3 check: chi_h = chi_c = 0.46 C/PWM, tau 80 s, u_b 50, 300 s, T(0) = T_room. "
                 f"Max |numerical - analytical| = {err:.1e} C", fontsize=9)
    print(f"300 s equal-chi check: amplitude {amp:.3f} C, lag {math.degrees(delta):.1f} deg "
          f"({delta / w:.1f} s), C = {Cc:.3f} C, max |num - ana| = {err:.2e} C")
    print(f"figure: {save(fig, 'sim_s18_3_analytic300.png')}")

    # (c) 100 s model-experiment comparison with the measured parameters
    d, t_inf = load_raw_window()
    dt = float(np.median(np.diff(d["time_s"])))
    room = t_inf["TEC_C"]
    crow = control_rows()["s18_3_pwm100"]
    t0, t1 = float(crow["t_start_s"]), float(crow["t_end_s"])
    w = 2 * math.pi / 100
    out_rows = []
    m_meas = (d["time_s"] >= t0) & (d["time_s"] <= t1)
    tm = d["time_s"][m_meas]
    meas_T = al.fit_sine(tm - dt / 2, d["TEC_C"][m_meas], w)
    meas_U = al.fit_sine(tm + dt / 2, d["signed_pwm"][m_meas].astype(float), w)
    out_rows.append({"source": "measured", "mean_minus_room_C": meas_T.mean - float(crow["room_C"]),
                     "amp_C": meas_T.amp, "lag_deg": lag_deg(meas_U, meas_T),
                     "resid_amp_C": meas_T.half_pp_resid})
    sims = {}
    # one-lump: start from the measured TEC temperature at the start of the 100 s run
    s34 = d["step"] == RAW_STEPS[100]
    idx0 = np.argmax(s34)
    sub = {k: v[idx0:] for k, v in d.items()}
    for label, model, dd in [("one-lump (measured chi_h, chi_c, tau)", Model(one, None, dt), sub),
                             ("coupled TEC-rod (whole record replayed)", Model(coupled, rod, dt), d)]:
        rec = replay(model, dd, t_inf, room)
        mm = (rec["time_s"] >= t0) & (rec["time_s"] <= t1)
        T = al.fit_sine(rec["time_s"][mm] - dt / 2, rec["TEC_C"][mm], w)
        out_rows.append({"source": label, "mean_minus_room_C": T.mean - float(crow["room_C"]),
                         "amp_C": T.amp, "lag_deg": lag_deg(meas_U, T), "resid_amp_C": T.half_pp_resid})
        sims[label] = rec
    print("S18.3 100 s model-experiment comparison (same recorded PWM, measured room and initial T)")
    print("| Source | Mean - room (C) | Amplitude (C) | Lag behind PWM (deg) | Residual (C) |")
    print("| --- | ---: | ---: | ---: | ---: |")
    for r in out_rows:
        print(f"| {r['source']} | {r['mean_minus_room_C']:.2f} | {r['amp_C']:.2f} | {r['lag_deg']:.1f} | "
              f"{r['resid_amp_C']:.2f} |")
    print(f"table: {write_rows(ANALYSIS / 'sim_s18_3_pwm100_compare.csv', out_rows)}")
    fig, axes = plt.subplots(2, 1, figsize=(10, 5.5), sharex=True, gridspec_kw={"height_ratios": [1, 2]})
    m = (d["time_s"] >= t0 - 300) & (d["time_s"] <= t1)
    axes[0].plot(d["time_s"][m], d["signed_pwm"][m], color=al.REFERENCE_COLOR)
    axes[1].plot(d["time_s"][m], d["TEC_C"][m], color=al.INK, lw=1.5, label="TEC measured")
    styles = [al.SENSOR_COLORS["T2_C"], al.SENSOR_COLORS["T1_C"]]
    for (label, rec), col in zip(sims.items(), styles):
        mm = (rec["time_s"] >= t0 - 300) & (rec["time_s"] <= t1)
        axes[1].plot(rec["time_s"][mm], rec["TEC_C"][mm], color=col, ls="--", label=label)
    for ax in axes:
        ax.axvspan(t0, t1, color=al.GRID, alpha=0.5, lw=0)
        al.style_axes(ax)
    axes[0].set_ylabel("signed PWM")
    axes[1].set_ylabel("TEC temperature (C)")
    axes[1].set_xlabel("Arduino time (s); shaded: fitted interval")
    axes[1].legend(fontsize=7, frameon=False)
    axes[0].set_title("S18.3, 100 s zero-centered PWM: measured TEC vs one-lump and coupled models", fontsize=9)
    print(f"figure: {save(fig, 'sim_s18_3_pwm100_compare.png')}")


# ------------------------------------------------------------------
# S18.4
# ------------------------------------------------------------------
def cmd_s18_4(args):
    rod, coupled, one, p = load_parameters()
    rows_c = control_rows()
    room = float(rows_c["s18_4_pos100"]["room_C"])
    out, store = [], {}
    for drive, label in [((0.0, 50.0), "s18_3_pwm100"), ((55.0, 50.0), "s18_4_pos100")]:
        r = rows_c[label]
        kind = "zero-centered" if drive[0] == 0 else "positive-only"
        out.append({"drive": kind, "source": "measured", "u_a": float(r["u_a"]), "u_b": float(r["u_b"]),
                    "mean_C": float(r["TEC_mean_C"]), "amp_C": float(r["TEC_amp_C"]),
                    "resid_amp_C": float(r["TEC_resid_amp_C"])})
        for mlabel, model, cycles in [("one-lump", Model(one, None, 1.0), 20),
                                      ("coupled TEC-rod", Model(coupled, rod, 1.0), 160)]:
            ctrl = Controller("PWM_SINE", 100.0, u_a=drive[0], u_b=drive[1], integer=True)
            rec, fits, sl = periodic_run(model, room, ctrl, 100.0, cycles)
            T = fits["TEC_C"]
            out.append({"drive": kind, "source": mlabel, "u_a": drive[0], "u_b": drive[1],
                        "mean_C": T.mean, "amp_C": T.amp, "resid_amp_C": T.half_pp_resid})
            store[(kind, mlabel)] = (rec, fits, sl)
    print(f"S18.4 tables, 100 s (room {room:.2f} C)")
    print("| Drive | Source | u_a | u_b | Mean (C) | Fundamental amplitude (C) | Residual amplitude (C) |")
    print("| --- | --- | ---: | ---: | ---: | ---: | ---: |")
    for r in out:
        print(f"| {r['drive']} | {r['source']} | {r['u_a']:.0f} | {r['u_b']:.0f} | {r['mean_C']:.2f} | "
              f"{r['amp_C']:.2f} | {r['resid_amp_C']:.3f} |")
    print(f"table: {write_rows(ANALYSIS / 'sim_s18_4.csv', out)}")

    # aligned plots: measured and simulated, both drives, three settled cycles
    d, t_inf = load_raw_window()
    fig, axes = plt.subplots(3, 2, figsize=(11, 7.5), sharex="col")
    for j, (kind, label) in enumerate([("zero-centered", "s18_3_pwm100"), ("positive-only", "s18_4_pos100")]):
        r = rows_c[label]
        t1 = float(r["t_end_s"])
        m = (d["time_s"] > t1 - 300) & (d["time_s"] <= t1)
        tm = d["time_s"][m] - d["time_s"][m][0]
        meas = al.fit_sine(d["time_s"][m] - 0.5, d["TEC_C"][m], 2 * math.pi / 100)
        axes[0, j].plot(tm, d["signed_pwm"][m], color=al.REFERENCE_COLOR, label="measured PWM")
        axes[1, j].plot(tm, d["TEC_C"][m] - meas.mean, color=al.INK, lw=1.6, label="measured")
        axes[2, j].plot(meas.t - meas.t[0], meas.resid, color=al.INK, lw=1.6, label="measured")
        for mlabel, col in [("one-lump", al.SENSOR_COLORS["T2_C"]), ("coupled TEC-rod", al.SENSOR_COLORS["T1_C"])]:
            rec, fits, sl = store[(kind, mlabel)]
            T, U = fits["TEC_C"], fits["signed_pwm"]
            # align the simulated cycle to the measured PWM phase
            shift = (al.wrap(U.phase - al.fit_sine(d["time_s"][m] + 0.5, d["signed_pwm"][m].astype(float),
                                                    2 * math.pi / 100).phase)) / (2 * math.pi / 100)
            ts = rec["time_s"][sl] + shift
            tt = ts - ts[0] + (d["time_s"][m][0] - d["time_s"][m][0])
            axes[1, j].plot(tt, rec["TEC_C"][sl] - T.mean, color=col, ls="--", label=mlabel)
            axes[2, j].plot(T.t - T.t[0] + shift, T.resid, color=col, ls="--", label=mlabel)
        axes[0, j].set_title(f"{kind} 100 s drive", fontsize=9)
        axes[2, j].set_xlabel("time in a three-cycle settled window (s)")
    axes[0, 0].set_ylabel("signed PWM")
    axes[1, 0].set_ylabel("T - fitted mean (C)")
    axes[2, 0].set_ylabel("residual (C)")
    for ax in axes.flat:
        al.style_axes(ax)
    axes[1, 1].legend(fontsize=7, frameon=False)
    print(f"figure: {save(fig, 'sim_s18_4_aligned.png')}")


# ------------------------------------------------------------------
# FOCUSED CONTROL SIMULATION
# ------------------------------------------------------------------
def control_metrics(rec, fits, sl, ctrl, u_max=255.0):
    T, S = fits["TEC_C"], fits["setpoint_C"]
    err = rec["TEC_C"][sl] - rec["setpoint_C"][sl]
    return {"mean_error_C": T.mean - ctrl.Ta, "amp_ratio": T.amp / ctrl.Tb,
            "lag_deg": lag_deg(S, T), "rms_tracking_C": float(np.sqrt(np.mean(err**2))),
            "saturated_fraction": float(np.mean(np.abs(rec["signed_pwm"][sl]) >= u_max - 0.5)),
            "pwm_min": float(np.min(rec["signed_pwm"][sl])), "pwm_max": float(np.max(rec["signed_pwm"][sl])),
            "T1_amp_C": fits["T1_C"].amp, "T1_mean_C": fits["T1_C"].mean}


def cmd_control(args):
    rod, coupled, one, p = load_parameters()
    rows_c = control_rows()
    cases = [("P only", "P", 100, 6.0, 0.0, "s18_7_p100"),
             ("PI, same Kp", "PI", 100, 6.0, 0.04, "s18_7_pi100"),
             ("PI, targeted adjustment", "PI", 100, 30.0, 0.3, "s18_7_pi100_t2"),
             ("PI selected, 200 s", "PI", 200, 30.0, 0.3, "s18_7_pi200_sel"),
             ("PI selected, 300 s", "PI", 300, 30.0, 0.3, "s18_7_pi300_sel")]
    out, store = [], {}
    for name, mode, P, kp, ki, label in cases:
        mr = rows_c[label]
        room = float(mr["room_C"])
        ctrl = Controller(mode, P, Ta=40, Tb=16, kp=kp, ki=ki, T_room_set=room)
        model = Model(coupled, rod, 1.0)
        rec, fits, sl = periodic_run(model, room, ctrl, P, max(40, int(4 * 3600 / P)), n_fit=5)
        sm = control_metrics(rec, fits, sl, ctrl)
        out.append({"case": name, "mode": mode, "period_s": P, "Kp": kp, "Ki": ki, "source": "simulation",
                    "room_C": room, **sm})
        out.append({"case": name, "mode": mode, "period_s": P, "Kp": kp, "Ki": ki,
                    "source": f"measured {label}", "room_C": room,
                    "mean_error_C": float(mr["mean_error_C"]), "amp_ratio": float(mr["amp_ratio"]),
                    "lag_deg": float(mr["TEC_lag_deg"]), "rms_tracking_C": float(mr["rms_tracking_C"]),
                    "saturated_fraction": float(mr["saturated_fraction"]),
                    "T1_amp_C": float(mr["T1_amp_C"]), "T1_mean_C": float(mr["T1_mean_C"])})
        store[name] = (rec, fits, sl)
    print("Focused control simulation (coupled TEC-rod model, feedback from the TEC node, "
          "Ta 40 C, Tb 16 C) vs the measured runs")
    print("| Case | Period (s) | Kp (PWM/C) | Ki (PWM/(C s)) | Source | Mean error (C) | Amplitude ratio | "
          "Lag (deg) | RMS error (C) | Saturated | T1 amplitude (C) |")
    print("| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in out:
        print(f"| {r['case']} | {r['period_s']} | {r['Kp']:g} | {r['Ki']:g} | {r['source']} | "
              f"{r['mean_error_C']:+.2f} | {r['amp_ratio']:.3f} | {r['lag_deg']:.1f} | "
              f"{r['rms_tracking_C']:.2f} | {100 * r['saturated_fraction']:.0f} % | {r['T1_amp_C']:.2f} |")
    print(f"table: {write_rows(ANALYSIS / 'sim_control.csv', out)}")

    fig, axes = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    colors = [al.SENSOR_COLORS["T2_C"], al.SENSOR_COLORS["T1_C"], al.SENSOR_COLORS["T3_C"]]
    first = True
    for (name, *_), col in zip(cases[:3], colors):
        rec, fits, sl = store[name]
        t = rec["time_s"][sl] - rec["time_s"][sl][0]
        if first:
            axes[0].plot(t, rec["setpoint_C"][sl], color=al.REFERENCE_COLOR, lw=1.2, label="setpoint 40 + 16 sin")
            first = False
        axes[0].plot(t, rec["TEC_C"][sl], color=col, label=name)
        axes[1].plot(t, rec["signed_pwm"][sl], color=col)
    axes[0].set_ylabel("TEC temperature (C)")
    axes[1].set_ylabel("signed PWM")
    axes[1].set_xlabel("time in the settled window (s)")
    axes[0].set_ylim(top=63)
    axes[0].legend(fontsize=7, frameon=False, ncol=4, loc="upper center")
    axes[0].set_title("Focused 100 s simulation: P (Kp 6) vs PI (Kp 6, Ki 0.04) and the adjusted "
                      "PI (Kp 30, Ki 0.3)", fontsize=9)
    for ax in axes:
        al.style_axes(ax)
    print(f"figure: {save(fig, 'sim_control_p_pi100.png')}")


# ------------------------------------------------------------------
# INTEGRATED COMPARISON (Stage 6)
# ------------------------------------------------------------------
def cmd_rod(args):
    rod, coupled, one, p = load_parameters()
    res_path = ANALYSIS / f"angstrom_{args.name}.csv"
    rows = {r["label"]: r for r in csv.DictReader(res_path.open())}
    if args.label not in rows:
        raise SystemExit(f"{args.label} not in {res_path.name} (have {', '.join(rows)})")
    r = rows[args.label]
    iv = al.Interval(args.label, "data/module_08/" + r["file"], step=int(float(r["step"])),
                     t_start=float(r["t_start_s"]), t_end=float(r["t_end_s"]))
    seg = al.select_segment(iv)
    st = seg.settings()
    P, w = seg.period, seg.omega
    Ta, Tb = float(st.get("Ta_C", 40)), float(st.get("Tb_C", 16))
    kp, ki = float(st.get("Kp", 30)), float(st.get("Ki", 0.3))
    # room: the interval's own baseline in intervals_<name>.csv, else --baseline
    spec = args.baseline
    iv_path = DATA / f"intervals_{args.name}.csv"
    if iv_path.exists():
        for other in al.read_intervals(iv_path):
            if other.label == args.label and other.baseline:
                spec = other.baseline
    t_inf = al.baseline_temperatures(spec) if spec else None
    room = (float(np.mean([t_inf[c][0] for c in ROD_COLUMNS])) if t_inf
            else float(st.get("T_room_C", 22)))
    kappa, nu = float(r["kappa_used"]), float(r["nu_used"])
    source = f"kappa, nu fitted to {args.label} itself"
    if args.independent:
        kappa, _, nu, _, used = kappa_nu_from(args.name, exclude=(args.label,))
        source = f"kappa, nu from {', '.join(used)} only (independent of {args.label})"
    rod = Rod(kappa, nu, rod.contact_m, rod.x1_m)
    print(f"room {room:.2f} C from {spec}; {source}: kappa {kappa:.4e} m^2/s, nu {nu:.3e} 1/s")
    dd = seg.data
    t = dd["time_s"]
    sel, n = al.complete_cycles(t, P)
    dt = float(np.median(np.diff(t[sel])))
    meas = {c: al.fit_sine(t[sel] - dt / 2, dd[c][sel], w) for c in TEMPERATURE_COLUMNS}
    meas_S = al.fit_sine(t[sel], dd["setpoint_C"][sel], w)

    ctrl = Controller(seg.mode, P, Ta=Ta, Tb=Tb, kp=kp, ki=ki, T_room_set=room, dt=dt)
    model = Model(coupled, rod, dt)
    rec, fits, sl = periodic_run(model, room, ctrl, P, max(30, int(5 * 3600 / P)), n_fit=n)
    S = fits["setpoint_C"]
    table = []
    for c in TEMPERATURE_COLUMNS:
        M, F = meas[c], fits[c]
        table.append({"label": args.label, "period_s": P, "sensor": SENSOR_LABELS[c],
                      "meas_mean_C": M.mean, "sim_mean_C": F.mean,
                      "meas_amp_C": M.amp, "meas_amp_sigma": M.sigma_amp, "sim_amp_C": F.amp,
                      "meas_lag_deg": lag_deg(meas_S, M), "meas_lag_sigma_deg": math.degrees(M.sigma_phase),
                      "sim_lag_deg": lag_deg(S, F)})
    print(f"Integrated comparison, {args.label} ({P:.0f} s, {seg.mode} Kp {kp:g} Ki {ki:g}, room {room:.2f} C; "
          f"lags relative to the requested setpoint)")
    print("| Sensor | Mean meas / sim (C) | Amplitude meas / sim (C) | Lag meas / sim (deg) |")
    print("| --- | ---: | ---: | ---: |")
    for row in table:
        print(f"| {row['sensor']} | {row['meas_mean_C']:.2f} / {row['sim_mean_C']:.2f} | "
              f"{row['meas_amp_C']:.2f} / {row['sim_amp_C']:.2f} | {row['meas_lag_deg']:.1f} / "
              f"{row['sim_lag_deg']:.1f} |")
    tag = f"{args.label}_independent" if args.independent else args.label
    for row in table:
        row["kappa_nu_source"] = source
    print(f"table: {write_rows(ANALYSIS / f'sim_rod_{tag}.csv', table)}")

    # analytical rod prediction from the measured T1 fit and the same kappa, nu
    pos = al.positions_for(seg.meta)
    # time axis: phase of the requested setpoint = 0 at t = 0 for both
    t_meas = t[sel] - (-meas_S.phase) / w
    t0m = t_meas[0]
    t_sim = rec["time_s"][sl] - (-S.phase) / w
    shift_sim = math.floor((t_sim[0] - t0m) / P) * P
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1.6, 1]})
    ax = axes[0]
    ax.plot(t_meas - t0m, dd["setpoint_C"][sel], color=al.REFERENCE_COLOR, lw=1, label="setpoint")
    ax.plot(t_meas - t0m, dd["TEC_C"][sel], color=al.INK, lw=1.5, label="TEC measured")
    ax.plot(t_sim - shift_sim - t0m, rec["TEC_C"][sl], color=al.INK, lw=1, ls="--", label="TEC simulated")
    ax.plot(t_meas - t0m, dd["T1_C"][sel], color=al.SENSOR_COLORS["T1_C"], lw=1.5, label="T1 measured")
    ax.plot(t_sim - shift_sim - t0m, rec["T1_C"][sl], color=al.SENSOR_COLORS["T1_C"], lw=1, ls="--",
            label="T1 simulated")
    ax.set_ylabel("T (C)")
    ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.set_title(f"{args.label}: boundary. Simulated rod base vs measured T1 first "
                 f"(both aligned on the requested setpoint phase); {source}", fontsize=8)
    ax = axes[1]
    for c in ["T1_C", "T3_C", "T5_C"]:
        col = al.SENSOR_COLORS[c]
        ax.plot(t_meas - t0m, dd[c][sel], color=col, lw=1.6, label=f"{SENSOR_LABELS[c]} measured")
        ax.plot(t_sim - shift_sim - t0m, rec[c][sl], color=col, lw=1, ls="--",
                label=f"{SENSOR_LABELS[c]} numerical model")
        if t_inf is not None and c != "T1_C":
            x = pos[c] / 100
            base = meas["T1_C"]
            pred = al.predict_rod(x, t[sel] - dt / 2, t_inf[c][0], base, kappa, nu, w)
            ax.plot(t_meas - t0m, pred, color=col, lw=1, ls=":",
                    label=f"{SENSOR_LABELS[c]} analytical (measured T1, m = sqrt(nu/kappa))")
    ax.set_ylabel("T (C)")
    ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    # residuals: measured minus model, the periodic simulation interpolated onto the
    # measured times
    ax = axes[2]
    xs = t_sim - shift_sim - t0m
    xm = t_meas - t0m
    for c in ["TEC_C", "T1_C", "T3_C", "T5_C"]:
        sim_at = np.interp(xm, xs, rec[c][sl], period=P)
        ax.plot(xm, dd[c][sel] - sim_at, color=al.SENSOR_COLORS[c], lw=1,
                label=f"{SENSOR_LABELS[c]} measured - numerical")
        if t_inf is not None and c in ("T3_C", "T5_C"):
            pred = al.predict_rod(pos[c] / 100, t[sel] - dt / 2, t_inf[c][0], meas["T1_C"],
                                  kappa, nu, w)
            ax.plot(xm, dd[c][sel] - pred, color=al.SENSOR_COLORS[c], lw=1, ls=":",
                    label=f"{SENSOR_LABELS[c]} measured - analytical")
    ax.axhline(0, color=al.MUTED, lw=0.8)
    ax.set_ylabel("residual (C)")
    ax.set_xlabel("time from a setpoint zero crossing (s); last three fitted cycles")
    ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    t_end = float(t_meas[-1] - t0m)
    axes[1].set_xlim(t_end - 3 * P, t_end)
    for a in axes:
        al.style_axes(a)
    print(f"figure: {save(fig, f'sim_rod_{tag}.png')}")


# ------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", nargs="?", default="all",
                    choices=["all", "calibrate", "s18-3", "s18-4", "control", "rod"])
    ap.add_argument("label", nargs="?", help="interval label for 'rod' (e.g. p300)")
    ap.add_argument("--name", default="production", help="angstrom_<name>.csv with kappa, nu")
    ap.add_argument("--baseline", default=None,
                    help="FILE:START,END room baseline for 'rod' if the interval has none in "
                         "intervals_<name>.csv")
    ap.add_argument("--independent", action="store_true",
                    help="rod: take kappa and nu from the other intervals only")
    args = ap.parse_args(argv)
    if args.command in ("all", "calibrate"):
        cmd_calibrate(args)
    if args.command in ("all", "s18-3"):
        cmd_s18_3(args)
    if args.command in ("all", "s18-4"):
        cmd_s18_4(args)
    if args.command in ("all", "control"):
        cmd_control(args)
    if args.command == "rod" or (args.command == "all" and args.label):
        if not args.label:
            raise SystemExit("rod needs an interval label, e.g. p300")
        cmd_rod(args)


if __name__ == "__main__":
    main()
