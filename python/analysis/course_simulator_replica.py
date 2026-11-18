# Phys 39 Module 8: a Python replica of the course's six-sensor TEC-rod simulator
#
# A line-by-line port of the model in
#   https://sethfraden.github.io/Phys39F26-course/simulators/angstrom-six-sensor/
# (page source checked on 18 Nov 2026; SHA-256 in data/module_08/course_simulator/
# course_sim_meta.json), so the course simulations can be run and plotted fully
# automatically, without a browser:
#   - coupled model: TEC1-12706 Peltier plate (Seebeck, Joule, conduction, heat sink),
#     contact to the rod, explicit finite-difference rod with natural convection,
#     six thermistors with a 2 s first-order lag; 0.01 s physics substeps
#   - PID on T0 every 0.1 s (% duty, 0.8 s derivative filter, conditional anti-windup,
#     integral clamp +-120 %, output clamp +-100 %), Angstrom sinusoid or constant setpoint
#   - open-loop PWM (constant or sinusoidal, clipped at +-255 counts)
#   - one-lump mode (S18.3): RK4 of tau dT/dt = Troom + chi(u) u - T, analytical overlay
#   - the page's chart buffer (one point per 0.25 s, three periods or 2000 s), its
#     "Settled response" readout and settled residual view, and its Save CSV format
# Same constants, same order of operations, same float64 arithmetic. The frost layer is
# not ported: on the page it only changes the drawing, never a temperature.
#
# The CLI runs the Module 8 cases (S18.3 table and startup check, S18.4, focused S19 P/PI,
# Stage 6 at 100-500 s), writes CSVs, draws figures in the page's style, and checks the
# replica against the browser records when they exist.
#
# Usage (from the repository root):
#   python python/analysis/course_simulator_replica.py              # everything, about 2 min:
#       runs, CSVs (data/module_08/course_simulator_replica/), page-style figures
#       (docs/figures/module_08/course_replica/), the check against the browser runs,
#       and the comparison tables and figures of module_08_course_sim.py (tag "replica")
#   python python/analysis/course_simulator_replica.py --quick      # shorter coupled runs (test)

import argparse
import json
import math
import sys
from dataclasses import dataclass, field
from multiprocessing import Pool
from pathlib import Path

import numpy as np

# ------------------------------------------------------------------
# CONSTANTS (as on the page)
# ------------------------------------------------------------------
MATERIALS = {
    "aluminum": dict(name="Aluminum 6061", k=167, rho=2700, c=896),
    "copper": dict(name="Copper", k=390, rho=8960, c=385),
    "brass": dict(name="Brass", k=110, rho=8500, c=380),
    "steel": dict(name="Mild steel", k=50, rho=7850, c=490),
    "stainless": dict(name="Stainless 304", k=16, rho=8000, c=500),
    "acrylic": dict(name="Acrylic", k=0.19, rho=1180, c=1470),
}
TAMB, DEW = 22.0, 11.0
TEC = dict(S=0.0513, R=1.861, Kt=0.560, Vsup=12.0, Imax=8.0,
           Cc=12.0, Ch=180.0, Rhs=0.55, Gleak=0.15, hContact=15000)
H_CONV = 10.0
SUB_DT = 0.01
PID_DT = 0.1
SENSOR_FRACS = [0, 0.12, 0.28, 0.48, 0.72, 0.95]
TAU_SENSOR = 2.0


def node_count(Lmm):
    return max(24, min(100, int(math.floor(Lmm / 2.5 + 0.5))))     # Math.round


def js_round(x):
    return math.floor(x + 0.5)


@dataclass
class Drive:                    # open-loop PWM
    manual: bool = False
    raw: float = 0.0
    rev: bool = False
    wave: str = "constant"
    offset: float = 0.0
    amplitude: float = 50.0
    period: float = 300.0


@dataclass
class Lump:
    mode: str = "coupled"
    room: float = 22.0
    initial: float = 22.0
    tau: float = 80.0
    heat: float = 0.46
    cool: float = 0.46


@dataclass
class Pid:
    sp: float = 10.0
    kp: float = 8.0
    ki: float = 0.15
    kd: float = 15.0
    pOn: bool = True
    iOn: bool = True
    dOn: bool = False
    integ: float = 0.0
    dFilt: float = 0.0
    prevPV: float = None
    u: float = 0.0


@dataclass
class SetpointDrive:
    mode: str = "sine"
    constant: float = 10.0
    Ta: float = 40.0
    Tb: float = 16.0
    period: float = 200.0


@dataclass
class Chart:
    t: list = field(default_factory=list)
    temps: list = field(default_factory=lambda: [[] for _ in SENSOR_FRACS])
    sp: list = field(default_factory=list)
    u: list = field(default_factory=list)
    tec: list = field(default_factory=list)
    exact: list = field(default_factory=list)
    max: int = 8000


class Simulator:
    """State and update rules of the page (sim, lump, drive, pid, setpointDrive, chart)."""

    def __init__(self, material="aluminum", d_mm=12.7, L_mm=300):
        self.matKey = material
        self.d = d_mm / 1000
        self.L = L_mm / 1000
        self.time = 0.0
        self.Tc = TAMB
        self.Th = TAMB
        self.Ts = np.full(6, TAMB)
        self.duty = 0.0
        self.V = self.I = self.W = 0.0
        self.pidAcc = 0.0
        self.chartAcc = 0.0
        self.lump = Lump()
        self.drive = Drive()
        self.pid = Pid()
        self.spd = SetpointDrive()
        self.chart = Chart()
        self.rebuild_rod()
        self.Ts[:] = TAMB

    # --- geometry -------------------------------------------------
    def rebuild_rod(self):
        self.N = node_count(self.L * 1000)
        self.T = np.full(self.N, TAMB)
        self.dx = self.L / self.N
        # sensor interpolation (rodTempAt) for the six fixed fractions
        self._sj, self._sf = self._interp_index(np.array(SENSOR_FRACS, float))

    def _interp_index(self, fracs):
        x = fracs * (self.N - 1)
        j = np.minimum(self.N - 2, np.floor(x)).astype(int)
        return j, x - j

    def rod_temp_at(self, frac):
        x = frac * (self.N - 1)
        j = min(self.N - 2, math.floor(x))
        f = x - j
        return self.T[j] * (1 - f) + self.T[j + 1] * f

    def rod_temps_at(self, fracs):
        j, f = self._interp_index(np.asarray(fracs, float))
        return self.T[j] * (1 - f) + self.T[j + 1] * f

    # --- drives ---------------------------------------------------
    def pwm_at(self, t):
        dr = self.drive
        if dr.wave == "sine":
            command = dr.offset + dr.amplitude * math.sin(2 * math.pi * t / dr.period)
        else:
            command = (-1 if dr.rev else 1) * dr.raw
        return max(-255.0, min(255.0, command))

    def analytical_at(self, t):
        lp, dr = self.lump, self.drive
        if (lp.mode != "lumped" or not dr.manual or dr.wave != "sine" or lp.heat != lp.cool
                or abs(dr.offset) + dr.amplitude > 255):
            return float("nan")
        r = 2 * math.pi * lp.tau / dr.period
        mean = lp.room + lp.heat * dr.offset
        a = lp.heat * dr.amplitude / math.sqrt(1 + r * r)
        phase = math.atan(r)
        return (mean + a * math.sin(2 * math.pi * t / dr.period - phase)
                + (lp.initial - mean + a * math.sin(phase)) * math.exp(-t / lp.tau))

    def update_setpoint(self):
        if self.spd.mode == "sine":
            omega = 2 * math.pi / self.spd.period
            self.pid.sp = self.spd.Ta + self.spd.Tb * math.sin(omega * self.time)
        else:
            self.pid.sp = self.spd.constant

    def pid_step(self, dt):
        p = self.pid
        self.update_setpoint()
        pv = self.Ts[0]
        e = p.sp - pv
        dpv = 0.0
        if p.prevPV is not None:
            raw = (pv - p.prevPV) / dt
            a = dt / (0.8 + dt)
            p.dFilt += a * (raw - p.dFilt)
            dpv = p.dFilt
        p.prevPV = pv
        P = p.kp * e if p.pOn else 0.0
        D = -p.kd * dpv if p.dOn else 0.0
        if p.iOn:
            u_unsat = P + p.integ + D
            if not ((u_unsat >= 100 and e > 0) or (u_unsat <= -100 and e < 0)):
                p.integ += p.ki * e * dt
            p.integ = max(-120.0, min(120.0, p.integ))
        Iterm = p.integ if p.iOn else 0.0
        p.u = max(-100.0, min(100.0, P + Iterm + D))
        return p.u

    def current_duty(self):
        if self.drive.manual:
            return self.pwm_at(self.time) / 255
        return self.pid.u / 100

    # --- physics --------------------------------------------------
    def lump_step(self, h):
        lp = self

        def f(t, T):
            u = self.pwm_at(t)
            return (self.lump.room + (self.lump.heat if u >= 0 else self.lump.cool) * u - T) / self.lump.tau
        t, T = lp.time, lp.Tc
        k1 = f(t, T)
        k2 = f(t + h / 2, T + h * k1 / 2)
        k3 = f(t + h / 2, T + h * k2 / 2)
        k4 = f(t + h, T + h * k3)
        lp.Tc += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        lp.Ts[0] = lp.Tc

    def physics_step(self, h):
        m = MATERIALS[self.matKey]
        S, R, Kt, Vsup, Imax = TEC["S"], TEC["R"], TEC["Kt"], TEC["Vsup"], TEC["Imax"]
        Cc, Ch, Rhs, Gleak, hContact = TEC["Cc"], TEC["Ch"], TEC["Rhs"], TEC["Gleak"], TEC["hContact"]
        A = math.pi * self.d * self.d / 4
        Gc = hContact * A
        alpha = m["k"] / (m["rho"] * m["c"])
        convC = 4 * H_CONV / (m["rho"] * m["c"] * self.d)
        dx2 = self.dx * self.dx
        T = self.T
        cap = m["rho"] * m["c"] * A * self.dx
        I = V = 0.0
        if abs(self.duty) > 0.002:
            V = self.duty * Vsup
            I = (V - (-S * (self.Th - self.Tc))) / R
            I = max(-Imax, min(Imax, I))
        TcK, ThK = self.Tc + 273.15, self.Th + 273.15
        Ic = -I
        Qc = S * Ic * TcK - 0.5 * Ic * Ic * R - Kt * (self.Th - self.Tc)
        Qh = S * Ic * ThK + 0.5 * Ic * Ic * R - Kt * (self.Th - self.Tc)
        self.Tc += h * (-Qc + Gc * (T[0] - self.Tc) + Gleak * (TAMB - self.Tc)) / Cc
        self.Th += h * (Qh - (self.Th - TAMB) / Rhs) / Ch
        self.V, self.I, self.W = V, I, abs(Qh - Qc)
        # rod conduction (explicit); Tl = T[0] at i = 0 and Tr = T[N-1] at i = N-1
        lap = np.empty_like(T)
        lap[1:-1] = T[:-2] - 2 * T[1:-1] + T[2:]
        lap[0] = T[0] - 2 * T[0] + T[1]
        lap[-1] = T[-2] - 2 * T[-1] + T[-1]
        dT = alpha * lap / dx2 - convC * (T - TAMB)
        dT[0] += Gc * (self.Tc - T[0]) / cap
        self.T = T + h * dT
        # six thermistors, first-order lag
        vals = self.T[self._sj] * (1 - self._sf) + self.T[self._sj + 1] * self._sf
        self.Ts += h * (vals - self.Ts) / TAU_SENSOR

    # --- chart ----------------------------------------------------
    def chart_window(self):
        if self.drive.manual:
            return 3 * self.drive.period if self.drive.wave == "sine" else 2000
        return 3 * self.spd.period if self.spd.mode == "sine" else 2000

    def chart_push(self):
        c = self.chart
        c.t.append(self.time)
        for i, a in enumerate(c.temps):
            a.append(float(self.Ts[i]))
        c.tec.append(self.Tc)
        c.exact.append(self.analytical_at(self.time))
        c.sp.append(self.pid.sp)
        c.u.append(self.duty * 100)
        c.max = math.ceil(max(2000, self.chart_window()) / 0.25) + 1
        if len(c.t) > c.max:
            for a in [c.tec, c.exact, c.t, c.sp, c.u] + c.temps:
                a.pop(0)

    # --- reset and main loop -------------------------------------
    def reset(self):
        self.time = 0.0
        self.Tc = self.Th = TAMB
        self.Ts[:] = TAMB
        self.duty = 0.0
        self.V = self.I = self.W = 0.0
        if self.lump.mode == "lumped":
            self.Tc = self.lump.initial
            self.Ts[:] = self.lump.initial
        self.pidAcc = self.chartAcc = 0.0
        self.rebuild_rod()
        self.pid.u = 0.0
        self.pid.integ = 0.0
        self.pid.dFilt = 0.0
        self.pid.prevPV = None
        self.update_setpoint()
        self.chart = Chart()
        self.duty = self.current_duty()
        self.chart_push()

    def advance(self, sim_dt):
        remaining, guard = sim_dt, 0
        lumped = self.lump.mode == "lumped"
        while remaining > 1e-9 and guard < 12000:
            h = min(SUB_DT, remaining)
            self.pidAcc += h
            self.chartAcc += h
            if self.pidAcc >= PID_DT:
                dt = self.pidAcc
                self.pidAcc = 0.0
                if not self.drive.manual:
                    self.pid_step(dt)
                self.duty = self.current_duty()
            if self.drive.manual:
                self.duty = self.pwm_at(self.time + h / 2) / 255
            if lumped:
                self.lump_step(h)
            else:
                self.physics_step(h)
            self.time += h
            if self.drive.manual:
                self.duty = self.pwm_at(self.time) / 255
            if self.chartAcc >= 0.25:
                self.chartAcc = 0.0
                self.chart_push()
            remaining -= h
            guard += 1

    # --- page readouts --------------------------------------------
    def settled_readout(self):
        """updateSettledReadout(): one complete cycle of the chart's TEC trace, 360
        uniform samples; returns None while waiting, else a dict with the page's numbers
        and the samples for the residual view."""
        c, period = self.chart, self.drive.period
        end = math.floor(c.t[-1] / period) * period
        start = end - period
        if len(c.t) < 2 or start < 5 * self.lump.tau or c.t[0] > start:
            return None
        count = 360
        samples = []
        mean = sine = cosine = 0.0
        j = max(1, next(i for i, t in enumerate(c.t) if t >= start))
        for i in range(count + 1):
            t = start + i * period / count
            while c.t[j] < t:
                j += 1
            fraction = (t - c.t[j - 1]) / (c.t[j] - c.t[j - 1])
            temperature = c.tec[j - 1] + fraction * (c.tec[j] - c.tec[j - 1])
            angle = 2 * math.pi * i / count
            samples.append(temperature)
            if i == count:
                break
            mean += temperature
            sine += temperature * math.sin(angle)
            cosine += temperature * math.cos(angle)
        amplitude = 2 * math.hypot(sine, cosine) / count
        lag = math.degrees(math.atan2(-cosine, sine))
        fitted = [(mean + 2 * sine * math.sin(2 * math.pi * i / count)
                   + 2 * cosine * math.cos(2 * math.pi * i / count)) / count for i in range(count + 1)]
        resid = [s - f for s, f in zip(samples, fitted)]
        rms = math.sqrt(sum(r * r for r in resid[:-1]) / count)
        pwm = [self.pwm_at(start + i * period / count) for i in range(count + 1)]
        return dict(period=period, mean=mean / count, mean_minus_room=mean / count - self.lump.room,
                    amplitude=amplitude, lag_deg=lag, resid_amp=(max(resid) - min(resid)) / 2,
                    resid_rms=rms, start=start, samples=samples, fitted=fitted, resid=resid, pwm=pwm)

    def csv_text(self):
        """Save CSV, as the page writes it."""
        one = self.lump.mode == "lumped"
        if self.drive.manual:
            drive = vars(self.drive).copy()
        else:
            drive = {**vars(self.spd), **{k: v for k, v in vars(self.pid).items()}}
        meta = {"note": "Settings at export; coupled geometry and PID controls may have changed during the record.",
                "model": self.lump.mode, "drive": drive,
                "lump": vars(self.lump) if one else None,
                "rod": None if one else {"material": self.matKey, "length_m": self.L, "diameter_m": self.d,
                                         "sensor_fractions": SENSOR_FRACS},
                "ambient_C": self.lump.room if one else TAMB, "source": "Python replica"}
        rows = ["# " + json.dumps(meta),
                "time_s,pwm_signed_counts,duty_percent,TEC_C,analytical_TEC_C,setpoint_C,T0_C,T1_C,T2_C,T3_C,T4_C,T5_C"]
        c = self.chart
        for i in range(len(c.t)):
            ex = "" if not np.isfinite(c.exact[i]) else repr(c.exact[i])
            sp = "" if self.drive.manual else repr(c.sp[i])
            temps = ["" if one else repr(a[i]) for a in c.temps]
            rows.append(",".join([repr(c.t[i]), repr(c.u[i] * 255 / 100), repr(c.u[i]), repr(c.tec[i]), ex, sp] + temps))
        return "\n".join(rows)


# ------------------------------------------------------------------
# SET-UP HELPERS (what a student does on the page)
# ------------------------------------------------------------------
def coupled(period, kp_pct, ki_pct, L_mm=285, d_mm=12.7, Ta=40, Tb=16, kd_pct=15, d_on=False):
    s = Simulator(d_mm=d_mm, L_mm=L_mm)
    s.lump.mode = "coupled"
    s.drive.manual = False
    s.spd.mode = "sine"
    s.spd.Ta, s.spd.Tb, s.spd.period = Ta, Tb, period
    s.pid.kp, s.pid.ki, s.pid.kd = kp_pct, ki_pct, kd_pct
    s.pid.pOn, s.pid.iOn, s.pid.dOn = True, ki_pct > 0, d_on
    s.reset()
    return s


def lumped(period, ua, ub, room=22, initial=22, tau=80, heat=0.46, cool=0.46):
    s = Simulator()
    s.lump = Lump("lumped", room, initial, tau, heat, cool)
    s.drive = Drive(manual=True, wave="sine", offset=ua, amplitude=ub, period=period)
    s.reset()
    return s


# ------------------------------------------------------------------
# MODULE 8 CASES
# ------------------------------------------------------------------
REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "data" / "module_08" / "course_simulator_replica"
FIG = REPO / "docs" / "figures" / "module_08" / "course_replica"
BROWSER = REPO / "data" / "module_08" / "course_simulator"
PWM_PER_PERCENT = 2.55
L_MM = 285                                         # slider step 5 mm; our rod 286.6 mm
HOLES_MM = [17.6 + x for x in (0, 24, 50, 76, 102)]  # holes 1, 3, 5, 6, 7 from the TEC face
ONE = dict(tau=59.1, heat=0.4115, cool=0.1732, room=21.73)   # measured one-lump (notes 8.0)

COUPLED = [  # label, period, Kp (PWM/C), Ki (PWM/(C s)), simulated s
    ("s19_p100", 100, 6, 0, 14400), ("s19_pi100", 100, 6, 0.04, 14400),
    ("s19_pi100_adj", 100, 30, 0.3, 14400), ("pi200", 200, 30, 0.3, 14400),
    ("pi300", 300, 30, 0.3, 18000), ("pi400", 400, 30, 0.3, 18000), ("pi500", 500, 30, 0.3, 18000)]
LUMPED = [(f"s18_3_P{P}_chic{cc}", P, 0, 50, 22, 22, 80, 0.46, cc, 12 * P, True)
          for P in (100, 200, 300) for cc in (0.23, 0.46)]
LUMPED += [("s18_3_P300_equal_startup", 300, 0, 50, 22, 22, 80, 0.46, 0.46, 600, False),
           ("s18_4_zero100", 100, 0, 50, ONE["room"], ONE["room"], ONE["tau"], ONE["heat"], ONE["cool"], 2000, True),
           ("s18_4_pos100", 100, 55, 50, ONE["room"], ONE["room"], ONE["tau"], ONE["heat"], ONE["cool"], 2000, True)]


def run_coupled(args):
    label, P, kp, ki, total = args
    s = coupled(P, kp / PWM_PER_PERCENT, ki / PWM_PER_PERCENT, L_mm=L_MM)
    fr = [x / (s.L * 1000) for x in HOLES_MM]
    rec = np.empty((total, 10))
    snap = None
    for k in range(total):
        s.advance(1.0)
        rec[k] = [s.time, s.pid.sp, s.pid.u, s.Tc, s.Ts[0], *s.rod_temps_at(fr)]
    snap = dict(T=s.T.copy(), Tc=s.Tc, Th=s.Th, Ts=s.Ts.copy(), duty=s.duty, chart=s.chart,
                N=s.N, L=s.L, d=s.d, time=s.time, kp=s.pid.kp, ki=s.pid.ki, iOn=s.pid.iOn,
                Ta=s.spd.Ta, Tb=s.spd.Tb, period=s.spd.period, csv=s.csv_text())
    return label, rec, snap


def run_lumped(args):
    label, P, ua, ub, room, init, tau, heat, cool, total, residual = args
    s = lumped(P, ua, ub, room, init, tau, heat, cool)
    rec = []
    for k in range(int(round(total / 0.1))):
        s.advance(0.1)
        rec.append([s.time, s.pwm_at(s.time), s.Tc, s.analytical_at(s.time)])
    return label, np.array(rec), dict(readout=s.settled_readout(), chart=s.chart, csv=s.csv_text(),
                                      lump=vars(s.lump).copy(), drive=vars(s.drive).copy(), time=s.time,
                                      residual=residual)


# ------------------------------------------------------------------
# FIGURES IN THE PAGE'S STYLE
# ------------------------------------------------------------------
BG, PANEL, EDGE, TXT, DIM = "#14171c", "#1d2129", "#323a48", "#dde3ec", "#8b95a7"
SENSOR_COLORS = ['#4fb3ff', '#57d98a', '#ffd166', '#ff8c61', '#c792ea', '#5ac8fa']
IRON = [[0, 0, 10], [32, 0, 140], [145, 0, 156], [226, 43, 0], [255, 143, 0], [255, 217, 0], [255, 255, 255]]


def _style(ax):
    ax.set_facecolor(PANEL)
    for sp in ax.spines.values():
        sp.set_color(EDGE)
    ax.tick_params(colors=DIM, labelsize=8)
    ax.grid(color="#505a6e", alpha=0.25, lw=0.6)


def _fmt_time(s):
    return f"{s:.0f}s" if s < 120 else f"{s / 60:.1f}m"


def _ironbow():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("ironbow", [np.array(c) / 255 for c in IRON])


def figure_coupled(label, rec, snap, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    P = snap["period"]
    fig = plt.figure(figsize=(13, 9.2), facecolor=BG)
    gs = fig.add_gridspec(4, 1, height_ratios=[1.0, 0.32, 1.6, 1.25], hspace=0.42)
    fig.text(0.012, 0.975, "TEC Rod Thermal Lab  (Python replica of the course simulator)", color=TXT,
             fontsize=12, weight="bold")
    fig.text(0.012, 0.952, f"Coupled TEC + rod  ·  Aluminum 6061 ⌀{snap['d'] * 1000:.1f} mm × {snap['L'] * 1000:.0f} mm  ·  "
             f"PID on T0: Kp {snap['kp']:.3f} %/°C, Ki {snap['ki']:.4f} %/(°C·s)"
             f"{'' if snap['iOn'] else ' (I off)'}  ·  Tset = {snap['Ta']:.0f} + {snap['Tb']:.0f} sin(2πt/{P:.0f} s)"
             f"  ·  t = {snap['time'] / 60:.1f} min  ·  ambient {TAMB:.1f} °C", color=DIM, fontsize=9)
    # thermal-cam strip
    ax = fig.add_subplot(gs[0])
    ax.set_facecolor("#000000")
    lo, hi = min(TAMB, snap["Ta"] - snap["Tb"]), max(TAMB, snap["Ta"] + snap["Tb"])
    cmap = _ironbow()
    Lmm = snap["L"] * 1000
    ax.imshow(snap["T"][None, :], aspect="auto", cmap=cmap, vmin=lo, vmax=hi,
              extent=[0, Lmm, -6.35, 6.35])
    ax.add_patch(plt.Rectangle((-9, -14), 9, 28, color=cmap((snap["Tc"] - lo) / (hi - lo))))
    ax.add_patch(plt.Rectangle((-20, -22), 11, 44, color=cmap((snap["Th"] - lo) / (hi - lo))))
    for i, fr in enumerate(SENSOR_FRACS):
        x = fr * Lmm
        ax.plot([x], [0], "o", color="white", ms=4)
        ax.text(x, 9 + 4 * (i % 2), f"T{i}", color="white", ha="center", fontsize=8)
    for k, x in enumerate(HOLES_MM):
        ax.plot([x, x], [-6.35, -11], color="#ff6a4d", lw=1.2)
        ax.text(x, -16, f"H{[1, 3, 5, 6, 7][k]}", color="#ff6a4d", ha="center", fontsize=7)
    ax.set_xlim(-25, Lmm + 10)
    ax.set_ylim(-24, 24)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(-22, 19, f"IRONBOW  •  {lo:.0f}…{hi:.0f} °C   (red ticks: our holes 1, 3, 5, 6, 7)", color="white", fontsize=8)
    for sp in ax.spines.values():
        sp.set_color(EDGE)
    # sensor tiles
    ax = fig.add_subplot(gs[1])
    ax.set_axis_off()
    for i, fr in enumerate(SENSOR_FRACS):
        x0 = i / 6 + 0.003
        ax.add_patch(plt.Rectangle((x0, 0.05), 1 / 6 - 0.008, 0.9, transform=ax.transAxes, color=PANEL, ec=EDGE))
        ax.add_patch(plt.Rectangle((x0, 0.87), 1 / 6 - 0.008, 0.08, transform=ax.transAxes, color=SENSOR_COLORS[i]))
        ax.text(x0 + 0.01, 0.55, f"T{i}", color=DIM, fontsize=8, transform=ax.transAxes, weight="bold")
        ax.text(x0 + 1 / 6 - 0.02, 0.55, f"{fr * Lmm:.0f} mm", color=DIM, fontsize=8, ha="right", transform=ax.transAxes)
        ax.text(x0 + 0.01, 0.15, f"{snap['Ts'][i]:.1f} °C", color=TXT, fontsize=11, family="monospace",
                transform=ax.transAxes, weight="bold")
    # page chart: last three periods of the chart buffer
    c = snap["chart"]
    t = np.array(c.t)
    span = 3 * P
    t1 = max(span, t[-1])
    m = t >= t1 - span
    ax = fig.add_subplot(gs[2])
    _style(ax)
    for i, a in enumerate(c.temps):
        ax.plot(t[m], np.array(a)[m], color=SENSOR_COLORS[i], lw=1.6, label=f"T{i}")
    ax.plot(t[m], np.array(c.tec)[m], color="white", lw=1.6, label="TEC")
    ax.plot(t[m], np.array(c.sp)[m], color="#ff6a4d", lw=1.6, ls="--", label="SP")
    ax.axhline(TAMB, color="#a7c4d9", ls=":", lw=1.1)
    ax.text(t1 - span + 5, TAMB + 0.6, f"Room {TAMB:.1f} °C", color="#a7c4d9", fontsize=8)
    ax2 = ax.twinx()
    ax2.plot(t[m], np.array(c.u)[m], color="#ffd166", alpha=0.55, lw=1)
    ax2.set_ylim(-100, 100)
    ax2.tick_params(colors="#7a6a3a", labelsize=8)
    ax2.set_ylabel("PWM % (right)", color="#b7a36b", fontsize=8)
    ax.set_xlim(t1 - span, t1)
    from matplotlib.ticker import FuncFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: _fmt_time(v)))
    ax.set_ylabel("°C", color=DIM)
    ax.legend(ncol=9, fontsize=7.5, frameon=False, labelcolor=TXT, loc="upper left", bbox_to_anchor=(0, 1.16))
    # our hole positions, same window
    ax = fig.add_subplot(gs[3])
    _style(ax)
    tt = rec[:, 0]
    mm = tt >= rec[-1, 0] - span
    ax.plot(tt[mm], rec[mm, 1], color="#ff6a4d", ls="--", lw=1.2, label="SP")
    ax.plot(tt[mm], rec[mm, 4], color=SENSOR_COLORS[0], lw=1.4, label="T0 (feedback, = our TEC sensor)")
    hole_cols = ["#57d98a", "#ffd166", "#ff8c61", "#c792ea", "#5ac8fa"]
    for k in range(5):
        ax.plot(tt[mm], rec[mm, 5 + k], color=hole_cols[k], lw=1.4,
                label=f"hole {[1, 3, 5, 6, 7][k]} ({HOLES_MM[k]:.1f} mm)")
    ax.set_xlim(rec[-1, 0] - span, rec[-1, 0])
    from matplotlib.ticker import FuncFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: _fmt_time(v)))
    ax.set_ylabel("°C", color=DIM)
    ax.set_title("Rod temperature at our sensor positions (TEC + holes 1, 3, 5, 6, 7)", color=TXT, fontsize=9, loc="left")
    ax.legend(ncol=4, fontsize=7.5, frameon=False, labelcolor=TXT, loc="upper center",
              bbox_to_anchor=(0.5, -0.12))
    fig.savefig(path, dpi=150, facecolor=BG, bbox_inches="tight")
    plt.close(fig)


def figure_lumped(label, rec, snap, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lp, dr, ro = snap["lump"], snap["drive"], snap["readout"]
    fig = plt.figure(figsize=(12, 8.2), facecolor=BG)
    fig.text(0.012, 0.975, "TEC Rod Thermal Lab  (Python replica)  ·  One-lump TEC (S18.3)", color=TXT,
             fontsize=12, weight="bold")
    fig.text(0.012, 0.95, f"Room {lp['room']:g} °C · Initial {lp['initial']:g} °C · τ {lp['tau']:g} s · "
             f"χ heating {lp['heat']:g} · χ cooling {lp['cool']:g} °C/PWM · u(t) = {dr['offset']:g} + {dr['amplitude']:g} "
             f"sin(2πt / {dr['period']:g} s) · t = {snap['time']:.0f} s", color=DIM, fontsize=9)
    if ro is not None:
        fig.text(0.012, 0.925, f"Settled response ({ro['period']:g} s): mean = {ro['mean']:.3f} °C | mean − room = "
                 f"{ro['mean_minus_room']:.3f} °C | fundamental amplitude = {ro['amplitude']:.3f} °C | "
                 f"phase lag = {ro['lag_deg']:.1f}°", color=TXT, fontsize=9, family="monospace")
    else:
        fig.text(0.012, 0.925, "Settled response: waiting for a complete cycle after five thermal time constants.",
                 color=TXT, fontsize=9, family="monospace")
    if snap["residual"] and ro is not None:
        x = np.linspace(0, ro["period"], len(ro["samples"]))
        gs = fig.add_gridspec(3, 1, top=0.89, hspace=0.35)
        panels = [(ro["pwm"], None, "Signed PWM (counts)", "#ffc979"),
                  (ro["samples"], ro["fitted"], "TEC temperature (°C)   Dashed: fitted sinusoid", "#4fb3ff"),
                  (ro["resid"], None, f"Residual (°C) · Amplitude {ro['resid_amp']:.3f} °C · RMS {ro['resid_rms']:.3f} °C",
                   "#57d98a")]
        for k, (vals, fit, title, col) in enumerate(panels):
            ax = fig.add_subplot(gs[k])
            _style(ax)
            ax.plot(x, vals, color=col, lw=1.4)
            if fit is not None:
                ax.plot(x, fit, color="#ff9b67", ls="--", lw=1.2)
            ax.axhline(0 if k != 1 else lp["room"], color="#607084", lw=0.8)
            ax.set_title(title, color=col, fontsize=9, loc="left")
            ax.set_xlim(0, ro["period"])
        ax.set_xlabel("time in one settled drive cycle (s)", color=DIM)
    else:
        ax = fig.add_axes([0.06, 0.08, 0.86, 0.8])
        _style(ax)
        ax.plot(rec[:, 0], rec[:, 2], color="white", lw=1.8, label="TEC (numerical)")
        if np.isfinite(rec[:, 3]).any():
            ax.plot(rec[:, 0], rec[:, 3], color="#ff8c61", ls="--", lw=1.4, label="exact (analytical)")
        ax.axhline(lp["room"], color="#a7c4d9", ls=":", lw=1.1)
        ax2 = ax.twinx()
        ax2.plot(rec[:, 0], rec[:, 1] / 2.55, color="#ffd166", alpha=0.55, lw=1)
        ax2.set_ylim(-100, 100)
        ax2.tick_params(colors="#7a6a3a", labelsize=8)
        ax2.set_ylabel("PWM % (right)", color="#b7a36b", fontsize=8)
        err = np.nanmax(np.abs(rec[:, 2] - rec[:, 3]))
        ax.set_title(f"Startup and first cycles: max |numerical − analytical| = {err:.1e} °C", color=TXT, fontsize=9, loc="left")
        ax.legend(fontsize=8, frameon=False, labelcolor=TXT)
        ax.set_xlabel("time (s)", color=DIM)
        ax.set_ylabel("°C", color=DIM)
    fig.savefig(path, dpi=150, facecolor=BG, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------
def write_records(label, rec, kind):
    OUT.mkdir(parents=True, exist_ok=True)
    if kind == "coupled":
        header = "time_s,setpoint_C,u_percent,pwm,Tplate_C,T0_C,T1_C,T2_C,T3_C,T4_C,T5_C"
        data = np.column_stack([rec[:, :3], rec[:, 2] * PWM_PER_PERCENT, rec[:, 3:]])
        np.savetxt(OUT / f"course_sim_{label}.csv", data, delimiter=",", header=header, comments="", fmt="%.10g")
    else:
        header = "time_s,pwm,TEC_C,analytical_C"
        np.savetxt(OUT / f"course_lump_{label}.csv", rec, delimiter=",", header=header, comments="", fmt="%.10g")


def compare_with_browser(label, rec, kind):
    """Largest difference from the browser record of the same case (the browser CSVs are
    rounded to 4 (coupled) or 6 (one-lump) decimals)."""
    if kind == "coupled":
        path = BROWSER / f"course_sim_{label}.csv"
        if not path.exists():
            return None
        b = np.genfromtxt(path, delimiter=",", names=True)
        n = min(len(b), len(rec))
        cols = ["setpoint_C", "Tplate_C", "T0_C", "T1_C", "T2_C", "T3_C", "T4_C", "T5_C"]
        idx = [1, 3, 4, 5, 6, 7, 8, 9]
        return max(float(np.max(np.abs(rec[:n, i] - b[c][:n]))) for c, i in zip(cols, idx))
    path = BROWSER / f"course_lump_{label}.csv"
    if not path.exists():
        return None
    b = np.genfromtxt(path, delimiter=",", names=True)
    n = min(len(b), len(rec))
    return float(np.max(np.abs(rec[:n, 2] - b["TEC_C"][:n])))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Python replica of the course six-sensor simulator: run the Module 8 cases")
    ap.add_argument("--quick", action="store_true", help="coupled runs 1 h only (test)")
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args(argv)
    FIG.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    report = []
    print("one-lump cases (S18.3, S18.4)")
    for case in LUMPED:
        label, rec, snap = run_lumped(case)
        write_records(label, rec, "lumped")
        (OUT / f"{label}_page.csv").write_text(snap["csv"] + "\n")
        figure_lumped(label, rec, snap, FIG / f"{label}.png")
        ro = snap["readout"]
        diff = compare_with_browser(label, rec, "lumped")
        line = (f"  {label}: " + ("waiting (startup)" if ro is None else
                f"mean - room {ro['mean_minus_room']:.3f} C, amplitude {ro['amplitude']:.3f} C, "
                f"lag {ro['lag_deg']:.1f} deg, residual {ro['resid_amp']:.3f} C")
                + ("" if diff is None else f"; max |replica - browser| {diff:.1e} C"))
        print(line)
        report.append(line)
    cases = [(l, P, kp, ki, 3600 if args.quick else tot) for l, P, kp, ki, tot in COUPLED]
    print(f"coupled cases ({len(cases)}, {args.jobs} in parallel; about 7 s per simulated hour each)")
    with Pool(args.jobs) as pool:
        for label, rec, snap in pool.imap(run_coupled, cases):
            write_records(label, rec, "coupled")
            (OUT / f"{label}_page.csv").write_text(snap["csv"] + "\n")
            figure_coupled(label, rec, snap, FIG / f"{label}.png")
            diff = compare_with_browser(label, rec, "coupled")
            line = (f"  {label}: {rec[-1, 0]:.0f} s simulated"
                    + ("" if diff is None else f"; max |replica - browser| {diff:.1e} C over {len(rec)} records"))
            print(line)
            report.append(line)
    (OUT / "replica_check.txt").write_text(
        "Python replica (python/analysis/course_simulator_replica.py) against the browser runs of the course page\n"
        "(data/module_08/course_simulator/; browser CSVs rounded to 1e-4 C coupled, 1e-6 C one-lump)\n"
        + "\n".join(report) + "\n")
    print(f"records: {OUT.relative_to(REPO)}; figures: {FIG.relative_to(REPO)}")
    if not args.quick:
        # the same comparison as for the browser runs: tables and figures with the tag "replica"
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import module_08_course_sim as cmp
        print("\ncomparison with the measurements and our model (module_08_course_sim.py, replica records)")
        cmp.main(["--sim-dir", str(OUT), "--tag", "replica"])


if __name__ == "__main__":
    sys.exit(main())
