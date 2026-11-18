# Phys 39 Module 8: self-test of the analysis scripts with synthetic data
#
# Writes raw run files in the GUI's exact format (18-column CSV + _meta.txt) from
# models with known parameters, runs module_08_control_fit.py and
# module_08_angstrom_fit.py on them, and checks the recovered values.
#
#   Control run: one-lump TEC, tau = 80 s, chi = 0.46 C/PWM in both directions, driven
#     by zero-centered sinusoidal PWM (u_b = 50) at 300, 200, 100 s; then a PI-like
#     segment with a known amplitude ratio, lag, and PWM saturation.
#   Rod run: finite insulated-tip rod (L = 26.90 cm) with kappa = 6.90e-5 m^2/s and
#     nu = 1.30e-3 1/s, sensors at 0, 2.4, 5.0, 7.6, 10.2 cm, a TEC-off baseline with
#     per-channel offsets, then 100, 200, 500 s steps with startup transients, 0.01 C
#     noise, and 0.01 C rounding.
# Checks: tau and chi from the PWM_SINE fits; amplitude ratio, lag, mean error and
# saturation of the setpoint segment; the settled start found after the transient;
# kappa and nu at 100 and 200 s; at 500 s the semi-infinite nu biased by the finite tip,
# the finite-rod fit recovering it, and --model auto choosing it; the combined kappa
# and nu; the mean-profile m; output files. The numerical TEC-rod model
# (module_08_simulation.py) against the analytical one-lump and finite-rod solutions.
# Everything is written to a temporary folder.
#
# Usage (from the repository root):
#   python python/analysis/module_08_analysis_selftest.py

import contextlib
import io
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import angstrom_lib as al  # noqa: E402
import module_08_angstrom_fit as afit  # noqa: E402
import module_08_control_fit as cfit  # noqa: E402
from module_08_common import DATA_COLUMNS  # noqa: E402

T_ROOM = 22.0
KAPPA, NU = 6.90e-5, 1.30e-3
L = 0.2690
X_CM = [0.0, 2.4, 5.0, 7.6, 10.2]
OFFSETS = [0.05, 0.0, 0.3, -0.2, 0.4, 0.1]       # per-channel reading offsets (TEC, T1-T5)
TAU_TH, CHI = 80.0, 0.46
RNG = np.random.default_rng(7)


def f(x, d):
    return "nan" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{d}f}"


class RunWriter:
    def __init__(self, path: Path):
        self.path = path
        self.rows = []
        self.events = []
        self.step = 0

    def event(self, t_s, text):
        self.step += 1
        self.events.append((t_s, f"step {self.step}: {text}"))

    def record(self, t_s, state, mode, ta, tb, period, sp, ff, temps, ctrl, ratio, pwm):
        self.rows.append([str(int(round(t_s * 1000))), state, mode, str(self.step), f(ta, 2),
                          f(tb, 2), f(period, 1), f(sp, 2), f(ff, 1),
                          *[f(T, 2) for T in temps], f(ctrl, 1), f(ratio, 3), str(int(pwm))])

    def close(self):
        with open(self.path, "w") as out:
            out.write(",".join(DATA_COLUMNS) + "\n")
            for r in self.rows:
                out.write(",".join(r) + "\n")
        with open(self.path.with_name(self.path.stem + "_meta.txt"), "w") as m:
            m.write(f"Raw data file: {self.path.name}\nChannel map (x from the first hole):\n")
            m.write("  TEC_C: TEC, input A0, TEC base\n")
            for i, x in enumerate(X_CM):
                m.write(f"  T{i + 1}_C: T{i + 1}, input A{i + 1}, x = {x:.1f} cm\n")
            m.write("Events (Arduino time, wall-clock time, then the Arduino step number):\n")
            for t, text in self.events:
                m.write(f"  [{t:9.2f} s | 10:00:00] {text}\n")


def noisy(values):
    return [v + o + RNG.normal(0, 0.01) for v, o in zip(values, OFFSETS)]


def make_control_run(path):
    w = RunWriter(path)
    t, T = 1.0, T_ROOM
    expected = []
    for P in (300, 200, 100):
        w.event(t, f"START mode=PWM_SINE Ta_C=40.00 Tb_C=16.00 period_s={P:.1f} Kp=6.000 "
                   f"Ki=0.0400 u_max=255 r=2.820 chi_H=0.4850 T_room_C={T_ROOM:.2f} u_a=0.0 u_b=50.0")
        t0 = t
        u_hold = 0
        for _ in range(int(8 * P)):
            # as on the Arduino: the PWM set at the previous record is held for ~1 s,
            # and the record carries the temperature averaged over that second
            T_sum = 0.0
            for _ in range(10):                       # Euler, dt = 0.1 s
                T += 0.1 / TAU_TH * (T_ROOM + CHI * u_hold - T)
                t += 0.1
                T_sum += T
            T_avg = T_sum / 10
            u_hold = round(50 * math.sin(2 * math.pi * (t - t0) / P))
            w.record(t, "RUN", "PWM_SINE", None, None, P, None, None,
                     noisy([T_avg, T_avg - 1, T_avg - 2, T_avg - 3, T_avg - 4, T_avg - 5]),
                     None, None, u_hold)
        expected.append((P, w.step))
    # setpoint segment: TEC = Ta + 0.8 Tb sin(wt - 0.3), PWM saturates at 255
    w.event(t, "START mode=PI Ta_C=40.00 Tb_C=16.00 period_s=100.0 Kp=20.000 Ki=0.2000 "
               "u_max=255 r=2.820 chi_H=0.4850 T_room_C=22.00 u_a=0.0 u_b=50.0")
    t0, P = t, 100.0
    sat = 0
    n = 600
    for _ in range(n):
        t += 1.0
        ph = 2 * math.pi * (t - t0) / P
        sp = 40 + 16 * math.sin(ph)                   # setpoint at the record time
        ph_avg = 2 * math.pi * (t - 0.5 - t0) / P     # temperature: mean of the last second
        tec = 40 + 0.8 * 16 * math.sin(ph_avg - 0.3)
        u = min(255, 120 + 200 * math.sin(ph + 0.4))
        sat += u >= 255
        w.record(t, "RUN", "PI", 40, 16, P, sp, None, noisy([tec] * 6), u, 2.82, round(u))
    w.close()
    return expected, w.step, sat / n


def rod_temperatures(t, t_start, P, B0, theta0, tau_osc, mean_frac):
    omega = 2 * math.pi / P
    q, qp = al.q_qprime(KAPPA, NU, omega)
    m = math.sqrt(NU / KAPPA)
    x = np.array(X_CM) / 100
    A = al.finite_amplitude(x, q, qp, L)
    grow = 1 - math.exp(-(t - t_start) / tau_osc)
    mean = T_ROOM + mean_frac * al.finite_mean(x, theta0, m, L)
    osc = B0 * grow * np.abs(A) * np.cos(omega * (t - t_start) + np.angle(A))
    rod = mean + osc
    tec = rod[0] + 3 + 1.3 * B0 * grow * math.cos(omega * (t - t_start) + 0.2)
    return [tec, *rod]


def make_rod_run(path):
    w = RunWriter(path)
    t = 0.0
    for _ in range(300):                                     # TEC-off baseline
        t += 1.0
        w.record(t, "OFF", "MANUAL", None, None, None, None, None,
                 noisy([T_ROOM] * 6), None, None, 0)
    steps = {}
    theta0, tau_mean = 15.0, 400.0
    t_on = t
    for P, cycles in ((100, 35), (200, 12), (500, 8)):
        w.event(t, f"START mode=PI Ta_C=40.00 Tb_C=16.00 period_s={P:.1f} Kp=20.000 Ki=0.2000"
                if P == 100 else f"SET PERIOD {P:.4f} (was 100.0000)")
        steps[P] = (w.step, t)
        ts = t
        for _ in range(int(cycles * P)):
            t += 1.0
            frac = 1 - math.exp(-(t - t_on) / tau_mean)
            temps = rod_temperatures(t, ts, P, 8.0, theta0, 150.0, frac)
            sp = 40 + 16 * math.sin(2 * math.pi * (t - ts) / P)
            w.record(t, "RUN", "PI", 40, 16, P, sp, None, noisy(temps), 0.0, 2.82, 60)
    w.close()
    return steps


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'} {name}" + (f": {detail}" if detail else ""))
    return ok


def quiet(fn, *a):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a)


def main():
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        figs, outs = tmp / "figs", tmp / "out"

        print("Control analysis (S18.3-S18.7)")
        run = tmp / "control_run.csv"
        expected, pi_step, sat_frac = make_control_run(run)
        rows = ["label,file,part,step,t_start_s,t_end_s,settle,exclude,notes"]
        rows += [f"pwm{P},{run},,{step},,,yes,," for P, step in expected]
        rows += [f"pi100,{run},,{pi_step},,,no,,"]
        (tmp / "intervals_control.csv").write_text("\n".join(rows) + "\n")
        res = quiet(cfit.main, ["--intervals", str(tmp / "intervals_control.csv"), "--room",
                                str(T_ROOM), "--fig-dir", str(figs), "--out-dir", str(outs)])
        for r in res[:3]:
            P = r["period_s"]
            lag_true = math.degrees(math.atan(2 * math.pi / P * TAU_TH))
            results.append(check(
                f"{P:g} s: lag, tau_eff, chi_eff",
                abs(r["TEC_lag_deg"] - lag_true) < 1.0 and abs(r["tau_eff_s"] / TAU_TH - 1) < 0.03
                and abs(r["chi_eff_C_per_PWM"] / CHI - 1) < 0.03 and abs(r["u_b"] - 50) < 0.5,
                f"lag {r['TEC_lag_deg']:.1f} (true {lag_true:.1f}) deg, tau {r['tau_eff_s']:.1f} s, "
                f"chi {r['chi_eff_C_per_PWM']:.3f}"))
        r = res[3]
        results.append(check("setpoint segment: amplitude ratio, lag, mean error, saturation",
                             abs(r["amp_ratio"] - 0.8) < 0.01
                             and abs(r["TEC_lag_deg"] - math.degrees(0.3)) < 0.5
                             and abs(r["mean_error_C"] - OFFSETS[0]) < 0.02   # TEC reading offset
                             and abs(r["saturated_fraction"] - sat_frac) < 0.01
                             and r["Kp"] == "20.000",
                             f"ratio {r['amp_ratio']:.3f}, lag {r['TEC_lag_deg']:.1f} deg, "
                             f"error {r['mean_error_C']:+.3f} C, saturated {r['saturated_fraction']:.0%} "
                             f"(true {sat_frac:.0%}), Kp {r['Kp']}"))

        print("Angstrom analysis (S18.9-S18.10, Stage 6)")
        run = tmp / "rod_run.csv"
        steps = make_rod_run(run)
        rows = ["label,file,part,step,t_start_s,t_end_s,settle,exclude,notes"]
        rows += [f"p{P},{run},,{step},,,yes,," for P, (step, _) in steps.items()]
        (tmp / "intervals_rod.csv").write_text("\n".join(rows) + "\n")
        res, sensors, comb = quiet(afit.main, [
            "--intervals", str(tmp / "intervals_rod.csv"), "--baseline", f"{run}:50,290",
            "--fig-dir", str(figs), "--out-dir", str(outs), "--name", "selftest"])
        by_p = {r["period_s"]: r for r in res}

        r = by_p[100]
        t_step = steps[100][1]
        results.append(check("100 s: settled start found after the transient, >= 5 cycles",
                             r["settled"] == "auto" and r["t_start_s"] > t_step + 1000
                             and r["n_cycles"] >= 5,
                             f"fit from {r['t_start_s'] - t_step:.0f} s after the step, "
                             f"{r['n_cycles']} cycles"))
        for P in (100, 200):
            r = by_p[P]
            q, qp = al.q_qprime(KAPPA, NU, 2 * math.pi / P)
            ok_k = abs(r["kappa_m2_s"] - KAPPA) < max(3 * r["kappa_sigma"], 0.02 * KAPPA)
            ok_n = abs(r["nu_finite_per_s"] - NU) < 3 * r["nu_finite_sigma"]
            results.append(check(
                f"{P} s: q, q' and kappa recovered; finite-fit nu within 3 sigma",
                abs(r["q_per_m"] / q - 1) < 0.02 and abs(r["qp_per_m"] / qp - 1) < 0.02 and ok_k and ok_n,
                f"q {r['q_per_m']:.2f} ({q:.2f}), q' {r['qp_per_m']:.2f} ({qp:.2f}), kappa "
                f"{r['kappa_m2_s']:.3e} ± {r['kappa_sigma']:.1e}, nu semi {r['nu_per_s']:.2e} / finite "
                f"{r['nu_finite_per_s']:.2e} ± {r['nu_finite_sigma']:.1e}"))
        r = by_p[500]
        t5 = [s for s in sensors if s["label"] == "p500" and s["sensor"] == "T5"][0]
        results.append(check(
            "500 s: semi-infinite nu biased by the tip, flagged; finite fit recovers nu",
            r["nu_per_s"] > 1.15 * NU and t5["finite_corr_exceeds_uncertainty"]
            and abs(r["nu_finite_per_s"] - NU) < 3 * r["nu_finite_sigma"]
            and abs(r["kappa_finite_m2_s"] - KAPPA) < 3 * r["kappa_finite_sigma"],
            f"nu semi {r['nu_per_s']:.2e}, finite {r['nu_finite_per_s']:.2e} ± "
            f"{r['nu_finite_sigma']:.1e} (true {NU:.2e}); T5 amplitude correction "
            f"{t5['finite_amp_corr']:+.1%}"))
        results.append(check(
            "auto model: finite fit used at 500 s; combined kappa and nu within 3 sigma",
            r["model_used"] == "finite"
            and abs(comb["kappa"][0] - KAPPA) < 3 * comb["kappa"][1]
            and abs(comb["nu"][0] - NU) < 3 * comb["nu"][1],
            "models " + ", ".join(f"{q['period_s']:g} s {q['model_used']}" for q in res)
            + f"; kappa {comb['kappa'][0]:.3e} ± {comb['kappa'][1]:.1e}, nu {comb['nu'][0]:.3e} ± "
              f"{comb['nu'][1]:.1e}"))
        ms = [r["m_mean_fit_per_m"] for r in res]
        results.append(check("mean profile: m from baseline-referenced means",
                             all(abs(m - math.sqrt(NU / KAPPA)) < 3 * r["m_mean_fit_sigma"] + 0.05
                                 for m, r in zip(ms, res)),
                             ", ".join(f"{m:.2f}" for m in ms) + f" (true {math.sqrt(NU / KAPPA):.2f})"))
        results.append(check("phase lag grows downstream (all sensors retained at 100 s)",
                             by_p[100]["retained"] == "T1 T2 T3 T4 T5" and by_p[100]["qp_per_m"] > 0))
        results.append(check("no false glitches in clean synthetic data",
                             all(r["glitches"] == 0 for r in res) and
                             all(r["glitches"] == 0 for r in quiet(cfit.main, [
                                 "--intervals", str(tmp / "intervals_control.csv"), "--room",
                                 str(T_ROOM), "--fig-dir", str(figs), "--out-dir", str(outs)]))))
        d = {k: v.copy() for k, v in al.select_segment(al.read_intervals(tmp / "intervals_rod.csv")[0]).data.items()}
        clean = len(al.find_glitches(d))
        for c in ("TEC_C", "T1_C", "T2_C", "T3_C", "T4_C", "T5_C"):
            d[c][100] -= 2.0                      # all channels drop 2 C for one record
        d["T4_C"][300] += 1.0                     # one downstream channel jumps 1 C
        found = al.find_glitches(d)
        times = sorted({round(t - d["time_s"][0]) for t, _ in found})
        results.append(check("glitch finder: common-mode drop and single-channel jump found",
                             clean == 0 and any("common-mode" in g for _, g in found)
                             and any("T4" in g for _, g in found)
                             and set(times) <= {100, 101, 102, 300, 301, 302},
                             f"{len(found)} flags at {times} s after the start"))
        expected_figs = ["control_pwm300.png", "control_intervals_control_comparison.png",
                         "angstrom_p100_fits.png", "angstrom_p100_profile_theory.png",
                         "angstrom_p100_settling.png", "angstrom_selftest_lengths.png",
                         "angstrom_selftest_frequency.png"]
        missing = [n for n in expected_figs if not (figs / n).exists()]
        results.append(check("figures and result tables written",
                             not missing and (outs / "angstrom_selftest.csv").exists()
                             and (outs / "control_intervals_control.csv").exists(),
                             f"missing {missing}" if missing else ""))
        if "--keep" in sys.argv:
            import shutil
            keep = Path(sys.argv[sys.argv.index("--keep") + 1])
            shutil.copytree(figs, keep, dirs_exist_ok=True)
            print(f"figures copied to {keep}")

    print("Numerical TEC-rod model (module_08_simulation.py)")
    results += check_simulation()

    print("Python replica of the course simulator (course_simulator_replica.py)")
    results += check_replica()

    print("PASS" if all(results) else "FAIL")
    return 0 if all(results) else 1


def check_replica():
    import course_simulator_replica as R
    out = []
    s = R.lumped(300, 0, 50, tau=80, heat=0.46, cool=0.46)
    worst = 0.0
    for _ in range(3000):
        s.advance(0.1)
        worst = max(worst, abs(s.Tc - s.analytical_at(s.time)))
    for _ in range(3000):
        s.advance(1.0)
    ro = s.settled_readout()
    w = 2 * math.pi / 300
    amp = 0.46 * 50 / math.sqrt(1 + (w * 80) ** 2)
    out.append(check("replica one-lump: RK4 on the analytical curve, settled readout = theory",
                     worst < 1e-6 and abs(ro["amplitude"] - amp) < 1e-3
                     and abs(ro["lag_deg"] - math.degrees(math.atan(w * 80))) < 0.05,
                     f"max |num - ana| {worst:.1e} C, amplitude {ro['amplitude']:.4f} ({amp:.4f}) C, "
                     f"lag {ro['lag_deg']:.2f} deg"))
    s = R.coupled(200, 30 / 2.55, 0.3 / 2.55, L_mm=285)
    for _ in range(600):
        s.advance(1.0)
    out.append(check("replica coupled: 100 nodes for 285 mm, PID drives T0 toward the setpoint",
                     s.N == 100 and abs(s.Ts[0] - s.pid.sp) < 5 and s.T[0] > 30,
                     f"N {s.N}, T0 {s.Ts[0]:.2f} C, setpoint {s.pid.sp:.2f} C after 600 s"))
    return out


def check_simulation():
    import module_08_simulation as sim
    out = []
    # one-lump, equal chi: settled amplitude and lag against the boxed S18.3 solution
    P, chi, tau, ub = 200.0, 0.46, 80.0, 50.0
    model = sim.Model(sim.Lump(chi, chi, tau), None, dt=0.1)
    ctrl = sim.Controller("PWM_SINE", P, u_b=ub, integer=False, dt=0.1)
    rec, fits, _ = sim.periodic_run(model, 22.0, ctrl, P, 12)
    w = 2 * math.pi / P
    amp = chi * ub / math.sqrt(1 + (w * tau) ** 2)
    lag = math.degrees(math.atan(w * tau))
    got_lag = sim.lag_deg(fits["signed_pwm"], fits["TEC_C"])
    out.append(check("one-lump: amplitude and lag match the analytical solution",
                     abs(fits["TEC_C"].amp - amp) < 2e-3 * amp and abs(got_lag - lag) < 0.1
                     and abs(fits["TEC_C"].mean - 22.0) < 1e-3,
                     f"amp {fits['TEC_C'].amp:.4f} ({amp:.4f}) C, lag {got_lag:.2f} ({lag:.2f}) deg"))
    # unequal chi: mean shift (chi_h - chi_c) u_b / pi
    model = sim.Model(sim.Lump(0.46, 0.23, tau), None, dt=0.1)
    ctrl = sim.Controller("PWM_SINE", P, u_b=ub, integer=False, dt=0.1)
    rec, fits, _ = sim.periodic_run(model, 22.0, ctrl, P, 12)
    shift = (0.46 - 0.23) * ub / math.pi
    out.append(check("one-lump: unequal chi raises the mean by (chi_h - chi_c) u_b / pi",
                     abs(fits["TEC_C"].mean - 22.0 - shift) < 0.01,
                     f"{fits['TEC_C'].mean - 22.0:.3f} ({shift:.3f}) C"))
    # rod: a sinusoidal TEC (no rod load, g = 0) through the contact gives the analytical
    # finite-rod transfer at T1 and the finite-rod ratio between T1 and T5
    rod = sim.Rod(KAPPA, NU, 0.04, 0.02, dx_m=0.001)
    model = sim.Model(sim.Lump(0.46, 0.46, 1.0, 0.0), rod, dt=1.0)
    ctrl = sim.Controller("PWM_SINE", P, u_b=20, integer=False, dt=1.0)
    rec, fits, _ = sim.periodic_run(model, 22.0, ctrl, P, 60, n_fit=5)
    h_num = (fits["T1_C"].amp / fits["TEC_C"].amp,
             al.wrap(fits["TEC_C"].phase - fits["T1_C"].phase))
    h = sim.contact_transfer(w, KAPPA, NU, 0.04, 0.02)
    A = al.finite_amplitude(0.102, *al.q_qprime(KAPPA, NU, w), L=sim.ROD_OVERALL_LENGTH_M - 0.02)
    r15 = fits["T5_C"].amp / fits["T1_C"].amp
    out.append(check("rod: contact transfer and T5/T1 ratio match the analytical finite rod",
                     abs(h_num[0] / abs(h) - 1) < 0.01 and abs(h_num[1] + np.angle(h)) < 0.01
                     and abs(r15 / abs(A) - 1) < 0.01,
                     f"T1/TEC {h_num[0]:.4f} ({abs(h):.4f}), lag {h_num[1]:.4f} "
                     f"({-np.angle(h):.4f}) rad, T5/T1 {r15:.4f} ({abs(A):.4f})"))
    # controller: PI removes the mean error that P leaves
    lump = sim.Lump(0.5, 0.2, 80.0)
    res = {}
    for name, ki in (("P", 0.0), ("PI", 0.05)):
        ctrl = sim.Controller("PI" if ki else "P", 300, Ta=40, Tb=5, kp=6, ki=ki, T_room_set=22)
        rec, fits, sl = sim.periodic_run(sim.Model(lump, None, 1.0), 22.0, ctrl, 300, 30)
        res[name] = fits["TEC_C"].mean - 40
    out.append(check("controller port: P leaves a droop, PI removes it",
                     res["P"] < -1 and abs(res["PI"]) < 0.05,
                     f"P {res['P']:+.2f} C, PI {res['PI']:+.3f} C"))
    return out


if __name__ == "__main__":
    sys.exit(main())
