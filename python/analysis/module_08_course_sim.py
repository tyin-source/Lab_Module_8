# Phys 39 Module 8: compare the course's six-sensor simulator with the measurements and
# with our TEC-rod model
#
# Input: data/module_08/course_simulator/course_sim_<label>.csv, written by
#   node tools/course_simulator_driver.js data/module_08/course_simulator
# (the course page run headlessly: coupled TEC + rod, aluminum 6061, d 12.7 mm, L 285 mm,
# Angstrom sinusoid 40 + 16 sin, PID on the simulator's T0 at x = 0; our gains converted
# to % duty per C by / 2.55; rod temperatures sampled at our sensor positions, TEC face +
# holes 1, 3, 5, 6, 7).
#
# 1. Focused S19 simulation at 100 s: P (Kp 6) vs PI (Kp 6, Ki 0.04) and the adjusted PI
#    (Kp 30, Ki 0.3): mean error, amplitude ratio, lag, RMS error, saturation of the
#    feedback temperature, next to the measured runs and our model (sim_control.csv).
# 2. Stage 6 integrated comparison at every production period (matched run 300 s):
#    amplitudes and lags (relative to the requested setpoint) of the feedback sensor and
#    T1-T5 for the measurement, the course simulator, and our model
#    (sim_rod_<label>_independent.csv); then the same normalized to T1, which removes the
#    boundary and isolates heat transport in the rod.
#
# 0. S18.3 and S18.4 from the simulator's one-lump mode (tools/course_simulator_lumped.js),
#    next to our Python one-lump values.
# Outputs: data/module_08/analysis/course_sim_lumped.csv, course_sim_s19.csv, course_sim_rod.csv;
#          docs/figures/module_08/course_sim_s19_p_pi100.png, course_sim_p300_comparison.png
#
# Usage (from the repository root):
#   python python/analysis/module_08_course_sim.py

import csv
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import angstrom_lib as al  # noqa: E402

SIM_DIR = al.MODULE8_DATA_DIR / "course_simulator"
TAG = ""                         # "" for the browser runs, "replica" for the Python replica


def out_name(base, ext):
    return f"{base}_{TAG}.{ext}" if TAG else f"{base}.{ext}"
ANALYSIS = al.ANALYSIS_DIR
FIGS = al.REPO_ROOT / "docs" / "figures" / "module_08"
ROD = ["T1_C", "T2_C", "T3_C", "T4_C", "T5_C"]
LABELS = {"T0_C": "TEC/T0", "T1_C": "T1", "T2_C": "T2", "T3_C": "T3", "T4_C": "T4", "T5_C": "T5"}


def load(label):
    return np.genfromtxt(SIM_DIR / f"course_sim_{label}.csv", delimiter=",", names=True)


def fit_last(d, period, cycles=5):
    t = d["time_s"]
    n = int(round(cycles * period))
    s = slice(len(t) - n, len(t))
    w = 2 * math.pi / period
    fits = {c: al.fit_sine(t[s], d[c][s], w) for c in ["setpoint_C", "T0_C", "Tplate_C"] + ROD}
    return fits, s


def lag_deg(ref, resp):
    """Lag behind ref in degrees, in [-90, 270): sensors far down the rod lag by more
    than 180 deg at 100 s."""
    return (math.degrees(al.wrap(ref.phase - resp.phase)) + 90) % 360 - 90


def lag_range(v):
    return (v + 90) % 360 - 90


def rows_of(path):
    return list(csv.DictReader(open(path)))


def write(path, rows):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    return path.relative_to(al.REPO_ROOT)


def s19():
    ours = {(r["case"], r["source"].split()[0]): r for r in rows_of(ANALYSIS / "sim_control.csv")}
    cases = [("P only", "s19_p100", 6, 0), ("PI, same Kp", "s19_pi100", 6, 0.04),
             ("PI, targeted adjustment", "s19_pi100_adj", 30, 0.3)]
    out, store = [], {}
    for case, label, kp, ki in cases:
        d = load(label)
        f, s = fit_last(d, 100)
        T0, S = f["T0_C"], f["setpoint_C"]
        err = d["T0_C"][s] - d["setpoint_C"][s]
        out.append({"case": case, "Kp_PWM_per_C": kp, "Ki_PWM_per_Cs": ki,
                    "Kp_percent_per_C": kp / 2.55, "Ki_percent_per_Cs": ki / 2.55,
                    "source": "course simulator (feedback T0)", "mean_error_C": T0.mean - 40,
                    "amp_ratio": T0.amp / 16, "lag_deg": lag_deg(S, T0),
                    "rms_tracking_C": float(np.sqrt(np.mean(err**2))),
                    "saturated_fraction": float(np.mean(np.abs(d["u_percent"][s]) >= 99.99)),
                    "pwm_min": float(d["pwm"][s].min()), "pwm_max": float(d["pwm"][s].max()),
                    "T1_amp_C": f["T1_C"].amp})
        store[case] = (d, s)
        for src, name in (("measured", "measured (feedback TEC)"), ("simulation", "our model (feedback TEC)")):
            r = ours[(case, src)]
            out.append({"case": case, "Kp_PWM_per_C": kp, "Ki_PWM_per_Cs": ki, "source": name,
                        "mean_error_C": float(r["mean_error_C"]), "amp_ratio": float(r["amp_ratio"]),
                        "lag_deg": float(r["lag_deg"]), "rms_tracking_C": float(r["rms_tracking_C"]),
                        "saturated_fraction": float(r["saturated_fraction"]),
                        "T1_amp_C": float(r["T1_amp_C"])})
    print("Focused S19 simulation, 100 s, Ta 40 C, Tb 16 C")
    print("| Case | Kp (PWM/C) | Ki (PWM/(C s)) | Source | Mean error (C) | Amplitude ratio | Lag (deg) | "
          "RMS error (C) | Saturated | T1 amplitude (C) |")
    print("| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in out:
        print(f"| {r['case']} | {r['Kp_PWM_per_C']:g} | {r['Ki_PWM_per_Cs']:g} | {r['source']} | "
              f"{r['mean_error_C']:+.2f} | {r['amp_ratio']:.3f} | {r['lag_deg']:.1f} | "
              f"{r['rms_tracking_C']:.2f} | {100 * r['saturated_fraction']:.0f} % | {r['T1_amp_C']:.2f} |")
    print(f"table: {write(ANALYSIS / out_name('course_sim_s19', 'csv'), out)}")

    fig, axes = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    cols = [al.SENSOR_COLORS["T2_C"], al.SENSOR_COLORS["T1_C"], al.SENSOR_COLORS["T3_C"]]
    first = True
    for (case, *_), col in zip(cases, cols):
        d, s = store[case]
        t = d["time_s"][s] - d["time_s"][s][0]
        if first:
            axes[0].plot(t, d["setpoint_C"][s], color=al.REFERENCE_COLOR, lw=1.2, label="setpoint 40 + 16 sin")
            first = False
        axes[0].plot(t, d["T0_C"][s], color=col, label=f"{case} (T0)")
        axes[1].plot(t, d["pwm"][s], color=col)
    axes[0].set_ylim(top=63)
    axes[0].legend(fontsize=7, frameon=False, ncol=4, loc="upper center")
    axes[0].set_ylabel("feedback temperature T0 (C)")
    axes[1].set_ylabel("PWM (2.55 x duty %)")
    axes[1].set_xlabel("time in the last five cycles (s)")
    axes[0].set_title("Course six-sensor simulator, 100 s: P (Kp 6 PWM/C) vs PI (Ki 0.04) and the "
                      "adjusted PI (Kp 30, Ki 0.3), gains / 2.55 in %/C", fontsize=9)
    for ax in axes:
        al.style_axes(ax)
    FIGS.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGS / out_name("course_sim_s19_p_pi100", "png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"figure: docs/figures/module_08/{out_name('course_sim_s19_p_pi100', 'png')}")


def rod():
    out = []
    for P in (100, 200, 300, 400, 500):
        label = "s19_pi100_adj" if P == 100 else f"pi{P}"
        f, _ = fit_last(load(label), P)
        S = f["setpoint_C"]
        ours = {r["sensor"]: r for r in rows_of(ANALYSIS / f"sim_rod_p{P}_independent.csv")}
        for c in ["T0_C"] + ROD:
            key = "TEC" if c == "T0_C" else c[:-2]
            o = ours[key]
            out.append({"period_s": P, "sensor": LABELS[c],
                        "meas_mean_C": float(o["meas_mean_C"]), "meas_amp_C": float(o["meas_amp_C"]),
                        "meas_lag_deg": lag_range(float(o["meas_lag_deg"])),
                        "course_mean_C": f[c].mean, "course_amp_C": f[c].amp, "course_lag_deg": lag_deg(S, f[c]),
                        "our_model_mean_C": float(o["sim_mean_C"]), "our_model_amp_C": float(o["sim_amp_C"]),
                        "our_model_lag_deg": lag_range(float(o["sim_lag_deg"]))})
    # normalized to T1
    for P in (100, 200, 300, 400, 500):
        rows = [r for r in out if r["period_s"] == P]
        t1 = next(r for r in rows if r["sensor"] == "T1")
        for r in rows:
            for src in ("meas", "course", "our_model"):
                r[f"{src}_amp_over_T1"] = r[f"{src}_amp_C"] / t1[f"{src}_amp_C"]
                r[f"{src}_lag_from_T1_deg"] = r[f"{src}_lag_deg"] - t1[f"{src}_lag_deg"]
    print("\nStage 6: measurement / course simulator / our model (lags behind the requested setpoint)")
    print("| Period (s) | Sensor | Amplitude (C) | Lag (deg) | B/B1 | Lag from T1 (deg) |")
    print("| ---: | --- | ---: | ---: | ---: | ---: |")
    for r in out:
        print(f"| {r['period_s']} | {r['sensor']} | {r['meas_amp_C']:.2f} / {r['course_amp_C']:.2f} / "
              f"{r['our_model_amp_C']:.2f} | {r['meas_lag_deg']:.1f} / {r['course_lag_deg']:.1f} / "
              f"{r['our_model_lag_deg']:.1f} | {r['meas_amp_over_T1']:.3f} / {r['course_amp_over_T1']:.3f} / "
              f"{r['our_model_amp_over_T1']:.3f} | {r['meas_lag_from_T1_deg']:.1f} / "
              f"{r['course_lag_from_T1_deg']:.1f} / {r['our_model_lag_from_T1_deg']:.1f} |")
    print(f"table: {write(ANALYSIS / out_name('course_sim_rod', 'csv'), out)}")
    return out


def figure_p300(out):
    # measured p300 record, aligned on the setpoint phase like the simulations
    res = {r["label"]: r for r in rows_of(ANALYSIS / "angstrom_production.csv")}["p300"]
    iv = al.Interval("p300", "data/module_08/" + res["file"], step=int(float(res["step"])),
                     t_start=float(res["t_start_s"]), t_end=float(res["t_end_s"]))
    seg = al.select_segment(iv)
    P, w = seg.period, seg.omega
    tm = seg.data["time_s"]
    Sm = al.fit_sine(tm, seg.data["setpoint_C"], w)
    xm = (tm + Sm.phase / w)            # setpoint phase zero (sine convention offset is common)
    d = load("pi300")
    f, s = fit_last(d, P)
    ts = d["time_s"][s]
    xs = ts + f["setpoint_C"].phase / w
    shift = round((xs[-1] - xm[-1]) / P) * P     # end the simulated cycles with the measured ones
    xs = xs - shift
    fig = plt.figure(figsize=(11, 10))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.3, 1.6, 1.2], hspace=0.35)
    a0 = fig.add_subplot(gs[0])
    axes = [a0, fig.add_subplot(gs[1], sharex=a0), fig.add_subplot(gs[2])]
    ax = axes[0]
    ax.plot(xm, seg.data["setpoint_C"], color=al.REFERENCE_COLOR, lw=1, label="setpoint")
    ax.plot(xm, seg.data["TEC_C"], color=al.INK, lw=1.5, label="TEC measured (feedback)")
    ax.plot(xs, d["T0_C"][s], color=al.INK, lw=1, ls="--", label="course sim T0 (feedback)")
    ax.plot(xm, seg.data["T1_C"], color=al.SENSOR_COLORS["T1_C"], lw=1.5, label="T1 measured")
    ax.plot(xs, d["T1_C"][s], color=al.SENSOR_COLORS["T1_C"], lw=1, ls="--", label="T1 course sim")
    ax.set_ylabel("T (C)")
    ax.set_title("300 s: course six-sensor simulator vs measurement, aligned on the setpoint phase. "
                 "Boundary first: the simulator's T1 swings far more", fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    ax = axes[1]
    for c in ["T1_C", "T3_C", "T5_C"]:
        col = al.SENSOR_COLORS[c]
        ax.plot(xm, seg.data[c], color=col, lw=1.6, label=f"{LABELS[c]} measured")
        ax.plot(xs, d[c][s], color=col, lw=1, ls="--", label=f"{LABELS[c]} course sim")
    ax.set_ylabel("T (C)")
    ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    ax = axes[2]
    rows = [r for r in out if r["period_s"] == 300 and r["sensor"] != "TEC/T0"]
    x = np.arange(len(rows))
    for k, (src, col, mk) in enumerate([("meas", al.INK, "o"), ("course", al.SENSOR_COLORS["T2_C"], "s"),
                                        ("our_model", al.SENSOR_COLORS["T1_C"], "^")]):
        ax.plot(x, [r[f"{src}_amp_over_T1"] for r in rows], marker=mk, color=col, lw=1,
                label={"meas": "measured", "course": "course simulator", "our_model": "our model"}[src])
    ax.set_xticks(x, [r["sensor"] for r in rows])
    ax.set_ylabel("amplitude / T1 amplitude")
    ax.legend(fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.set_title("Normalized to T1 (rod transport only): the three agree", fontsize=9)
    axes[1].set_xlabel("time from a setpoint zero crossing (s); last three cycles")
    for a in axes:
        al.style_axes(a)
    t_end = float(xm[-1])
    axes[0].set_xlim(t_end - 3 * P, t_end)
    fig.savefig(FIGS / out_name("course_sim_p300_comparison", "png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"figure: docs/figures/module_08/{out_name('course_sim_p300_comparison', 'png')}")


def lumped():
    """S18.3 course table and S18.4 simulation table from the course simulator's own
    one-lump mode (tools/course_simulator_lumped.js), next to our Python one-lump values."""
    ours3 = {(int(float(r["period_s"])), round(float(r["chi_c"]), 2)): r
             for r in rows_of(ANALYSIS / "sim_s18_3_course.csv")}
    out = []
    print("\nS18.3 table from the course simulator's one-lump mode (tau 80 s, chi_h 0.46, u_b 50, room 22 C); "
          "ours in brackets")
    print("| Period (s) | chi_c | Mean - room (C) | Amplitude (C) | Phase lag (deg) | Residual (C) | max abs(num - analytical) (C) |")
    print("| ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for P in (100, 200, 300):
        for cc in (0.23, 0.46):
            d = np.genfromtxt(SIM_DIR / f"course_lump_s18_3_P{P}_chic{cc}.csv", delimiter=",", names=True)
            t = d["time_s"]
            m = t > t[-1] - 3 * P + 1e-6
            w = 2 * math.pi / P
            T = al.fit_sine(t[m], d["TEC_C"][m], w)
            U = al.fit_sine(t[m], d["pwm"][m], w)
            ana = d["analytical_C"]
            err = float(np.nanmax(np.abs(d["TEC_C"] - ana))) if np.isfinite(ana).any() else float("nan")
            o = ours3[(P, cc)]
            r = {"case": "S18.3", "period_s": P, "chi_c": cc, "mean_minus_room_C": T.mean - 22,
                 "amp_C": T.amp, "lag_deg": lag_deg(U, T), "resid_amp_C": T.half_pp_resid,
                 "max_abs_num_minus_analytical_C": err,
                 "ours_mean_minus_room_C": float(o["mean_minus_room_C"]), "ours_amp_C": float(o["amp_C"]),
                 "ours_lag_deg": float(o["lag_deg"]), "ours_resid_amp_C": float(o["resid_amp_C"])}
            out.append(r)
            print(f"| {P} | {cc:.2f} | {r['mean_minus_room_C']:.2f} ({r['ours_mean_minus_room_C']:.2f}) | "
                  f"{r['amp_C']:.2f} ({r['ours_amp_C']:.2f}) | {r['lag_deg']:.1f} ({r['ours_lag_deg']:.1f}) | "
                  f"{r['resid_amp_C']:.3f} ({r['ours_resid_amp_C']:.3f}) | "
                  + ("—" if math.isnan(err) else f"{err:.1e}") + " |")
    ours4 = {(r["drive"], r["source"]): r for r in rows_of(ANALYSIS / "sim_s18_4.csv")}
    print("\nS18.4 table from the course simulator's one-lump mode (measured chi_h 0.412, chi_c 0.173, tau 59.1 s, "
          "room 21.73 C); measured and our one-lump in brackets")
    print("| 100 s drive | u_a | u_b | Fundamental amplitude (C) | Residual amplitude (C) |")
    print("| --- | ---: | ---: | ---: | ---: |")
    for label, kind, ua in (("s18_4_zero100", "zero-centered", 0), ("s18_4_pos100", "positive-only", 55)):
        d = np.genfromtxt(SIM_DIR / f"course_lump_{label}.csv", delimiter=",", names=True)
        t = d["time_s"]
        m = t > t[-1] - 300 + 1e-6
        T = al.fit_sine(t[m], d["TEC_C"][m], 2 * math.pi / 100)
        me, ol = ours4[(kind, "measured")], ours4[(kind, "one-lump")]
        out.append({"case": "S18.4", "drive": kind, "u_a": ua, "u_b": 50, "mean_C": T.mean,
                    "amp_C": T.amp, "resid_amp_C": T.half_pp_resid,
                    "measured_amp_C": float(me["amp_C"]), "measured_resid_amp_C": float(me["resid_amp_C"]),
                    "ours_amp_C": float(ol["amp_C"]), "ours_resid_amp_C": float(ol["resid_amp_C"])})
        print(f"| {kind} | {ua} | 50 | {T.amp:.2f} (meas {float(me['amp_C']):.2f}, ours {float(ol['amp_C']):.2f}) | "
              f"{T.half_pp_resid:.3f} (meas {float(me['resid_amp_C']):.3f}, ours {float(ol['resid_amp_C']):.3f}) |")
    print(f"table: {write(ANALYSIS / out_name('course_sim_lumped', 'csv'), out)}")


def main(argv=None):
    import argparse
    global SIM_DIR, TAG
    ap = argparse.ArgumentParser(description="course simulator vs measurement vs our model")
    ap.add_argument("--sim-dir", help="records directory (default: the browser runs)")
    ap.add_argument("--tag", default="", help="suffix for the output tables and figures")
    args = ap.parse_args(argv)
    if args.sim_dir:
        SIM_DIR = Path(args.sim_dir)
    TAG = args.tag
    lumped()
    s19()
    out = rod()
    figure_p300(out)


if __name__ == "__main__":
    main()
