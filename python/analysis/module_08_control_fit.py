# Phys 39 Module 8, S18.3-S18.8: fit the TEC boundary response of one or more intervals
#
# For each interval (one control mode and one forcing period, run_state RUN):
#   - optional settling test on the TEC thermistor (one-cycle windows, see
#     angstrom_lib.settling_test); the fit then uses the whole cycles after it
#   - time alignment: a record's temperatures are averages over the ~1 s before
#     time_ms (centred at time_ms - dt/2), the setpoint is the value at time_ms, and
#     the PWM computed at time_ms is held until the next record (centred at
#     time_ms + dt/2), dt = median record interval. The fits use these times, which
#     removes a spurious lag of dt/2 (setpoint) or dt (PWM): 1.8 or 3.6 deg at 100 s.
#   - common time fit T(t) = a + B cos(wt + phi) of TEC, T1, and the reference
#       PWM_SINE (S18.3, S18.4): reference = signed PWM, giving u_a and u_b
#       FF, P, PI (S18.5-S18.7):  reference = requested setpoint
#   - phase lag of TEC (and T1) behind the reference, residual amplitude (half the
#     residual peak-to-peak) and RMS residual
#   - FF/P/PI: mean error a - Ta, amplitude ratio B/Tb, RMS tracking error
#     sqrt(mean((T_TEC - T_set)^2)), PWM saturation in the interval and over the
#     whole step including startup
#   - PWM_SINE: effective one-lump parameters of the TEC with the rod mounted,
#     tau = tan(lag)/omega and chi = B / (u_b cos(lag)) (equal-susceptibility model;
#     with chi_h != chi_c this chi is the average of the two)
# Prints a table, writes data/module_08/analysis/control_<name>.csv, and saves one
# figure per interval plus a comparison of the settled cycles to docs/figures/module_08/.
#
# Usage (from the repository root):
#   python python/analysis/module_08_control_fit.py data/module_08/<run>.csv --step 4 --settle
#   python python/analysis/module_08_control_fit.py data/module_08/<run>.csv --t 1200,1800
#   python python/analysis/module_08_control_fit.py --intervals data/module_08/intervals_control.csv
# Options: --room T (room temperature, C; default: T_room from the run's settings),
#          --label NAME, --name NAME (output file names), --umax (default from settings).

import argparse
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

SETPOINT_MODES = ("FF", "P", "PI")
CYCLE_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
                "#6250d6", "#e34948"]


def fnum(x, fmt=".2f"):
    return "" if x is None or (isinstance(x, float) and math.isnan(x)) else format(x, fmt)


def analyse(seg: al.Segment, room=None, umax=None):
    iv, d = seg.interval, seg.data
    settings = seg.settings()
    if room is None and settings.get("T_room_C"):
        room, room_src = float(settings["T_room_C"]), "T_room setting"
    else:
        room_src = "--room" if room is not None else "none"
    if umax is None:
        umax = float(settings.get("u_max", 255))

    settling = None
    t = d["time_s"]
    start = None
    if iv.settle:
        settling = al.settling_test(d, ["TEC_C"], seg.period)
        start = settling.settled_start
    sel, n_cycles = al.complete_cycles(t, seg.period, start)
    if n_cycles < 1:
        raise ValueError(f"{iv.label}: less than one complete cycle after settling")
    dd = {k: v[sel] for k, v in d.items()}
    tt = dd["time_s"]
    w = seg.omega
    dt = float(np.median(np.diff(tt)))
    t_temp, t_pwm = tt - dt / 2, tt + dt / 2      # see the header: time alignment

    tec = al.fit_sine(t_temp, dd["TEC_C"], w)
    t1 = al.fit_sine(t_temp, dd["T1_C"], w)
    r = {"label": iv.label, "file": seg.path.name, "step": iv.step, "mode": seg.mode,
         "period_s": seg.period, "t_start_s": float(tt[0]), "t_end_s": float(tt[-1]),
         "n_cycles": n_cycles, "settled": ("auto" if iv.settle else "selected")
         if (not iv.settle or settling.settled_start is not None) else "NOT SETTLED",
         "room_C": room, "room_source": room_src,
         "TEC_mean_C": tec.mean, "TEC_mean_sigma": tec.sigma_mean,
         "TEC_mean_minus_room_C": tec.mean - room if room is not None else float("nan"),
         "TEC_amp_C": tec.amp, "TEC_amp_sigma": tec.sigma_amp,
         "TEC_resid_amp_C": tec.half_pp_resid, "TEC_resid_rms_C": tec.rms_resid,
         "T1_mean_C": t1.mean, "T1_amp_C": t1.amp, "T1_amp_sigma": t1.sigma_amp,
         "T1_resid_amp_C": t1.half_pp_resid,
         "Kp": settings.get("Kp", ""), "Ki": settings.get("Ki", ""), "u_max": umax,
         "sample_interval_s": dt, "notes": iv.notes}
    glitches = al.find_glitches(dd)
    r["glitches"] = len(glitches)
    if glitches:
        r["notes"] = ((r["notes"] + "; " if r["notes"] else "")
                      + f"{len(glitches)} glitch(es) in the fit interval, first at "
                      f"{glitches[0][0]:.0f} s: {glitches[0][1]}")

    if seg.mode == "PWM_SINE":
        ref = al.fit_sine(t_pwm, dd["signed_pwm"].astype(float), w)
        lag_tec, s_lag = al.lag(ref, tec)
        lag_t1, s_lag1 = al.lag(ref, t1)
        r.update({"reference": "signed PWM", "u_a": ref.mean, "u_b": ref.amp,
                  "TEC_lag_deg": math.degrees(lag_tec), "TEC_lag_sigma_deg": math.degrees(s_lag),
                  "T1_lag_deg": math.degrees(lag_t1), "T1_lag_sigma_deg": math.degrees(s_lag1)})
        if 0 < lag_tec < math.pi / 2:
            tau = math.tan(lag_tec) / w
            s_tau = s_lag / (w * math.cos(lag_tec) ** 2)
            chi = tec.amp / (ref.amp * math.cos(lag_tec))
            s_chi = chi * math.sqrt((tec.sigma_amp / tec.amp) ** 2 + (ref.sigma_amp / ref.amp) ** 2
                                    + (math.tan(lag_tec) * s_lag) ** 2)
            r.update({"tau_eff_s": tau, "tau_eff_sigma": s_tau,
                      "chi_eff_C_per_PWM": chi, "chi_eff_sigma": s_chi})
        else:
            r.update({"tau_eff_s": float("nan"), "chi_eff_C_per_PWM": float("nan")})
    elif seg.mode in SETPOINT_MODES:
        ref = al.fit_sine(tt, dd["setpoint_C"], w)
        lag_tec, s_lag = al.lag(ref, tec)
        lag_t1, s_lag1 = al.lag(ref, t1)
        ta = float(np.nanmedian(dd["Ta_C"]))
        tb = float(np.nanmedian(dd["Tb_C"]))
        err = dd["TEC_C"] - dd["setpoint_C"]
        sat_iv = np.abs(dd["signed_pwm"]) >= umax - 0.5
        sd = seg.step_data
        sat_step = np.abs(sd["signed_pwm"]) >= umax - 0.5
        first_sat = (float(sd["time_s"][sat_step][0] - sd["time_s"][0]) if sat_step.any()
                     else float("nan"))
        r.update({"reference": "setpoint", "Ta_C": ta, "Tb_C": tb,
                  "mean_error_C": tec.mean - ta, "amp_ratio": tec.amp / tb if tb else float("nan"),
                  "amp_ratio_sigma": tec.sigma_amp / tb if tb else float("nan"),
                  "TEC_lag_deg": math.degrees(lag_tec), "TEC_lag_sigma_deg": math.degrees(s_lag),
                  "T1_lag_deg": math.degrees(lag_t1), "T1_lag_sigma_deg": math.degrees(s_lag1),
                  "rms_tracking_C": float(np.sqrt(np.nanmean(err**2))),
                  "saturated_fraction": float(sat_iv.mean()),
                  "saturated_fraction_step": float(sat_step.mean()),
                  "first_saturation_after_step_start_s": first_sat})
    else:
        raise ValueError(f"{iv.label}: mode {seg.mode} has no forcing waveform to fit")
    return r, {"seg": seg, "sel_data": dd, "tec": tec, "t1": t1, "ref": ref, "dt": dt,
               "settling": settling}


def plot_interval(r, fits, fig_dir: Path):
    seg, dd = fits["seg"], fits["sel_data"]
    sd = seg.step_data if seg.interval.step is not None else seg.data
    w = seg.omega
    fig, axes = plt.subplots(3, 1, figsize=(12, 7.5), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1.6, 1.2]})
    ax = axes[0]
    ax.axvspan(dd["time_s"][0], dd["time_s"][-1], color=al.GRID, alpha=0.6, lw=0,
               label="fit interval")
    ax.plot(sd["time_s"], sd["TEC_C"], color=al.SENSOR_COLORS["TEC_C"], lw=1.6, label="TEC")
    ax.plot(sd["time_s"], sd["T1_C"], color=al.SENSOR_COLORS["T1_C"], lw=1.4, label="T1 (x = 0)")
    if seg.mode in SETPOINT_MODES:
        ax.plot(sd["time_s"], sd["setpoint_C"], color=al.REFERENCE_COLOR, lw=1.4, ls="--",
                label="requested setpoint")
    tt = dd["time_s"]
    ax.plot(tt, fits["tec"].model(tt - fits["dt"] / 2, w), color=al.MUTED, lw=1.0, ls=":",
            label="TEC fit")
    ax.set_ylabel("Temperature (C)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8, frameon=False)
    title = (f"{r['label']}: {seg.mode}, period {seg.period:g} s, {r['n_cycles']} cycles fitted"
             f" ({r['settled']})")
    ax.set_title(title, fontsize=10, loc="left")

    ax = axes[1]
    ax.plot(sd["time_s"], sd["signed_pwm"], color=al.INK, lw=1.4, label="applied PWM")
    if np.isfinite(sd["feedforward_pwm"]).any():
        ax.plot(sd["time_s"], sd["feedforward_pwm"], color=al.SENSOR_COLORS["T1_C"], lw=1.2,
                ls="--", label="feedforward (before x r)")
    if np.isfinite(sd["controller_pwm"]).any():
        ax.plot(sd["time_s"], sd["controller_pwm"], color=al.SENSOR_COLORS["T2_C"], lw=1.2,
                label="controller Kp e + Ki q")
    ax.axhline(0, color=al.MUTED, lw=0.8)
    ax.set_ylabel("Signed PWM")
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8, frameon=False)

    ax = axes[2]
    ax.plot(fits["tec"].t, fits["tec"].resid, color=al.INK, lw=1.0, label="TEC")
    ax.plot(fits["t1"].t, fits["t1"].resid, color=al.SENSOR_COLORS["T1_C"], lw=1.0, label="T1")
    ax.axhline(0, color=al.MUTED, lw=0.8)
    ax.set_ylabel("Fit residual (C)")
    ax.set_xlabel("Arduino time (s)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8, frameon=False)
    for a in axes:
        al.style_axes(a)
    fig.tight_layout()
    out = fig_dir / f"control_{r['label']}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_comparison(results, fitsets, fig_dir: Path, name: str):
    """Last three fitted cycles of every interval, aligned at a zero phase of its
    reference and shown as TEC minus its fitted mean."""
    fig, ax = plt.subplots(figsize=(9, 4.8))
    shown_tb = set()
    for i, (r, f) in enumerate(zip(results, fitsets)):
        seg, ref, tec = f["seg"], f["ref"], f["tec"]
        w, P = seg.omega, seg.period
        t = f["sel_data"]["time_s"] - f["dt"] / 2      # temperature times
        # latest time t0 <= t_end - 3P with w t0 + phi_ref = 0 (mod 2 pi)
        k = math.floor(((t[-1] - 3 * P) * w + ref.phase) / (2 * math.pi))
        t0 = (2 * math.pi * k - ref.phase) / w
        m = (t >= t0) & (t <= t0 + 3 * P)
        color = CYCLE_COLORS[i % len(CYCLE_COLORS)]
        ax.plot((t[m] - t0) / P, f["sel_data"]["TEC_C"][m] - tec.mean, color=color, lw=1.6,
                marker=al.SENSOR_MARKERS[al.ROD_COLUMNS[i % 5]], markevery=max(1, int(P / 8)),
                ms=5, label=f"{r['label']} ({r['mode']}, {r['period_s']:g} s)")
        if r.get("reference") == "setpoint" and round(r["Tb_C"], 2) not in shown_tb:
            x = np.linspace(0, 3, 400)
            ax.plot(x, r["Tb_C"] * np.cos(2 * np.pi * x), color=al.REFERENCE_COLOR, lw=1.2,
                    ls="--", label=f"requested, Tb = {r['Tb_C']:g} C")
            shown_tb.add(round(r["Tb_C"], 2))
    ax.axhline(0, color=al.MUTED, lw=0.8)
    ax.set_xlabel("Cycles after a zero phase of the reference (setpoint or PWM)")
    ax.set_ylabel("TEC minus its fitted mean (C)")
    ax.set_title("Settled TEC response: last three fitted cycles", fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    al.style_axes(ax)
    fig.tight_layout()
    out = fig_dir / f"control_{name}_comparison.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def print_table(results):
    print("\n| Interval | Mode | Period (s) | Mean TEC (C) | Mean - room (C) | Amp (C) | "
          "Lag vs ref (deg) | Resid amp (C) | Amp ratio | Mean err (C) | RMS track (C) | "
          "Saturation | T1 amp (C) |")
    print("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in results:
        print(f"| {r['label']} | {r['mode']} | {r['period_s']:g} | "
              f"{r['TEC_mean_C']:.2f} ± {r['TEC_mean_sigma']:.2f} | {fnum(r['TEC_mean_minus_room_C'])} | "
              f"{r['TEC_amp_C']:.2f} ± {r['TEC_amp_sigma']:.2f} | "
              f"{r['TEC_lag_deg']:.1f} ± {r['TEC_lag_sigma_deg']:.1f} | {r['TEC_resid_amp_C']:.2f} | "
              f"{fnum(r.get('amp_ratio'), '.3f')} | {fnum(r.get('mean_error_C'))} | "
              f"{fnum(r.get('rms_tracking_C'))} | "
              f"{fnum(r.get('saturated_fraction'), '.0%')} | {r['T1_amp_C']:.2f} |")
    pwm = [r for r in results if r["mode"] == "PWM_SINE"]
    if pwm:
        print("\nPWM_SINE: effective one-lump TEC parameters with the rod mounted "
              "(tau = tan(lag)/omega, chi = B/(u_b cos lag))")
        print("| Interval | Period (s) | u_a | u_b | tau_eff (s) | chi_eff (C/PWM) |")
        print("| --- | ---: | ---: | ---: | ---: | ---: |")
        for r in pwm:
            print(f"| {r['label']} | {r['period_s']:g} | {r['u_a']:.1f} | {r['u_b']:.1f} | "
                  f"{fnum(r.get('tau_eff_s'), '.1f')} ± {fnum(r.get('tau_eff_sigma'), '.1f')} | "
                  f"{fnum(r.get('chi_eff_C_per_PWM'), '.3f')} ± {fnum(r.get('chi_eff_sigma'), '.3f')} |")


def write_csv(results, path: Path):
    keys = []
    for r in results:
        keys += [k for k in r if k not in keys]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in results:
            w.writerow(r)


def intervals_from_args(args):
    if args.intervals:
        path = al.resolve_path(args.intervals)
        return al.read_intervals(path), path.stem
    if not args.run:
        raise SystemExit("give a run CSV or --intervals FILE")
    t0 = t1 = None
    if args.t:
        t0, t1 = (float(v) for v in args.t.split(","))
    label = args.label or (Path(args.run).stem + (f"_step{args.step}" if args.step is not None else ""))
    return [al.Interval(label=label, file=args.run, part=args.part, step=args.step,
                        t_start=t0, t_end=t1, settle=args.settle)], label


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run", nargs="?", help="raw run CSV")
    ap.add_argument("--intervals", help="intervals CSV (label,file,part,step,t_start_s,t_end_s,settle,exclude,notes)")
    ap.add_argument("--step", type=int)
    ap.add_argument("--t", help="START,END in Arduino seconds")
    ap.add_argument("--part", type=int, help="part of the file after an Arduino restart")
    ap.add_argument("--settle", action="store_true", help="find the settled start automatically")
    ap.add_argument("--label")
    ap.add_argument("--name", help="name for the summary CSV and comparison figure")
    ap.add_argument("--room", type=float, help="room temperature (C)")
    ap.add_argument("--umax", type=float)
    ap.add_argument("--fig-dir", default=str(al.MODULE8_FIG_DIR))
    ap.add_argument("--out-dir", default=str(al.ANALYSIS_DIR))
    args = ap.parse_args(argv)

    intervals, default_name = intervals_from_args(args)
    name = args.name or default_name
    fig_dir, out_dir = Path(args.fig_dir), Path(args.out_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    results, fitsets = [], []
    for iv in intervals:
        seg = al.select_segment(iv)
        r, f = analyse(seg, args.room, args.umax)
        if f["settling"] is not None:
            print(f"[{iv.label}] settling: {f['settling'].criterion}; settled from "
                  f"{fnum(f['settling'].settled_start, '.1f') or 'NEVER'} s")
        if r["glitches"]:
            print(f"[{iv.label}] WARNING: {r['notes']}")
        if r["n_cycles"] < 3:
            print(f"[{iv.label}] WARNING: only {r['n_cycles']} complete cycles (S18.8 asks for >= 3)")
        out = plot_interval(r, f, fig_dir)
        print(f"[{iv.label}] figure {out}")
        results.append(r)
        fitsets.append(f)
    print_table(results)
    csv_path = out_dir / f"control_{name}.csv"
    write_csv(results, csv_path)
    print(f"\nresults: {csv_path}")
    if len(results) > 1:
        print(f"comparison figure: {plot_comparison(results, fitsets, fig_dir, name)}")
    return results


if __name__ == "__main__":
    main()
