# Phys 39 Module 8, S18.9-S18.10 and Stages 5-6: Angstrom analysis of the rod sensors
#
# For each interval (one forcing period, run_state RUN):
#   1. settling test on T1 and the rod sensors (one-cycle windows; optional), then the
#      whole cycles after the settled start (>= 5 required for production)
#   2. common time fit T_i = a_i + B_i cos(wt + phi_i) for TEC and T1-T5
#   3. retained sensors: not excluded and measurable, B_i >= max(--min-amp, --min-snr sigma_B)
#   4. spatial fits  ln B_i = c - q x_i  and unwrapped  phi_i = c' - q' x_i  (phases
#      relative to the sensor nearest the base), weighted, with position uncertainty
#   5. kappa = omega/(2 q q'), nu = kappa (q^2 - q'^2); k = kappa rho c, h = nu rho c d/4
#      with the course rho and c (only for comparison)
#   6. checks: finite insulated-tip correction to amplitude and phase at every retained
#      sensor, compared with its measurement uncertainty; Bi_perp = nu d^2/(16 kappa);
#      q'R (radial uniformity of the oscillation)
#   7. mean profile: theta_i = a_i - T_inf,i fitted to theta0 cosh[m(L-x)]/cosh(mL);
#      m compared with sqrt(nu/kappa) and the Module 7 value (needs --room or --baseline)
#   8. the same amplitudes and phases fitted to the finite insulated-tip rod (model
#      check). nu = kappa (q^2 - q'^2) is a small difference of large numbers, so a
#      tip correction of 1-3 % in amplitude can move nu by tens of percent. --model
#      auto (default) uses the finite fit for an interval when the tip correction at
#      any retained sensor exceeds its measurement uncertainty, else the semi-infinite
#      fit; --model semi or finite forces one. The choice is reported per interval.
#   9. finite-rod prediction of T(x,t) at the retained sensors from the measured T1
#      fit and the kappa, nu used, overlaid on the data
# With several intervals: kappa and nu vs period, their weighted means and a chi^2
# consistency test, and the characteristic-length plot (1/m, 1/q, 1/q' vs period,
# 100-500 s) from the combined kappa and nu with the fitted lengths overlaid.
#
# Outputs: docs/figures/module_08/angstrom_<label>_*.png and angstrom_<name>_*.png,
# data/module_08/analysis/angstrom_<name>.csv (one row per interval) and
# angstrom_<name>_sensors.csv (one row per sensor and interval).
#
# Usage (from the repository root):
#   python python/analysis/module_08_angstrom_fit.py data/module_08/<run>.csv --step 7 --settle \
#       --baseline data/module_08/<baseline>.csv:60,600
#   python python/analysis/module_08_angstrom_fit.py --intervals data/module_08/intervals_production.csv \
#       --room 22.3
# Options: --exclude T5, --positions-cm 0,2.4,5.0,7.6,10.2, --length-cm 26.90,
#          --min-amp 0.05, --min-snr 3, --theory LABEL (interval for the theory overlay).

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
from module_08_common import ROD_COLUMNS, SENSOR_LABELS, TEMPERATURE_COLUMNS  # noqa: E402



def fnum(x, fmt=".3g"):
    return "" if x is None or (isinstance(x, float) and math.isnan(x)) else format(x, fmt)


def analyse(seg: al.Segment, args, t_inf):
    iv, d = seg.interval, seg.data
    exclude = {e.upper() for e in iv.exclude} | {e.upper() for e in (args.exclude or [])}
    candidates = [c for c in ROD_COLUMNS
                  if SENSOR_LABELS[c] not in exclude and np.isfinite(d[c]).mean() > 0.95]
    L = args.length_cm / 100
    pos_cm = al.positions_for(seg.meta, args.positions_cm)

    settling = None
    start = None
    if iv.settle:
        settling = al.settling_test(d, candidates, seg.period, k_sigma=args.k_sigma,
                                    floor_mean=args.floor_mean, floor_amp=args.floor_amp,
                                    floor_phase=args.floor_phase)
        start = settling.settled_start
    sel, n_cycles = al.complete_cycles(d["time_s"], seg.period, start)
    dd = {k: v[sel] for k, v in d.items()}
    tt = dd["time_s"]
    w = seg.omega

    fits = {}
    for c in TEMPERATURE_COLUMNS:
        if np.isfinite(dd[c]).sum() >= 8:
            fits[c] = al.fit_sine(tt, dd[c], w)
    retained = [c for c in candidates if c in fits
                and fits[c].amp >= max(args.min_amp, args.min_snr * fits[c].sigma_amp)]
    if len(retained) < 2:
        raise ValueError(f"{iv.label}: fewer than two measurable rod sensors ({retained})")

    x_m = np.array([pos_cm[c] / 100 for c in retained])
    sx = al.POSITION_UNCERTAINTY_MM / 1000
    sp = al.spatial_fit(retained, x_m, [fits[c] for c in retained], sx)
    kappa, s_kappa, nu, s_nu = al.kappa_nu(w, sp.q, sp.sigma_q, sp.qp, sp.sigma_qp)
    k, s_k, h, s_h = al.material(kappa, s_kappa, nu, s_nu)
    bi, qpR = al.biot_numbers(kappa, nu, sp.qp)
    amp_corr, ph_corr, bound = al.finite_correction(sp.x, sp.q, sp.qp, L)
    cplx_corr = al.finite_complex_correction(sp.x, sp.q, sp.qp, L)
    fsp = None
    if len(retained) >= 3:
        fsp = al.finite_spatial_fit(sp, L, sx)
        fk, s_fk, fn, s_fn = al.kappa_nu_cov(w, fsp.q, fsp.qp, fsp.cov)

    sensors = []
    for c in ROD_COLUMNS:
        if c not in fits:
            continue
        f = fits[c]
        row = {"label": iv.label, "period_s": seg.period, "sensor": SENSOR_LABELS[c],
               "x_cm": pos_cm[c], "mean_C": f.mean, "mean_sigma": f.sigma_mean,
               "amp_C": f.amp, "amp_sigma": f.sigma_amp, "phase_rad": f.phase,
               "phase_sigma": f.sigma_phase, "resid_rms_C": f.rms_resid,
               "resid_amp_C": f.half_pp_resid, "resid_drift_C": f.resid_drift,
               "status": ("retained" if c in retained else
                          "excluded" if SENSOR_LABELS[c] in exclude else "not measurable")}
        if c in sp.columns:
            i = sp.columns.index(c)
            row.update({"phase_unwrapped_rel_rad": sp.phase[i],
                        "finite_amp_corr": amp_corr[i], "finite_phase_corr_rad": ph_corr[i],
                        "finite_bound": bound[i], "finite_complex_corr": cplx_corr[i],
                        "sigma_lnB": sp.sigma_lnB[i], "sigma_phase_fit": sp.sigma_phase[i],
                        "finite_corr_exceeds_uncertainty":
                            bool(abs(amp_corr[i]) > sp.sigma_lnB[i]
                                 or abs(ph_corr[i]) > sp.sigma_phase[i])})
        sensors.append(row)

    # mean profile
    meanfit = None
    if t_inf is not None:
        mcols = [c for c in candidates if c in fits]
        theta = np.array([fits[c].mean - t_inf[c][0] for c in mcols])
        s_theta = np.array([math.sqrt(fits[c].sigma_mean**2 + t_inf[c][1] ** 2
                                      + args.sigma_floor**2) for c in mcols])
        xm = np.array([pos_cm[c] / 100 for c in mcols])
        if len(mcols) >= 3:
            meanfit = al.mean_profile_fit(xm, theta, s_theta, L, sigma_L=args.length_sigma_cm / 100)
            meanfit.x, meanfit.theta, meanfit.sigma, meanfit.cols = xm, theta, s_theta, mcols

    r = {"label": iv.label, "file": seg.path.name, "step": iv.step, "mode": seg.mode,
         "period_s": seg.period, "omega_rad_s": w, "t_start_s": float(tt[0]),
         "t_end_s": float(tt[-1]), "n_cycles": n_cycles,
         "settled": ("auto" if iv.settle else "selected")
         if (not iv.settle or settling.settled_start is not None) else "NOT SETTLED",
         "retained": " ".join(SENSOR_LABELS[c] for c in sp.columns),
         "T1_amp_C": fits["T1_C"].amp if "T1_C" in fits else float("nan"),
         "q_per_m": sp.q, "q_sigma": sp.sigma_q, "q_chi2_dof": _chi2(sp.amp_fit),
         "qp_per_m": sp.qp, "qp_sigma": sp.sigma_qp, "qp_chi2_dof": _chi2(sp.phase_fit),
         "l_amp_cm": 100 / sp.q, "l_amp_sigma_cm": 100 * sp.sigma_q / sp.q**2,
         "l_phase_cm": 100 / sp.qp, "l_phase_sigma_cm": 100 * sp.sigma_qp / sp.qp**2,
         "kappa_m2_s": kappa, "kappa_sigma": s_kappa, "nu_per_s": nu, "nu_sigma": s_nu,
         "k_W_mK": k, "k_sigma": s_k, "h_W_m2K": h, "h_sigma": s_h,
         "Bi_perp": bi, "qprime_R": qpR,
         "max_finite_amp_corr": float(np.max(np.abs(amp_corr))),
         "max_finite_phase_corr_rad": float(np.max(np.abs(ph_corr))),
         "notes": iv.notes}
    if fsp is not None:
        r.update({"q_finite_per_m": fsp.q, "q_finite_sigma": fsp.sigma_q,
                  "qp_finite_per_m": fsp.qp, "qp_finite_sigma": fsp.sigma_qp,
                  "kappa_finite_m2_s": fk, "kappa_finite_sigma": s_fk,
                  "nu_finite_per_s": fn, "nu_finite_sigma": s_fn,
                  "finite_fit_chi2_dof": fsp.chi2 / fsp.dof if fsp.dof > 0 else float("nan")})
    if meanfit is not None:
        r.update({"L_cm": 100 * meanfit.L, "L_sigma_cm": 100 * meanfit.sigma_L,
                  "mu_mean_fit": meanfit.mu, "mu_mean_fit_sigma": meanfit.sigma_mu,
                  "m_mean_fit_sigma_from_mu": meanfit.sigma_m_fit,
                  "m_mean_fit_per_m": meanfit.m, "m_mean_fit_sigma": meanfit.sigma_m,
                  "theta0_K": meanfit.theta0, "mean_fit_chi2_dof":
                      meanfit.chi2 / meanfit.dof if meanfit.dof > 0 else float("nan")})
    # model used for this interval
    exceeds = any(s.get("finite_corr_exceeds_uncertainty") for s in sensors)
    use_finite = fsp is not None and (args.model == "finite" or (args.model == "auto" and exceeds))
    if use_finite:
        r.update({"model_used": "finite", "q_used": fsp.q, "q_used_sigma": fsp.sigma_q,
                  "qp_used": fsp.qp, "qp_used_sigma": fsp.sigma_qp,
                  "kappa_used": fk, "kappa_used_sigma": s_fk, "nu_used": fn, "nu_used_sigma": s_fn})
    else:
        r.update({"model_used": "semi-infinite", "q_used": sp.q, "q_used_sigma": sp.sigma_q,
                  "qp_used": sp.qp, "qp_used_sigma": sp.sigma_qp,
                  "kappa_used": kappa, "kappa_used_sigma": s_kappa, "nu_used": nu,
                  "nu_used_sigma": s_nu})
    ku, s_ku, hu, s_hu = al.material(r["kappa_used"], r["kappa_used_sigma"], r["nu_used"],
                                     r["nu_used_sigma"])
    r.update({"k_used_W_mK": ku, "k_used_sigma": s_ku, "h_used_W_m2K": hu, "h_used_sigma": s_hu})
    r["m_angstrom_per_m"] = (math.sqrt(r["nu_used"] / r["kappa_used"]) if r["nu_used"] > 0
                             else float("nan"))
    r.update({"l_amp_used_cm": 100 / r["q_used"],
              "l_amp_used_sigma_cm": 100 * r["q_used_sigma"] / r["q_used"] ** 2,
              "l_phase_used_cm": 100 / r["qp_used"],
              "l_phase_used_sigma_cm": 100 * r["qp_used_sigma"] / r["qp_used"] ** 2})
    drift_cols = [c for c in ["T1_C"] + retained if c in fits]
    r["max_resid_drift_C"] = max(abs(fits[c].resid_drift) for c in drift_cols)
    if r["max_resid_drift_C"] > args.floor_mean:
        r["notes"] = ((r["notes"] + "; " if r["notes"] else "")
                      + f"residual drift up to {r['max_resid_drift_C']:.2f} C across the fit")
    glitches = al.find_glitches(dd)
    r["glitches"] = len(glitches)
    if glitches:
        r["notes"] = ((r["notes"] + "; " if r["notes"] else "")
                      + f"{len(glitches)} glitch(es) in the fit interval, first at "
                      f"{glitches[0][0]:.0f} s: {glitches[0][1]}")
    if r["n_cycles"] < 5:
        r["notes"] = (r["notes"] + "; " if r["notes"] else "") + f"only {n_cycles} cycles (< 5)"
    if len(retained) < 3:
        r["notes"] = (r["notes"] + "; " if r["notes"] else "") + "fewer than 3 measurable sensors"
    if nu <= 0 or sp.q <= sp.qp:
        r["notes"] = (r["notes"] + "; " if r["notes"] else "") + "q <= q' (nu <= 0): investigate"
    return r, sensors, {"seg": seg, "dd": dd, "fits": fits, "sp": sp, "fsp": fsp,
                        "settling": settling,
                        "meanfit": meanfit, "pos_cm": pos_cm, "retained": retained,
                        "candidates": candidates, "kappa": r["kappa_used"],
                        "nu": r["nu_used"], "L": L}


def _chi2(lf):
    return lf.chi2 / lf.dof if lf.dof > 0 else float("nan")


# ------------------------------------------------------------------
# FIGURES
# ------------------------------------------------------------------
def plot_fits(r, st, fig_dir):
    seg, dd, fits, sp = st["seg"], st["dd"], st["fits"], st["sp"]
    w, P = seg.omega, seg.period
    fig = plt.figure(figsize=(11, 8))
    gs = fig.add_gridspec(3, 2, height_ratios=[2.2, 1, 2])
    ax = fig.add_subplot(gs[0, :])
    t = dd["time_s"]
    m = t >= t[-1] - 3 * P
    for c in TEMPERATURE_COLUMNS:
        if c not in fits:
            continue
        kept = c == "TEC_C" or c in st["retained"]
        color = al.SENSOR_COLORS[c] if kept else al.MUTED
        lab = SENSOR_LABELS[c] + ("" if c == "TEC_C" else f" (x = {st['pos_cm'][c]:.1f} cm)")
        if not kept:
            lab += ", not used"
        ax.plot(t[m], dd[c][m], ls="none", marker=al.SENSOR_MARKERS[c], ms=3, color=color,
                alpha=0.55, markevery=3)
        ax.plot(t[m], fits[c].model(t[m], w), color=color, lw=1.6, label=lab)
    ax.set_ylabel("Temperature (C)")
    ax.set_xlabel("Arduino time (s)")
    ax.set_title(f"{r['label']}: period {P:g} s, {r['n_cycles']} cycles fitted ({r['settled']}); "
                 "markers = data, lines = a + B cos(wt + phi); last 3 cycles shown",
                 fontsize=10, loc="left")
    ax.legend(fontsize=8, frameon=False, ncol=3, loc="upper left", bbox_to_anchor=(1.0, 1.0))

    ax = fig.add_subplot(gs[1, :])
    for c in TEMPERATURE_COLUMNS:
        if c in fits and (c == "TEC_C" or c in st["retained"]):
            f = fits[c]
            ax.plot(f.t, f.resid, color=al.SENSOR_COLORS[c], lw=0.9, label=SENSOR_LABELS[c])
    ax.axhline(0, color=al.MUTED, lw=0.8)
    ax.set_ylabel("Residual (C)")
    ax.set_xlabel("Arduino time (s)")
    ax.legend(fontsize=8, frameon=False, ncol=2, loc="upper left", bbox_to_anchor=(1.0, 1.0))

    xs = np.linspace(0, sp.x.max() * 1.1, 50)
    ax = fig.add_subplot(gs[2, 0])
    for i, c in enumerate(sp.columns):
        ax.errorbar(100 * sp.x[i], sp.lnB[i], yerr=sp.sigma_lnB[i], xerr=0.13, ls="none",
                    marker=al.SENSOR_MARKERS[c], ms=7, color=al.SENSOR_COLORS[c], capsize=2)
    ax.plot(100 * xs, sp.amp_fit.intercept + sp.amp_fit.slope * xs, color=al.INK_2, lw=1.2,
            label="semi-infinite (straight line)")
    fsp = st["fsp"]
    if fsp is not None:
        A = al.finite_amplitude(xs, fsp.q, fsp.qp, st["L"])
        ax.plot(100 * xs, fsp.lnB0 + np.log(np.abs(A)), color=al.REFERENCE_COLOR, lw=1.2,
                ls="--", label="finite insulated tip")
    ax.legend(fontsize=8, frameon=False)
    ax.set_xlabel("x from the first hole (cm)")
    ax.set_ylabel("ln(B / 1 C)")
    ax.set_title(f"q = {sp.q:.2f} ± {sp.sigma_q:.2f} 1/m  (1/q = {100 / sp.q:.2f} cm), "
                 f"chi2/dof = {r['q_chi2_dof']:.2g}", fontsize=9, loc="left")

    ax = fig.add_subplot(gs[2, 1])
    for i, c in enumerate(sp.columns):
        ax.errorbar(100 * sp.x[i], sp.phase[i], yerr=sp.sigma_phase[i], xerr=0.13, ls="none",
                    marker=al.SENSOR_MARKERS[c], ms=7, color=al.SENSOR_COLORS[c], capsize=2,
                    label=SENSOR_LABELS[c])
    ax.plot(100 * xs, sp.phase_fit.intercept + sp.phase_fit.slope * xs, color=al.INK_2, lw=1.2)
    if fsp is not None:
        A = al.finite_amplitude(xs, fsp.q, fsp.qp, st["L"])
        ax.plot(100 * xs, fsp.phi0 + np.unwrap(np.angle(A)), color=al.REFERENCE_COLOR, lw=1.2,
                ls="--")
    ax.set_xlabel("x from the first hole (cm)")
    ax.set_ylabel(f"phi - phi({SENSOR_LABELS[sp.columns[0]]}) (rad), unwrapped")
    ax.set_title(f"q' = {sp.qp:.2f} ± {sp.sigma_qp:.2f} 1/m  (1/q' = {100 / sp.qp:.2f} cm), "
                 f"chi2/dof = {r['qp_chi2_dof']:.2g}", fontsize=9, loc="left")
    ax.legend(fontsize=8, frameon=False)
    for a in fig.axes:
        al.style_axes(a)
    fig.tight_layout()
    out = fig_dir / f"angstrom_{r['label']}_fits.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_profile_theory(r, st, t_inf, fig_dir):
    seg, dd, fits, mf = st["seg"], st["dd"], st["fits"], st["meanfit"]
    w, P, L = seg.omega, seg.period, st["L"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"width_ratios": [1, 1.6]})
    ax = axes[0]
    if mf is not None:
        xs = np.linspace(0, L, 200)
        for i, c in enumerate(mf.cols):
            ax.errorbar(100 * mf.x[i], mf.theta[i], yerr=mf.sigma[i], ls="none",
                        marker=al.SENSOR_MARKERS[c], ms=7, color=al.SENSOR_COLORS[c], capsize=2,
                        label=SENSOR_LABELS[c])
        ax.plot(100 * xs, al.finite_mean(xs, mf.theta0, mf.m, L), color=al.INK_2, lw=1.4,
                label=f"finite fit, m = {mf.m:.2f} ± {mf.sigma_m:.2f} 1/m")
        if r["kappa_used"] > 0 and r["nu_used"] > 0:
            ma = r["m_angstrom_per_m"]
            ax.plot(100 * xs, al.finite_mean(xs, mf.theta0, ma, L), color=al.REFERENCE_COLOR,
                    lw=1.2, ls="--", label=f"m = sqrt(nu/kappa) = {ma:.2f} 1/m")
        ax.axvline(100 * L, color=al.MUTED, lw=0.8, ls=":")
        ax.set_xlabel("x from the first hole (cm)")
        ax.set_ylabel("Mean excess temperature a_i - T_inf (K)")
        ax.legend(fontsize=8, frameon=False)
        ax.set_title("Mean profile (tip at the dotted line)", fontsize=10, loc="left")
    else:
        ax.text(0.5, 0.5, "mean profile needs --room or --baseline", ha="center",
                va="center", transform=ax.transAxes, color=al.INK_2)
        ax.set_axis_off()

    ax = axes[1]
    t = dd["time_s"]
    m = t >= t[-1] - 2 * P
    base = fits["T1_C"]
    t1_inf = t_inf["T1_C"][0] if t_inf is not None else None
    for c in ["T1_C"] + [c for c in st["retained"] if c != "T1_C"]:
        x = st["pos_cm"][c] / 100
        color = al.SENSOR_COLORS[c]
        ax.plot(t[m], dd[c][m], ls="none", marker=al.SENSOR_MARKERS[c], ms=3, color=color,
                alpha=0.6, markevery=2)
        if t1_inf is not None:
            pred = al.predict_rod(x, t[m], t1_inf, base, st["kappa"], st["nu"], w, L)
        else:   # no room temperature: oscillation only, about the measured mean
            q, qp = al.q_qprime(st["kappa"], st["nu"], w)
            A = al.finite_amplitude(x, q, qp, L)
            pred = fits[c].mean + base.amp * abs(A) * np.cos(w * t[m] + base.phase + np.angle(A))
        ax.plot(t[m], pred, color=color, lw=1.6,
                label=f"{SENSOR_LABELS[c]}, x = {100 * x:.1f} cm")
    ax.set_xlabel("Arduino time (s)")
    ax.set_ylabel("Temperature (C)")
    src = "mean + oscillation" if t1_inf is not None else "oscillation about measured means"
    ax.set_title(f"Markers: measured; lines: finite-rod prediction ({src}) from the T1 fit,\n"
                 f"kappa = {st['kappa']:.3g} m^2/s, nu = {st['nu']:.3g} 1/s ({r['model_used']} fit)",
                 fontsize=9, loc="left")
    ax.legend(fontsize=8, frameon=False, ncol=2)
    for a in fig.axes:
        al.style_axes(a)
    fig.tight_layout()
    out = fig_dir / f"angstrom_{r['label']}_profile_theory.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_settling(r, st, fig_dir):
    s = st["settling"]
    fig, axes = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
    for c, wins in s.windows.items():
        ok = [i for i, f in enumerate(wins) if f is not None]
        x = s.window_starts[ok]
        f = [wins[i] for i in ok]
        kw = dict(color=al.SENSOR_COLORS[c], marker=al.SENSOR_MARKERS[c], ms=4, lw=1.2,
                  label=SENSOR_LABELS[c])
        axes[0].plot(x, [g.mean for g in f], **kw)
        axes[1].plot(x, [g.amp for g in f], **kw)
        axes[2].plot(x, np.unwrap([g.phase for g in f]), **kw)
    for a, lab in zip(axes, ["Window mean (C)", "Window amplitude (C)", "Window phase (rad)"]):
        a.set_ylabel(lab)
        if s.settled_start is not None:
            a.axvline(s.settled_start, color=al.INK, lw=1, ls="--")
        al.style_axes(a)
    axes[-1].set_xlabel("Start of the one-cycle window (Arduino time, s)")
    axes[0].legend(fontsize=8, frameon=False, ncol=5)
    when = (f"settled from {s.settled_start:.0f} s (dashed)" if s.settled_start is not None
            else "NOT settled within the interval")
    axes[0].set_title(f"{r['label']}: {when}", fontsize=10, loc="left")
    fig.tight_layout()
    out = fig_dir / f"angstrom_{r['label']}_settling.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_frequency(results, comb, fig_dir, name):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    P = np.array([r["period_s"] for r in results])
    for ax, key, unit, scale in ((axes[0], "kappa", "kappa (1e-5 m^2/s)", 1e5),
                                 (axes[1], "nu", "nu (1e-3 1/s)", 1e3)):
        v = np.array([r[f"{key}_used"] for r in results]) * scale
        s = np.array([r[f"{key}_used_sigma"] for r in results]) * scale
        mean, sig = comb[key][0] * scale, comb[key][1] * scale
        ax.axhspan(mean - sig, mean + sig, color=al.GRID, lw=0)
        ax.axhline(mean, color=al.INK_2, lw=1)
        for model, mk, face in (("semi-infinite", "o", None), ("finite", "s", "white")):
            sel = np.array([r["model_used"] == model for r in results])
            if sel.any():
                ax.errorbar(P[sel], v[sel], yerr=s[sel], ls="none", marker=mk, ms=7,
                            mfc=face, color="#2a78d6", capsize=3, label=f"{model} fit")
        ax.legend(fontsize=8, frameon=False)
        ax.set_xlabel("Forcing period (s)")
        ax.set_ylabel(unit)
        ax.set_xticks(sorted(set(P)))
        ax.set_title(f"weighted mean {mean:.3g} ± {sig:.2g}; chi2/dof = "
                     f"{comb[key][2]:.2g}/{comb[key][3]}, p = {comb[key][4]:.2g}", fontsize=9,
                     loc="left")
        al.style_axes(ax)
    fig.tight_layout()
    out = fig_dir / f"angstrom_{name}_frequency.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_lengths(results, comb, fig_dir, name, length_cm):
    kappa, s_kappa = comb["kappa"][:2]
    nu, s_nu = comb["nu"][:2]
    periods = np.linspace(100, 500, 200)
    central, bands = al.length_curves(periods, kappa, s_kappa, nu, s_nu)
    fig, ax = plt.subplots(figsize=(11, 5.4))
    styles = {"mean": ("#1baf7a", "-", "l_mean = 1/m = sqrt(kappa/nu)"),
              "amp": ("#2a78d6", "-", "l_amp = 1/q (amplitude falls by e)"),
              "phase": ("#eb6834", "--", "l_phase = 1/q' (one radian of lag)")}
    for key, (color, ls, lab) in styles.items():
        ax.fill_between(periods, 100 * bands[key][0], 100 * bands[key][1], color=color,
                        alpha=0.18, lw=0)
        ax.plot(periods, 100 * central[key], color=color, ls=ls, lw=2, label=lab)
    P = np.array([r["period_s"] for r in results])
    ax.errorbar(P, [r["l_amp_used_cm"] for r in results],
                yerr=[r["l_amp_used_sigma_cm"] for r in results],
                ls="none", marker="o", ms=7, color="#2a78d6", capsize=3, label="1/q fitted")
    ax.errorbar(P, [r["l_phase_used_cm"] for r in results],
                yerr=[r["l_phase_used_sigma_cm"] for r in results], ls="none", marker="^", ms=7,
                color="#eb6834", capsize=3, label="1/q' fitted")
    mm = [r for r in results if not math.isnan(r.get("m_mean_fit_per_m", float("nan")))]
    if mm:
        m_mean, m_sig, *_ = al.weighted_mean([r["m_mean_fit_per_m"] for r in mm],
                                              [r["m_mean_fit_sigma"] for r in mm])
        ax.errorbar([112], [100 / m_mean], yerr=[100 * m_sig / m_mean**2],
                    ls="none", marker="D", ms=7, color="#1baf7a", capsize=3,
                    label="1/m from the Module 8 mean profile")
    ax.errorbar([106], [100 / al.MODULE7_M], yerr=[100 * al.MODULE7_M_SIGMA / al.MODULE7_M**2],
                ls="none", marker="s", ms=7, mfc="white", color="#1baf7a", capsize=3,
                label="1/m, Module 7 stationary profile")
    ax.axhline(length_cm, color=al.MUTED, lw=1, ls=":")
    ax.text(498, length_cm, f"analysis length L = {length_cm:.1f} cm", ha="right", va="bottom",
            fontsize=8, color=al.INK_2)
    ax.set_xticks([100, 200, 300, 400, 500])
    ax.set_xlim(90, 510)
    ax.set_ylim(0, max(length_cm, 100 * central["mean"][0]) * 1.25)
    ax.set_xlabel("Forcing period (s)")
    ax.set_ylabel("Length (cm)")
    ax.set_title(f"Characteristic lengths from the measured kappa = {kappa:.3g} ± {s_kappa:.2g} "
                 f"m^2/s, nu = {nu:.3g} ± {s_nu:.2g} 1/s (bands: 16-84 %)", fontsize=9, loc="left")
    ax.legend(fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    al.style_axes(ax)
    fig.tight_layout()
    out = fig_dir / f"angstrom_{name}_lengths.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


# ------------------------------------------------------------------
# REPORT
# ------------------------------------------------------------------
def print_interval(r, sensors, st):
    print(f"\n=== {r['label']}: {r['file']}, step {r['step']}, period {r['period_s']:g} s, "
          f"{r['t_start_s']:.0f}-{r['t_end_s']:.0f} s, {r['n_cycles']} cycles ({r['settled']})")
    if st["settling"] is not None:
        print(f"settling criterion: {st['settling'].criterion}")
    print("| Sensor | x (cm) | a (C) | B (C) | phi (rad) | resid amp (C) | status | finite amp / phase corr |")
    print("| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |")
    for s in sensors:
        corr = ""
        if "finite_amp_corr" in s:
            corr = (f"{s['finite_amp_corr']:+.1%} / {s['finite_phase_corr_rad']:+.3f} rad"
                    + (" (exceeds uncertainty)" if s["finite_corr_exceeds_uncertainty"] else ""))
        print(f"| {s['sensor']} | {s['x_cm']:.1f} | {s['mean_C']:.2f} ± {s['mean_sigma']:.2f} | "
              f"{s['amp_C']:.3f} ± {s['amp_sigma']:.3f} | {s['phase_rad']:+.3f} ± {s['phase_sigma']:.3f} | "
              f"{s['resid_amp_C']:.3f} | {s['status']} | {corr} |")
    print(f"q  = {r['q_per_m']:.3f} ± {r['q_sigma']:.3f} 1/m (chi2/dof {r['q_chi2_dof']:.2g}), "
          f"1/q  = {r['l_amp_cm']:.2f} ± {r['l_amp_sigma_cm']:.2f} cm")
    print(f"q' = {r['qp_per_m']:.3f} ± {r['qp_sigma']:.3f} 1/m (chi2/dof {r['qp_chi2_dof']:.2g}), "
          f"1/q' = {r['l_phase_cm']:.2f} ± {r['l_phase_sigma_cm']:.2f} cm")
    print(f"kappa = {r['kappa_m2_s']:.4g} ± {r['kappa_sigma']:.2g} m^2/s,  "
          f"nu = {r['nu_per_s']:.4g} ± {r['nu_sigma']:.2g} 1/s")
    if "kappa_finite_m2_s" in r:
        print(f"finite insulated-tip spatial fit (model check): q = {r['q_finite_per_m']:.3f} ± "
              f"{r['q_finite_sigma']:.3f}, q' = {r['qp_finite_per_m']:.3f} ± {r['qp_finite_sigma']:.3f} 1/m; "
              f"kappa = {r['kappa_finite_m2_s']:.4g} ± {r['kappa_finite_sigma']:.2g} m^2/s, "
              f"nu = {r['nu_finite_per_s']:.4g} ± {r['nu_finite_sigma']:.2g} 1/s "
              f"(chi2/dof {r['finite_fit_chi2_dof']:.2g})")
    print(f"(with rho = {al.RHO:g} kg/m^3, c = {al.C_P:g} J/(kg K): k = {r['k_W_mK']:.4g} ± "
          f"{r['k_sigma']:.2g} W/(m K), h = {r['h_W_m2K']:.3g} ± {r['h_sigma']:.2g} W/(m^2 K))")
    print(f"checks: Bi_perp = {r['Bi_perp']:.2e}, q'R = {r['qprime_R']:.3f}, largest finite-rod "
          f"correction {r['max_finite_amp_corr']:.2%} in amplitude, "
          f"{r['max_finite_phase_corr_rad']:.4f} rad in phase")
    print(f"used (k, h from these): k = {r['k_used_W_mK']:.4g} ± {r['k_used_sigma']:.2g} W/(m K), "
          f"h = {r['h_used_W_m2K']:.3g} ± {r['h_used_sigma']:.2g} W/(m^2 K)")
    print(f"used: {r['model_used']} fit, kappa = {r['kappa_used']:.4g} ± {r['kappa_used_sigma']:.2g} "
          f"m^2/s, nu = {r['nu_used']:.4g} ± {r['nu_used_sigma']:.2g} 1/s; largest residual drift "
          f"{r['max_resid_drift_C']:.3f} C")
    if "m_mean_fit_per_m" in r:
        print(f"mean profile: mu = mL = {r['mu_mean_fit']:.3f} ± {r['mu_mean_fit_sigma']:.3f} "
              f"(L = {r['L_cm']:.2f} ± {r['L_sigma_cm']:.2f} cm) -> "
              f"m = mu/L = {r['m_mean_fit_per_m']:.3f} ± {r['m_mean_fit_sigma']:.3f} 1/m "
              f"(chi2/dof {r['mean_fit_chi2_dof']:.2g}); sqrt(nu/kappa) = "
              f"{r['m_angstrom_per_m']:.3f} 1/m; Module 7: {al.MODULE7_M} ± {al.MODULE7_M_SIGMA} 1/m")
    if r["notes"]:
        print(f"notes: {r['notes']}")


def print_finite_table(sensors):
    """Stage 6 finite-length check at every retained sensor: the bound
    (e^{-2q(L-x)} + e^{-2qL}) / (1 - e^{-2qL}), the actual |H_fin/H_semi - 1| with the
    fitted q, q', and the amplitude and phase corrections against the measurement
    uncertainty of ln B and phi."""
    rows = [s for s in sensors if "finite_bound" in s]
    if not rows:
        return
    print("\nFinite-length check (semi-infinite q, q' of each interval):")
    print("| Interval | Period (s) | Sensor | x (cm) | Bound | abs(H_fin/H_semi - 1) | "
          "Amplitude corr | sigma(ln B) | Phase corr (rad) | sigma(phi) (rad) | Semi-infinite OK? |")
    print("| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |")
    for s in rows:
        ok = not s["finite_corr_exceeds_uncertainty"]
        print(f"| {s['label']} | {s['period_s']:g} | {s['sensor']} | {s['x_cm']:.1f} | "
              f"{s['finite_bound']:.2e} | {s['finite_complex_corr']:.2e} | "
              f"{s['finite_amp_corr']:+.2e} | {s['sigma_lnB']:.1e} | {s['finite_phase_corr_rad']:+.1e} | "
              f"{s['sigma_phase_fit']:.1e} | {'yes' if ok else 'no: correction > uncertainty'} |")


def write_rows(rows, path):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


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
    ap = argparse.ArgumentParser(description="Module 8 Angstrom analysis")
    ap.add_argument("run", nargs="?", help="raw run CSV")
    ap.add_argument("--intervals", help="intervals CSV (label,file,part,step,t_start_s,t_end_s,settle,exclude,notes)")
    ap.add_argument("--step", type=int)
    ap.add_argument("--t", help="START,END in Arduino seconds")
    ap.add_argument("--part", type=int)
    ap.add_argument("--settle", action="store_true")
    ap.add_argument("--label")
    ap.add_argument("--name")
    ap.add_argument("--exclude", type=lambda s: [v.strip().upper() for v in s.split(",") if v.strip()])
    ap.add_argument("--positions-cm", type=lambda s: [float(v) for v in s.split(",")])
    ap.add_argument("--length-cm", type=float, default=100 * al.ROD_ANALYSIS_LENGTH_M)
    ap.add_argument("--length-sigma-cm", type=float, default=0.14,
                    help="uncertainty of L: ruler 0.05 cm combined with the T1 thermistor's "
                         "1.3 mm placement in its hole")
    ap.add_argument("--room", type=float, help="room temperature T_inf (C) for the mean profile")
    ap.add_argument("--baseline", help="FILE:START,END with the TEC off; per-channel T_inf")
    ap.add_argument("--min-amp", type=float, default=0.05, help="smallest usable amplitude (C)")
    ap.add_argument("--min-snr", type=float, default=3.0)
    ap.add_argument("--sigma-floor", type=float, default=0.05, help="floor on mean uncertainty (K)")
    ap.add_argument("--k-sigma", type=float, default=3.0)
    ap.add_argument("--floor-mean", type=float, default=0.05)
    ap.add_argument("--floor-amp", type=float, default=0.05)
    ap.add_argument("--floor-phase", type=float, default=0.03)
    ap.add_argument("--theory", help="label of the interval for the theory figure (default: all)")
    ap.add_argument("--model", choices=("auto", "semi", "finite"), default="auto",
                    help="spatial model for the values used (see header)")
    ap.add_argument("--fig-dir", default=str(al.MODULE8_FIG_DIR))
    ap.add_argument("--out-dir", default=str(al.ANALYSIS_DIR))
    args = ap.parse_args(argv)
    if args.positions_cm is not None and len(args.positions_cm) != 5:
        raise SystemExit("--positions-cm needs five values (T1-T5)")

    intervals, default_name = intervals_from_args(args)
    name = args.name or default_name
    fig_dir, out_dir = Path(args.fig_dir), Path(args.out_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    if args.baseline:
        t_inf = al.baseline_temperatures(args.baseline)
        print(f"T_inf per channel from {args.baseline}")
    elif args.room is not None:
        t_inf = {c: (args.room, 0.0) for c in TEMPERATURE_COLUMNS}
    else:
        t_inf = None
        print("no --room or --baseline: the mean-profile fit is skipped")

    results, all_sensors = [], []
    for iv in intervals:
        seg = al.select_segment(iv)
        iv_inf = t_inf
        if iv.baseline:
            iv_inf = None if iv.baseline.lower() == "none" else al.baseline_temperatures(iv.baseline)
            print(f"\n{iv.label}: T_inf per channel from {iv.baseline}")
        r, sensors, st = analyse(seg, args, iv_inf)
        print_interval(r, sensors, st)
        print(f"figures: {plot_fits(r, st, fig_dir)}")
        if args.theory is None or args.theory == r["label"]:
            print(f"         {plot_profile_theory(r, st, iv_inf, fig_dir)}")
        if st["settling"] is not None:
            print(f"         {plot_settling(r, st, fig_dir)}")
        results.append(r)
        all_sensors += sensors

    good = [r for r in results if r["nu_used"] > 0 and r["kappa_used"] > 0]
    comb = {}
    if good:
        comb["kappa"] = al.weighted_mean([r["kappa_used"] for r in good],
                                         [r["kappa_used_sigma"] for r in good])
        comb["nu"] = al.weighted_mean([r["nu_used"] for r in good],
                                      [r["nu_used_sigma"] for r in good])
    print("\nSemi-infinite two-stage fit (primary) and the values used (--model "
          f"{args.model}):")
    print("| Interval | Period (s) | Cycles | Sensors | q (1/m) | q' (1/m) | kappa (1e-5 m^2/s) | "
          "nu (1e-3 1/s) | model used | kappa used | nu used | m mean fit (1/m) |")
    print("| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: |")
    for r in results:
        print(f"| {r['label']} | {r['period_s']:g} | {r['n_cycles']} | {r['retained']} | "
              f"{r['q_per_m']:.2f} ± {r['q_sigma']:.2f} | {r['qp_per_m']:.2f} ± {r['qp_sigma']:.2f} | "
              f"{1e5 * r['kappa_m2_s']:.3f} ± {1e5 * r['kappa_sigma']:.3f} | "
              f"{1e3 * r['nu_per_s']:.3f} ± {1e3 * r['nu_sigma']:.3f} | {r['model_used']} | "
              f"{1e5 * r['kappa_used']:.3f} ± {1e5 * r['kappa_used_sigma']:.3f} | "
              f"{1e3 * r['nu_used']:.3f} ± {1e3 * r['nu_used_sigma']:.3f} | "
              f"{fnum(r.get('m_mean_fit_per_m'), '.2f')} |")
    if len(good) > 1:
        for key in ("kappa", "nu"):
            mean, sig, chi2, dof, p = comb[key]
            print(f"weighted mean {key} = {mean:.4g} ± {sig:.2g}; chi2 = {chi2:.2f} for {dof} dof, "
                  f"p = {p:.3g}" + ("  (periods DISAGREE beyond uncertainty)" if p < 0.05 else ""))
        print(f"figure: {plot_frequency(good, comb, fig_dir, name)}")
    if good:
        print(f"figure: {plot_lengths(good, comb, fig_dir, name, args.length_cm)}")
    print_finite_table(all_sensors)
    write_rows(results, out_dir / f"angstrom_{name}.csv")
    write_rows(all_sensors, out_dir / f"angstrom_{name}_sensors.csv")
    print(f"results: {out_dir / f'angstrom_{name}.csv'}, {out_dir / f'angstrom_{name}_sensors.csv'}")
    return results, all_sensors, comb


if __name__ == "__main__":
    main()
