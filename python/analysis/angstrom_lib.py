# Phys 39 Module 8 - analysis library shared by the control and Angstrom scripts
#
# Loading
#   load_parts(csv)        raw run CSV, split where the Arduino clock restarts
#   read_metadata(csv)     sensor positions, events, and settings from <run>_meta.txt
#   read_intervals(file)   list of analysis intervals (data/module_08/intervals.csv)
#   select_segment(iv)     the records of one interval, with its period and mode
#
# Common time fit (course page, "Common Time Fit For The Stage Analysis")
#   T_i(t) = a_i + B_i cos(omega t + phi_i)
#   fit_sine() solves the linear form a + c cos(wt) + s sin(wt) by least squares:
#   B = sqrt(c^2 + s^2), phi = atan2(-s, c). A downstream lag gives a smaller phi.
#   The regression also carries the 2nd and 3rd harmonics as nuisance terms. Over whole
#   cycles they are orthogonal to the fundamental, so a, B, and phi are unchanged, but a
#   distorted waveform (zero-centered PWM with unequal heating and cooling) no longer
#   counts as noise: the uncertainties come from the residual of the full model, scaled
#   and inflated by (1 + rho)/(1 - rho), rho = lag-1 autocorrelation of that residual, so
#   slow drift or correlated noise is not counted as independent samples. The reported
#   residual (resid, half_pp_resid, rms_resid) is the course's: data minus the fitted
#   mean and fundamental, i.e. it includes the harmonic distortion.
#
# Settling (S18.9)
#   settling_test() fits successive one-cycle windows. Two windows agree when their
#   mean, amplitude, and phase differ by at most max(k sigma_difference, floor) for
#   every tested channel. sigma is the larger of the fit uncertainty and the measured
#   cycle-to-cycle scatter (robust spread of neighbour differences over the second half
#   of the windows): air movement makes the rod amplitudes scatter by 1-3 % from cycle
#   to cycle, several times the fit uncertainty, and that is noise, not drift. The response is settled from the first window after which
#   every window agrees with its neighbour (the course criterion) and with the last
#   window (so a slow drift that stays below the per-cycle tolerance is still caught).
#
# Spatial fit and parameters (Stage 6)
#   ln B_i = c - q x_i,   phi_i (unwrapped along x) = c' - q' x_i
#   weighted straight-line fits with the effective-variance method for the position
#   uncertainty; the slope uncertainty is multiplied by sqrt(chi2/dof) when that
#   exceeds 1.
#   kappa = omega / (2 q q'),   nu = kappa (q^2 - q'^2)
#
# Glitches (loose breadboard contacts)
#   find_glitches()        records where the second difference T[k] - 2T[k-1] + T[k-2]
#                          exceeds GLITCH_COMMON_C on >= 3 channels with the same sign,
#                          or GLITCH_SINGLE_C on a downstream rod channel (T2-T5).
#                          A smooth record stays far below both: at most ~0.06 C for
#                          the TEC at 100 s and 16 C amplitude, 0.15 C for one ADC
#                          count of flicker.
#
# Model checks
#   finite_correction()    finite insulated-tip vs semi-infinite complex amplitude
#   biot_numbers()         Bi_perp = h d / (4k) = nu d^2 / (16 kappa) and q' R
#   mean_profile_fit()     theta(x) = theta0 cosh[m(L - x)] / cosh(mL)
#   predict_rod()          T(x, t) from the measured T1 fit and kappa, nu
#   length_curves()        1/m, 1/q, 1/q' vs period with Monte Carlo bands

import csv
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2 as chi2_dist

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
# MODULE8_FIG_DIR and POSITION_UNCERTAINTY_MM are re-exported for the analysis scripts.
from module_08_common import (MODULE8_DATA_DIR, MODULE8_FIG_DIR, POSITION_UNCERTAINTY_MM,  # noqa: E402
                              REPO_ROOT, ROD_ANALYSIS_LENGTH_M, ROD_COLUMNS, ROD_DIAMETER_M,
                              SENSOR_LABELS, SENSOR_POSITIONS_CM, TEMPERATURE_COLUMNS,
                              read_run_csv)

ANALYSIS_DIR = MODULE8_DATA_DIR / "analysis"

# 6061-T6 properties used on the course page; needed only to turn kappa and nu into
# k and h, not for kappa and nu themselves.
RHO = 2700.0          # kg/m^3
C_P = 896.0           # J/(kg K)

# Module 7 mean-profile result (module_07_stationary_fin.md), for comparison.
MODULE7_M = 4.90      # 1/m
MODULE7_M_SIGMA = 0.04

# Figure styling: categorical slots of the reference palette (validated for CVD
# separation on a light surface); each sensor also has its own marker, so identity
# never rests on color alone.
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
SENSOR_COLORS = {"TEC_C": INK, "T1_C": "#2a78d6", "T2_C": "#eb6834", "T3_C": "#1baf7a",
                 "T4_C": "#eda100", "T5_C": "#e87ba4"}
SENSOR_MARKERS = {"TEC_C": "s", "T1_C": "o", "T2_C": "^", "T3_C": "D", "T4_C": "v", "T5_C": "P"}
REFERENCE_COLOR = "#6250d6"     # setpoint / PWM reference (slot 7, violet)


def style_axes(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=INK_2, labelsize=8)
    ax.xaxis.label.set_color(INK_2)
    ax.yaxis.label.set_color(INK_2)
    ax.title.set_color(INK)


# ------------------------------------------------------------------
# LOADING
# ------------------------------------------------------------------
def resolve_path(name: str, base: Path = None) -> Path:
    """A data file given absolute, relative to the working directory, to `base`, to
    the repository root, or to data/module_08/."""
    p = Path(name)
    candidates = [p] if p.is_absolute() else [Path.cwd() / p]
    if base is not None:
        candidates.append(base / p)
    candidates += [REPO_ROOT / p, MODULE8_DATA_DIR / p]
    for c in candidates:
        if c.exists():
            return c.resolve()
    raise FileNotFoundError(f"{name}: not found (looked in {', '.join(str(c) for c in candidates)})")


def load_parts(path: Path):
    """The run's records, split into parts wherever time_ms decreases (Arduino restart)."""
    data = read_run_csv(path)
    t = data["time_ms"]
    breaks = [0] + [i for i in range(1, len(t)) if t[i] < t[i - 1]] + [len(t)]
    parts = []
    for a, b in zip(breaks[:-1], breaks[1:]):
        parts.append({k: v[a:b] for k, v in data.items()})
    return parts


_POS_RE = re.compile(r"^\s+(T\d_C):.*x = ([-\d.]+) cm")
_EVENT_RE = re.compile(r"^\s+\[\s*([-\d.]+) s \|[^\]]*\]\s*(.*)$")
_KV_RE = re.compile(r"(\w+)=([-\w.]+)")
_SET_RE = re.compile(r"\bSET (\w+) ([-\d.]+)")
# SET command names -> the names used in START/CONFIG lines
SET_NAMES = {"MEAN": "Ta_C", "AMP": "Tb_C", "PERIOD": "period_s", "KP": "Kp", "KI": "Ki",
             "UMAX": "u_max", "RATIO": "r", "CHIH": "chi_H", "TROOM": "T_room_C",
             "PWMAMP": "u_b", "PWMOFF": "u_a", "MANUAL": "manual_pwm", "FFADD": "ff_add"}


def read_metadata(csv_path: Path):
    """Sensor positions (cm) and events [(part, time_s, text)] from <run>_meta.txt.
    `part` counts Arduino restarts, as in load_parts()."""
    meta_path = csv_path.with_name(csv_path.stem + "_meta.txt")
    positions, events = {}, []
    if not meta_path.exists():
        return {"positions_cm": positions, "events": events, "path": None}
    part = 0
    for line in meta_path.read_text().splitlines():
        m = _POS_RE.match(line)
        if m:
            positions[m.group(1)] = float(m.group(2))
            continue
        m = _EVENT_RE.match(line)
        if m:
            text = m.group(2)
            if "Arduino restarted" in text:
                part += 1
            events.append((part, float(m.group(1)), text))
    return {"positions_cm": positions, "events": events, "path": meta_path}


def settings_at(meta, part, t):
    """Arduino settings in force at Arduino time t of the given part, from the START,
    CONFIG, SET, and MODE events logged up to then (strings; empty if unknown)."""
    out = {}
    for p, te, text in meta.get("events", []):
        if (p, te) > (part, t):
            break
        if "START " in text or "CONFIG " in text:
            out.update(dict(_KV_RE.findall(text)))
        m = _SET_RE.search(text)
        if m and m.group(1) in SET_NAMES:
            out[SET_NAMES[m.group(1)]] = m.group(2)
        m = re.search(r"\bMODE (\w+)", text)
        if m:
            out["mode"] = m.group(1)
    return out


@dataclass
class Interval:
    label: str
    file: str
    part: int = None
    step: int = None
    t_start: float = None
    t_end: float = None
    settle: bool = False
    exclude: tuple = ()
    notes: str = ""
    base: Path = None
    baseline: str = None    # room baseline for this interval (see baseline_temperatures)


INTERVAL_COLUMNS = ["label", "file", "part", "step", "t_start_s", "t_end_s", "settle",
                    "exclude", "notes", "baseline"]   # baseline is optional


def _opt(text, kind):
    text = (text or "").strip()
    return kind(text) if text else None


def read_intervals(path: Path):
    """Intervals file: one row per analysis interval (columns INTERVAL_COLUMNS).
    Empty cells mean 'not used'; settle is yes/no; exclude lists sensors like T5 or T4;T5."""
    out = []
    with open(path, newline="") as f:
        rows = [r for r in csv.DictReader(f) if (r.get("label") or "").strip()
                and not r["label"].strip().startswith("#")]
    for r in rows:
        excl = tuple(s.strip().upper() for s in re.split(r"[;\s]+", r.get("exclude") or "")
                     if s.strip())
        out.append(Interval(
            label=r["label"].strip(), file=r["file"].strip(),
            part=_opt(r.get("part"), int), step=_opt(r.get("step"), int),
            t_start=_opt(r.get("t_start_s"), float), t_end=_opt(r.get("t_end_s"), float),
            settle=(r.get("settle") or "").strip().lower() in ("yes", "y", "true", "1"),
            exclude=excl, notes=(r.get("notes") or "").strip(), base=path.parent,
            baseline=(r.get("baseline") or "").strip() or None))
    return out


@dataclass
class Segment:
    interval: Interval
    path: Path
    data: dict              # arrays restricted to the interval
    step_data: dict         # the whole step(s) of the interval, for startup checks
    period: float
    omega: float
    mode: str
    meta: dict
    part: int

    def settings(self):
        return settings_at(self.meta, self.part, float(self.data["time_s"][0]))


def select_segment(iv: Interval) -> Segment:
    path = resolve_path(iv.file, iv.base)
    parts = load_parts(path)

    def mask_of(p, use_time):
        m = np.ones(len(p["time_ms"]), bool)
        if iv.step is not None:
            m &= p["step"] == iv.step
        if use_time and iv.t_start is not None:
            m &= p["time_s"] >= iv.t_start
        if use_time and iv.t_end is not None:
            m &= p["time_s"] <= iv.t_end
        return m

    if iv.part is not None:
        if not 0 <= iv.part < len(parts):
            raise ValueError(f"{iv.label}: part {iv.part} does not exist ({len(parts)} parts)")
        chosen = iv.part
    else:
        matches = [i for i, p in enumerate(parts) if mask_of(p, True).any()]
        if not matches:
            raise ValueError(f"{iv.label}: no records match step={iv.step}, "
                             f"t={iv.t_start}..{iv.t_end} s in {path.name}")
        if len(matches) > 1:
            raise ValueError(f"{iv.label}: the Arduino restarted in {path.name}; records in "
                             f"parts {matches} match. Give the 'part' column.")
        chosen = matches[0]
    p = parts[chosen]
    m = mask_of(p, True)
    if not m.any():
        raise ValueError(f"{iv.label}: the selection is empty")
    data = {k: v[m] for k, v in p.items()}
    step_mask = mask_of(p, False) if iv.step is not None else m
    step_data = {k: v[step_mask] for k, v in p.items()}

    periods = np.unique(data["period_s"][np.isfinite(data["period_s"])])
    if len(periods) != 1:
        raise ValueError(f"{iv.label}: the interval must have one forcing period; found "
                         f"{list(periods)} (select one step)")
    modes = np.unique(data["control_mode"])
    if len(modes) != 1:
        raise ValueError(f"{iv.label}: the interval must have one control mode; found {list(modes)}")
    states = set(np.unique(data["run_state"]))
    if states != {"RUN"}:
        raise ValueError(f"{iv.label}: the interval includes run_state {sorted(states)}; "
                         "select records with run_state RUN only")
    period = float(periods[0])
    return Segment(iv, path, data, step_data, period, 2 * math.pi / period, str(modes[0]),
                   read_metadata(path), chosen)


def complete_cycles(t, period, start=None):
    """Mask of the largest whole number of periods from `start` (default: first time).
    A span of n periods sampled once per record ends one record interval short of
    t0 + n P, so that interval is allowed for."""
    t0 = t[0] if start is None else start
    dt = float(np.median(np.diff(t))) if len(t) > 1 else 0.0
    n = int(math.floor((t[-1] - t0 + dt) / period + 1e-6))
    return (t >= t0) & (t < t0 + n * period - 1e-9), n


# ------------------------------------------------------------------
# COMMON TIME FIT
# ------------------------------------------------------------------
@dataclass
class SineFit:
    mean: float
    amp: float
    phase: float            # rad, in (-pi, pi]
    sigma_mean: float
    sigma_amp: float
    sigma_phase: float
    rms_resid: float
    half_pp_resid: float    # residual amplitude: half the residual peak-to-peak
    n: int
    inflation: float        # variance inflation for correlated residuals
    t: np.ndarray = field(repr=False)
    resid: np.ndarray = field(repr=False)

    @property
    def resid_drift(self):
        """Change of the residual's straight-line trend across the fitted interval (C):
        drift left in the fit window."""
        slope = np.polyfit(self.t, self.resid, 1)[0]
        return float(slope * (self.t[-1] - self.t[0]))

    def model(self, t, omega):
        return self.mean + self.amp * np.cos(omega * t + self.phase)


def fit_sine(t, y, omega, harmonics=3) -> SineFit:
    ok = np.isfinite(t) & np.isfinite(y)
    t, y = t[ok], y[ok]
    cols = [np.ones_like(t), np.cos(omega * t), np.sin(omega * t)]
    for k in range(2, harmonics + 1):
        cols += [np.cos(k * omega * t), np.sin(k * omega * t)]
    if len(t) < len(cols) + 5:
        raise ValueError("too few points for a sinusoid fit")
    X = np.column_stack(cols)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    noise = y - X @ beta                          # residual of the full model
    resid = y - X[:, :3] @ beta[:3]               # data - mean - fundamental (course residual)
    dof = len(t) - X.shape[1]
    s2 = float(noise @ noise) / dof
    if len(noise) > 2 and np.std(noise) > 0:
        rho = float(np.corrcoef(noise[:-1], noise[1:])[0, 1])
        rho = min(max(rho, 0.0), 0.95)
    else:
        rho = 0.0
    inflation = (1 + rho) / (1 - rho)
    cov = (s2 * inflation * np.linalg.inv(X.T @ X))[:3, :3]
    a, c, s = beta[:3]
    B = math.hypot(c, s)
    phi = math.atan2(-s, c)
    jB = np.array([0.0, c / B, s / B])
    jphi = np.array([0.0, s / B**2, -c / B**2])
    return SineFit(mean=float(a), amp=B, phase=phi,
                   sigma_mean=math.sqrt(cov[0, 0]), sigma_amp=math.sqrt(jB @ cov @ jB),
                   sigma_phase=math.sqrt(jphi @ cov @ jphi),
                   rms_resid=math.sqrt(float(np.mean(resid**2))),
                   half_pp_resid=0.5 * float(np.ptp(resid)), n=len(t), inflation=inflation,
                   t=t, resid=resid)


def wrap(angle):
    """Wrap to (-pi, pi]."""
    return -((-angle + math.pi) % (2 * math.pi) - math.pi)


def lag(ref: SineFit, resp: SineFit):
    """Phase lag of resp behind ref (rad, wrapped to (-pi, pi]) and its uncertainty."""
    return wrap(ref.phase - resp.phase), math.hypot(ref.sigma_phase, resp.sigma_phase)


# ------------------------------------------------------------------
# GLITCHES
# ------------------------------------------------------------------
GLITCH_COMMON_C = 0.25
GLITCH_SINGLE_C = 0.5
GLITCH_SINGLE_COLUMNS = ["T2_C", "T3_C", "T4_C", "T5_C"]


def glitch_at(prev2, prev1, now):
    """Glitch test for the newest of three consecutive records (dicts of the six
    temperatures). Returns a short description, or None."""
    d2 = {c: now[c] - 2 * prev1[c] + prev2[c] for c in TEMPERATURE_COLUMNS}
    d2 = {c: v for c, v in d2.items() if np.isfinite(v)}
    up = [c for c, v in d2.items() if v > GLITCH_COMMON_C]
    down = [c for c, v in d2.items() if v < -GLITCH_COMMON_C]
    for group in (up, down):
        if len(group) >= 3:
            return ("common-mode jump on " + ", ".join(SENSOR_LABELS[c] for c in group)
                    + f" ({max(abs(d2[c]) for c in group):.2f} C)")
    single = [c for c in GLITCH_SINGLE_COLUMNS if abs(d2.get(c, 0.0)) > GLITCH_SINGLE_C]
    if single:
        return "jump on " + ", ".join(f"{SENSOR_LABELS[c]} ({d2[c]:+.2f} C)" for c in single)
    return None


def find_glitches(data):
    """[(time_s, description)] for every glitch in a run segment."""
    out = []
    n = len(data["time_s"])
    for k in range(2, n):
        rec = [{c: data[c][j] for c in TEMPERATURE_COLUMNS} for j in (k - 2, k - 1, k)]
        g = glitch_at(*rec)
        if g:
            out.append((float(data["time_s"][k]), g))
    return out


# ------------------------------------------------------------------
# SETTLING
# ------------------------------------------------------------------
@dataclass
class SettlingResult:
    settled_start: float          # s, or None if never settled
    window_starts: np.ndarray
    windows: dict                 # column -> list of SineFit (one per window)
    passes: np.ndarray            # passes[k]: comparison of window k with window k-1 (k >= 1)
    criterion: str


def settling_test(data, columns, period, k_sigma=3.0, floor_mean=0.05, floor_amp=0.05,
                  floor_phase=0.03):
    omega = 2 * math.pi / period
    t = data["time_s"]
    n_win = int(math.floor((t[-1] - t[0]) / period + 1e-6))
    starts = t[0] + period * np.arange(n_win)
    windows = {c: [] for c in columns}
    for w0 in starts:
        m = (t >= w0) & (t < w0 + period)
        for c in columns:
            try:
                windows[c].append(fit_sine(t[m], data[c][m], omega))
            except ValueError:
                windows[c].append(None)
    # measured cycle-to-cycle scatter per channel and parameter (sigma of one window)
    scatter = {}
    for c in columns:
        fits = [f for f in windows[c][n_win // 2:] if f is not None]
        sc = {}
        for name, get in (("mean", lambda f: f.mean), ("amp", lambda f: f.amp),
                          ("phase", lambda f: f.phase)):
            v = np.unwrap([get(f) for f in fits]) if name == "phase" else np.array([get(f) for f in fits])
            d = np.diff(v)
            sc[name] = 1.4826 * float(np.median(np.abs(d - np.median(d)))) / math.sqrt(2) if len(d) >= 4 else 0.0
        scatter[c] = sc

    def agree(c, a, b):
        if a is None or b is None:
            return False
        sc = scatter[c]
        sm = math.hypot(max(a.sigma_mean, sc["mean"]), max(b.sigma_mean, sc["mean"]))
        sa = math.hypot(max(a.sigma_amp, sc["amp"]), max(b.sigma_amp, sc["amp"]))
        sp = math.hypot(max(a.sigma_phase, sc["phase"]), max(b.sigma_phase, sc["phase"]))
        return (abs(b.mean - a.mean) <= max(k_sigma * sm, floor_mean)
                and abs(b.amp - a.amp) <= max(k_sigma * sa, floor_amp)
                and abs(wrap(b.phase - a.phase)) <= max(k_sigma * sp, floor_phase))

    passes = np.zeros(n_win, bool)        # window k agrees with window k - 1
    with_last = np.zeros(n_win, bool)     # window k agrees with the last window
    for k in range(n_win):
        if k > 0:
            passes[k] = all(agree(c, windows[c][k - 1], windows[c][k]) for c in columns)
        with_last[k] = all(agree(c, windows[c][k], windows[c][-1]) for c in columns)
    settled = None
    for k in range(n_win - 1):
        if passes[k + 1:].all() and with_last[k:].all():
            settled = float(starts[k])
            break
    crit = (f"every later one-cycle window agrees with its neighbour and with the last window "
            f"in mean, amplitude, and phase within max({k_sigma:g} sigma, floor), sigma = larger "
            f"of fit uncertainty and measured cycle-to-cycle scatter, for "
            f"{', '.join(SENSOR_LABELS.get(c, c) for c in columns)}; floors {floor_mean:g} C, "
            f"{floor_amp:g} C, {floor_phase:g} rad")
    return SettlingResult(settled, starts, windows, passes, crit)


# ------------------------------------------------------------------
# SPATIAL FIT, KAPPA, NU
# ------------------------------------------------------------------
@dataclass
class LineFit:
    slope: float
    intercept: float
    sigma_slope: float       # statistical, from the weights
    sigma_intercept: float
    chi2: float
    dof: int
    scale: float             # max(1, sqrt(chi2/dof))
    resid: np.ndarray = field(repr=False)

    @property
    def sigma_slope_scaled(self):
        return self.sigma_slope * self.scale


def weighted_line(x, y, sy, sx=None, iterations=5) -> LineFit:
    x, y, sy = map(np.asarray, (x, y, sy))
    sx = np.zeros_like(x) if sx is None else np.broadcast_to(np.asarray(sx, float), x.shape)
    slope = np.polyfit(x, y, 1)[0] if len(x) > 1 else 0.0
    for _ in range(iterations):
        s_eff = np.sqrt(sy**2 + (slope * sx) ** 2)
        w = 1 / s_eff**2
        A = np.column_stack([x, np.ones_like(x)])
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        slope, intercept = cov @ (A.T @ (w * y))
    resid = y - (slope * x + intercept)
    chi2 = float(np.sum(resid**2 * w))
    dof = len(x) - 2
    scale = max(1.0, math.sqrt(chi2 / dof)) if dof > 0 else 1.0
    return LineFit(float(slope), float(intercept), math.sqrt(cov[0, 0]), math.sqrt(cov[1, 1]),
                   chi2, dof, scale, resid)


@dataclass
class SpatialFit:
    columns: list
    x: np.ndarray            # m
    lnB: np.ndarray
    sigma_lnB: np.ndarray
    phase: np.ndarray        # unwrapped, relative to the nearest-base sensor (rad)
    sigma_phase: np.ndarray
    amp_fit: LineFit
    phase_fit: LineFit
    q: float
    sigma_q: float
    qp: float
    sigma_qp: float


def spatial_fit(columns, x_m, fits, sigma_x_m):
    order = np.argsort(x_m)
    cols = [columns[i] for i in order]
    x = np.asarray(x_m, float)[order]
    f = [fits[i] for i in order]
    B = np.array([g.amp for g in f])
    lnB = np.log(B)
    s_lnB = np.array([g.sigma_amp / g.amp for g in f])
    ph = np.unwrap([g.phase for g in f])
    ph = ph - ph[0]
    s_ph = np.array([g.sigma_phase for g in f])
    af = weighted_line(x, lnB, s_lnB, sigma_x_m)
    pf = weighted_line(x, ph, s_ph, sigma_x_m)
    return SpatialFit(cols, x, lnB, s_lnB, ph, s_ph, af, pf,
                      -af.slope, af.sigma_slope_scaled, -pf.slope, pf.sigma_slope_scaled)


@dataclass
class FiniteSpatialFit:
    q: float
    qp: float
    cov: np.ndarray          # 2x2 covariance of (q, q'), scaled by max(1, chi2/dof)
    chi2: float
    dof: int
    lnB0: float
    phi0: float

    @property
    def sigma_q(self):
        return math.sqrt(self.cov[0, 0])

    @property
    def sigma_qp(self):
        return math.sqrt(self.cov[1, 1])


def finite_spatial_fit(sp: SpatialFit, L=ROD_ANALYSIS_LENGTH_M, sigma_x_m=0.0):
    """Fit the same amplitudes and unwrapped phases to the finite insulated-tip rod,
    ln B_i = ln B_0 + ln|A(x_i)|, phi_i = phi_0 + arg A(x_i), A = cosh[s(L-x)]/cosh(sL),
    s = q + i q' (four parameters). Compared with the semi-infinite straight-line fits,
    this shows how much the finite tip moves q, q', and especially nu."""
    from scipy.optimize import least_squares
    x = sp.x
    s_ln = np.sqrt(sp.sigma_lnB**2 + (sp.q * sigma_x_m) ** 2)
    s_ph = np.sqrt(sp.sigma_phase**2 + (sp.qp * sigma_x_m) ** 2)

    def model(p):
        lnB0, phi0, q, qp = p
        A = finite_amplitude(x, q, qp, L)
        return lnB0 + np.log(np.abs(A)), phi0 + np.unwrap(np.angle(A))

    def resid(p):
        mlnB, mph = model(p)
        return np.concatenate([(sp.lnB - mlnB) / s_ln, (sp.phase - mph) / s_ph])

    p0 = [sp.amp_fit.intercept, sp.phase_fit.intercept, sp.q, sp.qp]
    res = least_squares(resid, p0, x_scale="jac")
    J = res.jac
    chi2 = float(np.sum(res.fun**2))
    dof = 2 * len(x) - 4
    cov = np.linalg.pinv(J.T @ J)
    if dof > 0:
        cov = cov * max(1.0, chi2 / dof)
    return FiniteSpatialFit(float(res.x[2]), float(res.x[3]), cov[2:, 2:], chi2, dof,
                            float(res.x[0]), float(res.x[1]))


def kappa_nu_cov(omega, q, qp, cov):
    """kappa and nu with uncertainties from a full (q, q') covariance matrix."""
    kappa = omega / (2 * q * qp)
    nu = kappa * (q**2 - qp**2)
    gk = np.array([-kappa / q, -kappa / qp])
    gn = 0.5 * omega * np.array([1 / qp + qp / q**2, -(q / qp**2 + 1 / q)])
    return kappa, math.sqrt(gk @ cov @ gk), nu, math.sqrt(gn @ cov @ gn)


def kappa_nu(omega, q, sq, qp, sqp):
    kappa = omega / (2 * q * qp)
    s_kappa = kappa * math.hypot(sq / q, sqp / qp)
    nu = kappa * (q**2 - qp**2)
    dnu_dq = 0.5 * omega * (1 / qp + qp / q**2)
    dnu_dqp = -0.5 * omega * (q / qp**2 + 1 / q)
    s_nu = math.hypot(dnu_dq * sq, dnu_dqp * sqp)
    return kappa, s_kappa, nu, s_nu


def q_qprime(kappa, nu, omega):
    r = np.hypot(nu, omega)
    return np.sqrt((nu + r) / (2 * kappa)), np.sqrt((r - nu) / (2 * kappa))


def material(kappa, s_kappa, nu, s_nu, d=ROD_DIAMETER_M):
    """k = kappa rho c and h = nu rho c d / 4, assuming the course rho and c."""
    k = kappa * RHO * C_P
    h = nu * RHO * C_P * d / 4
    return k, s_kappa * RHO * C_P, h, s_nu * RHO * C_P * d / 4


def weighted_mean(values, sigmas):
    v, s = np.asarray(values, float), np.asarray(sigmas, float)
    w = 1 / s**2
    mean = float(np.sum(w * v) / np.sum(w))
    sigma = float(1 / math.sqrt(np.sum(w)))
    chi2 = float(np.sum(w * (v - mean) ** 2))
    dof = len(v) - 1
    p = float(chi2_dist.sf(chi2, dof)) if dof > 0 else float("nan")
    return mean, sigma, chi2, dof, p


# ------------------------------------------------------------------
# MODEL CHECKS AND PREDICTIONS
# ------------------------------------------------------------------
def finite_amplitude(x, q, qp, L=ROD_ANALYSIS_LENGTH_M):
    """Complex amplitude relative to x = 0 for the finite insulated-tip rod."""
    s = q + 1j * qp
    return np.cosh(s * (L - np.asarray(x))) / np.cosh(s * L)


def finite_correction(x, q, qp, L=ROD_ANALYSIS_LENGTH_M):
    """Relative amplitude change |A_fin/A_semi| - 1 and phase-lag change (rad) of the
    finite rod relative to the semi-infinite e^{-(q + i q') x}, and the Stage 6 bound
    on the complex correction |H_fin/H_semi - 1| <=
    (e^{-2q(L-x)} + e^{-2qL}) / (1 - e^{-2qL})."""
    x = np.asarray(x, float)
    s = q + 1j * qp
    R = finite_amplitude(x, q, qp, L) / np.exp(-s * x)
    bound = (np.exp(-2 * q * (L - x)) + np.exp(-2 * q * L)) / (1 - np.exp(-2 * q * L))
    return np.abs(R) - 1, -np.angle(R), bound


def finite_complex_correction(x, q, qp, L=ROD_ANALYSIS_LENGTH_M):
    """|H_fin/H_semi - 1|, the quantity the Stage 6 bound limits."""
    s = q + 1j * qp
    x = np.asarray(x, float)
    return np.abs(finite_amplitude(x, q, qp, L) / np.exp(-s * x) - 1)


def biot_numbers(kappa, nu, qp, d=ROD_DIAMETER_M):
    """Bi_perp = h d/(4k) = nu d^2/(16 kappa) (Module 7 transverse check) and q' R,
    the radius in units of the one-radian phase length (radial uniformity of the
    oscillation)."""
    return nu * d**2 / (16 * kappa), qp * d / 2


def finite_mean(x, theta0, m, L=ROD_ANALYSIS_LENGTH_M):
    return theta0 * np.cosh(m * (L - np.asarray(x))) / np.cosh(m * L)


@dataclass
class MeanFit:
    theta0: float
    sigma_theta0: float
    m: float                # mu / L
    sigma_m: float          # from sigma_mu and sigma_L
    chi2: float
    dof: int
    resid: np.ndarray = field(repr=False)
    mu: float = float("nan")            # mL, the fitted dimensionless parameter
    sigma_mu: float = float("nan")
    sigma_m_fit: float = float("nan")   # from sigma_mu alone
    L: float = ROD_ANALYSIS_LENGTH_M
    sigma_L: float = 0.0


def mean_profile_fit(x_m, theta, sigma, L=ROD_ANALYSIS_LENGTH_M, m0=4.9, sigma_L=0.0) -> MeanFit:
    """Stage 6 mean profile: fit (xi_i, theta_i), xi = x/L, to
    theta0 cosh[mu (1 - xi)] / cosh(mu); then m = mu/L with
    sigma_m^2 = (sigma_mu/L)^2 + (mu sigma_L / L^2)^2."""
    x_m, theta, sigma = map(np.asarray, (x_m, theta, sigma))
    xi = x_m / L
    popt, pcov = curve_fit(lambda z, a, mu: a * np.cosh(mu * (1 - z)) / np.cosh(mu), xi, theta,
                           p0=(theta[np.argmin(x_m)], m0 * L), sigma=sigma, absolute_sigma=True)
    theta0, mu = popt
    s_theta0, s_mu = math.sqrt(pcov[0, 0]), math.sqrt(pcov[1, 1])
    resid = theta - finite_mean(x_m, theta0, mu / L, L)
    chi2 = float(np.sum((resid / sigma) ** 2))
    m = mu / L
    s_m = math.hypot(s_mu / L, mu * sigma_L / L**2)
    return MeanFit(theta0, s_theta0, m, s_m, chi2, len(x_m) - 2, resid, mu, s_mu, s_mu / L,
                   L, sigma_L)


def predict_rod(x_m, t, t_inf, base_fit: SineFit, kappa, nu, omega, L=ROD_ANALYSIS_LENGTH_M):
    """Finite-rod prediction at position x from the measured T1 fit (mean and
    oscillation) and kappa, nu: T = T_inf + (a1 - T_inf) cosh[m(L-x)]/cosh(mL)
    + B1 |A(x)| cos(wt + phi1 + arg A(x))."""
    q, qp = q_qprime(kappa, nu, omega)
    m = math.sqrt(nu / kappa)
    A = finite_amplitude(x_m, q, qp, L)
    mean = t_inf + finite_mean(x_m, base_fit.mean - t_inf, m, L)
    return mean + base_fit.amp * abs(A) * np.cos(omega * t + base_fit.phase + np.angle(A))


def length_curves(periods, kappa, s_kappa, nu, s_nu, n_mc=4000, seed=1):
    """Central 1/m, 1/q, 1/q' (m) vs period and 16-84 % bands from sampling kappa and nu."""
    periods = np.asarray(periods, float)
    omega = 2 * np.pi / periods
    q, qp = q_qprime(kappa, nu, omega)
    central = {"mean": np.full_like(periods, math.sqrt(kappa / nu)), "amp": 1 / q, "phase": 1 / qp}
    rng = np.random.default_rng(seed)
    k = rng.normal(kappa, s_kappa, n_mc)
    n = rng.normal(nu, s_nu, n_mc)
    ok = (k > 0) & (n > 0)
    k, n = k[ok, None], n[ok, None]
    qs, qps = q_qprime(k, n, omega[None, :])
    samples = {"mean": np.sqrt(k / n) * np.ones_like(qs), "amp": 1 / qs, "phase": 1 / qps}
    bands = {key: np.percentile(v, [16, 84], axis=0) for key, v in samples.items()}
    return central, bands


def positions_for(meta, override_cm=None):
    """Sensor positions (cm) from --positions-cm, else the metadata file, else the
    module_08_common defaults."""
    if override_cm is not None:
        return dict(zip(ROD_COLUMNS, override_cm))
    pos = dict(SENSOR_POSITIONS_CM)
    pos.update(meta.get("positions_cm", {}))
    return pos


EXTRAPOLATION_FLOOR_C = 0.1


def baseline_temperatures(spec: str):
    """Per-channel room temperature T_inf (C) and its uncertainty, by channel.
      FILE:START,END              mean over that Arduino-time interval (TEC off, rod at room)
      FILE:START,END:extrapolate  the rod was still relaxing: fit T = T_inf + A exp(-t/tau)
                                  per channel and use the asymptote (uncertainty: fit error
                                  combined with EXTRAPOLATION_FLOOR_C)"""
    extrapolate = spec.endswith(":extrapolate")
    if extrapolate:
        spec = spec[: -len(":extrapolate")]
    name, _, rng = spec.rpartition(":")
    start, end = (float(v) for v in rng.split(","))
    path = resolve_path(name)
    parts = load_parts(path)
    for p in parts:
        m = (p["time_s"] >= start) & (p["time_s"] <= end)
        if m.sum() >= 10:
            if not extrapolate:
                return {c: (float(np.nanmean(p[c][m])),
                            float(np.nanstd(p[c][m]) / math.sqrt(m.sum())))
                        for c in TEMPERATURE_COLUMNS}
            out = {}
            t = p["time_s"][m] - start
            for c in TEMPERATURE_COLUMNS:
                y = p[c][m]
                popt, pcov = curve_fit(lambda t, Tr, A, tau: Tr + A * np.exp(-t / tau), t, y,
                                       p0=(y[-1] - 0.3, y[0] - y[-1] + 0.3, 300.0), maxfev=20000)
                out[c] = (float(popt[0]), math.hypot(math.sqrt(pcov[0, 0]), EXTRAPOLATION_FLOOR_C))
            return out
    raise ValueError(f"baseline {spec}: fewer than 10 records in that interval")
