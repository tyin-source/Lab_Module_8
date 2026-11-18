# Phys 39 Module 8 - automatic run schedule for the Angstrom GUI
#
# A schedule is a CSV file with one row per segment (see schedules/README.md). The GUI
# (python/module_08_angstrom_gui.py, "Automatic schedule") applies each row's settings
# through its normal controls, starts the run if needed, and calls
# ScheduleRunner.on_record() with every accepted DATA record. The runner decides when
# the segment is complete:
#   OFF / MANUAL rows: after min_s seconds.
#   PWM_SINE / FF / P / PI rows: when the settling test (angstrom_lib.settling_test,
#     one-cycle windows of the settle_on channels) finds a settled start and at least
#     cycles_after_settle whole periods have been recorded after it; or, if that has
#     not happened by max_s, the segment ends anyway and is marked NOT SETTLED.
#   TEST_T1_LIMIT / TEST_TEC_LIMIT rows (automatic S18.1 safety check): the GUI sets
#     the T1 (or TEC) limit in force to the present reading + limit_rise_C, heats with
#     MANUAL manual_pwm (default 40), and the row passes when the Arduino stops the
#     run with run_state SAFETY and that channel at the lowered limit. No shutdown
#     within max_s, or a shutdown by the other channel, aborts the schedule. The GUI
#     restores 56 C and 60 C before the next row and whenever a schedule ends.
#   max_spread_C on an OFF row: the six channel means must agree within this
#     (Module 7 criterion 1 C), else the schedule stops.
#   kind "check": the segment is not written to an interval file.
# Glitches (angstrom_lib.glitch_at: several channels jumping together, or one
# downstream channel jumping, as a loose breadboard contact does) are logged and
# restart the current segment's collection and settling test, so no interval written
# to the interval files contains one. Limit-test rows are not restarted.
# A segment's data are the records with the Arduino step of its last setting change;
# if the step changes during the segment (someone changed a setting), collection
# restarts. Each completed segment is appended to data/module_08/intervals_auto_<kind>.csv
# (the format the analysis scripts read), and a period change during a run keeps the
# drive phase continuous because the Arduino only receives SET PERIOD.
# The runner aborts on run_state SAFETY or TIMEOUT once the row has reached its run
# state (a SAFETY or TIMEOUT left over from before the row was applied is ignored
# while the row's START takes effect), or if the Arduino does not reach the expected
# run state within 15 s of a row being applied (for example a refused START). It never
# restarts a run by itself after an abort.

import csv
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
from module_08_common import REPO_ROOT, TEMPERATURE_COLUMNS  # noqa: E402

TEST_MODES = {"TEST_T1_LIMIT": ("T1_C", "T1LIMIT", 56.0),
              "TEST_TEC_LIMIT": ("TEC_C", "TECLIMIT", 60.0)}   # channel, setting, normal limit
MODES = ("OFF", "MANUAL", "PWM_SINE", "FF", "P", "PI", *TEST_MODES)
KINDS = ("check", "baseline", "control", "pilot", "production")
TEST_PWM_DEFAULT, TEST_PWM_MAX = 40.0, 80.0
TEST_RISE_DEFAULT = 2.0
# schedule column -> GUI setting (SET command name)
SETTING_COLUMNS = {"Ta_C": "MEAN", "Tb_C": "AMP", "period_s": "PERIOD", "Kp": "KP",
                   "Ki": "KI", "u_a": "PWMOFF", "u_b": "PWMAMP", "manual_pwm": "MANUAL",
                   "T_room_C": "TROOM"}
SCHEDULE_COLUMNS = (["label", "kind", "mode"] + list(SETTING_COLUMNS)
                    + ["settle_on", "cycles_after_settle", "min_s", "max_s", "limit_rise_C",
                       "max_spread_C", "notes"])
INTERVAL_COLUMNS = ["label", "file", "part", "step", "t_start_s", "t_end_s", "settle",
                    "exclude", "notes"]
START_TIMEOUT_S = 15.0
SETTLE_CHECK_EVERY = 5          # records between settling tests

# limits checked when a schedule is loaded (the sketch checks them again)
MIN_PERIOD_S = 100.0
SETPOINT_MIN_C, SETPOINT_MAX_C = 10.0, 56.0


@dataclass
class ScheduleRow:
    label: str
    kind: str
    mode: str
    settings: dict                 # SET name -> float (only the columns given)
    settle_on: list = field(default_factory=lambda: ["T1_C"])
    cycles_after_settle: int = 5
    min_s: float = 0.0
    max_s: float = None
    notes: str = ""
    limit_rise: float = TEST_RISE_DEFAULT
    max_spread: float = None

    @property
    def is_test(self):
        return self.mode in TEST_MODES

    def describe(self):
        parts = [f"{self.label} [{self.kind}] {self.mode}"]
        parts += [f"{k}={v:g}" for k, v in self.settings.items()]
        if self.is_test:
            ch = TEST_MODES[self.mode][0][:-2]
            parts.append(f"limit = {ch} now + {self.limit_rise:g} C, heat until SAFETY "
                         f"(max {self.max_s:g} s)")
        elif self.mode in ("OFF", "MANUAL"):
            parts.append(f"for {self.min_s:g} s")
            if self.max_spread is not None:
                parts.append(f"channel spread < {self.max_spread:g} C")
        else:
            parts.append(f"settled on {'/'.join(c[:-2] for c in self.settle_on)} + "
                         f"{self.cycles_after_settle} cycles (min {self.min_s:g} s, max {self.max_s:g} s)")
        return ", ".join(parts)


def _float(row, key, label):
    text = (row.get(key) or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"{label}: {key} = {text!r} is not a number")


def read_schedule(path: Path):
    """Parse and check a schedule CSV. Rows whose label starts with # are ignored."""
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        missing = [c for c in ("label", "kind", "mode") if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path.name}: missing columns {missing}")
        for raw in reader:
            label = (raw.get("label") or "").strip()
            if not label or label.startswith("#"):
                continue
            kind = (raw.get("kind") or "").strip().lower()
            mode = (raw.get("mode") or "").strip().upper()
            if kind not in KINDS:
                raise ValueError(f"{label}: kind must be one of {KINDS}")
            if mode not in MODES:
                raise ValueError(f"{label}: mode must be one of {MODES}")
            settings = {}
            for col, name in SETTING_COLUMNS.items():
                v = _float(raw, col, label)
                if v is not None:
                    settings[name] = v
            settle = [s.strip().upper() for s in (raw.get("settle_on") or "T1").replace(";", " ").split()]
            settle_cols = []
            for s in settle:
                col = s if s.endswith("_C") else f"{s}_C"
                if col not in TEMPERATURE_COLUMNS:
                    raise ValueError(f"{label}: settle_on {s!r} is not TEC or T1-T5")
                settle_cols.append(col)
            cycles = _float(raw, "cycles_after_settle", label)
            min_s = _float(raw, "min_s", label) or 0.0
            max_s = _float(raw, "max_s", label)
            rise = _float(raw, "limit_rise_C", label)
            spread = _float(raw, "max_spread_C", label)
            r = ScheduleRow(label, kind, mode, settings, settle_cols,
                            int(cycles) if cycles is not None else 5, min_s, max_s,
                            (raw.get("notes") or "").strip(),
                            rise if rise is not None else TEST_RISE_DEFAULT, spread)
            if r.is_test:
                r.settings.setdefault("MANUAL", TEST_PWM_DEFAULT)
                if not 0 < r.settings["MANUAL"] <= TEST_PWM_MAX:
                    raise ValueError(f"{label}: manual_pwm for a limit test must be 1-{TEST_PWM_MAX:g}")
                if not 0.5 <= r.limit_rise <= 5:
                    raise ValueError(f"{label}: limit_rise_C must be 0.5-5 C")
                if r.max_s is None:
                    r.max_s = 600.0
            elif mode in ("OFF", "MANUAL"):
                if min_s <= 0:
                    raise ValueError(f"{label}: {mode} rows need min_s (the duration)")
            elif max_s is None or max_s <= min_s:
                raise ValueError(f"{label}: give max_s larger than min_s")
            if r.cycles_after_settle < 1:
                raise ValueError(f"{label}: cycles_after_settle must be at least 1")
            rows.append(r)
    if not rows:
        raise ValueError(f"{path.name}: no schedule rows")
    check_limits(rows)
    return rows


def check_limits(rows):
    """Follow the settings row by row (blank = unchanged) and check the course limits
    where the values are known."""
    known = {}
    for r in rows:
        known.update(r.settings)
        if "PERIOD" in r.settings and r.settings["PERIOD"] < MIN_PERIOD_S:
            raise ValueError(f"{r.label}: period below {MIN_PERIOD_S:g} s")
        if r.mode in ("FF", "P", "PI") and "MEAN" in known and "AMP" in known:
            lo, hi = known["MEAN"] - known["AMP"], known["MEAN"] + known["AMP"]
            if lo < SETPOINT_MIN_C or hi > SETPOINT_MAX_C:
                raise ValueError(f"{r.label}: setpoint {lo:g}-{hi:g} C outside "
                                 f"{SETPOINT_MIN_C:g}-{SETPOINT_MAX_C:g} C")
        if r.mode == "PWM_SINE":
            ua, ub = known.get("PWMOFF", 0.0), known.get("PWMAMP", 0.0)
            if abs(ua) + ub > 255:
                raise ValueError(f"{r.label}: |u_a| + u_b above 255")


@dataclass
class Action:
    kind: str                     # "apply", "status", "note", "finished", "abort"
    text: str = ""
    row: ScheduleRow = None


class ScheduleRunner:
    def __init__(self, rows, csv_path_fn, intervals_dir: Path):
        self.rows = rows
        self.csv_path_fn = csv_path_fn          # returns the current raw CSV path
        self.intervals_dir = intervals_dir
        self.index = -1
        self.active = False
        self.results = []                       # one dict per completed row

    # ---- control ----
    def start(self):
        self.active = True
        self.index = 0
        return Action("apply", f"row 1/{len(self.rows)}: {self.rows[0].describe()}", self.rows[0])

    def abort(self):
        self.active = False

    def row_applied(self, t_now, test_limit=None):
        """Called by the GUI right after it sent a row's commands. test_limit: the
        lowered limit (C) a TEST_* row set."""
        self.applied_at = t_now
        self.test_limit = test_limit
        self.reached = False                    # the row's run state has been seen
        self.row_t0 = None                      # first record in the row's run state
        self._reset_segment()

    def _reset_segment(self):
        self.prev = []                          # last two records' temperatures
        self.seg_step = None
        self.seg_t0 = None
        self.data = {c: [] for c in ["time_s", *TEMPERATURE_COLUMNS]}
        self.period = None
        self.n_since_check = 0
        self.settled_start = None

    @property
    def row(self):
        return self.rows[self.index]

    # ---- one record ----
    def on_record(self, rec):
        if not self.active:
            return None
        t = rec["time_ms"] / 1000.0
        state = rec["run_state"]
        r = self.row
        if r.is_test and state == "SAFETY" and self.reached:
            ch = TEST_MODES[r.mode][0]
            v = rec[ch]
            if self.test_limit is not None and v >= self.test_limit - 0.15:
                self._append(rec, t)
                return self._finish(True, self.seg_t0 if self.seg_t0 is not None else t, t,
                                    f"shutdown with {ch[:-2]} = {v:.2f} C at the lowered "
                                    f"limit {self.test_limit:.2f} C: PASS")
            self.active = False
            return Action("abort", f"{r.label}: SAFETY, but {ch[:-2]} = {v:.2f} C is below the "
                                   f"test limit {self.test_limit:.2f} C (another limit stopped "
                                   "the run?)")
        if state in ("SAFETY", "TIMEOUT") and self.reached:
            self.active = False
            return Action("abort", f"run_state {state} during {r.label}: schedule stopped")
        want = "OFF" if r.mode == "OFF" else "RUN"
        want_mode = "MANUAL" if r.is_test else r.mode
        if state != want or (want == "RUN" and rec["control_mode"] != want_mode):
            if self.applied_at is not None and t - self.applied_at > START_TIMEOUT_S:
                self.active = False
                return Action("abort", f"{r.label}: the Arduino is {state} "
                                       f"{rec['control_mode']}, not {want} {want_mode}, "
                                       f"{START_TIMEOUT_S:g} s after the row was applied "
                                       "(START refused or a setting rejected?)")
            return Action("status", f"{r.label}: waiting for {want} {want_mode}")
        self.reached = True
        if self.row_t0 is None:
            self.row_t0 = t
        row_elapsed = t - self.row_t0

        if self.seg_step != rec["step"]:
            restarted = self.seg_step is not None
            self._reset_segment()
            self.seg_step, self.seg_t0 = rec["step"], t
            if restarted:
                self._append(rec, t)
                return Action("status", f"{r.label}: step changed to {rec['step']}; "
                                        "segment restarted")
        glitch = self._glitch(rec)
        if glitch and not r.is_test:
            self.glitches = getattr(self, "glitches", 0) + 1
            self._reset_segment()
            self.seg_step, self.seg_t0 = rec["step"], t
            self._append(rec, t)
            self.prev = [{c: rec[c] for c in TEMPERATURE_COLUMNS}]
            return Action("note", f"{r.label}: glitch at {t:.0f} s ({glitch}); segment restarted")
        self._append(rec, t)
        elapsed = t - self.seg_t0

        if r.is_test:
            ch = TEST_MODES[r.mode][0]
            if row_elapsed >= r.max_s:
                self.active = False
                return Action("abort", f"{r.label}: no shutdown within {r.max_s:g} s "
                                       f"({ch[:-2]} = {rec[ch]:.2f} C, limit {self.test_limit:.2f} C): "
                                       "limit test FAILED")
            return Action("status", f"{r.label}: heating, {ch[:-2]} = {rec[ch]:.2f} C, "
                                    f"lowered limit {self.test_limit:.2f} C")

        if r.mode in ("OFF", "MANUAL"):
            if elapsed >= r.min_s:
                return self._finish(True, self.seg_t0, t)
            if row_elapsed > max(3 * r.min_s, r.min_s + 600):
                self.active = False
                return Action("abort", f"{r.label}: glitches kept restarting the segment for "
                                       f"{row_elapsed:.0f} s; fix the breadboard contacts")
            return Action("status", f"{r.label}: {elapsed:.0f} / {r.min_s:g} s")

        P = rec["period_s"]
        self.period = P
        self.n_since_check += 1
        if elapsed >= max(r.min_s, 3 * P) and self.n_since_check >= SETTLE_CHECK_EVERY:
            self.n_since_check = 0
            self.settled_start = self._settled_start(P)
        if self.settled_start is not None:
            done = (t - self.settled_start) / P
            if done >= r.cycles_after_settle + 0.02:
                return self._finish(True, self.settled_start, t)
            text = (f"{r.label}: settled from {self.settled_start - self.seg_t0:.0f} s; "
                    f"{done:.1f} / {r.cycles_after_settle} cycles after settling")
        else:
            text = f"{r.label}: {elapsed:.0f} s, not settled yet ({row_elapsed:.0f} / max {r.max_s:g} s)"
        if row_elapsed >= r.max_s:
            return self._finish(False, max(self.seg_t0, t - r.cycles_after_settle * P), t)
        return Action("status", text)

    # ---- helpers ----
    def _glitch(self, rec):
        import angstrom_lib as al
        now = {c: rec[c] for c in TEMPERATURE_COLUMNS}
        g = al.glitch_at(self.prev[0], self.prev[1], now) if len(self.prev) == 2 else None
        self.prev = (self.prev + [now])[-2:]
        return g

    def _append(self, rec, t):
        self.data["time_s"].append(t)
        for c in TEMPERATURE_COLUMNS:
            self.data[c].append(rec[c])

    def _settled_start(self, period):
        import angstrom_lib as al
        d = {k: np.asarray(v, float) for k, v in self.data.items()}
        try:
            res = al.settling_test(d, self.row.settle_on, period)
        except ValueError:
            return None
        return res.settled_start

    def channel_spread(self):
        """Largest minus smallest channel mean over the current segment (C)."""
        means = [np.nanmean(self.data[c]) for c in TEMPERATURE_COLUMNS if self.data[c]]
        return float(np.max(means) - np.min(means)) if means else math.nan

    def _finish(self, settled, t_start, t_end, extra=""):
        r = self.row
        if r.mode in ("OFF", "MANUAL"):
            spread = self.channel_spread()
            extra = f"channel spread {spread:.2f} C"
            if r.max_spread is not None and not spread <= r.max_spread:
                self.active = False
                return Action("abort", f"{r.label}: channel spread {spread:.2f} C exceeds "
                                       f"{r.max_spread:g} C (check the thermistors and wiring)")
        self._write_interval(r, settled, t_start, t_end)
        self.results.append({"label": r.label, "kind": r.kind, "mode": r.mode,
                             "settled": settled, "t_start": t_start, "t_end": t_end,
                             "step": self.seg_step, "file": str(self.csv_path_fn())})
        text = f"{r.label} done: step {self.seg_step}, {t_start:.0f}-{t_end:.0f} s"
        if r.mode in ("PWM_SINE", "FF", "P", "PI"):
            text += ", settled" if settled else ", NOT SETTLED (max_s reached)"
        if extra:
            text += f", {extra}"
        self.index += 1
        if self.index >= len(self.rows):
            self.active = False
            return Action("finished", text + "; schedule complete")
        nxt = self.rows[self.index]
        return Action("apply", text + f"; row {self.index + 1}/{len(self.rows)}: {nxt.describe()}", nxt)

    def baseline_mean(self):
        """Mean of the six channels over the segment just finished (room temperature)."""
        vals = [np.nanmean(self.data[c]) for c in TEMPERATURE_COLUMNS if self.data[c]]
        return float(np.mean(vals)) if vals else math.nan

    def _write_interval(self, r, settled, t_start, t_end):
        if r.kind == "check":
            return
        path = self.intervals_dir / f"intervals_auto_{r.kind}.csv"
        new = not path.exists()
        csv_path = Path(self.csv_path_fn())
        try:
            name = str(csv_path.resolve().relative_to(REPO_ROOT))
        except ValueError:
            name = str(csv_path)
        notes = f"auto schedule; {r.mode}"
        if r.mode not in ("OFF", "MANUAL"):
            notes += "; settled" if settled else "; NOT SETTLED (max_s reached)"
        if r.notes:
            notes += f"; {r.notes}"
        with open(path, "a", newline="") as f:
            w = csv.writer(f)
            if new:
                w.writerow(INTERVAL_COLUMNS)
            w.writerow([r.label, name, "", self.seg_step, f"{t_start:.1f}", f"{t_end:.1f}",
                        "no", "", notes])
