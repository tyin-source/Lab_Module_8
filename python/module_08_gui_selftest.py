# Phys 39 Module 8: self-test of the Angstrom GUI without hardware
#
# Runs module_08_angstrom_gui.py off screen against a simulated Arduino that speaks
# the Module 8 protocol (DATA, EVENT, CONFIG, ERROR lines; HB, START, STOP, MODE,
# SET, ZERO, GET commands) and drives a one-lump TEC and a five-node rod. The
# control law itself is tested against the real sketch in
# tests/firmware_host/firmware_host_test.cpp; here the simulator only needs to
# answer like the Arduino. Checks:
#   1. start-up text and the header are ignored, malformed DATA lines are counted
#      as rejected, and the program keeps running
#   2. on the first record the GUI sends every setting, then GET; after every
#      record it sends HB
#   3. "Send settings" sends only the changed values (a period change during a run
#      is a single SET PERIOD), and Start sends START
#   4. the raw CSV has exactly the 18 DATA columns and one row per accepted
#      record, identical to the fields the Arduino sent (nan included)
#   5. EVENT, CONFIG, ERROR lines and notes reach the metadata file
#   6. an ERROR reply makes the next "Send settings" resend everything
#   7. a SAFETY record shows the red banner and the GUI does not restart the run
#   8. an Arduino restart (time goes back) is logged and the settings are resent
#   9. "New data file" is refused during a run, works when stopped, and a file
#      created in the same second gets a new name instead of overwriting
#  10. closing the GUI sends STOP
#  11. automatic schedule: a baseline, PI at 100 s, then 200 s; each segment ends after
#      the live settling test plus the requested cycles; T_room is set from the
#      baseline; the period change is a single SET PERIOD; the segments are written to
#      intervals_auto_*.csv and module_08_control_fit.py reads them; STOP at the end
#  12. a SAFETY record aborts a running schedule and sends STOP
#  13. command pacing: with COMMAND_SPACING_MS > 0 settings and START are queued and
#      sent one at a time in order; STOP goes out at once and discards the queue
#  14. an echoed SET value different from the one sent, or a SET with no echo within
#      ECHO_TIMEOUT_S, is treated as an ERROR
#  15. a schedule started while the Arduino is still in SAFETY waits for its START
#      instead of aborting
#  16. automatic S18.1 checks: channel spread within 1 C, TEST_TEC_LIMIT and
#      TEST_T1_LIMIT lower the limit to the reading + 2 C, heat, pass on SAFETY, and
#      the normal 56/60 C limits are restored; a channel spread above 1 C stops the run
#  17. autorun: the schedule starts by itself once records arrive, and the control
#      analysis runs afterwards on this run's segments
#  18. the serial thread answers each DATA line with HB while the window has handled a
#      record within HEARTBEAT_WATCHDOG_S, and stops when the window has stalled longer
# Data files go to a temporary folder, not data/module_08/.
#
# Usage (from the repository root):
#   python python/module_08_gui_selftest.py

import contextlib
import csv
import io
import math
import os
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6 import QtWidgets
import pyqtgraph as pg

import module_08_angstrom_gui as gui_module
from module_08_common import DATA_COLUMNS

T_AMB = 22.0
TAU_TEC = 90.0
CHI_H, CHI_C = 0.40, 0.14


class SimulatedArduino:
    """Fake serial link and simulated Module 8 sketch + TEC + rod."""

    def __init__(self):
        self.sent = []                 # every command the GUI wrote
        self.out = []                  # lines waiting to be read by the test
        self.reset()

    def reset(self):
        self.t_ms = 0
        self.state, self.mode, self.step = "OFF", "MANUAL", 0
        self.cfg = {"MEAN": 40.0, "AMP": 16.0, "PERIOD": 100.0, "KP": 6.0, "KI": 0.04,
                    "UMAX": 255.0, "RATIO": 2.82, "CHIH": 0.485, "TROOM": 22.0,
                    "PWMAMP": 50.0, "PWMOFF": 0.0, "MANUAL": 0.0, "FFADD": 0.0,
                    "T1LIMIT": 56.0, "TECLIMIT": 60.0}
        self.phase, self.q = 0.0, 0.0
        self.temps = [T_AMB] * 6
        self.force_unsafe = False
        self.out += ["Module 8 Angstrom rod controller (simulated)",
                     "DATA," + ",".join(DATA_COLUMNS)]

    # ---- serial link interface used by the GUI ----
    def readline(self):
        time.sleep(0.05)
        return b""

    def close(self):
        pass

    def event(self, text):
        self.step += 1
        self.out.append(f"EVENT,{self.t_ms},{self.step},{text}")

    def write_line(self, text):
        self.sent.append(text)
        words = text.strip().upper().split()
        if words == ["HB"]:
            return
        if words == ["GET"]:
            self.out.append(f"CONFIG,{self.t_ms},{self.step},mode={self.mode} "
                            + " ".join(f"{k}={v:g}" for k, v in self.cfg.items()))
        elif words == ["START"]:
            if (self.force_unsafe or self.temps[1] > self.cfg["T1LIMIT"]
                    or self.temps[0] > self.cfg["TECLIMIT"]):
                self.out.append("ERROR: START refused: T1 above the rod-base limit")
            else:
                self.state, self.phase, self.q = "RUN", 0.0, 0.0
                self.event(f"START mode={self.mode}")
        elif words == ["STOP"]:
            self.state = "OFF"
            self.event("STOP (output 0)")
        elif words == ["ZERO"]:
            self.q = 0.0
            self.event("ZERO integral")
        elif len(words) == 2 and words[0] == "MODE":
            self.mode, self.q = words[1], 0.0
            self.event(f"MODE {words[1]}")
        elif len(words) == 3 and words[0] == "SET" and words[1] in self.cfg:
            v = float(words[2])
            mean = v if words[1] == "MEAN" else self.cfg["MEAN"]
            amp = v if words[1] == "AMP" else self.cfg["AMP"]
            if words[1] == "PERIOD" and v < 100:
                self.out.append(f"ERROR: period outside the allowed range: {text}")
            elif (words[1] == "T1LIMIT" and not 0 <= v <= 56) or (words[1] == "TECLIMIT"
                                                                 and not 0 <= v <= 60):
                self.out.append(f"ERROR: limit outside the allowed range: {text}")
            elif mean - amp < 10 or mean + amp > 56:
                self.out.append(f"ERROR: Ta - Tb or Ta + Tb outside the setpoint limits: {text}")
            else:
                self.cfg[words[1]] = v
                self.event(f"SET {words[1]} {v:.4f}")
        else:
            self.out.append(f"ERROR: unknown command: {text}")

    # ---- one ~1 s average, then one DATA record ----
    def record(self, state_override=None):
        c = self.cfg
        self.t_ms += 1000
        sp = ff = ctrl = ratio = math.nan
        u = 0.0
        if self.state == "RUN":
            if self.mode != "MANUAL":
                self.phase += 2 * math.pi / c["PERIOD"]
            s = math.sin(self.phase)
            if self.mode == "MANUAL":
                u = c["MANUAL"]
            elif self.mode == "PWM_SINE":
                u = c["PWMOFF"] + c["PWMAMP"] * s
            else:
                sp = c["MEAN"] + c["AMP"] * s
                ratio = c["RATIO"]
                uff = (sp - c["TROOM"]) / c["CHIH"]
                if self.mode == "FF":
                    ff, u = uff, uff
                else:
                    e = sp - self.temps[0]
                    if self.mode == "PI":
                        self.q += e
                    ctrl = c["KP"] * e + c["KI"] * self.q
                    u = ctrl
                if u < 0:
                    u *= c["RATIO"]
            u = max(-c["UMAX"], min(c["UMAX"], u))
        u = round(u)
        if self.state == "RUN" and (self.temps[1] > c["T1LIMIT"] or self.temps[0] > c["TECLIMIT"]):
            self.state, u = "SAFETY", 0
            self.event("SAFETY SHUTDOWN: a lowered limit was exceeded")
        for _ in range(10):
            chi = CHI_H if u >= 0 else CHI_C
            self.temps[0] += 0.1 / TAU_TEC * (T_AMB + chi * u - self.temps[0])
            prev = self.temps[0]
            for i in range(1, 6):
                new = self.temps[i] + 0.1 * ((prev - self.temps[i]) / (25 if i == 1 else 45)
                                             - (self.temps[i] - T_AMB) / 900)
                prev, self.temps[i] = self.temps[i], new

        def f(x, d):
            return "nan" if math.isnan(x) else f"{x:.{d}f}"
        uses_sp = self.mode in ("FF", "P", "PI")
        fields = [str(self.t_ms), state_override or self.state, self.mode, str(self.step),
                  f(c["MEAN"] if uses_sp else math.nan, 2), f(c["AMP"] if uses_sp else math.nan, 2),
                  f(c["PERIOD"] if self.mode != "MANUAL" else math.nan, 1),
                  f(sp, 2), f(ff, 1), *[f(t, 2) for t in self.temps],
                  f(ctrl, 1), f(ratio, 3), str(u)]
        self.out.append("DATA," + ",".join(fields))


def pump(app, gui, arduino):
    """Deliver every queued Arduino line to the GUI; return the DATA records."""
    records = []
    with contextlib.redirect_stdout(io.StringIO()):
        while arduino.out:
            line = arduino.out.pop(0)
            gui.on_line(line)
            if line.startswith("DATA,") and not line.startswith("DATA,time_ms"):
                records.append(line)
            app.processEvents()
    return records


def run_records(app, gui, arduino, n, **kw):
    records = []
    for _ in range(n):
        arduino.record(**kw)
        records += pump(app, gui, arduino)
    return records


def quiet(fn, *args):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args)


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'} {name}" + (f": {detail}" if detail else ""))
    return ok


def main():
    app = QtWidgets.QApplication(sys.argv)
    pg.setConfigOptions(background="w", foreground="k")
    results = []
    gui_module.COMMAND_SPACING_MS = 0   # send at once; pacing is tested in 13
    warnings = []                    # message boxes would block off screen
    QtWidgets.QMessageBox.warning = staticmethod(lambda *a, **k: warnings.append(a[1:]))
    with tempfile.TemporaryDirectory() as tmp:
        gui_module.DATA_DIR = Path(tmp)
        arduino = SimulatedArduino()
        gui = quiet(gui_module.AngstromGUI, arduino)
        first_csv = gui.files.csv_path

        # 1-2. start-up, first record, sync, heartbeat, malformed lines
        good = pump(app, gui, arduino)
        good += run_records(app, gui, arduino, 1)
        sync = [c for c in arduino.sent if c != "HB"]
        results.append(check("first record: every setting sent, then GET",
                             len(sync) == len(gui_module.SETTINGS) + 3 and "MODE MANUAL" in sync and sync[-1] == "GET",
                             f"{len(sync)} commands"))
        good += pump(app, gui, arduino)            # EVENT and CONFIG replies
        good += run_records(app, gui, arduino, 20)
        bad = ["DATA,99,OFF,MANUAL,0,nan",                                         # too few fields
               "DATA,99,WALK,MANUAL,0,nan,nan,nan,nan,nan,22,22,22,22,22,22,nan,nan,0",  # bad state
               "DATA,99,OFF,MANUAL,0,nan,nan,nan,nan,nan,abc,22,22,22,22,22,nan,nan,0",  # not a number
               "DATA,99,OFF,MANUAL,0,nan,nan,nan,nan,nan,22,22,22,22,22,22,nan,nan,999"]  # PWM range
        for line in bad:
            quiet(gui.on_line, line)
        results.append(check("accepted and rejected counts",
                             gui.n_accepted == 21 and gui.n_rejected == 4,
                             f"accepted {gui.n_accepted}, rejected {gui.n_rejected}"))
        results.append(check("HB after every record", arduino.sent.count("HB") == 21))

        # 3. send only changes, then start
        n0 = len(arduino.sent)
        gui.mode_combo.setCurrentText("PI")
        gui.spins["AMP"].setValue(4.0)
        quiet(gui.start_button.click)
        cmds = [c for c in arduino.sent[n0:] if c != "HB"]
        results.append(check("Start sends the two changes and START",
                             cmds == ["MODE PI", "SET AMP 4.0000", "START"], str(cmds)))
        good += run_records(app, gui, arduino, 150)
        results.append(check("RUN state shown", gui.lbl_state.text() == "State: RUN"))
        n0 = len(arduino.sent)
        gui.spins["PERIOD"].setValue(200)
        quiet(gui.send_button.click)
        cmds = [c for c in arduino.sent[n0:] if c != "HB"]
        results.append(check("period change during the run is one SET PERIOD",
                             cmds == ["SET PERIOD 200.0000"], str(cmds)))
        good += run_records(app, gui, arduino, 20)

        # 9a. new data file refused during the run
        with contextlib.redirect_stdout(io.StringIO()):
            gui.new_data_file()
            app.processEvents()
        results.append(check("new data file refused while running",
                             gui.files.csv_path == first_csv and not gui.new_file_button.isEnabled()
                             and warnings and warnings[-1][0] == "Run in progress"))

        # 6. ERROR reply (Ta + Tb = 58 C above the 56 C limit) makes the next
        #    "Send settings" resend everything
        gui.spins["MEAN"].setValue(54.0)
        quiet(gui.send_button.click)
        good += pump(app, gui, arduino)
        results.append(check("ERROR reply: shown, and every setting marked unsent",
                             "ERROR" in gui.lbl_message.text()
                             and f"{len(gui_module.SETTINGS) + 2} setting(s) not sent yet" in gui.lbl_pending.text(),
                             gui.lbl_pending.text()))
        gui.spins["MEAN"].setValue(40.0)

        # 7. safety record
        n0 = len(arduino.sent)
        good += run_records(app, gui, arduino, 1, state_override="SAFETY")
        cmds = [c for c in arduino.sent[n0:] if c != "HB"]
        results.append(check("SAFETY shown in red, no automatic restart",
                             gui.lbl_state.text() == "State: SAFETY"
                             and "red" in gui.lbl_state.styleSheet() and cmds == []))
        quiet(gui.stop_button.click)
        good += pump(app, gui, arduino)
        good += run_records(app, gui, arduino, 3)

        # 8. Arduino restart
        arduino.reset()
        n0 = len(arduino.sent)
        good += run_records(app, gui, arduino, 2)
        cmds = [c for c in arduino.sent[n0:] if c != "HB"]
        results.append(check("restart logged and every setting resent",
                             len(cmds) == len(gui_module.SETTINGS) + 3 and cmds[-1] == "GET", f"{len(cmds)} commands"))

        # 4-5. raw CSV and metadata
        gui.note_edit.setText("self-test note")
        quiet(gui.note_button.click)
        quiet(gui.new_data_file)
        second_csv = gui.files.csv_path
        with open(first_csv, newline="") as f:
            rows = list(csv.reader(f))
        expected = [line.split(",")[1:] for line in good]
        results.append(check("CSV header is the 18 DATA columns", rows[0] == DATA_COLUMNS))
        results.append(check("CSV rows identical to the accepted records",
                             rows[1:] == expected, f"{len(rows) - 1} rows"))
        meta = first_csv.with_name(first_csv.stem + "_meta.txt").read_text()
        wanted = ("T3_C: T3, input A3, x = 5.0 cm", "START mode=PI", "SET PERIOD 200.0000",
                  "CONFIG mode=MANUAL", "ERROR: Ta - Tb or Ta + Tb", "Arduino restarted",
                  "NOTE: self-test note", second_csv.name, "analysis length L")
        missing = [s for s in wanted if s not in meta]
        results.append(check("metadata has channel map, events, config, errors, notes",
                             not missing, f"missing {missing}" if missing else ""))

        # 9b. new file receives records, no overwrite; 10. close sends STOP
        good2 = pump(app, gui, arduino) + run_records(app, gui, arduino, 3)
        third = gui_module.RunFiles(Path(tmp))
        third.close()
        n0 = len(arduino.sent)
        quiet(gui.close)
        with open(second_csv, newline="") as f:
            rows2 = list(csv.reader(f))
        results.append(check("new data file receives the next records",
                             second_csv != first_csv and len(rows2) == 1 + len(good2)))
        results.append(check("a file created in the same second is not overwritten",
                             third.csv_path not in (first_csv, second_csv), third.csv_path.name))
        results.append(check("closing the GUI sends STOP", "STOP" in arduino.sent[n0:]))

        # 11. automatic schedule
        sched_path = Path(tmp) / "schedule.csv"
        sched_path.write_text(
            "label,kind,mode,period_s,Ta_C,Tb_C,Kp,Ki,settle_on,cycles_after_settle,min_s,max_s,notes\n"
            "# comment row,,,,,,,,,,,,\n"
            "base,baseline,OFF,,,,,,,,60,,\n"
            "pi100,production,PI,100,40,10,6,0.04,TEC,3,300,3000,\n"
            "pi200,production,PI,200,,,,,TEC,3,600,4000,\n")
        arduino = SimulatedArduino()
        gui = quiet(gui_module.AngstromGUI, arduino)
        pump(app, gui, arduino)
        run_records(app, gui, arduino, 3)
        ok_load = quiet(gui.load_schedule, sched_path)
        gui.spins["TROOM"].setValue(25.0)       # the baseline (22 C) must replace this
        quiet(gui.send_settings)
        pump(app, gui, arduino)
        n0 = len(arduino.sent)
        quiet(gui.start_schedule, False)
        steps_seen = []
        for _ in range(20000):
            run_records(app, gui, arduino, 1)
            steps_seen.append(arduino.step)
            if not gui.runner.active:
                break
        cmds = [c for c in arduino.sent[n0:] if c != "HB"]
        i100 = cmds.index("START") if "START" in cmds else -1
        after = cmds[i100 + 1:]
        results.append(check("schedule loads and runs to completion",
                             ok_load and not gui.runner.active and cmds[0] == "STOP"
                             and "Schedule complete" in gui.lbl_schedule.text(),
                             gui.lbl_schedule.text()[:60]))
        results.append(check("PI row: mode, settings, T_room from the baseline, then START",
                             "MODE PI" in cmds[:i100] and "SET AMP 10.0000" in cmds[:i100]
                             and any(c.startswith("SET TROOM 22.") for c in cmds[:i100]),
                             str(cmds[:i100])))
        results.append(check("period change during the schedule is one SET PERIOD, then STOP",
                             after == ["SET PERIOD 200.0000", "STOP"], str(after)))
        rows_b = list(csv.DictReader(open(Path(tmp) / "intervals_auto_baseline.csv")))
        rows_p = list(csv.DictReader(open(Path(tmp) / "intervals_auto_production.csv")))
        spans = [(float(r["t_end_s"]) - float(r["t_start_s"])) for r in rows_p]
        results.append(check("segments written to intervals_auto_*.csv",
                             len(rows_b) == 1 and len(rows_p) == 2
                             and all("; settled" in r["notes"] for r in rows_p)
                             and spans[0] >= 300 and spans[1] >= 600
                             and rows_p[0]["step"] != rows_p[1]["step"],
                             f"production spans {spans} s, steps "
                             f"{[r['step'] for r in rows_p]}"))
        sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
        import module_08_control_fit as cfit
        res = quiet(cfit.main, ["--intervals", str(Path(tmp) / "intervals_auto_production.csv"),
                                "--fig-dir", str(Path(tmp) / "figs"), "--out-dir", tmp])
        results.append(check("control analysis reads the automatic intervals",
                             len(res) == 2 and [r["period_s"] for r in res] == [100.0, 200.0]
                             and all(r["n_cycles"] >= 3 for r in res),
                             ", ".join(f"{r['label']}: {r['n_cycles']} cycles, amp ratio "
                                       f"{r['amp_ratio']:.2f}" for r in res)))

        # 11b. a glitch (all channels drop 2 C for one record) restarts the segment and
        #      stays out of the interval written
        sched_path.write_text(
            "label,kind,mode,period_s,Ta_C,Tb_C,settle_on,cycles_after_settle,min_s,max_s\n"
            "pi100g,production,PI,100,40,10,TEC,3,300,3000\n")
        quiet(gui.load_schedule, sched_path)
        quiet(gui.start_schedule, False)
        run_records(app, gui, arduino, 40)
        arduino.temps = [v - 2.0 for v in arduino.temps]
        run_records(app, gui, arduino, 1)
        t_glitch = arduino.t_ms / 1000.0
        arduino.temps = [a + 2.0 for a in arduino.temps]   # contact restored
        for _ in range(4000):
            run_records(app, gui, arduino, 1)
            if not gui.runner.active:
                break
        rows_g = [r for r in csv.DictReader(open(Path(tmp) / "intervals_auto_production.csv"))
                  if r["label"] == "pi100g"]
        meta = gui.files.meta_path.read_text()
        results.append(check("glitch logged, segment restarted, interval starts after it",
                             "glitch at" in meta and "common-mode" in meta and len(rows_g) == 1
                             and float(rows_g[0]["t_start_s"]) > t_glitch + 1,
                             f"glitch at {t_glitch:.0f} s, interval from {rows_g[0]['t_start_s'] if rows_g else '?'} s"))

        # 12. SAFETY aborts a running schedule
        quiet(gui.start_schedule, False)
        run_records(app, gui, arduino, 5)
        n0 = len(arduino.sent)
        run_records(app, gui, arduino, 1, state_override="SAFETY")
        cmds = [c for c in arduino.sent[n0:] if c != "HB"]
        results.append(check("SAFETY aborts the schedule and sends STOP",
                             not gui.runner.active and cmds == ["STOP"]
                             and "SAFETY" in gui.lbl_schedule.text(), str(cmds)))
        quiet(gui.close)

        # 13. command pacing
        gui_module.COMMAND_SPACING_MS = 60
        arduino = SimulatedArduino()
        gui = quiet(gui_module.AngstromGUI, arduino)
        pump(app, gui, arduino)
        run_records(app, gui, arduino, 1)          # sync settings: queued
        n_queued = len(gui.cmd_queue)
        while gui.cmd_queue:
            quiet(gui.send_next_command)
        pump(app, gui, arduino)
        n0 = len(arduino.sent)
        gui.mode_combo.setCurrentText("PI")
        quiet(gui.start_button.click)
        queued = [c for c, _ in gui.cmd_queue]
        quiet(gui.send_next_command)                # only MODE PI goes out
        first = [c for c in arduino.sent[n0:] if c != "HB"]
        quiet(gui.stop_button.click)
        after_stop = [c for c in arduino.sent[n0:] if c != "HB"]
        results.append(check("pacing: commands queued in order, STOP immediate and clears the queue",
                             n_queued == len(gui_module.SETTINGS) + 3 and queued == ["MODE PI", "START"]
                             and first == ["MODE PI"] and after_stop == ["MODE PI", "STOP"]
                             and not gui.cmd_queue,
                             f"sync queued {n_queued}; then {queued}; sent {after_stop}"))

        # 14. echo checks
        pump(app, gui, arduino)
        gui.spins["PERIOD"].setValue(200)
        quiet(gui.send_settings)
        quiet(gui.send_next_command)                # write SET PERIOD; the simulator echoes it
        arduino.out = [x for x in arduino.out if "SET PERIOD" not in x]   # lose that echo
        quiet(gui.on_line, "EVENT,99000,40,SET PERIOD 150.0000 (was 100.0000)")
        bad_echo = "GUI check" in gui.lbl_message.text() and not gui.pending_echo
        gui.lbl_message.setText("")
        gui.pending_echo = {"KP": [["7.0000", time.monotonic() - 10]]}
        quiet(gui.check_missing_echoes)
        missing = "no reply to SET KP" in gui.lbl_message.text() and not gui.pending_echo
        results.append(check("wrong echo and missing echo both flagged as errors",
                             bad_echo and missing))

        # 15. schedule started while the Arduino is in SAFETY
        arduino.state = "SAFETY"
        pump(app, gui, arduino)
        run_records(app, gui, arduino, 2)
        sched_path.write_text(
            "label,kind,mode,period_s,Ta_C,Tb_C,settle_on,cycles_after_settle,min_s,max_s\n"
            "pi100,production,PI,100,40,10,TEC,2,200,2000\n")
        quiet(gui.load_schedule, sched_path)
        quiet(gui.start_schedule, False)
        run_records(app, gui, arduino, 2)           # START still queued: records say SAFETY
        still = gui.runner.active
        while gui.cmd_queue:
            quiet(gui.send_next_command)
        for _ in range(60):
            run_records(app, gui, arduino, 1)
        results.append(check("schedule waits through a leftover SAFETY until its START",
                             still and gui.runner.active and arduino.state == "RUN",
                             gui.lbl_schedule.text()[:70]))
        quiet(gui.abort_schedule, "end of test")
        quiet(gui.close)
        gui_module.COMMAND_SPACING_MS = 0

        # 16. automatic S18.1 checks
        def run_until_done(gui, arduino, limit=5000):
            for _ in range(limit):
                run_records(app, gui, arduino, 1)
                if not gui.runner.active:
                    break

        arduino = SimulatedArduino()
        gui = quiet(gui_module.AngstromGUI, arduino)
        pump(app, gui, arduino)
        run_records(app, gui, arduino, 3)
        sched_path.write_text(
            "label,kind,mode,manual_pwm,min_s,max_s,limit_rise_C,max_spread_C\n"
            "check_channels,check,OFF,,30,,,1.0\n"
            "test_tec,check,TEST_TEC_LIMIT,40,,300,2,\n"
            "test_t1,check,TEST_T1_LIMIT,40,,600,2,\n"
            "cool,check,OFF,,20,,,\n")
        quiet(gui.load_schedule, sched_path)
        n0 = len(arduino.sent)
        quiet(gui.start_schedule, False)
        run_until_done(gui, arduino)
        cmds = [c for c in arduino.sent[n0:] if c != "HB"]
        meta = gui.files.meta_path.read_text()
        lowered = [c for c in cmds if c.startswith(("SET T1LIMIT", "SET TECLIMIT"))]
        results.append(check("limit tests: lowered to reading + 2 C, SAFETY, PASS, limits restored",
                             "Schedule complete" in gui.lbl_schedule.text()
                             and meta.count(": PASS") == 2
                             and lowered[0].startswith("SET TECLIMIT 24.")
                             and any(c.startswith("SET T1LIMIT 2") and c != "SET T1LIMIT 56.0000"
                                     for c in lowered)
                             and lowered[-2:] in (["SET T1LIMIT 56.0000", "SET TECLIMIT 60.0000"],
                                                  ["SET TECLIMIT 60.0000", "SET T1LIMIT 56.0000"])
                             and arduino.cfg["T1LIMIT"] == 56.0 and arduino.cfg["TECLIMIT"] == 60.0
                             and "channel spread 0.0" in meta,
                             f"{lowered}"))
        arduino.temps[4] += 2.5                    # T4 reads 2.5 C high
        quiet(gui.start_schedule, False)
        run_until_done(gui, arduino)
        results.append(check("channel spread above 1 C stops the schedule",
                             "spread" in gui.lbl_schedule.text() and "exceeds" in gui.lbl_schedule.text(),
                             gui.lbl_schedule.text()[:70]))
        quiet(gui.close)

        # 17. autorun with the analysis afterwards
        arduino = SimulatedArduino()
        gui = quiet(gui_module.AngstromGUI, arduino)
        sched_path.write_text(
            "label,kind,mode,period_s,Ta_C,Tb_C,settle_on,cycles_after_settle,min_s,max_s\n"
            "room,baseline,OFF,,,,,,30,\n"
            "pi100,control,PI,100,40,10,TEC,3,300,3000\n")
        gui.autorun = {"schedule": sched_path, "yes": True, "analyse": True}
        pump(app, gui, arduino)
        for _ in range(4000):
            run_records(app, gui, arduino, 1)
            if gui.runner is not None and not gui.runner.active:
                break
        procs = list(gui.analysis_procs)
        for proc in procs:
            proc.wait(timeout=300)
        logs = sorted((Path(tmp) / "analysis").glob("*_control.log"))
        log_text = logs[0].read_text() if logs else ""
        results.append(check("autorun: schedule started by itself, control analysis ran afterwards",
                             gui.autorun.get("started") and "Schedule complete" in gui.lbl_schedule.text()
                             and len(procs) == 1 and procs[0].returncode == 0
                             and "| pi100 | PI | 100 |" in log_text,
                             f"{len(procs)} analysis process(es), log {logs[0].name if logs else 'missing'}"))
        quiet(gui.close)

    # 18. heartbeat from the serial thread, with the watchdog
    class Link:
        def __init__(self):
            self.sent = []

        def write_line(self, text):
            self.sent.append(text)

    with tempfile.TemporaryDirectory() as tmp2:
        gui_module.DATA_DIR = Path(tmp2)
        arduino = SimulatedArduino()
        gui = quiet(gui_module.AngstromGUI, arduino)
        link = Link()
        reader = gui_module.SerialReader(link, gui.reader.heartbeat_ok)
        data_line = b"DATA,1000,RUN,PI,5,40,16,100,41,nan,40,35,30,28,26,25,3,2.82,20\n"
        reader.handle_raw(data_line)
        reader.handle_raw(b"EVENT,1000,6,SET KP 6.0000 (was 6.0000)\n")
        fresh = list(link.sent)
        gui.last_processed -= gui_module.HEARTBEAT_WATCHDOG_S + 1     # window stalled 31 s
        reader.handle_raw(data_line)
        results.append(check("serial thread answers DATA with HB; stops when the window stalls > 30 s",
                             fresh == ["HB"] and link.sent == ["HB"], f"{fresh} then {link.sent}"))
        quiet(gui.close)

    passed = all(results)
    print("PASS" if passed else "FAIL")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
