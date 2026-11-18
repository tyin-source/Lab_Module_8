# Phys 39 Module 8 - six-thermistor Angstrom GUI
# (extends the Module 7 GUI, python/module_07_rod_acquisition_gui.py in the Module 7
# repository, which stays unchanged as the working checkpoint)
#
# Use with arduino/module_08_angstrom_controller/module_08_angstrom_controller.ino
# (115200 baud).
#
# Changes from Module 7:
#   - The drive and the P/PI law now run on the Arduino (see the sketch header).
#     This program sends the settings, START/STOP, and a heartbeat (HB) after
#     every record, and displays and saves what the Arduino reports. It no longer
#     computes the control law, so the saved record is exactly what the Arduino
#     applied.
#   - Parses the 18-field Module 8 record (S18.2):
#       DATA,time_ms,run_state,control_mode,step,Ta_C,Tb_C,period_s,setpoint_C,
#            feedforward_pwm,TEC_C,T1_C,T2_C,T3_C,T4_C,T5_C,controller_pwm,
#            heating_to_cooling_ratio,signed_pwm
#     and saves every accepted record unchanged to data/module_08/angstrom_<stamp>.csv
#     (flushed after every record; an existing file is never overwritten).
#   - Copies the Arduino's EVENT lines (start, stop, every setting, mode or period
#     change, safety, timeout; each with its step number) and CONFIG lines into the
#     run's metadata file, together with the GUI's own notes.
#   - Five modes: MANUAL, PWM_SINE (raw sinusoidal PWM, zero-centered or
#     positive-only), FF, P, PI. "Send settings" sends only the values that
#     changed, so a period change during a run is a single SET PERIOD command and
#     the Arduino continues the waveform without a phase jump.
#   - Live display of all six temperatures, the requested setpoint, run state,
#     mode, step, and the applied, feedforward, and controller PWM.
#   - Commands are paced (COMMAND_SPACING_MS) so the Uno's 63-byte receive buffer
#     never overflows; HB and STOP go out at once, and STOP discards anything queued.
#     Every SET is checked against the Arduino's echo: a different value, or no echo
#     within ECHO_TIMEOUT_S, is handled like an ERROR (everything is resent next time,
#     and a running schedule stops).
#
#   - Automatic schedule (optional): "Load schedule..." reads a CSV from schedules/
#     and "Start schedule" runs its rows in order through the same controls: each
#     row's settings are sent, the run is started if needed, and the row ends when the
#     live settling test has found a settled start and the requested whole cycles have
#     been recorded after it (or after a fixed time for OFF/MANUAL rows). Completed
#     segments are added to data/module_08/intervals_auto_<kind>.csv for the analysis
#     scripts. See python/module_08_schedule.py and schedules/README.md.
#
# The GUI keeps the Arduino running only while it is receiving records: if the
# GUI freezes or the USB link drops, the heartbeat stops and the Arduino ends the
# run after 3 s (run_state TIMEOUT). A schedule stops at once (and sends STOP) on
# SAFETY, TIMEOUT, any ERROR reply, the Stop button, or closing the window; it never
# restarts a run by itself.

import csv
import math
import re
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path

import serial
from PySide6 import QtCore, QtWidgets
import pyqtgraph as pg

import module_08_schedule as sched

from module_08_common import (ARDUINO_INPUTS, CONTROL_MODES, DATA_COLUMNS, INT_COLUMNS,
                              MODULE8_DATA_DIR, POSITION_UNCERTAINTY_MM, REPO_ROOT,
                              ROD_ANALYSIS_LENGTH_M, ROD_DIAMETER_M, ROD_OVERALL_LENGTH_M,
                              RUN_STATES, SENSOR_LABELS, SENSOR_POSITIONS_CM,
                              TEMPERATURE_COLUMNS, TEXT_COLUMNS)


# ------------------------------------------------------------------
# USER CONFIGURATION  (edit these near the top)
# ------------------------------------------------------------------
SERIAL_PORT = "/dev/cu.usbmodem1101"          # Mac. Windows example: "COM5"
BAUD_RATE = 115200

PLOT_UPDATE_MS = 200          # how often the plots redraw (ms)
# Commands other than HB and STOP are queued and sent one every COMMAND_SPACING_MS.
# The Uno's serial receive buffer holds 63 bytes, and the sketch prints an EVENT line
# (~50 characters, ~4 ms at 115200 baud) for each setting it accepts; a burst of
# settings would overflow the buffer and lose characters (a lost digit can turn
# "SET PWMAMP 50" into a valid "SET PWMAMP 5"). 0 sends at once (used by the self-test).
COMMAND_SPACING_MS = 60
ECHO_TIMEOUT_S = 3.0          # an accepted SET is echoed at once; none after this = lost
# The serial reader thread answers every DATA line with HB itself, so a run does not
# stop when the window's event loop is slowed for a few seconds (macOS App Nap with the
# display asleep or the window covered stopped the first S18.5 run after 10 min). It
# does so only while the window has handled a record within HEARTBEAT_WATCHDOG_S: a GUI
# that is really frozen stops sending HB, and the Arduino stops the run 3 s later.
HEARTBEAT_WATCHDOG_S = 30.0
TEMP_MIN_Y_SPAN = 0.2         # smallest temperature-axis span (C)

DATA_DIR = MODULE8_DATA_DIR

# Settings sent to the Arduino: (command name, label, minimum, maximum, decimals,
# step, default). Ranges match the checks in the sketch, which rejects anything
# outside them (Ta - Tb >= 10 C and Ta + Tb <= 56 C are checked there). Defaults:
# S18.5 setpoint 40 + 16 sin(wt) C at 100 s, S18.3 u_b = 50, Module 4 chi_H and r,
# Module 7 PI gains.
SETTINGS = [
    ("MEAN",   "Ta, mean (C)",             10.0, 56.0,   2, 0.5,   40.0),
    ("AMP",    "Tb, amplitude (C)",         0.0, 23.0,   2, 0.5,   16.0),
    ("PERIOD", "Period (s)",              100.0, 3600.0, 1, 10.0, 100.0),
    ("KP",     "Kp (PWM/C)",                0.0, 100.0,  3, 0.5,    6.0),
    ("KI",     "Ki (PWM/(C s))",            0.0, 5.0,    4, 0.01,   0.04),
    ("UMAX",   "u_max (PWM)",               0.0, 255.0,  0, 5.0,  255.0),
    ("RATIO",  "r = chi_H/chi_C",           1.0, 10.0,   3, 0.05,   2.82),
    ("CHIH",   "chi_H (C/PWM)",             0.01, 5.0,   4, 0.005,  0.485),
    ("TROOM",  "T_room (C)",                0.0, 40.0,   2, 0.1,   22.0),
    ("PWMOFF", "PWM_SINE u_a (offset)",  -255.0, 255.0,  0, 1.0,    0.0),
    ("PWMAMP", "PWM_SINE u_b (amplitude)",  0.0, 255.0,  0, 1.0,   50.0),
    ("MANUAL", "MANUAL PWM (signed)",    -255.0, 255.0,  0, 1.0,    0.0),
    ("T1LIMIT", "T1 limit (C)",             0.0, 56.0,   2, 0.5,   56.0),
    ("TECLIMIT", "TEC limit (C)",           0.0, 60.0,   2, 0.5,   60.0),
]
PERIOD_PRESETS_S = (100, 200, 300, 400, 500)

# Plot style for each channel: (color, line width). TEC is thicker.
TRACE_STYLE = {
    "TEC_C": ("k", 3.0),
    "T1_C": ((214, 39, 40), 1.6),
    "T2_C": ((255, 127, 14), 1.6),
    "T3_C": ((44, 160, 44), 1.6),
    "T4_C": ((31, 119, 180), 1.6),
    "T5_C": ((148, 103, 189), 1.6),
}
STATE_STYLE = {
    "OFF": "",
    "RUN": "color: white; background-color: green; font-weight: bold; padding: 2px;",
    "SAFETY": "color: white; background-color: red; font-weight: bold; padding: 2px;",
    "TIMEOUT": "color: white; background-color: darkorange; font-weight: bold; padding: 2px;",
}
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# SERIAL LINK  (unchanged from Modules 6 and 7)
# ------------------------------------------------------------------
class SerialLink:
    def __init__(self, port: str, baud: int):
        self.ser = serial.Serial(port, baud, timeout=0.5)
        self.write_lock = threading.Lock()

    def readline(self) -> bytes:
        return self.ser.readline()

    def write_line(self, text: str):
        with self.write_lock:
            self.ser.write((text + "\n").encode("ascii"))
            self.ser.flush()

    def close(self):
        with self.write_lock:
            try:
                self.ser.close()
            except Exception:
                pass


# ------------------------------------------------------------------
# LINE PARSER
# ------------------------------------------------------------------
# classify_line() sorts every serial line into:
#   ("record", raw_fields, values)  a complete DATA line
#   ("header", None, None)          the DATA,time_ms,... header
#   ("event", None, (time_ms, step, text))
#   ("config", None, (time_ms, step, text))   reply to GET
#   ("rejected", None, None)        starts with DATA but cannot be read
#   ("other", None, None)           start-up text, ERROR lines, anything else
def _to_float(text: str) -> float:
    value = float(text)                 # accepts "nan"
    if math.isinf(value):
        raise ValueError("infinite value")
    return value


def classify_line(line: str):
    for tag, kind in (("EVENT,", "event"), ("CONFIG,", "config")):
        if line.startswith(tag):
            parts = line.split(",", 3)
            try:
                return kind, None, (int(parts[1]), int(parts[2]), parts[3])
            except (IndexError, ValueError):
                return "other", None, None
    if not line.startswith("DATA,"):
        return "other", None, None
    fields = [f.strip() for f in line.split(",")[1:]]
    if fields == DATA_COLUMNS:
        return "header", None, None
    if len(fields) != len(DATA_COLUMNS):
        return "rejected", None, None
    values = {}
    try:
        for name, text in zip(DATA_COLUMNS, fields):
            if name in TEXT_COLUMNS:
                values[name] = text
            elif name in INT_COLUMNS:
                values[name] = int(text)
            else:
                values[name] = _to_float(text)
    except ValueError:
        return "rejected", None, None
    if (values["run_state"] not in RUN_STATES or values["control_mode"] not in CONTROL_MODES
            or not -255 <= values["signed_pwm"] <= 255 or values["time_ms"] < 0):
        return "rejected", None, None
    return "record", fields, values


class SerialReader(QtCore.QThread):
    line_received = QtCore.Signal(str)
    error = QtCore.Signal(str)

    def __init__(self, link, heartbeat_ok=None):
        super().__init__()
        self.link = link
        self.heartbeat_ok = heartbeat_ok      # callable: may this thread send HB now?
        self._stop = False

    def handle_raw(self, raw):
        line = raw.decode("utf-8", errors="replace").strip()
        if not line:
            return
        if line.startswith("DATA,") and self.heartbeat_ok is not None and self.heartbeat_ok():
            try:
                self.link.write_line("HB")
            except Exception:
                pass                          # reported by the reader loop if the port is gone
        self.line_received.emit(line)

    def run(self):
        try:
            while not self._stop:
                raw = self.link.readline()
                if raw:
                    self.handle_raw(raw)
        except Exception as e:
            self.error.emit(str(e))

    def stop(self):
        self._stop = True
        self.wait(1000)


# ------------------------------------------------------------------
# RAW DATA FILE + METADATA
# ------------------------------------------------------------------
def git_commit() -> str:
    """Short commit hash of the repository, with '+uncommitted changes' if dirty."""
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT,
                                capture_output=True, text=True, timeout=5).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=REPO_ROOT,
                               capture_output=True, text=True, timeout=5).stdout.strip()
        return (commit or "unknown") + (" +uncommitted changes" if dirty else "")
    except Exception:
        return "unknown"


class RunFiles:
    """One raw CSV and its metadata file. Opened with mode 'x', so an existing
    file is never overwritten."""

    def __init__(self, data_dir: Path):
        data_dir.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        suffix = 0
        while True:
            name = f"angstrom_{stamp}" + (f"_{suffix}" if suffix else "")
            self.csv_path = data_dir / f"{name}.csv"
            self.meta_path = data_dir / f"{name}_meta.txt"
            try:
                self.csv_file = open(self.csv_path, "x", newline="")
                break
            except FileExistsError:
                suffix += 1
        self.meta_file = open(self.meta_path, "x")
        self.rows = 0

        self.csv_file.write(",".join(DATA_COLUMNS) + "\n")
        self.csv_file.flush()

        m = self.meta_file
        m.write(f"Raw data file: {self.csv_path.name}\n")
        m.write(f"Created: {datetime.now().isoformat(timespec='seconds')}\n")
        m.write(f"Serial port: {SERIAL_PORT} at {BAUD_RATE} baud\n")
        m.write("Arduino sketch: arduino/module_08_angstrom_controller/"
                "module_08_angstrom_controller.ino\n")
        m.write("Python GUI: python/module_08_angstrom_gui.py\n")
        m.write(f"Git commit: {git_commit()}\n")
        m.write("Record format and units: see data/module_08/README.md. time_ms is the "
                "Arduino clock since its last reset; *_C in degrees C; PWM in counts, "
                "cooling negative; nan = field not used in this control mode.\n")
        m.write("Temperatures are raw readings (no calibration corrections).\n")
        m.write(f"Rod: 6061-T6, d = {ROD_DIAMETER_M * 1e3:.2f} mm, overall length "
                f"{ROD_OVERALL_LENGTH_M * 100:.2f} cm, analysis length L (first hole to tip) "
                f"{ROD_ANALYSIS_LENGTH_M * 100:.2f} cm (Module 7).\n")
        m.write(f"Channel map (x from the first thermistor hole, +/- {POSITION_UNCERTAINTY_MM} mm):\n")
        for col in TEMPERATURE_COLUMNS:
            where = ("TEC base" if col == "TEC_C"
                     else f"x = {SENSOR_POSITIONS_CM[col]:.1f} cm")
            m.write(f"  {col}: {SENSOR_LABELS[col]}, input {ARDUINO_INPUTS[col]}, {where}\n")
        m.write("Events (Arduino time, wall-clock time, then the Arduino step number):\n")
        m.flush()

    def write_record(self, raw_fields):
        self.csv_file.write(",".join(raw_fields) + "\n")
        self.csv_file.flush()
        self.rows += 1

    def event(self, time_ms, text: str):
        t = "      ---  " if time_ms is None else f"{time_ms / 1000.0:9.2f} s"
        self.meta_file.write(f"  [{t} | {datetime.now().strftime('%H:%M:%S')}] {text}\n")
        self.meta_file.flush()

    def close(self):
        for f in (self.csv_file, self.meta_file):
            try:
                f.close()
            except Exception:
                pass


def start_keep_awake():
    """macOS: keep the display and system awake while the GUI runs (caffeinate, ends
    with this process) and, if pyobjc is installed, opt out of App Nap. Elsewhere: no-op.
    Returns handles for stop_keep_awake()."""
    handles = {}
    if sys.platform != "darwin":
        return handles
    import os
    import shutil
    if shutil.which("caffeinate"):
        try:
            handles["caffeinate"] = subprocess.Popen(
                ["caffeinate", "-dims", "-w", str(os.getpid())],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print("[power] caffeinate: display and system kept awake while the GUI runs")
        except OSError:
            pass
    try:
        from Foundation import NSProcessInfo          # pyobjc, optional
        options = 0x00FFFFFF | 0xFF00000000           # user initiated, latency critical
        handles["activity"] = NSProcessInfo.processInfo().beginActivityWithOptions_reason_(
            options, "Module 8 data acquisition")
        print("[power] App Nap disabled for this process")
    except Exception:
        print("[power] pyobjc not installed: App Nap may still slow the window; the "
              "heartbeat comes from the serial thread, so the run continues")
    return handles


def stop_keep_awake(handles):
    p = handles.get("caffeinate")
    if p is not None and p.poll() is None:
        p.terminate()


def format_command_value(value: float) -> str:
    """A number the sketch's parser accepts: optional '-', digits, one '.'."""
    text = f"{value:.4f}"
    return "0.0000" if text == "-0.0000" else text


# ------------------------------------------------------------------
# GUI + PLOTS
# ------------------------------------------------------------------
class AngstromGUI(QtWidgets.QMainWindow):
    def __init__(self, link=None):
        super().__init__()
        self.setWindowTitle("Module 8 - Angstrom Method: Six-Thermistor Record and TEC Drive")

        self.last_time_ms = None      # Arduino time of the previous record
        self.t_offset = 0.0           # plot time = Arduino time + offset (grows if the Arduino resets)
        self.last_state = None        # run_state of the previous record
        self.sent = {}                # last value sent for each setting (and "MODE")
        self.pending_echo = {}        # setting -> [value, time written or None] awaiting echo
        self.synced = False           # settings sent since the Arduino (re)started
        self.schedule_rows = None     # loaded schedule (list of ScheduleRow)
        self.schedule_path = None
        self.runner = None            # ScheduleRunner while a schedule runs
        self.last_rec = None          # values of the latest accepted record
        self.autorun = None           # set by module_08_autorun.py: schedule to start
        self.analysis_procs = []      # analysis started after an automatic run

        self.n_accepted = 0
        self.n_rejected = 0

        # whole-session plot data
        self.ts = []
        self.temps = {col: [] for col in TEMPERATURE_COLUMNS}
        self.setpoints = []
        self.signed = []
        self.feedforward = []
        self.controller = []
        self.errors = []

        # open the serial link first, so a wrong port leaves no empty file
        self.link = link if link is not None else SerialLink(SERIAL_PORT, BAUD_RATE)
        self.files = RunFiles(DATA_DIR)
        print(f"[data] recording to {self.files.csv_path}")

        self._build_ui()
        self.update_file_label()

        self.last_processed = time.monotonic()   # last record handled by the window
        self.keep_awake = start_keep_awake()
        self.reader = SerialReader(
            self.link, lambda: time.monotonic() - self.last_processed < HEARTBEAT_WATCHDOG_S)
        self.reader.line_received.connect(self.on_line)
        self.reader.error.connect(self.on_serial_error)
        self.reader.start()

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_plots)
        self.timer.start(PLOT_UPDATE_MS)

        self.cmd_queue = []           # commands waiting to be sent (see COMMAND_SPACING_MS)
        self.cmd_timer = QtCore.QTimer(self)
        self.cmd_timer.timeout.connect(self.send_next_command)
        if COMMAND_SPACING_MS > 0:
            self.cmd_timer.start(COMMAND_SPACING_MS)

    # ---------------- build the GUI ----------------
    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)

        # --- settings ---
        box = QtWidgets.QGroupBox(
            "Drive (computed on the Arduino):  Tset = Ta + Tb sin(phase);  "
            "PWM_SINE: u = u_a + u_b sin(phase);  FF: u = (Tset - T_room)/chi_H;  "
            "P/PI: u = Kp e + Ki q, e = Tset - T_TEC;  cooling (u < 0) x r in FF, P, PI")
        grid = QtWidgets.QGridLayout(box)

        grid.addWidget(QtWidgets.QLabel("Mode:"), 0, 0)
        self.mode_combo = QtWidgets.QComboBox()
        self.mode_combo.addItems(CONTROL_MODES)
        grid.addWidget(self.mode_combo, 0, 1)
        self.ffadd_check = QtWidgets.QCheckBox("Add feedforward to P/PI")
        grid.addWidget(self.ffadd_check, 0, 2, 1, 2)

        grid.addWidget(QtWidgets.QLabel("Period presets:"), 0, 4)
        preset_row = QtWidgets.QHBoxLayout()
        for p in PERIOD_PRESETS_S:
            b = QtWidgets.QPushButton(f"{p} s")
            b.clicked.connect(lambda _=False, p=p: self.spins["PERIOD"].setValue(p))
            preset_row.addWidget(b)
        grid.addLayout(preset_row, 0, 5, 1, 3)

        self.spins = {}
        for i, (name, label, low, high, decimals, step, default) in enumerate(SETTINGS):
            row, col = 1 + i // 4, (i % 4) * 2
            grid.addWidget(QtWidgets.QLabel(label + ":"), row, col)
            spin = QtWidgets.QDoubleSpinBox()
            spin.setRange(low, high)
            spin.setDecimals(decimals)
            spin.setSingleStep(step)
            spin.setValue(default)
            self.spins[name] = spin
            grid.addWidget(spin, row, col + 1)

        buttons = QtWidgets.QHBoxLayout()
        self.send_button = QtWidgets.QPushButton("Send settings")
        self.send_button.clicked.connect(self.send_settings)
        self.start_button = QtWidgets.QPushButton("Start")
        self.start_button.clicked.connect(self.start_run)
        self.stop_button = QtWidgets.QPushButton("Stop (PWM 0)")
        self.stop_button.clicked.connect(self.stop_run)
        self.zero_button = QtWidgets.QPushButton("Zero integral")
        self.zero_button.clicked.connect(lambda: self.send_line("ZERO"))
        self.lbl_pending = QtWidgets.QLabel("")
        for w in (self.send_button, self.start_button, self.stop_button, self.zero_button):
            buttons.addWidget(w)
        buttons.addWidget(self.lbl_pending)
        buttons.addStretch(1)
        grid.addLayout(buttons, 1 + math.ceil(len(SETTINGS) / 4), 0, 1, 8)   # below the settings
        layout.addWidget(box)

        # --- automatic schedule ---
        sbox = QtWidgets.QGroupBox("Automatic schedule (stay with the apparatus while it runs)")
        srow = QtWidgets.QHBoxLayout(sbox)
        self.load_schedule_button = QtWidgets.QPushButton("Load schedule...")
        self.load_schedule_button.clicked.connect(lambda: self.load_schedule())
        self.start_schedule_button = QtWidgets.QPushButton("Start schedule")
        self.start_schedule_button.setEnabled(False)
        self.start_schedule_button.clicked.connect(lambda: self.start_schedule())
        self.abort_schedule_button = QtWidgets.QPushButton("Abort schedule (PWM 0)")
        self.abort_schedule_button.setEnabled(False)
        self.abort_schedule_button.clicked.connect(lambda: self.abort_schedule("aborted by user"))
        self.lbl_schedule = QtWidgets.QLabel("No schedule loaded")
        for w in (self.load_schedule_button, self.start_schedule_button,
                  self.abort_schedule_button):
            srow.addWidget(w)
        srow.addWidget(self.lbl_schedule, stretch=1)
        layout.addWidget(sbox)

        for spin in self.spins.values():
            spin.valueChanged.connect(self.update_pending)
        self.mode_combo.currentIndexChanged.connect(self.update_pending)
        self.ffadd_check.toggled.connect(self.update_pending)

        # --- data file, counters, notes ---
        file_box = QtWidgets.QGroupBox("Raw data")
        frow = QtWidgets.QHBoxLayout(file_box)
        self.lbl_file = QtWidgets.QLabel()
        frow.addWidget(self.lbl_file)
        self.lbl_counts = QtWidgets.QLabel("DATA accepted: 0   rejected: 0")
        frow.addWidget(self.lbl_counts)
        frow.addStretch(1)
        self.note_edit = QtWidgets.QLineEdit()
        self.note_edit.setPlaceholderText("Note for the metadata file (e.g. 'settled by eye at 100 s')")
        self.note_edit.setMinimumWidth(320)
        self.note_edit.returnPressed.connect(self.add_note)
        frow.addWidget(self.note_edit)
        self.note_button = QtWidgets.QPushButton("Add note")
        self.note_button.clicked.connect(self.add_note)
        frow.addWidget(self.note_button)
        self.new_file_button = QtWidgets.QPushButton("New data file")
        self.new_file_button.clicked.connect(self.new_data_file)
        frow.addWidget(self.new_file_button)
        layout.addWidget(file_box)

        # --- live temperatures ---
        row_t = QtWidgets.QHBoxLayout()
        self.lbl_temps = {}
        self.show_checks = {}
        for col in TEMPERATURE_COLUMNS:
            color, _ = TRACE_STYLE[col]
            check = QtWidgets.QCheckBox()
            check.setChecked(True)
            check.toggled.connect(self.update_trace_visibility)
            self.show_checks[col] = check
            lbl = QtWidgets.QLabel(f"{SENSOR_LABELS[col]}: --- C")
            css = color if isinstance(color, str) else "rgb({},{},{})".format(*color)
            lbl.setStyleSheet(f"color: {'black' if css == 'k' else css}; font-weight: bold;")
            self.lbl_temps[col] = lbl
            row_t.addWidget(check)
            row_t.addWidget(lbl)
            row_t.addSpacing(12)
        row_t.addStretch(1)
        layout.addLayout(row_t)

        # --- live state reported by the Arduino ---
        row_s = QtWidgets.QHBoxLayout()
        self.lbl_state = QtWidgets.QLabel("State: ---")
        self.lbl_mode = QtWidgets.QLabel("Mode: ---")
        self.lbl_step = QtWidgets.QLabel("Step: ---")
        self.lbl_wave = QtWidgets.QLabel("Ta/Tb/period: ---")
        self.lbl_tset = QtWidgets.QLabel("Tset: ---")
        self.lbl_pwm = QtWidgets.QLabel("PWM: ---")
        self.lbl_parts = QtWidgets.QLabel("ff / controller: ---")
        self.lbl_time = QtWidgets.QLabel("Time: --- s")
        for lbl in (self.lbl_state, self.lbl_mode, self.lbl_step, self.lbl_wave,
                    self.lbl_tset, self.lbl_pwm, self.lbl_parts, self.lbl_time):
            row_s.addWidget(lbl)
            row_s.addSpacing(10)
        row_s.addStretch(1)
        layout.addLayout(row_s)

        row_m = QtWidgets.QHBoxLayout()
        self.lbl_message = QtWidgets.QLabel("Arduino: ---")
        row_m.addWidget(self.lbl_message)
        row_m.addStretch(1)
        row_m.addWidget(QtWidgets.QLabel(
            "Plots: wheel = zoom (over an axis: that axis only), drag = pan, right-click = menu"))
        self.show_all_button = QtWidgets.QPushButton("Show all")
        self.show_all_button.clicked.connect(self.show_all)
        row_m.addWidget(self.show_all_button)
        layout.addLayout(row_m)

        # --- temperature plot ---
        self.temp_plot = pg.PlotWidget()
        self.temp_plot.setLabel("left", "Temperature (C)")
        self.temp_plot.addLegend(offset=(10, 10))
        self.temp_curves = {}
        for col in TEMPERATURE_COLUMNS:
            color, width = TRACE_STYLE[col]
            self.temp_curves[col] = self.temp_plot.plot(
                pen=pg.mkPen(color, width=width), connect="finite", name=SENSOR_LABELS[col])
        self.setpoint_curve = self.temp_plot.plot(
            pen=pg.mkPen("g", width=2, style=QtCore.Qt.DashLine), connect="finite", name="Tset")
        self.temp_plot.getViewBox().setLimits(minYRange=TEMP_MIN_Y_SPAN)
        layout.addWidget(self.temp_plot, stretch=3)

        # --- PWM plot ---
        self.cmd_plot = pg.PlotWidget()
        self.cmd_plot.setLabel("left", "Signed PWM")
        self.cmd_plot.addLegend(offset=(10, 10))
        self.pwm_curve = self.cmd_plot.plot(pen=pg.mkPen("k", width=2),
                                            name="applied (cooling negative)")
        self.ff_curve = self.cmd_plot.plot(pen=pg.mkPen("b", width=1.5), connect="finite",
                                           name="feedforward (before x r)")
        self.ctrl_curve = self.cmd_plot.plot(pen=pg.mkPen((255, 140, 0), width=1.5),
                                             connect="finite", name="controller Kp e + Ki q")
        self.cmd_plot.addLine(y=0, pen=pg.mkPen((160, 160, 160), width=1)).setZValue(-10)
        layout.addWidget(self.cmd_plot, stretch=1)

        # --- tracking error plot ---
        self.error_plot = pg.PlotWidget()
        self.error_plot.setLabel("bottom", "Arduino time (s)")
        self.error_plot.setLabel("left", "Tset - T_TEC (C)")
        self.error_plot.addLine(y=0, pen=pg.mkPen("k", width=1)).setZValue(-10)
        self.error_curve = self.error_plot.plot(pen=pg.mkPen("m", width=2), connect="finite")
        layout.addWidget(self.error_plot, stretch=1)

        self.plots = (self.temp_plot, self.cmd_plot, self.error_plot)
        for p in self.plots:
            p.showGrid(x=True, y=True, alpha=0.3)
            if p is not self.temp_plot:
                p.setXLink(self.temp_plot)

        self.resize(1400, 1050)

    # ---------------- settings ----------------
    def wanted_settings(self):
        """Settings shown in the GUI, as {command: value text}."""
        wanted = {"MODE": self.mode_combo.currentText(),
                  "FFADD": "1" if self.ffadd_check.isChecked() else "0"}
        for name, spin in self.spins.items():
            wanted[name] = format_command_value(spin.value())
        return wanted

    def changed_settings(self):
        return {k: v for k, v in self.wanted_settings().items() if self.sent.get(k) != v}

    def update_pending(self):
        n = len(self.changed_settings())
        self.lbl_pending.setText(f"{n} setting(s) not sent yet" if n else "")
        self.lbl_pending.setStyleSheet("color: darkorange; font-weight: bold;" if n else "")

    def send_settings(self):
        """Send the settings that differ from what the Arduino last accepted.

        MEAN and AMP are checked together on the Arduino (Ta - Tb and Ta + Tb must
        stay within the setpoint limits), so when both change the one that keeps
        the waveform inside the limits at every intermediate step is sent first."""
        changes = self.changed_settings()
        order = list(changes)
        if "MEAN" in changes and "AMP" in changes:
            new_amp = float(changes["AMP"])
            old_amp = float(self.sent.get("AMP", "inf"))
            first, second = ("AMP", "MEAN") if new_amp < old_amp else ("MEAN", "AMP")
            order = [k for k in order if k not in ("MEAN", "AMP")]
            order[0:0] = [first, second]
        for key in order:
            value = changes[key]
            self.send_line(f"MODE {value}" if key == "MODE" else f"SET {key} {value}")
            self.sent[key] = value
            if key != "MODE":
                self.pending_echo.setdefault(key, []).append([value, None])
        self.update_pending()
        return order

    def sync_settings(self):
        """After the Arduino (re)starts: send every setting, then ask for a CONFIG line."""
        self.sent = {}
        self.pending_echo = {}
        self.send_settings()
        self.send_line("GET")
        self.synced = True

    def start_run(self):
        self.send_settings()
        self.send_line("START")

    def stop_run(self):
        if self.runner is not None and self.runner.active:
            self.abort_schedule("Stop pressed")
        else:
            self.send_line("STOP")

    # ---------------- automatic schedule ----------------
    def load_schedule(self, path=None):
        if path is None:
            start = str(sched.REPO_ROOT / "schedules")
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load schedule", start,
                                                            "CSV files (*.csv)")
            if not path:
                return False
        path = Path(path)
        try:
            rows = sched.read_schedule(path)
        except (OSError, ValueError) as e:
            QtWidgets.QMessageBox.warning(self, "Schedule not loaded", str(e))
            return False
        self.schedule_rows, self.schedule_path = rows, path
        self.lbl_schedule.setText(f"{path.name}: {len(rows)} rows loaded")
        self.lbl_schedule.setStyleSheet("")
        self.start_schedule_button.setEnabled(True)
        print(f"[schedule] loaded {path}")
        for r in rows:
            print(f"   {r.describe()}")
        return True

    def start_schedule(self, confirm=True):
        if not self.schedule_rows or (self.runner is not None and self.runner.active):
            return False
        if self.last_time_ms is None:
            QtWidgets.QMessageBox.warning(self, "No data yet",
                                          "Wait for the first records from the Arduino.")
            return False
        if confirm:
            text = ("The schedule will drive the TEC automatically:\n\n"
                    + "\n".join(r.describe() for r in self.schedule_rows)
                    + "\n\nStay with the apparatus. The run stops on SAFETY, TIMEOUT, any "
                      "Arduino ERROR, Stop, or closing this window. Start now?")
            answer = QtWidgets.QMessageBox.question(self, "Start schedule", text)
            if answer != QtWidgets.QMessageBox.Yes:
                return False
        self.runner = sched.ScheduleRunner(self.schedule_rows, lambda: self.files.csv_path,
                                           DATA_DIR)
        self.files.event(self.last_time_ms, f"SCHEDULE START {self.schedule_path.name} "
                                            f"({len(self.schedule_rows)} rows)")
        for w in (self.send_button, self.start_button, self.load_schedule_button,
                  self.start_schedule_button, self.new_file_button):
            w.setEnabled(False)
        self.abort_schedule_button.setEnabled(True)
        self.handle_schedule_action(self.runner.start())
        return True

    def apply_schedule_row(self, row):
        prev = self.runner.rows[self.runner.index - 1] if self.runner.index > 0 else None
        if prev is not None and prev.kind == "baseline" and "TROOM" not in row.settings:
            room = self.runner.baseline_mean()
            if not math.isnan(room):
                self.spins["TROOM"].setValue(room)
                self.files.event(self.last_time_ms,
                                 f"SCHEDULE T_room set to {room:.2f} C (mean of the six "
                                 f"channels during {prev.label})")
        for name, value in row.settings.items():
            self.spins[name].setValue(value)
        test_limit = None
        if row.is_test:
            channel, setting, _ = sched.TEST_MODES[row.mode]
            test_limit = round(self.last_rec[channel] + row.limit_rise, 2)
            self.restore_limits()
            self.spins[setting].setValue(test_limit)
            self.files.event(self.last_time_ms, f"SCHEDULE {row.label}: {setting} lowered to "
                                                f"{test_limit:.2f} C for the shutdown test")
        else:
            self.restore_limits()
        if row.mode == "OFF":
            self.send_line("STOP")
            self.send_settings()                    # e.g. limits restored after a test
        else:
            self.mode_combo.setCurrentText("MANUAL" if row.is_test else row.mode)
            self.send_settings()
            if self.last_state != "RUN":
                self.send_line("START")
        self.runner.row_applied(self.last_time_ms / 1000.0, test_limit)

    def restore_limits(self):
        """Back to the normal safety limits (56 C for T1, 60 C for the TEC); the
        caller sends them with send_settings()."""
        for _, setting, normal in sched.TEST_MODES.values():
            self.spins[setting].setValue(normal)

    def handle_schedule_action(self, action):
        if action is None:
            return
        if action.kind == "status":
            self.lbl_schedule.setText(f"[{self.runner.index + 1}/{len(self.runner.rows)}] {action.text}")
            return
        if action.kind == "note":
            self.files.event(self.last_time_ms, f"SCHEDULE {action.text}")
            print(f"[schedule] {action.text}")
            self.lbl_message.setText(f"Glitch: {action.text}")
            self.lbl_message.setStyleSheet("color: darkorange; font-weight: bold;")
            return
        self.files.event(self.last_time_ms, f"SCHEDULE {action.text}")
        print(f"[schedule] {action.text}")
        if action.kind == "apply":
            self.lbl_schedule.setText(action.text)
            self.lbl_schedule.setStyleSheet("color: darkgreen; font-weight: bold;")
            self.apply_schedule_row(action.row)
        elif action.kind == "finished":
            self.send_line("STOP")
            self.restore_limits()
            self.send_settings()
            self.lbl_schedule.setText(f"Schedule complete; intervals added to "
                                      f"{DATA_DIR / 'intervals_auto_*.csv'}")
            self.lbl_schedule.setStyleSheet("color: darkgreen; font-weight: bold;")
            self._schedule_buttons_idle()
            if self.autorun and self.autorun.get("analyse"):
                self.launch_analysis()
        elif action.kind == "abort":
            self.abort_schedule(action.text)

    def abort_schedule(self, reason):
        if self.runner is not None:
            self.runner.abort()
        self.send_line("STOP")
        self.restore_limits()
        self.send_settings()
        self.files.event(self.last_time_ms, f"SCHEDULE ABORTED: {reason}")
        print(f"[schedule] ABORTED: {reason}")
        self.lbl_schedule.setText(f"Schedule stopped: {reason}")
        self.lbl_schedule.setStyleSheet("color: white; background-color: red; font-weight: bold;")
        self._schedule_buttons_idle()

    def launch_analysis(self):
        """After an automatic run: write this run's segments to interval files and
        start the matching analysis script for each kind in the background (log and
        figures in data/module_08/analysis/ and docs/figures/module_08/)."""
        results = self.runner.results
        run_name = self.files.csv_path.stem
        out = DATA_DIR / "analysis"
        out.mkdir(parents=True, exist_ok=True)
        baselines = [r for r in results if r["kind"] == "baseline"]
        scripts = {"control": "module_08_control_fit.py", "pilot": "module_08_angstrom_fit.py",
                   "production": "module_08_angstrom_fit.py"}
        launched = []
        for kind, script in scripts.items():
            rows = [r for r in results if r["kind"] == kind]
            if not rows:
                continue
            ivfile = out / f"intervals_{run_name}_{kind}.csv"
            with open(ivfile, "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(sched.INTERVAL_COLUMNS)
                for r in rows:
                    w.writerow([r["label"], r["file"], "", r["step"], f"{r['t_start']:.1f}",
                                f"{r['t_end']:.1f}", "no", "",
                                "auto run" + ("" if r["settled"] else "; NOT SETTLED")])
            fig_dir = (REPO_ROOT / "docs" / "figures" / "module_08" if DATA_DIR == MODULE8_DATA_DIR
                       else DATA_DIR / "figures")       # a self-test keeps its files out of the repo
            args = [sys.executable, str(REPO_ROOT / "python" / "analysis" / script),
                    "--intervals", str(ivfile), "--name", f"{run_name}_{kind}",
                    "--out-dir", str(out), "--fig-dir", str(fig_dir)]
            if kind != "control" and baselines:
                b = baselines[-1]
                args += ["--baseline", f"{b['file']}:{b['t_start']:.1f},{b['t_end']:.1f}"]
            log = out / f"{run_name}_{kind}.log"
            with open(log, "w") as lf:
                self.analysis_procs.append(subprocess.Popen(args, cwd=REPO_ROOT, stdout=lf,
                                                            stderr=subprocess.STDOUT))
            launched.append(log.name)
            self.files.event(self.last_time_ms, f"ANALYSIS started: {' '.join(args[1:])} > {log.name}")
        if launched:
            self.lbl_schedule.setText(self.lbl_schedule.text() + f"; analysis running: {', '.join(launched)}")
            print(f"[analysis] started; logs in {out}: {', '.join(launched)}")
        return launched

    def maybe_autostart(self):
        """module_08_autorun.py: once records arrive and the settings are synced, load
        the schedule and start it (with the confirmation dialog unless --yes)."""
        a = self.autorun
        if not a or a.get("started") or self.n_accepted < 3 or self.cmd_queue:
            return
        a["started"] = True
        if self.load_schedule(a["schedule"]):
            self.start_schedule(confirm=not a.get("yes"))

    def _schedule_buttons_idle(self):
        for w in (self.send_button, self.start_button, self.load_schedule_button):
            w.setEnabled(True)
        self.start_schedule_button.setEnabled(self.schedule_rows is not None)
        self.abort_schedule_button.setEnabled(False)

    def send_line(self, cmd: str, echo: bool = True):
        """HB and STOP go out at once (STOP also discards queued commands, so nothing
        queued, such as a START, can follow it); everything else is queued."""
        if cmd == "STOP":
            if self.cmd_queue:
                print(f"[queue] STOP: discarded {self.cmd_queue}")
            self.cmd_queue.clear()
            # STOP goes out at once and always succeeds: treat the Arduino as stopped
            # from now on, so a START decided before the next record is not skipped.
            self.last_state = "OFF"
        if cmd in ("HB", "STOP") or COMMAND_SPACING_MS <= 0:
            self._write(cmd, echo)
        else:
            self.cmd_queue.append((cmd, echo))

    def send_next_command(self):
        if self.cmd_queue:
            cmd, echo = self.cmd_queue.pop(0)
            self._write(cmd, echo)

    def _write(self, cmd: str, echo: bool):
        if cmd.startswith("SET "):
            parts = cmd.split()
            for entry in self.pending_echo.get(parts[1], []):
                if entry[1] is None:
                    entry[1] = time.monotonic()
                    break
        try:
            self.link.write_line(cmd)
            if echo:
                print(f"[sent] {cmd}")
        except Exception as e:
            print(f"[send error] {e}")

    def check_set_echo(self, text):
        """The Arduino echoes every accepted setting ("SET PERIOD 200.0000 (was ...)").
        A value different from the one sent means a corrupted command: treat it like
        an ERROR (resend everything next time; stop a running schedule)."""
        m = re.match(r"SET (\w+) ([-\d.]+)", text)
        if not m or not self.pending_echo.get(m.group(1)):
            return
        expected = self.pending_echo[m.group(1)].pop(0)[0]
        try:
            mismatch = abs(float(m.group(2)) - float(expected)) > 1e-3
        except ValueError:
            return
        if mismatch:
            self.on_line(f"ERROR (GUI check): Arduino accepted SET {m.group(1)} {m.group(2)}, "
                         f"but {expected} was sent")

    def check_missing_echoes(self):
        """A SET written more than ECHO_TIMEOUT_S ago without an echo or an ERROR was
        lost on the way: treat it like an ERROR."""
        now = time.monotonic()
        for key, entries in self.pending_echo.items():
            for value, written in entries:
                if written is not None and now - written > ECHO_TIMEOUT_S:
                    self.on_line(f"ERROR (GUI check): no reply to SET {key} {value} "
                                 f"within {ECHO_TIMEOUT_S:g} s (command lost)")
                    return

    # ---------------- data file, notes, counters ----------------
    def update_file_label(self):
        self.lbl_file.setText(f"File: {self.files.csv_path.name}  ({self.files.rows} rows)")

    def update_counts(self):
        self.lbl_counts.setText(f"DATA accepted: {self.n_accepted}   rejected: {self.n_rejected}")

    def add_note(self):
        text = self.note_edit.text().strip()
        if not text:
            return
        self.files.event(self.last_time_ms, f"NOTE: {text}")
        self.note_edit.clear()
        print(f"[note] {text}")

    def new_data_file(self):
        if self.last_state == "RUN":
            QtWidgets.QMessageBox.warning(self, "Run in progress",
                                          "Stop the run before starting a new data file.")
            return
        old = self.files
        self.files = RunFiles(DATA_DIR)
        old.event(self.last_time_ms, f"closed; continued in {self.files.csv_path.name}")
        old.close()
        self.files.event(self.last_time_ms, f"continues {old.csv_path.name}")
        self.send_line("GET")          # current settings into the new metadata file
        self.update_file_label()
        print(f"[data] saved {old.csv_path.name} ({old.rows} rows); "
              f"now recording to {self.files.csv_path}")

    # ---------------- one line received ----------------
    @QtCore.Slot(str)
    def on_line(self, line: str):
        kind, raw_fields, rec = classify_line(line)
        if kind == "header":
            return
        if kind == "event":
            time_ms, step, text = rec
            self.files.event(time_ms, f"step {step}: {text}")
            self.lbl_message.setText(f"Arduino: step {step}: {text}")
            print(f"[event] step {step}: {text}")
            self.check_set_echo(text)
            return
        if kind == "config":
            time_ms, step, text = rec
            self.files.event(time_ms, f"step {step}: CONFIG {text}")
            print(f"[config] step {step}: {text}")
            return
        if kind == "other":
            print(f"[arduino] {line}")
            if line.startswith("ERROR"):
                # The Arduino kept its previous value: forget what was sent, so the
                # next "Send settings" sends everything again.
                self.sent = {}
                self.pending_echo = {}
                self.update_pending()
                self.files.event(self.last_time_ms, f"Arduino: {line}")
                self.lbl_message.setText(f"Arduino: {line}")
                self.lbl_message.setStyleSheet("color: red; font-weight: bold;")
                if self.runner is not None and self.runner.active:
                    self.abort_schedule(f"Arduino {line}")
            return
        if kind == "rejected":
            self.n_rejected += 1
            self.update_counts()
            print(f"[rejected] {line}")
            return

        # ---- an accepted DATA record: save it exactly as received ----
        self.n_accepted += 1
        self.last_processed = time.monotonic()
        self.files.write_record(raw_fields)
        self.update_counts()
        self.update_file_label()

        # Heartbeat: keeps a run going only while this program receives records.
        self.send_line("HB", echo=False)
        self.check_missing_echoes()

        t_ms = rec["time_ms"]
        if self.last_time_ms is not None and t_ms < self.last_time_ms:
            print(f"[warning] Arduino time went back ({self.last_time_ms} -> {t_ms} ms): "
                  "the Arduino restarted")
            self.files.event(t_ms, f"Arduino restarted (time went back from {self.last_time_ms} ms)")
            self.synced = False
            if self.ts:
                self.t_offset = self.ts[-1] + 1.0 - t_ms / 1000.0
        self.last_time_ms = t_ms
        state = rec["run_state"]
        self.last_rec = rec
        self.last_state = state
        if not self.synced:
            self.sync_settings()

        # plots
        self.ts.append(t_ms / 1000.0 + self.t_offset)
        for col in TEMPERATURE_COLUMNS:
            self.temps[col].append(rec[col])
        self.setpoints.append(rec["setpoint_C"])
        self.signed.append(rec["signed_pwm"])
        self.feedforward.append(rec["feedforward_pwm"])
        self.controller.append(rec["controller_pwm"])
        self.errors.append(rec["setpoint_C"] - rec["TEC_C"])

        # labels
        for col in TEMPERATURE_COLUMNS:
            value = rec[col]
            text = "nan (check sensor)" if math.isnan(value) else f"{value:.2f} C"
            self.lbl_temps[col].setText(f"{SENSOR_LABELS[col]}: {text}")
        self.lbl_state.setText(f"State: {state}")
        self.lbl_state.setStyleSheet(STATE_STYLE.get(state, ""))
        self.lbl_mode.setText(f"Mode: {rec['control_mode']}")
        self.lbl_step.setText(f"Step: {rec['step']}")
        self.lbl_wave.setText(f"Ta/Tb/period: {rec['Ta_C']:g} / {rec['Tb_C']:g} C / "
                              f"{rec['period_s']:g} s")
        self.lbl_tset.setText(f"Tset: {rec['setpoint_C']:.2f} C")
        self.lbl_pwm.setText(f"PWM: {rec['signed_pwm']:+d}")
        self.lbl_parts.setText(f"ff / controller: {rec['feedforward_pwm']:.1f} / "
                               f"{rec['controller_pwm']:.1f}  (r = {rec['heating_to_cooling_ratio']:g})")
        self.lbl_time.setText(f"Time: {t_ms / 1000.0:.2f} s")
        schedule_on = self.runner is not None and self.runner.active
        self.new_file_button.setEnabled(state != "RUN" and not schedule_on)
        if schedule_on:
            self.handle_schedule_action(self.runner.on_record(rec))
        elif self.autorun:
            self.maybe_autostart()

    @QtCore.Slot(str)
    def on_serial_error(self, msg: str):
        print(f"[serial error] {msg}")
        self.files.event(self.last_time_ms, f"serial error: {msg}")
        self.lbl_message.setText(f"Serial error: {msg} (the Arduino stops the run after 3 s)")
        self.lbl_message.setStyleSheet("color: red; font-weight: bold;")

    # ---------------- plots ----------------
    def update_trace_visibility(self):
        for col, check in self.show_checks.items():
            self.temp_curves[col].setVisible(check.isChecked())

    def update_plots(self):
        if not self.ts:
            return
        xs = self.ts
        for col in TEMPERATURE_COLUMNS:
            self.temp_curves[col].setData(xs, self.temps[col])
        self.setpoint_curve.setData(xs, self.setpoints)
        self.pwm_curve.setData(xs, self.signed)
        self.ff_curve.setData(xs, self.feedforward)
        self.ctrl_curve.setData(xs, self.controller)
        self.error_curve.setData(xs, self.errors)

    def show_all(self):
        for p in self.plots:
            p.enableAutoRange()

    # ---------------- clean shutdown ----------------
    def closeEvent(self, event):
        stop_keep_awake(self.keep_awake)
        self.cmd_timer.stop()
        if self.runner is not None and self.runner.active:
            self.runner.abort()
            self.files.event(self.last_time_ms, "SCHEDULE ABORTED: GUI closed")
        self.send_line("STOP")
        self.timer.stop()
        self.reader.stop()
        self.link.close()
        self.files.event(self.last_time_ms, "GUI closed (STOP sent)")
        self.files.close()
        print(f"[data] saved {self.files.csv_path} ({self.files.rows} rows)")
        super().closeEvent(event)


# ------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------
def main():
    app = QtWidgets.QApplication(sys.argv)
    pg.setConfigOptions(background="w", foreground="k")
    try:
        window = AngstromGUI()
    except serial.SerialException as e:
        from serial.tools import list_ports
        ports = [p.device for p in list_ports.comports()]
        print(f"\nCould not open SERIAL_PORT = {SERIAL_PORT!r}:\n  {e}")
        print("Serial ports on this computer: "
              + (", ".join(ports) if ports else "none (is the Arduino plugged in?)"))
        print("Set SERIAL_PORT at the top of this file to the Arduino's port "
              "(a Mac name starts with /dev/cu.), and close the Arduino Serial Monitor.")
        sys.exit(1)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
