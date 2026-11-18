# Phys 39 Module 8 - run a schedule automatically
#
# Opens the Angstrom GUI on the Arduino's serial port (found automatically unless
# --port is given), waits for the first records, loads the schedule, and starts it.
# The usual confirmation dialog appears unless --yes is given. The GUI displays and
# saves every record as in a manual run. When the schedule completes, the matching
# analysis scripts start in the background (unless --no-analysis):
#   control rows                -> module_08_control_fit.py
#   pilot and production rows   -> module_08_angstrom_fit.py, with the run's last
#                                  baseline row as --baseline
# Logs, tables:  data/module_08/analysis/<run>_<kind>.log, *.csv
# Figures:       docs/figures/module_08/
#
# Usage (from the repository root):
#   python python/module_08_autorun.py schedules/s18_3_4_raw_pwm.csv
#   python python/module_08_autorun.py schedules/production.csv --port /dev/cu.usbmodem1101
#   python python/module_08_autorun.py --list-ports
#   python python/module_08_autorun.py            (no schedule: choose one in a dialog,
#                                                   e.g. when started with an IDE's Run button)
#
# Stay with the apparatus while it runs. Stop with the GUI's "Abort schedule" or
# "Stop" button, or by closing the window: each sends STOP at once.

import argparse
import sys
from pathlib import Path

import serial
from serial.tools import list_ports
from PySide6 import QtWidgets
import pyqtgraph as pg

import module_08_angstrom_gui as gui_module
import module_08_schedule as sched


def arduino_ports():
    """Serial ports that look like an Arduino (USB modem / ACM / Arduino in the name)."""
    found = []
    for p in list_ports.comports():
        text = " ".join(str(v) for v in (p.device, p.description, p.manufacturer) if v).lower()
        if any(k in text for k in ("arduino", "usbmodem", "ttyacm", "wchusbserial", "ch340")):
            found.append(p.device)
    return found


def main(argv=None):
    ap = argparse.ArgumentParser(description="Run a Module 8 schedule automatically")
    ap.add_argument("schedule", nargs="?", help="schedule CSV, e.g. schedules/s18_3_4_raw_pwm.csv")
    ap.add_argument("--port", help="serial port (default: the only Arduino-like port, "
                                   "else SERIAL_PORT in module_08_angstrom_gui.py)")
    ap.add_argument("--yes", action="store_true", help="start without the confirmation dialog")
    ap.add_argument("--no-analysis", action="store_true", help="do not start the analysis afterwards")
    ap.add_argument("--list-ports", action="store_true", help="list serial ports and exit")
    args = ap.parse_args(argv)

    if args.list_ports:
        for p in list_ports.comports():
            print(f"{p.device}\t{p.description}")
        return 0

    app = QtWidgets.QApplication(sys.argv[:1])
    pg.setConfigOptions(background="w", foreground="k")
    if not args.schedule:
        chosen, _ = QtWidgets.QFileDialog.getOpenFileName(
            None, "Choose a schedule to run", str(sched.REPO_ROOT / "schedules"),
            "Schedules (*.csv)")
        if not chosen:
            print("No schedule chosen (see schedules/README.md).")
            return 1
        args.schedule = chosen

    path = Path(args.schedule)
    if not path.is_absolute():
        path = (Path.cwd() / path) if (Path.cwd() / path).exists() else sched.REPO_ROOT / path
    try:
        rows = sched.read_schedule(path)          # check it before opening the port
    except (OSError, ValueError) as e:
        print(f"Schedule not usable: {e}")
        return 1
    print(f"Schedule {path.name}:")
    for r in rows:
        print(f"   {r.describe()}")

    port = args.port
    if port is None:
        ports = arduino_ports()
        if len(ports) == 1:
            port = ports[0]
        elif not ports:
            port = gui_module.SERIAL_PORT
            print(f"No Arduino-like port found; trying SERIAL_PORT = {port}")
        else:
            print(f"Several Arduino-like ports: {', '.join(ports)}. Choose one with --port.")
            return 1
    gui_module.SERIAL_PORT = port
    print(f"Serial port: {port} at {gui_module.BAUD_RATE} baud")

    try:
        window = gui_module.AngstromGUI()
    except serial.SerialException as e:
        print(f"Could not open {port}: {e}\nClose the Arduino Serial Monitor (only one program "
              "can use the port) or choose the port with --port (see --list-ports).")
        return 1
    window.autorun = {"schedule": path, "yes": args.yes, "analyse": not args.no_analysis}
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
