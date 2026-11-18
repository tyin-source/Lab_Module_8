# Phys 39 Module 8 - shared definitions for the Angstrom GUI and analysis
#
#   - repository paths for data and figures
#   - the DATA record fields sent by arduino/module_08_angstrom_controller
#   - the sensor map of the oscillatory (Angstrom) configuration
#   - reading one raw run CSV written by python/module_08_angstrom_gui.py

import csv
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE8_DATA_DIR = REPO_ROOT / "data" / "module_08"
MODULE8_FIG_DIR = REPO_ROOT / "docs" / "figures" / "module_08"

# Fields of one Arduino DATA record, in order (S18.2). These are also the CSV columns.
DATA_COLUMNS = [
    "time_ms", "run_state", "control_mode", "step", "Ta_C", "Tb_C", "period_s",
    "setpoint_C", "feedforward_pwm", "TEC_C", "T1_C", "T2_C", "T3_C", "T4_C", "T5_C",
    "controller_pwm", "heating_to_cooling_ratio", "signed_pwm",
]
TEXT_COLUMNS = ["run_state", "control_mode"]
INT_COLUMNS = ["time_ms", "step", "signed_pwm"]
TEMPERATURE_COLUMNS = ["TEC_C", "T1_C", "T2_C", "T3_C", "T4_C", "T5_C"]
ROD_COLUMNS = TEMPERATURE_COLUMNS[1:]
SENSOR_LABELS = {"TEC_C": "TEC", "T1_C": "T1", "T2_C": "T2",
                 "T3_C": "T3", "T4_C": "T4", "T5_C": "T5"}
ARDUINO_INPUTS = {"TEC_C": "A0", "T1_C": "A1", "T2_C": "A2",
                  "T3_C": "A3", "T4_C": "A4", "T5_C": "A5"}

RUN_STATES = ("OFF", "RUN", "SAFETY", "TIMEOUT")
CONTROL_MODES = ("MANUAL", "PWM_SINE", "FF", "P", "PI")

# Oscillatory sensor map (S18.1 starting layout): T1 in the first hole (x = 0) and
# T2-T5 about 1 inch apart, using the existing holes at 0, 2.4, 5.0, 7.6, 10.2 cm
# (Module 7 hole positions). Holes are numbered 1-9 from the TEC end: Module 8 uses
# holes 1, 3, 5, 6, 7 for T1-T5 (Module 7 used holes 1, 5, 7, 8, 9). Update to the MEASURED positions before a run; the GUI
# writes these values to every metadata file. The Module 7 full-length map
# (0, 5.0, 10.2, 15.3, 20.4 cm) belongs to the Module 7 data only.
SENSOR_POSITIONS_CM = {"T1_C": 0.0, "T2_C": 2.4, "T3_C": 5.0, "T4_C": 7.6, "T5_C": 10.2}
POSITION_UNCERTAINTY_MM = 1.3

# Rod geometry measured in Module 7 (module_07_rod_instrument.md, section 1).
ROD_DIAMETER_M = 12.70e-3          # +/- 0.02 mm
ROD_OVERALL_LENGTH_M = 28.66e-2    # +/- 0.05 cm
ROD_ANALYSIS_LENGTH_M = 26.90e-2   # first hole to tip, +/- 0.05 cm


def read_run_csv(path: Path):
    """Read one raw run CSV into numpy arrays (nan for placeholders).

    run_state and control_mode are returned as string arrays; every other column
    as float (time_ms, step, and signed_pwm are whole numbers).
    """
    cols = {name: [] for name in DATA_COLUMNS}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        missing = [c for c in DATA_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        for row in reader:
            for name in DATA_COLUMNS:
                cols[name].append(row[name])
    out = {}
    for name, values in cols.items():
        if name in TEXT_COLUMNS:
            out[name] = np.array(values)
        else:
            out[name] = np.array([float(v) for v in values])
    out["time_s"] = out["time_ms"] / 1000.0
    return out
