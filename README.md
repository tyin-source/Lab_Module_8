# Phys 39 Lab: The Angstrom Method (Module 8)

Module 8 drives the base of the aluminum rod with a sinusoidal TEC temperature and
records the TEC and five rod thermistors. The oscillation shrinks and lags along the
rod; its amplitude decay and phase lag separate axial diffusion ($\kappa$) from side
heat loss ($\nu$). This repository extends the Module 7 rod instrument
(`tyin-source/Lab_Module_7`), which stays unchanged as the working checkpoint.

## Contents

| Path | What it is |
| --- | --- |
| `arduino/module_08_angstrom_controller/module_08_angstrom_controller.ino` | Six-thermistor acquisition plus the drive: sinusoidal setpoint, raw sinusoidal PWM, feedforward, P, and PI, all computed on the Arduino; 18-field `DATA` record; `EVENT` line at every transition; TEC and rod-base (T1) safety limits; 3 s heartbeat timeout |
| `python/module_08_angstrom_gui.py` | Sends settings, START/STOP, and the heartbeat; live display of all six temperatures, setpoint, mode, and PWM; raw CSV and metadata file per run |
| `python/module_08_autorun.py` | One-command run: port detection, GUI, schedule start, analysis afterwards |
| `python/module_08_schedule.py`, `schedules/` | Automatic schedule: the GUI runs the rows of a schedule CSV in order, ends each segment once the live settling test passes and the requested cycles are recorded, and writes the segments to `data/module_08/intervals_auto_<kind>.csv` (see `schedules/README.md`) |
| `python/module_08_common.py` | Record fields, sensor map, rod geometry, raw-CSV reader |
| `python/module_08_gui_selftest.py` | Off-screen test of the GUI against a simulated Arduino |
| `tests/firmware_host/` | Compiles the sketch on a computer with a stand-in Arduino core and tests it against a simulated TEC and rod |
| `python/analysis/angstrom_lib.py` | Shared analysis: common time fit, settling test, spatial fits, $\kappa$ and $\nu$, finite-rod and Biot checks, mean-profile fit, predictions, characteristic lengths |
| `python/analysis/module_08_control_fit.py` | S18.3–S18.8: TEC and T1 fits against the PWM or setpoint; lag, residual amplitude, amplitude ratio, mean and RMS tracking error, saturation; effective $\tau$ and $\chi$ with the rod mounted |
| `python/analysis/module_08_angstrom_fit.py` | S18.9–S18.10 and Stage 6: settling, per-sensor fits, $q$, $q'$, $\kappa$, $\nu$, model checks, mean profile, theory overlay, frequency consistency, characteristic-length plot |
| `python/analysis/module_08_analysis_selftest.py` | Runs both analysis scripts on synthetic data with known parameters |
| `data/module_08/intervals_*.csv` | The analysis intervals (file, step, times) behind each result and figure |
| `data/module_08/` | Raw runs; the record format is in its README |
| `docs/module_notes/module_08_angstrom.md` | Parameters carried forward, sensor map, safety limits, S18 evidence table, predictions, run log |
| `docs/module_notes/module_08_modeling_app.md` | Formative modeling-app note |
| `docs/figures/module_08/` | Figures |

## How the drive works

Every ~1 s the Arduino averages the six thermistors, checks the safety limits, updates
the drive, applies the PWM, and prints one record. The phase advances by
$2\pi\,\Delta t/P$, so a period change during a run continues the waveform without a
jump. With $T_\mathrm{set}=T_a+T_b\sin\phi$:

| Mode | Command $u$ (signed PWM, cooling negative) | Used in |
| --- | --- | --- |
| `MANUAL` | constant | baselines |
| `PWM_SINE` | $u_a + u_b\sin\phi$ ($u_a=0$ zero-centered, $u_a>u_b$ positive-only), no cooling compensation | S18.3, S18.4 |
| `FF` | $(T_\mathrm{set}-T_\mathrm{room})/\chi_H$ | S18.5, S18.6 |
| `P` | $K_p e$, $e=T_\mathrm{set}-T_\mathrm{TEC}$ | S18.7 |
| `PI` | $K_p e + K_i q$, $q=\sum e\,\Delta t$, conditional-integration anti-windup | S18.7 onward |

In `FF`, `P`, and `PI` a negative command is multiplied by $r=\chi_H/\chi_C$ (static
compensation of the weaker cooling), then clamped to $\pm u_\max$. The GUI sends `HB`
after every record it receives. If it stops (GUI frozen, USB unplugged), the Arduino
ends the run after 3 s. The GUI sends settings one every 60 ms (a burst would
overflow the Uno's 63-byte serial receive buffer and lose characters) and checks
every setting against the Arduino's echo. The course limits are in the sketch: rod-base T1 at or below
56 °C, TEC 60 °C, forcing period at least 100 s, requested setpoint 10–56 °C. Verify
the T1 limit before the first run (`docs/module_notes/module_08_angstrom.md`, section 3).

## How to run

```bash
pip install -r requirements.txt
python python/module_08_gui_selftest.py          # GUI check, no hardware needed
python python/module_08_angstrom_gui.py          # set SERIAL_PORT at the top first
```

Upload the sketch with the Arduino IDE and close the Serial Monitor before starting the
GUI. The serial speed is now 115200 baud.

Automatic data collection, one command (finds the port, opens the GUI, starts the
schedule, runs the analysis afterwards):

```bash
python python/module_08_autorun.py schedules/s18_3_4_raw_pwm.csv
```

Or in the GUI: "Load schedule..." → a file in `schedules/` → "Start schedule". Each completed segment is added to
`data/module_08/intervals_auto_<kind>.csv`, which the analysis scripts read directly.
Stay with the apparatus while a schedule runs.

Analysis (list the intervals in `data/module_08/intervals_control.csv` and
`intervals_production.csv`, or give one run and step on the command line):

```bash
python python/analysis/module_08_analysis_selftest.py     # check on synthetic data
python python/analysis/module_08_control_fit.py --intervals data/module_08/intervals_control.csv --room 22.0
python python/analysis/module_08_angstrom_fit.py --intervals data/module_08/intervals_production.csv --name production   # room baseline per interval (column baseline)
python python/analysis/module_08_angstrom_fit.py data/module_08/<run>.csv --step 7 --settle --room 22.0
python python/analysis/module_08_simulation.py           # TEC-rod model: S18.3/S18.4 tables, P/PI simulation
python python/analysis/module_08_simulation.py rod p300 --independent  # integrated theory-experiment-simulation comparison
node tools/course_simulator_driver.js data/module_08/course_simulator  # course six-sensor simulator, headless
node tools/course_simulator_lumped.js data/module_08/course_simulator  # its one-lump mode: S18.3, S18.4
node tools/course_simulator_screens.js docs/figures/module_08/course_ui data/module_08/course_simulator/ui_csv  # screenshots + page CSVs
python python/analysis/module_08_course_sim.py   # course simulator vs measurement vs our model
python python/analysis/course_simulator_replica.py   # Python replica of the course simulator: all cases, figures, check, comparison (no browser, ~1.5 min)
```

Figures go to `docs/figures/module_08/`, tables to `data/module_08/analysis/`.

One caution from the self-test: $\nu=\kappa(q^2-q'^2)$ is a small difference of two
large numbers. A finite-tip correction of 1–3 % in amplitude at the farthest
sensor (300–500 s with the 0–10.2 cm layout) moves the semi-infinite $\nu$ by
25–35 %, while $\kappa$ barely changes. The Angstrom script therefore also fits the
finite insulated-tip rod and, with `--model auto` (default), uses that fit for any
interval where the tip correction exceeds the measurement uncertainty. Each interval
reports which model it used.

Firmware logic test (optional, needs a C++ compiler):

```bash
g++ -std=c++17 -O2 -o fw_test tests/firmware_host/firmware_host_test.cpp && ./fw_test
```
