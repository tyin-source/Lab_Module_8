# Module 8 data

Raw files written by `python/module_08_angstrom_gui.py`. Never edit or overwrite a raw
file; exclusions, corrections, and fit results belong in the analysis code and its
outputs.

- `angstrom_<YYYYmmdd-HHMMSS>.csv`: every accepted Arduino `DATA` record, exactly as
  sent. A new file starts when the GUI starts and with the "New data file" button
  (only while no run is active).
- `angstrom_<YYYYmmdd-HHMMSS>_meta.txt`: the matching metadata. It holds the serial
  settings, code files and Git commit, rod geometry, and the channel map with sensor
  positions. It also lists time-stamped events: every Arduino `EVENT` line (start,
  stop, mode, setting or period change, safety shutdown, timeout, each with its
  `step`), `CONFIG` replies, `ERROR` replies, Arduino restarts, and typed notes.

## Record format (S18.2)

One serial line per ~1 s average, header printed once at start-up:

```text
DATA,time_ms,run_state,control_mode,step,Ta_C,Tb_C,period_s,setpoint_C,feedforward_pwm,TEC_C,T1_C,T2_C,T3_C,T4_C,T5_C,controller_pwm,heating_to_cooling_ratio,signed_pwm
```

The CSV has the same 18 columns, without the `DATA` tag. The number and order of
fields never change. A field that does not apply in the current mode holds the
placeholder `nan`.

| Column | Meaning | Unit | Placeholder |
| --- | --- | --- | --- |
| `time_ms` | Arduino clock since its last reset | ms | never |
| `run_state` | `OFF`, `RUN`, `SAFETY` (a limit stopped the run), `TIMEOUT` (no heartbeat for 3 s) | — | never |
| `control_mode` | `MANUAL`, `PWM_SINE`, `FF`, `P`, `PI` | — | never |
| `step` | run-segment counter, +1 at every `EVENT` (start, stop, any setting, mode, or period change, safety, timeout) | — | never |
| `Ta_C` | requested mean TEC temperature $T_a$ | °C | `nan` in MANUAL and PWM_SINE |
| `Tb_C` | requested amplitude $T_b$ | °C | `nan` in MANUAL and PWM_SINE |
| `period_s` | forcing period | s | `nan` in MANUAL |
| `setpoint_C` | requested TEC temperature $T_a+T_b\sin\phi$ at this record | °C | `nan` unless FF, P, or PI is running |
| `feedforward_pwm` | $(T_\mathrm{set}-T_\mathrm{room})/\chi_H$, before the cooling scaling | PWM | `nan` unless FF is running, or P/PI with FFADD 1 |
| `TEC_C` | thermistor embedded in the TEC base (A0): actuator, feedback, and safety sensor | °C | `nan` if unplugged |
| `T1_C` | rod thermistor at the first hole, $x=0$ (A1): the Angstrom boundary and rod-base safety sensor | °C | `nan` if unplugged |
| `T2_C`…`T5_C` | downstream rod thermistors (A2–A5); positions in the metadata file | °C | `nan` if unplugged |
| `controller_pwm` | $K_p e + K_i q$, $e=T_\mathrm{set}-T_\mathrm{TEC}$, before the cooling scaling | PWM | `nan` unless P or PI is running |
| `heating_to_cooling_ratio` | $r=\chi_H/\chi_C$ applied to cooling commands | — | `nan` unless FF, P, or PI is running |
| `signed_pwm` | PWM applied to the TEC after scaling, clamping, and safety; heating positive, cooling negative | −255…255 | 0 when not running |

Temperatures are raw readings: no calibration corrections are applied in the file.

Settings that are not in the record are in the metadata file: $K_p$ and $K_i$
(PWM/°C, PWM/(°C s)), $u_\max$, $\chi_H$, $T_\mathrm{room}$, the PWM_SINE offset $u_a$
and amplitude $u_b$, the MANUAL PWM, and FFADD. They appear in the `START` event, in each
`SET` event, and in `CONFIG` lines.

## Analysis inputs and outputs

- `intervals_control.csv`, `intervals_production.csv`: one row per analysis interval
  (`label,file,part,step,t_start_s,t_end_s,settle,exclude,notes`, optionally
  `baseline` = `FILE:START,END` for the room temperature of that interval's run, with
  `:extrapolate` when the rod was still relaxing). They record which
  records each result and figure uses (A5 Part 1, figure sources). Rows whose label
  starts with `#` are ignored. `part` is needed only if the Arduino restarted within
  the file.
- `intervals_auto_<kind>.csv` (kind = baseline, control, pilot, production): written by
  the GUI's automatic schedule, one row per completed segment, with the settled time
  range. Rows marked NOT SETTLED reached their time limit without passing the
  settling test.
- `analysis/`: tables written by the analysis scripts (`control_<name>.csv`,
  `angstrom_<name>.csv`, `angstrom_<name>_sensors.csv`) and by the numerical model
  `python/analysis/module_08_simulation.py`:
  - `sim_parameters.csv`: calibrated model parameters, with their sources;
  - `sim_s18_3_course.csv`, `sim_s18_3_pwm100_compare.csv`: S18.3;
  - `sim_s18_4.csv`: S18.4;
  - `sim_control.csv`: the focused P/PI simulation;
  - `sim_rod_<label>.csv`: the integrated comparison.

  These are derived results; the raw files are never changed.

## Which file is which

| File | Role |
| --- | --- |
| `angstrom_20261116-201051.csv` | GUI opened and closed after 20 s; no records (not used) |
| `angstrom_20261116-204536.csv` | S18.2 GUI check: six channels live with the TEC off, 1226 records |
| `angstrom_20261116-211726.csv` | Automatic run of `schedules/s18_3_4_raw_pwm.csv`: S18.1 checks and shutdown tests, room baseline (step 28), S18.3 at 300/200/100 s (steps 32–34), S18.4 (step 35) |
| `angstrom_20261117-000509.csv` | First attempt at `schedules/s18_5_7_control.csv`: room baseline (step 17, 21.91 °C), then the feedforward run (step 20) stopped after 5 min by the heartbeat timeout (macOS App Nap with the display asleep slowed the GUI); idle afterwards. Not used for S18.6 |
| `angstrom_20261117-083725.csv` | Repeat of `schedules/s18_5_7_control.csv`: room baseline (step 17, 23.10 °C), S18.6 feedforward (step 20), S18.7 P (step 22) and PI $K_p$ 6, $K_i$ 0.04 at 100/200/300 s (steps 24–26) |
| `angstrom_20261117-114802.csv` | `schedules/s18_7_tuning.csv`: PI tuning at 100 s ($K_p$ 15 and 30, $T_b$ 16 and 10 °C) and the selected $K_p$ 30, $K_i$ 0.3 at 200 and 300 s (steps 20–26) |
| `angstrom_20261117-132324.csv` | `schedules/s18_9_10_pilots.csv`: baseline (step 17; the rod was still cooling from the tuning run, room extrapolated to 23.34–23.50 °C), pilots at 100 s (step 22) and 200 s (step 23); both used as production intervals |
| `angstrom_20261117-161210.csv` | `schedules/production.csv`: room baseline (step 17, 22.81–22.96 °C), production at 300 s (step 23), 400 s (step 24), 500 s (step 25), then 25 min of cool-down with the TEC off (step 26, room drift) |

The files `analysis/*angstrom_20261117-161210_production*` and the matching figures
were written by the GUI's automatic analysis at the end of the run, with the live
intervals and the code of that moment. The authoritative five-period results are
`analysis/angstrom_production*.csv`, from `intervals_production.csv`.

`course_simulator/`: records from the course's six-sensor simulator, run headlessly by
`tools/course_simulator_driver.js`:
- `course_sim_<label>.csv`, one row per simulated second: setpoint, duty and PWM, TEC
  plate, T0 (feedback), and the rod at our five sensor positions;
- `course_sim_meta.json`: page URL, fetch time, SHA-256, settings.
- `course_lump_<label>.csv`, `course_lump_meta.json`: the one-lump mode (S18.3, S18.4),
  run by `tools/course_simulator_lumped.js`, one row per 0.1 s.

These are simulation outputs, not measurements.

`course_simulator_replica/`: the same cases from the Python replica of the course page
(`python/analysis/course_simulator_replica.py`):
- `course_sim_<label>.csv` and `course_lump_<label>.csv`, with the same columns as the
  browser records;
- `<label>_page.csv`, in the page's Save CSV format;
- `replica_check.txt`, the agreement with the browser runs.

These are simulation outputs.
