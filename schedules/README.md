# Run schedules

Quickest start (finds the Arduino's port, opens the GUI, starts the schedule, and runs
the analysis when it finishes):

```bash
python python/module_08_autorun.py schedules/s18_3_4_raw_pwm.csv
```

Add `--port <port>` if several ports are found (`--list-ports` lists them), `--yes` to skip
the confirmation dialog, `--no-analysis` to skip the analysis.

Each CSV here is a sequence of run segments for the GUI's **Automatic schedule**
(`python/module_08_angstrom_gui.py`, logic in `python/module_08_schedule.py`). Load one
with "Load schedule...", read the summary, and press "Start schedule". The GUI displays
and saves every record as usual. The schedule only drives the controls you would use
by hand.

| File | Steps | Notes |
| --- | --- | --- |
| `s18_3_4_raw_pwm.csv` | S18.1 automatic checks (channel spread < 1 °C, TEC and T1 shutdown tests), 20 min cool-down, room baseline, S18.3 zero-centered PWM at 300, 200, 100 s ($u_b$ = 50), S18.4 positive-only at 100 s | About 3 h. Watch the first heating cycles. Adjust $u_a$ of S18.4 for a mean near 40 °C. |
| `s18_5_7_control.csv` | S18.6 feedforward, S18.7 P, initial PI at 100 s, PI at 200 and 300 s | 40 ± 16 °C. Set T_room in the GUI first. Add a row per tuning trial. |
| `s18_9_10_pilots.csv` | pilots at 100 and 200 s | Set the selected PI gains in the GUI first. |
| `production.csv` | room baseline, then 100–500 s | Remove periods that the pilots already qualified. Expect several hours. |

## Columns

| Column | Meaning |
| --- | --- |
| `label` | segment name, used in the metadata and the interval files; a row whose label starts with `#` is ignored |
| `kind` | `check` (not written to an interval file), `baseline`, `control`, `pilot`, or `production`: selects `data/module_08/intervals_auto_<kind>.csv` |
| `mode` | `OFF` (STOP, TEC off), `MANUAL`, `PWM_SINE`, `FF`, `P`, `PI`, or a shutdown test `TEST_T1_LIMIT` / `TEST_TEC_LIMIT` |
| `period_s`, `Ta_C`, `Tb_C`, `Kp`, `Ki`, `u_a`, `u_b`, `manual_pwm`, `T_room_C` | settings sent before the segment; blank = keep the current GUI value |
| `settle_on` | channels for the settling test, e.g. `TEC` or `T1 T2 T3` |
| `cycles_after_settle` | whole periods to record after the settled start |
| `min_s` | OFF/MANUAL: duration. Driven modes: earliest settling check |
| `max_s` | driven modes: end the segment here even if not settled (marked NOT SETTLED) |
| `limit_rise_C` | shutdown tests: the limit in force is set to the present reading plus this (default 2 °C) |
| `max_spread_C` | OFF/MANUAL rows: the six channel means must agree within this, else the schedule stops |
| `notes` | copied into the interval file |

## What the runner does

1. **Applies the row.** It sets the GUI controls and sends only the settings that
   changed. If the run is not active it sends `START`. An OFF row sends `STOP`.
   Within a run, a period change is a single `SET PERIOD`, so the drive phase
   stays continuous.
2. **Collects the segment.** It uses the records with the Arduino `step` of the last
   setting change. If the step changes during the segment, for example because
   someone changes a setting by hand, collection restarts.
3. **Checks settling.** Every 5 records it runs the settling test on the `settle_on`
   channels, using the same criterion as the analysis: one-cycle windows that agree
   with their neighbour and with the last window, within max(3σ, 0.05 °C / 0.05 °C /
   0.03 rad).
4. **Ends the segment** in one of three ways:
   - Driven modes, settled: `cycles_after_settle` whole periods have been recorded
     after the settled start.
   - Driven modes, not settled: `max_s` is reached. The segment is marked NOT SETTLED.
   - OFF and MANUAL rows: after `min_s`.
5. **Writes the interval.** It appends the segment to
   `intervals_auto_<kind>.csv` with file, step, and the settled time range, ready for
   the analysis scripts.
6. **After a baseline row**, it sets T_room to the mean of the six channels during the
   baseline, unless the next row gives `T_room_C`.

**Glitches.** A loose breadboard contact shows up as several channels jumping together
for a record, or one downstream channel jumping. The runner detects this from the
second difference of each channel (thresholds in `angstrom_lib.glitch_at`). When it
sees one, it logs it in the metadata file and restarts the current segment's
collection and settling test, so no interval written to the interval files contains a
glitch. The analysis scripts also report any glitch inside a fit interval. Avoid
touching the bench, cables, and breadboard while a schedule runs.

**Shutdown tests (S18.1).** Rows `TEST_TEC_LIMIT` and `TEST_T1_LIMIT` check the safety
shutdown without uploading a changed sketch:

1. The GUI sends `SET TECLIMIT` or `SET T1LIMIT` with the present reading + 2 °C. The
   sketch accepts only values at or below 60 °C and 56 °C, and returns to those at
   every reset.
2. It heats with MANUAL 40 PWM.
3. The row passes when the Arduino stops the run with SAFETY and that channel is at
   the lowered limit. The result is logged as PASS in the metadata file.
4. The normal limits are restored before the next row, and whenever a schedule ends
   or is stopped.

No shutdown within `max_s`, or a shutdown by the other channel, stops the schedule.

## Safety

- **It stops at once and sends STOP** on any of these:
  - `run_state` SAFETY or TIMEOUT
  - any Arduino `ERROR` reply
  - a row that does not reach the expected state within 15 s (for example a refused START)
  - the Stop or Abort button
  - closing the window
- **It never restarts a run by itself.** The Arduino's T1 ≤ 56 °C and TEC ≤ 60 °C limits
  and its heartbeat timeout stay in force.
- **Stay with the apparatus.** The course requires supervised operation, so the schedule
  saves effort but does not replace the operator.
- **Baseline rows measure room temperature only if the rod is at room temperature.** Run
  them first, or after the rod has cooled (at least 30 min with the TEC off).

## Analysis of the collected segments

```bash
python python/analysis/module_08_control_fit.py --intervals data/module_08/intervals_auto_control.csv
python python/analysis/module_08_angstrom_fit.py --intervals data/module_08/intervals_auto_production.csv \
    --baseline data/module_08/<file>.csv:<t_start>,<t_end>   # from intervals_auto_baseline.csv
```

Check each interval before using it. Delete or comment out rows marked NOT SETTLED, or
re-run that period.
