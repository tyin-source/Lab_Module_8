# Module 8, Part I: Angstrom Method — Team Notes

Working record for S18–S22: sensor map, settings, the S18 evidence table, predictions,
and run log. A5 links to this file; it is not rewritten for A5.

## 1. Parameters carried forward

| Quantity | Value | Source |
| --- | --- | --- |
| Rod | 6061-T6, $d$ = 12.70 ± 0.02 mm, overall length 28.66 ± 0.05 cm | Module 7, `module_07_rod_instrument.md` §1 |
| Analysis length $L$ (first hole to tip) | 26.90 ± 0.05 cm | Module 7 §1 |
| Hole positions ($x$ from the first hole) | 0, 1.2, 2.4, 3.7, 5.0, 7.6, 10.2, 15.3, 20.4 cm, ±1.3 mm | Module 7 §1 |
| Mean-profile decay, finite insulated-tip fit | $m$ = 4.90 ± 0.04 m⁻¹, $1/m$ = 20.4 cm, $mL$ = 1.32 | Module 7, `module_07_stationary_fin.md` (Run 2) |
| Implied $\nu/\kappa = m^2$ | 24.0 ± 0.4 m⁻² | from $m$ |
| Transverse Biot number | $\mathrm{Bi}_\perp=(md/4)^2$ = (2.42 ± 0.04) × 10⁻⁴ | Module 7 |
| Channel agreement | spread 0.94 °C in room air, no corrections; per-channel offsets cancel in amplitude and phase | Module 7 §5 |
| Heating susceptibility $\chi_H$ (TEC alone) | 0.485 °C/PWM | Module 4 |
| Cooling susceptibility $\chi_C$ (TEC alone) | 0.172 °C/PWM | Module 4 |
| $r=\chi_H/\chi_C$ | 2.82 | Module 4 |
| $\chi_H$ with the rod mounted | 0.44 (Run 1), 0.39 (Run 2) °C/PWM | Module 7 steady PWM and TEC rise at 45 °C |
| One-lump time constant $\tau$ (TEC alone) | 70 ± 3 s | Module 6, Part 3 |
| PI gains | $K_p$ = 6 PWM/°C, $K_i$ = 0.04 PWM/(°C s) (Module 7 with rod); 0.08 was the Module 6 best without the rod | Modules 6, 7 |
| Course planning values | $\kappa$ = 6.90 × 10⁻⁵ m²/s, $\nu$ = 1.30 × 10⁻³ s⁻¹ ($k$ = 167 W/(m K), $\rho$ = 2700 kg/m³, $c$ = 896 J/(kg K), $h$ = 10 W/(m² K)); $m$ = 4.34 m⁻¹ | Module 8 page, popup on dimensionless lengths |

### Course settings (Module 8 page)

| Item | Value | Where |
| --- | --- | --- |
| Requested TEC temperature | $T_\mathrm{set}=T_a+T_b\sin(\omega t)$ with $T_a$ = 40 °C, $T_b$ = 16 °C (24–56 °C); try $T_b$ = 10 °C if 16 °C is poorly tracked | S18 objective, S18.5, Stage 4 |
| Control periods | 100, 200, 300 s | S18.7 |
| Production periods | 100, 200, 300, 400, 500 s, ≥ 5 settled cycles, ≥ 3 rod sensors measurable | Stage 5 |
| S18.3 raw PWM | $u=u_a+u_b\sin(\omega t)$, $u_a$ = 0, start $u_b$ = 50 at 300 s, then 200 and 100 s with the same $u_b$ | S18.3 |
| S18.3 simulation | one-lump model, room and initial 22 °C, $\tau_\mathrm{th}$ = 80 s, $\chi_h$ = 0.46 °C/PWM, two values of $\chi_c$ per period | S18.3 |
| S18.4 positive-only | 100 s, same $u_b$, $u_a$ chosen for a mean near 40 °C with $0<u_a-u_b$ and $u_a+u_b\le255$ | S18.4 |
| Feedforward | $u_a=(T_a-T_\mathrm{room})/\chi_h$, $u_b=T_b/\chi_h$; stays on the heating branch for 40 ± 16 °C | S18.5 |
| Common time fit | $T_i(t)=a_i+B_i\cos(\omega t+\phi_i)$; unwrapped $\phi$ vs $x$ has slope $-q'$; lag = $\phi_\mathrm{ref}-\phi_\mathrm{resp}$ | Common time fit |
| Pilot targets | farthest retained sensor $B_\mathrm{far}\approx B_1/e$; unwrapped lag T1 → farthest 0.5–2 rad | Stage 4 |

Note: on the page $\tau$ is the forcing period, not the Module 6 relaxation time.

$\chi$ and $\tau$ were measured on the TEC alone. With the rod mounted, the Module 7
steady states give $\chi_H$ 10–20 % lower. The feedforward uses the Module 4 values
(S18.5), and feedback corrects the remainder. The S18.3 zero-centered runs give the
effective values with the rod: $\tau=\tan\varphi/\omega$ and
$\chi=\Delta T\sqrt{1+\omega^2\tau^2}/A$ at each period. The Angstrom results
($q$, $q'$, $\kappa$, $\nu$) use T1–T5 only and do not depend on $\chi$ or $\tau$.

## 2. Oscillatory sensor map (S18.1)

Starting layout: T1 in the first hole and T2–T5 about 1 inch apart, using existing
holes. The positions are the hole positions measured with a ruler in Module 7
(`module_07_rod_instrument.md`, section 1; ±1.3 mm, half the hole diameter); they match
`SENSOR_POSITIONS_CM` in `python/module_08_common.py`.

Holes are numbered 1–9 from the TEC end. Module 7 used holes 1, 5, 7, 8, 9; Module 8
uses holes 1, 3, 5, 6, 7 (T1 stays; T2–T5 move toward the base with their leads,
dividers, and Arduino inputs).

| Lead label | Arduino input | CSV column | Module 7 hole | Module 8 hole | $x$ (cm), measured in Module 7 |
| :---: | --- | --- | :---: | :---: | ---: |
| 1 | A0 | `TEC_C` | TEC base | TEC base | — |
| 2 | A1 | `T1_C` | 1 | 1 | 0.0 ± 0.13 |
| 6 | A2 | `T2_C` | 5 | 3 | 2.4 ± 0.13 |
| 3 | A3 | `T3_C` | 7 | 5 | 5.0 ± 0.13 |
| 4 | A4 | `T4_C` | 8 | 6 | 7.6 ± 0.13 |
| 5 | A5 | `T5_C` | 9 | 7 | 10.2 ± 0.13 |

Channel identity after the move: checked by unplugging each thermistor in turn; only its
own column read `nan`.

Planning check of this layout with the course planning values (semi-infinite model):

| Period (s) | $1/q$ (cm) | $1/q'$ (cm) | $B/B_1$ at 7.6 cm | $B/B_1$ at 10.2 cm | Lag at 10.2 cm (rad) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 4.6 | 4.7 | 0.19 | 0.11 | 2.15 |
| 200 | 6.5 | 6.8 | 0.31 | 0.21 | 1.51 |
| 300 | 7.9 | 8.4 | 0.38 | 0.27 | 1.22 |
| 400 | 9.0 | 9.8 | 0.43 | 0.32 | 1.04 |
| 500 | 10.0 | 11.0 | 0.47 | 0.36 | 0.92 |

The layout meets the Stage 4 targets ($B_\mathrm{far}\approx B_1/e=0.37\,B_1$, lag
0.5–2 rad) at 200–500 s. At 100 s T5 is somewhat beyond them (0.11 $B_1$, 2.15 rad),
with T4 close to the targets. The pilots decide; these are predictions, not results.

The Module 7 full-length map (0, 5.0, 10.2, 15.3, 20.4 cm) applies only to the Module 7
data and is kept there. If the pilots move a sensor, add a new table with the date and
the first file it applies to.

## 3. Safety limits and settings

| Item | Value | Where it is set |
| --- | --- | --- |
| TEC thermistor software limit | 60 °C | `TEC_LIMIT_C` in the sketch |
| Rod-base T1 limit | 56 °C (S18.1: at or below 56 °C, with margin for overshoot) | `ROD_BASE_LIMIT_C` in the sketch |
| Setpoint range $T_a\pm T_b$ | 10–56 °C (allows the S18.5 40 ± 16 °C) | `SETPOINT_MIN_C`, `SETPOINT_MAX_C` |
| Shortest forcing period | 100 s (S18.1) | `MIN_PERIOD_S` |
| Heartbeat timeout | 3 s | `COMMAND_TIMEOUT_MS` |
| Hardware thermal switch | about 70 °C | in series with the TEC |
| Record interval | about 1.0 s (500 × 2 ms) | `ADC_SAMPLES`, `SAMPLE_INTERVAL_US` |

With $T_a+T_b$ = 56 °C the requested peak equals the T1 limit. T1 sits below the TEC
(Module 7: 4.5 K lower at 45 °C), so the 60 °C TEC limit is likely to act first; overshoot
may still stop a run. If it does, reduce $T_b$ (Stage 4 suggests 10 °C) and record the change.

Verification of the limits: the first rows of `schedules/s18_3_4_raw_pwm.csv`
(`test_tec_limit`, `test_t1_limit`) lower the limit in force to the present reading + 2 °C
with `SET TECLIMIT` / `SET T1LIMIT`, heat with MANUAL 40, and pass when the Arduino stops
with SAFETY at that limit; 56 and 60 °C are restored afterwards (and at every reset). The
metadata file records each result as PASS. No upload of a changed sketch is needed. Record the result
in the S18.1 row.

Before long runs: with the TEC off and the GUI recording, gently move the USB cable,
the Arduino 5 V/GND jumpers, and the TEC supply leads, and check that no channel jumps.
This checks the loose contact behind the Module 7 power dropouts.

## 4. S18 evidence table

One row per action. Link the saved file, figure, code checkpoint, or note, state the
principal result, and record any departure from the procedure.

| Step | Action | Evidence (file, figure, commit) | Principal result | Departures |
| --- | --- | --- | --- | --- |
| S18.1 | Channels, units, positions, safety limits | Sensor map (section 2); `data/module_08/angstrom_20261116-211726_meta.txt` (rows `check_channels`, `test_tec_limit`, `test_t1_limit`, `room_before` of `schedules/s18_3_4_raw_pwm.csv`) | Channels identified by unplugging each thermistor (its column reads `nan`); the first heating shows T1–T5 responding in hole order (amplitude falling, lag growing). Channel spread with the TEC off: 0.94 °C at the start (rod still slightly warmer than the TEC), 0.23 °C in the room baseline; below 1 °C, no corrections. Positions: Module 7 hole positions, ±1.3 mm. Shutdown tests PASS: TEC limit lowered to 23.67 °C, stop at TEC 23.89 °C after 12 s; T1 limit lowered to 24.37 °C, stop at T1 24.42 °C after 47 s; limits restored to 60 and 56 °C. Room temperature 21.73 °C. Period ≥ 100 s enforced by the sketch | Positions from the Module 7 hole measurement; limits tested by lowering them by command instead of warming T1 by hand |
| S18.2 | Live display and saved record | `data/module_08/angstrom_20261116-204536.csv` (1226 records) and `angstrom_20261116-211726.csv` (9237 records), each with `_meta.txt`; sketch and GUI at the commits named in the metadata | The GUI showed all six channels, setpoint, mode, and PWM live; every record saved with the 18 documented fields, 0 rejected; placeholders `nan` documented in `data/module_08/README.md`; every start, setting and period change, and safety event logged with its step | — |
| S18.3 | Zero-centered sinusoidal PWM, 300/200/100 s | `data/module_08/intervals_control.csv`; `data/module_08/analysis/control_s18_3_7.csv`; figures `docs/figures/module_08/control_s18_3_pwm{300,200,100}.png`, `control_s18_3_7_comparison.png`; derivation check, course simulation table, and 100 s model comparison: section 8.1 (`sim_s18_3_*`) | $u_b$ = 50 at all periods. Mean − room 3.66, 3.91, 3.81 °C; fundamental amplitude 7.24, 5.59, 3.28 °C; lag behind PWM 51.3, 61.0, 75.2°; residual amplitude 1.17, 0.62, 0.56 °C (period 300, 200, 100 s). Same $u_b$ gives a smaller amplitude and larger lag at shorter period; the mean is raised because heating is stronger than cooling. One-lump fit: $\tau_\mathrm{eff}$ = 59.5, 57.5, 60.3 s; $\chi_\mathrm{eff}$ (average of heating and cooling) 0.231, 0.231, 0.257 °C/PWM | The schedule's first settling test was too lenient (fixed in commit d78e827); re-checked: at 300 s and 100 s amplitude and phase were steady after the first cycle (excluded) but the mean still rose 0.2 and 0.07 °C per cycle as the rod warmed. Not repeated; the mean values carry this drift |
| S18.4 | Positive-only sinusoidal PWM, 100 s | same files; figure `control_s18_4_pos100.png`; simulation table and aligned plots: section 8.2 (`sim_s18_4.csv`, `sim_s18_4_aligned.png`) | $u_a$ = 55, $u_b$ = 50: mean 44.36 °C, amplitude 4.77 °C, lag 75.3°, residual amplitude 0.06 °C against 0.56 °C for zero-centered PWM at 100 s: the positive-only response is nearly sinusoidal. Static heating susceptibility with the rod $(44.36-21.73)/55$ = 0.41 °C/PWM; with the zero-centered mean shift $(\chi_h-\chi_c)u_b/\pi$ = 3.81 °C this gives $\chi_c$ ≈ 0.17 °C/PWM, r ≈ 2.4 (Module 4: 0.485, 0.172, 2.82) | A mean near 40 °C is not reachable with $u_b$ = 50 and positive PWM ($u_a > u_b$ gives at least about 42 °C); $u_a$ = 55 kept, mean 44.4 °C. Settled after about 16 cycles as the rod warmed |
| S18.5 | Setpoint and feedforward construction | `schedules/s18_5_7_control.csv`; `data/module_08/angstrom_20261117-083725_meta.txt` | $T_\mathrm{set}=40+16\sin(2\pi t/100)$ °C. Room baseline 23.10 °C (spread 0.13 °C). $u_a=(40-23.10)/0.485$ = 34.8, $u_b=16/0.485$ = 33.0 PWM: 1.8–67.8, all on the heating branch, within 0–255 | Module 4 $\chi_h$ used as the course asks; with the rod the static $\chi_h$ is about 0.41 °C/PWM (S18.4) |
| S18.6 | Feedforward applied | `intervals_control.csv` row `s18_6_ff100`; `analysis/control_s18_3_7.csv`; `docs/figures/module_08/control_s18_6_ff100.png` | Prediction from the Module 4 calibration: mean 40.0 °C, amplitude 16 °C, lag 0. Measured: mean 37.12 °C (error −2.88), amplitude 3.19 °C (ratio 0.20), lag 76.9°, RMS tracking error 11.5 °C, residual 0.07 °C. The one-lump model with the rod ($\chi_h$ 0.41 °C/PWM, $\tau$ 60 s from S18.3–S18.4) predicts mean 37.4 °C, amplitude 3.5 °C, lag 75°: open loop cannot correct the lower susceptibility with the rod, nor the $1/\sqrt{1+\omega^2\tau^2}$ attenuation and $\tan^{-1}\omega\tau$ lag | First attempt (`angstrom_20261117-000509`) stopped by the heartbeat timeout; repeated |
| S18.7 | P and PI, tuning at 100/200/300 s | rows `s18_7_p100`, `s18_7_pi100`, `s18_7_pi200`, `s18_7_pi300` of `intervals_control.csv`; `control_s18_3_7.csv`; figures `control_s18_7_*.png`; focused 100 s P/PI simulation: section 8.3 (`sim_control.csv`, `sim_control_p_pi100.png`) | Gains in PWM/°C and PWM/(°C s); sample interval 1.0 s; no saturation. P ($K_p$ 6): mean 35.05 °C (error −4.95), ratio 0.48, lag 47.0°, RMS 10.1 °C. PI ($K_p$ 6, $K_i$ 0.04): mean error −0.01 °C at every period (integral action removes the droop); ratio 0.50, 0.63, 0.70 and lag 52.4, 35.9, 28.8° at 100, 200, 300 s. One-lump model predicts P: 0.48, 47.5°, mean 35.1 °C; PI: 0.50/51.5°, 0.68/33.5°, 0.74/25.5°. Model-based tuning (`schedules/s18_7_tuning.csv`, run `angstrom_20261117-114802`): at 100 s $K_p$ 15, $K_i$ 0.15: ratio 0.825, lag 29.5° (model 0.86, 32°); $K_p$ 30, $K_i$ 0.3: ratio 0.937, lag 15.4°, mean error +0.08 °C, RMS 3.4 °C (model 0.96, 18°); with $T_b$ = 10 °C: 0.937, 15.9°. The same $K_p$ 30, $K_i$ 0.3 at 200 s: 0.953, 9.2°, RMS 2.0 °C; at 300 s: 0.964, 6.9°, RMS 1.5 °C (model 0.97, 9°; 0.97, 6.5°). No sustained saturation: at 100 s with $K_p$ 30 the cooling PWM touched −255 in 0.4 % of the fitted records (two brief touches); at 200 and 300 s it stayed above −121 and −75; no retuning needed between periods. **Selected: PI, $K_p$ = 30 PWM/°C, $K_i$ = 0.3 PWM/(°C s) at 100, 200, 300 s.** Measured rod boundary T1 amplitude 4.6, 6.2, 7.2 °C | Initial gains from Module 7; gains chosen with the one-lump model fitted to S18.3–S18.4; one glitch (T3, −0.53 °C) restarted the first tuning segment; the 200 s fit starts one cycle after the period change |
| S18.8 | Boundary-response comparison | `data/module_08/analysis/control_s18_3_7.csv`; `docs/figures/module_08/control_s18_3_7_comparison.png` and per-run figures | Acceptance criteria for the requested 40 ± 16 °C: \|mean error\| ≤ 0.2 °C, amplitude ratio ≥ 0.90, lag ≤ 20°, no sustained PWM saturation. Raw PWM (S18.3–S18.4) has no temperature target: zero-centered PWM is distorted and raises the mean by 3.7–3.9 °C because heating is stronger than cooling; positive-only PWM is nearly sinusoidal (residual 0.06 °C) but its mean is set by the calibration. Feedforward at 100 s: error −2.9 °C, ratio 0.20, lag 77° (fails: open loop, rod lowers $\chi_h$, one-lump attenuation and lag). P ($K_p$ 6): error −5.0 °C (proportional droop), ratio 0.48, lag 47° (fails). PI ($K_p$ 6, $K_i$ 0.04): error −0.01 °C but ratio 0.50–0.70, lag 29–52° (fails: loop gain too low). PI ($K_p$ 30, $K_i$ 0.3): ratio 0.937/0.953/0.964, lag 15.4/9.2/6.9°, error ≤ 0.08 °C at 100/200/300 s (passes all criteria at every period). The one-lump model ($\chi$ 0.41 °C/PWM, $\tau$ 60 s) predicted every case within a few per cent, so the improvement comes from loop gain: integral action removes the mean error, higher $K_p$ and $K_i$ raise the amplitude ratio and cut the lag, and the cooling PWM limit bounds how far this can go at 100 s | Within the tested settings the full target was reached; amplitude ratio stays below 1 and lag above 0, as the model predicts for finite gain |
| S18.9 | 100 s pilot: settling, fits, layout | `schedules/s18_9_10_pilots.csv`, run `angstrom_20261117-132324` (step 22); `data/module_08/intervals_production.csv` row `p100`; `analysis/angstrom_production.csv`, `angstrom_production_sensors.csv`; figures `angstrom_p100_fits.png`, `angstrom_p100_settling.png`, `angstrom_p100_profile_theory.png` | PI $K_p$ 30, $K_i$ 0.3, 40 ± 16 °C. Settling criterion: successive one-cycle windows of T1–T5 agree with each other and with the last window in mean, amplitude, and phase within max(3σ, 0.05 °C / 0.05 °C / 0.03 rad), σ = larger of the fit uncertainty and the measured cycle-to-cycle scatter (amplitudes scatter 1–1.5 % per cycle at 100 s, random, from room air). Settled 3500 s after the start (rod warming); 13 cycles fitted. All five sensors measurable: B = 4.58, 2.74, 1.61, 0.94, 0.54 °C; $B_5/B_1$ = 0.12, T1→T5 lag 2.02 rad; residual amplitude ≤ 0.18 °C. q = 20.93 ± 0.38 m⁻¹, q′ = 19.88 ± 0.36 m⁻¹ (χ²/dof 0.07, 0.04); finite-tip correction ≤ 0.1 %. κ = (7.55 ± 0.19) × 10⁻⁵ m²/s, ν = (3.2 ± 1.6) × 10⁻³ s⁻¹ | The schedule's first settling test used only the fit uncertainty, did not pass within 4800 s, and ended the segment NOT SETTLED; the criterion now includes the measured scatter (commit after b41eb33) and the step was re-analysed from the raw file. Layout kept: $B_5/B_1$ below the 1/e target but 0.54 °C is far above the noise |
| S18.10 | 200 s pilot | same run (step 23); row `p200`; figures `angstrom_p200_*.png` | Settled 2800 s after the period change; 5 cycles. B = 6.22, 4.27, 2.94, 2.01, 1.35 °C; $B_5/B_1$ = 0.22, T1→T5 lag 1.41 rad. q = 14.86 ± 0.25, q′ = 13.87 ± 0.23 m⁻¹; T5 finite-tip phase correction 0.007 rad exceeds its uncertainty, so the finite-rod fit is used: κ = (7.58 ± 0.18) × 10⁻⁵ m²/s, ν = (2.0 ± 0.7) × 10⁻³ s⁻¹. κ agrees with 100 s (χ² 0.02, p = 0.89). Production plan: same layout and settings; 100 and 200 s qualified by the pilots; `schedules/production.csv` runs 300, 400, 500 s after a 10 min room baseline, ≥ 6 settled cycles each, 1 record per s, TEC and T1 limits as in section 3 | — |

## 5. Predictions before S18 (P3)

The numbers come from the planning table in section 2 (course values $\kappa$ = 6.90 × 10⁻⁵
m²/s, $\nu$ = 1.30 × 10⁻³ s⁻¹, semi-infinite rod), committed on 13 Nov before any
Module 8 run (commit 2897083). The wording was added on 17 Nov; the measured values in
brackets came later and were not used for the predictions.

- **$T(t)$ at the base and two downstream positions.** All three follow sinusoids at the
  forcing period. Downstream the oscillation is smaller and reaches its peak later, and
  it oscillates about a lower mean (the mean decays as $\cosh[m(L-x)]/\cosh(mL)$ with
  $1/m\approx$ 23 cm). At 200 s: $B(5.0\text{ cm})/B_1\approx e^{-5.0/6.5}$ = 0.46,
  peak later by $5.0/6.8$ = 0.74 rad (23 s); at 10.2 cm 0.21 and 1.51 rad (48 s).
  [Measured, 200 s pilot: 0.47 and 0.70 rad; 0.22 and 1.41 rad.]
- **Amplitude and phase with distance.** $\ln B$ falls linearly and the unwrapped phase
  lag grows linearly, with slopes $q\gtrsim q'$ (nearly equal, because
  $\Omega=\omega/\nu\gg1$ at these periods). Small departures are expected near the tip
  only; at most 2 % at 10.2 cm for 500 s.
- **Effect of a shorter period.** $q$ and $q'$ grow roughly as $\sqrt{\omega}$, so
  amplitude drops faster and lag accumulates faster. At 100 s, T5 at 10.2 cm keeps
  only about 11 % of $B_1$, about 0.5 °C for $B_1\approx$ 5 °C, with 2.2 rad of lag. The
  TEC itself follows less well at 100 s (one-lump low-pass, $\omega\tau\approx$ 4.4 with
  $\tau$ 70 s), so the controller must work harder there. [Measured: 0.12 and 2.02 rad;
  PI amplitude ratio 0.94 at 100 s, 0.96 at 300 s.]
- **Pilot settings.** 100 s and 200 s; PI control of the TEC thermistor to
  $T_a$ = 40 °C, $T_b$ = 16 °C, starting gains $K_p$ 6, $K_i$ 0.04 from Module 7 (final:
  $K_p$ 30, $K_i$ 0.3 after S18.7); TEC limit 60 °C, T1 limit 56 °C; one record per
  second; run until the settling test passes, then at least 5 cycles (about 25 min at
  100 s, 50 min at 200 s, plus 30–60 min of rod warming).

## 6. Run log

| File | Steps | Mode and settings | Interval used | Notes |
| --- | --- | --- | --- | --- |
| `angstrom_20261116-204536` | — | MANUAL, TEC off | — | S18.2 GUI check |
| `angstrom_20261116-211726` | 17–35 | `s18_3_4_raw_pwm.csv` | `intervals_control.csv` | S18.1 checks, S18.3, S18.4 |
| `angstrom_20261117-000509` | 17–22 | `s18_5_7_control.csv` | baseline only | Stopped at 10 min by TIMEOUT: the GUI's event loop was slowed (macOS App Nap, display asleep) and the heartbeat paused for 3 s. Fixed: the serial thread now sends the heartbeat (watchdog 30 s) and the GUI keeps macOS awake. Repeat the schedule |
| `angstrom_20261117-083725` | 17–26 | `s18_5_7_control.csv` (repeat) | `intervals_control.csv` | Complete: room 23.10 °C, S18.6 feedforward, S18.7 P and initial PI at 100/200/300 s; no timeout, no glitch |
| `angstrom_20261117-114802` | 20–26 | `s18_7_tuning.csv` | `intervals_control.csv` | PI tuning; selected $K_p$ 30, $K_i$ 0.3; one glitch (T3) restarted the first segment |
| `angstrom_20261117-132324` | 17–23 | `s18_9_10_pilots.csv` | `intervals_production.csv` (p100, p200) | Pilots at 100 and 200 s, both qualified for production. The "room" baseline (step 17, 23.90 °C, spread 0.82 °C) was taken while the rod was still cooling from the tuning run (T5 falling 0.9 °C per 5 min); the room temperature for the mean profile is the extrapolated asymptote, 23.34–23.50 °C (section 10) |
| `angstrom_20261117-161210` | 17–26 | `production.csv` | `intervals_production.csv` (p300, p400, p500) | Room baseline 22.81–22.96 °C (spread 0.14 °C, rod at room); 300, 400, 500 s settled after 1500, 1200, 2000 s, 6 cycles each, no glitch; 25 min cool-down recorded after the schedule (room drift, section 10) |

## 7. Theory: finite insulated-tip rod and the four checks

Rod equation for $\theta=T-T_\infty$: $\partial_t\theta=\kappa\,\partial_x^2\theta-\nu\theta$,
with $x=0$ at T1 (the measured boundary) and the insulated tip at $x=L$ = 26.90 cm.
Boundary conditions $T(0,t)=T_a+T_b\cos\omega t$ and $\partial_xT(L,t)=0$.

**Mean.** $\bar\theta''-m^2\bar\theta=0$, $m^2=\nu/\kappa$,
$\bar\theta(0)=T_a-T_\infty$, $\bar\theta'(L)=0$. The general solution
$A\cosh[m(L-x)]+B\sinh[m(L-x)]$ has slope $-m\{A\sinh[m(L-x)]+B\cosh[m(L-x)]\}$, which is
$-mB$ at $x=L$, so $B=0$; then $A\cosh(mL)=T_a-T_\infty$.

**Oscillation.** $\tilde\theta=F(x)e^{i\omega t}$ gives $i\omega F=\kappa F''-\nu F$, so
$F''=s^2F$ with $s^2=(\nu+i\omega)/\kappa$. Writing $s=q+iq'$ gives
$q^2-q'^2=\nu/\kappa$ and $2qq'=\omega/\kappa$. The same steps with $F(0)=T_b$ and
$F'(L)=0$ give $F=T_b\cosh[s(L-x)]/\cosh(sL)$. Hence

$$T(x,t)=T_\infty+(T_a-T_\infty)\frac{\cosh[m(L-x)]}{\cosh(mL)}
+\mathrm{Re}\left\{T_b\frac{\cosh[s(L-x)]}{\cosh(sL)}e^{i\omega t}\right\}.$$

**Check 1 ($x=0$).** Both ratios equal $\cosh(mL)/\cosh(mL)=1$ and
$\cosh(sL)/\cosh(sL)=1$. So $T(0,t)=T_\infty+(T_a-T_\infty)+\mathrm{Re}\{T_be^{i\omega t}\}
=T_a+T_b\cos\omega t$. ✓

**Check 2 (slope at $x=L$).** $\partial_x\cosh[m(L-x)]=-m\sinh[m(L-x)]$ and
$\partial_x\cosh[s(L-x)]=-s\sinh[s(L-x)]$. Both vanish at $x=L$ because $\sinh 0=0$.
The real part of a zero complex amplitude is zero, so the oscillation has zero slope at
the tip at every $t$. ✓

**Check 3 (long rod).** Write
$\dfrac{\cosh[m(L-x)]}{\cosh(mL)}=e^{-mx}\,\dfrac{1+e^{-2m(L-x)}}{1+e^{-2mL}}$.
When $m(L-x)\gg1$ (hence also $mL\gg1$), the fraction tends to 1 and the profile tends
to $e^{-mx}$. The same factorization with $s$ gives
$\dfrac{\cosh[s(L-x)]}{\cosh(sL)}=e^{-sx}\,\dfrac{1+e^{-2s(L-x)}}{1+e^{-2sL}}$ with
$|e^{-2s(L-x)}|=e^{-2q(L-x)}$. When $q(L-x)\gg1$ this is $e^{-sx}=e^{-qx}e^{-iq'x}$:
amplitude $T_be^{-qx}$, lag $q'x$. That is the semi-infinite solution. ✓

**Check 4 (why the tip region is most sensitive).** The finite and semi-infinite
solutions differ only by the factor
$(1+e^{-2s(L-x)})/(1+e^{-2sL})$. Its numerator is the wave reflected from the insulated
tip, which has travelled $2(L-x)$ and has decayed by $e^{-2q(L-x)}$. That term is
largest where $L-x$ is smallest, and at $x=L$ it equals 1: the tip temperature is the
sum of the incident and reflected waves, which doubles the amplitude relative to
$e^{-qL}$. So sensors near the tip carry the largest model difference, and only they
can tell the two models apart. For the mean, $m(L-x)$ at T5 is 5.5 × 0.167 = 0.92 and
$e^{-1.84}$ = 0.16. The mean profile therefore needs the finite form even at T5, which
is why the Stage 6 mean fit is finite and the oscillation fit can be semi-infinite.
With our $\kappa$, $\nu$ (7.57 × 10⁻⁵ m²/s, 2.24 × 10⁻³ s⁻¹) the oscillation's tip
factor at T5 (10.2 cm) grows with period: bound
$(e^{-2q(L-x)}+e^{-2qL})/(1-e^{-2qL})$ = 0.10 %, 0.72 %, 1.7 %, 2.9 %, 4.1 % at 100–500 s.
The predicted T5 corrections are amplitude +0.09 %, −0.06 %, −1.4 %, −2.7 %, −3.3 %
and phase +0.000, −0.007, −0.009, +0.002, +0.018 rad. Measured uncertainties at T5
are about 1–2 % in amplitude and 0.007–0.02 rad in phase, so the correction matters
from 200 s on. Because $\nu=\kappa(q^2-q'^2)$ is a small difference, `--model auto`
then uses the finite fit (section 9).

## 8. Simulations (S18.3, S18.4, focused S19 control, integrated comparison)

All simulations are in `python/analysis/module_08_simulation.py`; one command
reproduces every table and figure below (`python python/analysis/module_08_simulation.py`).
The script replaces the course's one-lump app and six-sensor simulator with the same
equations. It uses only our measured parameters, and its controller is a line-by-line
port of the sketch's `updateControl()`.

### 8.0 Model and calibration

- **TEC, one lump.** $\tau\dot\theta_T=\chi(u)u-\theta_T-g(\theta_T-\theta_\mathrm{rod}(0))$,
  with $\chi=\chi_h$ for heating and $\chi_c$ for cooling. $g=G/H$ is the TEC-to-rod
  conductance over the TEC-to-room conductance. With $g=0$ and no rod this is the
  Module 6 one-lump model.
- **Rod.** The rod equation above on the full 28.66 cm, from the TEC face to the tip.
  - Base: a contact boundary $-\partial_X\theta=(\theta_T-\theta)/\ell_c$. Tip: insulated.
  - T1 sits at $X_1$ and the other sensors at $X_1+x_i$.
  - Finite differences with 2 mm cells, integrated exactly over each 1.008 s record
    with the matrix exponential.
  - The PWM is held between records, and each record is the interval average (as the
    Arduino does).
- **Controller.** The sketch's PWM_SINE, FF, P, and PI (conditional integration,
  cooling × r = 2.82, clamp ±255, integer PWM), fed back from the TEC node like the
  apparatus. The course simulator feeds back from the rod base T0. Ours uses the same
  sensor as the experiment, so the gains need no conversion.

| Parameter | Value | How it was measured (`data/module_08/analysis/sim_parameters.csv`) |
| --- | --- | --- |
| $\kappa$, $\nu$ | (7.566 ± 0.085) × 10⁻⁵ m²/s, (1.91 ± 0.22) × 10⁻³ s⁻¹ | weighted five-period Angstrom results (`angstrom_production.csv`, section 10) |
| $\ell_c$, $X_1$ | 4.57 ± 0.14 cm, 2.05 ± 0.07 cm (geometric TEC face to hole 1: 1.76 cm) | fit of the TEC → T1 transfer $\cosh[s(L_r-X_1)]/(\cosh sL_r+\ell_cs\sinh sL_r)$ to the PI runs `s18_7_pi100_t2`, `pi200_sel`, `pi300_sel`: T1/TEC 0.308, 0.409, 0.468 (model 0.305, 0.411, 0.473); lag 0.862, 0.640, 0.539 rad (model 0.851, 0.650, 0.547) |
| $\chi_h$, $\chi_c$, $\tau$, $g$ (coupled) | 0.495 ± 0.017, 0.198 ± 0.013 °C/PWM, 77.7 ± 2.9 s, 1.2 ± 0.2 | the recorded PWM of run `angstrom_20261116-211726` (S18.3 at 300/200/100 s, S18.4, 40 min cool-down; 7560 records) replayed into the model and fitted to the TEC channel; rms residual 0.23 °C |
| one-lump $\chi_h$, $\chi_c$, $\tau$ | 0.412, 0.173 °C/PWM, 59.1 s | measured with the rod: S18.4 static rise, S18.3 mean shift $(\chi_h-\chi_c)u_b/\pi$, mean S18.3 $\tau_\mathrm{eff}$ |

The coupled lump's own values are close to the Module 4 and 6 values for the TEC
alone ($\chi_H$ 0.485, $\chi_C$ 0.172 °C/PWM, $\tau$ 70 ± 3 s). The rod is a second load
that lowers the effective $\chi$ and $\tau$ to the one-lump values measured with the
rod mounted. The model predicts the static mean ratio $(T_1-T_\infty)/(T_\mathrm{TEC}-T_\infty)$
as 0.758; the 300 s production interval measured (35.81 − 22.90)/(40.00 − 22.90) = 0.755.
This is an independent check, since the contact was fitted to oscillation data only.

Validation over the whole raw-PWM run (`docs/figures/module_08/sim_raw_pwm_replay.png`):
- **Fitted.** Only the TEC channel was fitted.
- **Predicted T1–T5.** The oscillations of T1–T5 follow the measurement. RMS residuals
  over all 7560 records are 0.44, 0.46, 0.64, 0.76, 0.86 °C for T1–T5.
- **Discrepancy.** Most of that residual is in the slow mean: during S18.4 and the
  cool-down the measured rod is up to about 1.5 °C warmer at T5 than the model, and
  cools more slowly. The model has no thermal mass for the TEC heat sink or the mount,
  so its slow heat exchange is too fast. With the pilot-only
  $\nu$ = 2.24 × 10⁻³ s⁻¹ the T5 residual was 1.26 °C. With the five-period
  $\nu$ = 1.91 × 10⁻³ s⁻¹ it dropped to 0.86 °C, which supports the smaller $\nu$. This
  discrepancy concerns the mean and the slow transients, not the oscillation used for
  $\kappa$ and $\nu$.

### 8.1 S18.3: check the derivation and the simulation table

**Substitution.** Put $\theta=a\sin\omega t+b\cos\omega t+Ce^{-t/\tau}$ into
$\tau\dot\theta+\theta=\chi u_b\sin\omega t$. The exponential gives
$\tau(-C/\tau)e^{-t/\tau}+Ce^{-t/\tau}=0$, so it solves the homogeneous equation for
any $C$. The sinusoidal part gives $a-\omega\tau b=\chi u_b$ and $b+\omega\tau a=0$,
the coefficients on the page. At $t=0$, $\theta(0)=b+C$. With
$C=\chi u_b\omega\tau/(1+\omega^2\tau^2)=-b$, this is 0, so $T(0)=T_\mathrm{room}$. ✓

**Startup.** The particular solution alone would start at $b<0$, below room. The
exponential cancels that offset at $t=0$ and decays with $\tau$. At 300 s, $\tau$ 80 s,
$\chi$ 0.46, $u_b$ 50: $C$ = 10.12 °C. The first half-cycle rides above the settled
sinusoid, and after $3\tau$ = 240 s (less than one period) under 5 % of $C$ remains.

**Analytical vs numerical.** The one-lump integration at 0.1 s and the boxed solution
agree to $3.8\times10^{-6}$ °C over five cycles, startup included
(`sim_s18_3_analytic300.png`). Settled amplitude 11.79 °C; lag 59.2° (49 s).

**Course simulator check.** The course page's own one-lump mode
(`angstrom-six-sensor/?s18=3`, run headlessly by `tools/course_simulator_lumped.js`;
`course_sim_lumped.csv`) gives the same six rows to every digit shown. Its numerical and
analytical curves agree to 1 × 10⁻⁶ °C.

**Course table.** One-lump, room and initial 22 °C, $\tau_\mathrm{th}$ 80 s, $\chi_h$ 0.46 °C/PWM,
$u_a$ 0, $u_b$ 50 (`data/module_08/analysis/sim_s18_3_course.csv`):

| Period (s) | $\chi_c$ (°C/PWM) | (1) Mean − room (°C) | (2) Amplitude (°C) | (3) Phase lag (°) | (4) Residual amplitude (°C) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 100 | 0.23 | 3.66 | 3.37 | 78.7 | 0.240 |
| 100 | 0.46 | 0.00 | 4.49 | 78.7 | 0.000 |
| 200 | 0.23 | 3.66 | 6.38 | 68.3 | 0.473 |
| 200 | 0.46 | 0.00 | 8.50 | 68.3 | 0.000 |
| 300 | 0.23 | 3.66 | 8.84 | 59.2 | 0.695 |
| 300 | 0.46 | 0.00 | 11.79 | 59.2 | 0.000 |

- **Equal susceptibilities.** These rows are the boxed solution:
  - zero mean shift;
  - amplitude $\chi u_b/\sqrt{1+\omega^2\tau^2}$ = 23/5.13, 23/2.71, 23/1.95;
  - lag $\tan^{-1}\omega\tau$ = 78.7°, 68.3°, 59.2°;
  - residual zero up to numerical error.
- **Changing $\chi_c$ at fixed period.** Lowering $\chi_c$ from 0.46 to 0.23 °C/PWM:
  - adds the $|\sin|$ term, which raises the mean by $(\chi_h-\chi_c)u_b/\pi$ = 3.66 °C at
    every period (it is a DC term, so the low-pass does not filter it);
  - scales the fundamental by $(\chi_h+\chi_c)/2\chi_h$ = 0.75;
  - leaves the lag unchanged;
  - adds a residual at $2\omega$.
- **Changing period at fixed $\chi_c$.** At longer period the amplitude grows and the
  lag shrinks, as the low-pass predicts. The residual grows from 0.24 to 0.70 °C,
  because the filter attenuates the $2\omega$ term less at low frequency.
- **Residual shape.** Its period is half the PWM period (`sim_s18_3_course_100s.png`).

**Model–experiment comparison at 100 s.** Both models are driven by the recorded
PWM, starting from the measured room and initial temperatures
(`sim_s18_3_pwm100_compare.csv`, `.png`):

| Source | Mean − room (°C) | Amplitude (°C) | Lag behind PWM (°) | Residual amplitude (°C) |
| --- | ---: | ---: | ---: | ---: |
| measured (`s18_3_pwm100`) | 3.81 | 3.28 | 75.2 | 0.56 |
| one-lump, measured $\chi_h$ 0.41, $\chi_c$ 0.17, $\tau$ 59 s | 3.64 | 3.80 | 75.0 | 0.39 |
| coupled TEC–rod model | 3.70 | 3.21 | 72.2 | 0.50 |

**What the unequal-susceptibility one-lump model reproduces:**
- the mean shift (from $\chi_h-\chi_c$);
- the lag;
- the twice-per-cycle residual.

**Where it falls short:**
- **Amplitude.** It overpredicts the amplitude by 16 %: with the rod attached, the
  susceptibility at 100 s is smaller than the static one, because the rod is a
  frequency-dependent heat load.
- **Slow warming.** It has no slow rod warming. The measured mean was still rising
  0.07 °C per cycle, which also adds to the measured residual.

The coupled model fixes the amplitude but lags 3° less than measured.

### 8.2 S18.4: simulation table

100 s, room 21.73 °C, the same measured parameters for both drives
(`sim_s18_4.csv`, aligned PWM/temperature/residual plots in `sim_s18_4_aligned.png`):

| 100 s drive | Source | $u_a$ | $u_b$ | Mean (°C) | Fundamental amplitude (°C) | Residual amplitude (°C) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| zero-centered | measured | 0 | 50 | 25.54 | 3.28 | 0.561 |
| zero-centered | one-lump | 0 | 50 | 25.50 | 3.79 | 0.339 |
| zero-centered | coupled TEC–rod | 0 | 50 | 25.66 | 3.21 | 0.314 |
| positive-only | measured | 55 | 50 | 44.36 | 4.77 | 0.062 |
| positive-only | one-lump | 55 | 50 | 44.36 | 5.34 | 0.011 |
| positive-only | coupled TEC–rod | 55 | 50 | 44.43 | 4.59 | 0.010 |

The course simulator's one-lump mode, with the same measured parameters, gives
3.80 °C / 0.335 °C (zero-centered) and 5.35 °C / 0.000 °C (positive-only)
(`course_sim_lumped.csv`). This matches our one-lump row. The 0.011 °C residual in ours
comes only from integer PWM; the course model uses continuous PWM.

**Positive-only PWM gives a more sinusoidal temperature**, in experiment and in both
models:
- the residual falls from 17 % of the fundamental to 1.3 % (measured) or 0.2 %
  (models);
- the whole cycle stays on the heating branch, so the $|\sin|$ term and its $2\omega$
  harmonic vanish.

**Remaining distortion.** The measured 0.06 °C is mostly not waveform distortion:
- it is slow drift and noise (the rod was still warming);
- the 0.01 °C simulated comes from integer PWM steps.

**The models do not capture all of the zero-centered residual** (0.31–0.34 against
0.56 °C):
- part of the measured value is the rising mean;
- the 2nd harmonic is not fully reproduced by a constant $\chi_c$.

### 8.3 Focused control simulation at 100 s (S19)

Coupled model, $T_a$ 40 °C, $T_b$ 16 °C, room as in each measured run, P and PI with
the same $K_p$, then one targeted adjustment, then the selected gains at 200 and 300 s
(`sim_control.csv`, `sim_control_p_pi100.png`). Gains in PWM/°C and PWM/(°C s).

| Case | Period (s) | $K_p$ | $K_i$ | Source | Mean error (°C) | Amplitude ratio | Lag (°) | RMS error (°C) | Saturated |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| P only | 100 | 6 | 0 | simulation | −5.03 | 0.448 | 47.8 | 10.23 | 0 % |
| P only | 100 | 6 | 0 | measured `s18_7_p100` | −4.95 | 0.482 | 47.0 | 10.06 | 0 % |
| PI, same $K_p$ | 100 | 6 | 0.04 | simulation | +0.00 | 0.466 | 52.2 | 9.26 | 0 % |
| PI, same $K_p$ | 100 | 6 | 0.04 | measured `s18_7_pi100` | −0.01 | 0.501 | 52.4 | 9.22 | 0 % |
| PI, adjusted | 100 | 30 | 0.3 | simulation | +0.84 | 0.903 | 16.4 | 3.67 | 13 % |
| PI, adjusted | 100 | 30 | 0.3 | measured `s18_7_pi100_t2` | +0.08 | 0.937 | 15.4 | 3.37 | 0.4 % |
| PI selected | 200 | 30 | 0.3 | simulation | −0.00 | 0.949 | 9.5 | 2.07 | 0 % |
| PI selected | 200 | 30 | 0.3 | measured `s18_7_pi200_sel` | −0.00 | 0.953 | 9.2 | 2.01 | 0 % |
| PI selected | 300 | 30 | 0.3 | simulation | +0.00 | 0.962 | 7.0 | 1.54 | 0 % |
| PI selected | 300 | 30 | 0.3 | measured `s18_7_pi300_sel` | −0.00 | 0.964 | 6.9 | 1.51 | 0 % |

**What integral action improves.** Adding $K_i$ = 0.04 at the same $K_p$ removes the
mean error, −5.0 °C → 0.0 °C, in simulation and in experiment. It barely changes the
oscillation: ratio 0.45 → 0.47 and lag 48° → 52°. The reason is that at 100 s,
$K_i/\omega$ = 0.64 PWM/°C, small next to $K_p$ = 6. The integral acts at low
frequency (the mean), not at the forcing frequency, and the low-frequency phase lag it
adds makes the lag slightly worse.

**What is still short, and the adjustment.** Amplitude and phase fall short because
the loop gain $K_p\chi$ ≈ 6 × 0.4 ≈ 2.4 is too low against the TEC's low-pass
($\omega\tau$ ≈ 4–5). The one adjustment, both gains × 5 ($K_p$ 30, $K_i$ 0.3), raises
the ratio to 0.90 and cuts the lag to 16° (measured 0.937 and 15.4°).

**Discrepancy at 100 s with the adjusted gains.** The model needs more cooling than
the apparatus:
- the model is at the −255 limit 13 % of the time, which gives it a +0.84 °C mean error;
- the apparatus touched −255 in only 0.4 % of records.

The cooling susceptibility was fitted with |u| ≤ 50 near room temperature, and a
constant $\chi_c$ does not describe the TEC at large cooling currents around a 40 °C
mean. At 200 and 300 s there is no saturation, and model and measurement agree within
0.6 % in ratio and 0.3° in lag. That supports the selected gains without retuning at
each period.

The simulation feeds back from the TEC node, exactly like the apparatus, so the gains
are used unchanged. (The course simulator's T0 feedback would need gains converted by
the TEC → T1 ratio of 0.3–0.5.)

### 8.4 Integrated comparison: theory, experiment, simulation (Stage 6 item 5)

`python python/analysis/module_08_simulation.py rod <label> --independent` takes a
production interval and simulates the same PI settings to a periodic steady state, from
room temperature.
- **Inputs.**
  - Room: the interval's own baseline in `intervals_production.csv`.
  - $\kappa$ and $\nu$: from the other four periods only.
  - Contact and TEC lump: as in 8.0.
- **Figure** (`sim_rod_<label>_independent.png`):
  - simulation and measurement aligned on the phase of the requested setpoint;
  - the simulated rod base against the measured T1 first, then T3 and T5;
  - the analytical rod prediction driven by the measured T1, with
    $m=\sqrt{\nu/\kappa}$;
  - a residual panel (measured − model).

**Matched run: 300 s** (`sim_rod_p300_independent.png`, `.csv`). Room 22.90 °C;
$\kappa$ = 7.578 × 10⁻⁵ m²/s and $\nu$ = 1.954 × 10⁻³ s⁻¹ from p100, p200, p400,
p500. Lags are relative to the requested setpoint:

| Sensor | Mean meas / sim (°C) | Amplitude meas / sim (°C) | Lag meas / sim (°) |
| --- | ---: | ---: | ---: |
| TEC | 40.00 / 40.00 | 15.43 / 15.39 | 6.8 / 7.0 |
| T1 | 35.81 / 35.81 | 7.26 / 7.28 | 37.7 / 38.3 |
| T2 | 34.30 / 34.53 | 5.29 / 5.42 | 53.2 / 53.6 |
| T3 | 33.29 / 33.33 | 3.88 / 3.93 | 70.0 / 70.2 |
| T4 | 32.36 / 32.32 | 2.82 / 2.84 | 86.6 / 86.7 |
| T5 | 31.50 / 31.47 | 2.03 / 2.04 | 103.3 / 103.4 |

**Inputs: measured, fitted, assumed.**
- *Measured:* positions, room, and the controller settings.
- *Fitted to other data:*
  - the contact (PI control runs);
  - the TEC lump (raw-PWM run);
  - $\kappa$ and $\nu$ (the other four periods).
- *Assumed:* the course $\rho c$ does not enter, and the rod is one-dimensional.

So this is an independent test of the 300 s record. With $\kappa$, $\nu$ fitted to
p300 itself (`sim_rod_p300.csv`) the amplitudes and lags are the same to within 1 %,
but the means come out 0.3–0.6 °C high. That is because p300's own $\nu$
(1.73 × 10⁻³ s⁻¹) is at the low end of the five values.

**Agreement.** The model gets the boundary right and carries it down the rod:
- **Boundary.** The simulated TEC matches the measured TEC within 0.3 % in amplitude
  and 0.2° in lag. T1 matches within 0.3 % and 0.6°, so the PI loop and the contact
  are reproduced.
- **Downstream.** T3–T5 then follow within 1.3 % in amplitude, 0.2° in lag, and
  0.04 °C in mean.
- **Residuals.** They stay within ±0.15 °C over three cycles and are mostly common to
  all channels (slow room-air wander), not shaped like a missing harmonic.

**Discrepancy at T2.** T2's simulated amplitude is 2.5 % too large (5.42 against
5.29 °C) and its mean 0.23 °C too high. The same pattern appears at every period:

| Period (s) | 100 | 200 | 300 | 400 | 500 |
| --- | ---: | ---: | ---: | ---: | ---: |
| T2 amplitude sim/meas − 1 (%) | −2.6 | +2.3 | +2.5 | +2.5 | +2.9 |

At 100 s the simulated boundary itself is 4 % low (see below). Relative to T1, T2 is
1.6–2.2 % high at every period, while its lag is within 0.3–0.6°. A sensor displaced
downstream would show up in both amplitude and phase. A small reading-gain or
thermal-contact error at T2 would show up in amplitude only, which fits better. It
is within what the T2 position uncertainty (±1.3 mm, 1.6 % in amplitude at 300 s)
plus a small calibration error allow. The spatial fits carry it, since the
uncertainty of each fitted point includes the position error.

**Boundary-control vs rod-model errors.** At 100 s the model's TEC saturates in
cooling (8.3), so its TEC runs 0.8 °C warm and 4 % low in amplitude. T1–T5 inherit
the same 4 % deficit and +0.5–0.9 °C means. That is a boundary-control error (the
cooling model), not a rod-model error: at 200–500 s, with no saturation, every
downstream sensor except T2 agrees within 1.5 %.

**Analytical traces.** The traces from the measured T1 agree with the numerical
model and the data in amplitude and phase. Their means follow
$m=\sqrt{\nu/\kappa}$ = 5.03 m⁻¹, which matches the measured mean profile
(section 10).

Other periods, all independent (`sim_rod_p*_independent.csv`):
- **100 s:** amplitudes 2–4 % low and means 0.5–0.9 °C high, from the boundary
  saturation above.
- **200 s:** amplitudes within 1.5 % (T2 2.3 %), lags within 0.6°.
- **400 s:** amplitudes within 1.3 % (T2 2.5 %), lags within 0.4°.
- **500 s:** amplitudes within 1.6 % (T2 2.9 %), lags within 0.9°; means 0.1–0.4 °C
  high.

### 8.5 The course's six-sensor simulator (S19 focused simulation and Stage 6)

The course page's simulator
(<https://sethfraden.github.io/Phys39F26-course/simulators/angstrom-six-sensor/>) was run
headlessly with `node tools/course_simulator_driver.js data/module_08/course_simulator`.
The driver:
- loads the page;
- sets its own controls: coupled TEC + rod, aluminum 6061, $d$ 12.7 mm, $L$ 285 mm
  (slider step 5 mm; our rod is 286.6 mm), Angstrom sinusoid 40 + 16 sin, Kd off;
- presses Apply / Reset;
- advances the page's own `advance()` in 1 s steps and records every second.

Nothing in the model code is changed. Page SHA-256 and settings are in
`course_simulator/course_sim_meta.json`. The S18.3 and S18.4 tables (8.1, 8.2) were
also run in its one-lump mode (`tools/course_simulator_lumped.js`). So every simulation
the course asks for has been done in the course's own simulator, and our Python model is
the additional, calibrated comparison.

**Screenshots in the course simulator's own format** (`docs/figures/module_08/course_ui/`,
made by `node tools/course_simulator_screens.js docs/figures/module_08/course_ui
data/module_08/course_simulator/ui_csv`). Each case is set up through the page's
controls, Apply / Reset is pressed, and the page's model is run to the time shown. Each
case also has the page's own **Save CSV** export in `data/module_08/course_simulator/ui_csv/`.

| File | What it shows | Use |
| --- | --- | --- |
| `s18_3_P{100,200,300}_chic{0.23,0.46}.png` | One-lump mode, settled readout, and "Show settled PWM, temperature, and residual" | S18.3 table; aligned plots for both $\chi_c$ at 100 s |
| `s18_3_P300_equal_startup.png` | 300 s, equal $\chi$, first 600 s: numerical (white) on the dashed analytical curve, startup included | S18.3 "Check the derivation" |
| `s18_4_zero100.png`, `s18_4_pos100.png` | One-lump mode with our measured $\chi_h$ 0.412, $\chi_c$ 0.173, $\tau$ 59.1 s, room 21.73 °C; residual view | S18.4 simulation table |
| `s19_p100_thermal.png`, `s19_pi100_thermal.png`, `s19_pi100_adj_thermal.png` | Coupled model at 100 s after 4 h: P (Kp 2.353 %/°C), PI (Ki 0.0157), adjusted PI (11.76, 0.1176) | Focused S19 simulation |
| `stage6_pi300_visual.png`, `stage6_pi300_thermal.png` | Coupled model at 300 s after 5 h, $L$ 285 mm, final PI gains | Stage 6 matched simulation (A5 Q4) |

The settled readouts printed by the page reproduce the tables above. For example,
100 s, $\chi_c$ 0.23: mean − room 3.661 °C, amplitude 3.366 °C, lag 78.7°, residual
0.240 °C. In the coupled screenshots, the page's own sensor tiles show the simulator's
positions (0, 34, 80, 137, 205, 271 mm on 285 mm); the comparison with our holes uses
the rod temperatures at our positions (`course_sim_rod.csv`).

**Python replica of the course simulator (fully automatic).**
`python/analysis/course_simulator_replica.py` is a line-by-line port of the page's model:
- the same constants;
- the TEC Peltier/Joule/heat-sink equations;
- the explicit rod with $N$ = 100 nodes at 285 mm;
- 0.01 s substeps and the 0.1 s PID with its filter, anti-windup and clamps;
- the 2 s sensor lag and the one-lump RK4;
- the chart buffer, the settled readout and residual view, and the Save CSV format.

It does not port the frost layer, which on the page changes only the drawing.

`python python/analysis/course_simulator_replica.py` runs every case above in about
1.5 min, with no browser. It writes the records and page-format CSVs to
`data/module_08/course_simulator_replica/` and the page-style figures to
`docs/figures/module_08/course_replica/`, then reruns the comparison tables with the tag
`replica`.

**Check against the browser runs** (`course_simulator_replica/replica_check.txt`):
- **Records.** Every record agrees within the browser CSV's rounding: 5 × 10⁻⁵ °C over
  14400–18000 s for the coupled cases, 5 × 10⁻⁷ °C for the one-lump cases.
- **Settled readouts.** Identical to the page's (e.g. 3.661 / 3.366 °C / 78.7° /
  0.240 °C).
- **Comparison tables.** The same to 6 significant figures.

So the replica can stand in for the page wherever a figure or rerun is needed.

**Sensor map.** Simulator T0 (its PID feedback, at the TEC face) = our TEC thermistor.
T1–T5 = holes 1, 3, 5, 6, 7, sampled from the simulated rod at 17.6, 41.6, 67.6, 93.6,
119.6 mm from the TEC face (hole 1 = 28.66 − 26.90 cm; then 2.4, 5.0, 7.6, 10.2 cm).

**Gains.** The simulator's controller outputs % duty (±100 %) and has no cooling scale
$r$. So Kp[%/°C] = Kp[PWM/°C]/2.55: Kp 6 → 2.353, Ki 0.04 → 0.0157; Kp 30 → 11.76,
Ki 0.3 → 0.1176. The feedback point is the same on both (TEC face), so no T0 → T1
conversion is needed. Analysis: `python python/analysis/module_08_course_sim.py`.

**What the course simulator fixes, unlike our model:**
- the rod is the handbook 6061 ($k$ = 167 W/(m K), $h$ = 10 W/(m² K));
- the TEC is a physical TEC1-12706 model (12 V);
- the TEC–rod contact is good (15000 W/(m² K)).

**Focused S19 simulation at 100 s** (`course_sim_s19.csv`, `course_sim_s19_p_pi100.png`):

| Case | Source | Mean error (°C) | Amplitude ratio | Lag (°) | RMS error (°C) | Saturated |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| P, Kp 6 | course simulator | −6.03 | 0.560 | 24.7 | 8.62 | 0 % |
| P, Kp 6 | measured | −4.95 | 0.482 | 47.0 | 10.06 | 0 % |
| PI, Kp 6, Ki 0.04 | course simulator | +0.00 | 0.590 | 26.9 | 6.16 | 0 % |
| PI, Kp 6, Ki 0.04 | measured | −0.01 | 0.501 | 52.4 | 9.22 | 0 % |
| PI, Kp 30, Ki 0.3 | course simulator | −0.00 | 0.904 | 8.1 | 1.88 | 0 % |
| PI, Kp 30, Ki 0.3 | measured | +0.08 | 0.937 | 15.4 | 3.37 | 0.4 % |

**What carries over and what doesn't:**
- **The same physics.** At the same Kp, integral action removes the droop (−6.0 → 0.0 °C)
  but barely changes amplitude or phase. Only higher gains raise the ratio above 0.9.
- **The simulated TEC is faster.** Its lag is about half the measured one, because the
  model TEC is stronger and faster than ours: it used at most 50 % duty (128 PWM)
  where the apparatus reached −255.
- **Our calibrated model matches the apparatus better.** It gives 47.8° and 52.2°.

**Stage 6 comparison at all five periods** (`course_sim_rod.csv`;
`course_sim_p300_comparison.png` for the matched 300 s run):

| 300 s | Amplitude meas / course / ours (°C) | Lag behind setpoint (°) | $B/B_1$ | Lag from T1 (°) |
| --- | ---: | ---: | ---: | ---: |
| TEC/T0 | 15.43 / 14.81 / 15.39 | 6.8 / 5.0 / 7.0 | — | — |
| T1 | 7.26 / 11.89 / 7.28 | 37.7 / 14.5 / 38.3 | 1 | 0 |
| T2 | 5.29 / 8.79 / 5.42 | 53.2 / 30.7 / 53.6 | 0.729 / 0.740 / 0.744 | 15.5 / 16.2 / 15.3 |
| T3 | 3.88 / 6.34 / 3.93 | 70.0 / 48.2 / 70.2 | 0.535 / 0.534 / 0.540 | 32.3 / 33.6 / 31.9 |
| T4 | 2.82 / 4.56 / 2.84 | 86.6 / 65.5 / 86.7 | 0.389 / 0.384 / 0.390 | 48.9 / 51.0 / 48.4 |
| T5 | 2.03 / 3.26 / 2.04 | 103.3 / 83.0 / 103.4 | 0.279 / 0.274 / 0.280 | 65.6 / 68.4 / 65.1 |

**Boundary first.** The feedback temperature tracks equally well in all three. But the
simulator's T1 swings 11.9 °C against the measured 7.3 °C and leads by 23°. Its good
contact puts the rod end almost at the TEC temperature, whereas in the apparatus the
TEC–rod contact ($\ell_c$ ≈ 4.6 cm, section 8.0) loses about half the amplitude at 300 s.
That boundary difference carries straight down the rod, so every absolute amplitude is
about 60 % too large.

**Rod transport.** Normalized to T1, the three agree. At 100–500 s:
- **Amplitude.** The simulator's $B_i/B_1$ matches the measurement within 0.014 at every
  sensor.
- **Lag.** Its lag from T1 to T5 is 2.8–6.4° (4.4–5.7 %) larger than measured.

This is a rod-model difference, and it has the size expected from the handbook $k$.
The phase slope scales as $q'\propto\kappa^{-1/2}$, and the simulator's
$\kappa$ = 6.90 × 10⁻⁵ against our 7.57 × 10⁻⁵ predicts about 4.7 % more lag. The
observed excess is 4.4–5.7 % (300 s: 68.4/65.6). So the course simulator, run with
handbook 6061, independently shows that our rod conducts heat about 9 % better than
handbook 6061-T6. That fits the measured $k$ = 183 W/(m K) (section 10).

**Answers to Stage 6 item 5 (A5 Q4):**
- **Feature reproduced:** the decay and lag of the oscillation along the rod, once the
  boundary is matched.
- **Meaningful discrepancy:** the boundary. The simulator's TEC–rod contact is much
  better than the apparatus's, so its T1 has 1.6× the amplitude and 23° less lag. This
  is a boundary/actuator error, distinct from the 4 % lag difference that comes from
  the material.
- **Inputs:** the simulator's parameters are the course's handbook values, not fitted
  to our data, so this is an independent test.

## 9. Stage 6 checks: mean profile and finite length

**Mean profile.** Means minus the per-channel room temperature of each run (column
`baseline` in `intervals_production.csv`), with $\xi_i=x_i/L$, fitted to
$\theta_0\cosh[\mu(1-\xi)]/\cosh\mu$. Then $m=\mu/L$ with
$\sigma_m^2=(\sigma_\mu/L)^2+(\mu\sigma_L/L^2)^2$. $L$ = 26.90 ± 0.14 cm: the ruler's
0.05 cm combined with the T1 thermistor's ±1.3 mm placement in its hole. The room
temperatures come from different runs:
- **Production run** (p300–p500): the TEC-off baseline at the start (rod at room,
  spread 0.14 °C).
- **Pilot run** (p100, p200): the asymptote extrapolated from that run's baseline,
  because the rod was still cooling (section 10).

| Interval | $\mu=mL$ | $m$ (m⁻¹) | $\chi^2/\nu_\mathrm{dof}$ | $\sqrt{\nu/\kappa}$ this period (m⁻¹) |
| --- | ---: | ---: | ---: | ---: |
| p100 | 1.300 ± 0.031 | 4.83 ± 0.12 | 0.5 | 6.55 |
| p200 | 1.353 ± 0.032 | 5.03 ± 0.12 | 0.3 | 5.17 |
| p300 | 1.331 ± 0.014 | 4.95 ± 0.06 | 1.9 | 4.79 |
| p400 | 1.325 ± 0.014 | 4.93 ± 0.06 | 2.0 | 4.92 |
| p500 | 1.340 ± 0.014 | 4.98 ± 0.06 | 2.1 | 5.08 |

**Agreement.** The five mean-profile values agree with each other, as they should:
the mean does not depend on the forcing period. Weighted mean
$m$ = 4.95 ± 0.03 m⁻¹ ($\mu$ = 1.33). It agrees with:
- the Module 7 full-length value, 4.90 ± 0.04 m⁻¹;
- the five-period Angstrom value $\sqrt{\nu/\kappa}$ = 5.03 ± 0.29 m⁻¹.

So the stationary profile and the oscillation imply the same heat loss.

**Uncertainty budget.** $\sigma_L$ adds 0.03 m⁻¹ to the 0.05 m⁻¹ from $\sigma_\mu$.
The room drift during the production run (section 10) is a larger systematic, of
−0.04 to −0.15 m⁻¹.

**Earlier pilot values.** Before the baseline correction, the pilot intervals gave
5.49 and 5.72 m⁻¹. That baseline was too warm, which shrinks every $\theta_i$ by the
same amount and so steepens the apparent decay.

**Finite-length check** at every retained sensor, with the semi-infinite $q$, $q'$ of
each interval. Printed by `module_08_angstrom_fit.py`; columns `finite_*` and `sigma_*`
in `angstrom_production_sensors.csv`:

| Interval | Sensor | $x$ (cm) | Bound | $\lvert H_\mathrm{fin}/H_\mathrm{semi}-1\rvert$ | Amplitude corr. | $\sigma(\ln B)$ | Phase corr. (rad) | $\sigma(\phi)$ (rad) | Semi-infinite OK? |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| p100 | T1 | 0.0 | 2.6e-05 | 0 | 0 | 2.8e-03 | 0 | 2.8e-03 | yes |
| p100 | T2 | 2.4 | 4.8e-05 | 3.0e-05 | −3.0e-05 | 5.0e-03 | +1.3e-06 | 5.0e-03 | yes |
| p100 | T3 | 5.0 | 1.2e-04 | 1.1e-04 | −7.5e-05 | 6.7e-03 | +8.1e-05 | 6.7e-03 | yes |
| p100 | T4 | 7.6 | 3.2e-04 | 3.2e-04 | +5.9e-05 | 1.4e-02 | +3.2e-04 | 1.4e-02 | yes |
| p100 | T5 | 10.2 | 9.3e-04 | 9.3e-04 | +8.7e-04 | 2.2e-02 | +3.3e-04 | 2.2e-02 | yes |
| p200 | T2 | 2.4 | 1.0e-03 | 4.7e-04 | +4.7e-04 | 2.2e-03 | +2.5e-05 | 2.2e-03 | yes |
| p200 | T3 | 5.0 | 1.8e-03 | 1.5e-03 | +1.3e-03 | 3.0e-03 | −6.2e-04 | 3.0e-03 | yes |
| p200 | T4 | 7.6 | 3.6e-03 | 3.4e-03 | +1.8e-03 | 5.6e-03 | −2.9e-03 | 5.6e-03 | yes |
| p200 | T5 | 10.2 | 7.3e-03 | 7.3e-03 | −6.7e-04 | 7.2e-03 | −7.3e-03 | 7.2e-03 | no: phase correction > uncertainty |
| p300 | T2 | 2.4 | 3.5e-03 | 1.4e-03 | +3.9e-04 | 1.8e-03 | −1.3e-03 | 1.8e-03 | yes |
| p300 | T3 | 5.0 | 5.6e-03 | 4.0e-03 | −3.4e-04 | 1.6e-03 | −3.9e-03 | 1.6e-03 | no |
| p300 | T4 | 7.6 | 9.5e-03 | 8.5e-03 | −4.3e-03 | 3.2e-03 | −7.4e-03 | 3.2e-03 | no |
| p300 | T5 | 10.2 | 1.7e-02 | 1.7e-02 | −1.4e-02 | 2.5e-03 | −8.8e-03 | 2.5e-03 | no |
| p400 | T2 | 2.4 | 7.0e-03 | 2.4e-03 | −1.2e-03 | 1.7e-03 | −2.1e-03 | 1.7e-03 | no |
| p400 | T3 | 5.0 | 1.1e-02 | 6.7e-03 | −4.9e-03 | 1.9e-03 | −4.6e-03 | 1.9e-03 | no |
| p400 | T4 | 7.6 | 1.7e-02 | 1.4e-02 | −1.3e-02 | 2.9e-03 | −5.2e-03 | 2.9e-03 | no |
| p400 | T5 | 10.2 | 2.8e-02 | 2.6e-02 | −2.6e-02 | 3.5e-03 | +6.7e-04 | 3.5e-03 | no |
| p500 | T2 | 2.4 | 1.1e-02 | 3.5e-03 | −3.0e-03 | 1.4e-03 | −1.8e-03 | 1.4e-03 | no |
| p500 | T3 | 5.0 | 1.6e-02 | 9.5e-03 | −9.1e-03 | 1.8e-03 | −2.5e-03 | 1.8e-03 | no |
| p500 | T4 | 7.6 | 2.4e-02 | 1.9e-02 | −1.9e-02 | 2.8e-03 | +9.4e-04 | 2.8e-03 | no |
| p500 | T5 | 10.2 | 3.8e-02 | 3.5e-02 | −3.2e-02 | 3.4e-03 | +1.4e-02 | 3.4e-03 | no |

(T1 at $x$ = 0 has zero correction at every period and is omitted after p100.) The
course bound is tight for the sensors farthest from the base.

**Where the semi-infinite model holds:**
- **100 s:** at all five sensors.
- **200 s:** at T1–T4. At T5 the phase correction equals its uncertainty.
- **300–500 s:** the corrections reach 1–3 % in amplitude at T4–T5 and exceed the
  measurement uncertainty from T3 (300 s) or T2 (400–500 s) on.

So `--model auto` uses the finite-rod fit from 200 s on, and the semi-infinite
interpretation is limited to 100 s for this layout.

**Effect of ignoring the tip.** The semi-infinite fits give $\kappa$ low by 2.4 %
(400 s) and 5.1 % (500 s), and $\nu$ high by 22 % and 10 %. Those numbers are the
size of the limit (section 10).

## 10. Stages 5–6: five-period results

**Production record (Stage 5).** All five intervals are in
`data/module_08/intervals_production.csv`:
- **Setup.** The same layout and sensor map throughout; PI $K_p$ 30, $K_i$ 0.3,
  $T_a$ 40 °C, $T_b$ 16 °C.
- **100 and 200 s.** The pilots (S18.9–S18.10).
- **300, 400, 500 s.** Run `angstrom_20261117-161210` (`schedules/production.csv`),
  settled ranges taken from the live settling test. Re-running the settling test
  offline on the whole steps returns exactly the same intervals: 2108–3907, 5115–7515,
  9525–12525 s.
- **Cycles.** 13, 5, 6, 6, 6.
- **Sensors and glitches.** All five sensors retained, no glitch in any interval.
- **Residual drift.** At most 0.16 °C (500 s).

**Room temperature and drift.**
- **Production run.**
  - Before: TEC-off baseline 22.81–22.96 °C (channel spread 0.14 °C).
  - After: the 25 min recorded after the schedule ended. Its cool-down extrapolates
    to 22.45–22.61 °C (fit $T_\infty+Ae^{-t/\tau}$, $\tau$ = 260–340 s).
  - The room fell 0.39 ± 0.14 °C over 3.8 h, about 0.1 °C/h.
  - It does not affect $q$, $q'$, $\kappa$, $\nu$, which come from the oscillation;
    each fit window is 0.5–0.8 h long and the residual drift stays ≤ 0.16 °C.
  - It does affect the mean profile. Interpolating the room linearly between the two
    baselines changes $m$ by −0.04 (300 s), −0.09 (400 s), and −0.15 m⁻¹ (500 s).
- **Pilot run.** Its baseline was taken right after the tuning run, with the rod still
  cooling (T1 −0.46, T5 −0.93 °C per 5 min).
  - The same exponential fit, per channel, gives asymptotes of 23.34–23.50 °C.
  - These keep the per-channel offsets of the clean production baseline (T1 highest,
    by about 0.1 °C). This supports the extrapolation.
  - They are used as $T_\infty$ for p100 and p200 (`…:extrapolate` in the intervals
    file; uncertainty includes a 0.1 °C floor).

**Results (Stage 6).** The table gives the values used: the semi-infinite fit at 100 s,
the finite-rod fit from 200 s on, chosen by `--model auto`. $k=\kappa\rho c$ and
$h=\nu\rho c\,d/4$ use the course $\rho$ = 2700 kg/m³ and $c$ = 896 J/(kg K)
(their uncertainty is Part II). Sources: `angstrom_production.csv`, figures
`angstrom_p*_fits.png`, `angstrom_production_frequency.png`,
`angstrom_production_lengths.png`.

| Period (s) | $B_1$ (°C) | $q$ (m⁻¹) | $q'$ (m⁻¹) | $\kappa$ (10⁻⁵ m²/s) | $\nu$ (10⁻³ s⁻¹) | $k$ (W/(m K)) | $h$ (W/(m² K)) | $\mathrm{Bi}_\perp$ | $q'R$ | model |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 100 | 4.58 | 20.93 ± 0.38 | 19.88 ± 0.36 | 7.55 ± 0.19 | 3.2 ± 1.6 | 183 ± 5 | 25 ± 12 | 4.3e-4 | 0.126 | semi-infinite |
| 200 | 6.22 | 14.86 ± 0.25 | 13.94 ± 0.23 | 7.58 ± 0.18 | 2.0 ± 0.7 | 183 ± 4 | 16 ± 6 | 2.9e-4 | 0.088 | finite |
| 300 | 7.26 | 12.29 ± 0.21 | 11.32 ± 0.19 | 7.53 ± 0.18 | 1.73 ± 0.51 | 182 ± 4 | 13 ± 4 | 2.9e-4 | 0.071 | finite |
| 400 | 7.95 | 10.80 ± 0.20 | 9.61 ± 0.17 | 7.57 ± 0.20 | 1.83 ± 0.41 | 183 ± 5 | 14 ± 3 | 3.1e-4 | 0.061 | finite |
| 500 | 8.43 | 9.82 ± 0.19 | 8.40 ± 0.16 | 7.61 ± 0.21 | 1.97 ± 0.34 | 184 ± 5 | 15 ± 3 | 3.0e-4 | 0.054 | finite |

**Frequency consistency.**
- **$\kappa$.** Independent of period: weighted mean
  (7.566 ± 0.085) × 10⁻⁵ m²/s, $\chi^2$ = 0.12 for 4 dof.
- **$\nu$.** (1.91 ± 0.22) × 10⁻³ s⁻¹, $\chi^2$ = 0.90 for 4 dof.
- **Largest departure.** 100 s, in $\nu$: 3.2 against 1.9, which is 0.8σ, within
  uncertainty. At 100 s, $\nu$ is least constrained because $q\approx q'$ there.
- **What remains of ν.** $q^2-q'^2$ is only 10–27 % of $q^2$, and $\nu$ carries a
  12 % uncertainty after combining.

**Uncertainty caveat.** Each period's uncertainty is dominated by the ±1.3 mm sensor
positions, and those errors are the same in every period:
- **Position-free fits.** Without the position term, $\sigma_\kappa$ per period is
  0.8–2 %. The five values (7.47–7.54) then still agree: $\chi^2$ 0.69 for 4 dof for
  $\kappa$, and 7.5 for 4 dof ($p$ = 0.11) for $\nu$.
- **Combined uncertainty.** Because the position errors are correlated across
  periods, combining five periods cannot beat them down. A fair combined uncertainty
  is the statistical ±0.04 plus a common position term of about ±0.16 × 10⁻⁵ m²/s:
  $\kappa$ ≈ (7.55 ± 0.17) × 10⁻⁵ m²/s (2.2 %).
- **The small $\chi^2$ = 0.12** reflects this double counting, not luck.

Part II will propagate it properly (Monte Carlo with a common position draw).

**Physical numbers.** With the course $\rho c$:
- $k$ = 183 ± 4 W/(m K), 9 % above the 6061-T6 handbook value of 167 W/(m K). Part II
  will check $\rho$, $c$, and the positions as causes.
- $h$ = 14.7 ± 1.7 W/(m² K), in the range expected for natural convection plus
  radiation from a warm horizontal rod.
- The transverse Biot number is ≤ 4 × 10⁻⁴ and $q'R$ ≤ 0.13. The one-dimensional model
  is justified at all five periods.

**Characteristic lengths** (`angstrom_production_lengths.png`). The curves come from
the weighted $\kappa$, $\nu$, with the fitted points overlaid.
- **The two oscillation lengths.** $1/q$ grows from 4.8 to 10.2 cm and $1/q'$ from
  5.0 to 11.9 cm between 100 and 500 s, roughly as $\sqrt{\text{period}}$. They
  separate at long period as $\Omega=\omega/\nu$ falls from 33 to 6.6.
- **The mean length.** $\ell_\mathrm{mean}=1/m$ = 20 cm is independent of period. The
  mean-profile points (5 periods, 4.95 m⁻¹ → 20.2 cm) and Module 7 (20.4 cm) lie in
  its band.
- **For the layout.** Over the 10.2 cm sensor span, the oscillation falls to 12 %
  (100 s) or 35 % (500 s) of $B_1$ and lags 2.0 or 0.88 rad. One radian of lag is
  reached inside the sensor span at 100–300 s. At 400 s and 500 s the lag at T5 is
  0.98 and 0.88 rad, because $1/q'$ = 10.4 and 11.9 cm are just beyond T5. The
  lag is still resolved to ±0.003 rad.

**T2 diagnosis (from existing data, no new run).** For each sensor, the measured value
is compared with the finite-rod prediction from T1. The amplitude and phase use the
five-period $\kappa$, $\nu$. The mean uses a finite profile fitted to the other four
sensors.

| | 100 s | 200 s | 300 s | 400 s | 500 s | T3–T5, same periods |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| T2 amplitude, meas/pred − 1 | −1.5 % | −1.9 % | −2.0 % | −2.0 % | −2.0 % | −0.2 to −1.0 % (T5 at 100 s: −3.5 %, $B$ = 0.5 °C) |
| T2 phase, meas − pred (mrad) | −7 | −1 | +2 | +2 | +2 | −8 to +7 (T5 at 100 s: −27) |
| T2 mean, meas − pred (°C) | −0.16 | −0.13 | −0.14 | −0.14 | −0.15 | 0.00 to +0.05 |
| T2 mean / excess θ | −1.4 % | −1.2 % | −1.2 % | −1.2 % | −1.3 % | 0.0 to +0.6 % |

**T2 reads low by a fixed fraction of its excess over room**, about 1.2–1.4 % in the
mean and 1.5–2.0 % in the amplitude (about 1 % more than T3–T5). That fraction is the
same at every period, and the phase is unchanged.
- **Position error ruled out.** A position error large enough for −2 % in amplitude
  would be $\Delta x$ = 0.02/$q$ = 1.6 mm at 300 s. That would also lag the phase by
  $q'\Delta x$ ≈ 18 mrad, and it would grow with $q$ toward 100 s (3.4 %). Neither is seen.
- **Slow thermal contact ruled out.** A contact time constant would cut the amplitude
  more at short periods and add a lag. Neither is seen either.

**The error follows the probe.** In Module 7 the same probe and channel (lead 6, A2)
sat in hole 5 (5.0 cm). The stationary fit left it at −0.251 K at θ = 12.1 K, i.e.
−2.1 % (`module_07_stationary_fin.md` §4), in a different hole. So the error goes with
the probe and its channel, not with the hole.

**Candidate causes.** Two remain:
- a $\beta$ (gain) difference of that thermistor, which is not tested by a
  one-temperature check;
- a steady lead-conduction loss that pulls that probe toward room.

A two-temperature water bath with T2 and T3 taped together would tell them apart. It
is not needed for the results.

**Effect on the results.** Excluding T2 changes the weighted $\kappa$ to
7.528 × 10⁻⁵ m²/s (−0.5 %), $\nu$ to 2.01 × 10⁻³ s⁻¹ (+5 %), and $m$ by +0.03–0.05 m⁻¹.
All of these are well inside the uncertainties. T2 is therefore kept, and the T2
error is listed as a systematic for Part II.
