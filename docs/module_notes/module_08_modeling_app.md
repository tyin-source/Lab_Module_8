# Module 8: Formative Modeling-App Check

Team note for the progress check. The authoritative app lives under `python/analysis/`;
its code, plots, and residuals become part of C6 and F1.

## Run instructions

```bash
python python/analysis/module_08_analysis_selftest.py
python python/analysis/module_08_angstrom_fit.py --intervals data/module_08/intervals_production.csv \
    --name production          # room baseline per interval (column baseline)
python python/analysis/module_08_simulation.py            # numerical TEC-rod model (notes section 8)
python python/analysis/module_08_simulation.py rod p300 --independent   # matched comparison
```

The app is `python/analysis/module_08_angstrom_fit.py` with `python/analysis/angstrom_lib.py`.
The numerical model is `python/analysis/module_08_simulation.py`. It reads the Angstrom
results, so run it after the fit.

## Input format

Raw CSV from `python/module_08_angstrom_gui.py` (18 columns; see
`data/module_08/README.md`) and its `_meta.txt` file for sensor positions and events.

## Sequence of fits

1. Settled interval: one-cycle windows of T1 and the rod sensors; settled from the first
   window after which every window agrees with its neighbour and with the last window in
   mean, amplitude, and phase within max(3σ, floor) (floors 0.05 °C, 0.05 °C, 0.03 rad).
   Then the whole cycles after it (at least 5 for production).
2. Common time fit $T_i=a_i+B_i\cos(\omega t+\phi_i)$ for each sensor, with
   uncertainties inflated for correlated residuals; residual amplitude and residual drift.
3. Retained sensors: not excluded and $B_i\ge\max(0.05\ ^\circ\mathrm{C},3\sigma_B)$.
4. Weighted fits of $\ln B_i$ and unwrapped $\phi_i$ against $x_i$ (±1.3 mm included):
   $q$, $q'$.
5. $\kappa=\omega/(2qq')$, $\nu=\kappa(q^2-q'^2)$.
6. Model check: the same data fitted to the finite insulated-tip rod; used instead of
   the straight-line fits where the tip correction exceeds the measurement uncertainty.
7. Means minus the TEC-off baseline, with $\xi=x/L$, fitted to
   $\theta_0\cosh[\mu(1-\xi)]/\cosh\mu$; $m=\mu/L$ with the uncertainties of $\mu$ and
   $L$ propagated; $m$ compared with $\sqrt{\nu/\kappa}$ and Module 7.
8. Weighted means of $\kappa$ and $\nu$ over the periods with a $\chi^2$ consistency test;
   characteristic-length plot.

## Preliminary results

All five production periods (notes section 10):

| Period (s) | $q$ (m⁻¹) | $q'$ (m⁻¹) | $\kappa$ (10⁻⁵ m²/s) | $\nu$ (10⁻³ s⁻¹) | model |
| --- | --- | --- | --- | --- | --- |
| 100 | 20.93 ± 0.38 | 19.88 ± 0.36 | 7.55 ± 0.19 | 3.2 ± 1.6 | semi-infinite |
| 200 | 14.86 ± 0.25 | 13.94 ± 0.23 | 7.58 ± 0.18 | 2.0 ± 0.7 | finite |
| 300 | 12.29 ± 0.21 | 11.32 ± 0.19 | 7.53 ± 0.18 | 1.73 ± 0.51 | finite |
| 400 | 10.80 ± 0.20 | 9.61 ± 0.17 | 7.57 ± 0.20 | 1.83 ± 0.41 | finite |
| 500 | 9.82 ± 0.19 | 8.40 ± 0.16 | 7.61 ± 0.21 | 1.97 ± 0.34 | finite |

Weighted: $\kappa$ = (7.57 ± 0.09) × 10⁻⁵ m²/s, $\chi^2$ 0.12 for 4 dof. The ±1.3 mm
positions are common to all periods; with them counted once, the uncertainty is
±0.17 × 10⁻⁵ m²/s. $\nu$ = (1.91 ± 0.22) × 10⁻³ s⁻¹, $\chi^2$ 0.90 for 4 dof.
Mean profile: $m$ = 4.95 ± 0.03 m⁻¹, against $\sqrt{\nu/\kappa}$ = 5.03 ± 0.29 m⁻¹ and
Module 7's 4.90 ± 0.04 m⁻¹.

## Approximation checks

- Finite length: the bound $(e^{-2q(L-x)}+e^{-2qL})/(1-e^{-2qL})$ and the actual correction
  at every retained sensor (table in `module_08_angstrom.md` section 9). The
  semi-infinite model holds at all sensors at 100 s. From 200 s on, the tip correction
  exceeds the measurement uncertainty: at T5 at 200 s, and from T2–T3 on at
  300–500 s, reaching 3 % at T5. Those intervals use the finite-rod fit.
- Transverse Biot number: $\mathrm{Bi}_\perp=\nu d^2/(16\kappa)$ = (2.9–4.3) × 10⁻⁴
  at the five periods, and $q'R$ ≤ 0.13. The cross-section is uniform to better than
  0.1 % for the mean and 1 % for the oscillation.

## Exclusions

None. All five sensors are retained at all five periods (amplitudes ≥ 0.54 °C,
SNR > 40), and no glitch falls in any fit interval. The room temperature of the pilot run
is extrapolated, because its baseline was taken while the rod was still cooling
(notes section 10).

## AI assistance

Claude Code wrote:
- the analysis library, both scripts, and the synthetic-data self-test;
- the numerical TEC–rod model and its calibration (`module_08_simulation.py`).

See the README's AI use note. The team checked the outputs against the raw data and
the course equations.
