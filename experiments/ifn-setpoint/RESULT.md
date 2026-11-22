# Round 96 result: blood IFN as a therapy-lowered patient set point (queue Q10)
Pre-registered in PREREG.md (commit bdd507b) before any per-donor score was computed. 18 donors: 7 CD19 CAR-T, 11 CD3xCD19.
Mean relative error on the late sample (CAR-T post, CD3xCD19 M3), leave-one-donor-out:

| Predictor | Late error | D28 error |
|---|---|---|
| B0 carry-forward | 0.410 | 0.411 |
| B1 population post level | 0.914 | 1.567 |
| M1 own baseline x trial fall | 0.325 | 0.474 |
| M2 own baseline x one shared fall | 0.305 | - |

- H1 (set point persists) FAIL. Error part passes (M1 < 0.9 x B1) but Spearman(base, late) = 0.49 in CD3xCD19 (rule >= 0.5);
  CAR-T 0.82. B1's mean is inflated by one donor with baseline 0.06 (P2); median errors B1 0.347 vs M1 0.305.
- H2 (therapy lowers it) PASS: M1 21% below carry-forward.
- H3 (one shared fall across CAR-T and bispecific) PASS: M2 0.305 <= 1.1 x M1. Late multiplier 0.69 (CAR-T), 0.75 (CD3xCD19), pooled 0.73.
- H4 descriptive: median multiplier 0.98 at day 28 vs 0.75 at month 3. At day 28 carry-forward beats the set-point model.
Exploratory (not pre-registered, exploratory.py): on an absolute scale (error / trial SD) M2 0.49, M1 0.54, B1 0.73, B0 0.79;
dropping P2, M2 0.29 < M1 0.31 < B0 0.42 < B1 0.53. Within-trial rank Spearman pooled 0.61; CD3xCD19 without P2 0.32.
Model change: setpoint.py, IFN_i(t) = S_i x f(t); f = 1 to day 28, then -> 0.73. Lag and half-life are bounded, not identified.
