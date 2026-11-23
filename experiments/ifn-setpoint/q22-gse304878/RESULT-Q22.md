# Q22 result (round 119): frozen IFN multiplier 0.73 on GSE304878 (BCMA/CD19 CAR-T, 5 SLE patients)
PREREG-Q22.md committed first (971c199), download second (4d9d18d), analysis third. analyse.py, out/result.json.
A1 checked after scoring: POST = median 90 days post-infusion (He et al., Arthritis Rheumatol 2026, doi 10.1002/art.70117);
setpoint.py f(90) = 0.74, so the frozen 0.73 applies.

Monocyte IFN6, mean relative error (n 5): carry-forward B0 0.948, population mean B1 0.580, frozen M2 (base x 0.73) 0.894.
Median post/pre 0.43 (bootstrap 90% 0.22-3.10); Spearman(base, post) 0.40.
Per patient ratio: P1 3.10, P2 1.64, P3 0.36, P4 0.43, P5 0.22.
VERDICT (pre-registered): INCONCLUSIVE. M2 beats carry-forward by 6% (needed 10%) and loses to the population mean.
FAIL (rebound) condition not met: median ratio 0.43 < 0.90.
T-cell IFN6 (secondary): B0 0.583, B1 0.204, M2 0.529, median 0.42 -> INCONCLUSIVE. QC-residualised: INCONCLUSIVE (no flip).
Not pre-registered, descriptive: the two ratio>1 patients had the lowest baselines (P1 0.075, P2 0.30) and the fewest
monocytes (P1 192/92, P2 43/244 cells pre/post); the three high-baseline patients fell to 0.22-0.43, deeper than 0.73.
Pattern = convergence to a common post-treatment level (population mean wins), not a proportional multiplier.
Decision per PREREG: no change to setpoint.py. Label: Emerging (n 5).
