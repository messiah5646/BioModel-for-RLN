# PREREG Q22 (round 119): out-of-sample test of the frozen IFN set-point multiplier on GSE304878

Committed BEFORE any expression file of GSE304878 is downloaded (rule 12a, push 1 of 3).

**Decision this can change:** whether LUMEN's blood-IFN multiplier after deep B-lineage depletion stays flat at
F_LATE = 0.73 indefinitely (current `setpoint.py`: no rebound) or needs a rebound term once B cells reconstitute.
If M2 (baseline x 0.73) beats carry-forward and the population mean -> keep the flat plateau and allow it for
BCMA/CD19 dual CAR-T arms. If carry-forward wins (median post/pre ratio >= 0.90) -> add a reconstitution rebound to
`setpoint.py` (f returns toward 1 after B-cell return), labelled Emerging. Otherwise -> no model change.

## What was seen before this PREREG (disclosure)
- GEO series summary and sample titles only: 5 SLE patients, PBMC scRNA (+BCR/TCR), PRE vs POST, where POST is
  "after peripheral B-cell reconstitution following dual-target (BCMA/CD19) CAR-T therapy". One library per sample.
- The paper has NOT been opened. Post-sample day is unknown. `setpoint.py` gives f(t) = 0.73 (+/-0.01) for every
  t >= ~100 d, so the frozen prediction is time-agnostic provided POST is >= 3 months after infusion (assumption A1,
  to be checked against the paper after scoring and reported; it does not change the rule).
- No expression value, figure or IFN statement from this study has been seen.

## Score (frozen, identical to round 96 `scores.py` CAR-T branch)
Cells with total counts > 500; log1p(CP10K); IFN6 = mean over IFI27, IFI44L, IFIT1, ISG15, RSAD2, SIGLEC1;
monocytes = cells with log1p(CP10K) LYZ > 2. Per-sample score = mean IFN6 over its monocytes.
If a gene symbol is missing, stop and report (no substitution).

## Predictors for each patient's POST score (relative error = |pred - obs| / base, as round 96)
- B0 carry-forward: pred = base.
- B1 population mean: mean POST score of the other 4 patients (leave-one-out).
- M2 frozen set point: pred = base x 0.73 (no refit on these data).

## Rules (mean relative error over the 5 patients)
- **PASS (keep flat plateau):** M2 <= 0.9 x B0 AND M2 <= 0.9 x B1.
- **FAIL (rebound):** median POST/PRE ratio >= 0.90 AND B0 <= M2.
- Otherwise **inconclusive** (no model change).
Labels: n = 5 patients < 10, so any outcome is at most **Emerging** (rule 12).

## Reported, never decisive
Per-patient ratios; Spearman(base, post) (set-point persistence, round 117 rho ~0.8); median ratio with a
patient-bootstrap 90% interval; the same three predictors with T-cell IFN6 (TCR-contig-free: cells with CD3E log1p
> 1 and LYZ <= 2); a QC sensitivity with IFN6 residualised on log counts, log genes and mito fraction (per-cell
linear fit pooled over all monocytes) - if the verdict flips, report both and call it inconclusive.
