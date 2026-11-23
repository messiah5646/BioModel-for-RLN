# Audit round 120 (rule 10): rounds 116-119

Scheduled run 2026-11-23. Clean worktree at 985cad4; every dataset re-downloaded from its public source; checksums checked.

**Decision this can change:** which r116-r119 labels go into the preprint (P4), and whether rounds-70-72 is ready to merge.
If the numbers reproduce and the labels hold, the audit passes and the labels go in unchanged. If not, downgrade them and mark the branch not ready.

## Reproduction (committed code, fresh data)
| Round | Data re-downloaded | Checksums | Output vs committed |
|---|---|---|---|
| 116 (P3, SLE anchor) | CT.gov JSON in repo | n/a | sri4/bicla tables byte-identical; p3_result.txt identical |
| 117 (P5, ILLUMINATE) | GSE88885/6 matrices and preprocessed expression, GPL17586 | 4/4 OK | result.json 67/67 values identical; per_patient.csv identical |
| 118 (P6, reservoir) | Mackensen 2022 source data MOESM2-8 | 7/7 OK | ratios 0.013/0.028/0.083/0.66/0.114 reproduce (Fig 3F) |
| 119 (Q22, BCMA/CD19) | GSE304878, 30 files | 30/30 OK | result.json byte-identical |

## Findings
1. **r118 factual error. The decision still holds.** "Mackensen: 1 post-treatment point" is wrong. MOESM5 (Fig 4B) has clinical
   anti-dsDNA from months 0 to 15. Two patients have >= 3 post-treatment titres with the last one >= 90 d, so 2 curves qualify, not 0.
   The gate needs 8, so the verdict is still UNDERPOWERED and the reservoir is still not added. r118 also took its ratios from
   the multiplex autoantigen panel (Fig 3F). On the clinical assay the month-3/baseline ratios are 0.0009, 0.0015 and 0.0021
   (all at the assay floor) and 0.089. The <= 0.11 bound is still a valid upper bound, just a loose one. Corrected in RESULT-P6.md.
2. **r117 reproducibility gap (fixed).** download.sh did not fetch gpl_data.txt, which analyse_p5.py needs. Added.
3. **r117 wording.** "70-80% better" is measured against the arm-mean population model (M0). Plain carry-forward of the
   patient's own baseline (M1) is nearly as good: MSE 0.43 vs 0.39 (ILL-1 -> ILL-2) and 0.58 vs 0.55 (reverse). The
   arm shift and the 0.86 shrink add only 6-11%. The set-point claim is about persistence, which still stands. The n should be
   precise: persistence n 135 placebo patients (94 + 41), out-of-trial n 434. **Supported (n 135 placebo, 2 trials)** stays.
4. **r116 n wording.** The primary LOTO test is 23 trials / 45 arms. The 40-trial figure is the all-arm robustness variant
   (same direction, D3 18.8% worse than D1). The label should read "Supported, n 23 trials primary; 40 in the all-arm variant".
5. **r119 robust to metric (exploratory, not pre-registered; r119_sens.py).** The population mean beats both the multiplier
   and carry-forward under absolute error (0.12 vs 0.40 vs 0.58), under |log ratio| (0.33 vs 0.95 vs 1.01), and with P1 or with
   P1 and P2 dropped. M2 vs B0 clears the 10% margin under absolute error but not under the pre-registered metric, so INCONCLUSIVE
   stands. Caveat for Q60: the two patients whose IFN rose had the fewest monocytes (P2 pre: 43 cells), so regression to the
   mean from noisy baselines predicts the same "convergence". Q60 must beat a noise-aware baseline before shrink-to-floor
   replaces the multiplier.

## Tests counted (rounds 116-119)
r116: 1 primary + 6 robustness + 1 Q58 split = 8. r117: 2 primary (A, B) + 3 secondary. r118: 0 (underpowered).
r119: 1 primary + 2 secondary. Total 17 comparisons, of which 4 primary pre-registered tests. Audit r120 adds 5 exploratory
sensitivity metrics (r119) and 1 descriptive source-data re-read (r118).

## Claims on n < 10
r118 bound (n 3-4) and r119 verdict (n 5) are already Emerging. Q58 "early read-out fast-forwards" (n 8 trials) is already
Emerging (r115). Nothing new to downgrade.

## Verdict
**PASS on code** (4/4 rounds reproduce exactly). **One factual correction** (r118, decision unchanged), one reproducibility fix,
two wording fixes. Branch rounds-70-72 is ready to merge (Thomas reviews).
