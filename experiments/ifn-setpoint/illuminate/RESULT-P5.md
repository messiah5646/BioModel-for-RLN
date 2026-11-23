# P5 RESULT (round 117): blood IFN IS a patient set point — PASS in ~450 patients across two placebo-controlled trials

PREREG efd7a52 -> download 9905d91 -> analysis (this commit). ILLUMINATE-1 (GSE88885) and -2 (GSE88886), whole-blood HTA 2.0,
tabalumab 120 mg Q2W/Q4W vs placebo, baseline / week 16 / week 52. All 6 IFN genes present; 908 + 418 samples matched.

## Primary (pre-registered)
**A. Persistence in placebo: PASS.**
| Trial | n | Spearman baseline vs week 52 | ICC over 3 visits |
|---|---|---|---|
| ILLUMINATE-1 | 94 | 0.81 | 0.87 (n 92) |
| ILLUMINATE-2 | 41 | 0.82 | 0.93 (n 37) |

**B. Out-of-trial prediction of week-52 IFN6: set-point model WINS** (gain vs arm mean >= 10% both directions).
| Train -> test | n test | M0 arm mean | M1 carry-forward | M2 baseline + arm shift | M3 shrunk (slope) |
|---|---|---|---|---|---|
| ILL-1 -> ILL-2 | 137 | 1.98 | 0.43 | 0.42 | **0.39** (slope 0.85) |
| ILL-2 -> ILL-1 | 297 | 1.83 | 0.58 | 0.59 | **0.55** (slope 0.87) |
Set-point models are 70-80% better than the population model.

**DECISION: keep the set-point layer.** Use the shrunk form (week-52 = arm shift + ~0.86 x baseline).

## Secondary (reported)
- C. Tabalumab vs placebo change in IFN6 (log2 units): ILL-1 -0.19 (95% CI -0.38 to -0.00; n 203 vs 94); ILL-2 -0.17
  (-0.41 to +0.07; n 96 vs 41). A small fall (~12%) with BAFF blockade, as predicted for a partial B-lineage hit.
- D. Week 16 predicts week 52 better than baseline (MSE 0.47 vs 0.59; 0.31 vs 0.44): an early read-out is informative for IFN.
- E. Baseline anti-dsDNA and C3 do not predict the IFN change (|rho| <= 0.24), consistent with round 113.

## What this changes
- Round 111 had downgraded "blood IFN is a slow patient set point" to Emerging (n 6 and 18). It is now **Supported**
  (n ~450 patients, two independent trials, pre-registered, out-of-trial prediction).
- Round 89's "54% of IFN variance is patient identity" (n 6) was an underestimate: 87-93% in placebo here (whole blood,
  bulk; round 89 was monocyte scRNA).
- Round 96's persistence "H1 FAIL (Spearman 0.49)" came from a single-cell B-depletion cohort; here persistence is
  strong. The difference is the measurement (bulk whole blood vs per-cell-type scRNA) and the larger n, not a contradiction
  of the structure.
- Not shown: the size of the multiplier on round 96's scale (scales differ, see PREREG). Whole-blood IFN is partly a
  cell-mix readout (neutrophils, monocytes), so "set point" here means the patient's whole-blood IFN level.

Labels: set-point persistence and set-point prediction — **Supported** (n 135 placebo / 434 test patients). Tabalumab IFN
effect — Emerging (one trial significant). Week-16 > baseline for prediction — Supported (n 417, both trials).
