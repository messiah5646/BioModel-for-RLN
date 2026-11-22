# Q15 result (round 105)
Part A (patient-level): BLOCKED — no open dataset with per-patient month-1 exposure (MPA AUC/trough or CD19 count),
month-3 UPCR and 12-month CR. Searched: GEO sample index (62,446 samples; only GSE116006 has a 'drug exposure' field =
binary on/off in an anti-IL-6 SLE trial), PubMed/PMC MPA-exposure and B-cell-depletion LN cohorts (open full texts
report only group medians/ROC; per-patient tables only in closed-access n=10 study, below the n>=20 rule).
Part B (binormal bound, primary month-3 UPCR AUC 0.844 from Ren Fail 2026 n=317):
- MPA-AUC (ROC 0.70): B-CONTINGENT — dAUC +0.029 at rho 0, >=0.05 only at rho <= -0.2.
- MPA-AUC (meta-analysis mean/SD, 0.77): B-CONTINGENT — +0.052 at rho 0, +0.032 at rho 0.2.
- Complete B-cell depletion, LUNAR (0.62): B-FAIL — max +0.028.
Needed standalone exposure AUC: 0.77 (rho 0), 0.84 (rho 0.3). Script bound.py, numbers result.json.
Exploratory (not pre-registered): over a baseline-only model (AUC 0.746) MPA 0.70 adds +0.054 at rho 0, +0.017 at 0.2.
