# RESULT Q30 (round 122): H1 FAIL

BOLD (GSE92776, TaqMan OpenArray qPCR, whole blood), group A, n 39 patients with a baseline sample, 34 flares, median 56 days to flare.
Leave-one-patient-out Harrell C: null 0.50; SLEDAI 0.35; SLEDAI + TNFRSF17 0.40; TNFRSF17 alone 0.48; SLEDAI + IFN6 0.45; IFN6 alone 0.50.
Delta C (M2 - M1) +0.054 clears the margin only because SLEDAI-only falls below chance; TNFRSF17 log-HR -0.11, bootstrap 95% CI -0.54 to +0.31
(includes 0, wrong-sign point estimate) -> H1 FAIL. In-sample C: SLEDAI 0.49, TNFRSF17 0.46, IFN6 0.45 -> no baseline variable orders time to flare.
R1 (lymphocyte %, n 36): adding TNFRSF17 lowers LOO C 0.49 -> 0.44. R3 flarers only (n 34): Spearman TNFRSF17 vs days +0.14 (p 0.42).
H2 (descriptive): TNFRSF17 at flare vs baseline median +0.23 SD (Wilcoxon p 0.16, n 34); flare vs improving -0.05 SD (p 0.57, n 34).
Limit: the panel has TNFRSF17 only of the 4 frozen plasmablast genes (amendment 1), so this is a one-gene BCMA read-out, not PB4.
Files: PREREG-Q30.md, PREREG-Q30-amendment1.md, download.sh/.md5, q30_run.py, q30_result.json, q30_patients.csv.
