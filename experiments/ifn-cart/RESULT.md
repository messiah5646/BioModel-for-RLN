# Result: interferon signature after CD19 CAR-T (GSE263931), 2026-11-21
Primary (pre-registered): post < pre in 6 of 7 patients -> PASS. Monocyte-only: 7 of 7.
Exception: pat4 (all-cell +0.06; monocytes still fell 0.07).
Pre-score vs size of drop: Spearman 0.32 (n=7, uninformative).
Caveats: all-cell score is confounded by composition (B cells removed); the monocyte score is the cleaner one and agrees.
pat5 pre has only 428 cells passing QC. Monocyte gate is a crude LYZ threshold, not clustering.
Meaning for LUMEN Layer 1: the cell layer reproduces the published direction on real patient data; it does not yet
show the score predicts anything clinical (no renal outcomes in this dataset).
