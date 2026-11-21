# Result: interferon signature after CD19 CAR-T (GSE263931), 2026-11-21
Primary (pre-registered): post < pre in 6 of 7 patients -> PASS. Monocyte-only: 7 of 7.
Exception: pat4 (all-cell +0.06; monocytes still fell 0.07).
Pre-score vs size of drop: Spearman 0.32 (n=7, uninformative).
Caveats: all-cell score is confounded by composition (B cells removed); the monocyte score is the cleaner one and agrees.
pat5 pre has only 428 cells passing QC. Monocyte gate is a crude LYZ threshold, not clustering.
Meaning for LUMEN Layer 1: the cell layer reproduces the published direction on real patient data; it does not yet
show the score predicts anything clinical (no renal outcomes in this dataset).

## Replication: CD3xCD19 bispecific (11 SLE donors, monocytes), as pre-registered
- Month 3 vs baseline: fell in 8 of 11 -> PASS.
- Day 28 vs baseline: fell in 6 of 11 -> FALSIFIED (secondary).
- Reading: the IFN fall is real but slow; it is not present at 4 weeks in about half the donors and is by 3 months.
  Exploratory, not tested: at 1.2 ug/kg all 3 donors fell by M3; at 0.3 ug/kg 3 of 6.
- For LUMEN Layer 2: the link from B-cell depletion to interferon is a delayed one (weeks to months), which the
  six-state model currently has no term for. This is a concrete, data-backed structural gap.
