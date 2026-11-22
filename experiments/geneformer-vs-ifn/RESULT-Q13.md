# Round 98 result (Q13): QC-adjusted IFN vs PCA vs Geneformer, 2026-11-22
Pre-registered in PREREG-Q13-qc-adjusted.md (commit 893e1d0) before any adjusted feature was computed.
Same 2,200 monocytes as round 92 (Base vs M3 after CD3xCD19, 11 donors, LODO). Unadjusted numbers reproduced
exactly (IFN 0.680, PCA-50 0.817, Geneformer 0.772, QC-only 0.661).

| Feature | raw | A1 QC regressed out | A2 QC + feature | A3 QC-matched cells (8 donors) |
|---|---|---|---|---|
| QC only (3 covariates) | 0.661 | - | - | 0.489 |
| F1 IFN6 score | 0.680 | 0.617 | 0.713 | 0.653 |
| F2 PCA-50 | 0.817 | 0.780 | 0.818 | 0.688 |
| F3 Geneformer-V2 | 0.772 | 0.720 | 0.772 | 0.685 |

- H1 IFN survives QC: PASS. Residualised 0.617 >= 0.60; adds +0.051 over QC alone, up in 9/11 donors.
- H2 Geneformer after A1: still DROPPED (+0.104 over IFN, but beats IFN in only 6/11 donors and trails PCA by 0.06).
- H3 PCA edge technical? By the pre-registered A1 rule: NO (edge 0.164 after regression vs 0.137 raw).
  But the secondary matched analysis cuts the PCA edge over IFN to 0.035 and Geneformer's to 0.032, while matching
  removes the QC signal entirely (0.49). Linear regression on 3 covariates does not remove non-linear technical
  structure from 50 dimensions; matching does. Reading: most of the high-dimensional edge is technical; it is
  not established which.
- Donor P12 inverts for PCA (0.10-0.16) and Geneformer (0.20-0.29) under every adjustment; IFN does not (0.56-0.60).
Files: q13_lodo_auc.csv, q13_matched_auc.csv, q13_verdict.json, q13_result.txt; code evaluate_q13.py.
