# RESULT round 100 (Q35, C1): scVI does not beat plain PCA on A-319 baseline vs month 3. Both survive QC matching at donor level.

PREREG: PREREG-C1-round100.md (24bc78f), deviation logged before evaluation (e0b2fd8). Code: c1_prep.py, c1_aggregate.py, c1_latent.py, c1_eval.py. Output: c1_result.json.
Data: CELLxGENE 6c147b03 (A-319 CD3xCD19 phase 1), 10 donors x {Base, M3} = 20 samples (P7 dropped, see deviation), 155,771 cells.
QC: 2,850 cells removed (1.8%), 17 doublets; Base 1,519, M3 1,331.

| Representation | Donor-level LODO AUC | Paired accuracy | Sign-flip p (1,024) |
|---|---|---|---|
| B0 QC-only | 0.74 | 0.9 | 0.029 |
| B1 IFN6 (monocytes) | 0.61 | 0.7 | 0.070 |
| B2 CLR composition | 0.81 | 0.9 | 0.002 |
| B3 pseudobulk (6 types) | 0.57 | 0.5 | 0.24 |
| **B4 PCA-50 (fold-fit), per-type means** | **0.94** | 0.9 | 0.002 |
| D1 scVI (fold-fit, scArches query) | 0.91 | 0.9 | 0.002 |
| D1 + B0 | 0.91 | 0.9 | 0.004 |

Decision rules:
1. D1 - best simple (B4) = -0.03 (95% donor bootstrap -0.15 to +0.08). **FAIL.**
2. D1 - B0 = +0.17; D1+B0 > B0. Pass.
3. QC matching (14,240 matched cells): B0 falls to 0.43 (check met). B4 0.87, D1 0.93. Pass as coded. There was no rule-1 gain to "survive", though.
4. Not applicable (each library = one donor x timepoint).

**H1: FAIL.** scVI adds nothing over a fold-fit PCA of the same cells.

Notes:
- New vs round 98: once cells are summarised per donor and type, both latent spaces keep their signal after QC matching (PCA 0.94 -> 0.87, scVI 0.91 -> 0.93). In round 98 the PCA edge collapsed after matching, but that was cell-level classification in monocytes only. At donor level the latent signal is mostly not technical here. B0 alone did reach 0.74 before matching.
- The IFN score is weak (0.61), consistent with round 96 (shared multiplier 0.73, small fall).
- Pseudobulk loses (0.57). Likely p >> n: 12,000 features reduced to 10 PCs on 18 training samples.
- Protected-population check: the QC flag rate was 7.7% in B/plasma and 11% in platelet/other, vs 0.1% in T/NK. The PREREG said to keep coherent flagged populations, but the code removed all flagged cells. This is a DEVIATION; the 11-donor rerun will apply the rule.
- scVI hit the 200-epoch cap in every fold (early stopping never triggered).
- Limits: n = 10 pairs; CIs are wide. Lineage files are author-annotated; neutrophils were absent from some samples' latent summaries.
