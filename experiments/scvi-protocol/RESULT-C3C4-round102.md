# RESULT round 102 (Q36, C3 + C4): scVI shows no demonstrated gain on Perez flare-vs-treated or SLE-vs-control

PREREG: PREREG-C3C4-round102.md (d169548); C4 HVG deviation logged before evaluation (b229d05). Code: p_design.py, p_prep.py, p_latent.py, p_eval.py. Output: c3_result.json, c4_result.json.
Data: CELLxGENE 218acb0f (Perez 2022, 10x 3' v2), raw.X streamed. C3: cohort 3, 10 donors x {flare, treated}, 51,483 cells. C4: cohort 4, 52 SLE vs 44 controls, 375,261 cells.
QC: C3 removed 150 cells (0.3%), C4 1,196 (0.3%); protected-population rule kept 157 / 630 flagged PB, Prolif, Progen, cDC cells (flag rates 3-60x baseline). Latent/pseudobulk types meeting the per-sample minimum: T4, T8, cM only.

## C3 flare vs post-treatment (10 pairs, leave-one-donor-out)
| Set | AUC | Paired acc | Sign-flip p |
|---|---|---|---|
| B0 QC | 0.33 | 0.3 | 0.97 |
| B1 IFN6 (monocytes) | 0.26 | 0.3 | 0.97 |
| **B2 composition** | **0.75** | 0.6 | 0.025 |
| B3 pseudobulk | 0.48 | 0.4 | 0.72 |
| B4 PCA-50 | 0.31 | 0.4 | 0.95 |
| D1 scVI | 0.50 | 0.5 | 0.44 |
| D1+B0 | 0.55 | 0.8 | 0.47 |
Rule 1 FAIL (D1 - B2 = -0.25, CI -0.58/+0.11). Rule 2 pass as coded (+0.17 over B0, which is below chance). Rule 3 FAIL (B0 0.53 after matching, check met; D1 0.50). Rule 4 pass (label skill 0.23 vs library 0.04). **H1 FAIL.**
Descriptive: IFN6 higher at flare in 6/10 donors. Flare vs treated composition: B cells up in 7/10, proliferating cells up 9/10, plasmablasts up 8/10, CD8 T down 9/10.

## C4 SLE vs control (96 donors, 5-fold donor CV)
| Set | AUC seed-0 | 20-repeat mean (sd) | Perm p |
|---|---|---|---|
| B0 QC | 0.86 | 0.86 (0.01) | 0.001 |
| B1 IFN6 | 0.76 | 0.75 (0.01) | 0.001 |
| B2 composition | 0.85 | 0.85 (0.01) | 0.001 |
| **B3 pseudobulk** | **0.91** | 0.92 (0.01) | 0.005 |
| B4 PCA-50 | 0.87 | - | 0.005 |
| D1 scVI | 0.93 | - | 0.005 |
| D1+B0 | 0.94 | - | - |
Rule 1 FAIL (D1 - B3 = +0.014, CI -0.033/+0.065). Rule 2 pass (+0.07 over B0; D1+B0 0.94 > 0.86). Rule 3 FAIL: after matching on counts/genes/mito within library x type, B0 ROSE to 0.93 (mandatory check <= 0.6 not met); D1 0.92, B4 0.85. Rule 4 pass (label 0.34 vs library 0.02). **H1 FAIL.**
Strict QC (3 MAD): B0 0.86, B1 0.76, B2 0.86 (unchanged).

## Unplanned descriptive check (no pass rule): what makes the "QC-only" set work in C4
Univariate donor-level AUC (SLE positive): ribosomal-protein fraction 0.16 (i.e. 0.84 lower in SLE), also 0.17 within CD4 T cells alone; log genes 0.69; mito 0.39; IG/HB 0.50. Ribo fraction was not a matching variable, which is why matching did not remove B0's signal. In C3 (flare vs treated) no QC metric separates (all 0.42-0.55).
