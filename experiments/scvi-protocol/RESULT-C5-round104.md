# Round 104 result (Q37, contrast C5): frozen Perez flare model on GSE254176, plus the bridging-donor batch check

PREREG-C5-round104.md (commit 36b5ae0), amendments Q43 (B5) and Q45 (rule 3 on all B0 metrics) in the same commit; scripts 9936c5a.
Data: Perez 2022 CELLxGENE 218acb0f cohort 3 (10 donors x flare/treated, 51,483 cells; 11,236 subsampled), 21,578 genes shared with
GSE254176 (19 PBMC samples, 6 patients, 119,843 cells with >= 200 genes, 117,746 kept). Full numbers: c5_result.json.

## H1 (learned latent beats simple baselines on transfer): FAIL on all four rules

| Rep | AUC flare vs not (19 samples) | patient-bootstrap 95% CI | Flare > Before (of 6) | within-patient AUC | domain-standardised AUC |
|---|---|---|---|---|---|
| B0 QC only | 0.42 | 0.12-0.75 | 3 | 0.35 | 0.37 |
| B1 IFN6 (monocytes; not blind) | 0.47 | 0.19-0.70 | 2 | 0.39 | 0.47 |
| **B2 composition (CLR)** | **0.74** | **0.52-1.00** | 4 | **0.77** | **0.81 (6/6)** |
| B3 pseudobulk | 0.60 | 0.42-0.81 | 5 | 0.79 | 0.58 |
| B4 PCA-50 | 0.56 | 0.33-0.82 | 3 | 0.62 | 0.65 |
| B5 nearest-context copy (k 3) | 0.39 | 0.13-0.72 | 2 | 0.31 | - |
| D1 scVI + scArches | 0.53 | 0.35-0.86 | 2 | 0.36 | 0.58 |
| D2 sysVI (transductive) | 0.36 | 0.14-0.89 | 3 | 0.29 | 0.64 |
| D1 + B0 | 0.53 | 0.23-0.86 | 2 | 0.32 | - |

- Rule 1: D1 - B2 = -0.22 (CI -0.33 to 0.00). FAIL. D2 - B2 = -0.39 (CI -0.46 to -0.11). FAIL.
- Rule 2: D1 - B0 = +0.10 and D1+B0 (0.53) > B0 (0.42): passes, but only because B0 is below chance.
- Rule 3 (all 9 QC metrics, 57,041 of 117,746 cells matched): B0 rises to 0.63 (> 0.6), so by amendment Q45 rule 3 FAILS
  ("QC metrics carry biology or batch that matching cannot remove"). D1 0.67 vs best matched baseline B4 0.54 is reported, not claimed.
  Composition is not interpretable after matching (matching removes cells unevenly).
- Rule 4: scVI latent predicts patient (skill 0.28) better than flare (0.15). FAIL.
- Descriptive transfer bar (AUC >= 0.70 and >= 5/6 episodes): no representation passes at the pre-registered default. B2 passes only
  in the pre-declared domain-standardised secondary (0.81, 6/6).
- Training in-sample AUC (all 20 Perez samples): B2 0.92, D1 0.90, D2 0.86, B3 0.89; so the latents fit Perez and do not carry over.

## Bridging-donor check (descriptive): scVI keeps donor identity across Perez cohorts; sysVI over-corrects

122 donor x cohort samples of 57 donors (4 processing cohorts), types T4, T8, cDC, cM, ncM; 146 retrieval queries.

| Rep | Same donor ranked first in other cohort | median rank | cohort-probe skill |
|---|---|---|---|
| PCA-50 | 0.34 | 3 | 0.91 |
| scVI (library as batch) | **0.78** | 1 | 0.34 |
| sysVI (cohort as system) | 0.22 | 5 | 0.16 |
| composition | 0.60 | 1 | 0.00 |

## Declared implementation details and deviations
- scArches query fine-tuned on a label-blind subsample of <= 1,500 cells per GSM (26k cells), then all cells encoded with the query model.
- Cell types transferred by kNN (k 15) from Perez author labels; median kNN vote 0.93. Transferred CD8 T fractions (24-43%) and
  plasmablast fractions (0.6-6.9%) are high for PBMC, so label transfer across chemistry is imperfect; composition results inherit this.
- GSE protected types (flag rate > 3x): PB (74% flagged), Prolif; kept.
- Only T4, T8, cM met the per-sample minima for latent and pseudobulk features (as in round 102).
- Bridging HVG selection used the span 1.0 fallback (same fallback as round 102). sysVI hit max_epochs 200 without early stopping.
- A bug (C.sub) stopped the first sysVI call after the scVI/annotation outputs were saved; sysVI was rerun from the saved files
  (step c5sysvi), identical inputs.
