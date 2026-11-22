# Pre-registration, round 98 (queue Q13): does any feature survive QC adjustment? 2026-11-22

Written before any QC-adjusted feature or AUC has been computed. Known from round 92 (unadjusted): IFN 0.680,
Geneformer 0.772, PCA-50 0.817, QC-only (exploratory) 0.661; M3 cells have higher % mito.

## Data and cells
Identical to PREREG.md (round 92): CELLxGENE 3979a248 Monocytes, Base vs M3, up to 100 cells per sample, numpy
seed 0 (prep.py, unchanged), 11 donors, 2,200 cells. Task, classifier (StandardScaler + L2 logistic, C=1,
balanced) and leave-one-donor-out (LODO) AUC unchanged. Features F1 IFN6, F2 PCA-50 of 2,000 HVGs, F3 Geneformer-V2
frozen mean embedding, all built exactly as in round 92.
QC covariates Q = [log1p nCount_RNA, log1p nFeature_RNA, percent.mito] from CELLxGENE obs.

## Adjustments (all fitted on the 10 training donors only, applied to the held-out donor)
- A1 regress-out: each feature column is replaced by its residual from an OLS fit on Q (intercept + 3 QC terms).
  For F2, PCA is fit first (as in round 92), then the 50 PC scores are residualised.
- A2 incremental: model on [Q + feature] compared with model on Q alone (F4).
- A3 QC-matched subsample (secondary): within each donor, each M3 cell is matched 1:1 without replacement to its
  nearest Base cell on standardised Q (greedy, caliper 0.25 SD Euclidean); unmatched cells dropped; LODO on matched
  cells; donors with < 20 matched pairs reported but excluded from the mean.

## Hypotheses and pass rules
- H1 (IFN survives QC): PASS if residualised IFN (A1) has mean LODO AUC >= 0.60 AND [Q+IFN] beats Q alone (A2) by
  >= 0.03 mean AUC with improvement in >= 7 of 11 donors.
- H2 (Geneformer, same SBB rule, after A1): KEPT only if resid-F3 - resid-F1 >= 0.05, resid-F3 > resid-F1 in >= 8/11
  donors, and resid-F3 >= resid-F2 - 0.02. Otherwise stays DROPPED.
- H3 (was the PCA edge technical?): "PCA edge largely technical" if (F2 - F1) after A1 <= half of the unadjusted
  (F2 - F1) = 0.137, i.e. <= 0.069. Otherwise the high-dimensional edge survives QC adjustment.
- Trivial baselines reported for every comparison: chance 0.5 and Q alone (F4).

## Caveats stated in advance
Regressing out % mito and depth can remove real biology (activation changes RNA content). A1 is conservative for
every feature; a feature that fails A1 but passes A2 is reported as "partly confounded", not as noise.
A3 changes the cell set; it is secondary and has no pass rule.
Research only.
