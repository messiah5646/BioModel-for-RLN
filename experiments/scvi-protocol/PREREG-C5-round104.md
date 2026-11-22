# PREREG addendum, round 104 (Q37, contrast C5): frozen Perez cohort-3 flare model, single pass on GSE254176

Frozen 2026-11-22 before any GSE254176 expression value is loaded in this round. Parent: PREREG.md (sha 49ca653a81dc), plus the
round-104 amendments AMENDMENT-Q43-VCC-baseline-round104.md (B5) and AMENDMENT-Q45-rule3-allB0-round104.md (rule 3 on all B0 metrics).

## What has already been seen (declared, not blind)
- GSE254176 labels (19 scRNA samples, 6 patients: 6 Flare, 6 Before flare, 7 Remission; P5 and P6 have one Remission sample only)
  and matrix sizes (2,531-11,794 barcodes per sample; 36,601 features, Cell Ranger 2020-A reference, so likely 10x 3' v3; Perez is 3' v2).
- Round 89 already scored monocyte IFN6 and plasmablast fraction on this dataset (IFN Flare > Before in 2/6 episodes). B1 is therefore
  NOT blind on this dataset. All other representations are blind.
- Perez cohort 3 C3 results (round 102): within-Perez LODO AUC scVI 0.50, composition 0.75, IFN6 0.26, PCA 0.31, B0 reported there.
  The training side is therefore not blind either; only the transfer is new.

## Training (Perez cohort 3, frozen once)
- Same 20 samples (10 donors x flare/treated), same Stage-1 QC and protected-population rule as PREREG-C3C4-round102.md (code p_prep.py logic).
- Subsample <= 100 kept cells per sample x author_cell_type (seed 0). Genes = intersection of Perez feature_name and GSE254176 symbols.
- Classifier per representation: L2 logistic regression (C = 1), standardised on the 20 training samples, fitted ONCE on all 20.
  High-dimensional sets (B3, B4, D1, D2) reduced to 10 PCs fitted on the training samples.

## Test (GSE254176, single pass)
- Stage 1 QC per sample: 5 MAD on log counts, log genes, top-20 %; 3 MAD upper on logit mito AND mito > 15%; cells < 200 genes dropped first.
  Doublets not scored (as round 102, declared). Protected-population rule applied with transferred labels.
- Cell types: transferred from Perez author_cell_type by kNN (k = 15, majority) in the D1 scVI latent after scArches mapping
  (each GSM sample is a new batch). Annotation is frozen before any label-dependent output is computed.
- Subsample <= 100 kept cells per sample x type (seed 0) for latent features. Types used = those with >= 10 subsampled cells in every
  sample of BOTH datasets (latent) or >= 20 kept cells in every sample (pseudobulk); a type that fails is dropped from both sides.
- Representations (all exactly as round 102 unless stated): B0 QC medians; B1 IFN6 mean over cM+ncM; B2 CLR composition (pseudocount 0.5);
  B3 pseudobulk log-CPM, 2,000 top-variance genes per type chosen on Perez; B4 PCA-50 on 2,000 batch-aware HVGs fit on Perez cells,
  GSE cells projected; B5 nearest-context copy (amendment Q43, k = 3 on standardised B2);
  D1 scVI (n_latent 20, n_layers 2, layer norm, NB, early stopping, max_epochs 200, encode_covariates, dropout 0.2, batch_key library_uuid)
  trained on Perez C3 cells only, GSE cells mapped by scArches (100 epochs, weight_decay 0);
  D2 sysVI (scvi.external.SysVI, defaults, system = dataset, batch covariate = library/GSM) trained on Perez C3 + GSE254176 cells
  WITHOUT labels: this is TRANSDUCTIVE (labelled as such); its classifier still sees Perez labels only.
- Primary metric: AUC, Flare (6) vs non-Flare (13), over the 19 samples, score = frozen classifier probability.
- Secondary: (a) paired episodes, Flare > preceding Before-flare in k of 6 (P1 1->2, P2 1->2, P3 1->2, P3 4->5, P4 2->3, P4 5->6);
  (b) within-patient-centred AUC (score minus that patient's mean; P5, P6 excluded, 17 samples);
  (c) domain-standardised variant: each dataset's features z-scored on its own samples (label-blind, transductive) before the frozen classifier.
- CI: patient-cluster bootstrap (2,000 resamples of the 6 patients) for AUC and AUC differences.

## Decision rules (all four for H1; parent rules with C5 specifics)
1. D1 - max(B1..B5) >= 0.05 AUC and bootstrap CI excludes 0 (D2 reported against the same rule, as secondary).
2. D1 - B0 >= 0.05 and AUC(D1+B0) > AUC(B0).
3. QC matching (amendment Q45): within patient x type, each Flare cell matched 1:1 (nearest neighbour, without replacement) to a non-Flare
   cell of the same patient on standardised log counts, log genes, top-20 %, mito, ribo, stress, HB, platelet, IG; caliper 0.5 SD x sqrt(9).
   P5/P6 keep all cells. Recompute B0, B4, D1 test features (latents reused), apply the frozen classifiers. B0 must fall to <= 0.6 AUC;
   the D1 gain over max(B4_matched, best non-latent baseline) must stay >= 0.05.
4. Batch probe on GSE cells in the D1 latent: 5-fold kNN (k = 15); skill = (balanced accuracy - chance)/(1 - chance) for PATIENT vs for
   Flare label. (Library = sample here, so patient is the batch-like nuisance.) Pass if label skill >= patient skill.

Separate, descriptive transfer bar (no H1 implication): a representation "transfers" only if AUC >= 0.70 AND paired >= 5/6.

## Bridging-donor batch check (descriptive, no pass rule; moved here from round 102)
Perez donors sampled in > 1 Processing_Cohort: <= 30 kept cells per sample x type. Representations: PCA-50 (HVGs), scVI (batch_key
library_uuid), sysVI (system = Processing_Cohort, batch = library_uuid), B2 composition. Per sample: concatenated per-type mean latent.
Metric: for each bridging sample, rank of the same donor's sample from a different cohort among all samples of that other cohort
(cosine distance): top-1 accuracy and median rank; plus a cohort-probe skill (kNN on per-sample vectors). Higher retrieval = better
batch correction that keeps donor identity.

## Expected
FAIL. C3 failed inside Perez (scVI 0.50), the chemistry differs (v2 vs v3), "treated" (post-treatment) is not the same state as
"Before flare"/"Remission", and round 89 showed blood IFN is a patient set point. A transfer AUC near 0.5 is the prior.
