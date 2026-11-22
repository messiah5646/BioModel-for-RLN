# PREREG addendum, round 100 (Q35, contrast C1): A-319 CD3xCD19, baseline vs month 3

Frozen 2026-11-22 before any feature below was computed. Parent protocol: experiments/scvi-protocol/PREREG.md (sha 49ca653a81dc).
Seen so far: file shapes, obs column names, timepoint / cell_type counts per file (design only). No expression, score or classifier output.

## Scope of this round
C1 only. C2 (CAR-T GSE263931) needs de-novo cell annotation of raw 10x matrices and moves to its own round (Q35b).

## Data
CELLxGENE 6c147b03, 10x 5' v2, 11 donors, timepoints Base and M3 (D28 excluded). Files: T (00260ebe), B (65c2c3c1), Myeloid (0f689b40).
The Monocytes (d0bce338) and LDG (729e31d1) files are re-clusterings of Myeloid cells; they are used only if their barcodes are NOT
contained in Myeloid (checked by obs index overlap; reported). Counts = `raw.X`.
Major types (from `cell_type`): CD4_T, CD8_T, NK, other_T, B_plasma, monocyte, neutrophil, DC, platelet_other.

## Stage 1 QC (per sample)
5 MAD on log1p counts, log1p genes, % counts in top-20 genes (both sides); 3 MAD upper on logit mito AND mito > 10% (deposited ceiling).
Doublets: scanpy scrublet per sample, flag only (DEVIATION from scDblFinder, which is R-only; declared now).
Ambient: measure HB (HBA1/HBA2/HBB), platelet (PPBP/PF4) and IG (IGHG*/IGKC/IGLC*/JCHAIN) fractions in non-erythroid, non-B cells; measure only.
Stress (FOS, FOSB, JUN, JUNB, HSPA1A, HSPA1B, DUSP1, IER2) and cell cycle scored, never regressed.
Protected-population check: test enrichment of QC flags in plasmablast/plasma (MZB1/JCHAIN) and neutrophils; any coherent flagged population is kept.
Cells failing a MAD flag are removed for all representations; doublet-flagged cells are removed. Removal counts are reported by sample x timepoint.

## Representations (per donor-timepoint sample; n = 22)
- B0 QC-only: median log counts, log genes, mito %, ribo %, IG/HB ambient fractions, doublet score, stress score, log cells recovered.
- B1 6-gene IFN score (IFI27, IFI44L, IFIT1, ISG15, RSAD2, SIGLEC1; mean log-normalised expression) averaged over monocytes.
- B2 CLR composition over major types (pseudocount 0.5).
- B3 pseudobulk log-CPM per major type (types with >= 20 cells in every sample), 2,000 HVGs.
- B4 PCA-50 of log-normalised 2,000 HVGs (fit inside fold on training cells), sample mean per major type.
- D1 scVI (n_latent 20, n_layers 2, layer norm, NB, early stopping, max_epochs 200, batch_key donor_id, no continuous covariates, 2,000 HVGs from training cells), trained inside each fold on training donors; held-out donor mapped with scArches (query epochs 100, weight_decay 0). Sample mean latent per major type.
- D3 Geneformer: NOT run this round (DEVIATION, declared now: CPU-only and already lost to PCA twice, rounds 92/98).
Cells for B4/D1: at most 150 per sample x major type, seed 0 (within the parent cap of 300).

## Classifier and metrics
High-dimensional sets (B3, B4, D1) are reduced to 10 PCs fitted on training samples only. Logistic regression (L2, C = 1, standardised) for every set.
Leave-one-donor-out (train 20 samples, predict both samples of the held-out donor).
Primary metric: AUC over the 22 out-of-fold probabilities (M3 = positive). Also: paired accuracy (held-out donor's M3 probability > Base) with exact sign-flip null (2^11).
CI: donor bootstrap (2,000) of AUC differences.

## Decision rules (parent protocol, all four needed for H1)
1. D1 - max(B1..B4) >= 0.05 AUC and bootstrap 95% CI excludes 0.
2. D1 - B0 >= 0.05 and AUC(D1 + B0) > AUC(B0).
3. QC-matched rerun: within donor x major type, 1:1 nearest-neighbour match of M3 to Base cells on standardised (log counts, log genes, mito), caliper 0.5 SD; recompute B0, B4, D1 (D1 reuses the fold's latent). Mandatory check: B0 AUC falls to <= 0.6. The D1 gain must survive.
4. Batch probe: NOT APPLICABLE here (each library is exactly one donor x timepoint, so library cannot be separated from label); declared now.

## Expected (literature + rounds 92/98)
B2 composition near ceiling (B cells are depleted by design), so a deep gain >= 0.05 is unlikely. A FAIL is the expected outcome.

## Deviation logged 2026-11-22 12:15 UTC, BEFORE any evaluation output was seen
Sample S08_M3 (donor P7, month 3) was silently killed by the 7 GB memory limit during per-sample QC (largest myeloid sample), so no
QC/latent data exist for it. Donor P7 is therefore dropped from evaluation (10 donor pairs, 20 samples; sign-flip null 2^10 = 1,024).
P7's baseline cells stayed in the scVI/PCA training sets of the other folds. A full 11-donor rerun is queued.
