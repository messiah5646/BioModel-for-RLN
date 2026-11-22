# PREREG addendum, round 102 (Q36, contrasts C3 and C4): Perez 2022 (CELLxGENE 218acb0f)

Frozen 2026-11-22 before any expression-derived feature was computed. Parent protocol: PREREG.md (sha 49ca653a81dc).
Seen so far: obs only (cohort, disease_state, donor, library, author cell type counts; output of p_design.py). No expression, score or classifier output.

## Design facts (from obs)
- C3: cohort 3, 10 donors with both a flare and a post-treatment ("treated") sample. Flare and treated samples of the same donor are NEVER in the same library (0/10). Positive class = flare. n = 20 samples.
- C4: cohort 4, 52 SLE (managed) vs 44 controls; 24 libraries (pooled wells). Sample unit = donor (replicate samples of one donor inside cohort 4 are merged). Positive class = SLE.
- 57 donors appear in >1 processing cohort (bridging donors).

## Stage 1 QC (per library_uuid x donor), from raw.X
- 5 MAD on log1p counts, log1p genes, top-20 %; 3 MAD upper on logit mito AND mito > 15%.
- Doublets: NOT re-scored (DEVIATION, declared now). The open matrix is already post Freemuxlet + Scrublet; B0 therefore has no doublet score.
- Ambient (HB, platelet PF4/PPBP, IG) and stress (FOS/FOSB/JUN/JUNB/HSPA1A/HSPA1B/DUSP1/IER2) measured, never regressed.
- Protected-population rule (fixes the round-100 deviation): if an author cell type's flag rate is > 3x the cohort's overall flag rate, flagged cells of that type are KEPT (reported).
- Strict setting (3 MAD on all metrics) is reported for B0-B2 only (DEVIATION: D1/B4 not refitted under strict, CPU budget).

## Representations (parent table)
- Major types = author_cell_type (T4, T8, NK, B, PB, Prolif, cM, ncM, cDC, pDC, Progen); latent / pseudobulk use types with >= 10 subsampled cells (latent) or >= 20 kept cells (pseudobulk) in every sample.
- B0 QC medians (log counts, log genes, mito, ribo, IG and HB in non-B, stress, log cells). B1 IFN6 mean over cM+ncM. B2 CLR composition (pseudocount 0.5). B3 pseudobulk log-CPM, 2,000 top-variance genes per type chosen on training samples. B4 PCA-50 on 2,000 HVGs fit on training cells. D1 scVI exactly as round 100 (n_latent 20, n_layers 2, layer norm, NB, early stopping, max_epochs 200, encode_covariates, dropout 0.2) with batch_key = library_uuid; held-out samples mapped by scArches (100 epochs, weight_decay 0). D3 Geneformer not run (as round 100).
- Subsample (fixed seed 0) per sample x type: C3 at most 100 cells; C4 at most 30 cells (DEVIATION from <= 300 cap is in the allowed direction; CPU budget).
- High-dimensional sets (B3, B4, D1) reduced to 10 PCs fitted on training samples. Logistic regression L2 C=1, standardised.

## Validation
- C3: leave-one-donor-out (10 folds; scVI refitted per fold). AUC over 20 out-of-fold probabilities; paired accuracy; exact sign-flip null (2^10).
- C4: 5-fold donor-grouped stratified CV. Baselines B0-B3: 20 repeats (seeds 0-19). B4/D1: latents fitted for ONE split (seed 0, 5 folds) then the classifier is re-run over 20 repeats of the outer split *using fold-0 latents transductively* is NOT allowed; instead B4/D1 are reported on the seed-0 split only (DEVIATION: 1 repeat, not 20). Label-permutation null (1,000): latents are unsupervised (no labels seen), so permutation re-fits only the classifier on fixed fold latents; this is exact for the supervised step.
- CI: donor bootstrap (2,000) of AUC differences, on the seed-0 split.

## Decision rules (all four for H1, per contrast)
1. D1 - max(B1..B4) >= 0.05 AUC and bootstrap CI excludes 0.
2. D1 - B0 >= 0.05 and AUC(D1+B0) > AUC(B0).
3. QC matching: within donor-pair x type (C3: flare vs treated of the same donor; C4: within library x type, SLE vs control cells) 1:1 nearest-neighbour on standardised (log counts, log genes, mito), caliper 0.5 SD x sqrt(3); recompute B0/B4/D1 (latents reused). B0 must fall to <= 0.6 AUC; the D1 gain must survive (D1_matched - max(B4_matched, best non-latent baseline) >= 0.05).
4. Batch probe (C4, and C3 where computable): on training cells of fold 1, 5-fold kNN (k = 15) on the D1 latent; skill = (balanced accuracy - chance)/(1 - chance) for library vs for label. Pass if label skill >= library skill.

## Bridging-donor batch check (descriptive, no pass rule)
For the bridging donors present in cohort 3 or 4 and another cohort, NOT done this round (would need a third cohort's latent; DEVIATION, moved to Q37 where cross-cohort transfer is the question).

## Q43 interaction
Q43 (VCC "copy nearest measured response" baseline) targets response prediction, not classification; it is not applicable to C3/C4 and stays a prerequisite for C5 (Q37).

## Expected
C4: composition and IFN near ceiling (SLE vs control is strong in Perez), so a deep gain >= 0.05 is unlikely. C3: n = 10 pairs, flare/treated each in different libraries, so library effects can mimic the label; FAIL expected.

## Deviation logged 2026-11-22 ~14:20 UTC, BEFORE any evaluation output was seen
C4 fold f0: batch-aware seurat_v3 HVG selection failed (loess near-singularity in small per-library cell sets, 30 cells/type cap).
Fallback: same call with span = 1.0 when the default span raises. C3 latents were all computed with the default span.
