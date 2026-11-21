# Pre-registration: Geneformer-V2 embeddings vs the 6-gene IFN score (queue Q3), 2026-11-21

Written after reading obs metadata only (sample, donor, timepoint, cell type counts) and a speed benchmark on random
token ids. No expression-derived feature, score or label association has been computed.

## Data
CELLxGENE collection 6c147b03 (CD3xCD19 bispecific phase 1 in SLE), dataset 3979a248 "Monocytes"
(76,421 cells, 11 donors, 10x 5' v2). Cells: all annotated monocytes (IFI+, Classical, CD16+ Non-Classical);
granulocytes excluded. Samples: Base and M3 only (D28 and M13 excluded). 11 donors x 2 samples.
Random subsample of up to 100 cells per sample (numpy seed 0), the same cells for every feature set.

## Task
Per-cell binary classification: M3 (after B-cell depletion) = 1 vs Base = 0.
Leave-one-donor-out (11 folds). Classifier for every feature set: StandardScaler + LogisticRegression
(L2, C = 1.0, class_weight = balanced, max_iter 5000), fitted on the 10 training donors only.
Metric: AUC on the held-out donor's cells (each held-out donor has both classes). Summary = mean over 11 donors.

## Feature sets
- F0 chance: AUC 0.5 by definition (reported as the floor).
- F1 IFN score: mean of CELLxGENE log-normalised X over IFI27, IFI44L, IFIT1, ISG15, RSAD2, SIGLEC1 (1 feature;
  same score as rounds 89 and the ifn-cart test).
- F2 PCA-50: top 2,000 highly variable genes (by variance of log-normalised X on training donors), PCA with 50
  components fitted on training donors only, applied to the held-out donor.
- F3 Geneformer: Geneformer-V2-104M (Hugging Face ctheodoris/Geneformer), frozen, no fine-tuning. Tokenisation as in
  the Geneformer tokenizer: raw counts / cell total x 10,000 / gene median (gc104M dictionary), non-zero genes ranked
  descending, <cls> prepended, <eos> appended, truncated to 2,048 tokens (deviation from the 4,096 limit, for CPU
  time). Cell embedding = mean over gene tokens (excluding <cls>/<eos>) of the second-to-last hidden layer
  (Geneformer's default emb_layer = -1). 768 features.

## Pass rule (SBB rule: an embedding is kept only if it beats the simple score)
Geneformer is KEPT only if all three hold:
1. mean AUC(F3) - mean AUC(F1) >= 0.05;
2. AUC(F3) > AUC(F1) in >= 8 of 11 held-out donors;
3. mean AUC(F3) >= mean AUC(F2) - 0.02 (it must not be worse than a plain PCA of the same cells).
Otherwise Geneformer is DROPPED as a feature for this axis.
Secondary, reported, no threshold: F1 vs F0 (does the IFN score itself carry held-out signal; expected yes from the
8/11 donor-level fall), and F2 vs F1.

## Known confound, stated in advance
Each sample is its own library; Base and M3 were collected months apart. High-dimensional features (F2, F3) can
learn technical sample differences. Leave-one-donor-out stops donor-specific batch from leaking, but a systematic
Base-vs-M3 processing difference would favour F2 and F3 over F1. A F3 or F2 win is therefore not by itself evidence of
biology.
