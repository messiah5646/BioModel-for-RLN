# Result: Geneformer-V2 vs the 6-gene IFN score (queue Q3), 2026-11-21

Pre-registered in PREREG.md (commit 77f694f) before any feature was computed.
Task: classify monocytes as after (M3) or before (Base) CD3xCD19 B-cell depletion; leave-one-donor-out, 11 donors,
2,200 cells (100 per sample).

| Feature set | Mean held-out AUC |
|---|---|
| F0 chance | 0.500 |
| F1 6-gene IFN score (1 feature) | 0.680 |
| F3 Geneformer-V2-104M, frozen mean embedding (768) | 0.772 |
| F2 PCA-50 of 2,000 HVGs | 0.817 |
| F4 QC only: log counts, log genes, % mito (EXPLORATORY, not pre-registered) | 0.661 |

Pass rule: (1) F3 - F1 >= 0.05: 0.092, met. (2) F3 > F1 in >= 8/11 donors: 7/11, NOT met.
(3) F3 >= F2 - 0.02: 0.772 vs 0.817, NOT met.
Verdict: Geneformer is DROPPED for this axis. A plain PCA of the same cells beats it.

Reading:
- The IFN score carries held-out signal (0.68), consistent with the 8/11 donor-level fall.
- Three cheap technical covariates alone reach 0.66, about the same as the IFN score. M3 cells have higher % mito
  (median 1.55 vs 1.04). So much of what high-dimensional features learn here may be sample processing, not biology.
- P12 is reversed for both PCA (0.10) and Geneformer (0.20): both learned a direction that holds across the other
  donors and inverts in one, the signature of a technical axis.
- The question "does a foundation-model embedding add biology beyond the IFN score" is not answerable on this design;
  it needs within-sample labels (e.g. a perturbation applied in one library) or batch-balanced sampling.

Files: lodo_auc.csv, verdict.json, exploratory_qc_auc.csv; code prep.py, embed.py, evaluate.py, exploratory_qc.py.
Data: CELLxGENE 3979a248 (re-downloadable), Geneformer from Hugging Face ctheodoris/Geneformer.
