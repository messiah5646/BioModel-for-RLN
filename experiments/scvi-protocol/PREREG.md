# PREREG (rounds 99+): scRNA QC + scVI vs simple baselines on open lupus data

Frozen 2026-11-22 from the methods report (reports/scRNA QC and scVI for lupus.md, notes in research_notes/), BEFORE any of these representations were computed on C1-C5/K1. Thomas asked for 'both': methods report, then apply it. Seed 0; scVI config: n_latent 20, n_layers 2, layer norm, gene_likelihood nb, early stopping, no continuous covariates, 2000 batch-aware HVGs. Hardware note: 2 CPUs, 7 GB RAM, no GPU -> Perez is streamed per cohort/library (backed h5ad), cells subsampled to <=300 per donor x major cell type (fixed seed) and the subsample is reported.

The protocol tests one hypothesis family: **H1, a learned latent space (scVI, sysVI where justified, Geneformer as negative control) gives a donor-level classification gain over the best simple baseline that survives QC control**. The null is the pattern seen in round 98, where the gain vanishes after QC matching. The tables below fix the settings in advance. Every threshold is a pre-registered default, and any change after unblinding is reported as a deviation.

**Stage 1: QC (identical code for all datasets, run per sample or library)**

| Step | Setting | Perez 2022 | CAR-T GSE263931 | A-319 | GSE254176 | Kidney (KPMP ref + E-MTAB-13596) |
|---|---|---|---|---|---|---|
| Start matrix | raw counts | `raw.X` (already post-QC) | filtered MTX | `raw.X`, 5 lineage files concatenated | filtered MTX | `raw.X` / MTX |
| MAD flags | 5 MAD on log counts, log genes, top-20 %; 3 MAD upper on logit mito | per `library_uuid` × donor | per sample | per sample (already floored) | per sample | per library, then per broad compartment |
| Mito hard ceiling | flagged only if above 3 MAD **and** above this ceiling | 15% | 15% | 10% (deposited ceiling) | 15% | 50% for scRNA (KPMP); MT genes dropped for snRNA |
| Ambient | measure HB, platelet and IG fractions in non-erythroid, non-B clusters; SoupX reduced mode only if load exceeds a pre-set level | measure only | measure; SoupX-reduced sensitivity arm | measure only | measure; SoupX-reduced sensitivity arm | SoupX as in KPMP; CellBender not possible without raw droplets |
| Doublets | scDblFinder per capture, flag only | per library; freemuxlet already applied | per sample | per sample | per sample | per library; drop clusters mixing two lineages |
| Stress / cycle | score FOS/JUN/HSP module and S/G2M; never regress | score | score | score | score | score; keep injury states |
| Protected populations | before removal, test whether flags are enriched in MZB1+/JCHAIN+, LDG (FCGR3B+) and PT cells; keep any coherent marker-defined flagged population | yes | yes | yes (LDG file) | yes | PT cells |
| Report | cells removed per sample by reason, cross-tabulated by timepoint or disease state | | | | | |

**Stage 2: representations (all features computed per donor-timepoint)**

| ID | Representation | Notes |
|---|---|---|
| B0 | QC-only: median log counts, log genes, mito, ribo, ambient score, doublet score, stress score, cells recovered | Technical floor; reproduces the round 98 analysis at donor level |
| B1 | 6-gene IFN score, donor mean in monocytes | Project's biological baseline |
| B2 | CLR cell-type composition (ECODA-style) | Most batch-robust baseline per Halter 2026 |
| B3 | Pseudobulk log-CPM per major cell type, PCA fitted inside each fold | Most batch-exposed baseline |
| B4 | PCA-50 on log-normalised HVGs, donor mean per cell type | Round 98 comparator |
| D1 | scVI latent (n_latent 20, n_layers 2, layer norm, NB, early stopping, no continuous covariates), donor mean per cell type | Trained inside each fold on training donors; held-out donors mapped by scArches |
| D2 | sysVI latent | Only for cross-chemistry pooling (Perez with trials) or KPMP sc+sn |
| D3 | Geneformer zero-shot latent | Negative control, expected to lose |
| E1 | MrVI sample distances | Exploratory only; not scored against H1 |

**Stage 3: contrasts, batch keys and validation**

| Contrast | Data | Unit and n | scVI batch key | Validation |
|---|---|---|---|---|
| C1: pre vs month 3 | A-319 | 11 donor pairs | donor | Leave-one-donor-out; within-donor paired score difference; exact sign-flip null (2,048 flips) |
| C2: pre vs post | CAR-T | 7 pairs | patient | Same as C1 (128 flips). Interpretation is capped by lymphodepletion and reconstitution timing |
| C3: flare vs treated | Perez cohort 3 | 10 pairs | `library_uuid` | Paired sign-flip |
| C4: SLE vs control | Perez cohort 4 | 52 vs 44 donors | `library_uuid` | 5-fold donor-grouped CV × 20 repeats; label-permutation null (1,000) |
| C5: external transfer | Train on Perez cohort 3; test on GSE254176 flare vs before-flare/remission | 6 patients | dataset as system (scVI vs sysVI compared on the 57 Perez bridging donors) | Frozen model, single pass |
| K1: kidney mapping | KPMP reference; query is the E-MTAB-13596 LN biopsy | n=1 | `LibraryID`; sc vs sn as a system | Descriptive only: scANVI labels plus kNN uncertainty; no hypothesis test |

**Pre-registered decision rules and controls.** The primary metric is donor-level AUC, or paired accuracy for the paired contrasts, with a donor-bootstrap 95% CI. H1 holds for a contrast only if all four conditions below are met. If any one fails, the result is reported as "no demonstrated gain", which is the outcome the literature predicts.

1. D1 beats the best of B1–B4 by at least **0.05 AUC**, and the paired-bootstrap CI excludes 0.
2. D1 beats B0 by at least the same margin, and D1+B0 beats B0 alone.
3. The gain survives re-running on QC-matched cells, with nearest-neighbour matching on B0 metrics within cell type. The mandatory check is that B0 drops to about 0.5 after matching.
4. A batch probe trained on the embedding does not predict library or cohort better than it predicts the label.

The rules also cover leakage, robustness and reporting. Every unsupervised model is fitted inside the fold, and any transductive fit is labelled as such. Results are reported under the default and a strict QC setting (3 MAD on all metrics), because composition shifts with filtering. Linear QC residualisation is reported alongside matching. Within-cell-type DE runs as pseudobulk (edgeR or limma-voom) on raw counts, with donor pairing in the design. One analytic seed and one hyperparameter set are frozen before results are seen. Each round's prereg doc records the hash of the QC and model config.
