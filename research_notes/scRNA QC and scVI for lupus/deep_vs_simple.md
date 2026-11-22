# Deep generative / foundation-model latent spaces vs simple baselines for single-cell tasks, and technical confounding of embeddings (evidence 2023-2026)

Scope note: written for the decision "is an scVI (or foundation-model) latent space worth using over PCA on log-normalised HVGs for sample-level classification (pre vs post therapy, flare vs remission, LN vs non-LN donors)?". Sources were fetched in this session unless marked *[not re-fetched; cited from prior knowledge, verify]*. Many 2026 items are preprints and have not been peer reviewed.

## 1. Benchmarks where simple baselines match or beat deep models, and where scVI does win

### Takeaway
Across independent benchmarks from 2023 to 2026, zero-shot foundation-model embeddings (Geneformer, scGPT, UCE, scFoundation) generally fail to beat HVG+PCA, scVI or Harmony at clustering and batch integration. Deep perturbation predictors fail to beat mean, additive or linear baselines. scVI is a different case. It is a strong, dedicated integration method that reliably beats foundation models and usually beats raw PCA at batch mixing. But no fetched source shows scVI beating PCA or pseudobulk at *sample-level* classification.

### Cited Findings
**Zero-shot foundation models vs HVG / scVI / Harmony (cell-level)**
- Kedzierska, Crawford, Amini & Lu (Genome Biology 2025, 26:101; preprint 2023) evaluated Geneformer and scGPT zero-shot on Pancreas, PBMC 12k, PBMC 95k, Immune and Tabula Sapiens. For cell-type clustering, HVG, Harmony and scVI outperformed both foundation models, and "HVG beats both foundation models on all metrics". scGPT beat the baselines only on PBMC 12k. — [Kedzierska et al. 2025, Genome Biol (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12007350)
- Batch integration in the same paper:
  - Geneformer ranked last on batch-mixing metrics on most datasets.
  - HVG selection achieved the best batch-integration scores on all datasets.
  - scVI and Harmony beat scGPT on Pancreas and PBMC, which have only technical batches.
  - scGPT beat scVI and Harmony on Tabula Sapiens and Immune, which have technical and donor batches. Those datasets were in scGPT's pretraining corpus, a confound the authors flag.
  - Source: [Kedzierska et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12007350)
- Kedzierska et al. also report that Geneformer's embedding "loses cell type information" and that its clustering is "primarily driven by batch effects". The share of variance explained by batch was *higher* in the Geneformer embedding than in the input data. — [Kedzierska et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12007350)
- A 2026 PLoS Comput Biol study of Geneformer at several scales found that ranking-based foundation models "remain less competitive than dedicated integration methods such as HVG-based representations and scVI" in zero-shot batch integration. — [Chen et al. 2026, PLoS Comput Biol (via Research Monitor)](https://research-monitor.com/papers/42531328)
- A bioRxiv preprint (Dec 2025, "Batch Effects Remain a Fundamental Barrier to Universal Embeddings in Single-Cell Foundation Models"; authors not shown on the fetched page) tested 8 foundation models zero-shot: Geneformer, scGPT, scFoundation, GeneCompass, CellPLM, UCE, Nicheformer and scCello.
  - Every model except scCello underperformed Harmony and scVI on scIB, and some fell below PCA.
  - The gap widened as batch heterogeneity grew (assay, RNA source).
  - Logistic-regression probes recovered *batch* more accurately from the foundation-model embeddings than from PCA or traditional methods.
  - Subtracting per-batch means reduced the batch signal only partly.
  - Source: [bioRxiv 2025.12.19.695371](https://www.biorxiv.org/content/10.64898/2025.12.19.695371v1.full)
- A 2026 bioRxiv preprint on spatial data reports that frozen foundation-model embeddings "offer little advantage over HVG + PCA on raw counts". — [McConnell et al. 2026, bioRxiv (via Research Monitor)](https://research-monitor.com/papers/id:213464)
- Souza & Mehta (arXiv 2602.16696, 2026) compared scTOP, a parameter-free rank-normalise / z-score / projection method, with foundation models.
  - Cross-species annotation: scTOP 0.679 macro-F1 vs 0.639 for Geneformer, the best foundation model.
  - Tabula Sapiens 2.0: a scTOP+PCA+logistic-regression pipeline reached 0.887 vs 0.911 for TF-Sapiens. Here the foundation model was slightly better.
  - SARS-CoV-2-infected vs uninfected cells: 0.850 macro-F1, "exceeding most foundation models".
  - Source: [Souza & Mehta 2026, arXiv](https://arxiv.org/pdf/2602.16696)

**Fine-tuned foundation models vs logistic regression**
- Boiarsky et al. (bioRxiv 2023, "A Deep Dive into Single-Cell RNA Sequencing Foundation Models") ran logistic regression on log library-size-normalised expression against fine-tuned scBERT and scGPT for cell-type annotation.
  - LR beat scBERT on PBMC Zheng68K in accuracy and macro-F1. Removing scBERT's pretraining did not change PBMC performance significantly (p = 0.34 for accuracy).
  - LR beat scGPT on myeloid data and tied on pancreas.
  - scGPT beat LR only on multiple sclerosis, a dataset-shift setting, and its advantage grew as training data shrank.
  - The authors conclude that pretraining may help mainly under distribution shift.
  - Source: [Boiarsky et al. 2023, bioRxiv](https://www.biorxiv.org/content/10.1101/2023.10.19.563100.full.pdf)

**Perturbation prediction**
- Ahlmann-Eltze, Huber & Anders (Nature Methods 2025, 22:1657–1661) compared scGPT, scFoundation, GEARS, CPA, and repurposed scBERT, Geneformer and UCE against simple baselines. "None outperformed the baselines."
  - Double perturbations (Norman data): every model's error was "substantially higher" than the additive baseline.
  - Unseen single perturbations (Adamson, Replogle K562, Replogle RPE1): no deep model consistently beat the mean model or a linear model on PCA embeddings (K = 10, λ = 0.1).
  - A linear model using scGPT or scFoundation gene embeddings matched or beat those models' own decoders. Pretraining on atlas data helped only slightly over random embeddings.
  - Source: [Ahlmann-Eltze et al. 2025 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12328236/)
- Csendes, Sanz, Szalay & Szalai (BMC Genomics 2025, 26:393; bioRxiv 2024) scored models by Pearson delta.
  - The "train mean" baseline beat scGPT and scFoundation on all 4 datasets: Adamson 0.711 vs 0.641 / 0.552; Norman 0.557 vs 0.554 / 0.459; Replogle K562 0.373 vs 0.327 / 0.269; Replogle RPE1 0.628 vs 0.596 / 0.471.
  - A random forest on GO features beat everything (0.739, 0.586, 0.480, 0.648).
  - Source: [Csendes et al. 2025 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12016270)
- Research Monitor's foundation-model field map also records that independent zero-shot evaluations found Geneformer and scGPT "do not consistently outperform simpler approaches in clustering and batch integration". — [Liu et al. 2026, Biosensors (via Research Monitor)](https://research-monitor.com/papers/42930128)

**Where scVI does win**
- scIB benchmark (Luecken et al., Nat Methods 2022): scANVI, Scanorama, scVI and scGen performed best on complex atlas-level integration tasks with nested batches. Harmony was competitive on simpler tasks. — [Luecken et al. 2022](https://www.nature.com/articles/s41592-021-01336-8) *[not re-fetched; cited from prior knowledge, verify]*
- Integration across heterogeneous technical batches is the setting where scVI's advantage over PCA is best documented. In the 2025 batch-effect preprint, scVI and Harmony beat most foundation models and PCA on scIB. — [bioRxiv 2025.12.19.695371](https://www.biorxiv.org/content/10.64898/2025.12.19.695371v1.full)
- Disagreement: Kedzierska et al. found that plain HVG selection scored *best* on batch-integration metrics on all 5 of their datasets, ahead of scVI. scVI's integration edge therefore depends on the dataset and the metric. — [Kedzierska et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12007350)

### Inferences
- The foundation-model case is weak. For an SLE/LN cohort, zero-shot Geneformer, scGPT or UCE embeddings should be treated as a baseline to beat, not a default. This matches the project's own rounds 92 and 98 (below).
- scVI is a separate question from foundation models. It earns its place for *integration* (joint clustering and annotation across processing batches). That is not evidence it helps *sample-level classification*, where downstream aggregation and donor-level n dominate.
- The repeated pattern across tasks is that a simple predictor (mean, additive, linear, LR, HVG+PCA) captures most of the signal that benchmarks measure. Benchmark wins for deep models are usually small and tied to particular datasets.

### Gaps
- I found no head-to-head benchmark of scVI latent vs PCA latent with the *same* downstream aggregation (e.g. per-donor mean of latents, then a classifier) for patient-level labels with donor-level CV.
- Exact AvgBIO and ASW numbers for Kedzierska et al. are in figures and supplements that were not retrieved.
- Small-N denoising: I did not retrieve a 2023-2026 source quantifying when scVI denoising beats log-normalisation for low-depth or small datasets.

## 2. Sample- and patient-level prediction from single-cell data

### Takeaway
For patient-level labels, the 2026 evidence favours simple representations: CLR-transformed cell-type proportions and pseudobulk. These match or beat MrVI, scPoli, MOFA+, GloScope and similar methods at a tiny fraction of the compute. The one large-scale win for a learned model (PaSCient) used thousands of patients and study-level splits, and reported that pseudobulk beat the foundation models. Cell-level cross-validation inflates accuracy, and donor-level (grouped) validation is mandatory.

### Cited Findings
**Composition and pseudobulk vs learned sample representations**
- Halter, Andreatta & Carmona (bioRxiv 2026, "Cell type composition drives patient stratification in single-cell RNA-seq cohorts") benchmarked 11 cohorts with 697 individuals across malignant, autoimmune and infectious diseases.
  - Baselines: ECODA (CLR-transformed cell-type proportions), pseudobulk, GloProp, and a shuffled-label control.
  - Learned methods: MOFA+, scITD, GloScope, MrVI, PILOT and scPoli.
  - ECODA ranked first on ANOSIM, graph modularity and ARI, followed by GloProp and pseudobulk, then GloScope, MOFA+, MrVI and scPoli. PILOT and scITD scored lowest.
  - The baselines ran in seconds. MrVI, scPoli and GloScope took hours on GPU and needed downsampling.
  - Source: [Halter et al. 2026, bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1.full)
- Batch result from the same paper, on the Joanito colorectal cohort: pseudobulk's batch ANOSIM (0.706) was *higher* than its tissue-type ANOSIM (0.322). ECODA's was 0.041 for batch vs 0.640 for tissue. Pseudobulk is the representation most vulnerable to batch, and composition the least. — [Halter et al. 2026](https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1.full)
- Caveats on Halter et al.:
  - The evaluation is unsupervised separation of known groups, not supervised AUC under donor-level CV.
  - No lupus cohort is included.
  - The authors do not claim learned representations are inherently inferior.
  - Source: [Halter et al. 2026](https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1.full)
- PaSCient (Liu, De Brouwer, … Regev, Heimberg, Genentech; bioRxiv 2024, Cell Systems 2026) trained on about 24.3M cells from over 5,000 samples, 135 disease labels and 413 studies on CELLxGENE.
  - The split was 60/20/20 *by study*, not by donor.
  - Baselines were pseudobulk and cell-type proportions with kNN or MLP, and CellPLM and SCimilarity with kNN or MLP. CloudPred and ProtoCell4P were cited but not benchmarked. ScRAT appears only in the appendix.
  - PaSCient beat all baselines "by a significant margin". The only numeric comparison in the text is severity accuracy: 0.88 vs 0.79 for pseudobulk.
  - The authors state that "a simple pseudo-bulk approach outperforms more complicated foundation models in this task."
  - Sources: [PaSCient bioRxiv v2](https://www.biorxiv.org/content/10.1101/2024.11.18.624166v2.full.pdf); [Cell Systems 2026](https://www.cell.com/cell-systems/fulltext/S2405-4712(26)00052-9)
- Muhire (arXiv 2605.03281, May 2026; single-author student work) built donor-aware benchmarks for IBD.
  - Datasets: SCP259, 30 donors (UC vs healthy); Kong 2023, 71 donors (CD vs healthy).
  - Compartment-stratified CLR composition with linear models reached AUROC 0.956 on SCP259.
  - A pseudobulk HVG XGBoost reached AUROC 1.000 at n = 30, which the author attributes to overfitting.
  - Compartment-level scVI latent features reached 0.931. That scVI model was trained on all donors before CV, so the estimate is effectively in-sample.
  - No between-method difference was statistically significant.
  - Source: [Muhire 2026, arXiv](https://arxiv.org/pdf/2605.03281)

**Leakage and pseudoreplication**
- Cang & Sun (bioRxiv, Sept 2026, "Cell-level random splits leak group-owned answers in single-cell benchmarks") show how leakage works.
  - Cells from the same patient, clone or batch share a label, so random cell-level splits let a model "memorise relatives".
  - They define two diagnostics: *exposure* (the fraction of test cells with relatives in training) and *retrievability* (whether nearest-neighbour search finds those relatives). Exposure decides whether leakage happens; retrievability decides how much it inflates scores.
  - In one patient cohort, the inflation was large enough to "overturn a clinical conclusion". Exact magnitudes were not on the fetched page.
  - Their tool is `leakcheck`.
  - Source: [Cang & Sun 2026](https://sciety.org/articles/activity/10.64898/2026.09.09.750484)
- Muhire 2026 asserts that cell-level splitting "can produce AUROC values approaching 1.0" but did not measure this head to head. — [Muhire 2026](https://arxiv.org/pdf/2605.03281)
- Pseudobulk DE methods outperform cell-level DE methods. Cell-level tests generate false discoveries because they treat cells as independent replicates. — [Squair et al. 2021, Nat Commun](https://www.nature.com/articles/s41467-021-25960-2); [Zimmerman et al. 2021, Nat Commun](https://www.nature.com/articles/s41467-021-21038-1) *[not re-fetched; cited from prior knowledge, verify]*

**Project-internal evidence (LUMEN project docs, not external literature)**
- Round 98 (2026-11-22) used 2,200 monocytes from 11 SLE donors before and 3 months after CD3×CD19 bispecific therapy (CELLxGENE dataset 3979a248), with leave-one-donor-out CV and cell-level AUC.
  - Raw AUCs: technical QC alone (depth, genes detected, % mito) 0.66; 6-gene IFN score 0.68; PCA-50 0.82; Geneformer 0.77.
  - On cells matched on QC (where QC alone fell to 0.49), PCA's lead over IFN shrank from 0.14 to 0.035, and Geneformer's to 0.03.
  - In donor P12, PCA and Geneformer predicted backwards. The IFN score did not.
  - Source: project doc `claude/round-98-the-IFN-score-survives-QC-adjustment-the-embedding-edge-mostly-does-not.md`

### Inferences
- In LN cohorts the number of donors is typically 10 to 100, so the effective n is donors (or donor-timepoints), not cells. With that n, a learned sample embedding (MrVI, scPoli, PaSCient-style) has no demonstrated advantage over CLR composition or pseudobulk with a regularised linear model. The large-n PaSCient result does not transfer.
- The project's round 98 still scores *cells* (cell-level AUC under LODO). Donor-level aggregation, e.g. per-donor mean score and one prediction per donor-timepoint, plus a permutation null over donor labels, would be the stricter test.
- Pseudobulk carries batch strongly (Halter: batch ANOSIM 0.71). In pre/post designs where each timepoint is its own library, pseudobulk and embeddings both risk classifying *library*. Composition is more robust to this, but is itself sensitive to QC filtering thresholds.

### Gaps
- No fetched benchmark reports supervised, donor-level-CV AUCs for CloudPred, ProtoCell4P or scPoli sample embeddings vs pseudobulk+LR on the same cohorts. CloudPred and ProtoCell4P papers were not retrieved.
- No lupus or LN cohort appears in Halter 2026, PaSCient or the IBD benchmark.
- Cang & Sun's quantitative inflation numbers were not retrieved.

## 3. Technical confounding of latent spaces (QC metrics, ambient RNA, batch)

### Takeaway
Foundation-model embeddings measurably retain batch information, often more than the input data does. Low-depth or high-dropout cells form artefactual structure. In paired pre/post designs where each timepoint is a separate library, nominal "biological" separation can be largely technical. The project's own matched analysis shows this. Practical controls: a QC-only baseline classifier, QC-matched cells, residualisation (which misses non-linear leakage), batch probes, donor-label permutation, and multiplexing designs.

### Cited Findings
- Geneformer embeddings had a higher share of variance explained by batch than the original data, and their clustering was "primarily driven by batch effects". On Pancreas, scGPT's main structure was also batch-driven. — [Kedzierska et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12007350)
- Linear probes predicted batch (assay, RNA source, donor, tissue) more accurately from foundation-model embeddings than from PCA, scVI or Harmony. Per-batch mean-centering only partly reduced this and could overcorrect.
  - The preprint did *not* test library size, % mito or processing date explicitly.
  - Recommended controls: report batch probes alongside label probes, benchmark against Harmony, scVI and PCA, use batch-centering as a post-hoc control, and monitor overcorrection.
  - Source: [bioRxiv 2025.12.19.695371](https://www.biorxiv.org/content/10.64898/2025.12.19.695371v1.full)
- Souza & Mehta used the number of detected genes as a dropout proxy. A high-dropout cluster of heart smooth-muscle cells was misclassified because of technical artefact. They argue that much of the curvature in Tabula Sapiens comes from technical noise and batch effects, which per-cell normalisation reduces. — [Souza & Mehta 2026](https://arxiv.org/pdf/2602.16696)
- In a slide-level spatial study, batch structure limited supervision of transcriptomic foundation models, and frozen embeddings added little over HVG+PCA. — [McConnell et al. 2026 (via Research Monitor)](https://research-monitor.com/papers/id:213464)
- Pseudobulk sample representations separated batch more strongly than biology in one cohort (batch ANOSIM 0.706 vs tissue 0.322), while CLR composition did not (0.041 vs 0.640). — [Halter et al. 2026](https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1.full)
- scVI models library size explicitly, as an observed size factor or a latent library variable, and takes batch as a conditioning covariate. It does *not* by default condition on % mito, ribosomal % or ambient profile. Those can be passed as continuous covariates, which also removes any biology correlated with them. — [scvi-tools scVI user guide](https://docs.scvi-tools.org/en/stable/user_guide/models/scvi.html) *[search result only; page content not fetched, verify]*
- Ambient RNA contaminates droplets with highly expressed genes from other cell types. SoupX and FastCAR were developed because it distorts differential expression between conditions. — [SoupX, Young & Behjati 2020 (PMC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7763177/); [FastCAR, BMC Genomics 2023](https://link.springer.com/10.1186/s12864-023-09822-3) *[search results only; quantitative content not fetched]*
- Project-internal evidence (round 98):
  - Three QC numbers alone gave AUC 0.66 for pre vs post in SLE monocytes, about the same as the IFN score (0.68).
  - Linear residualisation on the 3 QC numbers left PCA's edge intact.
  - QC-matching removed most of it.
  - So linear regression of QC covariates does not remove non-linear technical structure spread across 50 latent dimensions.
  - Source: project doc `claude/round-98-…md`

### Inferences
- Minimum control panel for any embedding-based sample classifier in this project:
  - a QC-only classifier (depth, genes, % mito, % ribo, ambient score, doublet score) under the same CV;
  - an "embedding + QC vs QC alone" increment;
  - QC-matched or QC-stratified re-evaluation;
  - a batch/date probe on the embedding;
  - a donor-label permutation null;
  - donor-level aggregation before scoring.
- Residualising % mito can remove real activation or stress biology. Matching is the cleaner diagnostic and residualisation the conservative one. Report both.
- Design fix: the definitive test is multiplexed (hashtag or genotype) libraries in which pre and post, or flare and remission, cells share a library. Only that breaks the library-timepoint confound. No retrospective method fully does.

### Gaps
- I found no 2023-2026 paper that quantifies specifically how much variance in scVI latent dimensions is explained by % mito, ribosomal % or ambient RNA in a case-control cohort.
- I did not retrieve a documented published case in which a disease or therapy classifier was later shown to separate on processing date. Cang & Sun describe a "conclusion-overturning" leakage case, but the details were not retrieved.

## 4. MrVI for sample-level heterogeneity in a lupus cohort

### Takeaway
MrVI (Boyeau et al., Nature Methods 2025, 22:2264–2274) is the most principled deep option for cohort designs.
- It separates a sample-invariant latent *u* from a sample-aware latent *z*, conditioning on one nuisance covariate (e.g. processing batch).
- It produces per-cell sample-by-sample distance matrices, plus counterfactual DE and DA without predefined clusters.
- Its published validation is mostly qualitative or figure-only.
- In a 2026 independent unsupervised benchmark it ranked below CLR composition and pseudobulk.
- The published version handles only one target and one nuisance covariate, with no paired or repeated-measures structure.

Use it as an exploratory tool for finding *which cell states* vary across patients. It is not a classifier that beats pseudobulk.

### Cited Findings
- What MrVI provides:
  - *u* is invariant to sample and batch. *z* adds sample-specific effects while correcting nuisance covariates such as sequencing run or processing day.
  - For each cell it computes counterfactual *z* under every sample, which yields a cell-specific sample-sample distance matrix. Clustering these matrices finds cell subsets with distinct patient stratifications.
  - DE: counterfactuals are regressed on sample covariates with χ² tests, giving per-cell log fold-changes.
  - DA: log-ratio of aggregated posteriors between groups.
  - Citation: Boyeau, Hong, Gayoso et al., Nat Methods 2025, 22(11):2264–2274, doi:10.1038/s41592-025-02808-x.
  - Source: [scvi-tools MrVI guide](https://docs.scvi-tools.org/en/stable/user_guide/models/mrvi.html)
- Validation in the preprint:
  - Semi-synthetic PBMC with 32 synthetic subjects: MrVI's DA was more accurate than Milo, and its DE correlated better with pseudobulk references than miloDE did.
  - COVID-19 PBMC, 419K cells: the *u* space mixed samples better than PCA and scVI, and the patient split was significant (Mann-Whitney p < 0.05). The stratification did not perfectly match infection timelines.
  - Crohn's, 46 patients and 25 controls: default MrVI mixed samples better than scVI but "slightly lagged in preserving cell-state distinctions".
  - sci-Plex: hyperparameters were tuned on the same metrics used for evaluation.
  - Source: [Boyeau et al. bioRxiv 2022.10.04.510898](https://www.biorxiv.org/content/10.1101/2022.10.04.510898.full.pdf)
- Stated scope limits: the version described handles one target and one nuisance covariate. Multiple covariates, continuous covariates and other modalities are future work. No minimum sample number is given. — [Boyeau et al. bioRxiv](https://www.biorxiv.org/content/10.1101/2022.10.04.510898.full.pdf); [scvi-tools MrVI guide](https://docs.scvi-tools.org/en/stable/user_guide/models/mrvi.html)
- Independent benchmark: in Halter et al. 2026 (11 cohorts, 697 individuals), MrVI ranked below ECODA, GloProp and pseudobulk for unsupervised patient stratification and needed hours of GPU time. — [Halter et al. 2026](https://www.biorxiv.org/content/10.64898/2026.03.27.714811v1.full)

### Inferences
- What MrVI offers a lupus cohort:
  1. Set target = sample (donor × timepoint) and nuisance = processing batch or 10x run.
  2. Cell-state-resolved sample distances can show *which* compartments, e.g. ISG-high monocytes or plasmablasts, separate flare from remission or LN from non-LN, without predefined clusters.
  3. Counterfactual DE gives per-cell fold-changes conditioned on batch.
- Hazards:
  - If processing batch is collinear with timepoint (pre and post in separate runs), declaring batch as nuisance will remove the therapy effect, and leaving it out lets the effect be confounded. MrVI cannot fix a confounded design.
  - Repeated measures per patient are not modelled explicitly. Donor pairing would have to be handled downstream.
  - With about 10 to 30 donors, the DSS-style tests have limited power. Validate any MrVI-derived stratification against CLR composition and pseudobulk baselines, with donor-level permutation.

### Gaps
- Quantitative MrVI benchmark numbers (scIB scores, DA AUCs) are in figures that were not retrieved.
- No application of MrVI to SLE or LN data was found.
- No published guidance was found on the minimum number of samples for MrVI's sample-distance tests.
