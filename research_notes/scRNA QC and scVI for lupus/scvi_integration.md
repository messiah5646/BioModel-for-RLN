# scVI-family models for single-cell integration and batch correction (2022-2026)

Scope note: 15 tool calls were used. Pages actually fetched this session are cited plainly. Standard references that were not fetched this session are marked "(not fetched this session)". Treat those as pointers, not verified quotes.

## Q1. scVI/scANVI recommended settings in scvi-tools 1.x

### Takeaway
The scVI defaults are n_latent=10, n_layers=1, n_hidden=128, zinb likelihood and gene-level dispersion. Treat them as a starting point. For reference mapping with scArches the docs recommend different settings: n_layers=2, dropout 0.2, layer norm on, batch norm off, encode_covariates=True. The docs say covariate keys are for nuisance factors only. A scverse developer advises leaving continuous covariates out unless a technical gradient is clearly visible.

### Cited Findings
- `scvi.model.SCVI` constructor defaults:
  - n_hidden=128, n_latent=10, n_layers=1, dropout_rate=0.1.
  - dispersion='gene'. Options are 'gene-batch', 'gene-label' and 'gene-cell'.
  - gene_likelihood='zinb'. Options are 'nb', 'poisson' and 'normal' (marked experimental).
  - latent_distribution='normal', or 'ln' (logistic normal).
  - Source: [scvi-tools 1.3.3 SCVI API](https://docs.scvi-tools.org/en/1.3.3/api/reference/scvi.model.SCVI.html)
- `setup_anndata` arguments:
  - batch_key: the obs column holding batch. If None, all cells get one batch.
  - labels_key: the obs column holding labels.
  - size_factor_key: replaces library size as the likelihood offset, on a linear scale.
  - categorical_covariate_keys / continuous_covariate_keys: "added alongside batch and treated as nuisance factors, which the model tries to remove from the latent space". The docs say they "should not be used for biologically-relevant factors".
  - Source: [SCVI API](https://docs.scvi-tools.org/en/1.3.3/api/reference/scvi.model.SCVI.html)
- `train()` options:
  - max_epochs=None uses `get_max_epochs_heuristic()`, which scales with the number of cells.
  - early_stopping defaults to False and must be turned on explicitly.
  - Source: [SCVI API](https://docs.scvi-tools.org/en/1.3.3/api/reference/scvi.model.SCVI.html)
- `get_normalized_expression(transform_batch=...)`: None uses the observed batch. An int or a batch name decodes every cell as if it came from that batch, which is a counterfactual decode. Source: [SCVI API](https://docs.scvi-tools.org/en/1.3.3/api/reference/scvi.model.SCVI.html)
- `differential_expression` defaults:
  - mode='vanilla', with 'change' as the alternative.
  - delta=0.25 (LFC threshold in change mode).
  - fdr_target=0.05.
  - batch_correction=False by default.
  - Source: [SCVI API](https://docs.scvi-tools.org/en/1.3.3/api/reference/scvi.model.SCVI.html)
- Settings the scvi-tools team found to work with scArches. For the reference SCVI model:
  - use_layer_norm="both", use_batch_norm="none"
  - encode_covariates=True
  - dropout_rate=0.2, n_layers=2
  - HVGs: 2000 batch-aware genes for SCVI, chosen on the reference only. totalVI uses 4000 seurat_v3 genes.
  - Source: [scvi-tools reference mapping tutorial](https://docs.scvi-tools.org/en/1.0.3/tutorials/notebooks/scrna/scarches_scvi_tools.html)
- Developer advice on covariates (scverse forum, user cane11):
  - Correlated or nested categorical covariates are fine.
  - "adding continuous covariates can make the latent space significantly worse".
  - Train without them first. Add them only if a strong technical gradient remains.
  - Check the result against known cell types or trajectories, or with scib-metrics.
  - Source: [scverse discourse: nested batch effects with scVI](https://discourse.scverse.org/t/nested-batch-effects-with-scvi/1910)
- scANVI should be initialised from a pre-trained scVI model with `SCANVI.from_scvi_model(..., unlabeled_category, labels_key)`. The tutorial trains about 20 more epochs with n_samples_per_label=100. Source: [reference mapping tutorial](https://docs.scvi-tools.org/en/1.0.3/tutorials/notebooks/scrna/scarches_scvi_tools.html)
- Primary method papers, not fetched this session:
  - scVI: Lopez et al. 2018, Nat Methods 15:1053, [doi:10.1038/s41592-018-0229-2](https://doi.org/10.1038/s41592-018-0229-2)
  - scANVI: Xu et al. 2021, Mol Syst Biol 17:e9620, [doi:10.15252/msb.20209620](https://doi.org/10.15252/msb.20209620)
  - scvi-tools: Gayoso et al. 2022, Nat Biotechnol 40:163, [doi:10.1038/s41587-021-01206-w](https://doi.org/10.1038/s41587-021-01206-w)

### Inferences
- A practical default for PBMC or kidney atlases:
  - Raw counts in a counts layer and 2000-4000 batch-aware HVGs.
  - batch_key = the strongest technical unit (sample, 10x lane or pool). Put other technical factors in categorical_covariate_keys.
  - n_latent 10-30 and n_layers 1-2. Use n_layers=2 with layer norm if the model will serve as a scArches reference.
  - gene_likelihood: 'nb' is usually enough for UMI data, and 'zinb' is the default. I found no head-to-head numbers; see Gaps.
  - early_stopping=True.
- Continuous covariates such as percent_mito and n_counts get the same treatment as batch, so their effect is removed from the latent space. If a QC metric correlates with biology, this risks removing biology (Q4).

### Gaps
- I did not reach the current 1.4/1.5 docs on max_epochs (heuristic of about min(round(20000/n_cells*400), 400)) or on nb vs zinb guidance. Both are from memory and unverified. I found no quantitative benchmark of nb vs zinb.
- I found no official guidance on choosing n_latent beyond the defaults.

## Q2. sysVI for strong batch effects (cross-system)

### Takeaway
sysVI is a cVAE with a VampPrior and a cycle-consistency loss, published by Hrovatin et al. in BMC Genomics in 2025 (bioRxiv 2023). It is scvi-tools' recommended model when batch effects are substantial: cross-species, organoid vs tissue, and sc vs sn. It gives clearly stronger batch correction than scVI while keeping comparable biological preservation. For weak, within-lab batch effects scVI preserves biology slightly better.

### Cited Findings
- Citation: Hrovatin K, Moinfar AA, Zappia L, et al., Theis FJ. "Integrating single-cell RNA-seq datasets with substantial batch effects." BMC Genomics 26:974, published 30 Oct 2025, [doi:10.1186/s12864-025-12126-3](https://doi.org/10.1186/s12864-025-12126-3). The bioRxiv preprint is [10.1101/2023.11.03.565463](https://doi.org/10.1101/2023.11.03.565463). Sources: [MIT DSpace record](https://dspace.mit.edu/handle/1721.1/163744); [scvi-tools sysVI guide](https://docs.scvi-tools.org/en/stable/user_guide/models/sysvi.html)
- Mechanism:
  - Cycle-consistency: embed a cell, decode it as another system, re-embed it, then minimise the distance between the two embeddings. Only same-cell comparisons are made, which helps preserve biology.
  - VampPrior: a multimodal learned prior that is more robust to the number of modes than a GMM prior.
  - Class: `scvi.external.SysVI`.
  - Input should be roughly normal, i.e. size-factor normalised and log1p. Raw counts are not used.
  - Source: [scvi-tools sysVI guide](https://docs.scvi-tools.org/en/stable/user_guide/models/sysvi.html)
- When to use which model, per the docs:
  - sysVI for substantial effects (cross-species, organoid-tissue), "where other models often fail".
  - scVI for weak effects such as multiple subjects from one lab, where scVI has "slightly higher biological preservation".
  - The docs advise comparing several runs over cycle weights and seeds.
  - Source: [sysVI guide](https://docs.scvi-tools.org/en/stable/user_guide/models/sysvi.html)
- Benchmark datasets:
  - Mouse vs human islets: 263k mouse cells and 192k human cells.
  - Organoid vs tissue retina: 43.5k and 54.5k cells.
  - Adipose sc vs sn: 28k cells and 58k nuclei.
  - Retina atlas sc vs sn: 266k cells and 1.78M nuclei.
  - Mouse vs human skin.
  - Source: [Hrovatin 2025 BMC Genomics](https://link.springer.com/article/10.1186/s12864-025-12126-3)
- Metrics: iLISI for batch mixing, NMI for cell types, and Moran's I for within-cell-type variation. Source: [Hrovatin 2025](https://link.springer.com/article/10.1186/s12864-025-12126-3)
- Results:
  - Raising the KL weight "does not improve integration". It collapses latent dimensions, so the apparent batch gain comes from fewer effective dimensions.
  - Adversarial methods (ADV, GLUE) mix unrelated cell types, such as acinar with immune cells and astrocytes with Müller glia. GLUE had the highest iLISI but significantly lower Moran's I.
  - scVI kept biology well but corrected batch weakly: beta and ductal cells stayed split by system.
  - Harmony-py kept within-cell-type structure well but corrected batch weakly.
  - Seurat had low within-cell-type preservation and took more than 20 hours on large data.
  - VAMP+CYC (sysVI) had significantly higher batch correction than the cVAE baseline on most datasets (adj. p<0.1), with comparable biological preservation.
  - Runtime: about 500k cells in under 50 minutes, about twice scVI's time. It scaled to the 2M-cell retina atlas.
  - Source: [Hrovatin 2025](https://link.springer.com/article/10.1186/s12864-025-12126-3)
- Hyperparameters:
  - Cycle weight of about 5 in the tutorial. The pancreas runs selected 2 or 5.
  - An automatic choice of 50 on the unbalanced skin data was cautioned against.
  - About 5 VampPrior components is enough.
  - Source: [Hrovatin 2025](https://link.springer.com/article/10.1186/s12864-025-12126-3)

### Inferences
- In a lupus or LN project, the "system" split that fits sysVI is scRNA-seq PBMC or kidney vs snRNA-seq kidney, e.g. KPMP mixes sc and sn data. Within-cohort sample or lane batches are weak effects, and scVI is the better choice for those.
- A sysVI model trains on normalised log data, so its decoded output is not a count model. Use it for the embedding, not for DE.

### Gaps
- I did not extract exact iLISI, NMI or Moran's I values; they are in the figures and supplement. I found no independent third-party benchmark of sysVI.

## Q3. Benchmarks: scIB 2022, scIB-E and newer; where Harmony matches scVI; over-correction

### Takeaway
In scIB 2022, scANVI, Scanorama, scVI and BBKNN did best on complex atlas-scale tasks. Harmony and Seurat did well on simple tasks and simulations, but Harmony usually ranked outside the top third on real data. Newer work (scIB-E 2024, sysVI 2025) shows label-level metrics miss within-cell-type biology. Strong correction, whether adversarial, high KL or high supervision, erodes that biology. That is the over-correction risk.

### Cited Findings
- scIB design:
  - 10 methods and 38 method-preprocessing combinations.
  - 9 tasks (2 simulations, 5 RNA, 2 ATAC) totalling 1.2M cells, 77 batches and 23 publications.
  - 14 metrics. The overall score weights batch removal 40% and bio-conservation 60%.
  - Source: [Luecken et al. scIB preprint full text](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2.full-text); published as Nat Methods 19:41 (2022), [doi:10.1038/s41592-021-01336-8](https://link.springer.com/article/10.1038/s41592-021-01336-8)
- scIB rankings:
  - "BBKNN, Scanorama, and scVI perform well, particularly on complex integration tasks."
  - Seurat v3 and Harmony do well on simulations.
  - Harmony usually ranks outside the top third on real data but is very usable, and is recommended for limited hardware or large data.
  - Seurat v3 failed above 100k cells.
  - LIGER performed poorly.
  - Source: [scIB full text](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2.full-text)
  - The published version ranks scANVI (with labels), Scanorama, scVI and scGen top on complex tasks (from memory; that version was not fetched this session).
- HVG selection and scaling in scIB:
  - HVG selection raised overall scores in 72% of comparisons.
  - Full gene sets were better for trajectory and cell-cycle conservation.
  - Scaling helped batch removal in 63% of comparisons but hurt bio-conservation in 72%.
  - Source: [scIB full text](https://www.biorxiv.org/content/10.1101/2020.05.22.111161v2.full-text)
- scIB-E (bioRxiv Dec 2024, "Benchmarking deep learning methods for biologically conserved single-cell integration"):
  - Splits bio-conservation into inter-cell-type and intra-cell-type parts, adding PCR-comparison-batch and Jaccard metrics.
  - Tests 16 VAE loss variants at three supervision levels, with scVI and scANVI as baselines.
  - Adding a Corr-MSE loss raised intra-cell-type conservation by 6.61±5.49% and cut batch correction by 4.02±1.84%.
  - High CellSupcon weight raised inter-cell-type conservation but eroded intra-cell-type information, which is over-correction.
  - On the fetal lung atlas, plain scVI scored highest on intra-cell-type conservation.
  - Source: [bioRxiv 2024.12.09.627450](https://www.biorxiv.org/content/10.1101/2024.12.09.627450.full.pdf)
- Over-correction mechanisms in sysVI 2025:
  - Adversarial alignment merges unrelated cell types.
  - KL up-weighting removes information by collapsing dimensions.
  - Harmony keeps within-type structure but under-corrects strong effects.
  - Source: [Hrovatin 2025](https://link.springer.com/article/10.1186/s12864-025-12126-3)
- A 2025 Genome Biology integration benchmark exists ([s13059-025-03869-z](https://link.springer.com/article/10.1186/s13059-025-03869-z)), but the fetch failed and its content is unverified.

### Inferences
- Harmony matches or beats scVI when batch effects are moderate and within-cell-type structure matters more than full mixing. It is also fast. scVI and scANVI pull ahead at atlas scale with many complex, nested batches.
- Disease confounded with batch, e.g. all SLE samples on one chemistry or site: any method given batch_key will partly remove the disease signal. Supervision (scANVI) protects cell-type labels, not disease states inside a cell type.
  - Design mitigation: batch_key = processing batch, not donor, when disease is donor-level.
  - Check that a known disease programme, such as the IFN score in SLE monocytes, survives integration.
- scIB metrics reward cell-type separation. A high scIB score does not show that disease-relevant within-type variation survived. Add Moran's I or scIB-E-style intra-type metrics.

### Gaps
- I did not verify the exact published scIB rank order, including scANVI's position, or the 2025 Genome Biology benchmark's findings.
- I found no benchmark specifically on SLE or LN data.

## Q4. Regressing QC metrics with continuous_covariate_keys

### Takeaway
This is supported mechanically: the covariate is fed to the decoder and its effect is removed from the latent space. The scvi-tools developers advise against using it by default because it can degrade the latent space. It will remove any biology correlated with the QC metric.

### Cited Findings
- Continuous covariates are nuisance factors the model removes from the latent space. Do not use them for biologically relevant factors. Source: [SCVI API](https://docs.scvi-tools.org/en/1.3.3/api/reference/scvi.model.SCVI.html)
- "adding continuous covariates can make the latent space significantly worse." The advice is to add one only if a strong technical gradient remains after a model without it. Source: [scverse discourse](https://discourse.scverse.org/t/nested-batch-effects-with-scvi/1910)
- library size is already modelled as an offset, either observed or via size_factor_key. Source: [SCVI API](https://docs.scvi-tools.org/en/1.3.3/api/reference/scvi.model.SCVI.html)

### Inferences
- n_counts is already handled by the library-size offset, so adding it as a covariate is redundant and can over-correct.
- percent_mito tracks real biology: metabolic state, tubular epithelial cells in kidney (PT cells are mito-rich) and stressed or apoptotic cells. Regressing it in kidney data risks removing PT injury states.
- Preferred order:
  1. Filter with MAD-based QC.
  2. Train without the covariate.
  3. Check whether latent dimensions correlate with percent_mito within cell types.
  4. Only then compare a run with the covariate, using scIB plus a within-type biology check.

### Gaps
- I found no published quantitative study of scVI continuous covariates removing biology. The evidence is developer advice plus the mechanism.

## Q5. Reference mapping (scArches, scANVI label transfer) onto atlases

### Takeaway
scArches-style mapping freezes the reference weights and trains only the new-batch parameters on the query. The reference needs scArches-friendly settings, and the query genes must be aligned to the reference genes. scANVI then transfers labels via `predict()`. The tutorial reported about 97% accuracy on its example.

### Cited Findings
- Workflow:
  1. `SCVI.prepare_query_anndata(query, ref_model)` reorders the genes and zero-pads missing ones.
  2. `load_query_data` builds the query model.
  3. Train the query model: SCVI max_epochs=200 or SCANVI 100, with plan_kwargs weight_decay=0.0. This keeps the reference embeddings unchanged.
  4. scANVI `predict()` gives the label transfer, about 97% accuracy in the tutorial.
  - Source: [scvi-tools reference mapping tutorial](https://docs.scvi-tools.org/en/1.0.3/tutorials/notebooks/scrna/scarches_scvi_tools.html)
- Human Lung Cell Atlas mapping and label-transfer tutorials exist: an scArches tutorial and an scvi-tools kNN-uncertainty tutorial. Sources: [scArches HLCA tutorial](https://docs.scarches.org/en/stable/hlca_map_classify.html); [scvi-tools query_hlca_knn](https://docs.scvi-tools.org/en/0.20.1/tutorials/notebooks/query_hlca_knn.html)
- scArches paper (not fetched this session): Lotfollahi et al. 2022, Nat Biotechnol 40:121, [doi:10.1038/s41587-021-01001-7](https://doi.org/10.1038/s41587-021-01001-7)
- Perez et al. 2022 lupus PBMC atlas (not fetched this session): Science 376:eabf1970, [doi:10.1126/science.abf1970](https://doi.org/10.1126/science.abf1970). It is a multiplexed scRNA-seq dataset of about 1.2M PBMCs from SLE patients and controls.
  - Recalled from memory and not verified: 162 SLE and 99 control donors, and the data may be on CELLxGENE.
  - I found no published pretrained scVI or scANVI model for it.

### Inferences
- To map onto Perez 2022, or onto KPMP or HCA kidney data, you usually need to train your own scArches-compatible reference on the atlas counts. Most atlases ship only embeddings, not scvi-tools models.
- Check query gene overlap with the reference HVGs. Heavy zero-padding degrades mapping.
- Use kNN or scANVI uncertainty to flag unseen states, such as LN-specific infiltrating cells, rather than forcing labels.
- Disease states absent from a healthy-heavy reference may be mapped onto the nearest healthy type. scArches does not create new axes, so check where the residual variation goes.

### Gaps
- I did not verify the existence or location of pretrained scvi-tools models for KPMP or HCA kidney, or for Perez 2022.
- I did not fetch scArches' quantitative benchmarks.

## Q6. DE and counterfactuals in scVI latent space; corrected values vs pseudobulk

### Takeaway
scVI offers Bayesian DE and counterfactual decoding via transform_batch. Cell-level tests treat cells as independent and inflate false discoveries. Pseudobulk across biological replicates is the validated standard, so DE should not be run on batch-corrected or denoised values as if they were data.

### Cited Findings
- Squair et al. 2021, Nat Commun 12:5692:
  - Used 18 ground-truth datasets pairing bulk and sc data from the same purified cells.
  - Of 14 DE methods, the 6 best were all pseudobulk.
  - Single-cell methods are biased toward highly expressed genes and call spike-ins DE.
  - They gave hundreds of false positives in null replicate splits.
  - RNAscope validated 5/6 pseudobulk hits vs 3/13 Wilcoxon hits.
  - Source: [Squair 2021 Nat Commun](https://www.nature.com/articles/s41467-021-25960-2)
- scVI `differential_expression` compares posterior expression between groups (vanilla or change mode, delta=0.25, fdr_target=0.05). Batch correction is optional (batch_correction=False by default). Source: [SCVI API](https://docs.scvi-tools.org/en/1.3.3/api/reference/scvi.model.SCVI.html)
- `transform_batch` decodes cells as if they came from a specified batch, which is a counterfactual expression decode. Source: [SCVI API](https://docs.scvi-tools.org/en/1.3.3/api/reference/scvi.model.SCVI.html)

### Inferences
- Use the scVI latent space for clustering, annotation and neighbourhood analysis.
  - Run DE as pseudobulk on raw counts per donor × cell type (edgeR, DESeq2 or limma-voom), with batch as a design covariate.
  - Corrected or denoised values are smoothed and have false certainty, and cells within a donor are not independent replicates.
- Counterfactual decodes (transform_batch, or decoding with a changed covariate) only change the nuisance axes the model was conditioned on. They are not causal treatment simulations. An scVI counterfactual of "SLE vs healthy" is not possible unless disease was a covariate, and that would remove it from the latent space.
- For perturbation counterfactuals, scvi-tools has other models such as contrastiveVI and CPA-like ones. That is outside this scope and was not researched.

### Gaps
- I found no benchmark comparing scVI's built-in Bayesian DE against pseudobulk on ground truth. The Squair 2021 benchmark did not, to my knowledge, include scVI DE (unverified).
