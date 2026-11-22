# scRNA-seq QC best practice (human 10x) for lupus PBMC and kidney biopsy data

Scope: cell-level QC thresholds, ambient RNA, doublets, stress artefacts and cell cycle, and the risk of over-filtering. Sources checked 2026-11-22. "Classic" marks pre-2023 work. Where a number comes only from a fetched-page summary and not from the primary PDF tables, I say so.

## 1. MAD-based adaptive thresholds: metrics, multipliers, mito caps, per-sample vs pooled, tissue-specific mito

### Takeaway
The scverse book's default: flag outliers at 5 MAD on log1p(total counts), log1p(n genes) and % counts in the top 20 genes, and at 3 MAD (upper side only) on mito fraction, combined with a hard mito floor of 8%. Run it per sample and keep it permissive. Lupus and kidney studies have used much looser fixed mito caps: 8% for 10x blood, 25% for kidney cells (Arazi 2019), and under 50% for KPMP kidney scRNA-seq. A global human 10% mito rule would remove proximal tubule cells and other naturally mito-high cells.

### Cited Findings
- scverse "Single-cell best practices" QC chapter: the outlier metrics are `log1p_total_counts` (5 MAD, both sides), `log1p_n_genes_by_counts` (5 MAD, both sides) and `pct_counts_in_top_20_genes` (5 MAD, both sides). The current page version scores mito as `logit_mt_fraction` with fraction = (mt+1)/(total+2) and flags at 3 MAD, upper side only. — [sc-best-practices QC](https://www.sc-best-practices.org/preprocessing-visualization/quality-control/)
- A cell counts as a mito outlier only if it is above 3 MAD **and** has `pct_counts_mt > 8`. The book warns that 8% "is not universal since mitochondrial content differs between cell types and tissues." — [sc-best-practices QC](https://www.sc-best-practices.org/preprocessing-visualization/quality-control/)
- The book calls 5 MAD "a relatively permissive filtering strategy." It advises being as permissive as possible and re-assessing filtering after cell annotation. Ribosomal (RPS/RPL) and haemoglobin (HB*) fractions are computed as diagnostics but do not drive filtering. In its example, filtering took 16,934 cells to 15,998 (about 5.5% removed). Gene filter used: `min_cells=20`. — [sc-best-practices QC](https://www.sc-best-practices.org/preprocessing-visualization/quality-control/)
- Heumos et al., Nat Rev Genet 2023: "Quality control is performed at the sample level as thresholds can vary substantially between samples." Metrics "are considered jointly to prevent the misinterpretation of cellular signals." Authors recommend "setting permissive thresholds initially and potentially removing more cells as necessary." — [Heumos 2023](https://www.nature.com/articles/s41576-023-00586-w.pdf) / [PMC10066026](https://pmc.ncbi.nlm.nih.gov/articles/PMC10066026)
- Osorio & Cai, Bioinformatics 2021 (classic): meta-analysis of 5,530,106 cells from 1,349 PanglaoDB datasets (91.3% 10x). They "suggest a standardized mtDNA% threshold of 10% for scRNA-seq QC of human samples" and 5% for mouse. In mouse, whole kidney and the proximal and distal tubule cell types are among the few with mean mito significantly above 5%. In human, neutrophils and macrophages are among the cell types above 5%. — [bioRxiv preprint text via SciSpace](https://scispace.com/pdf/systematic-determination-of-the-mitochondrial-proportion-in-rc82n0wqis.pdf); [bioRxiv 2020.03.02.973925](https://www.biorxiv.org/content/biorxiv/early/2020/03/04/2020.03.02.973925.full.pdf). Exact per-tissue means are in supplementary tables S1–S4, which I did not read.
- miQC (Hippen et al., PLoS Comput Biol 2021, classic) is a data-driven alternative that fits a mixture model of mito% against n genes instead of a fixed cap. On human menstrual blood 10x data, a fixed 10% cutoff removed 41% of cells and miQC removed 18.5%. — [miQC bioRxiv](https://www.biorxiv.org/content/10.1101/2021.03.03.433798v1.full-text)
- KPMP/HuBMAP kidney atlas (Lake et al., Nature 2023), 10x scRNA-seq (scCv3): mito "<50%", genes per cell >500 and <5,000 (the upper cap is to reduce multiplets), ambient correction with SoupX v1.5.0. snRNA-seq (snCv3): genes 400–7,500, MT-* genes removed, DoubletDetection v2.4.0. Altered or injured states (low genes, high mito/ER transcripts, CST3/HAVCR1/LCN2 up) were **flagged and kept**, not filtered. Clusters that mixed major classes were removed as ambiguous or low quality. — [Lake 2023 Nature](https://www.nature.com/articles/s41586-023-05769-3.pdf)
- AMP lupus nephritis Phase 1 (Arazi et al., Nat Immunol 2019, classic, CEL-Seq2 plate-based): kidney cells needed ≥1,000 genes and mito <25%, with <5,000 genes to remove putative multiplets. Urine cells used ≥500 genes. 10x cells used ≥250 genes, and 10x blood used mito <8%. Regressing on plate, UMI and mito was tested and "had a negligible effect," so it was not applied. Epithelial cells (CD10+CD45−, 10% of each sort) had lower quality, which the authors attributed to lower viability, and were kept as a single cluster. Two dividing subclusters were removed for high mito and stress genes. — [Arazi 2019 PDF](https://hacohenlab.mgh.harvard.edu/wp-content/uploads/2023/06/Arazi-et-al-Nat-Immun-2019.pdf)
- 10x neutrophil tutorial: example mito cutoff of 15% for blood. — [10x Cell Ranger neutrophil tutorial](https://www.10xgenomics.com/support/software/cell-ranger/tutorials/cr-tutorial-neutrophils)

### Inferences
- PBMC: use the scverse 5 MAD rules per sample, with a mito upper bound of 3 MAD capped at about 8–15%. Published lupus blood work used 8% (Arazi 10x blood) and 10x uses 15% in its example.
- Kidney biopsy, whole-cell 10x: a 10% mito cap is inappropriate. Use per-sample MAD on logit mito with a hard ceiling of roughly 25% (AMP) to 50% (KPMP). Better still, apply thresholds per broad compartment (epithelial vs immune) or re-check after annotation, because proximal tubule cells are mito-high by biology.
- For snRNA-seq kidney data, mito reads should be near zero. High mito in nuclei suggests cytoplasmic contamination; KPMP removes MT genes altogether.

### Gaps
- I did not obtain the Osorio & Cai per-tissue numbers (human kidney, PBMC means) from the supplementary tables.
- I found no published head-to-head comparison of per-sample vs pooled MAD thresholds with numbers. The per-sample recommendation rests on expert consensus (Heumos 2023, scverse).
- AMP Phase 2 (10x) lupus nephritis QC thresholds were not retrieved.

## 2. Ambient RNA removal: SoupX vs CellBender vs DecontX

### Takeaway
On genotype ground truth from mouse kidney, CellBender gives the most accurate ambient estimates and the largest marker-specificity gains. SoupX underestimates contamination in whole-cell data. Clustering and classification are robust unless contamination is high, and correction can distort fine structure at low contamination. A 2026 benchmark of 7 tools ranks CellBender first, DecontX (full mode) second and SoupX (reduced mode) third. That benchmark advises SoupX when contamination is not evident or only filtered matrices exist, because CellBender and DecontX can over-correct.

### Cited Findings
- Janssen et al., Genome Biol 2023: mouse kidney data, three scRNA-seq and two snRNA-seq replicates, each pooling CAST/EiJ, C57BL/6J and 129S1 mice so that SNPs give ground truth. Background was about 3–35% of UMIs per cell. Replicate medians were about 3%, 11% (whole cell), 17% and 35% (nuclei). Under 10% of background came from barcode swapping; the rest resembled empty-droplet (ambient) profiles. — [Janssen 2023 Genome Biol](https://link.springer.com/10.1186/s13059-023-02978-x)
- Janssen 2023, method accuracy: **CellBender** was most accurate and captured cell-to-cell variation, but cost about 45 CPU hours per dataset. **SoupX** underestimated background in scRNA-seq and did not capture per-cell variation. **DecontX** was comparable to CellBender on whole cells, worse on nuclei, and failed on default settings for one nuclei sample. — [Janssen 2023](https://link.springer.com/10.1186/s13059-023-02978-x)
- Kidney marker example from Janssen 2023: the PT marker Slc34a1 was detected in about 60% of non-PT cells before correction. After correction it was 54% (SoupX), 7% (CellBender), 9% (DecontX with empty-droplet profile) and 2% (DecontX default). In nuc2, DecontX left 87% and CellBender reduced it to 4%. — [Janssen 2023](https://link.springer.com/10.1186/s13059-023-02978-x)
- Janssen 2023, downstream effects: correction changed under 1% of whole-cell labels and up to 1.3% in nuc2. Recommendation: always correct for marker or DE analysis, preferably with CellBender. For clustering, classification or pseudotime, correct only when background is high. — [Janssen 2023](https://link.springer.com/10.1186/s13059-023-02978-x)
- Cargnelli, Nielsen & Madsen, bioRxiv Jan 2026 (not peer-reviewed): 7 methods (CellBender, DecontX, FastCAR, scAR, scCDC, SoupX, CellClear) on 49 real samples (97,357 cells and nuclei) plus simulations. 81.56% of cells contained ambient RNA. Mean ambient fraction was 11.80% in nuclei vs 7.79% in whole cells. Ambient load correlated with empty-droplet expression (R = 0.857). The top 10% most frequent features contributed 82.30% of ambient RNA. — [Cargnelli 2026 bioRxiv](https://www.biorxiv.org/content/10.64898/2026.01.13.699237.full.pdf)
- Cargnelli 2026, scores: scAR removed the most ambient RNA (0.679) but also removed 0.445 of endogenous signal, and it worsened integration. CellBender removed 0.604 and was the only method that improved clustering of the high-contamination sample. Overall ranking: CellBender, then DecontX (full), then SoupX (reduced); "no method is universally best." CellBender needs a GPU. Only SoupX reduced mode works with filtered matrices alone. — [Cargnelli 2026](https://www.biorxiv.org/content/10.64898/2026.01.13.699237.full.pdf)
- scverse book: ambient correction must be run separately for every sample. It is optional and depends on contamination level. The book's notebook uses SoupX. — [sc-best-practices QC](https://www.sc-best-practices.org/preprocessing-visualization/quality-control/). Heumos 2023 says one should "consider removing ambient RNA as an initial analysis step" for many tissues, naming SoupX and CellBender. — [Heumos 2023](https://www.nature.com/articles/s41576-023-00586-w.pdf)
- KPMP used SoupX v1.5.0 on 10x kidney whole-cell data. — [Lake 2023](https://www.nature.com/articles/s41586-023-05769-3.pdf)
- 10x suggests CellBender when neutrophils do not separate from background. — [10x neutrophil tutorial](https://www.10xgenomics.com/support/software/cell-ranger/tutorials/cr-tutorial-neutrophils)

### Inferences
- Kidney biopsy digests (enzymatic digestion, fragile tubule cells) are the high-contamination case. Ambient PT transcripts (e.g. SLC34A1, LRP2) and haemoglobin will bleed into immune cells. CellBender, run on the raw `raw_feature_bc_matrix`, is favoured, especially before marker or DE analysis of infiltrating immune cells.
- In lupus PBMC, the likely soup components are haemoglobin (HBB/HBA from RBC lysis), platelet genes (PPBP/PF4) and, in active SLE with plasmablast expansion, immunoglobulin transcripts (IGHG/IGKC/IGLC). Ig soup could create false "Ig+" signal in T, NK or myeloid cells, and could distort B-cell isotype or plasmablast scores. This is my inference: none of the benchmarks I read tested Ig soup directly.
- For PBMC with low contamination, SoupX (or no correction) avoids CellBender or DecontX over-correction (Cargnelli 2026). Measure first: check the HB and Ig fraction in non-erythroid, non-B clusters.

### Gaps
- I found no benchmark specific to human PBMC or human kidney biopsy ambient RNA; Janssen 2023 is mouse kidney.
- I found no published quantification of Ig ambient RNA in lupus PBMC.

## 3. Doublet detection: scDblFinder vs Scrublet vs DoubletFinder; multiplexed data

### Takeaway
Xi & Li 2021 compared 9 methods on 16 real datasets. scDblFinder was not one of them. DoubletFinder had the best mean AUPRC, and it was only 0.537. scDblFinder, added in the benchmark's later update and in Germain et al. (F1000Research 2021/2022), has the top mean AUPRC on those same datasets, which its parameters were tuned on. The scverse book and Heumos 2023 recommend scDblFinder, run per sample (per 10x capture lane). In genotype-multiplexed lupus PBMC (Perez 2022), donor demultiplexing (freemuxlet) removed a mean 22% doublet rate, and Scrublet then removed a further 67,969 droplets.

### Cited Findings
- Xi & Li, Cell Systems 2021 (classic): the 9 methods were doubletCells, Scrublet, cxds, bcds, hybrid, Solo, DoubletDetection, DoubletFinder and DoubletDecon, tested on 16 real datasets with experimentally annotated doublets plus synthetic data. "The highest mean AUPRC value, which was achieved by DoubletFinder, is only 0.537." Solo had the highest mean AUROC. The best per-dataset AUPRC ranged from 0.239 to 1.000. cxds was fastest and the only method robust to splitting data into batches. DoubletFinder, DoubletDetection and Scrublet did well for clustering and trajectories. "No method is universally the best." — [Xi & Li 2021](https://escholarship.org/content/qt6gr648np/qt6gr648np_noSplash_32134c0e501d9fdf306702f2a4dcc052.pdf)
- Germain et al., scDblFinder, F1000Research 2021 (v2 2022): scDblFinder has the highest AUPRC on most of the 16 Xi & Li datasets and the top two mean AUPRCs, at a fraction of the runtime. The authors note that parameters were tuned on these same datasets, which biases the comparison. At the automatic threshold, about 87% of heterotypic doublets are recovered, with a real FDR of about 17% in one sample. Processing samples separately was efficient and usually as accurate as processing them jointly; joint processing can bias against samples with different cell numbers or depth. Cluster-based and random artificial-doublet modes perform comparably. Doublet-origin calls are unreliable: both parents were correct for only 20% of doublets. — [Germain 2022 F1000Research](https://f1000research.com/articles/10-979/v2)
- scverse book: use scDblFinder, which reportedly had "the highest doublet detection accuracy and a good computational efficiency and stability" in the Xi & Li addendum. Detection "has to be run per sample" and "should not be run on aggregated scRNA-seq data representing multiple batches." The book suggests flagging doublets and inspecting them rather than deleting them immediately. — [sc-best-practices QC](https://www.sc-best-practices.org/preprocessing-visualization/quality-control/)
- Heumos 2023: scDblFinder outperformed other methods in several benchmarks. "It can be beneficial to apply multiple doublet detection methods and compare the results." — [Heumos 2023](https://www.nature.com/articles/s41576-023-00586-w.pdf)
- Perez et al., Science 2022 (classic, lupus PBMC): 10x 3′ v2 data from 264 samples (162 SLE including flares, 99 controls) plus 91 replicates, in 23 pools across 4 batches. Donors were demultiplexed by genotype with freemuxlet (mean doublet rate 22.12%). Scrublet removed 67,969 further droplets. Platelets, megakaryocytes and RBCs were removed by markers (112,805 cells). The final dataset had 1,263,676 cells, mean 3,560 singlets per sample. Cell filtering used "standard approaches in Scanpy 1.6"; numeric thresholds are in the supplement, which I did not read. — [Perez 2022 PMC9297655](https://pmc.ncbi.nlm.nih.gov/articles/PMC9297655)
- Doublet handling in kidney atlases: KPMP used DoubletDetection v2.4.0 (snRNA-seq) and a 5,000-gene cap (scRNA-seq). Arazi 2019 used FSC/SSC singlet gates in sorting plus a 5,000-gene cap. — [Lake 2023](https://www.nature.com/articles/s41586-023-05769-3.pdf); [Arazi 2019](https://hacohenlab.mgh.harvard.edu/wp-content/uploads/2023/06/Arazi-et-al-Nat-Immun-2019.pdf)

### Inferences
- Run scDblFinder per capture (GEM well). In genotype-pooled data, use the genotype-detected (between-donor) doublets as a calibration check, then run scDblFinder to catch within-donor heterotypic doublets. Perez 2022 did the same with Scrublet.
- Superloaded pools like Perez's (22% doublets) require genotype demultiplexing. For single-donor kidney biopsies, use scDblFinder, flag rather than delete, and then remove clusters that co-express markers of two lineages (for example PT plus T-cell markers). KPMP reasoned the same way about "ambiguous" clusters.
- Hard gene caps (5,000 genes) are crude. Their effect depends on chemistry and depth (v3 and GEM-X detect more genes), so do not copy them.

### Gaps
- I did not retrieve exact scDblFinder AUPRC values from the Xi & Li addendum.
- I found no 2023–2026 independent doublet benchmark that includes scDblFinder. One may exist; I did not search exhaustively.

## 4. Dissociation stress / ex vivo artefacts, cell cycle and sex: regress or not

### Takeaway
Dissociation and ex vivo handling induce an immediate-early and heat-shock programme (FOS, JUN, HSPA1A and others) that is tissue-specific. Recommended practice is to score and flag it, not regress it blindly. Wet-lab fixes (transcription inhibitors, cold proteases) reduce it at source. For cell cycle, the consensus is to score cells with built-in Scanpy or Seurat functions and to regress only if cycling obscures the question, because proliferation is itself biology. For a lupus study this includes cycling plasmablasts and proliferating T cells. Arazi 2019 found that regressing technical covariates had a negligible effect.

### Cited Findings
- A 2026 Communications Biology paper on neuroimmune tissue: enzymatic digestion, temperature shifts and mechanical forces induce "stress-altered transcription" (e.g. Fos, Jun, Hspa1a), and the signature is tissue-specific (brain differs from skull bone marrow). An actinomycin D + triptolide + anisomycin cocktail cut artefactually activated cells by about 95%. A random-forest classifier on tissue-specific panels (53 and 55 genes) reduced artefact signal by 60–64%. Caveat: the inhibitors may mask genuine rapid responses. — [Springer Nature "Behind the paper" blog on Commun Biol 2026](https://communities.springernature.com/posts/the-noise-that-wasn-t-biology-cleaning-up-single-cell-data). This is a secondary source; I did not read the primary paper.
- A Genome Biology 2019 paper on dissociation artefacts (cold-active protease vs warm collagenase in solid tissues) surfaced in search, but I did not open it. — [Genome Biol s13059-019-1830-0](https://link.springer.com/article/10.1186/s13059-019-1830-0)
- Arazi 2019 removed two dividing kidney subclusters for high mito and stress-gene expression. Regressing plate, UMI and mito "had a negligible effect" and was not applied. — [Arazi 2019](https://hacohenlab.mgh.harvard.edu/wp-content/uploads/2023/06/Arazi-et-al-Nat-Immun-2019.pdf)
- Heumos 2023 on cell cycle: effects can be removed, but "knowing whether cells are cycling may provide valuable insights into the underlying biology." A benchmark recommends Scanpy or Seurat built-in labelling and correction as a baseline, with Tricycle for highly heterogeneous data. — [Heumos 2023](https://www.nature.com/articles/s41576-023-00586-w.pdf)
- Heumos 2023 on normalisation: the shifted logarithm works well for variance stabilisation before dimensionality reduction. Do not apply the shifted log to CPM-scale data. — [Heumos 2023](https://www.nature.com/articles/s41576-023-00586-w.pdf)

### Inferences
- Lupus PBMC: cryopreserved and thawed PBMCs, and differences in time-to-processing between sites, can induce FOS/JUN/HSP and IFN-unrelated stress programmes. Score a stress module per cell and test it as a covariate by site, batch and disease. Do not regress it, because regression can strip real activation signal (e.g. NR4A1 or EGR in activated T cells).
- Kidney biopsy: tubule stress and injury programmes (HAVCR1, LCN2, SPP1) overlap with real lupus nephritis injury biology. KPMP kept these "altered states." Flag them; do not filter them.
- Sex: XIST and chrY genes (RPS4Y1, DDX3Y) should be used to verify donor sex and to check genotype-demultiplexing assignments. Do not regress them; exclude them from HVGs if sex confounds clustering. This is standard practice from general knowledge, not from a source fetched here.

### Gaps
- I did not retrieve the van den Brink 2017 (Nat Methods, classic) or Denisenko 2020 (Genome Biol, classic) dissociation signature gene lists.
- I found no source quantifying sex effects on scRNA-seq QC or giving a regress-vs-not recommendation for sex.
- I found no lupus-specific study of stress artefacts from cryopreservation in PBMC or kidney.

## 5. Low-quality cells vs biologically distinct states: over-filtering risks

### Takeaway
Uniform thresholds systematically remove real populations: low-RNA neutrophils, mito-high proximal tubule cells, injured or "altered" tubule states, low-quality but real epithelial cells in biopsies, and (by inference) high-RNA plasmablasts, which outlier caps on counts and top-20-gene % can flag. Use permissive per-sample MAD thresholds, then re-check every filtered cell against annotation, as both scverse and KPMP do.

### Cited Findings
- Neutrophils express only a few hundred genes and look like background, so default Cell Ranger cell calling drops them. In the 10x example, default calling gave 1,886 cells against about 8,000 expected, and neutrophils were later recovered as 3,343 cells. 10x recommends `--force-cells` (overestimate), `--include-introns` (default from v7) and a feature-count threshold between neutrophils and empties (log2 7.5 in one example, 6 for 5′). Confirm identity with FCGR3B. GEM-X chemistry may not need forcing. — [10x neutrophil tutorial](https://www.10xgenomics.com/support/software/cell-ranger/tutorials/cr-tutorial-neutrophils)
- Neutrophils and macrophages are among human cell types with mean mito above 5% (Osorio & Cai). — [Osorio & Cai preprint](https://scispace.com/pdf/systematic-determination-of-the-mitochondrial-proportion-in-rc82n0wqis.pdf)
- KPMP flagged and kept altered or injured states (low gene count, high mito or ER transcripts, injury markers). It removed only ambiguous clusters that mixed major classes. — [Lake 2023](https://www.nature.com/articles/s41586-023-05769-3.pdf)
- Arazi 2019 found kidney epithelial cells had lower quality than leukocytes and kept them as one cluster rather than discarding them. — [Arazi 2019](https://hacohenlab.mgh.harvard.edu/wp-content/uploads/2023/06/Arazi-et-al-Nat-Immun-2019.pdf)
- Heumos 2023: permissive thresholds first. Analogous warning from CITE-seq: removing cells with low counts "may remove cell types" with intrinsically low signal. — [Heumos 2023](https://www.nature.com/articles/s41576-023-00586-w.pdf)
- scverse: high mito fraction can reflect real respiratory biology, so mito should be considered jointly with counts and genes. — [sc-best-practices QC](https://www.sc-best-practices.org/preprocessing-visualization/quality-control/)
- A fixed 10% mito cut removed 41% of cells from a human blood-derived (menstrual) 10x dataset, against 18.5% with the adaptive miQC. — [miQC](https://www.biorxiv.org/content/10.1101/2021.03.03.433798v1.full-text)

### Inferences
- Plasmablasts and plasma cells have very high Ig transcript content. They can exceed the upper 5 MAD bound on total counts and score very high on `pct_counts_in_top_20_genes`, which IG genes dominate. Check that MAD outlier flags are not enriched in the MZB1+/JCHAIN+ cluster before removing them. Plasmablast expansion is a key SLE activity readout, so this is my inference and a priority check.
- PBMC from Ficoll preparations contains few neutrophils, but low-density granulocytes (LDGs), which are relevant in SLE, co-sediment with PBMC. They need the neutrophil-style cell-calling handling above. This is my inference from lupus biology; no source was fetched here.
- Recommended order: per-sample cell calling (EmptyDrops or Cell Ranger, forced for granulocytes), then ambient correction, then permissive MAD flags, then doublet flags, then annotation, then a re-check of flagged cells by cluster, then final removal. In the re-check, keep flagged cells that form coherent, marker-defined populations.

### Gaps
- I found no published quantification of plasmablast loss under MAD filtering.
- I found no lupus-specific study of LDG retention under 10x QC.
- AMP Phase 2 lupus nephritis 10x QC thresholds were not retrieved.
