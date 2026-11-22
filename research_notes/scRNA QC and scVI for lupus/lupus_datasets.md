# Public lupus single-cell datasets: design, author processing, confounders and open availability (for a scVI re-analysis)

Retrieval date: 2026-11-22. GEO pages blocked automated browsing (reCAPTCHA), so GEO facts come from NCBI E-utilities / GEO SOFT text (`acc.cgi?...&form=text`) and the GEO FTP `filelist.txt`. CELLxGENE facts come from the Discover curation API and from reading the `obs` table of the published h5ad files remotely over HTTP range requests (no full download). Counts marked "(obs)" were computed by me from those h5ad files.

## Q1. Perez et al. 2022 Science (1.26M PBMCs, 162 SLE / 99 controls): what is open, how was it processed, what confounds it?

### Takeaway
The full cell-level count matrix (raw counts in `raw.X`) is openly downloadable from CELLxGENE as one 12.2 GB h5ad. The open donor metadata are disease (SLE/normal), `disease_state` (managed/flare/treated/na), ancestry, sex, age, `Processing_Cohort` (4 batches), `library_uuid` (88 10x wells), sample/suspension IDs and author cell types. Nephritis, SLEDAI and medication data are not in the open file. Processing cohort is badly confounded with disease state: cohort 1 is all controls, and every flare and post-flare-treated sample is in cohort 3. A scVI model that uses `Processing_Cohort` as a batch key will therefore partly remove the flare signal.

### Cited Findings
**Design and technology**
- 1.2M+ PBMCs from 162 cases and 99 controls, profiled by multiplexed scRNA-seq ("mux-seq") — [GEO GSE174188 summary (E-utilities)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188)
- "The 264 samples corresponded to 162 SLE cases, including 19 disease flare cases and 10 matched samples post–flare treatment, along with 99 healthy controls". "The 264 samples and 91 replicates were profiled in 23 pools across four batches" — [Perez et al. 2022, PMC9297655](https://pmc.ncbi.nlm.nih.gov/articles/PMC9297655/) (DOI [10.1126/science.abf1970](https://doi.org/10.1126/science.abf1970))
- PBMCs were Ficoll-isolated and cryopreserved (FBS + 10% DMSO). Frozen PBMCs were pooled (Kang et al. demuxlet design) and run on 10x Chromium Single Cell 3' v2. Cells in batches 3 and 4 were also stained with AbSeq/CITE antibodies — [GEO GSE174188 SOFT, Sample_extract_protocol](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188)
- Surface protein was profiled in processing batches 3 (155,034 cells) and 4 (375,261 cells), with 16 and 99 DNA-conjugated antibodies respectively — [Perez 2022 PMC9297655](https://pmc.ncbi.nlm.nih.gov/articles/PMC9297655/)
- The GEO series has 88 GSM samples, one per 10x well (titles like `CLUES1-POOL01-1`, "10x well: 1–4"). 64 are RNA only and 24 are RNA+ADT. The genome build is hg19. GEO has no supplementary matrix (the GEO FTP `suppl/` directory returns 404); the processed data live on CELLxGENE/HCA — [GEO GSE174188 SOFT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188)

**Authors' QC, doublets and integration**
- Pipeline: Cell Ranger, then Freemuxlet to assign donor and remove doublets (mean doublet rate 22.12%), then Scrublet (removed 67,969 more droplets), then removal of platelets and red blood cells (112,805 cells). 1,263,676 cells remained — [Perez 2022 PMC9297655](https://pmc.ncbi.nlm.nih.gov/articles/PMC9297655/)
- Platelet contamination: cells with >1 count of PF4, SDPR, GNG11 or PPBP were removed. RBCs were removed as a Louvain community with high haemoglobin expression — [GEO GSE174188 SOFT data_processing](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188)
- "Technical variation was removed using COMBAT and regressing out nUMIs, and mitochondrial percent". ADTs were CLR-normalised and used only to confirm annotations — [Perez 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9297655/); [GEO SOFT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188)
- No explicit per-cell gene/UMI/mito thresholds are given in the GEO processing text or in the main-text methods I could access. No ambient-RNA correction (SoupX/CellBender) is mentioned — [GEO SOFT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188)
- Authors' claims: "Different pools and processing batches had no observable effects on the distribution of cells". Cell-type percentages were reproducible across biological replicates (median Pearson R 0.79 for cases, 0.85 for controls) — [Perez 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9297655/)
- Batch was a covariate in the eQTL models — [Perez 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9297655/)

**Open availability (CELLxGENE)**
- Collection `436154da-bcf1-4130-9c8b-120ff9a888f2`, dataset `218acb0f-9f2f-4f76-b90b-15a4b7c7f629`: 1,263,676 cells, 261 donors, 10x 3' v2, one H5AD asset of 12.22 GB. Raw counts are in `raw.X`. Linked resources: dbGaP phs002812 (genotypes), GEO GSE137029, an HCA project, Zenodo 4724043 and the UCSC Cell Browser `lupus-pbmc` — [CELLxGENE curation API, collection 436154da](https://cellxgene.cziscience.com/collections/436154da-bcf1-4130-9c8b-120ff9a888f2)
- Download URL: `https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad` — [CELLxGENE curation API](https://api.cellxgene.cziscience.com/curation/v1/collections/436154da-bcf1-4130-9c8b-120ff9a888f2)
- Open obs columns (obs): `Processing_Cohort`, `disease`, `disease_state`, `donor_id`/`ind_cov`, `library_uuid` (88 levels), `sample_uuid` (274), `suspension_uuid` (291), `self_reported_ethnicity`, `sex`, `development_stage` (age in years), `author_cell_type` (11 types: T4, cM, B, T8, NK, ncM, Prolif, cDC, pDC, PB, Progen), `ct_cov` (14 finer subtypes), `cell_state`. The `obsm` slots are X_pca and X_umap. Gene annotation is GENCODE 19 — [CELLxGENE h5ad c55dc602](https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad)
- Donor counts per `disease_state` (obs): na (controls) 99 donors / 486,418 cells; managed 146 / 696,626; flare 19 / 55,120; treated 10 / 25,512 — [CELLxGENE h5ad](https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad)
- Ancestry (obs): European American 149 donors, Asian 107, African American 3, Hispanic/Latin 2. Sex: 244 female, 17 male — [CELLxGENE h5ad](https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad)
- `Processing_Cohort` (obs): 1.0 = 47 donors / 175,273 cells / 12 libraries; 2.0 = 142 / 558,108 / 40; 3.0 = 41 / 155,034 / 12; 4.0 = 96 / 375,261 / 24. 57 donors appear in more than one cohort, which gives cross-batch bridging replicates. 11 donors have more than one `sample_uuid` — [CELLxGENE h5ad](https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad)
- Pooling (obs): 7–19 donors per library, mostly 16. Cells per library range 5,080–17,525 (median 14,577) — [CELLxGENE h5ad](https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad)

**Cohort × disease_state × ancestry, in donors (obs)**
- Cohort 1: controls only (46 European American, 1 Asian).
- Cohort 2: 120 managed SLE (63 Asian, 57 European American) and 22 controls.
- Cohort 3: all 19 flare and all 10 treated donors, plus 4 managed and 18 controls; this is the only cohort with African American and Hispanic donors.
- Cohort 4: 52 managed (26 Asian / 26 European American) and 44 controls (22/22), i.e. balanced.

Source: [CELLxGENE h5ad](https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad)

**Labels that are not open**
- Genotypes are in dbGaP phs002812.v1.p1 (controlled) — [Perez 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9297655/)
- A 2026 Commun Biol re-analysis split GSE174188 into "99 HCs, 68 SLE, and 76 LN samples" and thanked Perez "for sharing sequencing data and clinical information". This confirms that LN status was obtained by request and is not in the open file — [Commun Biol 2026, PMC12868680](https://pmc.ncbi.nlm.nih.gov/articles/PMC12868680/) (DOI [10.1038/s42003-025-09431-8](https://doi.org/10.1038/s42003-025-09431-8))

**Re-analyses of this dataset**
- SLECA (Comput Struct Biotechnol J 2026): integrated Perez with other SLE datasets using scVI v1.1.6 on raw UMIs, and noted that "residual batch effects may persist across studies" — [PMC13358165](https://pmc.ncbi.nlm.nih.gov/articles/PMC13358165/) (DOI [10.34133/csbj.0163](https://doi.org/10.34133/csbj.0163))
- "Transcriptional patient map of SLE" (bioRxiv 2026): reprocessed Perez "following the original publication", used it as the reference for Scanpy `ingest` label transfer of Nehar-Belaid GSE135779, and modelled the four Perez processing batches as separate batch groups (e.g. "Pérez 1.0, N = 27") — [PMC13142532](https://pmc.ncbi.nlm.nih.gov/articles/PMC13142532/) (DOI [10.64898/2026.04.28.721379](https://doi.org/10.64898/2026.04.28.721379))
- Commun Biol 2026: used Scanpy with Harmony and BBKNN for batch correction, with LN labels obtained from the authors — [PMC12868680](https://pmc.ncbi.nlm.nih.gov/articles/PMC12868680/)
- Paek et al. 2025 (EMM) used Perez as external validation — [PMC12411632](https://pmc.ncbi.nlm.nih.gov/articles/PMC12411632/)
- Other papers reusing GSE174188 (single-cell eQTL methods, PARAFAC2-RISE in Cell Systems 2025, scDesignPop, "Private information leakage from single-cell count matrices" in Cell 2024) appear in a PMC search for "GSE174188 AND lupus" (44 hits) — [PMC search via E-utilities](https://pmc.ncbi.nlm.nih.gov/search/?term=GSE174188+lupus)

### Inferences
- In Perez, healthy vs SLE contrasts are only partly separable from batch: cohort 1 has no cases, and cohort 2 is 85% cases. Flare vs managed is almost fully nested in cohort 3. Cohort 4 is the cleanest within-batch case/control comparison; the authors themselves used "processing batch 4 cases only" in one analysis.
- For scVI, a reasonable batch key is `library_uuid` (88 wells; donors are pooled within a well, so disease is not nested in a well), with `Processing_Cohort` as an extra categorical covariate. Batch-corrected latents should be checked for flare signal loss. The 57 bridging donors seen in more than one cohort can be used to measure how well batches are corrected.
- Perez already removed doublets, platelets and RBCs, so the open matrix is post-QC. Re-QC is limited to re-applying thresholds (for example MAD-based) on `raw.X`. Ambient-RNA correction is not possible without the raw droplet matrices.
- The 3' v2 chemistry and hg19/GENCODE 19 annotation differ from the 5' v2 GRCh38 datasets below. Genes should be harmonised by Ensembl ID, and scVI should be given assay/chemistry as a batch covariate.
- The "treated" donors are re-samples of flare patients (10 matched post-treatment samples), so `donor_id` spans states. 146 + 19 + 10 = 175 is more than 162 SLE donors, so some donors carry more than one state label; check this per donor before any donor-level split.

### Gaps
- I did not obtain the exact per-cell QC thresholds (genes, UMIs, mito %) from the supplementary methods; the GEO text and main methods do not state them.
- I could not confirm whether the open `raw.X` contains all genes or the filtered gene set.
- The CELLxGENE obs has no explicit replicate flag. Which of the 274 samples are the "91 replicates" vs the 264 primary samples has to be inferred from `sample_uuid`/`donor_id` duplication.

## Q2. CD19 CAR-T in SLE (GSE263931; Wilhelm, Krönke, Schett, Mackensen group, Erlangen)

### Takeaway
GSE263931 is open on GEO: per-patient 10x 5' v2 PBMC count matrices (MTX/TSV) for 7 patients before and after CD19 CAR-T, plus matched TCR libraries. Authors' QC thresholds were not reported numerically. Two things limit the dataset: post-treatment samples were taken after B-cell reconstitution (9–28 weeks), and patients received fludarabine/cyclophosphamide lymphodepletion, which the authors say they cannot rule out as a contributor.

### Cited Findings
- Title: "Selective CAR-T cell mediated B cell depletion suppresses interferon signature in SLE [PBMC]". 28 samples, released 2024-05-10, PMID 38722688. Contributors: Wilhelm, Chambers, Müller, Bozec, Grieshaber-Bouyer, Winkler, Mougiakakos, Mackensen, Schett, Krönke. BioProject PRJNA1100404 — [GEO GSE263931 (E-utilities/SOFT)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263931)
- Samples are Patient1–7 × pre/post × (scRNAseq, TCRseq): 14 GEX ("cell type: PBMCs") and 14 TCR ("cell type: T cells"). Treatment labels: "pre/post anti CD19 CART-cell therapy". SepMate/Ficoll PBMCs, Chromium Next GEM 5' v2, `cellranger multi` v7.0.0, hg38 — [GEO GSE263931 SOFT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263931)
- Files: `GSE263931_RAW.tar` (557 MB) with per-sample `*_barcodes.tsv.gz`, `*_features.tsv.gz`, `*_matrix.mtx.gz`. Per-sample `*_samples.csv.gz` files (about 145 bytes) were added on 2026-11-13 — [GEO FTP filelist.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263931/suppl/filelist.txt)
- QC: "we applied appropriate thresholds for UMIs, number of detected genes, and number of UMIs assigned to mitochondrial genes" (no numbers given). Seurat v4.3. "All single-patient data sets were integrated into 1 object". The integration method is not named beyond "the Seurat pipeline". About 40,700 PBMCs after integration. Clustree was used for resolution. scRepertoire v1.10.0 for TCR/BCR — [Wilhelm et al. 2024 JCI Insight, PMC11383166](https://pmc.ncbi.nlm.nih.gov/articles/PMC11383166/) (DOI [10.1172/jci.insight.179433](https://doi.org/10.1172/jci.insight.179433))
- No doublet or ambient-RNA step is described — [PMC11383166](https://pmc.ncbi.nlm.nih.gov/articles/PMC11383166/)
- Timing: "B cell repopulation was observed in the individual patients between 9 and 28 weeks after CD19 CAR T cell treatment". The post samples are "after early B cell reconstitution", and "Among the data sets obtained after B cell reconstitution, we could not detect CAR T cells" — [PMC11383166](https://pmc.ncbi.nlm.nih.gov/articles/PMC11383166/)
- Confounder in the authors' words: "we cannot exclude that the lymphotoxic conditioning therapy with fludarabine and cyclophosphamide ... contributes to the beneficial effects" — [PMC11383166](https://pmc.ncbi.nlm.nih.gov/articles/PMC11383166/)

### Inferences
- The raw MTX files are Cell Ranger filtered (cell-called) matrices, not raw droplet matrices, so SoupX/CellBender would need the unfiltered matrices, which were not deposited. Doublet detection (scDblFinder/Solo) can still run per sample.
- Patient is perfectly paired with timepoint, so patient should be the scVI batch key (or sample key with patient as a covariate). Pre vs post can then be tested within patient.
- Clinical labels (SLEDAI, LN) are only in the paper's Supplemental Table 1, not in GEO.

### Gaps
- Fresh vs cryopreserved handling is not stated in the GEO text I retrieved.
- I did not read the contents of the new `*_samples.csv.gz` files, which may be per-sample metadata.
- No independent re-analysis of GSE263931 was identified in this pass.

## Q3. CD3×CD19 bispecific A-319 phase 1 in SLE (Nat Med 2026; CELLxGENE 6c147b03)

### Takeaway
Cell-level data are openly on CELLxGENE as 5 lineage-split h5ad files (T, B, LDG/neutrophil progenitors, myeloid, monocytes; 10x 5' v2; 11 donors; raw counts in `raw.X`). Each has Baseline / Day 28 / Month 3 timepoints, dose, a composite disease label (including LN) and IFN module scores. There are no healthy controls, and raw reads are in China's GSA-human (HRA017319), which is controlled access.

### Cited Findings
- Paper: Xu J, Mei C, Guan X, et al. "A bispecific CD3×CD19 antibody for systemic lupus erythematosus: a phase 1 trial". Nat Med 2026;32(9):3372–3384, PMID 42722888. 12 patients with active SLE, IV arm. A-319 at 0.3–1.2 µg/kg three times weekly for 3 weeks after 1 week of priming (0.05 µg/kg). 52 weeks follow-up — [PubMed 42722888](https://pubmed.ncbi.nlm.nih.gov/42722888/) (DOI [10.1038/s41591-026-04572-7](https://doi.org/10.1038/s41591-026-04572-7))
- Collection `6c147b03-08e4-4bb7-afb7-2bcc0a814b65`, published 2026-09-21, 5 datasets: T cells 112,040 cells (1.08 GB); B cells 4,829; low-density granulocytes and neutrophil progenitors 77,112; myeloid 167,156 (1.1 GB); monocytes 76,421. All 10x 5' v2, blood, 11 donors, raw counts in `raw.X`. The only external link is GSA-human HRA017319 — [CELLxGENE collection 6c147b03](https://cellxgene.cziscience.com/collections/6c147b03-08e4-4bb7-afb7-2bcc0a814b65)
- T-cell obs (obs): `timepoint` ∈ {Base 27,103 cells; D28 47,122; M3 37,815}, 11 donors at each. `sample` IDs are S01_Base…S17_M3 (S14 has an M13 sample instead of M3). `dose` ∈ {0.3 µg/kg: 6 donors; 0.6: 2; 1.2: 3}. `author_disease` includes SLE, SLE/LN, SLE/LN/sSS, SLE/LN/NPSLE, SLE/sSS, SLE/RA and SLE/SSc; 3 donors carry LN. All donors are Asian and female, aged 20–54. The file has TCR clonotype columns (`raw_clonotype_id`, `cdr3s_nt`, `t.clone.timepoints`), many published IFN module scores (Banchereau, Kirou, Feng, etc.), `nCount_RNA`, `nFeature_RNA` and `percent.mito` — [CELLxGENE h5ad 00260ebe](https://datasets.cellxgene.cziscience.com/00260ebe-ef0d-400f-b937-8869459e2940.h5ad)
- Observed QC floors in the deposited T-cell object (obs): `nCount_RNA` min 501, median 4,450; `nFeature_RNA` min 228; `percent.mito` max 9.99 — [CELLxGENE h5ad 00260ebe](https://datasets.cellxgene.cziscience.com/00260ebe-ef0d-400f-b937-8869459e2940.h5ad)
- No batch_condition field is declared on the datasets — [CELLxGENE curation API](https://api.cellxgene.cziscience.com/curation/v1/collections/6c147b03-08e4-4bb7-afb7-2bcc0a814b65)

### Inferences
- The observed floors imply the authors filtered at about UMI > 500 and mito < 10% (Seurat-style column names). This is inferred from the data, not stated in the methods.
- The lineage split means the 5 files must be concatenated (on `observation_joinid`/barcodes) to rebuild whole-PBMC composition. The deposited "LDG" compartment shows granulocytes were retained, which is not standard for PBMC preps; this may reflect sample handling.
- 11 of 12 trial patients are deposited. Ethnicity is uniform (Asian), so the data cannot test ancestry effects. The Base→D28→M3 design gives within-patient trajectories similar to the CAR-T data but without lymphodepletion.

### Gaps
- The Nat Med article page could not be fetched, so the authors' stated QC, doublet, ambient and integration methods, and whether PBMCs were fresh or cryopreserved, are unconfirmed.
- Per-timepoint SLEDAI is not in the obs columns I read (only `author_disease`); it may be in the paper's supplement.

## Q4. GSE254176: 6 SLE patients (LN at onset), longitudinal before-flare / flare / remission

### Takeaway
Open on GEO: 19 PBMC 10x 5' v2 matrices (GEX + TCR + BCR) with per-sample disease state labels (Before flare 6, Flare 6, Remission 7) and sex. The 32 healthy controls come from AIDA (HCA-Asia) and are not in this GEO series. Authors filtered at >500 genes and <5% mito, used Scrublet, and integrated with Harmony.

### Cited Findings
- GEO: "Dynamic changes of immune cells in systemic lupus erythematosus with nephritis using longitudinal peripheral blood single-cell RNA-seq". 57 samples = 19 timepoints × (scRNA, scTCR, scBCR). Patients 1–4 have 2–6 timepoints each; patients 5 and 6 have 1 each. Released 2025-02-01. Platforms GPL18573 and GPL24676 (NextSeq 500 and NovaSeq 6000) — [GEO GSE254176 (E-utilities/SOFT)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254176)
- Characteristics counted over all 57 samples: disease state Remission 21, Flare 18, Before flare 18 (i.e. 7/6/6 timepoints); sex Female 48, Male 9 — [GEO GSE254176 SOFT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254176)
- Handling: frozen cells were thawed, filtered (40 µm) and counted with AO/PI; only samples with viability >70% were sequenced. Chromium Next GEM 5' v2, Chip K. Cell Ranger 6.0.0, GRCh38 — [GEO GSE254176 SOFT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254176)
- Files: `GSE254176_RAW.tar` (732 MB); per-sample MTX/TSV (GEX) and CSV (clonotypes). Sample library names look like `OM-HL-014-C8` and `OM-HL-014-C6-F` — [GEO FTP filelist.txt](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254176/suppl/filelist.txt)
- Authors' QC: cells with >500 genes and <5% mitochondrial content kept; platelets (PPBP, PF4, NRGN) removed; Rscrublet v0.1.0 plus manual marker-based doublet removal. TCR/IG V genes, ribosomal, mitochondrial and dissociation-stress genes were excluded from features. Harmony v0.1.0 was used on PCA "to eliminate technical discrepancies between samples". 153,483 cells after filtering — [Paek et al. 2025 Exp Mol Med, PMC12411632](https://pmc.ncbi.nlm.nih.gov/articles/PMC12411632/) (DOI [10.1038/s12276-025-01504-2](https://doi.org/10.1038/s12276-025-01504-2))
- Clinical: every patient had at least one flare with increased SLEDAI, and "All individuals initially presented with lupus nephritis". SLEDAI was 4–14 before flare and 9–19 at flare (Supplementary Table 1). Controls: "32 Korean participants ... from the Asian Immune Diversity Atlas" — [PMC12411632](https://pmc.ncbi.nlm.nih.gov/articles/PMC12411632/)
- Perez GSE174188 was re-processed by these authors "using a pipeline consistent with" theirs for validation — [PMC12411632](https://pmc.ncbi.nlm.nih.gov/articles/PMC12411632/)

### Inferences
- Without the AIDA controls this is a patients-only longitudinal set. To add controls, download the AIDA Korean donors from CELLxGENE/HCA; they are a different study, so their batch is confounded with disease.
- The C6/C8 suffix in library names may be a chip or lane identifier and should be tested as a batch covariate. This is a guess from naming only.
- Per-timepoint SLEDAI values are only in the supplement and must be hand-joined to GSM IDs.

### Gaps
- Whether the C6/C8 suffix means chip, lane or run is not stated.
- The time intervals between samples were not extracted.

## Q5. Lupus nephritis kidney single-cell data, and KPMP as a kidney reference

### Takeaway
The flagship LN kidney datasets are all controlled access and outside this plan: Arazi 2019 (ImmPort SDY997 / dbGaP phs001457), Der 2019 (AMP; no open repository found) and the AMP Phase 2 atlas (155 LN / 30 controls, ARK Portal/Synapse). Open kidney LN data are very thin. E-MTAB-13596 has only 1 LN biopsy (CD45+ immune cells; plus 5 ANCA and 1 nephrectomy control). GSE306378 is not kidney at all: it is PBMCs from 3 LN and 3 controls. KPMP v1.5 on CELLxGENE is open (AKI/CKD/healthy, scRNA + snRNA) but contains no LN.

### Cited Findings
**Arazi et al. 2019 Nat Immunol (AMP Phase 1)**
- 24 LN biopsies and 10 living-donor controls from 10 US sites. Intact tissue was cryopreserved, then processed in batches at a central site. Viable cells were sorted into 384-well plates for a modified CEL-Seq2 (2,736 leukocytes, 145 epithelial cells), plus a small droplet arm. Urine cells (577 high-quality cells from 8 patients) were also profiled. The authors state "the main sources of variability ... corresponded to cell types, rather than batch or technical factors" — [Arazi et al. 2019, PMC6726437](https://pmc.ncbi.nlm.nih.gov/articles/PMC6726437/) (DOI [10.1038/s41590-019-0398-x](https://doi.org/10.1038/s41590-019-0398-x))
- Two clusters with high mitochondrial and stress genes were excluded as low quality — [PMC6726437](https://pmc.ncbi.nlm.nih.gov/articles/PMC6726437/)
- Data: "deposited in the ImmPort repository (accession code SDY997)", with raw data in "dbGAP (accession code phs001457.v1.p1)". Both are controlled — [PMC6726437](https://pmc.ncbi.nlm.nih.gov/articles/PMC6726437/)
- There is an HCA Data Explorer project page for this paper (project 2d559a6e) — [HCA Data Explorer](https://explore.data.humancellatlas.org/projects/2d559a6e-7cd9-432f-9f6e-0e4df03b0888)

**Der et al. 2019 Nat Immunol (kidney + skin)**
- 21 LN renal biopsies, 17 with a matched non-lesional skin punch biopsy. The group's earlier study used the Fluidigm C1 96-well platform, which missed glomerular cells — [Der et al. 2019, PMC6584054](https://pmc.ncbi.nlm.nih.gov/articles/PMC6584054/) (DOI [10.1038/s41590-019-0386-1](https://doi.org/10.1038/s41590-019-0386-1))
- An HCA project titled "Single cell RNA sequencing to dissect the molecular heterogeneity in lupus nephritis" exists (project 4627f43e; probably Der 2017 JCI Insight) — [HCA Data Explorer](https://explore.data.humancellatlas.org/projects/4627f43e-a43f-44dd-8c4b-7efddb3f296d)

**AMP Phase 2 LN atlas (Sugiarto, Gurajala, Curtis, ..., Raychaudhuri; bioRxiv 2025, v3)**
- 538K scRNA and 140K snRNA profiles from biopsies of 155 LN patients and 30 pre-implantation transplant controls, plus 325K PBMC profiles. After QC: 155 LN / 30 controls with matched PBMC from 118 LN / 19 controls — [PMC13060338](https://pmc.ncbi.nlm.nih.gov/articles/PMC13060338/) (DOI [10.1101/2025.08.11.669754](https://doi.org/10.1101/2025.08.11.669754))
- QC removed cells with low UMIs, low gene counts and high mitochondrial reads (Table S3), plus Scrublet with marker-based doublet removal. Integration used Harmony on weighted PCs "to correct for processing batch, site, sample, and differences in scRNAseq and snRNAseq technologies". Cryopreserved-tissue enzymatic dissociation depleted non-immune cells in scRNA: 93.9% of immune profiles came from scRNA. Activity and chronicity explained cell states better than ISN class — [PMC13060338](https://pmc.ncbi.nlm.nih.gov/articles/PMC13060338/)
- Data: processed matrices are on the ARK Portal (Synapse): scRNA `syn54406421`, snRNA `syn54407831`. The preprint says they will be public on acceptance, under AMP policies — [PMC13060338](https://pmc.ncbi.nlm.nih.gov/articles/PMC13060338/)
- A related AMP study, "Deep profiling of lupus nephritis kidneys reveals dynamic changes in myeloid cells associated with disease progression" (Ann Rheum Dis 2026; bioRxiv 2025.12.17), has an ARK Portal dataset `syn74339161` — [ARK Portal](https://arkportal.synapse.org/Explore/Datasets/DetailsPage?id=syn74339161); [ARD](https://ard.eular.org/article/S0003-4967(26)00149-4/fulltext)

**E-MTAB-13596 (Amsterdam UMC; Vegting, Hilhorst)**
- "ScRNA-seq of human kidney immune cells of patients with ANCA-associated glomerulonephritis, Lupus Nephritis against a healthy nephrectomy control". 7 samples: 5 AGN, 1 LN, 1 nephrectomy control (papillary urothelial carcinoma). Fresh needle-core biopsy was processed immediately and FACS-sorted for live CD45+ cells. TotalSeq-C (130 antibodies) plus TCR/BCR in 2 samples (1 AGN and the control). Released 2024-06-28 — [BioStudies E-MTAB-13596](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13596)
- Open files: per-sample `*_matrix.mtx`, `*_barcodes.tsv` and `*_features.tsv` (AGN1–5, LN, NC), `feature_reference.csv`, SDRF/IDF. The LN matrix is 60 MB — [BioStudies FTP](https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/596/E-MTAB-13596/Files/)

**GSE306378 (Anhui University of Chinese Medicine; Cheng LL et al., BMC Med Genomics 2026)**
- This is PBMC, not kidney: 3 LN and 3 normal controls (titles SLE_rep1–3 and N_rep1–3), all "treatment: Immunosuppressant". 10x 3' v3, Cell Ranger 3.0.1, hg38. Seurat v3.0.2 QC: 200–6,000 genes per cell, genes detected in at least 3 cells, plus a mito filter. Integration with Seurat CCA anchors (FindIntegrationAnchors/IntegrateData). 69,069 cells. `GSE306378_RAW.tar` is 458 MB of MTX/TSV — [GEO GSE306378 SOFT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306378); [Cheng et al. 2026, PMC12895912](https://pmc.ncbi.nlm.nih.gov/articles/PMC12895912/) (DOI [10.1186/s12920-025-02277-3](https://doi.org/10.1186/s12920-025-02277-3))
- Data-processing text says reads were aligned to "the mouse genome, version hg38", a copy-paste error in the GEO record — [GEO GSE306378 SOFT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306378)

**KPMP / Lake et al. kidney atlas (open reference)**
- CELLxGENE collection `0f528c8a-a25c-4840-8fa3-d156fa11086f` (v1.5, 2024-05-13): scRNA 225,177 cells / 77 donors (10x 3' v3, batch_condition LibraryID, Run, 1.82 GB) and snRNA 304,989 nuclei / 80 donors (batch_condition library_id, 3.03 GB). Diseases: AKI, CKD, normal. Raw counts in `raw.X`. The original Lake et al. Nature 2023 collection is `bcb61471-2a44-4d00-a0af-ff085512674c`. KPMP repository: atlas.kpmp.org — [CELLxGENE KPMP v1.5](https://cellxgene.cziscience.com/collections/0f528c8a-a25c-4840-8fa3-d156fa11086f); [Lake et al. 2023](https://doi.org/10.1038/s41586-023-05769-3)
- A keyword scan of all CELLxGENE collection names and descriptions for "lupus" or "nephritis" (2026-11-22) returned only the Perez collection and the A-319 collection; there is no LN kidney collection on CELLxGENE — [CELLxGENE curation API collections list](https://api.cellxgene.cziscience.com/curation/v1/collections)

### Inferences
- Open LN kidney single-cell data at cell level is effectively 1 biopsy (E-MTAB-13596). A kidney arm in the scVI plan would need either a Synapse/ARK data-use application (AMP Phase 2; out of scope under the open-only constraint) or other GEO kidney series that this pass did not find.
- KPMP gives a healthy/AKI/CKD kidney reference. Its injured PT/TAL states (adaptive/maladaptive) are relevant to LN chronicity, but any LN-to-KPMP mapping would be cross-disease and cross-protocol: KPMP is mostly fresh dissociation (scRNA) plus frozen nuclei, while AMP used cryopreserved tissue.
- GSE306378 can be used as a small extra blood LN set (3' v3), not a kidney set.

### Gaps
- Der et al. 2019 data availability: the PMC text I parsed did not show the statement. It is probably AMP-controlled like Arazi, but this is unverified.
- I did not systematically search GEO for 2024–2026 LN kidney scRNA/snRNA series (e.g. childhood-onset LN spatial, Chinese LN biopsy cohorts). The bioRxiv spatial childhood-onset LN study (PMC10680641) appeared in search, but I did not check its data availability. A targeted GEO E-utilities query (`lupus nephritis[All] AND kidney AND single cell AND gse[ETYP]`) is the next step.
- I did not check whether any KPMP participants have LN (KPMP enrols AKI/CKD; LN is likely excluded, unverified).
