# PREREG addendum, round 103 (Q38, K1): kidney immune mapping, KPMP reference -> E-MTAB-13596 query

Frozen 2026-11-22 before any query expression value was read. Parent protocol: PREREG.md (K1 row: descriptive only, no hypothesis test).
Seen so far: KPMP obs (cell types, donors, libraries, disease); E-MTAB-13596 SDRF (7 samples: LN n=1 class II, ANCA-GN n=5, nephrectomy control n=1; CD45+ sorted, 10x 3' v3) and matrix headers (cell counts). No query gene values.

## Data (open, no account)
- Reference: KPMP adult kidney scRNA v1.5, CELLxGENE collection 0f528c8a, dataset 08648660 (225,177 cells, 77 donors, 10x 3' v3; normal/CKD/AKI; no LN).
- Query: E-MTAB-13596 (BioStudies, CC public), Cell Ranger filtered matrices for AGN1-5, LN, NC.

## Reference labels
- Immune cells (subclass.l1 in T, B, MYL) keep `author_cell_type`; labels with < 50 cells overall are merged into their l1 class.
- Non-immune cells get one label per subclass.l1 (PT, TAL, EC, ...), so contaminants are not forced into immune labels. Non-immune capped at 1,000 cells per l1 class (seed 0).
- Immune capped at <= 300 cells per donor x label (seed 0), per parent hardware note.

## Stage 1 QC (query, per library = sample; parent K1 column)
- 5 MAD on log1p counts, log1p genes, top-20 %; mito flagged only if > 3 MAD upper (logit) AND > 50 %. Doublets: scanpy scrublet per sample, flag and remove. Protected-population rule as round 102 (a cell label's flag rate > 3x the sample's overall rate -> keep its flagged cells). Ambient/stress scored, never regressed. Report removals per sample.

## Model (one seed, one config, frozen)
- 2,000 HVGs, seurat_v3 flavour, batch_key = LibraryID, on reference counts; query genes matched by Ensembl ID (missing genes zero-filled, count reported).
- scVI: n_latent 20, n_layers 2, layer norm, NB, dropout 0.2, encode_covariates, batch = LibraryID, max_epochs 200 with early stopping. scANVI from that scVI, max_epochs 20, 100 labelled samples per label. Query: scArches surgery (load_query_data, batch = sample), 100 epochs, weight_decay 0.

## Held-out sanity check (the trustworthiness gate; trivial baselines on held-out data)
- 20 % of reference donors held out (seed 0, donor-grouped). Models refitted on the other 80 %; held-out donors mapped as a query.
- Compare on held-out cells: (a) majority class; (b) kNN k = 15 on PCA-50 of log-normalised HVGs fitted on training cells; (c) scANVI. Metrics: accuracy and macro-F1 over immune labels.
- Gate: scANVI labels are used as the primary query labels only if scANVI macro-F1 >= kNN-PCA macro-F1 - 0.05 AND > majority. Otherwise kNN-PCA labels are primary and scANVI is reported secondary.

## Query report (descriptive, no pass rule, n = 1 LN)
- Per sample: label proportions among immune cells; fraction "unresolved" (kNN uncertainty = share of 15 nearest reference neighbours in the scANVI latent that disagree with the predicted label; > 0.5 = unresolved); median scANVI max probability.
- Marker cross-check: per predicted label, mean of canonical markers (CD3E, NKG7, MS4A1, JCHAIN, CD14, FCGR3A, C1QA, CD1C, LILRA4, TPSAB1); flags a label whose own marker is not its top marker.
- Round-91 interstitial non-responder programme (CD99, MAP4K4, CSF1R, PTPN6, TGFB1, HLA-DRA, IFI16, LILRB1, LILRB2, C1S) scored (scanpy score_genes) in query myeloid cells; LN vs NC vs AGN medians reported. Descriptive only: one LN patient, class II, cannot test the round-91 claim.
- IFN6 score (round 89 genes) in query myeloid and T cells, LN vs others, descriptive.

## Systematic open-data search (catalogue, no test)
- Query data/real/catalogues/geo_sample_index.tsv.gz and GEO E-utilities for 2024-2026 human LN kidney single-cell/single-nucleus series; record accession, n LN, labels (class, activity/chronicity, response), open matrix yes/no.
