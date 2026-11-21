# Q2 (Perez 2022 blood IFN vs nephritis beyond SLEDAI): blocked, 2026-11-21

No test was run, so no pre-registration was needed. The labels the test needs are not public.

What was checked:
- CELLxGENE dataset 218acb0f (1,263,676 cells; h5ad read remotely, obs only). Its obs fields are disease, disease_state
  (flare / managed / treated / na), Processing_Cohort, ancestry, sex, age and cell type. It has no nephritis,
  renal, SLEDAI, anti-dsDNA or complement field.
- GEO GSE174188 series matrix: pooled 10x wells only (no per-donor clinical characteristics).
- Zenodo 4724043 (authors' code): the clinical analysis reads `Clinical_Data.V2_updated.xlsx` (columns sledaiscore,
  acrrenalc, acrantidsdna, sledailowcomp, ...). That file is not deposited. Genotypes and phenotypes are in dbGaP
  phs002812 (controlled access).
- The PMC supplementary-tables zip returned an HTML challenge page instead of the file.

Published aggregate result (authors' notebook Figure_3_4 cell 87, not our test): membership of the
IFN-high "High" molecular group gave odds ratio 2.02 for ACR renal involvement (p = 0.14), against 2.68 for anti-dsDNA
(p = 0.02) and 5.50 for anti-Smith (p = 0.001). Not adjusted for SLEDAI; renal criterion is "ever", not current.

Unblock: dbGaP phs002812 authorised access (institutional signing official), or a public cohort with per-donor
nephritis labels and blood expression (candidate: GSE224705, queue Q11).
