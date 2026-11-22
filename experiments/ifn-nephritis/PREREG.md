# PREREG round 99 (Q14 + Q29 design): does blood IFN6 flag nephritis beyond extra-renal activity and serology?

Committed 2026-11-22 BEFORE any expression matrix of these series was downloaded. Only GEO sample-characteristic
keys/levels (catalogue `geo_sample_index.tsv.gz`) were inspected to define variables. No IFN score has been computed.

## Question (Q2, blocked on Perez 2022, re-asked in open data)
In SLE patients, does the 6-gene blood interferon score (IFN6) separate patients with ACTIVE nephritis from those without,
beyond (a) extra-renal SLEDAI and (b) serology (anti-dsDNA, low complement)?

## Cohorts (open GEO, whole blood, Illumina HT-12 or Affymetrix)
- A = GSE65391 (pediatric SLE, Banchereau 2016). SLE only. ONE sample per patient: the lowest `visit` number.
- B = GSE72326 (adult SLE, Chiche/Chaussabel, group=SLE). ONE sample per patient: patient = title before the final `V<k>`,
  lowest k kept. GSE49454 is the same 157 SLE samples (identical field counts) and is NOT used (duplicate).
- C = GSE99967 (adult SLE, biopsy-proven active LN `ALN` vs SLE without LN `ANLN`). Secondary only (no renal items or serology).

## Variables
- y (active nephritis): A `renal`==1 (renal SLEDAI domain active); B `renal`==Y; C `nephritis state`==ALN.
- Extra-renal SLEDAI (erSLEDAI): A sledai - (proteinuria+hematuria+pyuria+urinary_casts) [items coded 0/4].
  B sledai - 4*(haematuria Y + leukocyturia Y + uriry casts Y + [24h proteinuria > 0.5]); floor at 0. C: not computable.
- lowC (low complement, binary): A `low_complement`>0; B `low c3`==Y or `low c4`==Y.
- dsDNA (binary): A `increased_dna_binding`>0; B `a-dsd` >= 10 (pre-specified cut; sensitivity: within-cohort z of log1p).
- IFN6 = mean of per-gene z-scores (z within cohort across the analysed SLE patients) of IFI27, IFI44L, IFIT1, ISG15, RSAD2,
  SIGLEC1; probe -> gene from the GPL annotation, highest-mean probe per gene (same genes and rule as rounds 89-98).
- Complete cases on (y, erSLEDAI, lowC, dsDNA, IFN6) for A and B.

## Models (logistic regression, L2 C=1, standardised features)
- T1 trivial serology baseline: lowC + dsDNA.
- T0 activity baseline: erSLEDAI.
- BASE: erSLEDAI + lowC + dsDNA.
- FULL: BASE + IFN6.
- IFN-alone: IFN6.

## Evaluation (held-out)
1. Leave-one-cohort-out (LOCO): train on A test on B, and train on B test on A (coefficients transported; IFN6 is z-scored within
   cohort, a stated limitation).
2. Within-cohort: 5-fold stratified CV repeated 20x (seed 0..19), out-of-fold AUC.
Uncertainty: 2000 patient-level bootstrap resamples of the held-out predictions (paired, FULL vs BASE).

## Pass rules
- H1 (primary): FULL - BASE AUC >= 0.05 AND paired-bootstrap 95% CI lower bound > 0, in BOTH LOCO directions. Otherwise FAIL.
- H2: IFN-alone AUC > T1 serology AUC in BOTH LOCO-test cohorts (point estimates), and IFN-alone AUC >= 0.60 in both.
- H3 (secondary, cohort C): IFN6 AUC for ALN vs ANLN >= 0.65 with bootstrap 95% CI lower bound > 0.5.
Prior (stated now): H1 expected to FAIL — IFN tracks global activity (round 97: rho 0.56 with SLEDAI) and serology,
and kidney IFN falls with chronicity (round 95). H3 uncertain.

## Reported regardless
Within-cohort CV AUCs for all models; raw AUC of erSLEDAI, IFN6; Spearman IFN6 vs erSLEDAI; n per cohort.
Research and simulation only. Not clinical guidance.
