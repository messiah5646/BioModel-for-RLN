# Autonomous research queue (LUMEN / lupus nephritis)

Each run: take the FIRST item with status `todo`, set it to `doing`, do it, set it to `done` with a one-line result
and the commit hash, and add any new follow-ups at the end. One item per run unless it finishes early.

## Standing rules (never break)
1. Pre-register (commit a PREREG file) BEFORE looking at the outcome data of any test. State the pass rule.
2. Every model is compared with trivial baselines on held-out data. Report failures as failures.
3. Label claims Established / Supported / Emerging / Hypothetical. Human data first.
4. Write each round to the claude.ai project as `claude/round-NN-<finding>.md`; push code to branch `rounds-70-72`.
5. Never contact anyone, create accounts, submit forms, or spend money. Never delete data.
6. Research and simulation only. No clinical advice.

## Open-data-only policy (Thomas, 2026-11-22)
Use only data downloadable without an account or a data-use agreement (GEO, CELLxGENE, ArrayExpress/BioStudies, ClinicalTrials.gov,
Zenodo/figshare, paper supplements). ARK/AMP Synapse is ALL controlled (every file carries the ARK controlled AR), ImmPort and dbGaP
need accounts: do not plan around them. When a task needs such data, look for an open substitute or mark it [dropped: not open].

## Queue
- [done] Q0 OPEN-DATA HARVEST (live session 2026-11-22: 597 GEO series -> 62,446-sample characteristics index + per-series clinical label flags; 34 expression matrices of the 28 clinically labelled series (0.49 GB); 6 CELLxGENE lupus datasets; ArrayExpress non-GEO list; ARK/AMP confirmed all controlled. scripts/harvest/README.md, data/real/catalogues/) Thomas: "stick to open data but get all of them".
- [done] Q1 GSE254176 flare test (round 89: all 4 tests fail; IFN is a patient set point, not a flare marker): pre-register whether the monocyte IFN score rises BEFORE a flare (pre-flare vs
  post-flare samples, 6 patients), then run it. Data: /home/claude/data/geo/GSE254176 (re-download from GEO if missing).
- [blocked: nephritis/SLEDAI labels not public (dbGaP phs002812 controlled); see experiments/perez-ln/BLOCKED.md] Q2 Perez 2022 (CELLxGENE 436154da): stream a subset; pre-register whether the blood IFN score separates SLE
  with vs without nephritis beyond SLEDAI; compare with a trivial baseline (anti-dsDNA/complement if available).
- [done] Q3 (round 92: Geneformer DROPPED; LODO AUC IFN 0.68, Geneformer 0.77, PCA-50 0.82, QC-only 0.66; beats IFN in 7/11 donors, rule needed 8; commit 508aaac) Geneformer-V2 embeddings vs the 6-gene IFN score as features for pre/post B-cell-depletion classification
  (CD3xCD19 data); embeddings must beat the gene score to be kept (SBB rule).
- [done] Q4 (round 93: D7 tab written; non-response = 4 separable failures: under-exposure, incomplete B-cell hit, persistent interstitial myeloid activity, damage floor; baseline barely separates; commit 8477717) Research doc "SLE, LN and refractory LN: mechanistic reconstruction" (claude.ai artifact
  Mv7f8uVm8ELsioUZqncEBH): write D7, responsive vs refractory comparison, from primary literature.
- [done] Q5 (where-this-stands.md rewritten through round 90) Fold rounds 62-88 into project doc where-this-stands.md.
- [done] Q6 (round 90: frozen wk-104 forecasts for NOBILITY/TULIP-LN; PASS, MSE 0.055 vs carry-forward 0.242) Extend resp_table.csv with published later visits (NOBILITY wk76/104, BLISS-LN wk52 CRR, TULIP-LN
  wk104) and re-score lumen_forecast.py leave-one-trial-out; pre-register first.
- [done] Q7 (found GSE224705 blood responders n=448, GSE200306 serial biopsies; round 91 run on GSE200306) Search for patient-level refractory-LN switch data that is public (no access requests); catalogue only.
- [done] Q8 (round 95: D2 LN map tab written; biopsy activity and chronicity are two separate axes in human data, matching the six-state split; kidney immune response shifts myeloid -> lymphoid and kidney IFN falls as chronicity rises; commit f5a6a29) Research doc: D2 LN map (D1 onward moved to Q20).
- [weekly] W1 (last done 2026-11-21, round 94: no LN trial newly posted results; 13 unreported LN RCTs given frozen prior-only forecasts in prospective/w1/, commit 631d3fc; nothing scorable yet) Check ClinicalTrials.gov for LN trials that newly posted results; for each, FREEZE a lumen_forecast
  prediction in experiments/lumen-bench/prospective/ BEFORE reading results, then score it. This is the real test.
- [weekly] W2 Check Virtual Cell Challenge 2026 winners and new single-cell lupus datasets on GEO/CELLxGENE.
- [blocked: no open patient-level longitudinal urine data with flare labels; best set AMP Phase II urine proteomics needs Synapse login+clickwrap; catalogue + frozen PREREG in experiments/urine-flare, commit 73d8dc0] Q9 Add a flare early-warning test using urine markers (CD163/CD206) if any public longitudinal urine dataset exists (round 89 follow-up).
- [done] Q10 (round 96: IFN_i(t) = own baseline x shared therapy multiplier 0.73, no fall by day 28; H1 set-point persistence FAIL (CD3xCD19 Spearman 0.49 < 0.5), H2 beats carry-forward PASS, H3 one multiplier for CAR-T and CD3xCD19 PASS; experiments/ifn-setpoint, commit 5e44266) Put IFN in the patient-state layer as a therapy-lowered set point; re-fit the CAR-T and CD3xCD19 kinetics jointly (round 89 follow-up).
- [done] Q11 (round 97: H1 MMF FAIL, delta AUC -0.025 [-0.08,+0.02], UPCR+SLEDAI 0.87 vs +IFN/PB 0.84; H2 AZA FAIL; UPCR alone 0.84 carries the signal; commit 2b2469f) GSE224705 (448 blood samples, MMF/AZA responders vs non-responders, longitudinal): pre-register whether a blood IFN/plasmablast signature predicts MMF response beyond UPCR and SLEDAI.
- [blocked: no public independent LN kidney set with response labels; AMP Ph1 (ImmPort SDY997) and Ph2 (Synapse) need accounts; GEO checked, see experiments/kidney-validation/CATALOGUE.md] Q12 Validate the round-91 H3 interstitial myeloid/TGF-beta non-responder programme in an independent kidney dataset (AMP, KURDS-like, or GEO), pre-registered.
- [done] Q13 (round 98: IFN survives QC adjustment, H1 PASS, resid 0.62, +0.05 over QC in 9/11; Geneformer still DROPPED; PCA edge survives linear regression (H3 'not technical' by rule) but shrinks 0.14->0.035 on QC-matched cells; commit 73c81f1) (round 92 follow-up) Batch-aware embedding test: regress out or match on QC covariates (% mito, counts) and re-run IFN vs PCA vs Geneformer LODO; also test whether the IFN score alone survives QC adjustment. Pre-register.
- [todo] Q14 (Q2 follow-up) Find a public blood dataset with per-donor nephritis labels and SLEDAI (GSE224705 first, then GEO bulk SLE PBMC series with renal annotation) and run the Q2 question there. Pre-register.
- [todo] Q15 (round 93 follow-up) Exposure-gate test: find public LN data with month-1 drug exposure or B-cell counts (LUNAR/NOBILITY published tables, MPA PK cohorts) and pre-register whether exposure adds >=0.05 AUC over month-3 UPCR for 12-month CR.
- [todo] Q16 (round 93 follow-up) Damage-floor test: in AURORA biopsy sub-study or LNTN aggregates, check whether baseline chronicity >=4 caps proteinuria response (lower plateau, not slower rate) in the six-state model. Pre-register.
- [todo] Q17 (round 94 follow-up) Rebuild resp_table.csv (CRR counts per arm per visit) from data/real/ctgov/lupus_with_results.json.gz with a script, and commit it with git add -f so it is never lost again; confirm it reproduces result86.txt and the round-94 reconstructed arm table.
- [todo] Q18 (round 94 follow-up) Each W1 run: re-check the 13 trials in prospective/w1/w1_frozen_2026-11-21.json for posted results and score them under PREREG-W1; add early read-out updates only as separately frozen files, never edits.
- [todo] Q19 (round 95 follow-up) Pre-register and test in GSE200306 (or another public LN biopsy set with chronicity index): baseline kidney IFN score correlates NEGATIVELY with NIH chronicity index (AMP 2026 claim). Compare with a trivial baseline (chronicity vs proteinuria/eGFR alone).
- [todo] Q20 Research doc: D1 SLE map, then D4, D5, D6 ... in the doc's order (one deliverable per run).
- [todo] Q21 (round 95 follow-up) D2 gap: write a class V (membranous) LN sub-section from primary literature (PLA2R-negative antigens such as EXT1/EXT2, NCAM1; podocyte-facing mechanism).
- [todo] Q22 (round 96 follow-up) Out-of-sample test of the frozen set-point multiplier (setpoint.py, F_LATE 0.73, flat to day 28): find a third public B-cell-depletion lupus blood dataset (rituximab/obinutuzumab bulk or scRNA, another CAR-T cohort) and pre-register that baseline x 0.73 beats carry-forward and the population mean at the late visit.
- [todo] Q23 (round 96 follow-up) Identify the IFN lag and half-life: look for public lupus blood data with weekly or day-14/day-56 samples after B-cell depletion; until then the lag (28 d) and half-life are bounds, not estimates.
- [dropped 2026-11-22: Thomas chose open data only; ARK/AMP data are all controlled-access (DUC + institutional signing official), checked file by file] Q24 (Q9 follow-up) Run experiments/urine-flare/PREREG-AMP-urine.md on AMP Phase II urine proteomics (syn53124971) once the files are in Downloads/BioModel-data/amp-urine.
- [todo] Q25 (round 97 follow-up) GSE224705 longitudinal: does the within-patient IFN6 trajectory (first->last visit) differ between MMF responders and non-responders (KURDS claim: IFN falls only in responders)? Pre-register a trajectory-vs-carry-forward test; patient-level, labels are constant per patient.
- [todo] Q26 (round 97 follow-up) Exploratory AZA signal (higher baseline IFN6 with non-response, n 18) needs an independent AZA cohort; search GEO for azathioprine-treated SLE/LN blood with outcome labels; catalogue, then pre-register.
- [todo] Q27 (round 98 follow-up) The QC-matched result says the PCA/Geneformer edge over IFN is mostly technical. Pre-register a within-library test: a dataset where perturbed and control cells share one library (hashtag/multiplexed lupus PBMC stimulation, e.g. IFN-alpha in vitro on GEO), then rerun IFN vs PCA vs Geneformer.
- [todo] Q28 (round 98 follow-up) Donor P12 inverts for PCA and Geneformer under every adjustment but not for IFN: check P12's sample metadata (dose, chemistry, cell counts, cell-type mix) descriptively; no pass rule.
- [todo] Q29 (Q0 harvest; answers Q2/Q14) Blood IFN vs nephritis beyond activity in OPEN cohorts with per-sample renal labels: GSE65391 (pediatric, 204 pts, SLEDAI, renal items, UPCR, C3/C4), GSE72326 (24h proteinuria, renal flags), GSE99967 (SLE with/without LN, biopsy class), GSE49454. Pre-register: IFN6 + nephritis status, adjusted for SLEDAI-without-renal-items; trivial baseline = anti-dsDNA + C3/C4. Leave-one-cohort-out.
- [todo] Q30 (Q0 harvest; round 89 follow-up) Flare test in a large open cohort: BOLD GSE92776 (steroid withdrawal, 153 pts, longitudinal, flare outcome) and GSE121239 (65 pts, longitudinal SLEDAI). Pre-register whether baseline/pre-flare IFN6 or plasmablast score predicts flare beyond baseline SLEDAI; replicate-or-refute round 89 (n 6) with ~200 patients.
- [todo] Q31 (Q0 harvest; IFN set-point) Set-point persistence at scale: in GSE65391 (996 samples) and GSE121239, estimate the fraction of IFN6 variance explained by patient identity vs visit-level SLEDAI (round 89 said 54% in n 6). Pre-register the 54% +/- 15 band.
- [todo] Q32 (Q0 harvest; trial simulator) ILLUMINATE-1/2 (GSE88884-6, tabalumab vs placebo, 1,760 pts, baseline + week 16/52 blood): pre-register whether the round-96 shared-multiplier set-point model predicts week-52 IFN6 better than carry-forward and the arm mean, out of sample (BAFF blockade = partial B-cell hit; frozen multiplier 0.73 should OVER-predict the fall).
- [todo] Q33 (Q0 harvest) Interventional pharmacodynamic sets for IFN-axis drugs: IFN-alpha kinoid (GSE39088, GSE72747/54/98 incl. LN with GFR), anti-IFN-gamma AMG 811 (GSE78193, GSE100430/1). Catalogue arms/timepoints, then pre-register a mechanism-specific IFN6 response test (type I blockade lowers IFN6, type II blockade does not).
- [todo] Q34 (Q0 harvest; Q12/Q19 substitute) Open LN kidney sets with histology: GSE32591/GSE32592 (glomerular + tubulointerstitial microarray), GSE200306, GSE42648 (repeat biopsies), E-MTAB-12257 (whole kidney), E-MTAB-13596 (kidney immune scRNA). Check which carry activity/chronicity or response labels; then run Q19 (kidney IFN vs chronicity) and Q12 (interstitial myeloid/TGF-beta) where labels exist.

