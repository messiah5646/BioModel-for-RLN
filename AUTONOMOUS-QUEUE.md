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

## Queue
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
- [todo] Q9 Add a flare early-warning test using urine markers (CD163/CD206) if any public longitudinal urine dataset exists (round 89 follow-up).
- [todo] Q10 Put IFN in the patient-state layer as a therapy-lowered set point; re-fit the CAR-T and CD3xCD19 kinetics jointly (round 89 follow-up).
- [todo] Q11 GSE224705 (448 blood samples, MMF/AZA responders vs non-responders, longitudinal): pre-register whether a blood IFN/plasmablast signature predicts MMF response beyond UPCR and SLEDAI.
- [todo] Q12 Validate the round-91 H3 interstitial myeloid/TGF-beta non-responder programme in an independent kidney dataset (AMP, KURDS-like, or GEO), pre-registered.
- [todo] Q13 (round 92 follow-up) Batch-aware embedding test: regress out or match on QC covariates (% mito, counts) and re-run IFN vs PCA vs Geneformer LODO; also test whether the IFN score alone survives QC adjustment. Pre-register.
- [todo] Q14 (Q2 follow-up) Find a public blood dataset with per-donor nephritis labels and SLEDAI (GSE224705 first, then GEO bulk SLE PBMC series with renal annotation) and run the Q2 question there. Pre-register.
- [todo] Q15 (round 93 follow-up) Exposure-gate test: find public LN data with month-1 drug exposure or B-cell counts (LUNAR/NOBILITY published tables, MPA PK cohorts) and pre-register whether exposure adds >=0.05 AUC over month-3 UPCR for 12-month CR.
- [todo] Q16 (round 93 follow-up) Damage-floor test: in AURORA biopsy sub-study or LNTN aggregates, check whether baseline chronicity >=4 caps proteinuria response (lower plateau, not slower rate) in the six-state model. Pre-register.
- [todo] Q17 (round 94 follow-up) Rebuild resp_table.csv (CRR counts per arm per visit) from data/real/ctgov/lupus_with_results.json.gz with a script, and commit it with git add -f so it is never lost again; confirm it reproduces result86.txt and the round-94 reconstructed arm table.
- [todo] Q18 (round 94 follow-up) Each W1 run: re-check the 13 trials in prospective/w1/w1_frozen_2026-11-21.json for posted results and score them under PREREG-W1; add early read-out updates only as separately frozen files, never edits.
- [todo] Q19 (round 95 follow-up) Pre-register and test in GSE200306 (or another public LN biopsy set with chronicity index): baseline kidney IFN score correlates NEGATIVELY with NIH chronicity index (AMP 2026 claim). Compare with a trivial baseline (chronicity vs proteinuria/eGFR alone).
- [todo] Q20 Research doc: D1 SLE map, then D4, D5, D6 ... in the doc's order (one deliverable per run).
- [todo] Q21 (round 95 follow-up) D2 gap: write a class V (membranous) LN sub-section from primary literature (PLA2R-negative antigens such as EXT1/EXT2, NCAM1; podocyte-facing mechanism).
