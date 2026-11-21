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
- [done] Q3 (round 92: Geneformer DROPPED; LODO AUC IFN 0.68, Geneformer 0.77, PCA-50 0.82, QC-only 0.66; beats IFN in 7/11 donors, rule needed 8; commit RESULTHASH) Geneformer-V2 embeddings vs the 6-gene IFN score as features for pre/post B-cell-depletion classification
  (CD3xCD19 data); embeddings must beat the gene score to be kept (SBB rule).
- [todo] Q4 Research doc "SLE, LN and refractory LN: mechanistic reconstruction" (claude.ai artifact
  Mv7f8uVm8ELsioUZqncEBH): write D7, responsive vs refractory comparison, from primary literature.
- [done] Q5 (where-this-stands.md rewritten through round 90) Fold rounds 62-88 into project doc where-this-stands.md.
- [done] Q6 (round 90: frozen wk-104 forecasts for NOBILITY/TULIP-LN; PASS, MSE 0.055 vs carry-forward 0.242) Extend resp_table.csv with published later visits (NOBILITY wk76/104, BLISS-LN wk52 CRR, TULIP-LN
  wk104) and re-score lumen_forecast.py leave-one-trial-out; pre-register first.
- [done] Q7 (found GSE224705 blood responders n=448, GSE200306 serial biopsies; round 91 run on GSE200306) Search for patient-level refractory-LN switch data that is public (no access requests); catalogue only.
- [todo] Q8 Research doc: D2 LN map, then D1, D4, D5 ... in the doc's order.
- [weekly] W1 Check ClinicalTrials.gov for LN trials that newly posted results; for each, FREEZE a lumen_forecast
  prediction in experiments/lumen-bench/prospective/ BEFORE reading results, then score it. This is the real test.
- [weekly] W2 Check Virtual Cell Challenge 2026 winners and new single-cell lupus datasets on GEO/CELLxGENE.
- [todo] Q9 Add a flare early-warning test using urine markers (CD163/CD206) if any public longitudinal urine dataset exists (round 89 follow-up).
- [todo] Q10 Put IFN in the patient-state layer as a therapy-lowered set point; re-fit the CAR-T and CD3xCD19 kinetics jointly (round 89 follow-up).
- [todo] Q11 GSE224705 (448 blood samples, MMF/AZA responders vs non-responders, longitudinal): pre-register whether a blood IFN/plasmablast signature predicts MMF response beyond UPCR and SLEDAI.
- [todo] Q12 Validate the round-91 H3 interstitial myeloid/TGF-beta non-responder programme in an independent kidney dataset (AMP, KURDS-like, or GEO), pre-registered.
- [todo] Q13 (round 92 follow-up) Batch-aware embedding test: regress out or match on QC covariates (% mito, counts) and re-run IFN vs PCA vs Geneformer LODO; also test whether the IFN score alone survives QC adjustment. Pre-register.
- [todo] Q14 (Q2 follow-up) Find a public blood dataset with per-donor nephritis labels and SLEDAI (GSE224705 first, then GEO bulk SLE PBMC series with renal annotation) and run the Q2 question there. Pre-register.
