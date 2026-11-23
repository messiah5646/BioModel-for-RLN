# PREREG Q53 / P2 (round 112): does "carry the early read-out forward" survive the week-104 rows? (2026-11-23)

## Decision this can change
Whether LUMEN's forecaster keeps "update the prior with the drug's own early (<=26 wk) read-out, carried forward" as its
DEFAULT update at every horizon. If H1 and H2 pass -> keep it unchanged. If H1 passes but H2 fails -> the default
becomes raw carry-forward (the Bayesian shrink adds nothing). If H1 fails -> the early read-out is only trusted up to
~52 weeks; beyond that the forecaster reports the prior (molecule > target > pooled) and flags the early read-out
as uninformative, and the round-85-87 claim "early read-outs fast-forward trials" is downgraded to Emerging.

## Honesty note (what has already been seen)
The 3 added week-104 rows (NOBILITY 26/63 vs 14/62; TULIP-LN1 int 14/51, basic 8/45 vs 9/49) were read in round 90
and are visible in resp_table_plus_pub.csv. TULIP-LN1 has no visit <=26 wk in the table, so only the NOBILITY row can
enter Task C'. NOBILITY's week-24 log OR (-0.10) and week-104 log OR (~0.86) are therefore known before this test,
so H1 is NOT a blind test; it is a robustness check of a rule fixed in round 86. Scoring code is the round-86 code
unchanged (run86.py Task C'), applied to resp_table_plus_pub.csv. H2 and the robustness variants are fixed here,
before any number is computed on the extended table.

## Data
experiments/lumen-bench/resp_table_plus_pub.csv = resp_table.csv (rebuilt by build_resp_table.py from the committed
CT.gov snapshot data/real/ctgov/lupus_with_results.json.gz, round 107) + 3 published week-104 rows (add_publications.py).
Step 2 commit (data): rebuild both CSVs from the snapshot and confirm byte-identical content (or report the diff).

## Rows
Task C' rows as in run86.py: every arm with a CRR visit <= 26 wk (early = latest such visit) and every visit >= 48 wk
of that arm (one row per late visit; late visits <= 104 wk).

## Predictors (late log OR from early)
- C0 carry-forward: early log OR
- C1 null: 0
- C2 shrunk carry-forward: pooled IV mean + tau2/(tau2+se_early^2) * (early - pooled); pooled/tau2 from the primary-endpoint
  table (latest visit <= 104 per arm), exactly as run86.py
- F  LUMEN forecaster posterior, leave-one-TRIAL-out: lumen_forecast.forecast(molecule, target, early=<early counts>,
  R=<table without the held-out trial>, holdout=<trial>); drift_var also computed without the held-out trial
  (passing R without the trial does this). Prediction = post_logOR.

## Score
IV-weighted MSE (weights 1/se_late^2), as round 86. Also reported, not decisive: noise-corrected MSE (mean err^2 - se^2).

## Hypotheses and pass rules
- H1 (queue rule, primary): on the extended rows, C0 beats C1 by >= 10% IV-MSE (C0 <= 0.9 * C1).
- H1-robust (all three must also hold for H1 to count as "robust"):
  (a) one row per arm (latest late visit only);
  (b) arms with n_drug < 10 dropped (ALXN2050 fD_120);
  (c) trial-balanced weights (each trial's rows rescaled to equal total weight).
- H2: F IV-MSE <= C0 IV-MSE and F <= 0.9 * C1, on the extended rows.
- H3 (descriptive, no pass rule; n of long-horizon pairs is 2): drift variance (mean of (late-early)^2 - se_e^2 - se_l^2)
  for late visits 48-52 wk vs late visits 76-104 wk. If long-horizon drift >= 2x short-horizon drift, record as
  Emerging (n = 2 arms) and propose a horizon-scaled drift only via a new PREREG.

## Reporting
All numbers with the number of rows / arms / trials. Failures stated as failures. Same table also scored on the
original resp_table.csv to show exactly what the added rows change.
