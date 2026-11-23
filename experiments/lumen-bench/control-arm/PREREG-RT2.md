# PREREG RT-20261123-2 (round 128): control-arm response model for LN and SLE trials

Written 2026-11-23 (scheduled run) from `docs/PREREG-TEMPLATE.md`. Under rule 12b this run commits and pushes this file
and ENDS. Scoring runs in a later scheduled run, at least 1 hour after this commit. No download is needed: all inputs are
already committed (md5 below), so `fetch.sh` is not used; the result JSON is written through `scripts/provenance/stamp.py`.

## 1. Decision this can change (rule 7)
Whether LUMEN outputs ABSOLUTE arm response rates (control rate x relative effect) or relative effects only.
- PASS (below) -> LUMEN outputs absolute control and drug-arm rates using the meta-regression control model; the preprint
  may report absolute rates at the scored n.
- FAIL, INCONCLUSIVE or UNDERPOWERED -> LUMEN outputs relative effects (log OR) on a pooled control rate only, labelled
  "control rate = pooled, not modelled"; absolute-rate claims stay out of the preprint.
- LN guard: even on a combined PASS, LN absolute rates use the meta-regression only if the LN-only stratum (S2) is not a
  FAIL by the same rule; otherwise LN stays relative-only.

## 2. What has already been seen (disclosure)
The control-arm counts (r_ctrl/n_ctrl) entered every log OR in rounds 85-90, 112, 116 and 127, so the analyst has seen
them indirectly; control RATES have not been tabulated against posting year or definition stringency in any round.
resultsFirstPostDate was used as an ordering variable in round 127 (RT-1). The outcome-measure titles/descriptions have
NOT been read for this coding; the coding rule below is frozen before reading them. Structural facts checked for this
PREREG (no rates computed): LN table 20 trials, visits 24-104 wk; SLE SRI-4 table 40 trials, visits 4-52 wk.

## 3. Data (frozen, committed)
- `data/real/ctgov/lupus_with_results.json.gz`: resultsFirstPostDate per NCT; outcome-measure title + description.
- LN: `experiments/lumen-bench/resp_table.csv` (20 trials, CT.gov only; publication rows not used).
- SLE: `experiments/lumen-bench/sle/sri4_table.csv` (40 trials).
- Unit = one control arm per trial. Visit: LN latest week <= 104; SLE latest week <= 52 (the same targets as RT-1 / P3).
  If several rows of one trial at that visit carry different (r_ctrl, n_ctrl), take the row with the largest n_ctrl.
- Outcome y = logit control rate with 0.5 continuity correction: log((r+0.5)/(n-r+0.5)); se^2 = 1/(r+0.5) + 1/(n-r+0.5).

## 4. Covariates (exactly two, pre-named; frozen coding)
- YEAR: resultsFirstPostDate as decimal year, centred on the training-set mean.
- STRINGENT (binary): the outcome measure from which the control count was taken is located by re-running the committed
  selection (`resp_table/build_resp_table.py` for LN, `sle/build_sle_table.py` SRI path for SLE) and recording its title
  and description. STRINGENT = 1 if title+description matches, case-insensitive,
  `predniso|corticosteroid|glucocorticoid|steroid|rescue (medication|therap)`; else 0. Missing description -> title only.
  The coded table (nct, title, description[:300], STRINGENT) is written to `rt2_coding.csv` BEFORE any model is fitted
  and committed with the result. No hand edits to the coding.
- If STRINGENT is constant within a training set, the model drops it for that fold (year only) and the fold is flagged.

## 5. Models (fitted within stratum, on training trials only)
- B1 pooled IV mean of y (weights 1/se^2).  B2 unweighted mean of y.  (Trivial baselines.)
- M1 meta-regression: WLS of y on [1, YEAR, STRINGENT], weights 1/(se^2 + tau2), tau2 = DerSimonian-Laird tau2 of the
  intercept-only training fit.
- Strongest trivial baseline (rule 2b): min IV-MSE of B1, B2. Carry-forward is NOT applicable for the primary use case
  (control rate forecast before a trial reports anything, so no own earlier value exists); it is reported descriptively
  only (S6).

## 6. Held-out scheme and metric
Leave-one-trial-out within stratum (LN, SLE); combined LN+SLE = union of the within-stratum predictions (no cross-stratum
pooling). Primary metric: IV-weighted MSE of predicted vs observed logit control rate, weights 1/se_obs^2.

## 7. Minimum n
Combined >= 40 scored trials with LN >= 15 and SLE >= 25. Below that: UNDERPOWERED (not a pass).

## 8. Pass rule
PASS: IV-MSE(M1) <= 0.90 x min(IV-MSE(B1), IV-MSE(B2)) on the combined set.

## 9. FAIL rule (rule 2b, mandatory)
FAIL: IV-MSE(M1) >= 1.10 x min(IV-MSE(B1), IV-MSE(B2)). Losing to ANY trivial baseline by >= 10% is FAIL, not
inconclusive. Anything between = INCONCLUSIVE (decision treats it as "not shown").

## 10. Secondaries (labelled; none can turn a FAIL into a pass)
- S1 each covariate alone (YEAR only; STRINGENT only).
- S2 LN-only and SLE-only verdicts by the same pass/FAIL rule (LN-only drives the LN guard in section 1).
- S3 visit-matched: LN restricted to 48 <= week <= 56, SLE to week 52 (removes visit-length confounding).
- S4 95% predictive-interval coverage of observed y for B1 (mu +/- 1.96 sqrt(tau2 + se_mu^2 + se_obs^2)) and M1;
  if B1 coverage < 0.90, LUMEN must widen control-rate intervals regardless of the verdict.
- S5 time-forward version (train on strictly earlier-posted trials of the same stratum, burn-in >= 5).
- S6 descriptive: for trials with a control visit <= 26 wk and a later target, error of carrying the early control rate
  forward (shows how much control rates drift with visit length).

## 11. Labels (rule 12)
Supported only if the verdict rests on >= 10 trials per stratum used; n trials per stratum stated next to the label.
Verifiably pre-registered (timing only) if the analysis commit is >= 1 h after this commit; the analyst is not blind to
published trial results.

## 12. Output
`experiments/lumen-bench/control-arm/run_rt2.py` -> `rt2_coding.csv`, `rt2_result.json` (via stamp.py, prereg path this
file), `rt2_predictions.csv` (one row per trial x model).

Research and simulation only. Not clinical or drug-development guidance.

## Input md5 (frozen at PREREG)
- 1b31cd52ccc62dc6124cd611342ba5d0  data/real/ctgov/lupus_with_results.json.gz
- 3ce67bdb51c6597bce6666d914618778  experiments/lumen-bench/resp_table.csv
- b37dd7cd72c9e912ba079088532c3765  experiments/lumen-bench/sle/sri4_table.csv
- 3828f6ed31154f404c69db84be887d5f  experiments/lumen-bench/resp_table/build_resp_table.py
- ea42d9c7d6b37934fb5dc8ef74eb310d  experiments/lumen-bench/sle/build_sle_table.py
- ca14be70b213b98bbc2b830e2b092329  scripts/provenance/stamp.py
