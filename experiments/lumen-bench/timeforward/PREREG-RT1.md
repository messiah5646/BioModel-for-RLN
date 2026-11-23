# PREREG RT-20261123-1 (round 126): time-forward backtest of the shipped LN/SLE trial forecaster

Written 2026-11-23 (scheduled run) from `docs/PREREG-TEMPLATE.md`. Under rule 12b this run commits and pushes this file
and ENDS. Scoring runs in a later scheduled run, at least 1 hour after this commit. No download is needed: all inputs are
already committed (hashes below), so `fetch.sh` is not used; the result JSON is written through `scripts/provenance/stamp.py`.

## 1. Decision this can change (rule 7)
Whether LUMEN's default trial forecast includes the shrunk early read-out update (`lumen_forecast.forecast(early=...)`).
- If PASS (below) -> keep the early update as the default and the preprint may claim it, at the scored n.
- If FAIL, INCONCLUSIVE or UNDERPOWERED -> switch the early update off by default (pooled prior only; `--early` stays as an
  option labelled Emerging), and the preprint abstract claims the pooled prior only.
- Secondary decision (W1): if pooled-only 95% predictive intervals cover >= 90% of observed arms time-forward, W1 criterion
  (a) cannot tell models apart; PREREG-W1 scoring then reports criterion (b) as primary (amendment note, not a rule change
  to the frozen forecasts).

## 2. What has already been seen (disclosure)
Leave-one-trial-out results on these same tables are known: r86/r112 (LN carry-forward 0.286 vs null 0.592, rests on
voclosporin; without it null wins), r116 (SLE pooled beats anchor by 31%; SLE early carry-forward 0.031 vs null 0.245,
n 8 arms). Structural facts checked for this PREREG (no outcomes scored): every LN (20) and SLE (40) NCT has a
`resultsFirstPostDate` in the snapshot; arms with both an early and a late visit: LN 12 arms / 8 trials, SLE 8 arms / 4
trials (all-time, before burn-in and voclosporin exclusion). The analyst knows headline results of several trials.
The time-forward ORDER of scoring is new; the outcome numbers are not.

## 3. Data (frozen, committed)
- `data/real/ctgov/lupus_with_results.json.gz` -> `protocolSection.statusModule.resultsFirstPostDateStruct.date` per NCT.
- LN: `experiments/lumen-bench/resp_table.csv` (20 trials, CT.gov only; publication rows of resp_table_plus_pub.csv are
  NOT used because they have no CT.gov posting date of their own).
- SLE: `experiments/lumen-bench/sle/sri4_table.csv` (40 trials, SRI-4).
- Effect = log OR with 0.5 continuity correction, se as in `lumen_forecast.load`.

## 4. Visits
- LN late (target) per arm: latest visit with week <= 104 (as `lumen_forecast.forecast`). LN early: latest visit week <= 26.
  For the early-update subset the target is the latest visit with 48 <= week <= 104.
- SLE late per arm: latest visit week <= 52 (as P3 primary). SLE early: latest visit 12 <= week <= 26. Early-subset target:
  latest visit 48 <= week <= 52.

## 5. Time-forward scheme
- Trials are ordered by resultsFirstPostDate. For trial t, the training set = trials of the SAME stratum (LN or SLE) whose
  date is STRICTLY earlier (ties are not in each other's training set).
- Burn-in: trial t is scored only if its training set has >= 5 trials.
- Combined LN+SLE = union of the within-stratum forecasts (no cross-indication prior pooling).

## 6. Models (all fitted on the training set only)
- M0 null: log OR 0.
- M1 pooled-only: IV mean mu of training-set late arms; prior variance v0 = tau2_DL + se(mu)^2.
- M2 shipped default = M1 prior updated by the arm's own early read-out exactly as `lumen_forecast.forecast`:
  post = (mu/v0 + le/ve)/(1/v0 + 1/ve), ve = early se^2 + drift_var(training set). drift_var is computed from training-set
  arms with both early and late visits (same function); if fewer than 2 such training pairs exist the arm is NOT scored for
  M2 (it counts against the minimum n).
- C0 carry-forward (trivial baseline): the arm's own early log OR.
- M3 anchor (secondary): molecule > target > pooled, shrunk to mu with tau2 as `forecast(anchor=True)`.
- Strongest trivial baseline for the early subset (rule 2b): carry-forward C0 (the arm's own earlier value exists);
  null M0 and pooled M1 also count as baselines.

## 7. Metric
IV-weighted MSE on log OR, weights 1/se_obs^2 of the scored (target) arm, as rounds 86 and 116.

## 8. Primary analysis and minimum n
Early-update subset, combined LN+SLE, **voclosporin excluded** (molecule == voclosporin, i.e. AURA-LV and AURORA arms).
Minimum n: >= 6 scored arms from >= 4 trials. Below that: **UNDERPOWERED** (not a pass).

## 9. Pass / FAIL rules (rule 2b clauses included)
- **PASS**: IV-MSE(M2) <= 0.90 x IV-MSE(M1), AND M2 is not worse than C0 or M0 by >= 10%.
- **FAIL**: M2 loses to ANY trivial baseline (M0, M1, C0) by >= 10%, i.e. IV-MSE(M2) >= 1.10 x min(M0, M1, C0).
  Losing to a trivial baseline by >= 10% is FAIL, not inconclusive.
- **INCONCLUSIVE**: anything else. (Decision line treats it as "not shown".)

## 10. Secondaries (labelled; none can turn a FAIL into a pass)
- S1 same primary with voclosporin included; S2 LN-only and SLE-only early subsets (reported with n; each flagged
  underpowered under the same minimum n).
- S3 all-arm time-forward, strata LN, SLE, combined: M1 vs M0, M3 vs M1 (IV-MSE).
- S4 95% predictive-interval coverage of the observed log OR for M1 (interval mu +/- 1.96 sqrt(v0 + se_obs^2)) and M2
  (post +/- 1.96 sqrt(v_post + se_obs^2)), all-arm combined and early subset. Feeds the W1 decision in section 1.
- S5 leave-one-trial-out version of the same primary subset, for comparison with the time-forward numbers.

## 11. Labels (rule 12)
Supported only if PASS/FAIL rests on >= 10 trials in the primary subset; otherwise Emerging, with n trials and arms stated.
This test is the first under rule 12b, so if the analysis commit is >= 1 h after this commit it may be called verifiably
pre-registered (timing only; the analyst was not blind to published trial results).

## 12. Output
`experiments/lumen-bench/timeforward/run_rt1.py` -> `rt1_result.json` (via stamp.py, prereg path
`experiments/lumen-bench/timeforward/PREREG-RT1.md`) and `rt1_predictions.csv` (one row per scored arm and model).

Research and simulation only. Not clinical or drug-development guidance.

## Input md5 (frozen at PREREG)
- 1b31cd52ccc62dc6124cd611342ba5d0  data/real/ctgov/lupus_with_results.json.gz
- 3ce67bdb51c6597bce6666d914618778  experiments/lumen-bench/resp_table.csv
- b37dd7cd72c9e912ba079088532c3765  experiments/lumen-bench/sle/sri4_table.csv
- eae508f5fc57703d25699071ddce9bb3  experiments/lumen-bench/lumen_forecast.py
