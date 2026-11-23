# RESULT RT-20261123-1 (round 127): time-forward backtest — shrunk early update FAILS

PREREG c8e07af (2026-11-23T15:29:02Z); analysis run_utc 2026-11-23T16:29:42Z (60.7 min later; rule 12b timing met).
Inputs md5-identical to the PREREG. Script: run_rt1.py -> rt1_result.json, rt1_predictions.csv.

Primary (early subset, LN+SLE, voclosporin excluded, n 13 arms / 7 trials; IV-MSE on log OR):
null 0.367 | pooled-only 0.311 | shrunk update (shipped) 0.259 | carry-forward 0.146 | anchor 0.311
Verdict by rule: FAIL (update 0.259 >= 1.10 x carry-forward 0.146; it is 77% worse). Update beats pooled by 17% (not enough).
Label: Emerging (7 trials < 10).

Secondaries: S1 with voclosporin (16/9) FAIL, same pattern (0.276 vs CF 0.126). S2 LN-novoc (7/4) FAIL: null best 0.279,
CF 0.476, update 0.587. S2 SLE (6/3) UNDERPOWERED: CF 0.037 best. S3 all-arm: pooled beats null in LN (0.327 vs 0.427,
22 arms/15 trials), SLE (0.132 vs 0.190, 71/35) and combined (0.161 vs 0.224, 93/50); anchor ~ pooled (0.165).
S4 pooled 95% PI coverage 0.892 (n 93) < 0.90 -> W1 criterion (a) stays. S5 LOTO same subset (17/10): FAIL, update 0.191 vs CF 0.123.
Unscored for M2 (< 2 training drift pairs): CALIBRATE, NOBILITY, TULIP-1 x2.

Decision applied (pre-registered): early update OFF by default in the lumen_forecast CLI; headline = pooled prior;
--early prints the update under early_update_EMERGING; --use-early-update restores it. forecast() unchanged so frozen
forecasts and audits reproduce.

Exploratory (not pre-registered): drift_var re-estimated as 0 in 15 of 16 scored fits, so the shrinkage comes from a tight
pooled prior (small DL tau2), not from early->late drift. Raw carry-forward beats the update in SLE and in voclosporin
LN, but null beats carry-forward in LN without voclosporin (0.279 vs 0.476): no single early rule wins everywhere.
