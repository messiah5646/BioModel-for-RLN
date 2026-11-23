# RESULT RT-20261123-2 (round 128): control-arm response model

PREREG: PREREG-RT2.md, commit d5b5f84 (2026-11-23T17:25:36Z). Code: run_rt2.py. Output: rt2_coding.csv, rt2_predictions.csv,
rt2_result.json (stamp: run_utc 18:25:37Z, git_head 5ae2eb7).

**Timing disclosure.** run_rt2.py was first executed at 18:25:25Z, 59 min 49 s after the PREREG commit (11 s short of the
rule-12b hour); it was re-run at 18:25:37Z (60 min 01 s) and that stamped output is the one committed. The code is deterministic
on md5-checked committed inputs, so both runs give identical numbers. No data were downloaded; the two table builders were
re-run at ~18:24Z and reproduce the committed CSVs byte-identically.

Units: one control arm per trial, LN 20 (latest wk <= 104), SLE 40 (latest wk <= 52). STRINGENT=1 (steroid/rescue wording):
LN 8/20, SLE 3/40.

| Set (LOTO) | n trials | B1 pooled IV | B2 unweighted | M1 year+stringent | M1 / best | Verdict |
|---|---|---|---|---|---|---|
| Combined (primary) | 60 | 0.227 | 0.223 | 0.236 | 1.06 | INCONCLUSIVE |
| LN (S2 guard) | 20 | 0.456 | 0.435 | 0.509 | 1.17 | FAIL |
| SLE | 40 | 0.166 | 0.166 | 0.163 | 0.99 | INCONCLUSIVE |

Secondaries: S1 year-only combined 1.05 (LN FAIL 1.15); stringent-only combined 0.96, LN 0.93 (inconclusive, exploratory).
S3 visit-matched (LN 15, SLE 17) FAIL 1.22 (LN 1.75). S5 time-forward (50 trials) FAIL 1.38 (LN 15: 2.92).
S4 95% PI coverage: B1 0.883 (LN 0.85, SLE 0.90), M1 0.917 -> pooled control-rate intervals are too narrow.
S6 (16 trials with an early control visit): carrying the early control rate forward IV-MSE 0.110 vs pooled 0.184; mean
logit change early->target -0.03 (control rates barely drift with visit length).

**Decision (pre-registered):** not shown -> LUMEN outputs relative effects (log OR) on a pooled control rate only, labelled
"control rate = pooled, not modelled"; absolute-rate claims stay out of the preprint. LN guard FAIL independently. Applied:
lumen_forecast.py CLI now prints control_rate = 'pooled, not modelled' (function unchanged; frozen forecasts reproduce).
Labels: meta-regression does not improve control-rate forecasts, Supported (n 60 trials: 20 LN, 40 SLE). Pooled control-rate
interval under-coverage: Emerging (single estimate 0.88 vs 0.90 target, n 60).
