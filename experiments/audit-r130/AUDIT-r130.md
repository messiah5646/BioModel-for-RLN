# Audit round 130 (rule 10): rounds 126-129

Scheduled run 2026-11-23 ~19:25 UTC. Clean worktree at b1f57fd (branch rounds-70-72).

Decision this can change: whether the RT-1 / RT-2 verdicts (and the two LUMEN changes they triggered: pooled-prior
CLI headline; relative effects only) stand, whether any label is downgraded, and whether the branch is ready to merge.

## Reproduction (committed md5-frozen inputs)
| Check | Result |
|---|---|
| Input md5 vs PREREG-RT1 / PREREG-RT2 (snapshot 1b31cd52, resp_table 3ce67bdb, sri4 b37dd7cd, builders 3828f6ed / ea42d9c7, stamp.py ca14be70) | all match |
| `build_resp_table.py` and `build_sle_table.py` rebuilt from the snapshot | resp_table.csv, sri4_table.csv, bicla_table.csv byte-identical |
| `run_rt1.py` -> rt1_result.json, rt1_predictions.csv | identical except run_utc / git_head |
| `run_rt2.py` -> rt2_result.json, rt2_predictions.csv, rt2_coding.csv | identical except run_utc / git_head |
| PREREG-RT1.md / PREREG-RT2.md edited after their commits? | no (git diff c8e07af..HEAD, d5b5f84..HEAD empty) |
| `lumen_forecast.py` md5 now 04957a31, frozen eae508f5 in PREREG-RT1 | differs: docstring + CLI only (rounds 127, 129); shipped functions (load, iv, dl, drift_var, forecast) unchanged; RT-1 reproduces exactly with the current file. Future frozen-md5 checks should expect this. |
| CLI smoke test (`belimumab BLyS --early 24 30 100 20 100`) | headline = pooled prior, `early_update_EMERGING` block, `control_rate: pooled, not modelled` - matches the r127 and r129 decisions |

## Timing (rule 12b)
- RT-1: PREREG c8e07af 15:29:02Z; run_utc 16:29:42Z (60.7 min); analysis commit fcaa984 16:30:21Z. Meets 12b.
- RT-2: PREREG d5b5f84 17:25:36Z; committed run_utc 18:25:37Z (60 min 01 s); analysis commit a27840d 18:26:07Z.
  The round-129 disclosure of a first execution at 18:25:25Z (11 s early) cannot be checked from git (no artefact of
  that run exists) and needs none: the code is deterministic on md5-frozen inputs and reproduces the committed output
  exactly, and 12 s leaves no room for outcome-informed edits. Verdict: verifiably pre-registered on commit timing,
  slip on record. No label change.

## Tests counted, rounds 126-129
- Rounds 126 and 128: PREREG only, no tests.
- RT-1 (r127): 1 primary + S1, S2 (3 strata), S3 (3 strata), S4, S5 = 10 scored comparisons, + 1 exploratory (drift variance).
- RT-2 (r129): 1 primary (combined) + LN guard + SLE stratum + S1 (2 covariates x 3 strata) + S3 + S4 + S5 + S6 = 13.
- Total 23 comparisons, 2 primary pre-registered tests; both primaries negative (FAIL; INCONCLUSIVE with LN guard FAIL).

## Labels (rule 10: n < 10 or not pre-registered -> downgrade)
| Claim | n | Pre-reg | Label | Audit |
|---|---|---|---|---|
| r127 shrunk update does not beat carry-forward, time-forward | 7 trials / 13 arms | yes (12b) | Emerging | correct (n < 10) |
| r127 pooled beats null time-forward | 50 trials (LN 15, SLE 35) | secondary of 12b test | Supported | correct; numbers 0.327/0.427, 0.132/0.190 reproduce |
| r127 S2 LN-only: null best | 4 trials | secondary | not labelled, stated as Emerging in preprint | correct |
| r129 meta-regression does not beat pooled | 60 trials | yes (12b) | Supported | correct |
| r129 LN FAIL 1.17x | 20 trials | yes (guard stratum) | Supported | correct |
| r129 pooled control-rate coverage 0.88 | 60 trials | secondary | Emerging | correct |
| r129 stringency-only LN -7% | 20 trials | exploratory reading of S1 | Hypothetical | correct |
| r129 early control read-out beats pooled | 16 trials | descriptive secondary | Emerging | correct |

## Verdict
PASS on code, PASS on labels. No downgrades. Documents brought up to date: where-this-stands.md through round 130,
preprint v0.1.5 (RT-2 rows in 3.1, abstract (1), 2.2, 5, 6). Branch rounds-70-72 is ready to merge (Thomas reviews).
