# Audit round 155 (rule 10): rounds 151-154

Scheduled run 2026-11-27 ~10:20 UTC. Clean worktree at 99298b3 (branch rounds-70-72).

Decision this can change: whether the O1 (r151), O3 (r152) and O2 (r154) verdicts and their rule-8 settlements stand,
which labels enter the preprint, and whether rounds-70-72 is ready to merge. If a number does not reproduce, a PREREG or
run script was edited after commit, or the analysis ran < 1 h after its PREREG -> the verdict is withdrawn and the
rule-8 entry removed; if all check -> verdicts stand, preprint gets the three rows, branch ready.

## Reproduction (clean worktree, `python3 -I experiments/lumen-v3/run_o{1,3,2}.py`)
| Check | Result |
|---|---|
| `pytest tests` (17 lumen/ tests) | 17 passed |
| `python3 -m lumen.cli verify` | all 5 md5 true |
| r151 `run_o1.py` -> `o1_result.json`, `frozen_O1_2026-11-27.json` | identical except the stamp fields run_utc / git_head |
| r152 `run_o3.py` -> `o3_result.json`, `frozen_O3_2026-11-27.json` | identical except stamps |
| r154 `run_o2.py` -> `o2_result.json`, `frozen_O2_2026-11-27.json`, `o2_predictions.csv` | identical except stamps (csv md5 9fb6e189 unchanged) |
| Load-bearing numbers in round docs vs JSON | all match: O1a 0.836 / M1 0.857 / M0 1.232, diff -0.021 [-0.054, 0.004], LOTO +0.0005, LN -0.051 (n 15); O3 P1 0.1674 / M1 0.1606 / M0 0.2244 (x1.042), LN x1.11, b(ph3-ph2) LN +0.03 [-0.27, 0.32] n 20, SLE +0.10 [-0.02, 0.22] n 40, M1 cov95 ph3 0.75 (14 trials / 20 arms) vs ph2 0.926; O2 R1 0.923 / B1 0.858 / B2 0.876 / R0 0.835, LN 1.029 vs 0.883 (n 15), LOTO 0.913 vs 0.852 (n 60), D1 CF 0.584 vs R1-tf 1.031 (n 16) |
| PREREG-O1/O2/O3 and run_o1/o2/o3.py edited after first commit? | no (one commit each: 4e8975a, 7ae773f, 05328df; 815ea38, 220a413, fb692ae) |
| Run script committed before the result commit? | yes in all three (O1 33 s, O3 35 s, O2 69 s before) |

## Timing (rule 12b)
| Test | PREREG push | Result commit | Gap |
|---|---|---|---|
| O1 | 04:15:47Z | 06:25:47Z | 2 h 10 m |
| O3 | 04:26:01Z | 07:26:31Z | 3 h 00 m |
| O2 | 08:26:16Z | 09:27:20Z | 1 h 01 m |
All three are correctly labelled **12b-timed (committed data; earlier looks not excludable)**; none is verifiably
pre-registered (no outcome data fetched after the PREREG). O2's 60 control arms had already been analysed in RT-2
(r128-129); r154 discloses this.

## Tests counted, rounds 151-154
3 primary decision comparisons (O1a vs M1/M0; P1 vs M1/M0; R1 vs B1/B2), 1 pre-registered guard (O2 LN guard,
FAIL). About 20 sensitivities/secondaries (O1: LOTO, LN, SLE, O1b-d, cov80; O3: LOTO, LN, SLE, ph2/3 recode, RE,
by-phase, full fit; O2: LOTO, SLE, R0, R1-HN0.5, R1-w0.5, D1). 3 prospective freezes (descriptive). All three primaries
are negative (2 INCONCLUSIVE, 1 FAIL); with 3 primaries and none passing there is no multiplicity concern for a positive claim.

## Labels (two wording corrections, no n-based downgrades)
No Supported/Established claim in rounds 151-154 rests on n < 10. Two claims say more than the numbers do:
1. **r151** "Supported (n 50 trials): ... moves the predictive score by about 0.02 nats per trial — a direction
   consistent with the simulation". The primary CI crosses 0 and the LOTO difference is +0.0005 (CI -0.023 to 0.020),
   so no direction is shown. Corrected: **Supported (n 50 trials, 12b-timed): REML/RE/t did not beat M1 by the
   pre-registered margin (O1a - M1 -0.021, CI -0.054 to 0.004).** "Direction consistent with simulation" -> **Emerging**
   (time-forward only; LOTO shows none).
2. **r154** "R1 robust MAP worse than pooled for control-rate forecasting: Supported, n 50". R1 FAILS by the
   pre-registered point rule (+0.065 nats past a 0.05 margin), but the paired CI [-0.027, 0.146] includes 0, so "worse"
   is not shown. Corrected: **Supported (n 50 trials, 12b-timed): R1 does not beat the pooled control rate and fails
   rule 2b (+0.065 nats, CI -0.027 to 0.146); its 95% intervals over-cover (1.00).** The decision (stay
   relative-effects-only) is unchanged because it needed a PASS.
r152 labels are correct (Supported n 50 for "phase does not improve the forecast"; Emerging for direction and for the
phase-3 coverage of 0.75, n 14 trials).

## Observations (no action needed)
- Log scores are not comparable across O1 and O3: O1 averages arms within trial, then over trials (PREREG-O1 §7); O3's
  secondary log score is arm-averaged (`lumen.baselines.log_score`). M1 is 0.857 in O1 and 0.732 in O3 for that reason.
  Neither round compares them; noted so the preprint does not.
- The run scripts rewrite the committed `frozen_O*_2026-11-27.json` files (stamp fields only) when re-run. Harmless for
  the forecasts (numbers identical) but a re-run in the main tree would change the freeze timestamps; audits must re-run
  only in a scratch worktree (done here, worktree discarded).

## Verdict
**PASS on code; two label wordings corrected (r151, r154), no verdict changes.** O1 INCONCLUSIVE, O3 INCONCLUSIVE,
O2 FAIL stand, and so do their rule-8 settlements. Preprint v0.2.0 gets three §3.1 rows. Branch rounds-70-72 is ready to merge.
