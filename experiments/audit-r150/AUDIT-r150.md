# Audit round 150 (rule 10): rounds 146-149

Scheduled run 2026-11-27 ~05:25 UTC. Clean worktree at 5589f4d (branch rounds-70-72).

Decision this can change: whether the O1 and O3 analyses proceed under PREREG-O1 / PREREG-O3 unchanged, and whether
rounds-70-72 is ready to merge. If a PREREG was edited after push, a candidate was already scored on real trials, or a
pre-PREREG statement is wrong -> amendment (restarting the 1-h clock if it changes model/metric/pass rule/inclusion)
before the run; if all check -> O1 runs next, then O3, unchanged.

## Reproduction
| Check | Result |
|---|---|
| r148: `pytest tests` (17 tests: frozen inputs md5, r90 frozen, W1-SLE frozen, legacy CLI parity, `--early` no headline, gate rule 2b, contracts, candidates) | 17 passed (1.6 s) |
| r148: `python3 -m lumen.cli verify` | all 5 md5 true (lumen_forecast.py, resp_table, resp_table_plus_pub, sri4_table, bicla_table) |
| r148: `sim_o1.py` (synthetic, seed 0) | 9 result rows byte-identical to committed `sim_o1_output.txt` |
| r147: W1-SLE freeze | reproduced by `test_w1_sle_frozen` (shipped lumen_forecast.py, 8 trials) |
| r146: `catalogue_w1sle.py` on committed raw | `ctgov_screen.csv` and `candidates_w1sle.csv` identical; raw md5 688da1d7 = MANIFEST |
| r149: phase counts stated in PREREG-O3 §4 (no outcomes read) | reproduce: LN 20 = PHASE2 12 / PHASE3 7 / PHASE2-3 1; SLE 40 = PHASE2 26 / PHASE3 13 / PHASE1-2 1 |
| r149: rule-8 check "phase on the relative effect untested" | confirmed: no script in experiments/ uses phase as a covariate; project search finds only r148/r149 and the prior-art note |
| PREREG-O1, PREREG-O3, PREREG-W1-SLE, frozen_2026-11-26.json, lumen/candidates.py, run_rt1.py edited after commit? | no (one commit each: 4e8975a, 05328df, 5515355, a803568, 91876dc, fcaa984) |
| Any O1/O3 candidate run on real trials before its PREREG? | no: the only importer of `lumen.candidates` outside lumen/ and tests/ is `sim_o1.py` (synthetic); no result file in experiments/lumen-v3/ |
| O3 prospective freeze input (phase of the 8 W1-SLE trials) | present: `phase` column in `candidates_w1sle.csv` |

## Timing (rule 12b)
- PREREG-O1 4e8975a pushed 04:15:47Z; `lumen/candidates.py` (the O1a code) committed 26 s earlier in 91876dc. The 1-h
  clock runs from the PREREG; O1 analysis is due from 05:15:47Z. Both O1 and O3 are correctly labelled 12b-timed.
- PREREG-O3 05328df pushed 04:26:01Z; analysis due from 05:26:01Z.
- W1-SLE: catalogue -> PREREG gap 59 m 39 s, already disclosed in DISCLOSURE-run2-timing.md (catalogue is registry
  metadata only, so rule 12b's 1-h rule does not strictly apply to it).

## Tests counted, rounds 146-149
0 decision comparisons on real outcome data. 1 synthetic coverage simulation (9 cells x 4 models, descriptive).
8 prospective forecasts frozen (r147). 2 PREREGs written (O1, O3).

## Labels
No Supported/Established claims were added in rounds 146-149. The synthetic coverage numbers (M1 0.84-0.94 vs O1a
0.94-0.99 at tau >= 0.3) are described as synthetic in the RESEARCH note and PREREG-O1 §3; correct. Phase-2 optimism:
Hypothetical (r149); correct. **No downgrades.**

## Verdict
PASS on code and labels. O1 and O3 proceed under their PREREGs unchanged. Branch rounds-70-72 is ready to merge.
