# Audit round 135 (rule 10): rounds 131-134

Scheduled run 2026-11-24 ~15:25 UTC. Clean worktree at 1996a2a (branch rounds-70-72).

Decision this can change: whether RT-3 run 3 proceeds under PREREG-RT3 unchanged (and what its result may be called), whether
any label is downgraded, and whether the branch is ready to merge. If the frozen catalogue/PREREG were edited or outcome data
reached the repo before the PREREG -> amend labels before run 3; if clean -> run 3 proceeds as written.

## Reproduction
| Check | Result |
|---|---|
| r132: `run_rt1.py` in clean worktree -> rt1_result.json | identical except run_utc / git_head; M1 0.3114, M2 0.2586, C0 0.1461; S2 LN-novoc M1 0.6259 vs M0 0.2794 (the numbers the r132 note cites) |
| r132: CLI `lumen_forecast.py belimumab BLyS --early 24 30 100 20 100` | prints the r132 early_note ("NO early rule is validated") and `control_rate: pooled, not modelled` |
| r131: ALLEGORY NCT04963296 posted CT.gov results | yes, resultsFirstPostDate 2026-11-12 (API re-query 2026-11-24) |
| r133/r134: 10 input md5s frozen in PREREG-RT3 vs HEAD | all match |
| CATALOGUE-RT3.md, candidates_frozen.csv edited after 836f753? PREREG-RT3.md after 68900fb? | no (git diff empty) |
| Outcome files (sources/, new_arms.csv) committed before run 3? | no |
| Timing: catalogue 13:38:11Z -> PREREG 14:39:01Z (60.8 min) | meets rule 12b |

## Findings (both about what run 3 may be called; neither is an error in a reported result)
1. **N30 JASMINE early + late by-arm numbers were already committed** in round 131 (sweep_2026-11-23.md, b748458,
   2026-11-23 22:27Z): SRI-4 wk 24 and wk 52 percentages by arm. PREREG-RT3 section 2 says the analyst knows no early-visit
   by-arm counts; that is wrong for N30. N30 is a CT.gov-sourced candidate (no PMID), so run 3 may well include it.
2. **All 28 PubMed-sourced candidate abstracts (N01-N28) were committed before the PREREG** in raw/pubmed_efetch.xml
   (836f753; disclosed in CATALOGUE section 1 as downloaded, not parsed). Blind flag scan (abstract_flags.py, flags only):
   numeric content 28/28, an early-visit sentence with a number 5/28 (N01, N11, N12, N19, N24). PREREG-RT3 section 11 would
   let run 3 claim "verifiably pre-registered" from re-fetch times, which would mask this.
   Fix: AMENDMENT-A135-RT3.md (c181f31, committed before any run-3 fetch): N01-N28 and N30 are 12b-timed; "verifiably
   pre-registered" only if the decision is reached on the clean candidates (N29, N31-N36) alone with >= 10 trials (it cannot
   be: at most 7); reported sensitivity S7 without N30. Models, metric, minimum n and decision rules unchanged, so the
   amendment cannot move the decision.

## Tests counted, rounds 131-134
r131: descriptive sweep, 0 qualifying trials, 0 tests. r132: wording/CLI, no new analysis (numbers re-quoted from RT-1).
r133: wording + catalogue (no outcomes). r134: PREREG only. **0 new comparisons.**

## Labels
| Claim | n | Pre-reg | Label | Audit |
|---|---|---|---|---|
| r132: no early rule validated (pooled 0.311 vs update 0.259 vs CF 0.146) | 7 trials / 13 arms | 12b-timed (RT-1) | Emerging | correct |
| r132: LN-novoc null best 0.279 vs pooled 0.626 | 4 trials | secondary | Emerging | correct |
| r131: 0 qualifying early-before-late trials | sweep | descriptive | no label | correct; ALLEGORY claim reproduces |
| r133: RT-1/RT-2 = 12b-timed | - | - | wording | correct |
| r134: RT-3 "may be called verifiably pre-registered" | - | - | conditional | **amended**: cannot apply to N01-N28, N30 |

## Verdict
PASS on code (all reproduce), PASS on labels for reported results; one pre-run label amendment for RT-3 (c181f31).
No downgrades of existing claims. Branch rounds-70-72 is ready to merge (Thomas reviews). Run 3 may proceed under
PREREG-RT3 + AMENDMENT-A135-RT3.
