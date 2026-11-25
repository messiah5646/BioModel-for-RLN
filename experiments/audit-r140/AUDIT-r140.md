# Audit round 140 (rule 10): rounds 136-139

Scheduled run 2026-11-25 ~13:25 UTC. Clean worktree at fe4e977 (branch rounds-70-72).

Decision this can change: whether RT-4 run 3 proceeds under PREREG-RT4 + AMENDMENT-R20261125-RT4 unchanged (and what its
result may be called), whether any label is downgraded, and whether the branch is ready to merge. If outcome data for
RT-4 reached the repo before the PREREG -> amend labels/sensitivities before run 3; if clean -> run 3 proceeds as written.

## Reproduction
| Check | Result |
|---|---|
| r136: `newtrials/run_rt3.py` in clean worktree -> rt3_result.json | identical except run_utc / git_head; primary 0 trials (UNDERPOWERED); S3 n 3 trials / 7 arms M1 0.3492 vs M0 0.3822 vs M3 0.3420 |
| r136: exclusion counts | 36 rows = 33 excluded (rule 3: 18, rule 5: 8, rule 4: 3, rule 2/4: 1, rules 1/2/6: 1 each) + 3 late-only '5(primary)' -> S3; matches queue line |
| r136: fetch times vs PREREG-RT3 (68900fb 14:39:01Z) | 46 sources fetched 16:24:17Z-16:30:43Z (> 1 h 45 min after); extraction 5abd3f3 before scoring c7c4954 |
| r137: W1 frozen trials posted results? | CT.gov API re-query 2026-11-25 13:2xZ: 0/13 hasResults; NCT05268289 TERMINATED (reproduces) |
| r138: 12 input md5s frozen in PREREG-RT4 vs HEAD | all OK |
| r138: MUSE/LOTUS filtered by EXC regex, not EXCLUDED set | reproduces: MUSE titles hit 'Corticosteroid'/'Interferon'; LOTUS titles hit 'Composite Response' |
| PREREG-RT4, AMENDMENT-R20261125-RT4, CATALOGUE.md edited after their commits? | no (one commit each: 67f2d26, 414a21e, 36dd117) |
| regdocs/sources/ or any Table A source fetched before PREREG/amendments? | no regdocs/sources; r137b MANIFEST-harvest.tsv has 0 Table A PMC/EMA/FDA URLs |
| r139 amendment before any fetch | yes (12:25:39Z, no sources dir) |

## Findings (about what RT-4 run 3 may be called; none is an error in a reported result)
1. **LOTUS early + late by-arm SRI-4 are in a committed PREREG input** (lupus_with_results.json.gz: "SRI-4 Composite
   Response at Week 24 / Week 52"). Never tabulated, titles only read (r138 and this audit). PREREG-RT4's "early values
   new" cannot be claimed for LOTUS; MUSE's committed CT.gov values are OCS-conditioned (fail the plain-SRI-4 rule anyway).
2. **Telitacicept ph2b source PMC10958275 full text was committed in r136**, before PREREG-RT4. Masked-digit scan: text
   holds no by-arm early SRI-4 value (only SLEDAI-reduction wording and an atacicept citation). Disclosure only.
   Fix for both: AMENDMENT-A140-RT4.md (03df5c2, before any run-3 fetch): S9 clean-early sensitivity without LOTUS/MUSE,
   label "early values new for 8 of 10 Table A trials", LOTUS CT.gov counts usable as population cross-check.
   Models, metric, minimum n and decision rules unchanged, so the amendment cannot move the decision.

## Tests counted, rounds 136-139
r136: RT-3 primary empty (0 decision comparisons) + S3 descriptive (3 models, n 3 trials). r137/137b: catalogue + harvest,
0 tests. r138: PREREG only. r139: amendment only. **0 new decision comparisons; 1 descriptive set.**

## Labels
| Claim | n | Pre-reg | Label | Audit |
|---|---|---|---|---|
| r136: RT-3 UNDERPOWERED (0 of >= 10 new trials qualify) | 36 candidates | PREREG-RT3 + A135 amendment | verdict | correct; reproduces |
| r136: S3 pooled 0.349 vs null 0.382 | 3 trials / 7 arms | secondary | descriptive, no label | correct |
| r137: 0/13 W1 frozen trials have results | 13 trials | registry fact | Established | correct; reproduces 2026-11-25 |
| r137: Table A source locations | 10 trials | - | Emerging | correct (sub-agent compiled, 2 spot-checked) |
| r138/r139: RT-4 label "12b-timed (late outcomes already seen; early values new)" | - | - | conditional | **amended**: early values not new for LOTUS (MUSE) |

## Verdict
PASS on code (r136 reproduces exactly), PASS on labels for reported results; one pre-run label amendment for RT-4
(03df5c2). No downgrades. Branch rounds-70-72 is ready to merge (Thomas reviews). RT-4 run 3 may proceed under
PREREG-RT4 + AMENDMENT-R20261125-RT4 + AMENDMENT-A140-RT4 (round 141, due >= 14:27Z, 1 h after this amendment, the
same convention as r139).
