# Audit round 160 (rule 10): rounds 156-159

Scheduled run 2026-11-27 ~16:25 UTC. Clean scratch worktree at 2b41e1c (branch rounds-70-72), discarded after use.

Decision this can change: whether CAL1 run 2 (RT-20261127-1) proceeds under PREREG-CAL1 unchanged, and whether
rounds-70-72 is ready to merge. If PREREG-CAL1 was edited after its commit, a cal1 script/result predates it, a freeze
md5 differs from its first commit, or a write-once re-run dirties git -> CAL1 run 2 is withheld and the branch is not
ready; if all check -> CAL1 run 2 proceeds unchanged and the branch is ready.

## Reproduction (clean worktree)
| Check | Result |
|---|---|
| `pytest tests` | 19 passed (17 + 2 write-once tests from r158) |
| `python3 -m lumen.cli verify` | all 5 md5 true |
| `run_o1/o3/o2.py` (python3 -I) | exit 0; outputs go to `.reruns/` (git-ignored); o1/o2/o3_result + frozen_O1/O2/O3 identical to committed minus stamps |
| `w1/freeze.py` (from its dir), `w1-sle/freeze_w1sle.py` (from experiments/lumen-bench, as its docstring says) | re-runs identical to committed minus stamps (`_meta.git_head`, `_meta.frozen_utc`) |
| `git status --porcelain` after all re-runs | empty (r158 write-once works) |
| FROZEN-MD5.tsv (6 freezes) vs current files | 6/6 match |
| Each freeze: commits touching it / md5 at first commit | 1 commit each; first-commit md5 = TSV for 6/6 (r158 claim "unchanged since first commit" reproduces) |
| r156 D-BICLA: NCT02437890 mBICLA wk24/wk48 in committed CT.gov snapshot vs bicla_table.csv | 10/10 rows match (28,28,22,24,22 / 28,32,22,19,22) |
| r157: correction boxes present in project docs r151, r152, r154 (project_search) | yes; r151 and r154 contain "does not beat"; "not over-optimistic" survives only in r152 file name + correction box |
| PREREG-CAL1: commits | 1 (6c2dd0a, 2026-11-27T15:26:04Z) |
| Any cal1/CAL1 file other than the PREREG, in any commit on any branch | none |
| PREREG-CAL1 names strongest trivial baseline (rule 2b) and a FAIL branch at >= 10% | yes: H (time-forward historical fraction), Brier(P1) >= 1.10 x Brier(H) or x Brier(C05), or abs gap >= 0.10 -> DROP; INCONCLUSIVE and UNDERPOWERED (< 20 trials) branches pre-registered |

Minor note (no action): `freeze_w1sle.py` fails with ModuleNotFoundError if run from the repo root; it must be run from
experiments/lumen-bench as its docstring says. Not a defect.

## Tests counted, rounds 156-159
0 decision comparisons. r156 = 1 data check (registry match); r157 = wording only; r158 = provenance infrastructure
(+2 unit tests); r159 = PREREG only, no outcome data read.

## Labels
No new Supported/Established claims in rounds 156-159. r156 uses "Verified" for the data check (corrected in r157).
r158's "unchanged since first commit" is a provenance fact and reproduces. Nothing to downgrade.

## Verdict
**PASS; nothing wrong found.** CAL1 run 2 proceeds under PREREG-CAL1 unchanged (as round 161; earliest 16:26Z,
already past). Branch rounds-70-72 is ready to merge (Thomas's step).
