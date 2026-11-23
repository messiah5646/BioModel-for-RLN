# Audit round 115 (rule 10): rounds 111-114 and 116. Code PASSES; 2 label corrections

Run live 2026-11-23 in a clean `git worktree` of HEAD (c6e1cbd). Not independent: the auditor (Claude) wrote rounds 113,
114 and 116. The reproduction part is mechanical; the label part should be re-checked by the daily reviewer.

## A. Reproduction from committed code (clean checkout)
| Round | What was re-run | Result |
|---|---|---|
| 112 | `q53/run_q53.py` -> result.txt | identical |
| 113 | `dara-ifn/analyse.py` on GEO files (SHA256SUMS) -> final_R1, rho36 -0.30, median m36 1.103, R2, reliable rho -0.80 | identical |
| 113/Q57 | `kidney-compartments/regenerate_r91.py` (matrix) | identical |
| 114 | `check_A1_cd38.py`, `check_geneset.py` | identical |
| 116 | `build_sle_table.py` (both tables byte-identical), `run_p3.py` -> p3_result.json | identical |

## B. Tests counted
~24 score comparisons (r112: 4 models x 6 score sets), 7 (r113), 2 (r114), 5 (Q57), 13 (r116 incl. Q58 and robustness):
about 51. Few rest on p-values; none of the verdicts rests on a single p < 0.05.

## C. Label corrections (authoritative from this round)
| Round | Claim | Was | Now | Why |
|---|---|---|---|---|
| 112 | Forecaster's shrunk early update beats carry-forward and null at 48-104 wk | Supported (n 8 trials) | **Emerging** | n < 10 trials; the edge over carry-forward comes from one trial (NOBILITY), as the round itself reports |
| 114 | Assumption A1 (GEO patient numbers = paper SLE #) "now supported" | prose "supported" | **consistent (n 5, exact p 0.04)** | a 5-patient consistency check, not an evidence label |

Unchanged: r113 H1-fails (Emerging), r113 reliable subset (Emerging, not pre-registered), Q57 programme (Emerging),
r116 "anchor hierarchy does not beat the pooled mean in SLE" (Supported, n 40 trials, pre-registered; analyst not blind,
stated), r116 Q58 (Emerging, n 8 arms).

## D. Process findings
1. Round numbering collision (fixed): P3 was labelled "round 115" in its PREREG/data commits, which is an audit slot; the
   result is filed as round 116 and this audit takes 115.
2. r116 robustness "one arm per trial (highest dose)": for arms without an mg dose in the title (tabalumab q2w/q4w) the
   tie is broken by sort order, not by dose. Affects one robustness line only; verdict direction unchanged in every variant.
3. r113 PREREG trusted GEO sample titles; the TCR audit showed two are wrong. Lesson added for future single-cell rounds:
   verify sample identity (clonotypes or genotype) before demultiplexing on metadata.
4. r116 parser: an analyst edit briefly broke the placebo rule; caught by the pre-committed before/after row check. Keep
   that check for every table rebuild.

## Branch status
rounds-70-72 is a clean fast-forward of main; ready to merge as code. Label corrections live here and do not block it.
