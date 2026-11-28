# Audit round 170 (rule 10): rounds 166-169

Scheduled run 2026-11-28 ~14:25 UTC. Clean scratch worktree at 9168d4a (branch rounds-70-72), removed after use.

Decision this can change: whether the W3 first-batch freeze proceeds under PREREG-W3 (ec6a487) unchanged, whether any
r166-169 label is downgraded, and whether rounds-70-72 is ready to merge. If PREREG-W3 was edited, a W3 freeze/script
predates it, any of the 38 pooled trials has posted results, a freeze md5 changed, or REC1/POW1 do not reproduce ->
W3 freeze withheld and branch not ready; if all check -> W3 freeze proceeds unchanged, branch ready.

## Reproduction / checks
| Check | Result |
|---|---|
| `pytest tests` (clean worktree) | 21 passed |
| `python3 -m lumen.cli verify` | 5/5 md5 true |
| r166 `freeze_rec1.py --selftest` | PASS (w_hat 0.447 vs true 0.4) |
| r166 `freeze_rec1.py` real re-run | write_once refused overwrite, wrote .reruns/; rerun identical to committed frozen_REC1 except _meta.frozen_utc and _meta.git_head; committed md5 b59f83d9 = FROZEN-MD5.tsv |
| r166 timing | PREREG-REC1 08:32:25Z -> sources fetched 10:24:07Z+ (manifest) -> script 441ceac 10:25:53Z -> freeze 7f1e1af 10:26:21Z (1 h 54 min; fetch after PREREG) |
| r167 `pow1_sim.py` re-run | 49 s; pow1_result.json identical except run_utc/git_head (main-tree re-run would overwrite: R-20261128-5, still open) |
| r167 numbers | n 21: power to FAIL a gap 0.30 = 0.889 (MDD 0.30), gap 0.10 = 0.198; MDD n 50 0.20, n 100 0.15; frozen O1/O2/O3 8/5/8; REC1 P(any switch) max 0.0925 over competitor truths a/b/c and n 5-100 (claim "<= 10% even if a competitor is true" holds); KEEP-OPEN >= 0.88 at every n |
| r168 CLI | `lumen forecast ianalumab BAFFR --table sri4` prints OR 1.37 [0.98, 1.91], base rate 0.74 (80 arms/40 trials) labelled 'calibration not shown', calibration_note with +0.11 gap; no P_benefit |
| r168 Supported claim (H over-predicts, n 50) | r168 said "CI not computed for H". Computed here (experiments/audit-r170/h_gap_ci.py, same trial units and bootstrap scheme as run_cal1): gap +0.1136, 95% CI +0.005 to +0.231 (excludes 0). Label holds. |
| r168 propagation (rule 10 line) | project_search 'ranks trials better than the base rate on Brier': in the r161 doc only the correction box quotes the old phrase; the body reads 'has a lower Brier score'. Applied. |
| r169 `catalogue_w3.py` re-run | git status clean (w3_screen.csv, w3_candidates.csv byte-identical); 22 pass rules 1-7; manual_screen.csv excludes 5 -> 17 |
| r169 CT.gov completeness | live countTotal for the two committed queries = 59 and 200 = committed study counts (no truncated page) |
| r169 blinding | catalogue drops hasResults before reading other fields; only registry metadata committed |
| PREREG-W3 | 1 commit (ec6a487, 13:28:35Z); no freeze_w3.py, w3/frozen* or other W3 freeze in any commit on any branch |
| 38 pooled trials (21 W1/REC1 + 17 W3) | CT.gov v2 re-queried 14:3x UTC: 38/38 hasResults = false (2 COMPLETED and 1 TERMINATED among the 21, already known from the REC1 topline screen) |
| FROZEN-MD5.tsv (7 freezes) | 7/7 match; 1 commit each |
| main / deletions | origin/main last commit 2026-11-04; no file deletions in the window |

## Findings
1. **None that change a verdict or label.**
2. **Descriptive note (r167).** When no competitor is true (truth = P1), REC1's false SWITCH probability reaches 0.108 at n 10
   (0.098 at n 21). The round-167 sentence is about competitor truths and holds (max 0.0925); this adds the false-switch side.
3. **r168 CI filled (above).** The 'Supported (n 50)' over-prediction of H now has its CI (+0.005 to +0.231): as marginal as P1's.
4. Open reviewer minors R-20261128-3/-4/-5/-6 are confirmed still open (where-this-stands stale; labels on counts in r166/r167;
   pow1_sim bare write; F1 amendment). Not fixed here: they are their own queue items.

## Tests counted, rounds 166-169
0 decision comparisons on outcome data. r166 freeze (0 scored); r167 one simulation study (power tables; changed rule 8 re-open
conditions); r168 0 new tests (display/wording); r169 catalogue + PREREG (0 comparisons). This audit: 1 descriptive CI.

## Labels
One Supported claim in the window (r168, n 50 trials, 12b-timed): holds, CI now stated. Counts/procedures in r166/r167 carry
labels (R-20261128-4, open). No claim on n < 10. No downgrades.

## Verdict
**PASS on code and labels.** W3 first-batch freeze proceeds under PREREG-W3 unchanged (1-h clock passed 14:28Z; next run).
Branch rounds-70-72 is ready to merge (Thomas's step).
