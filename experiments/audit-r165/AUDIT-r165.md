# Audit round 165 (rule 10): rounds 161-164

Scheduled run 2026-11-28 ~09:25 UTC. Clean scratch worktree at 30e25f4 (branch rounds-70-72), removed after use.

Decision this can change: whether the REC1 freeze run proceeds under PREREG-REC1 (6f95104) unchanged, whether CURE2
proceeds on the CURE1 catalogue, and whether rounds-70-72 is ready to merge. If PREREG-REC1 was edited, a REC1
script/freeze predates it, any W1 trial has reported, a freeze md5 changed, or CAL1 does not reproduce -> REC1 freeze
withheld / CAL1 decision re-opened and the branch is not ready; if all check -> REC1 freeze proceeds unchanged,
CURE2 PREREG proceeds, branch ready.

## Reproduction / checks
| Check | Result |
|---|---|
| `pytest tests` (clean worktree) | 20 passed |
| `python3 -m lumen.cli verify` | 5/5 md5 true |
| `python3 -I experiments/lumen-v3/run_cal1.py` (real run) | exit 0; cal1_result.json identical to committed except run_utc/git_head/prereg_commit; cal1_units.csv byte-identical (write_once_text no-op); decision_a DROP, decision_b nothing added |
| `lumen forecast ianalumab BAFFR --table sri4` | OR 1.37 [0.98, 1.91], base_rate_positive 0.74 (80 arms, 40 trials), P_benefit not printed, CAL1 note present = round-161 doc |
| r161 timing | PREREG-CAL1 15:26:04Z -> run_cal1.py 16:30:14Z (1 h 04 min) -> result 16:31:36Z; 12b-timed holds |
| FROZEN-MD5.tsv (6 freezes) | 6/6 match; 1 commit each |
| r162 W1 status re-queried (CT.gov v2, 09:3x UTC) | 21/21 trials hasResults = false (r162 reproduces; REC1 window still open) |
| PREREG-REC1 | 1 commit (6f95104, 08:32:25Z); no freeze_rec1.py, prospective/rec1/ or other REC1 file in any commit on any branch |
| PREREG-REC1 rule 2b | strongest trivial baseline = B (base-rate print); FAIL at Brier >= 1.10 x B or gap >= B + 0.10; UNDERPOWERED < 5; SWITCH needs >= 10 and bootstrap bound |
| van Zwet constants in PREREG-REC1 s5b | match Table 1 of PMC9290572 (proportions 0.32/0.31/0.30/0.07; SNR SD 0.61/1.42/2.16/5.64; z SD 1.17/1.74/2.38/5.73 = sqrt(SD^2+1), consistent with the posterior formula) |
| CURE1 blinding | CATALOGUE.md and search_meta.json hold counts, follow-up durations from abstracts and header/row labels only; no per-patient values committed |
| CURE1 gate n | 40906188 n 11 (abstract regex '11','11'); 41612420 regex finds '12' and 'Nine' |

## Findings
1. **Minor (CURE2 must handle; does not change the gate).** For PMID 41612420 the abstract regex returns both '12'
   and 'Nine'; the predecessor 38777376 reported n 10. CURE1 took 12 from the table header level. If the extended
   follow-up set has fewer SLE patients with per-patient follow-up (e.g. 9), the gate still holds (11 + 9 = 20 >= 20)
   but the held-out set would be < 10 and CURE3 could not claim Supported. CURE2 PREREG must state: held-out n is
   counted at extraction; if < 10, the CURE3 label is capped at Emerging.
2. **Wording (r161).** "The pooled predictive ranks trials better than the base rate on Brier" — Brier is not a
   ranking metric. Corrected reading: "P1 has a lower Brier score than the base rate (0.210 vs 0.226, CI of the
   difference -0.034 to 0.004, crosses 0); the tertile reliability order is right". Label unchanged (Emerging).
3. **Fragility noted, no downgrade.** r161's Supported (n 50 trials) for "printed P_benefit over-confident" rests on
   the printed value (0.977 vs 0.688, gap CI 0.179-0.412) and on the pre-registered P1 gap (+0.110, CI 0.002-0.229,
   just past the 0.10 margin; LOTO 0.050). The printed-value part is robust; the P1 part is marginal, as r161 says.

## Tests counted, rounds 161-164
2 pre-registered comparisons (r161 (a) FAIL, (b) INCONCLUSIVE) + 8 descriptives; r162 status check (0 scorable);
r163 catalogue (0 outcomes); r164 PREREG only (3 primaries registered, 0 run).

## Labels
One Supported claim (r161, n 50 trials, 12b-timed) - n stated, pre-registered, holds. No claim on n < 10. No downgrades.

## Verdict
**PASS on code and labels; one wording fix, one CURE2 constraint.** REC1 freeze run proceeds under PREREG-REC1
unchanged (earliest 09:32Z, already past; next run). CURE2 PREREG proceeds with finding 1 added. Branch rounds-70-72 is
ready to merge (Thomas's step).
