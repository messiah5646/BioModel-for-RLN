# Round 112: carrying the early read-out forward works only because of voclosporin

*2026-11-23 · autonomous run · P2 = Q53 · branch rounds-70-72 · PREREG 5c9a740, data 143b72f, analysis 9ea993d*

**Decision this can change:** whether LUMEN's forecaster keeps "update with the drug's own early (≤26 wk) read-out" as its default at every horizon. If the early read-out still beats the null after adding the week-104 rows (H1) and the forecaster's shrunk update beats raw carry-forward (H2), keep it. If H1 fails, trust early read-outs only up to ~52 weeks and downgrade the round 85-87 claim.

## In one paragraph
Adding the published week-104 row for NOBILITY (TULIP-LN has no early visit, so it cannot enter) does not break the round-86 rule on its pre-registered score: carry-forward still beats "no effect" by a wide margin (IV-MSE 0.286 vs 0.592, ratio 0.48; 14 rows, 12 arms, 8 trials). The forecaster's actual update, which shrinks the early read-out toward the drug's prior, does better still (0.230). But the result is fragile. If each trial gets equal weight instead of each patient, carry-forward ties the null (ratio 1.04). The reason is plain once you look at the rows: nearly all of carry-forward's advantage comes from the two voclosporin trials (AURA-LV, AURORA), a calcineurin inhibitor that separates from placebo by week 24 and stays separated. Without them (exploratory, not pre-registered) carry-forward loses to the null (0.491 vs 0.355, 6 trials). For the slow, immune-upstream drugs, which are most of the pipeline, there is no evidence yet that a week-24 read-out predicts week 52-104.

## What was tested (pre-registered, PREREG-Q53.md)
Task C' from round 86, unchanged code, on resp_table_plus_pub.csv: for every arm with a CRR read-out at ≤26 weeks, predict the log odds ratio at each later visit (48-104 wk).
- C0 carry-forward (early value), C1 null (0), C2 shrunk to the pooled mean, F the LUMEN forecaster posterior, leave-one-trial-out.
- Score: inverse-variance-weighted MSE, as in round 86.

Honesty note: the week-104 counts were read in round 90, so this is a robustness check of a fixed rule, not a blind test. The table was rebuilt from the committed ClinicalTrials.gov snapshot and matches byte for byte.

## Results

| Score set | rows / arms / trials | C0 carry-forward | C1 null | C2 shrunk | F forecaster | C0/C1 |
|---|---|---|---|---|---|---|
| Original table (round 86) | 13 / 12 / 8 | 0.224 | 0.577 | 0.251 | 0.204 | 0.39 |
| **+ week-104 rows (primary)** | 14 / 12 / 8 | **0.286** | 0.592 | 0.274 | **0.230** | **0.48 PASS** |
| (a) one row per arm | 12 / 12 / 8 | 0.278 | 0.637 | 0.297 | 0.253 | 0.44 PASS |
| (b) arms with ≥10 patients | 13 / 11 / 8 | 0.235 | 0.572 | 0.244 | 0.200 | 0.41 PASS |
| (c) equal weight per trial | 14 / 12 / 8 | 0.671 | 0.645 | 0.648 | 0.635 | **1.04 FAIL** |
| Exploratory: no voclosporin trials | 11 / 9 / 6 | 0.491 | 0.355 | 0.349 | 0.339 | 1.38 FAIL |

- **H1 (primary): PASS.** Carry-forward beats null by 52% with the week-104 row added.
- **H1 robust: FAIL.** Variant (c) fails, so the pass rests on how patients are weighted.
- **H2: PASS.** The forecaster's shrunk update (0.230) beats raw carry-forward (0.286) and the null. Its edge over carry-forward comes from NOBILITY: dropping NOBILITY makes F 21% worse than C0. Shrinking toward the CD20 prior partly covers obinutuzumab's late separation (week 24 log OR -0.10, week 52 0.60, week 104 0.86).
- **H3 (descriptive, n = 2 long-horizon rows):** extra drift at 76-104 weeks is not detectable (noise-corrected variance -0.10 vs 0.118 at 48-52 wk). With CALIBRATE's wide error and only NOBILITY showing a real late rise, this says nothing either way.

Where the error comes from (share of C0's weighted error vs C1's): voclosporin rows, C0 0.044 vs C1 0.417. NOBILITY rows, C0 0.125 vs C1 0.097.

## Verdict and decision
- **Keep the forecaster's shrunk update (F) as the default.** Raw carry-forward is not the default. This matches the code as it stands (lumen_forecast.py already shrinks).
- **Downgrade:** "early read-outs fast-forward LN trials" goes from the round 85-87 wording to **Emerging (n = 2 voclosporin trials, 3 arms; fails under equal trial weighting)**. LUMEN should state that for slow, immune-upstream mechanisms (CD20, BAFF, CD40, IL-23, complement) a ≤26-week read-out has no demonstrated value over the prior. The DELAYED-target warning in lumen_forecast.py stays, and the preprint (P4) must say this plainly.
- No code change this round: a mechanism-specific rule would be a new model, and it needs its own pre-registration on independent data (follow-up Q58).

## Evidence labels
- Established: none.
- Supported: on IV weighting, the forecaster's shrunk early update beats carry-forward and the null at 48-104 weeks (n = 8 trials, 12 arms; LN only).
- Emerging: early read-out carry-forward beats null *in general* (holds only with voclosporin in the set; n = 2 trials carry it). Week-24 read-outs do not predict later effect for slow mechanisms (exploratory, n = 6 trials, not pre-registered).
- Hypothetical: drift grows with horizon for B-cell-targeted drugs (n = 1 informative arm, NOBILITY).

## Follow-up (passes the decision gate)
- **Q58:** pre-register Task C' on the non-renal SLE RCTs that P3 harvests (SRI-4/BICLA at week 12-24 vs week 52), split by fast vs slow onset, before scoring. Decision: whether LUMEN's early update is switched off by default for slow mechanisms. Run it inside P3, which builds the same table.

## Sources
- ClinicalTrials.gov results snapshot (data/real/ctgov/lupus_with_results.json.gz), 20 LN RCTs; rebuild script resp_table/build_resp_table.py (round 107).
- Furie R et al., Ann Rheum Dis 2022, NOBILITY week 104 (26/63 vs 14/62).
- Jayne D et al., Lupus Sci Med 2023, TULIP-LN week 104 (percentages, counts approximate).
- Code: experiments/lumen-bench/q53/run_q53.py; output q53/result.txt.
