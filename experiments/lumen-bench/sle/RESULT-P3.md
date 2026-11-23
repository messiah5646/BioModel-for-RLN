# P3 / Q58 RESULT (round 116): the LN anchor hierarchy does NOT generalise to non-renal SLE; early read-outs do carry forward

Numbering: the PREREG (64f5bd5) and data commit (34a85a9) say "round 115". Round numbers ending in 5 are audit rounds
(queue rule 10), so this result is round 116 and round 115 is the audit of rounds 111-114. Nothing else changed.

## Primary (pre-registered): leave-one-trial-out, SRI-4, arms with an anchor (23 trials, 45 arms)
| Model | IV-MSE | noise-corrected MSE |
|---|---|---|
| D0 null | 0.170 | 0.094 |
| D1 pooled mean of other trials | **0.082** | 0.050 |
| D3 anchored (molecule > target > pooled) | 0.108 | 0.075 |
D3 is 31% WORSE than D1 -> **FAIL** by the pre-registered rule.

| Robustness (pre-registered, not decisive) | D3 vs D1 |
|---|---|
| LOTO all arms (40 trials, 80 arms) | +19% worse |
| LOMO anchored (target-only anchors; 16 trials) | +16% worse |
| LOTO one arm per trial (highest dose; 23 trials) | +42% worse |
| BICLA LOTO anchored (14 trials) | +21% worse |
Every variant points the same way.

Why (descriptive): same-molecule trials disagree about as much as different molecules do (TULIP-1 anifrolumab log OR
-0.17 / -0.11 vs anifrolumab-Asia +1.19; SLE-BRAVE-I vs -II baricitinib 0.16-0.44 vs 0.03-0.06; ILLUMINATE-1 vs -2).
Between-trial tau^2 is small (~0.03) around a pooled log OR of ~0.31 (OR ~1.36), so the pooled mean is hard to beat and the
anchor (weight ~0.8 on a few trials) adds noise.

## Q58 (pre-registered): week 12-26 read-out -> week 48-52, 8 arms
C0 carry-forward 0.031, C2 shrunk 0.124, C1 null 0.245 -> **early read-out fast-forwards** (C0 87% better than null).
Fast vs slow split: fast n = 2 arms (TULIP-1), slow n = 6. Prediction "helps more for fast" **not supported** (C0 ties null
in fast; strongly beats null in slow). Too small to read as evidence either way; Emerging.

## Labels
- "The molecule -> target anchor hierarchy does not improve on the pooled mean for SLE SRI-4/BICLA effects":
  **Supported**, n = 40 trials (23 anchored), pre-registered, code committed. Caveat: analyst knew several published results.
- "Early (<= 26 wk) read-out carried forward predicts the 1-year effect": Emerging, n = 8 arms / 5 trials.

## Decision
- The hierarchy is an LN-specific heuristic. Combined with round 112 (LN signal rests on voclosporin's 2 trials), LUMEN's
  default prior for a new lupus drug becomes the **pooled mean** (all classes), not the molecule/target anchor. Anchoring is
  kept only as a reported sensitivity analysis.
- Keep carry-forward of an early read-out as the update step (supported here in SLE and in round 112 for LN, both small).
- For the preprint (P4): report P3 as a negative generalisation test; do not claim the anchor hierarchy as a method.
