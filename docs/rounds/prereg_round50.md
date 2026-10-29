# Pre-registration — Round 50: a 20-point holdout, and the route-3 repair
Written 2026-10-29 BEFORE any fit on this data. Seeds 0-11.

## Why this round exists
Rounds 48-49 scored against TWO held-out numbers and looked at them twice. Round 49
section 5 ruled that a third test tuned to fix its miss would not be a holdout.
`data/real/aurora1_per_visit.csv` is newly harvested: 12 visits x 2 arms x 2 endpoints.
Only week 52 was ever used before, so weeks 30-48 are genuinely unseen.

## Split, fixed now
- **FIT:** weeks 2-24, both arms, both endpoints = 28 points.
- **HELD OUT:** weeks 30, 36, 42, 48 = 16 points, never previously used in any round.
- Week 52 is reported separately and EXCLUDED from scoring, because rounds 48-49 used it.

## Models compared (same fit set, same free-parameter count where possible)
- **V1** activity-driven damage (round 48)
- **V2** activity-independent progression (round 49's winner on the old target)
- **V4 (new, the repair)** V2 plus route 3 entering the DAMAGE DRIVER, not just the
  observable. Round 49 section 4 found the drug's protein-leak reduction was applied
  post-hoc to the observation while the ODE's toxicity term carried no haemodynamic
  factor, so the drug lowered what was measured but not what damaged the kidney.

## Hypotheses
- **H1.** V4 beats both V1 and V2 on held-out RMSE over the 16 unseen points.
  *Falsified if either beats V4.*
- **H2.** V4's held-out RMSE beats both pre-declared trivial comparators:
  persistence (hold the week-24 value) and linear extrapolation of each arm's
  week 16->24 slope. *Falsified otherwise.* This is the bar that matters; this project
  has repeatedly found trivial rules tying mechanistic models.
- **H3.** V4 reproduces the SIGN of both arms' week 24->48 eGFR trend (treated falling,
  placebo falling more slowly). *Falsified if either sign is wrong.*

## Stopping rule
Any variant whose fit does not converge (cost not improving across 12 multistarts) is
reported as failed and its held-out numbers are not quoted.

## Declared in advance about the data
The treated arm's week-52 eGFR (-1.5) sits off the week 36-48 trend (-1.9, -2.8, -3.6).
That is noted now, before fitting, because rounds 48-49 fitted to week 52 alone. It is
NOT grounds for excluding week 52 from the figure; it is grounds for preferring a
trajectory fit over a single-endpoint fit, which is what this round does.

## What this cannot establish
Aggregate arm-level means; nothing per-patient. 16 held-out points from one trial with
strong serial correlation, so the effective sample size is well below 16 and no formal
p-value is claimed. Tests model dynamics, not any therapy.
