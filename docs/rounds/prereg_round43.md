# Pre-registration — Round 43: the placebo falsifier
Written 2026-10-28, BEFORE any test number exists. Seed 20261028.

## Why
`lessons.md` §8: "Placebo response is the live falsifier. Real placebo arms in lupus
trials improve 33–44%. A model with a stable disease attractor cannot reproduce that,
and this is the project's clearest outstanding falsification."

## Sharpening the test
A **deterministic** model yields a response *rate* of either 0% or 100%. A response rate
of 33–44% is a statement about a **population**, so the real question is:

> Does a plausible between-patient parameter distribution produce a control-arm complete
> renal response (CRR) rate in the target band, **while simultaneously** keeping the
> treated arm's rate higher, and without absurd parameter spreads?

Note the model is **not** claimed to have a stable disease attractor under control
conditions: `supp(t)` encodes background immunosuppression (MMF/steroids), which is what
a lupus "placebo" arm actually receives. So the falsifier as written in `lessons.md` is
partly already answered; what is NOT answered is whether the *rate* is reproducible.

## Endpoint
CRR at week 52, AURORA-1-like: **UPCR ≤ 0.5** AND **eGFR ≥ 60 or no >20% fall from
baseline**. Both criteria required.

## Hypotheses and thresholds
- **H1 (primary).** A virtual population with between-patient variability on the
  top-sensitivity parameters produces control-arm CRR in **20–45%**.
  *Falsified if* no spread in the pre-declared grid reaches that band.
- **H2.** At the spread satisfying H1, the treated arm's CRR exceeds control by
  **≥ 10 percentage points** (AURORA-1 reported roughly 41% vs 23%).
  *Falsified if* the separation is < 10 pp or inverted.
- **H3.** The spread required for H1 is "plausible", pre-declared as **CV ≤ 50%** on
  jittered parameters. *Falsified if* CV > 50% is needed.

## Method, declared in advance
1. Sensitivity screen FIRST; jitter the top-ranked parameters **by sensitivity, not
   alphabetically** (lessons §7, the R24 error).
2. Log-normal jitter, 2000 patients per arm, seed 20261028.
3. Report the **trivial comparator**: a 2-parameter logistic fitted to the control-arm
   mean trajectory, asked whether it can hit the same CRR band. If it ties, that goes in
   the headline.
4. Report binomial CI on every rate.

## Inclusion rule, declared before seeing results
Exclude a simulated patient only if the integrator fails or eGFR hits the clamp
(non-physiological). Declared now so it cannot be chosen after seeing the table.

## What the round cannot establish
Synthetic data, illustrative parameters, one model structure. A matching CRR rate is
**not** evidence the mechanism is right — round 24 showed near-identical trajectories
from different hypotheses. This tests *reachability of a rate*, nothing more.
