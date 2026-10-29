# Pre-registration — Round 46
Written 2026-10-29 BEFORE any test number. Seeds 2000-2007 declared now.

Three tests, in increasing order of value.

## T1 — falsify my own explanation of the squeeze
R44 §3 explained the variability-power squeeze by claiming the jittered parameters enter
MULTIPLICATIVELY (`UPCR = p0*A*N`), so spread reshapes rather than translates the outcome
distribution. **That explanation was offered after seeing the result and was never tested.**

- **H1.** Under ADDITIVE jitter of matched variance, the squeeze is absent or materially
  weaker (slope of Cohen's d vs spread at least 2x shallower than multiplicative).
  *If additive jitter shows the same squeeze, my mechanism explanation is WRONG* and
  R44 §3 must be corrected — the phenomenon would survive, the explanation would not.

## T2 — is the viable window empty?
R43/R44 left this open and it is the decisive structural question. A window requires
control CRR in 20-45% AND treated-control >= 10pp SIMULTANEOUSLY.
- **H2.** There exists a CV in 0.50-0.75 where both hold with the seed-SD interval
  excluding the boundary. *Falsified if no CV satisfies both* — which would mean the model
  cannot be simultaneously realistic about placebo response and treatment effect, a
  structural defect, and that becomes the headline over everything else.
- Fine grid: CV 0.50, 0.55, 0.60, 0.65, 0.70, 0.75. 8 seeds.

## T3 — can enrichment break the squeeze? (the one that matters)
If heterogeneity erodes effect size, a trial might recover it by enrolling a narrower
population. **The test only counts if enrichment uses BASELINE-MEASURABLE quantities.**
Enriching on a latent rate constant is not actionable; enriching on baseline UPCR or eGFR is.

- **H3.** At least one enrichment rule defined solely on baseline UPCR and/or baseline eGFR
  raises Cohen's d by >= 25% relative to the unenriched population at the same CV.
  *Falsified if no observable-only rule helps*, which would be its own useful result:
  the squeeze would be irreducible by enrollment criteria.
- Pre-declared rules (fixed now, no post-hoc additions): baseline-UPCR tertiles and
  quartiles (top/middle/bottom), baseline eGFR above/below 60, and the four 2x2
  combinations of UPCR-median x eGFR-60.
- Reported alongside: required N per arm within each enriched subpopulation, since a
  larger effect in a smaller eligible pool is not automatically a better trial.

## Method
- Virtual-twin design for T3: the SAME patient is simulated under control and treatment,
  so the measured quantity is the true average treatment effect in that subpopulation.
  Real trials randomise; this measures a property of the model population, not a trial.
  Required-N figures still use the two-proportion formula on independent arms.
- N=4000 at CV=0.70 for T3, drawn once; enrichment rules evaluated post-hoc on that draw.
- Checkpoint per cell. rtol 1e-7 (verified 2e-9 vs tight solve).
- Report seed SD separately from binomial CI.

## Inclusion rule, declared before results
Exclude a simulated patient only on integrator failure or eGFR clamp.

## What none of this can establish
Synthetic, one model structure, one operating point, a stylised drug. An enrichment rule
that works here is a hypothesis about this model, NOT a trial-design recommendation.
