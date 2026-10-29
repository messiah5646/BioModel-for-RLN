# Pre-registration — Round 49: fix the round-48 falsification
Written 2026-10-29 BEFORE any fit. Seeds 0-11.

## What is being fixed
Round 48's holdout FALSIFIED the model on the placebo arm: predicted -1.384 vs observed
-3.05 mL/min/yr (error 1.666, outside the pre-declared +-1.5). Diagnosis in R48 section 3:
the model's damage accumulation decelerates because dD/dt is proportional to activity A,
and A falls as background immunosuppression takes hold -- while the real arm declines
close to linearly for two further years.

## The declared model change
Add an ACTIVITY-INDEPENDENT progression term:   dD/dt = (k_dam*A + k_tox*UPCR + k_prog)*N
k_prog continues to accrue damage when activity is controlled. Declared before fitting.

## The test: held-out model selection, not a rescue
Three variants, EACH with exactly 2 free placebo parameters on the SAME 2 year-1 targets,
so each is exactly determined and gives a POINT prediction. No variant gets extra slack.

- **V1** {s_inf, k_dam}            activity-driven damage only (round 48's model)
- **V2** {s_inf, k_prog}           activity-independent progression only, k_dam = 0
- **V3** {s_inf, k_both}           both tied to one scale: k_dam = k_both, k_prog = 0.3*k_both*A_ref

## Hypotheses
- **H1.** V2's held-out error on the placebo slope is below **0.75** mL/min/yr. This is
  deliberately TIGHTER than round 48's +-1.5, which was too loose to bite. *Falsified if
  V2 is above 0.75.*
- **H2.** V2 beats V1 on held-out placebo error. *Falsified if V1 <= V2* -- which would
  mean the R48 section 3 diagnosis was wrong.
- **H3.** Under the best variant, the treated arm's admissible route-3 interval is
  NARROWER than round 48's 1.474 mL/min/yr span. *Falsified if it is not* -- this asks
  whether a better structural model buys identifiability, which R40 suggested it might not.

## Stopping rule, declared now
Any variant whose year-1 fit leaves a residual cost above 1e-10 is declared failed and its
prediction is not reported. Same rule caught a real bug in round 48.

## What this cannot establish
Held-out n = 1 per variant for the placebo slope. This can rank three variants and rule
one out; it cannot validate any of them. Aggregate arm-level published data. Tests model
dynamics, not any therapy; no clinical or drug-development conclusion follows.
