# Pre-registration — Round 48: a real out-of-sample temporal holdout
Written 2026-10-29 BEFORE any fit is run or any test number exists.

## Why this is not a repeat of round 40
Round 40 fitted the model against **seven** AURORA summary features and matched six.
The seventh — the month 12→36 eGFR slope — was the miss (predicted −2.07 / −3.38 against
observed −0.72 / −3.05). **That was a fit, not a holdout:** the slopes were in the
comparison set. This round fits on **year-1 data only** and predicts years 2–3, which
makes the slope genuinely held out.

## Data (published aggregate, arm-level; from round 35, EU Clinical Trials Register)
**Fitting set — AURORA 1, week 52 change from baseline (4 numbers):**
- placebo: ΔUPCR −1.88, ΔeGFR +1.5 ; baseline UPCR 3.87, eGFR 77.4
- treated: ΔUPCR −2.65, ΔeGFR −1.5 ; baseline UPCR 4.14, eGFR 78.3

**Held-out set — AURORA 2, eGFR slope month 12→36 (2 numbers):**
- placebo **−3.05** mL/min/yr, treated **−0.72** mL/min/yr

These two held-out numbers are NOT used in fitting, by construction (§Method).

## Model change required, declared now
Rounds 35/38 established the data demands a **third route**: a reversible
drug-haemodynamic eGFR term, distinct from monotone damage. `model.py` has no such term,
so one is added: `eGFR *= (1 - h_drug * u_drug(t))`, reversible, active only on drug.
Declared before fitting so it cannot be introduced later to rescue a miss.

## Free parameters: exactly 4, fitted to exactly 4 targets
`k_dam` (shared, disease), `cni_supp` (drug immunosuppression), `haemo` (drug UPCR
offset), `h_drug` (drug eGFR offset). Everything else fixed at reference values.
**Four parameters on four targets means the fit has zero slack, so the fit itself is
uninformative about model quality — the entire test is the held-out slope.**

## Hypotheses
- **H1 (ordering).** The model predicts the placebo arm declining faster than the treated
  arm over months 12→36. *Falsified if the ordering is wrong or equal.*
- **H2 (magnitude).** Both predicted slopes fall within **±1.5 mL/min/yr** of observed.
  *Falsified otherwise.* (Tolerance set now; round 40's miss was 1.35 off on the treated
  arm, so ±1.5 is a real but not trivial bar.)
- **H3 (baselines).** The model beats both trivial comparators on held-out RMSE across the
  two slopes:
  - **persistence:** ΔeGFR stays at its week-52 value → slope 0.0 for both arms
  - **linear extrapolation:** continue each arm's own baseline→week-52 slope

## Inclusion / stopping rules, declared before results
- Optimiser: least-squares on the 4 standardised residuals, multistart from 12 seeds,
  best final cost retained. If the best fit leaves any target residual > 5% of its
  target's scale, the fit is declared **failed** and no prediction is reported.
- `rtol=1e-9` for fitting (tighter than the population runs; this is a 2-arm problem).

## What this cannot establish
Arm-level aggregate means, not patient-level, so nothing about per-patient behaviour.
Two held-out numbers is a weak test in the statistical sense — it can falsify the model
but cannot validate it. The treated arm's identity is incidental: this is a test of the
model's dynamics, **not an evaluation of any therapy**, and no clinical or
drug-development conclusion follows from it either way.

*Research and simulation only. Not clinical or drug-development guidance.*
