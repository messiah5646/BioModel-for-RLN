# Pre-registration: LUMEN benchmark on harvested LN trials (round 85), 2026-11-21

## Honesty note
The registry renal values were viewed while building fr_table.csv (unavoidable: extraction). Mitigations:
model forms, class labels and the delay constant are fixed here before any model is scored; nothing is tuned after.
Class labels are assigned by mechanism only:
- fast (acts on filtration barrier / complement effector / haemodynamics): voclosporin, ALXN2050 (factor D), ravulizumab (C5)
- slow (acts on lymphocytes / costimulation / plasma cells / interferon): abatacept, iscalimab (CD40), zetomipzomib, anifrolumab
Delay constants: tau_fast = 1 wk; tau_slow = 6.9 wk, derived from CD3xCD19 monocyte IFN kinetics (ifn-cart), not renal data.
Exclusions by rule: SELUNE, ORCHID-LN (no baseline UPCR in registry or abstract); ixazomib (baseline imbalance, n<10/arm);
REGENCY, BLISS-LN (responder rates only).
Unit: y = log(fraction of baseline UPCR remaining). Every metric is leave-one-trial-out (LOTO), each trial weighted equally.

## Task A: control arm of an unseen trial
- A0: constant = mean y of other trials' control rows (ignores time)
- A1: y = -k t (one parameter, fit on others)
- A2: FR = f_inf + (1 - f_inf) exp(-t / tau) (two parameters, fit on others)
Claim "control arm is predictable" if A2 or A1 beats A0 by >= 10% RMSE.

## Task B: add-on effect d = y_drug - y_ctrl (same trial, same week)
- B0: d = 0 (no effect)
- B1: d = pooled mean of other trials
- M1: d = class mean from other trials
- M2: d = E_c (1 - exp(-t / tau_c)); E_c fit on other trials of class c
Claim "mechanism predicts effect" if M1 or M2 beats both B0 and B1 by >= 10% RMSE.

## Task C: fast-forwarding a trial (early read -> late effect, within trial)
Arms with an early (<= 26 wk) and a late (>= 36 wk) visit. Predict late d from the earliest-to-26-wk d.
- C0: carry forward (d_late = d_early)
- C1: d_late = 0
- C2: delay-scaled: d_late = d_early * g(t_late) / g(t_early), g = 1 - exp(-t / tau_c)
Claim "early read-out can be fast-forwarded" if C2 beats C0 and C1 by >= 10% RMSE.

## Content criterion for a "possible fix"
At least one task gives a claim that survives; otherwise the fix must come from data the benchmark lacks,
and the round names which.
