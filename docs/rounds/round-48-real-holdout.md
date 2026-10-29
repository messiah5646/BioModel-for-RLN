# Round 48 — a real out-of-sample test: beats the baselines, falsified on the placebo arm, and year-1 data cannot pin the long-term slope (2026-10-29)

Pre-registration: `docs/rounds/prereg_round48.md`, written before any fit was run.
Code: `experiments/rln-identifiability/r48_holdout.py`.
Data: published **aggregate, arm-level** endpoints from the EU Clinical Trials Register
(harvested in round 35).

**H1 supported. H2 FALSIFIED on the placebo arm. H3 supported — the model beats both
trivial comparators on held-out RMSE even in its worst case.** The most useful output is
§5: an identifiability result established against real trial data.

This is a test of model dynamics. It is **not** an evaluation of any therapy, and no
clinical or drug-development conclusion follows from it.

---

## 1. Why this is a holdout and round 40 was not

Round 40 compared the model against seven AURORA summary features and matched six; the
seventh — the month 12→36 eGFR slope — was the miss. **But those slopes were in the
comparison set, so that was a fit.** Here the model is fitted to **year-1 endpoints only**
and the month 12→36 slopes are never shown to the optimiser.

| | fitted | held out |
|---|---|---|
| placebo | ΔUPCR −1.88, ΔeGFR +1.5 at wk 52 | slope 12→36 = **−3.05** |
| treated | ΔUPCR −2.65, ΔeGFR −1.5 at wk 52 | slope 12→36 = **−0.72** |

## 2. A pre-registration error, and a bug the stopping rule caught

**The prereg's parameter accounting was wrong**, found during implementation and recorded
before any held-out number existed. It claimed "4 parameters on 4 targets, zero slack."
Parameters do not pool across arms:

- placebo: `{s_inf, k_dam}` — 2 free on 2 targets, **exactly determined**
- treated: `{cni_supp, haemo, h_drug}` — 3 free on 2 targets, **under-determined by one**

So the treated arm has no unique prediction. Rather than pick a point, the under-determined
direction is scanned and the prediction reported as an interval (§4).

**The pre-registered stopping rule then earned its place.** The first run's treated-arm fit
failed at every scanned value (cost ≈ 4.0, `cni_supp` pinned at its boundary), so under the
rule no prediction was reportable. The cause was a real bug: the route-3 haemodynamic
offset was applied from t = 0, scaling baseline *and* week 52 equally, so it could not move
the delta at all. Fixed by starting the offset at day 1, which is also the correct
physiology. **Had the rule not been pre-declared, the broken fit's numbers were reportable
and wrong.**

## 3. The placebo arm — exactly determined, and it MISSES

Year-1 fit is exact (residual cost 2.5×10⁻³², both targets hit to three decimals), so
there is no slack to hide in.

| | predicted | observed | error |
|---|---|---|---|
| placebo slope 12→36 | **−1.384** | **−3.05** | **1.666** |

Outside the pre-declared ±1.5 tolerance. **H2 is falsified on the placebo arm.**

**The direction of the miss is mechanistically informative.** The model under-predicts
progressive decline by more than half. Its damage accumulation decelerates, because
`dD/dt ∝ A` and activity `A` falls as background immunosuppression takes hold — whereas the
real arm declines close to linearly over two further years. That points to an irreversible
progression term that **does not track current activity**, which the model does not have.
A concrete, testable structural gap, and it falls straight out of the holdout.

## 4. The treated arm — consistent, but the prediction is weak

Admissible route-3 strength (those `h_drug` values admitting an exact year-1 fit):
**h_drug ∈ [0.0375, 0.0700]**, a factor of 1.87.

| | predicted | observed | |
|---|---|---|---|
| treated slope 12→36 | **[−1.354, +0.120]** | **−0.72** | **contains observed** |

So the model does not falsify here. But the interval is **1.474 mL/min/yr wide and changes
sign**, so it is a weak, non-committal prediction — and the whole interval sits inside the
pre-declared ±1.5 tolerance, meaning **that tolerance was too loose to test this arm at
all.** My own threshold failed to bite; recorded rather than reinterpreted.

**An overclaim avoided.** A coarse first scan converged at only one grid point, h_drug =
0.05, predicting **−0.782** against observed −0.72 — an apparent near-exact hit. Refining
the grid showed that point is one of thirteen admissible values spanning sign change.
**Reporting it would have been selecting a point from a continuum.**

## 5. The result worth keeping: an identifiability verdict on real data

Round 35 asked directly: *can route 2 (nephron loss) and route 3 (reversible
drug-haemodynamic) be separated from UPCR + eGFR time courses?*

> **No — not from year-1 two-endpoint data.** A 1.87× range of route-3 strength fits the
> year-1 endpoints *exactly and indistinguishably*, and across that range the predicted
> long-term filtration slope spans 1.474 mL/min/yr and changes sign.

This is the project's central thesis meeting real trial data: the year-1 readout does not
determine the long-term trajectory, and the width of that indeterminacy is now measured
rather than asserted. It also says what would fix it — the discriminating signal is the
*shape* of the early response (a step-then-plateau versus a sustained slope), so denser
sampling in the first weeks is worth more than a longer tail.

## 6. H3 — the model beats both trivial comparators

Held-out RMSE across the two slopes:

| | RMSE |
|---|---|
| **model** (treated = interval midpoint) | **1.180** |
| **model** (treated = worst case in interval) | **1.319** |
| persistence (slope 0 for both arms) | 2.216 |
| linear extrapolation of each arm's year-1 slope | 3.264 |

**Both comparators were pre-declared, and the model wins even taking its worst admissible
value.** Given this project's history — a tissue lookup beating a cell model, a 2-parameter
curve matching a 17-state mechanism — a clean win over trivial rules on genuinely held-out
real data is the strongest empirical contact the model has made.

Note what carries the win: mostly the *ordering* and the treated arm. The placebo arm
contributes most of the model's remaining error.

## 7. What this does not establish

- **Two held-out numbers.** This can falsify the model; it cannot validate it. No
  statistical confidence attaches to a win on n = 2.
- **Arm-level aggregate means, not patient-level.** Nothing here speaks to per-patient
  behaviour, which is where round 29's confusability result lives.
- The placebo arm is standard-of-care, not natural history, so its −3.05 is damage accrual
  on treatment.
- Fixed by assumption and not fitted: `k_tox = 0.2·k_dam`, `phi = 0`, and every humoral
  parameter. The holdout tests the kidney block, not the humoral layer.
- The two source studies report different week-52 ΔeGFR for the same arms (the extension
  is a subset), so levels were fitted from the randomised study and only the **slope** —
  less sensitive to that discrepancy — was predicted. Mixing the two levels would be an
  error.
- Round 47's retraction applies to the humoral headline, not here; these are different
  parameters.

## 8. Next

1. **Add an activity-independent progression term** and re-run this exact holdout. §3 names
   the gap and the test is now built, so this is cheap and decisive.
2. **Tighten the tolerance before re-testing** — ±1.5 was too loose for the treated arm.
   Pre-declare something that can fail.
3. Chase per-visit data for a non-CNI agent, which removes the route-3 confound entirely
   and would make the placebo-arm test cleaner.
4. The holdout harness now exists; every future model change should be scored through it.

*Research and simulation only. Not clinical or drug-development guidance.*
