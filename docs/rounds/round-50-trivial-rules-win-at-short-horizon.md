# Round 50 — on 16 genuinely unseen points, a two-point linear rule beats the model, and round 48's headline was horizon-specific (2026-10-29)

Pre-registration: `docs/rounds/prereg_round50.md`, written before any fit on this data.
Data: `data/real/aurora1_per_visit.csv` (newly harvested, now in-repo).
Code: `experiments/rln-identifiability/r50_trajectory_holdout.py`.

**H1 FALSIFIED. H2 FALSIFIED. H3 FALSIFIED for two of three variants.** The round-49
repair did not help, every mechanistic variant loses to both trivial comparators, and
reconciling this with round 48 shows **round 48's "beats the baselines" claim was
horizon-specific and should not have been stated unconditionally.**

Published aggregate arm-level data. Tests model dynamics, not any therapy.

---

## 1. The holdout erosion problem is now fixed

Round 49 §5 ruled that a third test against the same two numbers would be fitting, not
prediction. So the per-visit series was harvested: **12 visits × 2 arms × 2 endpoints**,
with SDs and N. Only week 52 had ever been used, so **weeks 30–48 are genuinely unseen**.

- **FIT:** weeks 2–24, both arms, both endpoints — 28 points
- **HELD OUT:** weeks 30, 36, 42, 48 — **16 points, never used in any round**
- Week 52 excluded from scoring because rounds 48–49 fitted to it

The data is saved in-repo. Round 35's copy was lost to the ephemeral container; that will
not happen again.

## 2. Result

All three variants carry **5 free parameters**, so the comparison is parameter-count fair.

| model | fit cost | held-out RMSE | UPCR | eGFR |
|---|---|---|---|---|
| V1 activity-driven damage | 82.45 | 2.3499 | 0.5453 | 3.2782 |
| **V2 activity-independent progression** | 80.88 | **1.7674** | 0.5588 | 2.4362 |
| V4 = V2 + route-3 repair | 82.67 | 1.9824 | 0.5548 | 2.7480 |
| **persistence** (hold week 24) | — | **1.3074** | 0.1275 | 1.8446 |
| **linear extrapolation** (wk 16→24 slope) | — | **1.1371** | **0.1194** | **1.6037** |

**H1 falsified: the repair made things worse.** V4 (1.98) is beaten by V2 (1.77). Round 49
§4 identified a genuine inconsistency — the drug's protein-leak reduction was cosmetic,
entering the observable but not the damage driver — and repairing it did **not** improve
held-out prediction. A correct criticism of the code does not imply a better model.

**H2 falsified, and this is the headline.** A **two-point linear extrapolation beats the
five-parameter mechanistic model by 1.55×** on genuinely held-out real data, and
persistence beats it too. On UPCR the gap is **4.7×** (0.119 vs 0.559).

This is the pattern this project keeps finding — a tissue lookup beating a cell model, a
2-parameter curve matching a 17-state mechanism — now reproduced on **real trial
trajectories** rather than synthetic data.

**H3 falsified for V1 and V4.** Observed week 24→48 eGFR change: treated −3.30, placebo
−1.70, both falling with the treated arm falling faster. V1 predicts both *rising* (+0.24,
+0.29) — the wrong sign. V4 gets the signs but the wrong ordering. **Only V2 is directionally
correct** (−0.78, −0.74), and even it understates the magnitude ~4×.

## 3. Reconciling with round 48 — the verdict depends on the horizon

Round 48 reported the model beating persistence (2.22) and linear extrapolation (3.26)
with 1.18–1.32. Round 50 reports the opposite ordering. Both are correct, and the
difference is the **prediction horizon**:

| test | horizon | winner |
|---|---|---|
| R50, weeks 30–48 (16 points) | ≤ 6 months | **trivial** — linear 1.14 vs best model 1.77 |
| R48, month 12→36 slope (2 points) | 2 years | **model** — 1.18–1.32 vs 2.22 and 3.26 |

The RMSEs are on different scales, so **only the ordering within each test is comparable.**
But the ordering reverses, and the mechanism is clear: over a few months the trajectory is
near-flat, so persistence is close to optimal and linear extrapolation is better still;
over two years linear extrapolation diverges while the model's saturating dynamics do not.

> **Correction to round 48.** Its claim that the model "beats both trivial comparators" was
> true only at the long horizon and only on a 2-point target. Stated unconditionally it was
> an overclaim. The defensible version: *this model has no short-horizon advantage over
> linear extrapolation, and whatever value it has is at horizons beyond a year.*

That scoped claim is more useful than the unconditional one, and it is testable: it
predicts a crossover, and the crossover point is a measurable quantity.

## 4. Two bugs, both caught by the pre-registered stopping rule

1. **The CSV reader** mis-parsed the commented header and silently produced wrong columns.
2. **Argument order.** `least_squares` calls `f(x, *args)` — x first — while the residual
   function was declared `(variant, x, weeks)`. Every evaluation threw and returned the
   penalty, so all three fits sat at cost 1.4×10⁷. **The first run produced a full,
   plausible-looking results table from entirely failed fits**, and the only thing that
   stopped it being reported was the pre-declared rule that a non-converged fit is not
   quotable.

That is now three rounds in a row where the stopping rule caught a bug that would
otherwise have produced publishable-looking numbers.

## 5. What this does not establish

- **16 held-out points from one trial with strong serial correlation**, so the effective
  sample size is well below 16 and no p-value is claimed.
- Arm-level aggregate means. Nothing per-patient.
- The horizon-crossover account in §3 is an explanation offered after seeing both results.
  **It is not separately tested.** The clean test is to score both model and trivial rules
  across a range of horizons on the same data and locate the crossover; that needs a
  longer per-visit series than this trial provides.
- Week 52's treated eGFR (−1.5) sits off the week 36–48 trend (−1.9, −2.8, −3.6). Flagged
  in the pre-registration before fitting. It is excluded here for holdout reasons, but it
  also means **rounds 48–49 fitted to a single visit that does not sit on its own trend
  line** — an additional reason to prefer trajectory fits.
- V4's failure does not vindicate the cosmetic implementation; the inconsistency is real
  and remains worth fixing for correctness, just not for predictive gain.

## 6. Next

1. **Locate the crossover horizon** — the one result here that generalises. Needs a trial
   with a longer per-visit series, or a second trial (15 more are available with results).
2. **Stop treating "beats a trivial rule" as a milestone to be announced** and start
   reporting it per horizon, as a curve.
3. V2 is the standing best variant. Its remaining defect is magnitude: it understates the
   observed week 24→48 decline roughly 4×.

*Research and simulation only. Not clinical or drug-development guidance.*
