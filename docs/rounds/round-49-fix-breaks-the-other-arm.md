# Round 49 — the fix works, and it breaks the other arm (2026-10-29)

Pre-registration: `docs/rounds/prereg_round49.md`, written before any fit.
Code: `experiments/rln-identifiability/r49_variants.py`.

**H1, H2 and H3 all supported — and the model is now falsified in a new place.** The
activity-independent progression term fixes the placebo arm exactly as round 48's
diagnosis predicted, narrows the identifiable interval, and **overshoots the treated arm.
No single variant fits both.**

Tests model dynamics on published aggregate arm-level data. Not an evaluation of any
therapy; no clinical or drug-development conclusion follows.

---

## 1. Held-out model selection

Three variants, **each exactly determined** (2 free placebo parameters on the same 2
year-1 targets), so each yields a point prediction and none gets extra slack. Held-out
target: placebo slope month 12→36, observed **−3.05** mL/min/yr.

| variant | structure | fit cost | predicted slope | held-out error |
|---|---|---|---|---|
| V1 | activity-driven damage only (round 48) | 2.2e−31 | −1.384 | **1.666** |
| **V2** | **activity-independent progression only** | 2.5e−32 | **−2.545** | **0.505** |
| V3 | both, fixed ratio | 2.2e−31 | −1.740 | 1.310 |

**H1 supported** — V2's error is 0.505, under the pre-declared 0.75 (deliberately tighter
than round 48's ±1.5, which was too loose to bite). **H2 supported** — V2 beats V1, so
round 48 §3's diagnosis was correct: damage that tracks current activity decelerates as
suppression takes hold, and the real arm does not.

## 2. H3 — a better structure *did* buy identifiability

Round 40 suggested structure might not help identifiability. Here it did.

| | admissible `h_drug` | predicted slope interval | span |
|---|---|---|---|
| round 48 (V1) | [0.0375, 0.0700] | [−1.354, +0.120] | **1.474** |
| round 49 (V2) | [0.0375, 0.0625] | [−2.504, −1.872] | **0.633** |

The span narrows **2.3×** and no longer changes sign. **H3 supported.** Worth recording
against round 40's expectation: replacing a wrong mechanism with a better one tightened
what the same data could determine, where adding observables had not.

## 3. The new falsification

The treated arm's interval is **[−2.504, −1.872]** against observed **−0.72**. It no
longer contains the observation — it misses by 1.15 to 1.78 mL/min/yr.

**So fixing one arm broke the other.** V1 under-predicted decline in the placebo arm; V2
over-predicts it in the treated arm. The data demands something neither has:

> progression accruing **independently enough of activity** to make the untreated arm
> decline near-linearly, yet **reducible by the drug** enough to leave the treated arm
> nearly flat.

## 4. A specific inconsistency, found by looking rather than guessing

Route 3 as implemented is **cosmetic**. The ODE's damage driver is
`k_tox · p0 · A · N`, with no haemodynamic factor, while the drug's proteinuria reduction
is applied post-hoc to the observation only:

```
model.py:85   (k_dam*A + k_tox*p0*A*N + k_prog) * N      <- no haemo term
r49:41        O[1:, 0] *= (1.0 - haemo)                   <- observation only
```

So in the model the drug lowers **what you measure** but not **what damages the kidney**.
Physiologically a haemodynamic reduction in protein leak should lower proteinuria-driven
tubular toxicity too. Repairing that would give the drug a protective route capable of
flattening the treated arm while progression still drives the untreated one — i.e. it
could satisfy both constraints in §3 without a new free parameter.

**That is the named next test. It is not run here, for the reason in §5.**

## 5. Holdout erosion — the honest cost of this round

These two held-out numbers have now been looked at **twice**: round 48 and round 49.

Round 49 was pre-registered and its variants were declared before fitting, so it is a
legitimate second look. But **a third test, designed after seeing §3's miss and tuned to
fix it, would not be a holdout at all** — it would be fitting with extra steps. Each look
spends holdout integrity, and this project's own rules make that explicit (never choose a
policy after seeing the comparison).

**Therefore: §4's repair must be scored against data not yet used.** The workaround is in
`docs/workarounds.md` — re-harvesting the per-visit series (13 visits per arm in year 1)
gives roughly an order of magnitude more held-out points from a source already confirmed
reachable, and converts this from an n = 2 test into a real one.

## 6. What this does not establish

- Held-out n = 1 per variant on the placebo slope. This ranks three variants and rules one
  out; it validates none.
- V3's mixing ratio (0.3) was fixed by assumption, not fitted, so V3 is one point in a
  family rather than the best mixture. A fitted mixture would likely beat all three and
  would also be under-determined.
- `k_tox = 0.2 · scale` throughout, fixed, not fitted.
- Aggregate arm-level means; nothing per-patient.
- §4's inconsistency is identified by reading the code, not by a test showing it matters.
  Its predicted consequence is a hypothesis.

## 7. Next, in order

1. **Fresh held-out data first** (`docs/workarounds.md` W1). Without it §4 cannot be tested
   honestly.
2. **Then** repair route 3 so the drug's protein-leak reduction enters the damage driver,
   and score it on the new holdout with a tolerance tight enough to fail.
3. A fitted activity-dependence fraction, reported as an interval rather than a point,
   once there is enough held-out data to support it.

*Research and simulation only. Not clinical or drug-development guidance.*
