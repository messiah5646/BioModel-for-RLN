# Round 53 — at 15–27 months ahead the trivial rules still win, and round 48's long-horizon result is falsified (2026-10-29)

Computed. Data: `data/real/aurora2_per_visit.csv` (newly harvested).
Code: `experiments/rln-identifiability/r53_long_horizon.py`. Reuses round 50's fits —
nothing re-fitted.

**The model has no predictive advantage over trivial rules at any horizon from 6 weeks to
27 months. Round 48's long-horizon win does not survive being scored on levels instead of
a slope.**

---

## 1. What this round could do that round 52 could not

Round 52 was capped at 28 weeks ahead and found the model losing everywhere, leaving the
crossover hypothesis alive beyond that window. The continuation study extends to month 36.
**Months 18, 24 and 30 were never used in any round** — rounds 48/49 used only the month
12→36 slope, i.e. the two endpoints.

Fit: AURORA 1 weeks 2–24. Held out: **12 points at 15–27 months ahead.**

## 2. Result

Error in units of each endpoint's own SD. Offset-corrected scores anchor the model at
month 12 exactly as the trivial rules are anchored, which isolates trajectory **shape**
from the level discrepancy in §4.

| method | raw | offset-corrected |
|---|---|---|
| V1 | 0.300 | 0.276 |
| V2 | 0.296 | **0.240** |
| V4 | 0.288 | 0.239 |
| **persistence** | **0.163** | **0.163** |
| linear extrapolation | 0.200 | 0.200 |

**Best model 0.239 against best trivial 0.163 — ratio 1.47. Trivial still wins.**

## 3. The gap narrows with horizon and never crosses

| horizon ahead | best model / best trivial |
|---|---|
| 6 weeks (R52) | 3.00 |
| 18 weeks (R52) | 1.87 |
| 28 weeks (R52) | 2.96 |
| **15–27 months (this round)** | **1.47** |

So round 50's horizon account was **directionally right** — the gap does close, 3.00 → 1.47
— but there is **no crossover within 27 months**, which is the entire span of available
data. The honest conclusion:

> **This model has no demonstrated predictive advantage over a persistence rule at any
> horizon that exists in the data.** The crossover is no longer an open hypothesis waiting
> on longer data; the longer data is now in hand and the crossover is not there.

**A reversal worth noting:** persistence (0.163) beats linear extrapolation (0.200) at long
horizon, the opposite of the short-horizon ordering where linear won. Over two years the
trajectory is near-flat, so holding the last value wins and extrapolating a slope diverges.

## 4. Round 48 is falsified, and the reason is instructive

Round 48 reported the model beating persistence (2.216) and linear extrapolation (3.264)
with 1.18–1.32 on the month 12→36 eGFR slope. That was this project's strongest empirical
claim.

**It does not survive.** Scoring the same models on the actual per-visit **levels** at
months 18/24/30 reverses the verdict. Both results are arithmetically correct:

- the model gets the **slope** roughly right (−2.5 against −3.05 observed)
- while being **wrong about the levels** at every intermediate visit

A slope is a two-point summary that discards the trajectory. Scoring it rewarded a model
that has the right overall tilt and the wrong position. **Round 48's choice of a
slope-based target flattered the model, and no level-based metric at any horizon favours
it.**

That makes three corrections now attached to round 48: its parameter accounting (round 49
§2), its unconditional baseline claim (round 50 §3), and now its target choice.

## 5. The confound, declared

The continuation study is a **continuer subset** (N=116/100 against 179/178), so a model
fitted to the full randomised population is predicting a survivorship subset. The level
discrepancy is visible directly — month-12 UPCR change −3.17 here against AURORA 1's
week-52 −2.65. Offset correction removes the level but not the selection: continuers may
also have different *trajectories*, which would flatter persistence. **Unquantified**, and
it is the main threat to this round's conclusion.

Per-visit N also declines (116→99, 100→87) and attrition is not modelled.

## 6. What this does not establish

- One trial family, aggregate arm means, 12 held-out points with heavy serial correlation.
- Round 50's fits were reused rather than refitted per horizon.
- It does **not** establish that no model could win here — only that these three variants,
  fitted this way, do not.
- It does **not** test the B-cell/autoantibody axis, which is `ROADMAP.md` Stage 1 and
  carries its own gate. **This round is about the renal model only.**

## 7. Consequence for the roadmap

The renal predictive model is finished as a predictor, and that is a result rather than a
defeat: it is now established with real data at four horizons rather than suspected.

- **Stage 0 (identifiability) is where the evidence is.** Nothing here touches those
  results — they concern what is measurable, not what is predictable.
- **Stage 1 remains the live question** and is now the only thing that can revive the
  predictive programme. It must be the B-cell/autoantibody axis, where the data is best,
  and it must clear a **level-based** out-of-sample bar. Not a slope.
- P1 harvest is now **2 of 16**.

*Research and simulation only. Not clinical or drug-development guidance.*
