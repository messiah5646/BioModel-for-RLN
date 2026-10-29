# Round 52 — the gap narrows with horizon but the model loses at every horizon tested (2026-10-29)

Computed. Code: `experiments/rln-identifiability/r52_crossover.py`. Reuses round 50's
fits on weeks 2–24, so nothing is re-fitted and no holdout is re-burned.

**Round 50's horizon account is weakly supported in DIRECTION and the crossover is NOT
demonstrated.** The model loses to linear extrapolation at every horizon available.

---

## Result

Error in units of each endpoint's own SD, so UPCR and eGFR are commensurable.

| horizon | ahead | best model | best trivial | ratio | trend |
|---|---|---|---|---|---|
| week 30 | 6 wk | 0.1253 | 0.0418 | **3.00** | — |
| week 36 | 12 wk | 0.1678 | 0.0620 | **2.71** | closing |
| week 42 | 18 wk | 0.1626 | 0.0871 | **1.87** | closing |
| week 48 | 24 wk | 0.2085 | 0.0999 | **2.09** | widening |
| week 52 | 28 wk | 0.1859 | 0.0629 | **2.96** | widening |

Linear extrapolation wins at **every** horizon. V2 is the best variant at every horizon,
confirming round 49's ranking.

## What this does and does not support

**Direction, weakly.** The ratio falls 3.00 → 1.87 from 6 to 18 weeks ahead, consistent
with round 50's account that trivial rules excel where the trajectory is near-flat.

**But the crossover is not demonstrated, and three things argue against over-reading it:**

1. **The trend is non-monotone** — it reverses at weeks 48 and 52.
2. **No crossover occurs within the data.** The ratio never approaches 1 in a 28-week-ahead
   window. Extrapolating the ratio to parity gives ~9 months ahead, but **that
   extrapolation is precisely the move round 50 showed fails**, so it locates a direction,
   not a point. Quoting "~9 months" as a crossover would repeat the error this round exists
   to avoid.
3. **Week 52 is the worst ratio of all** (2.96), and it is the visit rounds 48–49 fitted to.

**Consequence.** Round 48's long-horizon win stands as a single 2-point result on a
different scale, now with no supporting short-horizon trend to interpolate toward it. The
honest statement:

> **This model has no demonstrated predictive advantage over linear extrapolation at any
> horizon that has been measured.** Whether one exists beyond 28 weeks is an open question
> requiring per-visit data past week 52, which this trial does not provide.

## What this does not establish

- One trial, aggregate arm means, five horizons, heavy serial correlation between them.
- A single fit per variant; no refit per horizon and no seed variation.
- The SD-standardisation makes endpoints comparable but is not the only defensible choice.
- Absence of a crossover within 28 weeks is not evidence against one at 2 years.

## Next

1. Per-visit data past week 52 is the only thing that settles this — the extension study's
   per-visit series, or a longer trial from the 15 unharvested.
2. Until then, **the predictive programme does not clear its own bar** and should not be
   built on. See `docs/ROADMAP.md` Stage 1 gate.

*Research and simulation only. Not clinical or drug-development guidance.*
