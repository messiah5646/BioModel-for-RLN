# Roadmap — from identifiability verdicts to a validated subsystem

Written 2026-10-29, after round 52. Forward-looking companion to `positioning.md`
(which ranks *where we add value*) and `where-this-stands.md` (which records *what is
established*). This document says what the model is **for**, in stages, with a gate
between each.

**The stated goal** is a cure for lupus/RLN, or failing that a realistic simulator of the
immune system and body on which drugs and interventions can be tested. This roadmap takes
that goal seriously by refusing to treat it as a single step.

---

## The bar, stated honestly before anything else

- **No general-purpose immune simulator is quantitatively validated against measured human
  immune trajectories.** Not one, from anyone, at any budget (round 37). The Immune Digital
  Twin Blueprint takes dimensionless percent-chance inputs, so it cannot be initialised
  from patient measurements even in principle. The best-validated immune model fitted to
  real human data at scale is 17 ODEs with interferon readouts only, and its strength is
  its narrowness.
- **No cell foundation model beats a tissue lookup at perturbation prediction**
  (Ahlmann-Eltze 2025). Round 50 reproduced the same asymmetry on real trial trajectories.
- **Our own model has no demonstrated predictive advantage over linear extrapolation at any
  measured horizon** (round 52).

So the simulator is not a thing we are behind on. It is an open problem, and the only
domains where this class of problem has been solved — structure prediction above all —
were solved by being **narrow with abundant ground truth**. That is the shape of this plan.

**The binding constraint is not modelling.** No dataset gives molecular and clinical
measurements longitudinally in the same subjects. Any model bridging them does so by
assumption, and that assumption is its principal scientific risk. More rounds cannot fix a
missing dataset.

---

## Stage 0 — Identifiability verdicts · NOW → ~3 months · no new data required

**Deliverable.** The LN measurement-design verdict: what an LN trial should measure and
what it should not bother with. Five computed findings already in hand (round 38's biopsy
result, round 38's placebo-arm result, round 40's sampling-density result, round 39's
proteinuria-timing result, this session's severity-dominance result), plus round 48's
route-2/route-3 inseparability established on real data.

**Why this is upstream of the simulator, not a detour.** If a quantity is unmeasurable, no
simulator can ever learn it. Stage 0 tells Stage 1 which parameters are worth grounding.

**Gates before publication.** Gate A (prior-art check on trial-design reviews — six routes
failed, needs a browser or institutional access) and G0 (3 of 4 targets remaining).

**Kill criterion.** If trial statisticians already know all five findings, contribution A
is not a contribution. **This is checkable and unchecked**, and it is the single largest
risk to Stage 0.

---

## Stage 1 — One validated subsystem · ~6–12 months · THE decisive milestone

**Scope.** The B-cell / autoantibody axis only: `B`, `Ps`, `Pl`, `Ab`. Nothing renal,
nothing systemic.

**Why this axis.** It has the best data in all of lupus. CAR-T trials report B-cell counts,
autoantibody titres **and** vaccine IgG at matched visits — which is exactly the paired
measurement `lessons.md` already identified as the cheapest decisive one, and exactly what
separates the CD19⁺ renewable pool from the CD19⁻ archive.

**Why it matters.** This is the first point at which *"test an intervention in silico"*
means anything. And the interventions acting on this axis are the live ones: CD19
depletion, BCMA-directed therapy, CAAR-T, anti-BAFF.

**The gate, and it is hard.** The subsystem must **beat linear extrapolation and
persistence out-of-sample, at the horizon it claims**, on data never used in fitting.
Round 52 says the current model does not clear this at any horizon ≤ 28 weeks.

**Kill criterion.** If a model of the axis with the best available data cannot beat trivial
rules, **the predictive programme is dead and Stages 2–4 must not be attempted.** The
project then becomes the identifiability programme, full stop. That is a real possible
outcome and naming it now is what makes the gate meaningful.

---

## Stage 2 — Negative regulation and germinal centre · ~12–24 months

**Tregs first**, because their omission already put the model in the wrong dynamical class
and manufactured a false cooperativity requirement (round 27). Then the germinal centre,
where tolerance breakdown is actually computed in the literature.

Each addition needs **its own validation data and its own Stage-1-style gate**. Stacking
mechanisms is not monotonically better — round 24 found a stacked hypothesis scoring worse
than a simpler one.

---

## Stage 3 — The unimplemented transfer function · the scientific prize

Germinal-centre output → circulating autoantibody → renal immune-complex deposition →
proteinuria, end to end. Round 37 found this implemented by nobody. It is where
cure-relevant reasoning would live.

**Requires paired molecular and clinical longitudinal data in the same subjects**, which
does not currently exist in any public dataset. So this stage is gated on data access, not
on modelling effort.

---

## Stage 4 — Multi-system · requires an institution

Not a solo project, and saying so is not pessimism. It is what every comparable programme
has required.

---

## Standing processes — these run across sessions, not within a round

| | process | definition of done | status |
|---|---|---|---|
| **P1** | **Trial-data harvest.** All 16 EU CTR lupus-nephritis trials with results sections → per-visit series in `data/real/` | 16/16 harvested, each with its source citation | **1/16.** AURORA 1 done |
| **P2** | **Identifiability atlas.** Structural + practical identifiability for every published lupus/LN model against the observables trials actually collect | every located model has a verdict | not started; needs a structured search first |
| **P3** | **G0 clearance.** Read primary model equations, not reviews | 4/4 targets verified at source | **1/4.** Route proven to work |
| **P4** | **Labelled triage set.** 50–100 abstracts hand-labelled for the model-state-variable question | 50+ labelled, keyword baseline scored | not started; n=9 pilot exists |
| **P5** | **Institutional access.** ImmPort (free, same-day), Synapse syn64064827 (days) | one patient-level source in hand | not started — **the highest-leverage item on this page** |
| **P6** | **Crossover characterisation.** Per-visit data past week 52 | model-vs-trivial ratio measured beyond 1 year | blocked on P1 |

**P5 outranks everything else on this page.** It is not technical, it has not moved, and it
is what separates Stage 3 from being permanently out of reach.

---

## What would falsify the whole programme

1. **Stage 1's gate fails** — the predictive half dies, and that is a legitimate finding
   rather than a defeat.
2. **Stage 0's kill criterion fires** — the measurement-design findings turn out to be
   known, and the identifiability half loses its lead deliverable.
3. **Both** — the project reduces to a well-documented set of negative results, which is
   still more than most efforts in this area produce, and should be written up as such.

## The discipline that makes any of this worth believing

Seven retractions are recorded in `where-this-stands.md`, most found by attacking our own
best results. The pre-registered stopping rule has caught a bug in three consecutive
rounds, each of which would otherwise have produced publishable-looking numbers from failed
fits. **That capability is the project's most transferable asset** and it should survive any
change of direction this roadmap triggers.

*Research and simulation only. Not clinical or drug-development guidance.*
