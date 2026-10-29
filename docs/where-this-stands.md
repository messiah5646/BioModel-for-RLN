# Where this stands — after round 50 (2026-10-29)

Single source of truth. The previous version was written after round 40 and was nine
rounds stale; it also carried two claims since retracted and one environment fact since
made false. Superseded content is named in §Retractions rather than deleted.

## The one-paragraph version

**Strand A (methods, rounds 18–23):** a nonlinear forecaster with an unlabelled amplitude
gate beats linear baselines across four real corpora; JEPA/latent structure has never
earned its keep in-distribution. **Strand B (lupus nephritis, rounds 26–50):** a
mechanistic model used not to predict patients but to establish what measurement can and
cannot identify. The position after round 50 is sharper and less comfortable than after
round 40: the identifiability results are solid and now partly confirmed against real
trial data, **but on genuinely held-out short-horizon real trajectories a two-point linear
rule beats the five-parameter model.** The project's value is in impossibility and
identifiability results, not in predictive accuracy.

---

## What is established, computed and reproducible

- **The structural zero.** Under feed-forward coupling an observable upstream of a
  degeneracy carries exactly zero information about it. Reproduced independently in this
  repo; apparent molecular sensitivity scales as exactly 1/h, the signature of a true zero.
- **Feedback does not rescue it at realistic strength** (round 41 code). ~2× at implausible
  feedback, under 20% at plausible. **Severity is worth 437×.** Closes round 40's #1 item.
- **The CD19⁻ long-lived plasma cell compartment is load-bearing.** After depletion the
  protected pool dominates what remains of autoantibody production — 93% at a 5% baseline
  share, ≥99% above 40%. A B-cell count is a misleading surrogate.
- **A better structure can buy identifiability where more observables cannot** (round 49).
  Replacing activity-driven damage with activity-independent progression narrowed the
  treated arm's admissible interval 2.3×, from a span that changed sign to one that does not.
- **The placebo falsifier is survived, conditionally.** Realistic control-arm response
  needs between-patient CV ≈ 0.65–0.70, above the 0.50 pre-declared plausible. The window
  is **not empty** (round 46).
- **The variability–power squeeze is real and not a threshold artifact.** Cohen's d falls
  0.488 → 0.234 as CV rises 0.40 → 1.00 (p = 5.4e−6); required N rises 4.5×.
- **The squeeze is irreducible by enrollment criteria** on baseline renal function
  (round 46). Baseline UPCR moves control response 9.5% → 46.3% across quartiles while
  leaving effect size flat — strongly **prognostic, non-predictive**. No rule beats
  enrolling everybody once screening burden is counted.
- **Routes 2 and 3 are not separable from year-1 two-endpoint data** (round 48, real data).
  A 1.87× range of route-3 strength fits year-1 endpoints exactly and indistinguishably,
  over which the predicted long-term slope spans 1.474 mL/min/yr and changes sign.
- **G0 partially cleared at source** (round 49): the Bull Math Biol NETosis model's four
  equations contain no long-lived plasma cells, no CD19⁻ plasma cells and no Tregs, with
  B-cell involvement explicitly omitted by its authors. Three targets remain.

## Retractions and corrections

| claim | status | where |
|---|---|---|
| "deep B-cell depletion buys a ~7.9% autoantibody fall" | **magnitude retracted** — set by an assumed parameter, spans 1.2–75.7% | R47 |
| the pulse-schedule / NC3 dose-sparing finding was novel | **novelty retracted** — PWM optogenetics established since 2013 | R42 |
| the model "beats both trivial comparators" on held-out data | **corrected** — true only at a 2-year horizon on a 2-point target; reverses at ≤6 months | R50 |
| multiplicative parameter structure causes the squeeze | **downgraded** to "contributing factor, not established" — the squeeze appears under additive jitter too | R46 |
| CV 0.70 satisfies both placebo-realism and ≥10pp separation | **corrected** — single-seed artifact; 8 seeds put the window at 0.65–0.70 | R44 |
| cosine-VIF as the instrument for observable-set comparison | **retracted** — not monotone, reported impossible sub-1.0 gains | R41 code |
| round 48's long-horizon predictive win | **falsified (R53)** — it scored a *slope*; on per-visit levels at 15–27 months the model loses to persistence (0.239 vs 0.163) | R53 |
| route 3 as implemented | **inconsistency confirmed** — cosmetic, affects the observable not the damage driver. Repairing it did **not** improve held-out prediction | R49/R50 |

## The current live problem

**Round 53 settles the predictive question for the renal model: no advantage over a
persistence rule at any horizon from 6 weeks to 27 months, which is the full span of
available data. The gap narrows with horizon (3.00 → 1.47) but never crosses.** Stage 1 of
`ROADMAP.md` — the B-cell/autoantibody axis, scored on levels rather than slopes — is now
the only route that can revive the predictive programme.

Round 50 gave the first sign of this. On 16 genuinely unseen points (weeks 30–48, fit on
weeks 2–24), held-out RMSE was: linear extrapolation **1.14**, persistence **1.31**, best
mechanistic variant **1.77**. On UPCR the gap is 4.7×. Only one variant even gets the sign
of the week 24→48 trend right, and it understates the magnitude ~4×.

**The scoped claim that survives:** this model has no short-horizon advantage over linear
extrapolation; whatever predictive value it has appears beyond a year. That predicts a
crossover, and locating it is the highest-value next test.

## Where the project adds value (round 45 positioning, unchanged by round 50)

Identifiability work needs no data, which turns the binding constraint — no patient-level
access — into the niche. Ranked: (A) the LN measurement-design verdict, (B) an
identifiability atlas of published lupus models, (C) the refractory-state model.
**Round 50 strengthens this**: predictive accuracy is where the project loses; impossibility
and identifiability verdicts are where it wins.

## Standing holes

No Tregs — the highest-priority gap, since their omission previously changed the dynamical
class. No germinal centre, no T-cell help, no interferon. Refractoriness is not a modelled
state. Everything outside rounds 48–50 is synthetic at one operating point.

## Standing environment facts

- **There is now a repo.** `messiah5646/BioModel-for-RLN`. The container remains ephemeral,
  so commit as work is produced. This supersedes the previous "no repo" note.
- **Real data is in-repo**: `data/real/aurora1_per_visit.csv` (12 visits) and
  `aurora2_per_visit.csv` (months 12–36, continuer subset — read its header traps). Round 35's copy was lost to the container; this one will not be.
- **Test reachability per host AND per tool.** WebFetch and the shell take different
  routes; the EU Clinical Trials Register is unreachable by curl and reachable by WebFetch.
- **EU CTR is the reliable source** for per-arm per-visit endpoints: 64 lupus nephritis
  trials, **16 with results sections**, 15 still unharvested.
- `huggingface.co` and `api.typesafe.ai` unreachable. `pypi.org` and
  `raw.githubusercontent.com` reachable.
- Background jobs stop at a 2-hour limit; checkpoint per cell.
- **The pre-registered stopping rule has caught a bug in three consecutive rounds**, each
  of which would otherwise have produced publishable-looking numbers from failed fits.

## Blocked, with workarounds

See `docs/workarounds.md` — seven blockers, each with a tested status. Four dissolve once
the work is scoped so they are not dependencies. The two that genuinely bind are harvesting
jobs on a confirmed-reachable source.

## Next, in order

1. **Locate the crossover horizon** (round 50 §6). The one result here that generalises.
2. Harvest a second trial's per-visit series — 15 are available — for a whole-trial holdout.
3. Add Tregs before building anything else on the humoral layer.
4. Clear the remaining three G0 targets by reading primary equations.
5. Fit `g_l`/`g_s` to paired autoantibody and vaccine-IgG titres; it swings round 47's
   retracted headline across its whole range.

*Research and simulation only. Not clinical or drug-development guidance.*
