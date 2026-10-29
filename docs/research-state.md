# Research state

> Living document. Last updated: 2026-10-28.

This file tracks where the BioModel-for-RLN thinking actually stands — distinct from the
aspirational README. It is the in-repo mirror of the Notion project hub.

## Current stage

First code committed. `experiments/rln-identifiability/` holds the model and the scripts
that reproduce rounds 40–44. `docs/system-design.md` is the canonical architecture spec.

**The project's position has shifted** from "build a biology world model" to a sharper
claim: *establish what is knowable in lupus nephritis.* The reasoning is in
`docs/positioning.md`. The short version: identifiability work needs no data, which turns
the project's binding constraint — no patient-level access — into its niche.

## What is established (computed, reproducible)

- **The structural zero.** Under feed-forward coupling an observable upstream of a
  degeneracy carries exactly zero information about it. Reproduced independently here;
  the apparent sensitivity scales as exactly `1/h`.
- **Feedback does not rescue it at realistic strength.** ~2× at implausible feedback,
  under 20% at plausible. **Severity is worth 437×** — the information lives in the
  clinical channel and scales with accumulated damage. This closes round 40's #1 open item.
- **The CD19⁻ long-lived plasma cell compartment is load-bearing.** After depletion the
  protected pool dominates what remains of autoantibody production — robustly, 93% at a
  5% baseline share and ≥99% above 40%. A B-cell count is a misleading surrogate.
  **Round 47 retracted the magnitude** (an earlier "7.9% fall" figure): it is set by an
  assumed parameter and spans 1.2–75.7% as that parameter moves.
- **The placebo falsifier is survived conditionally.** A realistic control-arm response
  rate needs between-patient CV ≈ 0.6–0.7, above the 0.50 pre-declared plausible. Anyone
  rejecting that spread should read the falsifier as still standing.
- **The variability–power squeeze is real and is not a threshold artifact.** Cohen's d
  falls 0.488 → 0.234 as CV rises 0.40 → 1.00 (p = 5.4e-6); required N rises 4.5×.

## Open questions

- **State representation / perturbation encoding** — answered in `docs/system-design.md`
  §3 and §5. The load-bearing result is that an intervention module must carry **internal
  state**; a stateless gate is provably dose-equivalent (measured 0.00% spread) and is an
  output head in disguise.
- **Data** — unchanged and still the binding constraint. No dataset gives molecular and
  clinical modalities longitudinally in the same subjects, so any bridge between them is
  an assumption. The route to patient-level data is institutional, not technical.
- **Evaluation** — the validation ladder is in `docs/system-design.md` §9. L6 (reproducing
  a 33–44% placebo response) is now partly addressed; the rest is outstanding.
- **Novelty** — partly resolved, and one bet died. See `docs/rounds/round-42-gates.md`:
  the pulse-schedule claim was **retracted** (PWM optogenetics is established since 2013).
  Two gates remain open (a trial-design review that needs a browser, and G0).

## Known holes

No Tregs — the highest-priority gap, since their omission previously changed the model's
dynamical class. No germinal centre, no T-cell help, no interferon. Refractoriness is not
yet a modelled state. Everything is synthetic at one operating point.

## Notes

Round docs live in `docs/rounds/`, each with its pre-registration where one exists.
Retractions are named rather than deleted — see round 42 (§D) and round 44 (§4, §5), both
of which correct earlier claims made in this repo's own history.

*Research and simulation only. Not clinical or drug-development guidance.*
