# Where this project can actually add value — lupus and RLN

2026-10-28 · strategic positioning, written after rounds 37–41

**Nothing here is computed. This is a judgement document.** It reasons over computed
results from rounds 18–41 but is itself a set of bets, and should be argued with rather
than followed.

---

## 1. The premise: stop building what exists

Round 41 settled the supply question. The field already has, deployed and tested:
structure prediction, protein and binder design, docking, molecule generation, DNA and
protein sequence foundation models. Those are the strong end of the field, with abundant
objective ground truth, and **we should consume them, never rebuild them.**

What nobody has, on either the preprint or the deployed side: **a transition operator for
biology.** No action axis, no iterable rollout, no time index. Plus, specific to our
target, no model with long-lived plasma cells, no model with Tregs, no model where
refractoriness is a *state* rather than a trial inclusion criterion.

So the opening is real. The question is which part of it a single researcher with no wet
lab, no patient-level data, and two CPUs can actually take.

## 2. The structural advantage — and it comes from the constraint

The honest competitive position: we cannot win on compute, data, or wet-lab throughput.
Arc, Recursion, Isomorphic and NVIDIA win all three permanently.

**But identifiability work requires no data.** You can prove a quantity is unmeasurable
without measuring it. That is the one class of result where our binding constraint —
no access to patient-level trial data (round 37 §8: every formal route keys on an
institutional signing official) — stops being a handicap.

And the project is already unusually good at it. Rounds 29–40 produced a run of
identifiability verdicts, five of which killed our own prior claims. `lessons.md` §0 is
right that this self-correction rate is the rarer asset. The field systematically does not
do this: round 41 found that even commercial benchmark reports omit simple baselines.

**So the bet is: be the group that says what is knowable in lupus, not the group that
builds the biggest model of it.**

### One honesty correction

I checked whether identifiability methodology is itself unoccupied. **It is not.** There
is a mature QSP identifiability literature — a 2021 *CPT:PSP*-adjacent perspective on
parameter identifiability in QSP, a 2024 workflow-recommendation paper, profile-likelihood
methods dating to Raue 2009, and `StructuralIdentifiability.jl` tutorials. **The method is
standard. We must not claim it as novel.**

What the searches did *not* surface is any systematic application of it to lupus or LN
models. That is a weaker claim than "new method" but it is the honest one, and it is
better news operationally: we do not have to invent methodology, only apply mature
methodology rigorously to a domain that has not had it. That is fast and cheap.

**[HYPOTHESIS — failed search, not proof of absence.]** Must be verified at source before
any novelty claim, per the round-34 rule.

## 3. Five candidate contributions, ranked by value per unit cost

### A. The LN measurement-design verdict — **do this first**

Convert rounds 29–41 into one actionable artefact: *what an LN trial should measure, and
what it should not bother with.* We already hold five concrete, non-obvious, computed
findings:

1. A direct biopsy activity index **does not** resolve the damage split — `r` stays at
   −0.982, `k_dam` SD improves 2% against a required 40% (R38). Whole-slide imaging is
   being oversold for this purpose.
2. **The placebo arm does most of the identifiability work** — single-arm costs 5.2× on
   `D_final` and halves route-3 detection. Two-arm per-visit trial data outranks
   whole-slide images (R38).
3. **Fortnightly autoantibody around a depletion window** pins the long-lived pool's decay
   `d_Pl` to ~1.5×, but **never** its replenishment `f_l` (σ 0.85 vs 0.30 required).
   Quarterly sampling straddles the nadir and costs ~2× (R40). A durability study should
   be designed for the first quantity and must not expect the second.
4. A proteinuria-lowering intervention needs **≥75% reduction within ~2 weeks** to reach
   single-digit VIF; timing dominates magnitude (50% at day 14 beats 90% at day 364).
   RAAS blockade gives 30–40% and is usually background from entry (R39).
5. **The damage split is identifiable only in advanced disease** — σ falls 438× from
   severity alone, versus at most 2× from any observable set (this session). Precisely
   backwards from where it would be clinically useful.

Why it is unoccupied: unglamorous, produces no model to demo, and nobody is rewarded for
telling a field what not to measure. Why it is valuable: it is directly actionable by
anyone designing an LN trial, and several of these findings would save real money.

**Cost: low — it is writing, not computing.** Most of it is already computed.
**Killed by:** finding that trial statisticians already know all five. Check before writing.

### B. The identifiability atlas of published lupus/LN models — **do this second**

Run structural and practical identifiability on every published lupus/LN model against the
observables trials actually collect. Targets already named: Budu-Grajdeanu 2010 (4 ODEs,
the only one with a damage variable, and its damage is reversible), Park/Mugundu/Singh
2025 (Takeda, CAR-T, no kidney compartment, SLEDAI by linear regression),
Volkova/Sokolov anifrolumab IFN model (17 ODEs, 52 parameters, 254 patients), the new
integrative type-1-IFN SLE model (PMC12256560), and the Simulations Plus SLE model
(excludes LN explicitly).

Output: *"these published models contain N parameters that no feasible trial could ever
identify."* A field-wide service, cheap, needs no new data.

**Cost: low-moderate.** `StrikePy` installs; no Julia, so `StructuralIdentifiability.jl`
and DAISY are unavailable — a real tooling constraint to solve first.
**Killed by:** discovering such an atlas exists. Verify before starting.

### C. The refractory-state model — **the most scientifically interesting**

The one genuinely unoccupied *modelling* ground specific to our target:

- **Long-lived plasma cells.** This session showed why it is load-bearing: 89% B-cell and
  97% short-lived-PC depletion buys a **7.9%** autoantibody fall, because **99.7%** of
  production comes from the CD19⁻ pool. That is the mechanism of serological relapse and
  a demonstration that B-cell count is a misleading surrogate. Round 37 found LLPCs absent
  from every model whose cell list could be retrieved.
- **Tregs.** Also absent field-wide. Their omission previously changed our dynamical class
  and manufactured a false cooperativity requirement (R27). This is our own largest hole.
- **Refractoriness as a state**, not an inclusion criterion.

**Cost: moderate.** **Killed by:** G0 — the claim that published models omit LLPCs still
rests on one low-reliability survey and must be verified at source first.

### D. The action-axis methodology — ~~the most transferable~~ **NOVELTY RETRACTED (R42)**

> **Killed by round 42.** Pulse-width-modulated optogenetics is established since
> 2013, and bioRxiv 2023.12.19.572411 screens intensity × period × duty cycle and
> attributes the differences to photoswitch activation/reversion kinetics — the same
> mechanism, named. The schedule-dependence phenomenon is **not ours**. What survives
> is narrow: NC3-as-audit-protocol for world models, which is a framing contribution,
> not a discovery. **It must not anchor a paper.** The 0.00% stateless result is
> standard periodic averaging and should be cited as such, not claimed.

Original text retained below for the record:


The NC3 law from this session generalises well beyond lupus: *a module whose action enters
through a stateless gate is provably dose-equivalent (measured 0.00% spread) and is an
output head in disguise; passing NC3 requires internal state with a relaxation time
comparable to the schedule's period.* Plus the schedule-at-fixed-dose test, which is sharp
enough to fail and which anti-BAFF correctly fails.

This speaks to the whole world-model field, not just immunology, and it is the piece most
likely to interest the AI-for-biology community.

**Cost: already largely done.** **Killed by:** finding the equivalent result in the
control-theory or pulsed-dosing literature, where it may well exist in another vocabulary.
**Check this before claiming it.** Periodic-averaging and pulsed-stimulation results are a
long-established field.

### E. The unimplemented transfer function — **right, but not yet**

Germinal centre output → circulating autoantibody → renal immune-complex deposition →
proteinuria, end to end. Round 37 flagged that nobody implements it. It is the real
scientific prize and the thing a cure argument would eventually need.

**Cost: high, and data-hungry in exactly the way we are blocked.** Defer.

## 4. The tension worth naming

A and B are the most defensible science. D is the most likely to get attention from the
AI-for-biology world. They are not the same thing, and pretending otherwise leads to doing
neither well.

The resolution that does not compromise either: **do A and B for the science, write D as
the methodological frame around them.** A paper whose contribution is "here is what is
knowable in lupus nephritis, established with a method that also exposes why every
published biology world model lacks an action axis" is one piece of work, not two.

What would compromise it: letting D's generality pull the work away from the lupus target
until nothing decides anything about RLN. Rule 5 exists for this.

## 5. The non-technical bottleneck, restated

Round 37 §8's conclusion still stands and still outranks every technical item here:
**the highest-leverage single action is not technical.** Vivli, NIDDK-CR, ImmPort's
controlled tier and UK Biobank all key on an institutional affiliation or a signing
official. A nominal academic affiliation unlocks patient-level data that no amount of
modelling substitutes for.

Note the asymmetry: **contributions A, B and D need no such access.** C needs some.
E needs a lot. That ordering is not a coincidence — it is the ranking.

## 6. What not to do

- **Do not build a foundation model.** Lost on compute and data before starting.
- **Do not put Geneformer, or any cell foundation model, in a surrogate slot** without a
  paired win over tuned PCA/scVI (R41).
- **Do not claim identifiability methodology as novel.** It is mature. Claim the
  application and the specific verdicts.
- **Do not claim the NC3 law** before checking periodic-averaging and pulsed-dosing
  literature.
- **Do not add model capability before the placebo falsifier passes.** A model that cannot
  reproduce a 33–44% placebo response has an unvalidated foundation, and `lessons.md`
  names building on that as spinning.

## 6b. Post-R42 ranking

With D's novelty retracted, **A is the clear lead** — by elimination as much as by
merit, which is the honest statement. A's own gate (the *Kidney Int Rep* trial-design
review, PMID 42006208) is **unresolved**, blocked by reCAPTCHA and publisher errors.
B holds weakly. C is still gated on G0, untouched since round 40.

## 7. Next three rounds, concretely

1. **R42 — clear the gates.** Three verification searches that gate three novelty claims:
   does a lupus identifiability atlas exist (B); is the NC3 law already in the
   periodic-averaging literature (D); G0, do published lupus models really omit LLPCs (C).
   Cheap, and each can kill a bet before it is paid for.
2. **R43 — the placebo falsifier.** Can the §4 architecture produce a 33–44% placebo
   response? If not, that is the headline and it outranks everything above.
3. **R44 — write A.** The measurement-design verdict, assuming R42 shows it is not
   already known.

*Research and simulation only. Not clinical or drug-development guidance.*
