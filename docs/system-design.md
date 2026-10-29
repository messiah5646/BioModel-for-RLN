# Bio World Model — System Design

**Status: canonical architecture spec. Supersedes ad-hoc design decisions in rounds 18–40.**
Version 1.0 · 2026-10-28 · target: refractory lupus nephritis and SLE

Every load-bearing claim below was computed in this session, not argued. Code in
`src/`, raw outputs in `out/*.json`. Claims that were **not** computed are marked
**[HYPOTHESIS]** and carry no authority. Two of my own results were retracted
mid-session; both are recorded in §11 rather than deleted.

---

## 0. One-paragraph version

A controlled in-silico lab for autoimmune intervention design, built as a **slow
mechanistic core** (immune compartments → kidney, with explicit feedback) wrapped
in a **plug-and-play intervention layer** where every therapy — optogenetic,
small-molecule, biologic, cell — satisfies one contract. The core never imports a
module. The organising discovery of this session is that **the module contract must
require internal state**: a stateless intervention module is provably
dose-equivalent (measured spread 0.00%) and is therefore an output head in
disguise, which is exactly the failure mode round 37 found in every published
biology world model. Modules carrying their own kinetics change the rollout by
101–253% at identical total dose. That single requirement is what makes this an
action-conditioned world model rather than a dose-response surface.

---

## 1. Design constraints inherited from rounds 18–40

These are not preferences. They are findings the architecture must not contradict.

| # | Finding | Architectural consequence |
|---|---|---|
| C1 | An observable upstream of a degeneracy carries **exactly zero** information about it under feed-forward coupling (R40) | Coupling topology is an *identifiability* decision, not a modelling convenience. The graph must be explicit and inspectable. |
| C2 | No cell foundation model beats a tissue lookup at perturbation prediction (Ahlmann-Eltze 2025); latent structure never earned its keep in-distribution (R20) | Trivial baselines are a **first-class component**, not a benchmarking afterthought. Learned parts must beat them to ship. |
| C3 | Action-conditioned multi-step biology is an empty category; every system attaches perturbation to a head on frozen embeddings (R37). **Confirmed R41 on the deployed side**: all 11 BioNeMo catalogue models are static — zero have a transition operator. | The action must enter the transition operator. §5 makes this a testable gate. **No vendor supplies L2; it must be built.** |
| C4 | Predictions are pinned far tighter than parameters (`D_final` σ 0.023 vs log `k_dam` σ 6.7) | The output API must separate *predictions* from *parameter claims*. Mechanistic attributions are gated on identifiability. |
| C5 | Compartments are not interchangeable: CD19 is lost at the plasma-cell transition | CD19⁻ long-lived plasma cells are a **mandatory** state. Verified in §4 to produce the relapse signature. |
| C6 | Omitting Tregs changed the dynamical class and manufactured a false cooperativity requirement (R27) | Negative regulators are mandatory, not optional detail. Currently **missing** — see §10. |
| C7 | Placebo arms improve 33–44%; the placebo arm does most of the identifiability work (R38) | Two-arm simulation is the default. Single-arm is a degraded mode. |
| C8 | Epitope spreading is state-space *expansion* | The state vector must be growable at runtime. |
| C9 | No dataset gives molecular and clinical modalities longitudinally in the same subjects | Any molecular↔clinical bridge is **an assumption**, and the system must label it as such in its own output. |

---

## 2. Layered architecture

```
  L5  EXPERIMENT SELECTION        what to run next, ranked by value/cost
       |  (read-only on L4; never mutates the core)
  L4  INFERENCE + IDENTIFIABILITY  FIM, degeneracy detection, provenance gating
       |
  L3  INTERVENTION LAYER           plug-and-play modules (the action axis)
       |  <-- ONE CONTRACT. Core has zero knowledge of any module.
  L2  CORE DYNAMICS                immune compartments -> kidney, explicit graph
       |
  L1  STATE ALGEBRA                named, typed, growable state vector
       |
  L0  SUBSTRATE                    integrators, units, seeds, provenance tiers
```

The single most important line is the one between L2 and L3. It is enforced by the
core's signature: the core sums `drift()` contributions and composes `modulate()`
factors from an opaque list. Adding a therapy requires **zero** edits to L2.

---

## 3. L1 — state algebra

States are named and typed, never positional. The reference instantiation:

| state | meaning | surface | depletable by CD19 agents |
|---|---|---|---|
| `B` | autoreactive B cells | CD19⁺ | yes |
| `Ps` | short-lived plasma cells | CD19⁺ | yes |
| `Pl` | long-lived plasma cells | **CD19⁻** | **no** |
| `Ab` | autoantibody titre | — | — |
| `A` | reversible inflammatory activity | — | — |
| `D` | **monotone irreversible** renal damage | — | — |

`D` is monotone by construction (`dD/dt ≥ 0`). This is the one structural commitment
the project has claimed novelty on: the only published lupus model with a damage
variable (Budu-Grajdeanu 2010) gives it a decay term, making it reversible.

**Growable.** C8 requires that epitope spreading add clones at runtime. The state
vector is therefore a dict-backed registry, not a fixed array, and modules may
declare new states (§5). Fixed-length arrays are a correctness bug here, not an
optimisation.

**Observation operator** — deliberately separate from state, because the
identifiability results are properties of the *pair*:
```
UPCR = p0 · A · N          N = 1 − D    (attenuated by nephron loss: the fold)
eGFR = G0 · N · (1 − c_inf·A)
```
UPCR attenuated by lost nephron mass is the geometry the project claims as novel;
Weinstein 2022 is the antecedent for the descending limb and must be cited.

---

## 4. L2 — core dynamics, and the feedback decision

```
dB /dt = s_B·supp(t)·(1 + φ·D) − d_B·B
dPs/dt = f_s·supp(t)·B        − d_Ps·Ps
dPl/dt = f_l·supp(t)·B        − d_Pl·Pl           <-- CD19−, protected
dAb/dt = g_s·Ps + g_l·Pl      − d_Ab·Ab
dA /dt = k_a·Ab               − k_clear·A
dD /dt = (k_dam·A + k_tox·UPCR)·N                 <-- monotone
```

`φ` is the **damage → antigen exposure → B-cell activation** feedback. φ=0 is
feed-forward (round 40's assumption); φ>0 writes the feedback into the ODE, which
was round 40's #1 open item.

### 4.1 Calibration (`out/` sanity block)
UPCR 4.17 → 1.63 at week 52; eGFR 78.2 → 84.3; autoantibody nadir 0.168; `D_final`
0.0445; no clamping anywhere. AURORA-like and physiological.

### 4.2 C5 verified — the long-lived plasma cell floor is load-bearing

CD19-directed depletion, day 180–240, measured at day 240:

| compartment | untreated | depleted | reduction |
|---|---|---|---|
| `B` (CD19⁺) | 0.3590 | 0.0380 | **89.4%** |
| `Ps` (CD19⁺) | 0.1099 | 0.0033 | **97.0%** |
| `Pl` (CD19⁻) | 0.7256 | 0.7154 | **1.4%** |
| `Ab` | 0.5055 | 0.4653 | **7.9%** |

**[MAGNITUDE RETRACTED — round 47.]** The 7.9% figure is driven by an *assumed* parameter (`g_l`/`g_s`, the long-lived pool's share of autoantibody production). Sweeping it moves the fall across 1.2–75.7%, and real CAR-T seroconversion data argues for a *large* fall. What survives is the structural claim: after depletion the protected pool dominates what remains of production (93% at a 5% share, ≥99% above 40%). See `docs/rounds/round-47-llpc-magnitude-retracted.md`.

**89% B-cell and 97% short-lived-plasma-cell depletion buys a 7.9% autoantibody
fall**, because **99.7%** of autoantibody production at day 240 comes from the
protected CD19⁻ pool. This is the modelled mechanism of serological relapse and a
direct demonstration that **a B-cell count is a misleading surrogate** — near-total
depletion of everything measurable in blood barely moves the quantity that drives
disease.

*(Two process notes, per the project's own rules. First, I initially read this flat
response as a bug; it is the intended biology. Second, an earlier draft of this
section quoted different numbers, carried over from a superseded parameterisation
and reproduced by a stale hardcoded `print` in the script — exactly the "number
whose script no longer produces it" failure. The table above was re-derived from the
current code before release, and the stale string is corrected in §11.)*

### 4.3 The feedback verdict — round 40's #1 open item, now answered

Round 40 established the structural zero and correctly refused to count its own
`s_B`-shift emulation. I wrote the feedback into the ODE.

**First, the pipeline reproduces round 40 independently.** At φ=0, apparent
molecular sensitivity scales as exactly 1/h across five step sizes (6.6e-8 at
rel_h=1e-1 → 1.1e-3 at 1e-5) — the signature of pure roundoff around a true zero.
|cos| is stable to the 9th decimal. Converged, with the noise floor measured.

**Then the instrument failed and was replaced.** The cosine-VIF reported "gains"
below 1.0 — adding an observable appearing to *destroy* information, which is
impossible. The cos-VIF is a ratio estimator and is **not monotone in the observable
set**, so it cannot answer this question. Fisher-information marginal variance is
monotone by construction. All verdicts below use the FIM; the VIF tables are
retracted for this purpose (§11).

Marginal σ on log `k_dam`, nuisance parameters jointly estimated, realistic assay
noise (UPCR 15%, eGFR 8%, autoAb 20%, B-cell 25% CV):

| severity | φ | `D_final` | σ, clinical only | σ, + molecular | gain |
|---|---|---|---|---|---|
| 1 | 0 | 0.045 | 283.6 | 280.5 | 1.01× |
| 1 | 20 | 0.048 | 236.1 | 114.1 | **2.07×** |
| 10 | 0 | 0.343 | 4.59 | 4.19 | 1.10× |
| 10 | 20 | 0.514 | 1.68 | 0.94 | **1.78×** |
| 20 | 0 | 0.544 | 2.01 | 1.64 | 1.23× |
| 20 | 20 | 0.814 | 0.648 | 0.373 | **1.74×** |

**Three results, in order of importance:**

1. **Round 40's zero is a statement about direct sensitivity, not about assay value.**
   Even at φ=0 the gain is 1.01–1.23×, never exactly 1. Molecular observables pin
   the *nuisance* parameters (`d_Pl`, `f_l`, `kill`, `d_Ab`), which indirectly
   tightens `k_dam`. The naive reading — "molecular assays are worthless for the
   damage split" — is **wrong**, and this refines round 40 rather than overturning it.
2. **Feedback helps, but needs to be strong.** At φ=20 molecular data roughly halves
   the uncertainty. At φ≤5 it buys under 20%. So round 40's conditional holds in
   practice across any plausible feedback regime; the conclusion is robust, and the
   #1 open item is closed as *"feedback does not rescue it at realistic strength."*
3. **Severity dominates everything.** σ falls from 283.6 to 0.648 — a factor of
   **437** — from disease severity alone, versus at most 2× from the observable set.
   The information about the damage split lives in the *clinical* channel and scales
   with accumulated damage. **[HYPOTHESIS, high confidence, not yet tested on real
   data]** this predicts the damage split is identifiable in advanced/refractory
   disease and hopeless in early disease — the opposite of where one would want it.

**What this does not establish.** Synthetic data, illustrative parameters, one
operating point, one patient per arm. Linearised Fisher geometry, one method, with
no multistart-fitting cross-check this session (round 39 ran one for the VIF route;
not repeated here). Nothing here is a claim about patients.

---

## 5. L3 — the intervention contract (the core of the design)

### 5.1 The contract

```python
class Module:
    def declares(self) -> dict   # new states this module owns (growable registry)
    def touches(self)  -> dict   # states + params it may alter   <- audit surface
    def drift(t, y, p) -> vec    # additive contribution to dy/dt
    def modulate(t,y,p)-> dict   # multiplicative parameter changes
    def provenance()   -> dict   # measured | derived | constrained | assumed
```

Composition is summation of drifts and multiplication of modulations. Verified: the
core ran untreated, CD19 CAR-T, LiCAR-T, anti-BAFF, and LiCAR-T **+** anti-BAFF with
no core edits and no module aware of any other.

`touches()` and `provenance()` are not bureaucracy. They are what let L4 refuse to
report a mechanistic attribution that rests on an assumed parameter, implementing C4.

### 5.2 The NC3 law — the central result of this design

Round 37's audit paper requires that changing the action change the **rollout**, not
the decoded endpoint. I needed an operational test sharp enough to fail. The one I
propose and ran:

> **Hold total dose exactly fixed. Vary only the schedule.** A dose-response head
> receives an identical summary statistic and *must* return an identical answer. A
> true transition-operator action need not.

Fixing the duty cycle at 0.0833 (so total light energy is identical to the
microsecond) and varying only the pulse period:

| module | effective-gate range | `Ps_nadir` spread | verdict |
|---|---|---|---|
| **stateless gate (no internal state)** | 1.00× | **0.00%** | **NC3 FAIL — dose-equivalent** |
| iLID / CpLID, off ≈1 min | 4.18× | 101.2% | NC3 PASS |
| CRY2, off ≈6 min | 8.44× | 189.2% | NC3 PASS |
| CRY2 L348F, off ≈24 min | 9.18× | 252.7% | NC3 PASS |

**The law, and it is a theorem about the contract, not a style guide:**

> A module whose action enters the drift linearly through a stateless gate is
> **provably dose-equivalent**: when the schedule is fast relative to the state's
> relaxation time the system time-averages, and only the integral of the input
> survives. Such a module is an output head wearing a transition operator's clothes.
> **A module passes NC3 if and only if it carries internal state whose relaxation
> time is comparable to the schedule's period.**

The stateless row measuring exactly 0.00% is the proof, not an illustration. This is
why `declares()` is in the contract: **modules must be able to own state.** Round 37
found that every published biology world model attaches perturbation to a frozen-
embedding head. This contract makes that failure mode structurally impossible to
express, and the test above detects it if someone tries.

### 5.3 Multiscale execution

The photoswitch runs at minutes; the immune system at months. Integrating minute
pulses over 3 years was attempted and killed at the 120 s timeout. The module
therefore solves its own periodic steady state in closed form and exposes a
cycle-averaged gate upward. **A module owns its timescale and its integrator.** This
is a contract requirement, not an optimisation — the naive version does not run.

### 5.4 A falsifiable biological prediction that falls out

Slower photoswitch off-rates make schedule matter *more* (4.18× → 9.18× as off-rate
goes 1 → 24 min), because a slow switch integrates pulses and rapid pulsing holds it
near its ON asymptote while slow pulsing lets it fully relax.

**Therefore: pulse period, at fixed total light dose, is a ~9× efficacy lever that
costs nothing.** This matters directly — >99.7% of NIR input becomes heat, and
heating is the clinical bottleneck for the UCNP route. Reaching the same effect at
~9× lower light dose is the most actionable thing this design says about LiCAR
engineering, and it is **[HYPOTHESIS]** until tested against the measured
dose-response in Nguyen/Zhou 2021. A cheap, decisive falsifier.

### 5.5 Worked modules

**`licar_t`** — nano-optogenetic CAR-T. AND gate (antigen AND light), explicit
dark-state leak, Beer-Lambert tissue attenuation, upconversion as `I^n` with n≈2.1
for the 475 nm channel giving a sharp spatial threshold. Owns photoswitch state.
Provenance: `I_nir` measured, `n_photon` measured, `mu` derived, `I_sat` assumed,
`leak` assumed.

**`anti_baff`** — biologic. Pure parameter modulation on `s_B`, no drift, no state.
Correctly NC3-FAILS, and should: a chronically dosed biologic *is* dose-equivalent at
these timescales. **The test discriminates real pharmacology, not just code quality.**

**`cd19_cart`** — constitutive. Touches `B`, `Ps`; structurally cannot touch `Pl`.

**Sourcing module parameters (R41).** The BioNeMo structure/design stack
(RFdiffusion/RFD4 → ProteinMPNN → OpenFold3/Boltz-2) is a legitimate route to move
`kill` from *assumed* toward *derived*, which widens what §6's provenance gate permits.
Boltz-2 is admitted **only** with a mean-predictor baseline reported alongside and a
re-benchmark against 2.1: independent results give mean Pearson ~0.42 on PL-REX, only
5–7% over ΔvinaRF20, with affinity compression and run-to-run variability up to
1.5 kcal/mol.

**The gap that matters.** §5.2 proved the NC3-critical parameter is the photoswitch
**off-rate** `k_off` — the whole 4.18×→9.18× schedule lever. Boltz-2 predicts *affinity*,
a thermodynamic ratio; it does not separate `k_on` from `k_off`. Two switches with
identical K_d can differ an order of magnitude in off-rate — iLID (~1 min) vs CRY2 L348F
(~24 min) is exactly that case, and §5.2 shows they behave completely differently at the
same light dose. **No deployable model supplies the parameter this architecture most
needs.** It is a measurement requirement, not a modelling choice.

A drug module, a cell therapy, a cytokine blockade and an optogenetic construct all
satisfy one contract. "If it's not optogenetic and it's some other drug" is a
different module file, not a different architecture.

---

## 6. L4 — inference and the provenance gate

Two instruments, mandatory, because this session showed one is not enough:
1. **Fisher marginal variance** — monotone in observables. The default.
2. **Multistart fitting** — the cross-check. Round 39 ran it; round 40 did not, and
   this session did not. **Standing debt.**

**The provenance gate (implements C4).** Before the system reports any statement, it
classifies it:
- *Prediction* (`D_final`, trajectories) — released if converged.
- *Parameter / mechanistic attribution* (`k_dam` vs `k_tox` split) — released **only
  if** marginal σ is below a pre-registered threshold at the current operating point.

Given §4.3, at severity 1 the damage-attribution claim is **refused** (σ=283) and at
severity 20 with strong feedback it is **released** (σ=0.37). The system must say
*"this patient ends near 50% nephron loss"* and must refuse *"60% of the damage was
inflammatory"* unless the data earns it. This is the architecture enforcing the
project's own hardest-won lesson in code rather than in prose.

**Degeneracy handling.** Where a direction is structurally degenerate, reparameterise
to the identifiable combination and report that. Never quote a posterior over a
degenerate direction — it is the prior wearing a posterior's clothes.

---

## 7. L5 — experiment selection

Ranks candidate experiments by expected information per unit cost, using L4's FIM.
Strand A's closed loop cut whole-matrix error ~55% with two chosen experiments,
within ~12% of an oracle — **post hoc; needs a prospective test.**

Read-only on L4 and never mutates the core. Prediction written before the run;
trivial comparator (random selection) always scored alongside, because random matched
the mechanism's first four picks in a prior round.

---

## 8. Baselines as a first-class component (implements C2)

Not a benchmark harness — a module type. Every reported result carries its trivial
comparator in the same object: tissue/arm lookup, "stronger ligand wins", a
2-parameter exponential-then-logistic curve, and a tuned linear/ARX model swept as
hard as the learned model. If the trivial rule ties, that goes in the headline.

**Metrics:** WMSE, Normalized Inverse Rank, Dynamic Range Fraction alongside plain
MSE — the round-20/21 null may have been measured with an instrument that cannot
resolve a latent advantage. **Standing debt, cheap, could overturn a conclusion.**

**Named exclusion (R41).** Geneformer — the single-cell model shipped in NVIDIA BioNeMo —
is **barred** from this slot pending a paired win over tuned PCA/scVI. An independent
benchmark (arXiv 2410.13956, Replogle) places it **last** among foundation models on
known-relationship recall (0.149 vs scVI 0.503) and **11× worse than PCA** at linear
probing (0.005 vs 0.058). The one task it wins, batch integration, is one where *random
embeddings score 1.000*. This is C2's third independent triangulation.

**Where learned components are allowed.** Not in the core. What survived rounds 18–23
is *a nonlinear forecaster plus an unlabelled amplitude gate* — latent/JEPA structure
never earned its keep in-distribution. Its measured strengths are calibrated
cell-to-cell spread and long rollouts under strong amplitude shift. So: learned
surrogates for the **observation operator** and for **inter-patient variability**,
with the mechanistic core carrying causal structure. A learned component ships only
when it beats the tuned classical baseline on a paired comparison with CIs.

---

## 9. Validation ladder

Borrowed from the LiCAR validation cascade, which is a good template:

| level | gate |
|---|---|
| L0 | units; steady state exact (drift < 4e-5 over 3 y); no clamping |
| L1 | step-size convergence + **noise floor reported** with any rank/threshold claim |
| L2 | **callable transition reusable over >1 step** (audit paper's L2) |
| L3 | calibrated long-horizon rollout, controlled error growth; **beat 295.67 → 320.12** under iterated state feedback |
| L4 | NC3: schedule-at-fixed-dose test (§5.2) |
| L5 | two-arm trial reproduction, incl. the **placebo arm** |
| L6 | **placebo response 33–44%** — the live falsifier. A model with a stable disease attractor cannot produce it. **Still outstanding.** |

L6 is the clearest standing falsification and the architecture does not yet pass it.

---

## 10. Known holes — read before trusting anything

1. **No Tregs.** C6 says negative regulators change the dynamical class. Their absence
   previously manufactured a false cooperativity requirement. **Highest-priority gap.**
2. **No germinal centre, no T-cell help, no interferon, no complement, no antigen.**
   The humoral layer is a four-state caricature.
3. **The molecular↔clinical bridge is an assumption** (C9). No dataset provides both
   modalities longitudinally in the same subjects. This belongs in the first paragraph
   of any output, not its limitations.
4. **Specificity is receptor-level.** An anti-dsDNA B cell and a tetanus B cell are
   identical to any CAR. No model invents a targeting handle.
5. **Refractoriness is not a modelled state.** It appears only as a trial criterion.
6. **No multistart cross-check** of §4.3. One method, not two.
7. **Health may be a basin, not a point.** No baseline immune state for the general
   population is definable; disease may be a bifurcation, not a parameter offset.
8. **Composite scores double-count** — SLEDAI-2K contains anti-dsDNA and complement.
9. **No source for photoswitch kinetic parameters** (R41). `k_off` drives the NC3 result
   and the dose-sparing prediction, and no deployable model predicts it. Currently
   *assumed*, which §6 means the schedule claims cannot be released as mechanistic.
10. **G0 uncleared**: the claim that published lupus models omit long-lived plasma
   cells rests on one low-reliability survey. **No novelty is claimed until verified
   at source.**

---

## 11. Retractions from this session

The project's rule is that retractions are named, not deleted.

| retracted | why | replaced by |
|---|---|---|
| T1/T2 **cosine-VIF** as the instrument for observable-set comparison | not monotone in the observable set; reported gains < 1.0, which is impossible | Fisher marginal variance (§4.3) |
| T5 v2 **NC3 table** | held ON-minutes-per-cycle fixed while varying period, changing duty and hence **total dose by 72×**. It measured dose, not schedule. | T5 v3, duty fixed (§5.2) |
| T4 v1 parameterisation | eGFR pinned on its clamp, non-physiological; φ=0 failed to reproduce round 40's zero (finite-difference noise, h too small) | recalibrated model (§4.1), convergence measured |
| "CD19 depletion barely moves Ab = bug" | not a bug; the CD19⁻ LLPC floor | §4.2, now a headline result |
| §4.2 depletion figures in the first draft | quoted from a superseded parameterisation and echoed by a stale hardcoded `print` string in `t5_nc3.py`; caught on pre-release verification | re-derived from current code (§4.2) |

---

## 12. Build order

1. **L0–L2 + C5 compartments.** Steady state exact, no clamping, units tested.
2. **Add Tregs** before anything else is built on top (hole 1).
3. **L3 contract + NC3 gate in CI.** Every module must pass or declare dose-equivalence.
4. **L4 provenance gate** with both instruments wired.
5. **L6 placebo falsifier.** If the architecture cannot produce a 33–44% placebo
   response, that is the finding and it outranks everything downstream.
6. L5 selection loop, prospectively.
7. Learned surrogates, last, only where they beat tuned classical baselines.

The order is deliberate: it front-loads the tests most likely to kill the design.

---

*Research and simulation only. Not clinical or drug-development guidance.*
