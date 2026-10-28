# Round 41 — the BioNeMo stack: what it can and cannot supply (2026-10-28)

**Nothing here was computed. This is retrieval, so by rule 1 every claim below is a
hypothesis about the literature and the vendor catalogue, not a result.** It is
reconnaissance aimed at one decision: can any publicly deployable, tested model slot
into the §2 architecture, and at which layer.

Headline: **the catalogue confirms round 37's empty category from the production side,
and the one model NVIDIA ships for single-cell work is the one an independent benchmark
ranks last.** The stack is usable, but one layer lower than it first appears.

---

## 1. What is actually in the catalogue

Eleven models on `build.nvidia.com/explore/biology`, plus the framework's four
pretrained checkpoints. Sorted by what they *are*, not by what they are marketed as:

| category | models | temporal axis? | action axis? |
|---|---|---|---|
| structure prediction | OpenFold2, OpenFold3, Boltz-2, MSA-search | no | no |
| protein/binder design | RFdiffusion, RFD4-Proteina, ProteinMPNN | no | no |
| small-molecule | MolMIM, GenMol, DiffDock | no | no |
| sequence foundation | Evo2-40B (DNA), ESM-2, AMPLIFY (protein) | no | no |
| single-cell | Geneformer | no | no |

**Every one is a static predictor or generator.** Not one has a transition operator,
an iterable rollout, or a time index. Round 37 established that action-conditioned
multi-step biology is an empty category across the *preprint* literature; this says the
same of the *deployed, API-callable, commercially supported* literature. The category
is empty on both sides.

**Consequence for §2: none of this can be L2.** The core dynamics layer has no vendor
to buy it from. That is simultaneously the project's risk and its entire opportunity.

## 2. Geneformer fails C2 — third independent triangulation

C2 (no cell foundation model beats simple baselines at perturbation prediction) was
previously supported by Ahlmann-Eltze 2025 and by our own rounds 20/21. A third,
independent source now measures the *specific model NVIDIA ships*.

From "Benchmarking Transcriptomics Foundation Models for Perturbation Analysis: one PCA
still rules them all" (arXiv 2410.13956), Replogle dataset:

| task | best classical | Geneformer | other FMs |
|---|---|---|---|
| linear probing (top-5) | **PCA 0.058**, scVI 0.032 | **0.005** | UCE 0.008, cellPLM 0.009, scGPT 0.011 |
| known-relationship recall (CORUM) | **scVI 0.503**, PCA 0.450 | **0.149** | scGPT 0.389, UCE 0.366, cellPLM 0.284 |
| reconstruction (Spearman) | **scVI 0.265**, PCA 0.215 | — | scGPT 0.180, UCE 0.172, cellPLM 0.169 |
| batch integration (iLISI) | PCA 0.919 | 0.979 | scGPT 0.941 |

**Geneformer is last among the foundation models on CORUM recall — 0.149 against scVI's
0.503, a factor of 3.4 — and 11× worse than PCA at linear probing.** The single task it
wins, batch integration, is the one where **random embeddings score 1.000**, i.e. the
metric rewards destroying biological signal. Fine-tuning scGPT improved batch
integration while cutting CORUM recall 0.389 → 0.166, the same trade.

**This also partly discharges a standing debt.** Round 37 worried that our round-20/21
null "may have been measured with the wrong instrument." Here the null holds across six
metric families including recall and reconstruction, not just MSE. That is weak evidence
the instrument was not the problem — it does **not** remove the need to re-score
rounds 20/21 with WMSE/NIR/DRF, because this is a different corpus and a different task.

**Decision: Geneformer is barred from the L8 learned-surrogate slot** unless it beats a
tuned PCA/scVI baseline on our own paired comparison. On present evidence it will not.

## 3. Boltz-2 affinity: the marketing and the independent benchmarks disagree

Boltz-2 is the most tempting component, because §5.5's modules carry affinity-type
parameters and `lessons.md` already tracks binder K_d (obe-cel CAT19 116 nM vs FMC63
0.9 nM). Vendor framing is "unprecedented accuracy." Independent results collected at
rowansci.com/blog/boltz2-benchmarks:

- **PL-REX:** mean Pearson **~0.42**. Placed *second* to SQM 2.20, and only **5–7%**
  better than ΔvinaRF20 and classical docking scores (GlideSP, Gold ChemPLP), while
  running an order of magnitude slower than Vina or Glide.
- **Molecular glues (93 compounds vs OpenFE FEP):** "dramatically underperformed," with
  "generally poor or even negative correlations."
- **Uni-FEP (~350 proteins, ~5,800 ligands):** significantly lags FEP where buried water
  matters. **Affinity compression** — on 75 of 350 targets predictions fell within
  ~2 kcal/mol, including 21 targets whose experimental values span >4 kcal/mol.
- **TNKS2/3CLPro:** run-to-run variability up to **1.5 kcal/mol**, enough to obscure
  enrichment; authors judge it unsuitable to replace physics-based methods late-stage.
- **ASAP-Polaris antiviral challenge:** worst mean absolute error of the methods studied
  (comparators were fine-tuned).
- **"No simple baselines (such as a mean predictor) are reported"** — C2's exact
  complaint, in a different subfield.

**Affinity compression is the same pathology as control bias / signal dilution in the
perturbation benchmarks**: the model regresses toward the mean and the metric hides it.
That a 5–7% margin over ΔvinaRF20 is presented as a generational advance is the
baseline failure this project has hit repeatedly (R18, R24).

Caveat: these predate the June 2026 Boltz-2.1 release, which changed the affinity head.
Any use must re-benchmark against 2.1 **with a mean-predictor baseline included.**

## 4. What the stack *can* supply — and the gap that matters most

The usable role is **grounding intervention-module parameters at L3**, not supplying L2.
A concrete loop for the `licar_t` module of §5.5:

```
RFdiffusion / RFD4-Proteina  ->  candidate split-CAR or photoswitch backbone
ProteinMPNN                  ->  sequences for that backbone
OpenFold3 / Boltz-2          ->  structure, and a (weak) affinity estimate
                             ->  constrains module parameters currently 'assumed'
```

This matters because §5.5 lists `I_sat`, `leak` and `kill` as **assumed**, and §6's
provenance gate refuses mechanistic claims resting on assumed parameters. Moving any of
them to *derived* widens what the system is allowed to say.

**But the gap is sharp, and it is the one that matters.** §5.2 proved the NC3-critical
parameter is the photoswitch **off-rate** `k_off`: schedule sensitivity scaled 4.18× →
9.18× as off-rate went 1 → 24 min, and that is the entire ~9× dose-sparing lever.

> **Boltz-2 predicts binding affinity, a thermodynamic ratio. It does not predict
> `k_on` and `k_off` separately.** Two photoswitches with identical K_d can have off-rates
> differing by an order of magnitude — iLID (~1 min) and CRY2 L348F (~24 min) are exactly
> that case — and §5.2 shows they behave completely differently under the same light
> dose.

**So the BioNeMo stack cannot constrain the parameter this architecture most needs.**
Kinetics must come from the photoswitch literature (the toolbox is already compiled:
iLID, CpLID, CRY2 and its L348F/W349R variants, Magnets, BphP1, REDMAP/REDLIP) or from
dedicated measurement. No deployable model supplies it. **[HYPOTHESIS]** — I did not
find a kinetics-prediction model; absence of evidence here is weak, and a targeted search
for photoswitch-kinetics predictors is worth one round.

## 5. What changes in the system design

| § | change |
|---|---|
| §8 | Geneformer named and **barred** from the learned-surrogate slot pending a paired win over tuned PCA/scVI. |
| §5.5 | Structure/design stack named as the route to move `kill` from *assumed* to *derived*; Boltz-2 admitted **only** with a mean-predictor baseline and 2.1 re-benchmark. |
| §5.2 | Reinforced: `k_off` is NC3-critical **and** unavailable from any deployable model. It is a measurement requirement, not a modelling choice. |
| §10 | New hole: no source for photoswitch kinetic parameters. |
| §1/C3 | C3 strengthened — the empty category is confirmed on the deployed side, not only the preprint side. |

## 6. What this does not establish

- Retrieval only; nothing computed. No number above was reproduced by us.
- The PCA benchmark is one corpus (Replogle/L1000). It does **not** transfer
  automatically to our corpora, and does not discharge the WMSE/NIR/DRF re-scoring debt.
- Boltz-2 results predate 2.1. The criticism may be partly stale; the *baseline-reporting*
  criticism is not, since it is about method, not version.
- I did not test any BioNeMo endpoint. Reachability, licensing and cost are unexamined.
- Absence of a kinetics predictor is a failed search, not a proof of absence.

## 7. Next, ranked by value per unit cost

1. **One targeted search for photoswitch kinetics prediction.** Cheap; would close the
   §4 gap or confirm it as a measurement requirement.
2. **Re-score rounds 20/21 with WMSE/NIR/DRF.** Still the standing debt; §2 lowers the
   prior that the instrument was the problem but does not settle it.
3. Test §5.4's ~9× schedule prediction against the measured dose-response in
   Nguyen/Zhou 2021. Still the cheapest decisive falsifier in the design.
4. If any BioNeMo endpoint is used, benchmark it against a mean predictor **first**.

*Research and simulation only. Not clinical or drug-development guidance.*
