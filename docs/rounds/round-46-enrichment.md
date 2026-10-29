# Round 46 — the squeeze is irreducible by enrollment criteria, and my explanation of it was wrong (2026-10-29)

Pre-registration: `docs/rounds/prereg_round46.md`, written before any test number.
Code: `experiments/rln-identifiability/r46_enrichment.py`. Seeds 2000–2007.

**H1 passes the letter of the pre-registration and fails its substance. H2 resolved
favourably — the window is not empty. H3 FALSIFIED, and that is the round's most useful
result.**

---

## T1 — my own explanation of the squeeze does not survive

Round 44 §3 explained the variability–power squeeze by claiming the jittered parameters
enter **multiplicatively** (`UPCR = p0·A·N`), so spread reshapes rather than translates
the outcome distribution. That explanation was offered *after* seeing the result and was
never tested. H1 pre-declared that additive jitter of matched variance would show a slope
at least 2× shallower.

| mode | CV 0.40 | CV 0.70 | CV 1.00 | slope | relative decline |
|---|---|---|---|---|---|
| multiplicative | 0.535 ± 0.019 | 0.361 ± 0.013 | 0.301 ± 0.012 | −0.390 | **44%** |
| additive | 0.161 ± 0.024 | 0.063 ± 0.005 | 0.053 ± 0.003 | −0.181 | **67%** |

Slope ratio **0.46**, against a pre-declared threshold of 0.50. **H1 passes — by 0.04.**

**But the substance is the opposite of what H1 was testing for.** The squeeze appears under
*both* jitter structures. The additive slope is shallower in absolute terms only because
additive jitter crushes the baseline effect size (d 0.161 vs 0.535); in *relative* terms
the additive case erodes **more** (67% vs 44%).

> **R44 §3's explanation is downgraded from "the mechanism" to "a contributing factor, not
> established."** Multiplicativity is not the cause of the squeeze. The phenomenon
> survives; my account of why does not.

Recorded as a near-miss on a threshold I set myself: had I pre-registered the relative
measure, H1 would have failed outright. The pre-registration is what makes that visible
instead of arguable.

## T2 — the window is NOT empty

Round 43/44 left open whether any variability level gives a realistic placebo rate **and**
a detectable effect simultaneously. Requirement: control CRR in 20–45% **and** diff ≥ 10 pp,
with the seed-SD interval not crossing either boundary.

| CV | control % | diff (pp) | both, strictly? |
|---|---|---|---|
| 0.50 | 15.8 ± 1.75 | +13.5 ± 1.35 | no |
| 0.55 | 18.4 ± 1.71 | +12.5 ± 1.25 | no |
| 0.60 | 20.9 ± 1.49 | +12.3 ± 1.47 | mean only |
| **0.65** | **22.7 ± 1.68** | **+12.2 ± 1.82** | **YES** |
| **0.70** | **24.8 ± 1.67** | **+11.6 ± 1.46** | **YES** |
| 0.75 | 27.0 ± 1.77 | +10.9 ± 1.60 | mean only |

**The model can be simultaneously realistic about placebo response and about treatment
effect.** That closes the structural worry R44 raised.

**Why these diffs exceed R44's.** This round uses a **virtual-twin** design — the same
patient simulated under control and treatment — which removes between-arm sampling noise
and estimates the true average treatment effect. R44 used independent arms, as a real
trial does, and got +9.8 pp at CV 0.70 against +11.6 pp here. **Not a contradiction: a
different and more powerful estimator of a different quantity.** Required-N figures still
use independent arms.

## T3 — enrichment cannot break the squeeze (H3 FALSIFIED)

The actionable question: if heterogeneity erodes effect size, can a trial recover it by
enrolling a narrower population? The test only counts if the rule uses **baseline-measurable**
quantities — enriching on a latent rate constant is not something a trial can do.

N = 4000 at CV 0.70, rules pre-declared, evaluated post-hoc on one draw.

| enrichment rule | eligible | control % | diff pp | Cohen d | gain | req N/arm |
|---|---|---|---|---|---|---|
| baseline UPCR middle tertile | 33% | 22.6% | +11.7 | 0.408 | 1.15× | 231 |
| baseline UPCR top quartile | 25% | 9.5% | +7.6 | 0.397 | 1.12× | 313 |
| UPCR > median & eGFR ≥ 60 | 50% | 13.9% | +8.8 | 0.393 | 1.11× | 302 |
| baseline UPCR top tertile | 33% | 11.8% | +8.2 | 0.392 | 1.11× | 311 |
| baseline UPCR bottom quartile | 25% | 46.3% | +13.3 | 0.390 | 1.10× | 218 |
| **none (all patients)** | **100%** | **26.0%** | **+11.0** | **0.354** | **1.00×** | **274** |

Best observable-only rule: **1.15×**, against a pre-declared 1.25×. **H3 falsified.**

**And it is worse than the table suggests.** Accounting for screening burden — a bigger
effect in a smaller eligible pool is not a better trial — **no enrichment rule beats
enrolling everybody**:

| | must screen per arm |
|---|---|
| none (all patients) | **~274** |
| UPCR > median & eGFR ≥ 60 | ~604 |
| UPCR middle tertile | ~692 |
| UPCR top tertile | ~933 |
| UPCR top quartile | ~1250 |

### The finding underneath, which is the useful one

Baseline UPCR moves the **control response rate** enormously — from **9.5%** in the top
quartile to **46.3%** in the bottom — while leaving Cohen's d essentially flat at
**0.39 ± 0.01** across every stratum.

> **In this model, baseline UPCR is a strongly prognostic and non-predictive marker.** It
> tells you who will do well; it does not tell you who will benefit. That is exactly why
> enrichment cannot break the squeeze: breaking it needs a *predictive* marker, and
> baseline proteinuria is not one.
>
> **The squeeze is irreducible by enrollment criteria defined on baseline renal function.**

A negative result, and more decision-useful than a positive would have been: it says not to
spend a trial's eligibility budget on proteinuria stratification expecting a power gain.

## What this does not establish

- **The eGFR arm of T3 was vacuous.** "Baseline eGFR ≥ 60" captured 100% of simulated
  patients, because the model's baseline eGFR distribution (78.2 with jitter) barely
  straddles 60. So H3 was effectively tested on UPCR alone, and the eGFR rules carry no
  information. A real cohort has far more baseline eGFR spread. **This is a limitation of
  the model's baseline distribution, not a finding about eGFR.**
- Enrichment rules were evaluated on **one draw** of 4000, not across seeds. The d values
  differ by less than their likely sampling error, so the *ranking* among rules should not
  be read as meaningful — only the overall verdict that none reaches 1.25×.
- Synthetic throughout, one model structure, one operating point, a stylised drug whose
  effect is weaker than voclosporin's observed. **An enrichment rule's performance here is
  a hypothesis about this model, not a trial-design recommendation.**
- T1's additive jitter was variance-matched but floored at a small positive value to keep
  rate constants physical; at CV 1.00 that truncation is non-negligible and was not
  characterised.

## Next

1. **T1's falsifier was not the right one.** If multiplicativity is not the cause, find
   what is: the candidate is that `D` is monotone, so damage accumulated before treatment
   is irreversible and heterogeneity in the pre-treatment period partitions patients into
   states the drug can no longer change. Testable by varying when treatment starts.
2. **A predictive marker, if one exists in this model.** T3 says baseline renal function is
   not one. The model's own parameters are candidates — the question is whether any maps to
   something measurable.
3. Widen the baseline eGFR distribution so the eGFR arm of T3 is actually testable.

*Research and simulation only. Not clinical or drug-development guidance.*
