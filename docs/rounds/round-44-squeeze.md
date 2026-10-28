# Round 44 — the variability–power squeeze is real, and it is not a threshold artifact (2026-10-28)

Pre-registration: `out/PREREG_round44.md`, written before any test number.
Code: `src/t7_squeeze.py`, `src/t7_analyze.py`. Seeds 1000–1007, n=600/arm/cell, 48 cells.

**H1 supported. H2 supported — and H2 was the one that could have killed the finding.
H3 supported.** Round 43 §4 survives a deliberate attempt to destroy it.

**One overclaim was identified and killed before it was made (§5).**

---

## 1. What was attacked

R43 §4 reported, unpredicted and on **one seed per cell**, that as between-patient
variability rises the control response rate rises toward realism while the treated−control
gap shrinks. A single-seed observation is not a result.

**The alternative explanation that had to be ruled out:** CRR is a threshold on a
continuous quantity. As CV → ∞, the fraction of parameter space clearing UPCR ≤ 0.5 is
increasingly set by spread rather than by treatment, so *any* two arms must converge. If
that is all this is, the squeeze is arithmetic, carries no trial-design meaning, and R43 §4
should be retracted.

**Discriminator, declared in advance:** compute a **continuous** standardised effect size
(Cohen's d on log UPCR₅₂) on the same populations. Binary-only erosion ⇒ artifact.
Erosion in the continuous measure too ⇒ the distributions themselves are separating less.

## 2. Result

8 seeds per cell, **paired within seed** (both arms drawn from the same seed). `±` is
**seed SD**, reported separately from binomial CI per `lessons.md` §2.

| CV | control % | treated % | diff (pp) | **Cohen d** | median req N/arm |
|---|---|---|---|---|---|
| 0.40 | 11.0 ± 1.08 | 22.7 ± 1.62 | +11.7 ± 1.57 | **0.488 ± 0.066** | 149 |
| 0.50 | 16.6 ± 0.92 | 28.4 ± 1.95 | +11.8 ± 2.05 | **0.399 ± 0.069** | 213 |
| 0.60 | 21.5 ± 1.16 | 32.2 ± 2.17 | +10.7 ± 2.48 | **0.344 ± 0.072** | 259 |
| 0.70 | 25.7 ± 1.16 | 35.5 ± 2.26 | +9.8 ± 3.02 | **0.306 ± 0.076** | 406 |
| 0.80 | 29.1 ± 1.54 | 38.0 ± 2.63 | +8.8 ± 3.03 | **0.278 ± 0.073** | 484 |
| 1.00 | 34.5 ± 1.87 | 42.2 ± 2.98 | +7.7 ± 3.13 | **0.234 ± 0.076** | 671 |

**H1 — binary squeeze.** diff_pp falls +11.7 → +7.7 pp. Welch t = +3.22, **p = 8.8e-3**.
Supported, though note the seed SD *doubles* across the range (1.57 → 3.13), so the binary
endpoint also becomes noisier — itself part of the story.

**H2 — the discriminator.** Cohen's d falls **0.488 → 0.234**. Welch t = +7.16,
**p = 5.4e-6**; regression of d on CV gives slope **−0.404**, **p = 1.3e-9**, r² = 0.554.

> **The continuous effect size erodes as fast as the binary one. The squeeze is not a
> threshold artifact — the underlying distributions genuinely separate less as
> heterogeneity rises. R43 §4 survives.**

Note the continuous test is *stronger* evidence than the binary one (p = 5e-6 vs 9e-3),
which is the opposite of what the artifact hypothesis predicts.

**H3 — required sample size.** Monotonic: 149 → 213 → 259 → 406 → 484 → 671 per arm,
a **4.5× escalation** across the CV range.

## 3. What it means

Between-patient heterogeneity is not merely noise around a fixed treatment effect; in this
model it **shrinks the effect itself** at the distributional level. The mechanism is that
the jittered parameters enter multiplicatively (`UPCR = p0·A·N`), so spread does not
translate the outcome distribution — it reshapes it, and a fixed multiplicative treatment
effect buys progressively less separation.

Combined with R43: the variability needed to reproduce a realistic lupus placebo arm
(CV ≈ 0.70) is also the variability at which the detectable effect is **37% smaller**
(d 0.488 → 0.306) and required N is **2.7× larger** than at CV 0.40.

## 4. The window, now with error bars

R43 guessed the viable window at CV 0.6–0.8. With 8 seeds: control CRR enters the 20–45%
band at **CV 0.60** (21.5 ± 1.16) rather than 0.70, and the ≥10 pp separation is lost
between **CV 0.60 (+10.7) and CV 0.70 (+9.8)**. So the window is **narrower and lower than
R43 reported — roughly CV 0.60–0.65, and it may be empty**: at CV 0.60 the +10.7 pp mean
carries a seed SD of 2.48, so individual seeds fall below 10.

**This is a partial correction to R43 §2**, which reported CV 0.70 as the value satisfying
both criteria on one seed. With 8 seeds, CV 0.70 gives +9.8 pp and **fails H2 of R43**.

## 5. An overclaim identified and killed

The required-N table invites the headline *"LN trials are systematically underpowered once
realistic heterogeneity is accounted for."* **It does not follow, and it is wrong.**

| | control | treated | diff | required N/arm |
|---|---|---|---|---|
| our model, CV 0.70 | 25.7% | 35.5% | +9.8 pp | 343 |
| **AURORA-1 observed** | ~23% | ~41% | **+18 pp** | **102** |

AURORA-1 randomised ~357 total (~179/arm) — **comfortably powered for the effect it
actually observed.** Our large required-N comes from our *stylised treatment effect being
weaker than voclosporin's real one* (+9.8 vs +18 pp), not from any demonstration about real
trials. The 4.5× escalation **with CV** is real; the absolute numbers are a property of our
drug model and must not be quoted about real trials.

Recorded here because it is exactly the kind of claim this project has had to retract
before — plausible, quotable, and not supported by what was computed.

## 6. What this does not establish

- Synthetic throughout, one model structure, one operating point, one stylised drug.
- The erosion mechanism (multiplicative parameters reshaping the distribution) is an
  explanation offered **after** seeing the result and is **not separately tested**. A
  dedicated test would compare additive vs multiplicative jitter. **Not done.**
- Cohen's d on log UPCR₅₂ is one continuous summary. Others (E-distance, Wasserstein)
  might behave differently; untested.
- 8 seeds give seed SD, not a converged distribution. The CV 0.60 window verdict rests on
  means ± SD, not a formal test of whether any CV satisfies both criteria simultaneously.
- Required-N uses a two-proportion normal approximation with no interim analyses,
  dropout, or covariate adjustment.

## 7. Next

1. **Additive vs multiplicative jitter**, to test §3's proposed mechanism directly. Cheap,
   and it is the obvious falsifier of the explanation rather than the phenomenon.
2. **Formally test whether the R43/R44 window is empty** — the most decision-relevant open
   question, since an empty window would mean the model cannot be simultaneously realistic
   about placebo and about treatment effect.
3. Feed §2–§4 into contribution A; the heterogeneity-vs-power relation is a genuine
   trial-design output.
4. Unchanged: gate A (PMID 42006208, needs a browser), G0.

*Research and simulation only. Not clinical or drug-development guidance.*
