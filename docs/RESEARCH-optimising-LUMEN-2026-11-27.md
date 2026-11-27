# Where LUMEN can still be improved, and where it cannot (2026-11-27)

**The question.** Thomas asked for the architecture to be built and for research on how to optimise it. This note
ranks the candidate improvements by one test: could the result change what LUMEN outputs (queue rule 7)? It is a
targeted reading of the literature (PubMed, 9 papers), not a systematic review. One synthetic simulation was run; no
candidate has touched real trial outcomes.

## 1. Bigger learned models: no gain expected

- **Established elsewhere:** five foundation models and two other deep-learning models did not beat deliberately simple
  linear baselines at predicting transcriptome changes after perturbation (Ahlmann-Eltze, Huber & Anders,
  *Nat Methods* 2025, [doi:10.1038/s41592-025-02772-6](https://doi.org/10.1038/s41592-025-02772-6)).
- **Here:** Geneformer, scVI and sysVI lost to PCA or composition five times (rounds 92-104). The learned virtual cell
  lost to the population average, 0 of 11 wins (VC-1, round 142; Emerging, n 11).
- **Consequence for the design:** a learned cell-state or world-model module gets a slot (registry entry
  `cell_state.*`), but it can only be retried when open data cross a size threshold (>= 30 longitudinal patients
  under one therapy; W2 watches for it). Tuning the architecture before that is optimising noise.

## 2. The forecaster's weak point is its uncertainty, not its centre

What ships (M1) is a fixed-effect mean with a DerSimonian-Laird between-trial variance and a 1.96 interval.

- With few or unbalanced studies, DerSimonian-Laird underestimates heterogeneity and normal-approximation intervals
  under-cover; t-based adjustments and Bayesian models with a sensible heterogeneity prior do better
  (Seide, Röver & Friede 2019, [doi:10.1186/s12874-018-0618-3](https://doi.org/10.1186/s12874-018-0618-3);
  Partlett & Riley 2017, [doi:10.1002/sim.7140](https://doi.org/10.1002/sim.7140), who also warn that frequentist
  *prediction* intervals can under-cover when heterogeneity is low and study sizes vary).
- Empirical heterogeneity priors for log odds ratios exist: from 14,886 Cochrane meta-analyses (Turner et al. 2015,
  [doi:10.1002/sim.6381](https://doi.org/10.1002/sim.6381)), and for health technology assessment with very few
  studies (Lilienthal et al. 2024, [doi:10.1002/jrsm.1685](https://doi.org/10.1002/jrsm.1685)). Röver et al. 2021
  ([doi:10.1002/jrsm.1475](https://doi.org/10.1002/jrsm.1475)) give guidance on weakly informative choices.
- **Already seen in LUMEN:** pooled control-rate intervals covered 0.88 against 0.90 (round 129; Emerging).
- **Synthetic check (this round, `experiments/lumen-v3/sim_o1_output.txt`).** With a true between-trial SD of 0.3-0.6
  on the log OR scale, the shipped interval covered 0.84-0.94 against 0.95 nominal. REML + random-effects weights +
  a t quantile (O1a) covered 0.94-0.99. Bayesian variants sat in between. At LUMEN's pool sizes (>= 20 trials) the
  synthetic gap is only 1-3 points of coverage, so **the expected gain is modest**. The biggest gap is in the
  time-forward burn-in, where there are few training trials.
- **Action:** PREREG-O1 (commit 4e8975a) tests O1a against M1 and against a null with the same width. The primary
  metric is the log score, because a forecast is a distribution, and the FAIL rule is rule 2b.

## 3. Robust priors for genuinely new mechanisms

A pooled prior shrinks every new drug towards the historical mean. Robust meta-analytic-predictive priors add a vague
mixture component, so a conflicting new result is not over-shrunk (Schmidli et al. 2014,
[doi:10.1111/biom.12242](https://doi.org/10.1111/biom.12242)). In PREREG-O1 this is the secondary O1d. It cannot ship
from that test.

## 4. Absolute response rates (currently refused)

LUMEN outputs odds ratios only, because covariates did not predict control rates (RT-2, round 129). The standard
industry tool for borrowing historical controls is the same robust MAP prior. It does not need covariates: it needs
honest width. **Candidate O2 (gated, not pre-registered):** a robust MAP control-rate prior, scored time-forward on
coverage and log score against the pooled rate. Decision: whether LUMEN may print absolute arm rates. O2 should wait
for O1, because both depend on the same heterogeneity estimate.

## 5. Phase-2 optimism

Phase-2 effects tend to shrink in phase 3 (anifrolumab phase 2 → TULIP-1 is the case in `PRIOR-ART`, arXiv
2006.05795). Two of the eight frozen W1-SLE trials are phase 2. **Candidate O3:** phase as a covariate on the effect.
Before writing a PREREG, check that RT-2's negative result (posting year and definition, on control rates) does not
already cover it under rule 8. RT-2 was about control rates, so O3 looks untested, but verify this first.

## 6. Things that will not help

More mechanistic detail (Tregs, germinal centre) fills real biological gaps but has never improved prediction on
these endpoints: six harnesses lost to 0-1 parameter rules. Anchoring by molecule or target failed (round 116), and
early read-out rules are settled (rule 8). Compute is not a bottleneck: every forecast runs in under a second.

## Ranked plan

| # | Item | Can change | Status |
|---|---|---|---|
| 1 | W1 prospective scoring | whether the pooled prior is validated out of time | weekly, due 2026-11-28 |
| 2 | **O1** heterogeneity / interval model | shipped interval and P(benefit) | **PREREG pushed 4e8975a; analysis >= 1 h later** |
| 3 | O2 robust MAP control rate | whether absolute rates may be printed | gated, after O1 |
| 4 | O3 phase covariate | whether phase-2 forecasts are shrunk | check rule 8 first |
| 5 | learned cell-state retry | whether a cell-state layer exists | blocked by data size (W2) |

*Research and simulation only. Not clinical or drug-development guidance.*
