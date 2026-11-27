# PREREG-O2 — can LUMEN give an honest absolute control-arm response rate with a robust MAP prior?

1. **Item / round / date.** O2 (queue, from round 148 research note), round 153, 2026-11-27 UTC (scheduled run).
   Analysis runs in a LATER run, >= 1 h after this file is pushed (rule 12b). This run writes no analysis code that
   touches outcomes and scores nothing.

2. **Decision this can change (rule 7).** Whether LUMEN may print an absolute control-arm response rate (and hence an
   absolute drug-arm rate = expit(logit control + log OR)). If R1 PASSES -> a `control_rate` layer enters `lumen/` as a
   v3.1 candidate (registry status "candidate", CLI `--absolute` prints control rate + 95% interval and the implied
   drug-arm rate), labelled Supported (n trials, 12b-timed), pending the W1 prospective record; LN absolute rates only if
   the LN-only stratum is not a FAIL (LN guard, as in RT-2). If R1 FAILS or is INCONCLUSIVE -> LUMEN stays
   relative-effects-only (`control_rate: "pooled, not modelled"` unchanged) and "the predictive distribution for control
   rates (robust MAP / random-effects)" is added to rule 8's settled list, re-openable only by W1 prospective scoring.

3. **Rule-8 check.** RT-2 (round 129) tested *covariates* (posting year, definition stringency) on control-rate point
   forecasts. O1 (round 151) tested heterogeneity/interval models for the *drug effect* (log OR). No test has scored a
   predictive *distribution* (log score) for control rates, and no MAP prior has been run on any trial. Untested.

4. **What has already been seen (disclose plainly).** The unit table and every control count were analysed in RT-2
   (rounds 128-129) and entered every log OR since round 86. RT-2 already reported, for the same units:
   time-forward (S5, 50 scored trials) point IV-MSE on logit control rate B1 (fixed-effect IV pooled) 0.212,
   B2 (unweighted mean) 0.192 (LN 0.248 vs 0.165; SLE 0.203 vs 0.198); LOTO 95% coverage of the B1 interval 0.88
   (LN 0.85, SLE 0.90). So the analyst knows the unweighted mean is the strongest point baseline and the shipped pooled
   interval under-covers; this is why B2 is given an honest-width interval below and why the pass bar is set against the
   better of B1 and B2. No MAP, random-effects predictive or log score has been computed on these units. In this run
   only rt2_result.json summaries (above) and the W1-SLE frozen file's keys/endpoints were read; no outcome rows.
   Label: **12b-timed (committed data; earlier looks not excludable)**.

5. **Data (no download).** Exactly RT-2's units: `run_rt2.unit_table(lambda st, g: pick(g, 0, CAP[st]))` over
   `experiments/lumen-bench/resp_table.csv` (LN, md5 3ce67bdb51c6597bce6666d914618778) and `sle/sri4_table.csv` (SLE,
   md5 b37dd7cd72c9e912ba079088532c3765), dates from `data/real/ctgov/lupus_with_results.json.gz`
   (md5 1b31cd52ccc62dc6124cd611342ba5d0). `run_rt2.py` executes and writes files at import, so the analysis script
   copies `unit_table`, `pick`, `dl_tau2` verbatim and **asserts** (i) its unit table equals `rt2_coding.csv` on
   stratum/trial/nct/week, and (ii) its time-forward B1 and B2 means equal the `spec == 'S5'` rows of
   `rt2_predictions.csv` to 1e-10. Any assertion failure -> stop, no scoring. Outcome per trial: y = logit((r+0.5)/
   (n-r+0.5)), se2 = 1/(r+0.5) + 1/(n-r+0.5) (RT-2 definitions). One unit per trial.

6. **Models.** For each held-out trial, training set = trials of the same stratum with strictly earlier
   resultsFirstPostDate; scored only if >= 5 training trials (RT-2 S5 scheme). k = training trials.
   - **B1 (shipped "pooled"):** mean = fixed-effect IV mean (w = 1/se2), var = DL tau^2 + 1/sum(w); normal.
   - **B2 (strongest trivial point baseline, given honest width):** mean = unweighted mean of training y,
     var = s^2 (1 + 1/k), s^2 = sample variance (ddof 1) of training y; normal.
   - **R1 (only candidate that can ship): robust MAP** (Schmidli et al. 2014, doi:10.1111/biom.12242).
     Normal-normal hierarchy on logit scale with known se2; flat prior on mu; tau ~ half-normal(scale 1.0).
     Posterior of tau on a fixed grid tau = 0.0005, 0.0015, ..., 3.9995 (4000 midpoints, step 0.001):
     p(tau | y) ∝ HN(tau; 1) * |V_mu(tau)|^{1/2} * prod_i N(y_i; mu_hat(tau), se2_i + tau^2), with
     V_mu = 1/sum(1/(se2_i+tau^2)), mu_hat the RE-weighted mean. MAP predictive for the new trial's true logit:
     mixture over the grid of N(mu_hat(tau), V_mu(tau) + tau^2) with weights p(tau | y). **Robustification:** mixture
     0.8 x MAP + 0.2 x N(m_MAP, 2.0^2), m_MAP = MAP predictive mean (sd 2.0 ≈ unit-information sd for a binary logit).
     Predictive for the observed y: every normal component's variance + se2 of the held-out trial.
   - Secondary only (cannot ship): R0 = non-robust MAP (weight 0); R1-HN0.5 (tau prior scale 0.5); R1-w0.5.
   - Carry-forward is not applicable: the decision is a pre-trial forecast with no read-out of the trial itself.
     (RT-2 S6 found an early control read-out predicts the later one well; reported again here as descriptive D1.)

7. **Metric and scheme.** Time-forward as above, LN + SLE combined. **Primary: mean negative log predictive density
   of y per trial** (one unit per trial; lower is better), mixture densities computed exactly (sum of normal pdfs).
   Secondary: IV-MSE of the predictive mean (`lumen.baselines.iv_mse`), 95% and 80% central-interval coverage
   (mixture quantiles by bisection to 1e-8), LN-only and SLE-only, LOTO mode.

8. **Minimum n.** >= 20 scored trials combined and >= 10 in LN, else **UNDERPOWERED** (same consequence as
   INCONCLUSIVE; LN-only underpowered -> LN absolute rates refused).

9. **Pass rule.** R1's primary score is lower than **the better of B1 and B2** by >= 0.05 nats per trial, the paired
   trial-bootstrap (2000 resamples, seed 20261127) 95% CI of (R1 - that baseline) lies below 0, combined 95% coverage
   is in [0.88, 0.99], and the FAIL rule does not fire.

10. **FAIL rule (rule 2b).** R1 FAILS if its primary score is worse than B1 or B2 by >= 0.05 nats per trial, OR its
    IV-MSE is >= 10% worse than B1's or B2's. (The 0.05-nat margin replaces the template's 10% ratio for a log score,
    as in PREREG-O1.) Anything else is INCONCLUSIVE. **LN guard:** the same FAIL rule applied to LN-only scored trials;
    if it fires, LN absolute rates are refused even if the combined verdict is PASS.

11. **Sensitivities / descriptives (cannot turn FAIL into PASS).** (a) LOTO mode; (b) SLE-only; (c) R0, R1-HN0.5,
    R1-w0.5; (d) 80% coverage; (e) D1: on RT-2 S6's trials with an early (<= 26 wk) control visit, log score of the
    early control logit carried forward (var = its se2 + target se2) vs R1 on the same trials. **Prospective record:**
    freeze R1 control-rate forecasts (mean, 95% interval, as a rate) for the W1-SLE trials whose endpoint is SRI-4,
    trained on all 40 SLE units; BICLA-endpoint trials are listed as not frozen (no BICLA control pool in RT-2).
    Written as `prospective/w1-sle/frozen_O2_<date>.json`, descriptive; PREREG-W1-SLE is not changed.

12. **Labels.** n trials stated next to any Supported claim; a PASS on committed data is at most Supported with the
    12b-timed qualifier and is not called validated before the W1 prospective record; n < 10 subgroup -> Emerging at most.

13. Result JSON through `lumen.provenance.stamp(..., prereg='experiments/lumen-v3/PREREG-O2.md')`; script
    `experiments/lumen-v3/run_o2.py`, committed before it is first executed.

*Research and simulation only. Not clinical or drug-development guidance.*
