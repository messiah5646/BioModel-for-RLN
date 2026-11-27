# PREREG-O1 — does a better heterogeneity model improve LUMEN's predictive intervals?

1. **Item / round / date.** O1 (LUMEN v3 build, round 148, live session with Thomas), 2026-11-27 UTC.
   Analysis runs in a LATER run, >= 1 h after this file is pushed (rule 12b).

2. **Decision this can change (rule 7).** Which predictive distribution `lumen.layers.drug_effect` ships for a new
   trial's log OR. If O1a PASSES -> O1a becomes the v3.1 default (mean and 95% interval, P(benefit)); the frozen
   `lumen_forecast.py` and every frozen forecast stay as they are. If O1a FAILS or is INCONCLUSIVE -> M1 stays,
   `lumen/candidates.py` stays unshipped, and the interval question is closed for open data (add it to rule 8).

3. **What has already been seen.** Every table and the CT.gov date snapshot used here are committed and were
   analysed in RT-1 (round 127), RT-4 (round 141) and audits 130/135/140/145, incl. M1 coverage in RT-1 S4.
   No candidate in this file has been run on real trials. One row of `rt1_predictions.csv` was printed while writing
   the analysis plan (ALXN2050 fD_120). Label: **12b-timed (committed data; earlier looks not excludable)**.
   The only new information since is synthetic (`experiments/lumen-v3/sim_o1_output.txt`): with true tau >= 0.3,
   the shipped M1 interval covered 0.84-0.94 against 0.95 nominal; O1a covered 0.94-0.99.

4. **Data.** Exactly the inputs of `experiments/lumen-bench/timeforward/run_rt1.py`, imported unchanged
   (tables, `DATE` snapshot, `CFG`, `late_arms`, `score_stratum`); md5s checked by `lumen verify` before the run.

5. **Models.** For each held-out unit, using only its training trials:
   - **M1 (reference, shipped):** mean = fixed-effect IV mean, var = DL tau^2 + Var(mean); normal.
   - **M0 (strongest trivial for a distributional score):** mean 0, var = M1's var; normal.
   - **O1a (primary candidate):** `lumen.candidates.predictive_frequentist` with REML tau^2, random-effects weights,
     Student-t with df = k - 2 (k = training late arms); for scoring, the t scale is sqrt(sd^2 + se_obs^2).
   - Secondary candidates: O1b Bayes half-normal(0.5) tau prior; O1c Bayes log-normal tau^2 (-2.56, 1.74);
     O1d robust mixture of M1 (w 0.2, vague N(0, 1)). Secondary candidates cannot be shipped by this test.
   Carry-forward is not applicable (no early read-out in the 'all' subset).

6. **Metric and scheme.** Time-forward (train = trials whose results were first posted strictly earlier), strata LN
   and SLE combined, 'all' subset, burn-in >= 5 training trials, as in RT-1. **Primary: mean negative log
   predictive density per trial** (arms averaged within a trial, then averaged over trials; lower is better).
   Secondary: IV-MSE (`lumen.baselines.iv_mse`), 95% coverage, LOTO mode, LN-only and SLE-only.

7. **Minimum n.** >= 20 scored trials, else UNDERPOWERED.

8. **Pass rule.** O1a's primary score is lower than M1's by >= 0.05 nats per trial, the paired trial-bootstrap
   (2000 resamples, seed 20261127) 95% CI of (O1a - M1) lies below 0, combined 95% coverage is in [0.88, 0.99],
   and the FAIL rule does not fire.

9. **FAIL rule (rule 2b).** O1a FAILS if its primary score is worse than M0 or M1 by >= 0.05 nats per trial, OR its
   IV-MSE is >= 10% worse than M0's or M1's. (A 10% ratio is undefined for a log score that can be negative, so the
   distributional margin is fixed at 0.05 nats; this is a disclosed departure from the template's wording, not from
   its intent.) Anything else is INCONCLUSIVE.

10. **Sensitivities (cannot turn FAIL into PASS).** LOTO mode; LN-only; SLE-only; O1b-O1d; coverage at 80%.
    **Prospective record:** the same run freezes O1a and O1b forecasts for the 8 trials in
    `prospective/w1-sle/frozen_2026-11-26.json` (no outcomes read) so they are scored next to M1 when W1 scores;
    descriptive only, PREREG-W1-SLE is not changed.

11. **Labels.** n trials stated next to any Supported claim; a PASS on committed data is at most **Supported** with
    the 12b-timed qualifier, and needs the W1 prospective record before it is called validated.

12. Result JSON through `lumen.provenance.stamp(..., prereg='experiments/lumen-v3/PREREG-O1.md')`.

*Research and simulation only. Not clinical or drug-development guidance.*
