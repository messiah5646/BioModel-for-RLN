# PREREG-CAL1 — is the printed P_benefit calibrated? (display decision)

1. **Item / round / date.** RT-20261127-1 (reviewer route, claude/REVIEW-2026-11-27.md Route 1), round 159,
   2026-11-27 UTC, scheduled run. This run commits and pushes this file and ENDS (rule 12b). The analysis script
   (`experiments/lumen-v3/run_cal1.py`) is written, committed and run in a LATER run, >= 1 h after this push.

2. **Decision this can change (rule 7).** What the LUMEN display layer prints next to the pooled forecast.
   `lumen_forecast.py` (md5-frozen) and every frozen forecast file (write-once) are NOT changed in any branch;
   only `lumen/layers/drug_effect.py` output fields, the CLI text, where-this-stands and the preprint can change.
   - **DROP (FAIL branch):** the layer stops printing `P_benefit` as a headline; it prints instead
     `base_rate_positive` = the empirical fraction of past trial arms with observed log OR > 0 (with n), labelled
     "historical base rate, not a trial-specific probability"; the preprint says P_benefit was mis-calibrated.
   - **KEEP (PASS branch):** `P_benefit` stays, printed with its observed calibration (mean predicted vs observed
     fraction, n trials).
   - **INCONCLUSIVE:** `P_benefit` stays, with the observed calibration and "calibration not established (n)" printed.
   Rule 8: the predictive model (M1) is not re-opened; no model, interval or tau estimator changes in any branch.

3. **What has already been seen.** All inputs are committed and were analysed before: RT-1 (round 127; IV-MSE and
   95% coverage of M1 in `rt1_result.json` S4, combined coverage 0.89), O1/O2/O3 (rounds 151-154), audits
   130-155 and the reviewer re-run of 2026-11-27. The fraction of trials with observed log OR > 0 and the fraction
   significant are not computed in any committed script or result JSON (grep of experiments/ and lumen/ on 2026-11-27 for lor>0 / obs>0 / 1.96-se thresholds found none), but the tables are
   public and were read, so earlier looks cannot be excluded. While writing this file only the shape of
   `rt1_predictions.csv` (93 time-forward 'all' arm rows; 50 trials: 15 LN, 35 SLE) and its column names were printed;
   no outcome column was printed. Label: **12b-timed (committed data; earlier looks not excludable)**.

4. **Data.** Exactly the inputs of `experiments/lumen-bench/timeforward/run_rt1.py`, imported unchanged
   (`LN`, `SLE`, `DATE`, `CFG`, `late_arms`, `score_stratum`), md5s checked by `lumen verify` (provenance.verify) first.
   Units = the time-forward 'all' subset rows of `score_stratum(.., 'tf')` (burn-in >= 5 training trials, strata LN and
   SLE, the stratum's own training trials). The script asserts that M1, v_M1, obs and se_obs equal the committed
   `rt1_predictions.csv` (tf, all) to 1e-9. No download. BICLA (`sle/bicla_table.csv`) is used only in D7.

5. **Predictions and trivial baselines** (all per held-out arm, using only that unit's training trials):
   - **Event E1:** observed late log OR > 0 (the table row RT-1 scores: 0.5-corrected counts).
   - **P1 (the model, primary):** P(observed log OR > 0) under the shipped predictive incl. the arm's own sampling
     error: Phi( M1 / sqrt(v_M1 + se_obs^2) ). This is the observable counterpart of the printed P_benefit.
   - **H (strongest trivial baseline):** time-forward constant historical fraction = share of the training trials'
     late arms (`late_arms(train, stratum)`) with observed log OR > 0. No carry-forward exists ('all' subset).
   - **C05:** constant 0.5.
   - **D-latent (descriptive only):** the printed number itself, Phi( M1 / sqrt(v_M1) ), scored against E1.

6. **Metric and scheme.** Time-forward, LN+SLE combined. Unit = trial: arm scores are averaged within a trial, then
   over trials (as O1). **Primary metric: calibration-in-the-large gap** g = mean over trials of
   (mean arm P1 - mean arm E1). Co-primary for the FAIL rule: **Brier score** (mean over trials of the within-trial
   mean of (p - E1)^2) for P1, H and C05. Paired trial bootstrap, 2000 resamples, numpy default_rng seed 20261127,
   percentile 95% CI for g and for Brier(P1) - Brier(H).

7. **Minimum n.** >= 20 scored trials, else UNDERPOWERED (P_benefit kept, printed "calibration not scored").

8. **PASS (KEEP) rule.** |g| < 0.10 AND the bootstrap 95% CI of g lies inside (-0.10, 0.10) AND
   Brier(P1) <= Brier(H) AND Brier(P1) <= Brier(C05), and the FAIL rule does not fire.

9. **FAIL (DROP) rule, rule 2b.** |g| >= 0.10 (point estimate) OR Brier(P1) >= 1.10 x Brier(H) OR
   Brier(P1) >= 1.10 x Brier(C05). The strongest trivial baseline is H. Anything else = **INCONCLUSIVE**
   (expected to be likely at n 50: the CI half-width of g may exceed 0.10).

10. **Secondary (b): P(significant benefit).** Event E2: obs / se_obs > 1.96 (two-sided alpha 0.05 Wald test on the
    table row, favourable direction). Note: this is the table's late endpoint, not necessarily the sponsor's
    registered primary analysis. Model P2 = Phi( (M1 - 1.96 se_obs) / sqrt(v_M1 + se_obs^2) ); baseline H2 = share
    of training late arms with lor/se > 1.96; C05. Same gap/Brier/bootstrap and the same three-branch rules.
    Decision: (b) PASS -> the layer may add `P_sig_benefit` labelled Emerging with n; FAIL or INCONCLUSIVE -> nothing
    added. (b) cannot change the (a) verdict. Sensitivity (b'): either-direction significance |obs/se| > 1.96.

11. **Descriptives (cannot change any verdict).**
    D1 D-latent vs E1 (gap, Brier). D2 reliability table: trials split into tertiles of mean-arm P1, mean predicted
    vs observed fraction per tertile, plus the same for H. D3 LN-only and SLE-only gaps and Brier. D4 LOTO mode
    (`score_stratum(.., 'loto')`, H from all other trials of the stratum). D5 arm-level units (each arm weight 1).
    D6 95% interval coverage of M1 (repeat of RT-1 S4 / O1; settled, reported only). D7 BICLA: same time-forward
    scheme and rules on `bicla_table.csv` (CFG 'SLE'), descriptive; the BICLA P_benefit display follows the (a)
    verdict and is printed "calibration scored on LN CRR + SLE SRI-4 only".
    Molecule note (round 139 style): first line of the round doc states how many of the 50 trials are belimumab
    and anifrolumab.

12. **Labels.** Any Supported claim states n trials and "12b-timed". PASS on committed data = at most Supported
    and is re-checked when W1 / W1-SLE prospective trials score. n < 10 in any stratum -> Emerging for that stratum.

13. Result JSON `experiments/lumen-v3/cal1_result.json` and `cal1_units.csv` written through
    `lumen.provenance.stamp(..., prereg='experiments/lumen-v3/PREREG-CAL1.md')` and `provenance.write_once`.
    `run_cal1.py --selftest` (simulated units only) is run before the real run; script committed before the real run.

*Research and simulation only. Not clinical or drug-development guidance.*
