# PREREG-REC1 — prospective freeze of three recalibration competitors on W1 (printed-default decision)

1. **Item / round / date.** REC1 (PRIORITY, cross-field items, claude/RESEARCH-cross-field-methods-2026-11-28.md),
   round 164, 2026-11-28 UTC, scheduled run. This run commits and pushes this file and ENDS (rule 12b). The freeze
   script (`experiments/lumen-v3/freeze_rec1.py`), the CT.gov enrollment download and the write-once freeze happen in a
   LATER run, >= 1 h after this push. Scoring happens only in W1 weekly runs, when W1 / W1-SLE trials report.

2. **Decision this can change (rule 7).** What `lumen/layers/drug_effect.py` prints as its probability headline.
   Today (round 161, CAL1 FAIL) it prints `base_rate_positive` (B). `lumen_forecast.py` (md5-frozen), M1, and every
   existing frozen file stay unchanged in every branch.
   - **SWITCH:** a competitor X meets the switch rule (s9) at a look with >= 10 scored trials -> the layer prints X's
     probability as `P_positive_recalibrated` (label Emerging, n trials, "prospective, REC1"), B stays printed beside it;
     CAL1 is recorded as "reopened: prospective pass for X (n)".
   - **CLOSE:** all three competitors meet the FAIL rule (s9) -> B stays the print and the recalibration line is closed
     (rule 8 gains "REC1 recalibration competitors: settled FAIL, n").
   - **KEEP-OPEN (INCONCLUSIVE / UNDERPOWERED):** B stays the print; nothing changes; the line stays open for more W1
     outcomes or the rolling freeze POW1 may add. (This refines the queue text, which offered only switch/close: at
     W1's size most looks are expected to be inconclusive, and "inconclusive" must not be read as "close".)
   Rule 8: this is NOT a retrospective re-test of CAL1/O1. No competitor is scored on the 50 historical trials; they are
   used only to fit one parameter of competitor (a). No retrospective Brier, gap or coverage of any competitor is
   computed or reported, in this run or the freeze run.

3. **What has already been seen.** The 50 historical trials and their CAL1 result (round 161: P1 mean 0.799 vs observed
   0.688, gap +0.110; printed latent P_benefit 0.977) are known; the fitted value of (a) will therefore not be blind.
   That is acceptable because (a) is scored only on W1 outcomes, none of which exists: round 162 (2026-11-28) found 0/21
   W1 + W1-SLE trials scorable, CT.gov hasResults 0/21, no topline. The analyst model's training cutoff is June 2026.
   The van Zwet constants were read today from the paper (s5b); nothing else was opened.

4. **Trials and data.** The 21 trials of the existing freezes: 13 LN in `prospective/w1/w1_frozen_2026-11-21.json`
   (PREREG-W1-2026-11-21, endpoint CRR, table `resp_table.csv`) and 8 SLE in `prospective/w1-sle/frozen_2026-11-26.json`
   (PREREG-W1-SLE, SRI-4 or BICLA at the registered visit, tables `sle/sri4_table.csv` / `sle/bicla_table.csv`).
   Freeze run: fetch each trial's CT.gov v2 record via `scripts/provenance/fetch.sh` into
   `prospective/rec1/sources/` (+ `MANIFEST.tsv`). Any trial with hasResults = true or a public topline found by the
   W1 screen procedure is EXCLUDED in the freeze (listed, not frozen). Inputs per trial: enrollment count (actual if
   given, else anticipated), number of arms and number of control arms (`armGroups`, type PLACEBO_COMPARATOR /
   SHAM / NO_INTERVENTION, or the arm named placebo / standard of care).
   Fitting data for (a): `experiments/lumen-v3/cal1_units.csv` (93 time-forward arms, 50 trials; columns M1, v_M1,
   se_obs, E1, P1), md5 checked by the freeze script against the commit 2b11bca version.

5. **Forecasts frozen per trial** (all are P(observed late log OR > 0), the event CAL1 scored). Shared pieces:
   M1, v = mean and variance of `drug_effect.prior(molecule, target, table, anchor=False)` (the shipped pooled
   predictive; for W1-SLE it must equal the frozen M1_pooled prior_logOR / prior_sd to 1e-3, asserted).
   Anticipated sampling se: n_a = enrollment / n_arms (each arm); p0 = median over the table's trials of the control
   rate r_ctrl/n_ctrl (one value per trial, from the arms the core pools: week <= 104, last visit per arm; O2: the
   control rate is not modelled, this is only a nominal value for the anticipated se); p1 = expit(logit(p0) + M1); se_ant^2 = 1/(n_a p1) + 1/(n_a (1-p1)) + 1/(n_a p0) + 1/(n_a (1-p0)).
   P1 = Phi(M1 / sqrt(v + se_ant^2)).
   - **(a) Mixture-weight recalibration (Novartis-style null component; one parameter).** P_a = (1 - w) P1 + w 0.5.
     w is fitted once on cal1_units.csv by maximising the trial-weighted Bernoulli log-likelihood of E1 given
     (1 - w) P1 + w 0.5 (each arm weight 1/n_arms_in_trial), grid w = 0, 0.001, ..., 1. The fitted w is printed in the
     freeze. Descriptive only (not scored for any decision): temperature variant P_T = Phi(Phi^-1(P1) / T), T fitted the
     same way on grid 0.50..5.00 step 0.01.
   - **(b) van Zwet signal-to-noise shrinkage of the pooled mean.** z = M1 / sqrt(v). Prior on SNR: zero-mean normal
     mixture, weights (0.32, 0.31, 0.30, 0.07), SNR sds tau = (0.61, 1.42, 2.16, 5.64) (van Zwet, Schwab & Senn 2021,
     Stat Med, Table 1, 23 551 Cochrane RCTs, doi 10.1002/sim.9173). Posterior of SNR | z: component k has weight
     q_k ∝ w_k N(z; 0, tau_k^2 + 1), mean m_k = z tau_k^2/(tau_k^2 + 1), var s_k^2 = tau_k^2/(tau_k^2 + 1).
     Effect theta = sqrt(v) SNR. P_b = sum_k q_k Phi( sqrt(v) m_k / sqrt(v s_k^2 + se_ant^2) ).
     Stated stretch: the Cochrane prior describes single-trial SNRs; here it is applied to LUMEN's pooled mean, as the
     research doc proposed. This is the first out-of-sample use of the idea for next-trial forecasting known to us.
   - **(c) Equal-weight median ensemble.** P_c = median{P1, 0.5, P_cf}, where P_cf = Phi(M3 / sqrt(v3 + se_ant^2)) from
     `prior(..., anchor=True)` when the anchor found is the same molecule (carry-forward of the molecule's own earlier
     trials), else P_cf = B (the base rate below). No fitted weights.
   - **B (the current print, reference and strongest trivial baseline):** `base_rate_positive` of the trial's table as
     the layer prints it at freeze time (share of past arms with log OR > 0).
   - Also frozen for reference (cannot change the decision): P1 itself, D = the legacy printed latent P_benefit
     Phi(M1/sqrt(v)), C05 = 0.5.
   Each record carries `uncertainty_sources`: P1/(a)/(c) = "between-trial tau + anticipated sampling"; (b) adds
   "Cochrane SNR prior"; B = "none (base rate)".
   Interval forecasts for secondaries (log OR scale): P1 and (c)-member predictive N(M1, v + se_ant^2); (a) the mixture
   (1-w) N(M1, v + se_ant^2) + w N(0, se_ant^2); (b) the posterior mixture plus se_ant^2. 50% and 95% central intervals
   frozen.

6. **Outcome and scoring (W1 runs only).** Outcome extraction exactly as PREREG-W1-2026-11-21 (LN, CRR) and
   PREREG-W1-SLE (SLE, registered endpoint and visit): counts by arm, fetched via fetch.sh, extraction CSV committed
   before scoring. Per active arm vs its own control, e = 1[observed log OR (0.5 cc) > 0]. Unit = trial. Primary
   metric: **Brier** = mean over trials of the within-trial mean of (P - e)^2 (as CAL1). Calibration gap
   g = mean over trials of (P - mean e). Frozen P is scored as frozen (anticipated se); a re-score with realised se is
   descriptive only.

7. **Minimum n.** < 5 scored trials (LN + SLE pooled) = UNDERPOWERED, interim numbers only. The decision look is the
   first W1 run with >= 5 scored trials; SWITCH additionally needs >= 10 scored trials (at 5-9 a qualifying competitor
   is reported as "interim lead", no change). Looks continue as trials report (max 21 + any later REC1 additions); a
   decision, once taken, flips only at >= 10 more scored trials under the same rule.

8. **Paired uncertainty.** Paired trial bootstrap (2000 resamples, numpy default_rng seed 20261128) of
   Brier(X) - Brier(B) and of |g(X)| - |g(B)|, percentile 95% CI, reported for every competitor.

9. **Rules per competitor X in {a, b, c}, against B (rule 2b: B is the strongest trivial baseline; carry-forward is
   a member of (c), not a separate baseline, because no W1 trial has an early read-out).**
   - **FAIL:** Brier(X) >= 1.10 x Brier(B), OR |g(X)| >= |g(B)| + 0.10.
   - **SWITCH-eligible:** Brier(X) <= 0.90 x Brier(B) AND the bootstrap 95% upper bound of Brier(X) - Brier(B) < 0
     AND |g(X)| < |g(B)| AND Brier(X) < Brier(D) AND Brier(X) <= Brier(C05).
   - Otherwise INCONCLUSIVE for X.
   Decision: >= 10 scored trials and >= 1 SWITCH-eligible -> SWITCH to the eligible X with the lowest Brier (ties: fewer
   fitted parameters, order c, b, a). All three FAIL -> CLOSE. Else KEEP-OPEN. Three competitors are compared; this
   multiplicity is stated in every report and is why SWITCH needs the bootstrap bound, not the point estimate.

10. **Secondaries (descriptive; none changes the decision).** Brier skill 1 - Brier(X)/Brier(B); AUC (trial-level, mean
    e > 0.5 as positive; reported only if both classes have >= 2 trials); 50% / 95% interval coverage of observed arm
    log OR; PIT histogram (5 bins); reliability by halves of P; LN-only and SLE-only numbers; P_T (temperature);
    realised-se re-score. Comparisons: 3 primary (a, b, c vs B), the rest descriptive.

11. **Write-once and provenance.** Freeze file `prospective/rec1/frozen_REC1_<date>.json` written through
    `lumen.provenance.stamp(..., prereg='experiments/lumen-v3/PREREG-REC1.md')` and `provenance.write_once`, its md5
    appended to `prospective/FROZEN-MD5.tsv` in the same commit. Order of commits in the freeze run: (1) fetch +
    manifest, (2) `freeze_rec1.py` with `--selftest` (simulated inputs) passing, (3) the real freeze. Any amendment
    that changes a model, metric, rule or trial list restarts the 1-h clock and must precede any W1 trial scoring
    (POW1 may amend s7/s9 this way).

12. **Labels.** A SWITCH is Emerging (n trials, "prospective") until >= 20 scored trials. No Supported claim from REC1
    below 20 scored trials. Every report states n scored trials, n arms and the three-way multiplicity.

*Research and simulation only. Not clinical or drug-development guidance.*
