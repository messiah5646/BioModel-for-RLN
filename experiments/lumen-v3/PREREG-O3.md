# PREREG-O3 — should LUMEN's forecast for a new trial depend on whether it is phase 2 or phase 3?

1. **Item / round / date.** O3 (queue, from round 148 research note), round 149, 2026-11-27 UTC (scheduled run).
   Analysis runs in a LATER run, >= 1 h after this file is pushed (rule 12b).

2. **Decision this can change (rule 7).** Whether `lumen.layers.drug_effect` takes trial phase as a covariate on the
   log OR. If P1 PASSES -> phase enters drug_effect as a v3.1 candidate (CLI gets `--phase 2|3`; mean shifts by the
   fitted phase term), labelled Supported (n trials, 12b-timed), pending the W1 prospective record; frozen
   `lumen_forecast.py` and every frozen forecast stay unchanged. If P1 FAILS or is INCONCLUSIVE -> no phase term;
   "phase as a covariate on the relative effect" is added to rule 8's settled list; one row in preprint 3.1.

3. **Rule-8 check (done in this run, before writing this file).** Settled list and project docs searched
   (`grep phase` over experiments/, docs/, lumen/; project_search "phase 2 vs phase 3 covariate"). RT-2 (round 129)
   tested posting year + definition stringency on **control-arm rates** only; r86/r116/RT-1 tested anchor (molecule,
   target) and mechanism class on the effect. **Phase has never been a covariate on the relative effect: untested.**

4. **What has already been seen.** All outcome tables and the CT.gov snapshot are committed and were analysed in
   r86, r116, RT-1 (r127), RT-4 and audits; RT-1 S3 found the time-forward pooled mean beats null (LN 0.327 vs 0.427,
   SLE 0.132 vs 0.190). **No analysis has split effects by phase.** In this run, only the phase *distribution* of the
   trials was counted (no outcomes): LN 20 trials = PHASE2 12, PHASE3 7, PHASE2/PHASE3 1; SLE 40 = PHASE2 26,
   PHASE3 13, PHASE1/PHASE2 1. The analyst knows several public results that suggest phase-2 optimism (e.g.
   anifrolumab MUSE vs TULIP-1), which is why the question was asked. Label: **12b-timed (committed data; earlier
   looks not excludable)**.

5. **Data.** Inputs of `experiments/lumen-bench/timeforward/run_rt1.py`, imported unchanged (tables, `DATE`, `CFG`,
   `late_arms`, `score_stratum`, shipped `iv`/`dl`); phase from `designModule.phases` in the same committed snapshot
   `data/real/ctgov/lupus_with_results.json.gz`. md5s checked by `lumen verify` before the run. No download.
   **Phase coding (fixed):** ph3 = 1 if 'PHASE3' is in the trial's phases list, else 0 (so PHASE2/PHASE3 -> 1,
   PHASE1/PHASE2 -> 0). Any trial with no phase -> excluded from P1 and M1 scoring alike, listed.

6. **Models.** For each held-out arm, using only its training trials (same stratum, as in RT-1):
   - **M1 (reference, shipped; strongest trivial baseline):** fixed-effect IV pooled mean of training late arms
     (weights 1/se^2), var = DL tau^2 + Var(mean). Carry-forward is not applicable (no early read-out in 'all').
   - **M0:** 0.
   - **P1 (only candidate that can ship):** weighted least squares of training late-arm lor on [1, ph3] with the
     **same weights as M1** (1/se^2), so P1 with the phase column dropped reproduces M1 exactly (asserted to 1e-10 in
     the analysis). Mean = a + b * ph3(held-out trial). Var = tau^2 (same shipped DL tau^2 as M1) + x0'(X'WX)^-1 x0.
     **Fallback:** if the training set has < 3 trials in either phase class, P1 = M1 for that unit (flagged).
   - Secondary only: P1-RE (weights 1/(se^2 + tau^2)); cannot ship.

7. **Metric and scheme.** Time-forward (train = trials whose results were first posted strictly earlier), 'all'
   subset, burn-in >= 5 training trials, as in RT-1. **Primary: IV-MSE (`lumen.baselines.iv_mse`) over all scored
   arms, LN + SLE combined.** Verdict via `lumen.gate.verdict(P1, {M0, M1}, reference='M1', margin=0.10)`.

8. **Minimum n.** >= 20 scored trials in which P1 did not fall back, including >= 6 phase-3 and >= 6 phase-2 trials;
   else **UNDERPOWERED** (treated as not shipping, same consequence as INCONCLUSIVE).

9. **Pass rule.** P1 IV-MSE <= 0.90 x M1 and below 1.10 x M0 (gate PASS).

10. **FAIL rule (rule 2b).** P1 IV-MSE >= 1.10 x the best of M0, M1 = **FAIL**. Anything else = INCONCLUSIVE.

11. **Sensitivities / descriptives (cannot turn FAIL into PASS).** (a) LOTO mode; (b) LN-only, SLE-only;
    (c) PHASE2/PHASE3 recoded as phase 2; (d) P1-RE; (e) log score and 95% coverage for M1 and P1 (`lumen.baselines`);
    (f) descriptive: b with its 95% CI from the full-data fit per stratum ("how much bigger are phase-2 effects?");
    (g) IV-MSE on phase-2 held-out arms only and phase-3 only. **Prospective record:** freeze P1 forecasts for the 8
    trials in `prospective/w1-sle/frozen_2026-11-26.json` (phase from the W1 catalogue; no outcomes read),
    descriptive, scored next to M1 when W1 scores; PREREG-W1-SLE unchanged.

12. **Labels.** n trials stated next to any Supported claim; a PASS on committed data is at most Supported
    (12b-timed); n < 10 in any subgroup claim -> Emerging at most.

13. Result JSON through `lumen.provenance.stamp(..., prereg='experiments/lumen-v3/PREREG-O3.md')`.

*Research and simulation only. Not clinical or drug-development guidance.*
