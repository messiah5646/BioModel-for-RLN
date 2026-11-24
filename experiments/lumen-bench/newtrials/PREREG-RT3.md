# PREREG RT-20261124-1 run 2 (round 134): early read-out rules scored on NEW trials only

Written 2026-11-24 (scheduled run) from `docs/PREREG-TEMPLATE.md`, on top of the frozen catalogue
(`CATALOGUE-RT3.md`, `candidates_frozen.csv`, commit 836f753, 13:38 UTC). Under rule 12b this run commits and pushes this
file and ENDS. Run 3 (>= 1 h after this commit) downloads every source through `scripts/provenance/fetch.sh`, commits the
manifest, commits the extraction table, and only then scores. No outcome number of any candidate has been extracted or read
by this project before this commit.

## 1. Decision this can change (rule 7)
The CLI headline of `lumen_forecast.py` when `--early` is given (today: pooled prior M1, with the round-132 note that no
early rule is validated).
- **UNDERPOWERED** (< 10 new trials in the primary subset) -> park Q61-type early-read-out work; the CLI keeps pooled + the
  "no early rule validated" note; preprint unchanged.
- **CF-SWITCH** -> `--early` headline becomes raw carry-forward C0 (label at the scored n).
- **UPDATE-SWITCH** -> `--early` headline returns to the shrunk update M2 (label at the scored n).
- **POOLED-WINS** -> keep pooled and drop the early read-out from the CLI output entirely (no early number printed).
- **DEFAULT-FAIL** (pooled M1 loses to a trivial baseline by >= 10% and no switch rule is met) -> the `--early` headline is
  withdrawn: the CLI prints M0, M1, C0 side by side with "no validated rule; pooled lost to <baseline> on new trials".
- **INCONCLUSIVE** -> no CLI change; preprint adds the new-trial numbers.

## 2. What has already been seen (disclosure)
- RT-1 (round 127, 13 arms / 7 trials, 12b-timed on committed data): C0 0.146 < M2 0.259 < M1 0.311; LN without
  voclosporin null best (0.279). LOTO results r86/r112/r116. These are the 60 tabulated trials; none is scored here.
- Candidate trials: titles, registry metadata and years only (run 1). Analyst knows headline primary results of several
  candidates (EMBODY 1/2, EXPLORER, ROSE, BLISS-BELIEVE missed; telitacicept ph2b/ph3 met) but no early-visit by-arm
  counts. So the test is blind to the early-vs-late relationship of the new trials, not to their headline outcomes.

## 3. Data
- Frozen inputs (md5 below): candidates and eligibility rules (`CATALOGUE-RT3.md` section 3, applied unchanged),
  `resp_table.csv` (LN, 20 trials), `sle/sri4_table.csv` (SLE, 40 trials), `sle/bicla_table.csv` (BICLA secondary),
  `data/real/ctgov/lupus_with_results.json.gz` (resultsFirstPostDate of the 60 trials), `lumen_forecast.py`, `run_rt1.py`.
- New data: CT.gov results JSON or open full text / supplement of each candidate, fetched in run 3 via `fetch.sh` into
  `newtrials/sources/` with `MANIFEST.tsv`. Extraction -> `newtrials/new_arms.csv` (columns as resp_table: trial, arm,
  molecule, target, week, r_trt, n_trt, r_ctl, n_ctl, source, page/table), committed BEFORE the scoring commit, with
  `exclusions.csv` (candidate id, rule number, reason).
- Effect = log OR with 0.5 continuity correction; se as `lumen_forecast.load`.

## 4. Visits (as CATALOGUE section 3.5 and PREREG-RT1 section 4)
SLE: early = latest visit 12 <= wk <= 26; target = latest 48 <= wk <= 52. LN: early = latest 12 <= wk <= 26; target = latest
48 <= wk <= 104. An arm with no qualifying early visit is not in the primary subset (enters S3 only).

## 5. Time-forward scheme
Merged set = 60 tabulated trials + included new trials, ordered by date (tabulated: resultsFirstPostDate; new: CATALOGUE
rule 8). Each NEW trial is forecast from all trials (old and new) of its stratum dated strictly earlier; tabulated trials
are never scored. Burn-in: >= 5 training trials in the stratum, else the trial is not scored (counts against minimum n).
Combined = union of within-stratum forecasts (no cross-indication pooling).

## 6. Models (fitted on the training set only; code = `run_rt1.py` functions, unchanged except the data loader)
- M0 null: log OR 0.
- M1 pooled prior (current CLI headline): IV mean of training late arms; v0 = tau2_DL + se(mu)^2.
- M2 shrunk update: M1 updated by the arm's own early log OR, ve = early se^2 + drift_var(training pairs); needs >= 2
  training pairs, else the arm is not scored for M2 and the trial counts against minimum n.
- C0 carry-forward: the arm's own early log OR.
- M3 anchor (secondary only).
Strongest trivial baseline (rule 2b): C0. M0 and M1 are also trivial baselines. Rule 2b applies to whichever model the
decision ships, including the default M1.

## 7. Metric
IV-weighted MSE on the target-visit log OR, weights 1/se_obs^2 (as rounds 86, 116, 127).

## 8. Primary subset and minimum n
New trials with both an early and a target visit, SRI-4 (SLE) or CRR (LN), combined LN+SLE, voclosporin excluded (none
expected). **Minimum n: >= 10 new trials** (and >= 12 arms). Below -> UNDERPOWERED.

## 9. Decision rules (evaluated in this order; first match wins)
Let A = IV-MSE on the primary subset. LN clause: if the LN part of the subset has >= 3 new trials, require
A_LN(shipped model) < 1.10 x A_LN(M0); if < 3 LN trials, the switch ships for SLE only and LN keeps M1 with "not tested in LN".
1. UNDERPOWERED: section 8 not met.
2. CF-SWITCH: A(C0) <= 0.90 x A(M1) AND A(C0) < 1.10 x min(A(M0), A(M2)) AND LN clause for C0.
3. UPDATE-SWITCH: A(M2) <= 0.90 x A(M1) AND A(M2) < 1.10 x min(A(M0), A(C0)) AND LN clause for M2.
   (If 2 and 3 both hold, the lower A ships.)
4. DEFAULT-FAIL: A(M1) >= 1.10 x min(A(M0), A(C0)).
5. POOLED-WINS: A(M1) <= min(A(M0), A(M2), A(C0)).
6. INCONCLUSIVE: anything else.
Losing to any trivial baseline by >= 10% is FAIL for that model, never inconclusive.

## 10. Secondaries (labelled; none changes the decision)
- S1 equal trial weights (mean of per-trial IV-MSE). S2 LN-only and SLE-only primary subsets with n.
- S3 all new trials incl. those without an early visit: M1 vs M0 vs M3 at the late visit (does the pooled prior itself
  generalise to new trials?).
- S4 BICLA-only SLE trials, same four models.
- S5 trial-level bootstrap (2,000 resamples, seed 134) 95% interval of A(C0)/A(M1) and A(M2)/A(M1), descriptive.
- S6 new + RT-1 primary arms pooled (labelled: RT-1 part is 12b-timed only).

## 11. Labels (rule 12)
Supported only if the deciding comparison rests on >= 10 new trials; otherwise Emerging, n trials and arms stated. If every
source fetch_utc in the manifest is later than this commit and the extraction commit precedes the scoring commit, the test
may be called **verifiably pre-registered** (rule 12b), with the section-2 disclosure attached.

## 12. Output
`newtrials/run_rt3.py` -> `rt3_result.json` via `scripts/provenance/stamp.py` (prereg path
`experiments/lumen-bench/newtrials/PREREG-RT3.md`) and `rt3_predictions.csv`.

Research and simulation only. Not clinical or drug-development guidance.

## Input md5 (frozen at PREREG)
- dca7a1578e5f9029cfedf232197a006d  experiments/lumen-bench/newtrials/candidates_frozen.csv
- 8fc695bb2ad2d1f1e4126cdb4db55451  experiments/lumen-bench/newtrials/CATALOGUE-RT3.md
- 3ce67bdb51c6597bce6666d914618778  experiments/lumen-bench/resp_table.csv
- b37dd7cd72c9e912ba079088532c3765  experiments/lumen-bench/sle/sri4_table.csv
- 99eb279e207302973fef6899f0e86722  experiments/lumen-bench/sle/bicla_table.csv
- c8f6274a0fb8b273c107933fff80e1ef  experiments/lumen-bench/lumen_forecast.py
- 2880d4cd9ea28101bcc0d42a7f9ef59b  experiments/lumen-bench/timeforward/run_rt1.py
- 1b31cd52ccc62dc6124cd611342ba5d0  data/real/ctgov/lupus_with_results.json.gz
- 4165bab8f72dee66fd8929e3edff191d  scripts/provenance/fetch.sh
- ca14be70b213b98bbc2b830e2b092329  scripts/provenance/stamp.py
