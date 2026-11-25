# PREREG RT-4 run 2 (round 138): early read-outs from regulatory documents and open full texts

Written 2026-11-25 (scheduled run) from `docs/PREREG-TEMPLATE.md`, on top of the frozen catalogue
`experiments/unblock-2026-11-25/CATALOGUE.md` Table A (commit 36dd117, 10:30:31 UTC). This PREREG is pushed >= 1 h after that
commit. Under rule 12b this run commits and pushes this file and ENDS. Run 3 (>= 1 h after this commit) fetches every
source through `scripts/provenance/fetch.sh` into `experiments/lumen-bench/regdocs/sources/` with `MANIFEST.tsv`, commits
the manifest, commits the extraction table, and only then scores. No early-visit number from any Table A source has been
extracted or read by this project before this commit (round 137 recorded locations only).

## 1. Decision this can change (rule 7)
The CLI headline of `lumen_forecast.py` when `--early` is given (today: pooled prior M1, with the round-132 note that no
early rule is validated; round 136 RT-3 was UNDERPOWERED on new trials).
- **UNDERPOWERED** (section 8 not met) -> keep pooled + "no early rule validated" and **stop early-read-out work for good**
  (no further Q61/RT-3/RT-4-type routes are queued; W1 keeps only its late-result freeze).
- **CF-SWITCH** -> `--early` headline becomes raw carry-forward C0 (Emerging at the scored n; see section 11).
- **UPDATE-SWITCH** -> `--early` headline returns to the shrunk update M2 (Emerging at the scored n).
- **POOLED-WINS** -> keep pooled and drop the early read-out from the CLI output entirely (no early number printed).
- **DEFAULT-FAIL** (pooled M1 loses to a trivial baseline by >= 10% and no switch rule is met) -> the `--early` headline is
  withdrawn: the CLI prints M0, M1, C0 side by side with "no validated rule; pooled lost to <baseline>".
- **INCONCLUSIVE** -> no CLI change, early-read-out work stops (as UNDERPOWERED); preprint adds the numbers.

## 2. What has already been seen (disclosure)
- **Late outcomes of every Table A trial except MUSE and LOTUS are already in analysed tables** (r86, r112, r116, RT-1 r127;
  telitacicept in RT-3 r136). MUSE and LOTUS late SRI-4 values were never tabulated, but their headline results are
  known to the analyst (MUSE met its OCS-taper primary; LOTUS stopped for futility). Telitacicept ph2b wk48 counts were
  extracted in RT-3 (round 136, `newtrials/new_arms.csv`).
- RT-1 primary (13 arms / 7 trials, 12b-timed): C0 0.146 < M2 0.259 < M1 0.311 < M0 0.367; LN without voclosporin M0 best
  (0.279). Those 7 trials are re-scored here (they are part of the merged primary); this PREREG therefore cannot be
  blind to about half of its primary subset. Section 10 S1 reports the RT-4-new trials alone for that reason.
- Round 137 opened the source documents only to record table/figure locations. The analyst did not read early-visit values.
- **Label: "12b-timed (late outcomes already seen; early values new)"**. Never "verifiably pre-registered".
- Catalogue correction (titles only, checked 2026-11-25 for this PREREG): MUSE and LOTUS are NOT in `build_sle_table.py`'s
  `EXCLUDED` set as CATALOGUE.md says. They are filtered by its `EXC` regex: MUSE's CT.gov SRI-4 outcomes are titled
  "SRI(4) Response With Oral Corticosteroid ..." (OCS-taper condition, Day 169 and Day 365); LOTUS's are titled "SRI-4
  Composite Response" (Week 24 and Week 52; the word 'composite' trips the regex). TULIP-2 is in `bicla_table.csv` only.

## 3. Data
- Frozen inputs (md5 at the end): CATALOGUE.md Table A, `resp_table.csv` (LN 20 trials), `sle/sri4_table.csv` (SLE 40),
  `sle/bicla_table.csv`, `newtrials/new_arms.csv` (telitacicept late rows), `lupus_with_results.json.gz` (dates and MUSE/LOTUS
  CT.gov outcomes), `lumen_forecast.py`, `run_rt1.py`, `run_rt3.py`.
- New data: the documents in Table A, fetched in run 3 via `fetch.sh`. Source priority per trial = the order in Table A's
  "Open source" column (first listed; second only if the first is unreachable or prints no qualifying value).
- Extraction -> `regdocs/rt4_rows.csv` (columns: trial, nct, stratum, arm, molecule, target, week, r_drug, n_drug, r_ctrl,
  n_ctrl, kind {early, late}, source, page/table/figure, how {table, text, figure-label}) and `regdocs/rt4_exclusions.csv`
  (trial, rule, reason). Both committed BEFORE the scoring commit.
- Effect = log OR with 0.5 continuity correction; se as `lumen_forecast.load`.

## 4. Inclusion list and per-trial rules (frozen)
| Trial | Stratum | What run 3 extracts | Rows added |
|---|---|---|---|
| ACCESS (NCT00774852) | LN | wk52 CRR by arm (ITT) | late wk52 (wk24 row already in table becomes its early visit). Included only if overall by-arm wk52 CR counts with N are printed or are a sum of printed strata; else excluded. |
| BLISS-LN (NCT01639339) | LN | early CRR (not PERR) | early |
| BLISS-SC (NCT01484496) | SLE | early SRI-4 | early |
| BLISS-NEA (NCT01345253, table name belimumab-NEA) | SLE | early SRI-4 | early |
| BLISS-76 (NCT00410384) | SLE | early SRI-4, both dose arms | early |
| PLUTO (NCT01649765) | SLE | early SRI-4 | early |
| TULIP-2 (NCT02446899) | SLE | early BICLA (S4 only); SRI-4 at both visits if printed | early; primary only if SRI-4 at both visits is printed |
| MUSE (NCT01438489) | SLE | plain SRI-4 (no OCS-taper condition) at an early and a target visit, both dose arms | early + late; included only if BOTH visits are printed without the OCS condition; else excluded "definition" |
| LOTUS (NCT03517722) | SLE | SRI-4 at wk24 and wk52 | early + late; included only if the methods text defines its SRI-4 as the standard SRI-4 (SLEDAI-2K >= 4 reduction, no new BILAG A / <= 1 new B, PGA < 0.3 worsening, standard treatment-failure rules) with no added OCS condition; judged from the methods section before reading any results table |
| telitacicept ph2b (NCT02885610) | SLE | early SRI-4, three dose arms | early (late wk48 rows from RT-3 `new_arms.csv`) |
No other trial may be added in run 3 (cross-check sources for AURORA/AURA-LV/NOBILITY/TULIP-1/CALIBRATE are not used).

Extraction rules (RT-3 CATALOGUE section 3 rules 4-7, applied unchanged, plus):
- **Visit:** latest printed visit with 12 <= week <= 26 (months x 52/12, nearest reported visit); target windows as section 5.
- **Endpoint match:** the early value must use the same endpoint definition and analysis population as the trial's
  existing target-visit row (LN: CRR as in `build_resp_table.py`; SLE: SRI-4 responders, primary analysis, NRI as reported).
- **Figure values:** used only if the number is printed as a label in the figure and is unambiguously attributable to one
  arm and visit. Percentages are converted to counts with that visit's printed N, else the arm's randomised N; rounded
  half-up. Unlabelled bars/lines are NOT digitised; the trial is excluded "figure not labelled".
- **Population check:** where the same source also prints the target-visit value for a trial whose target row already
  exists, a difference of > 2 responders in any arm versus the table flags a population mismatch; that trial is excluded
  from the primary (it stays in S3).
- **Arms:** each active dose arm vs shared control, arm names mapped to the existing table labels by dose.
- **Dates:** CT.gov `resultsFirstPostDate` (MUSE 2016-08-15, LOTUS 2022-03-09); telitacicept 2023-12-21 (RT-3).

## 5. Visits and scoring scheme (as RT-1, `run_rt1.CFG`, unchanged)
SLE: early 12-26 wk, target latest 48-52 wk. LN: early <= 26 wk, target latest 48-104 wk; early week < target week.
Merged tables = 60 tabulated trials + RT-4 rows (+ telitacicept rows), each trial dated as section 4. Every trial is
forecast time-forward from all trials of its stratum dated strictly earlier (`run_rt1.score_stratum` mode 'tf', burn-in
>= 5 training trials). Adding rows changes training sets for later trials; that is intended and disclosed.

## 6. Models (`run_rt1` functions, unchanged except the loader appending section-4 rows)
- M0 null: log OR 0. M1 pooled prior (current CLI headline). M2 shrunk update (needs >= 2 training drift pairs, else the
  arm is not M2-scored and is outside the primary). C0 carry-forward (the arm's own early log OR). M3 anchor (secondary).
Strongest trivial baseline (rule 2b): C0; M0 and M1 are also trivial baselines. Rule 2b applies to whichever model the
decision ships, including the default M1.

## 7. Metric
IV-weighted MSE on the target-visit log OR, weights 1/se_obs^2 (as rounds 86, 116, 127, 136).

## 8. Primary subset and minimum n
Early-subset arms (both visits) of the merged tables, M2-scored, SRI-4 (SLE) or CRR (LN), LN+SLE combined, **voclosporin
excluded** (AURA-LV, AURORA). **Minimum n: >= 12 trials and >= 18 arms, of which >= 5 trials carry an early value
extracted in RT-4** (so the result is not an RT-1 re-run, rule 8). Below either -> UNDERPOWERED.

## 9. Decision rules (evaluated in this order; first match wins)
Let A = IV-MSE on the primary subset. LN clause: if the LN part has >= 3 trials, require A_LN(shipped model) <
1.10 x A_LN(M0); if < 3 LN trials, a switch ships for SLE only and LN keeps M1 with "not tested in LN".
1. UNDERPOWERED: section 8 not met.
2. CF-SWITCH: A(C0) <= 0.90 x A(M1) AND A(C0) < 1.10 x min(A(M0), A(M2)) AND LN clause for C0.
3. UPDATE-SWITCH: A(M2) <= 0.90 x A(M1) AND A(M2) < 1.10 x min(A(M0), A(C0)) AND LN clause for M2.
   (If 2 and 3 both hold, the lower A ships.)
4. DEFAULT-FAIL: A(M1) >= 1.10 x min(A(M0), A(C0)).
5. POOLED-WINS: A(M1) <= min(A(M0), A(M2), A(C0)).
6. INCONCLUSIVE: anything else.
Losing to any trivial baseline by >= 10% is FAIL for that model, never inconclusive.

## 10. Secondaries (labelled; none changes the decision)
- S1 RT-4-new trials only (trials whose early value is RT-4-extracted): same four models, n stated. If S1 and the primary
  point to different rules, the round doc says so in its first lines.
- S2 LN-only and SLE-only primary subsets with n. S2b voclosporin included.
- S3 all merged trials incl. excluded-from-primary ones, M0 vs M1 vs M3 at the late visit.
- S4 BICLA SLE trials (bicla_table + TULIP-2 early), same four models.
- S5 trial-level bootstrap (2,000 resamples, seed 138) 95% interval of A(C0)/A(M1) and A(M2)/A(M1), descriptive.
- S6 equal trial weights (mean of per-trial IV-MSE).
- S7 primary without trials whose early value came from a figure label (tables/text only).

## 11. Labels (rule 12)
Emerging at the scored n (trials and arms stated) whatever the outcome: about half the primary was seen in RT-1 and all
late outcomes were already analysed. Supported only if S1 alone (RT-4-new trials) has >= 10 trials and points to the same
rule as the primary. Label text: "12b-timed (late outcomes already seen; early values new)".

## 12. Output
`regdocs/run_rt4.py` -> `rt4_result.json` via `scripts/provenance/stamp.py` (prereg path
`experiments/lumen-bench/regdocs/PREREG-RT4.md`) and `rt4_predictions.csv`.

Research and simulation only. Not clinical or drug-development guidance.

## Input md5 (frozen at PREREG)
- df17c984754ba49a97d8469c4c273d24  experiments/unblock-2026-11-25/CATALOGUE.md
- 3ce67bdb51c6597bce6666d914618778  experiments/lumen-bench/resp_table.csv
- b37dd7cd72c9e912ba079088532c3765  experiments/lumen-bench/sle/sri4_table.csv
- 99eb279e207302973fef6899f0e86722  experiments/lumen-bench/sle/bicla_table.csv
- 7b5654cc67e22d37e352a98bd86b602f  experiments/lumen-bench/newtrials/new_arms.csv
- c8f6274a0fb8b273c107933fff80e1ef  experiments/lumen-bench/lumen_forecast.py
- 2880d4cd9ea28101bcc0d42a7f9ef59b  experiments/lumen-bench/timeforward/run_rt1.py
- e99d081ed5cfe9641ce23c23226703e9  experiments/lumen-bench/newtrials/run_rt3.py
- ea42d9c7d6b37934fb5dc8ef74eb310d  experiments/lumen-bench/sle/build_sle_table.py
- 1b31cd52ccc62dc6124cd611342ba5d0  data/real/ctgov/lupus_with_results.json.gz
- 4165bab8f72dee66fd8929e3edff191d  scripts/provenance/fetch.sh
- ca14be70b213b98bbc2b830e2b092329  scripts/provenance/stamp.py
