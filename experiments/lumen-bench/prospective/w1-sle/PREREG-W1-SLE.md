# PREREG-W1-SLE (RT-20261126-1 run 2, round 147, 2026-11-26 UTC): prospective test of the shipped pooled prior on unreported non-renal SLE trials

Written from `docs/PREREG-TEMPLATE.md` on top of the run-1 catalogue (`CATALOGUE-W1SLE.md`, commit 4ed6f82, pushed 16:27 UTC;
this file is written >= 1 h later, rule 12b). The frozen forecasts (`frozen_2026-11-26.json`) follow as a SEPARATE, LATER
commit in this same run (route spec). Scoring happens only in later W1 runs, when results appear.

## 1. Item / round / date
RT-20261126-1 REVIEWER ROUTE, run 2; round 147; 2026-11-26.

## 2. Decision this can change (rule 7)
The pooled-prior headline of `lumen_forecast.py` (CLI default since round 116) and the preprint's pooled-prior claim.
- **PASS** (pooled M1 beats null M0 by >= 10% on >= 5 prospectively scored trials) -> the headline gets its first
  prospective support; label Emerging until >= 10 trials scored, then Supported (n stated).
- **FAIL** (M1 loses to M0 by >= 10%) -> the CLI headline becomes null-centred intervals (OR 1, width = the pooled
  prior's), and the preprint drops the pooled-prior claim.
- **INCONCLUSIVE** -> no change; numbers reported with n.
- **UNDERPOWERED** (< 5 scored trials) -> interim numbers only, no decision.

## 3. What has already been seen (disclosure)
- The 40-trial SRI-4 and 26-trial BICLA SLE tables (`sle/`, round 116 P3) from which the pooled prior is built; RT-3/RT-4
  trials. None of the candidates below is in them (catalogue rule 2).
- For the candidates: registry metadata only (run 1), plus web-search result TITLES (run 1 and today); no outcome number
  of any candidate has been read. The analyst model's training cutoff is June 2026; none of the 8 included trials had
  reached primary completion by then (earliest: frexalimab 2026-09-04), so no topline could have been in training data.
- Run-2 re-screen (`rescreen_run2.csv`): CT.gov hasResults = false for all 10; **sirolimus NCT04582136 EXCLUDED**
  (a German news headline on mTOR inhibition improving lupus outcomes = possible public topline, unresolved) and
  **MiSLE NCT02633163 EXCLUDED** (Lupus Foundation news "MSC therapy shows promise" + an ACR abstract surfaced, may be
  MiSLE, unresolved). Both excluded conservatively without opening the pages.

## 4. Data / trials (frozen inclusion list: 8 trials, 6 molecules)
| NCT | trial | molecule | target (for M3 only) | endpoint | primary visit | table for prior |
|---|---|---|---|---|---|---|
| NCT05039840 | frexalimab APATURA ph2 | frexalimab | CD40L | SRI-4 | wk24 | sri4 |
| NCT06559163 | obexelimab SunStone ph2 | obexelimab | CD19/FcgRIIb | BICLA | wk24 | bicla |
| NCT05672576 | cenerimod OPUS-1 ph3 | cenerimod | S1P1 | SRI-4 | wk52 | sri4 |
| NCT05648500 | cenerimod OPUS-2 ph3 | cenerimod | S1P1 | SRI-4 | wk52 | sri4 |
| NCT05843643 | upadacitinib SELECT-SLE ph3 (2 studies) | upadacitinib | JAK | BICLA | wk52 | bicla |
| NCT05639114 | ianalumab SIRIUS-SLE 1 ph3 | ianalumab | BAFF-R | SRI-4 | wk60 | sri4 |
| NCT05624749 | ianalumab SIRIUS-SLE 2 ph3 | ianalumab | BAFF-R | SRI-4 | wk60 | sri4 |
| NCT06015230 | GR1603 ph1b/2 (phase 2 part) | GR1603 | unknown | BICLA | wk24 | bicla |
Exact target strings only (BAFF-R is not BAFF), as W1-LN. A trial whose topline becomes public between this push and
the freeze commit is excluded in the freeze commit (none expected within minutes).
Outcome data (later W1 runs): CT.gov posted results JSON, or a peer-reviewed paper / congress abstract / regulatory
document giving responders and n by arm, fetched via `scripts/provenance/fetch.sh` into `prospective/w1-sle/sources/`
with `MANIFEST.tsv`; extraction -> `w1sle_arms.csv` committed BEFORE the scoring commit. A press release giving % and n by
arm may be used (counts = round(% x n)), flagged `source=press`; S6 drops them.

## 5. Models compared (all frozen now, `frozen_2026-11-26.json`)
- **M1 pooled prior (the shipped CLI headline, under test):** `lumen_forecast.forecast(molecule, target, R=load(<table>),
  anchor=False)` with `lumen_forecast.py` at HEAD (md5 below), table = `sle/sri4_table.csv` for SRI-4 trials and
  `sle/bicla_table.csv` for BICLA trials. Output: log OR mean m and sd -> 95% interval m +/- 1.96 sd.
- **M0 null (trivial baseline, strongest available):** log OR 0. No early read-out of any trial exists, so carry-forward is
  not available for the headline test; rule 2b's strongest trivial baseline is M0. (M1 is itself a trivial baseline.)
- **M3 anchored (secondary only):** same call with anchor=True; whatever anchor level the shipped code finds is recorded
  in the freeze (cenerimod has a same-molecule SRI-4 anchor; the others are expected to fall back to pooled).
- **C0 carry-forward slots (reserved, descriptive):** ianalumab SIRIUS-SLE 1/2 register SRI-4 at wk36 and wk60. If both
  visits are reported, C0 = wk36 log OR carried to wk60 is scored against M1 and M0 at wk60. Rule 8: the early-read-out
  rule is settled and re-opens only with >= 5 frozen trials reporting both visits; 2 trials cannot re-open it.

## 6. Metric and scheme
Prospective (forecast frozen before any result exists). **Primary metric: IV-weighted MSE of the observed log OR (drug vs
placebo, 0.5 continuity correction, se as `lumen_forecast.load`) at the registered primary endpoint and visit**, over every
active dose arm vs its own control (SELECT-SLE: each study's arms vs that study's placebo), weights 1/se_obs^2.
The endpoint must be the registered one (SRI-4 or BICLA) at the registered visit; a trial reporting only another visit or
endpoint is unscorable (listed, counts against n).

## 7. Minimum n
>= 5 scored trials. Below -> UNDERPOWERED (interim report, no decision).

## 8-9. Decision rules (evaluated at the FIRST W1 run where >= 5 trials are scored = the decision look)
Let A(.) = primary IV-MSE.
1. UNDERPOWERED: < 5 scored trials.
2. **FAIL (rule 2b):** A(M1) >= 1.10 x A(M0).
3. **PASS:** A(M1) <= 0.90 x A(M0).
4. INCONCLUSIVE: otherwise.
Later looks (more trials) update n and the label; a decision flip after the decision look requires >= 10 scored trials
and the same rule. No other stopping.

## 10. Secondaries (labelled; none changes the decision)
- S1 equal trial weights (mean of per-trial IV-MSE).
- S2 one molecule per same-molecule pair (cenerimod, ianalumab): keep the trial whose results post first.
- S3 M3 anchored vs M1 vs M0.
- S4 SRI-4-only and BICLA-only subsets, with n.
- **S5 95% interval coverage of M1 (named in the route):** share of scored arms whose observed log OR lies inside M1's 95%
  interval; < 0.80 at the decision look -> flagged "pooled intervals too narrow" in the preprint (descriptive; no
  CLI change from this test alone).
- S6 drop press-release-sourced arms.
- S7 C0 slots (section 5).

## 11. Labels (rule 12)
Supported only at >= 10 scored trials (n stated); otherwise Emerging at most. "Verifiably pre-registered" applies if every
outcome source's fetch_utc in the manifest is later than this PREREG commit and the extraction commit precedes scoring.

## 12. Output
`prospective/w1-sle/score_w1sle.py` (written at the first scoring run, reading the frozen JSON unchanged) ->
`w1sle_result.json` via `scripts/provenance/stamp.py` (prereg path `experiments/lumen-bench/prospective/w1-sle/PREREG-W1-SLE.md`).

Research and simulation only. Not clinical or drug-development guidance.

## Input md5 (frozen at PREREG)
- 2c11b9a5aba5c047c0817140c7d5cbe6  experiments/lumen-bench/lumen_forecast.py (last changed a24cea8)
- b37dd7cd72c9e912ba079088532c3765  experiments/lumen-bench/sle/sri4_table.csv
- 99eb279e207302973fef6899f0e86722  experiments/lumen-bench/sle/bicla_table.csv
- 55385007c38b29ee11b52f57f698f332  experiments/lumen-bench/prospective/w1-sle/candidates_w1sle.csv
- d7b40ceed35e84074657c4f310d92901  experiments/lumen-bench/prospective/w1-sle/CATALOGUE-W1SLE.md
- 168e99dc3de80a884ebf2327c52f4cb5  experiments/lumen-bench/prospective/w1-sle/rescreen_run2.csv
- 4165bab8f72dee66fd8929e3edff191d  scripts/provenance/fetch.sh
- ca14be70b213b98bbc2b830e2b092329  scripts/provenance/stamp.py
