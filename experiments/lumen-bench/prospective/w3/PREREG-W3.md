# PREREG-W3 (round 169, 2026-11-28 UTC): standing rolling prospective freeze of every new SLE/LN RCT

Written from `docs/PREREG-TEMPLATE.md` after the run-1 catalogue (`CATALOGUE-W3.md`, commit e8205db, pushed 13:27 UTC).
Rule 12b: this run pushes the PREREG and ENDS. The first-batch freeze (`freeze_w3.py`, committed with a `--selftest`
before its first real execution) runs in a later scheduled run >= 1 h after this push. No outcome of any candidate
exists or was read (20 of the 22 registry-passers are not yet enrolled/completed; see the catalogue).

## 1. Item / round / date
Queue W3 (from POW1, round 167); round 169; 2026-11-28.

## 2. Decision this can change (rule 7)
Whether any rule-8 settlement (CAL1 FAIL, O1 INCONCLUSIVE, O2 FAIL, O3 INCONCLUSIVE, the early-read-out rule) can ever
reopen, and on what pooled prospective n. POW1 (round 167, `pow1_result.json` md5 below) showed W1's 21 trials cannot
reopen O1/O2/O3 and detect only gross (>= 0.30) calibration gaps.
- Catalogue gate met (17 eligible >= 10) -> this protocol is adopted; settlements stay closed until the pooled
  prospective set P reaches the n in section 8, and LUMEN prints the re-open n and the projected earliest year.
- If later W3 batches find accrual < 5 eligible trials/year for 2 consecutive years -> LUMEN prints the settled lines as
  'closed for the foreseeable future'.

## 3. What has already been seen (disclosure)
- Registry metadata of the 22 registry-passers (titles, arms, primary measures, dates, enrolment). No outcome.
- Everything in the historical tables (LN 20 trials `resp_table.csv`, SLE 40 SRI-4 / 26 BICLA) and all CAL1/O1/O2/O3/
  REC1/POW1 results. None of the 17 is in any table. The analyst model's training cutoff is June 2026; no candidate had
  reached registry primary completion by then.

## 4. Trials
**Pooled prospective set P = W1-LN (13, frozen 2026-11-21) + W1-SLE (8, frozen 2026-11-26/27) + every W3 batch.**
W1/W1-SLE/REC1 freezes are not touched; W3 adds the fields in section 5 only for W3 trials. For the 21 earlier trials,
the O1/O2/O3/REC1 records already frozen are the ones scored (fields not frozen for them are 'not available', not
back-filled after 2026-11-28; exception: assurance, section 5 item 8, may be added for the 21 in the first-batch freeze
because none has reported, per R-20261128-6).

**First batch (frozen inclusion list, 17 trials; `manual_screen.csv` 'eligible'):**
| NCT | molecule (registry) | endpoint | registered primary visit | table |
|---|---|---|---|---|
| NCT06791772 | CM313 | SRI-4 | wk24 | sri4 |
| NCT06859931 | TQB3702 | SRI-4 | wk24 | sri4 |
| NCT07447986 | SG301 | SRI-4 | wk25 | sri4 |
| NCT07509892 | TQH3906 | SRI-4 | wk32 | sri4 |
| NCT07299422 | SHR-2173 | SRI-4 | wk24 | sri4 |
| NCT07760155 | ESG206 | SRI-4 | wk25 | sri4 |
| NCT07311200 | D-2570 | SRI-4 | wk32 | sri4 |
| NCT07405970 | MIL62 | SRI-4 | wk52 | sri4 |
| NCT07438496 | nipocalimab | SRI-4 | wk52 | sri4 |
| NCT06535412 | IMC-002; IMM0306S (each arm vs placebo) | SRI-4 | wk24 | sri4 |
| NCT07409181 | BI 3000202 | SRI-4 | wk32 | sri4 |
| NCT06617325 | dapirolizumab pegol | BICLA | wk48 | bicla |
| NCT07515014 | E6742 | BICLA | wk24 | bicla |
| NCT05835310 | anifrolumab (paediatric, Part B) | BICLA | wk52 | bicla |
| NCT07937501 | IPG11406 | CRR | wk24 | ln |
| NCT03673748 | allogeneic MSC | CRR | wk24 | ln |
| NCT07678203 | dapagliflozin (add-on vs SoC monitoring) | KDIGO 2024 renal remission, scored as CRR | registry time frame | ln |

Targets are read from registry intervention text at freeze; unknown -> `unknown` (used only by anchored M3, secondary).
Exclusion at freeze: any trial whose CT.gov record shows hasResults=true, or for which a topline is found in the
freeze-run re-screen (web/PubMed titles only, yes/no recorded), is excluded and listed.

**Later batches (standing).** Every W1 weekly run re-pulls the two CT.gov queries via `fetch.sh`, re-runs
`catalogue_w3.py` unchanged, applies the same manual-screen reasons, and freezes newly eligible trials in a commit after
the catalogue commit, with `freeze_w3.py` unchanged. A trial is eligible only if no result can exist (not yet at
registry primary completion) or a topline screen finds none. Catalogue rule changes need an amendment (restarts the
1-h clock, rule 12b).

## 5. Frozen fields per W3 trial (write-once via `lumen.provenance.write_frozen`; md5 appended to FROZEN-MD5.tsv)
1. **M1** shipped drug_effect pooled prior: `lumen_forecast.forecast(molecule, target, R=load(table), anchor=False)`,
   `lumen_forecast.py` md5 below; mean m, sd, 95% interval.
2. **M0** null (log OR 0) - implicit.
3. **P1** printed-style P_benefit (as CAL1) and **B** printed base rate (as the shipped CLI for that table).
4. **REC1 competitors** a/b/c with REC1's frozen w = 0.26 and van Zwet constants (`freeze_rec1.py` functions, unchanged).
5. **O1a/O1b** candidate predictives (`lumen.candidates`, as `frozen_O1_2026-11-27.json`).
6. **O3 P1** phase-covariate predictive (as `frozen_O3_2026-11-27.json`, registry phase).
7. **O2 R1, B1, B2** control-rate predictives (as `run_o2.py`/`frozen_O2_2026-11-27.json`) for each table.
8. **Assurance** A = P(observed log OR / se_ant > 1.96) = Phi((m - 1.96 se_ant) / sqrt(sd^2 + se_ant^2)), se_ant from
   registry enrolment split equally across arms and the table's median control rate (`freeze_rec1.se_ant2`), (R-20261128-6).
9. **C0 slots:** registered earlier responder visits (registry secondary outcomes), for the early-read-out rule.

## 6. Scoring (later W1 runs, when a trial reports)
Same as PREREG-W1-SLE s4/s6: registered endpoint at registered visit, every active arm vs its own control, log OR with
0.5 continuity correction; sources via `fetch.sh` into `w3/sources/` + MANIFEST, extraction CSV committed before the
scoring commit; result JSONs via `lumen.provenance` write-once + `stamp.py`. A trial reporting another endpoint/visit is
unscorable (listed, counts against nothing). Press-release counts flagged, dropped in a sensitivity.

## 7. Metrics per settlement (each the ORIGINAL PREREG's primary metric and rule, unchanged)
CAL1: PREREG-CAL1 s8-s9 (gap g, Brier vs H and 0.5). O1: PREREG-O1 s8-s9 (log score, 0.05 nats). O2: PREREG-O2 s9-s10
(NLPD vs better of B1/B2, LN guard). O3: PREREG-O3 s9-s10 (IV-MSE, gate margin 0.10). Early-read-out rule: RT-1 rules on
trials reporting both visits. REC1: PREREG-REC1 unchanged (its own trials only). M1 vs M0 (pooled headline): PREREG-W1-SLE
s8-9 rule applied to all of P (secondary; W1-SLE's own 8-trial look keeps its own rule).
**Added for all margin rules (POW1):** a re-open (verdict different from the settled one) also needs the paired
trial-bootstrap 95% CI (2000 resamples, seed 20261128) of the margin statistic to exclude 0.

## 8. Re-open n (pre-stated from POW1; decision looks only at these n, scored trials in P that carry the needed fields)
| settlement | re-open look at n scored trials | 2nd look |
|---|---|---|
| O2 control-rate prior (FAIL) | 20 | 50 |
| CAL1 P_benefit calibration (FAIL) | 50 | 100 |
| O1 interval model (INCONCLUSIVE) | 50 | 100 |
| O3 phase covariate (INCONCLUSIVE) | 50 | 100 |
| early-read-out rule | 5 trials reporting both registered visits (rule 8) | 10 |
Below the look n: interim numbers only, labelled UNDERPOWERED, no decision. No other looks, no stopping.

## 9. Decision rules at a look
The original PREREG verdict is recomputed on P at the look. Settled verdict reproduced or INCONCLUSIVE -> settlement
stands (label updated with prospective n). Different verdict AND section-7 CI excludes 0 -> settlement reopens: the
original PREREG's PASS branch / FAIL branch action applies, through `lumen.gate`, at most Supported with n stated.

## 10. Labels
Supported only at a look with n >= the section-8 n (stated); "verifiably pre-registered" for W3 trials if every outcome
fetch_utc > the freeze commit time and > this PREREG's commit time.

## 11. Output
`w3/freeze_w3.py` (selftest-first) -> `w3/frozen_W3_batch1_<date>.json`; later batches `frozen_W3_batchK_<date>.json`;
scoring `w3/score_w3.py` (written at the first look) -> `w3_result.json`.

Research and simulation only. Not clinical or drug-development guidance.

## Input md5 (at PREREG)
- 2c11b9a5aba5c047c0817140c7d5cbe6  experiments/lumen-bench/lumen_forecast.py
- b37dd7cd72c9e912ba079088532c3765  experiments/lumen-bench/sle/sri4_table.csv
- 99eb279e207302973fef6899f0e86722  experiments/lumen-bench/sle/bicla_table.csv
- 3ce67bdb51c6597bce6666d914618778  experiments/lumen-bench/resp_table.csv
- 99f696da016fb51762cb668402e005b4  lumen/candidates.py
- 3d034506004f13456a7ed32d131f6436  lumen/provenance.py
- 8417a36106408161135506a53ef2fe63  experiments/lumen-v3/run_o1.py
- 831af043b12c559a30bf4f8c7dbc7314  experiments/lumen-v3/run_o2.py
- e011d90d0545bf8897e7d80f8701b404  experiments/lumen-v3/run_o3.py
- 2111dbb549398369e7b18d9654fb3b3f  experiments/lumen-v3/run_cal1.py
- 471998ce8246821e0377091da65a88ff  experiments/lumen-v3/freeze_rec1.py
- 56d7ee671bc78c53d9f9f34d01e9f91c  experiments/lumen-v3/pow1_result.json
- 2ad1a322cb8f5b46aa44681300185a99  experiments/lumen-bench/prospective/w3/catalogue_w3.py
- 4c7e3b10ec0061ff8f06636995d4d810  experiments/lumen-bench/prospective/w3/w3_candidates.csv
- e255fde51282fb90467c0c315975df69  experiments/lumen-bench/prospective/w3/manual_screen.csv
- 1482093b1e63c81cb59a18210ac9aa99  experiments/lumen-bench/prospective/w3/CATALOGUE-W3.md
- 4165bab8f72dee66fd8929e3edff191d  scripts/provenance/fetch.sh
- ca14be70b213b98bbc2b830e2b092329  scripts/provenance/stamp.py
