# RT-20261124-1, run 1: frozen catalogue and extraction rules (no outcomes extracted)

Written 2026-11-24 ~13:50 UTC (round 133, scheduled run). Reviewer route RT-20261124-1 (claude/REVIEW-2026-11-24.md).
This file and `candidates_frozen.csv` are frozen at this commit. Run 2 (>= 1 h later, rule 12b) writes the PREREG on
top of them; run 3 fetches the full texts and results through `scripts/provenance/fetch.sh`, extracts, then scores.

**Decision this can change** (from the route): the CLI headline when `--early` is given. If carry-forward beats pooled by
>= 10% and is not >= 10% worse than null in LN -> headline becomes carry-forward (Emerging/Supported at that n). If pooled
wins -> keep pooled and drop the early read-out from the CLI output. If fewer than 10 new trials qualify -> UNDERPOWERED:
park Q61-type work; early read-outs stay unvalidated.

## 1. What was looked at in this run (and what was not)
- ClinicalTrials.gov API v2, randomised interventional LN and SLE studies, completed / terminated / active-not-recruiting,
  **without posted results** (`raw/ctgov_noresults_*.json`; fields: IDs, titles, phase, enrolment, registered outcome
  measures and time frames, sponsor, countries, references). Script `catalogue_ctgov.py`. Registry entries without
  results contain no outcome numbers.
- PubMed E-utilities: query in `raw/MANIFEST.tsv` (lupus in title AND randomized controlled trial[pt] AND placebo /
  standard of care / standard therapy in title-abstract, 2000-2026; 296 records). Script `screen_pubmed.py` parses
  **titles, years and registry accession numbers only**; abstracts were downloaded inside the efetch XML but were not
  parsed, printed or read. Screen: 70 records map to a tabulated NCT, 96 titles are clearly secondary, 130 screened by
  title by hand -> `candidates_frozen.csv`.
- Not done: no full text, supplement, abstract or CT.gov results section of any candidate was opened. No outcome
  number was extracted.
- **Analyst knowledge (disclosure).** I know the headline primary results of several candidates from training-era
  literature (e.g. EMBODY 1/2, EXPLORER, ROSE, BLISS-BELIEVE missed their primary endpoints; telitacicept phase 2b and
  phase 3 met theirs). I do not know their early-visit by-arm counts. The time-forward early-vs-late comparison on these
  trials has not been computed by this project.

## 2. Candidate list
`candidates_frozen.csv`: 36 candidates (29 SLE, 7 LN): 28 from PubMed titles, 8 CT.gov studies with no posted results.
Excluded before listing (title or registry level): anything already in `resp_table.csv`, `sle/sri4_table.csv`,
`sle/bicla_table.csv`, the NCT maps and `EXCLUDED` set of `sle/build_sle_table.py` (this keeps MUSE, EMBRACE,
daxdilimab, LOTUS and brentuximab out: they were already seen and excluded by pre-registered definition rules), the 13
W1 frozen prospective trials (`prospective/w1/w1_frozen_2026-11-21.json`, e.g. IRIS NCT05138133, iptacopan LN
NCT05268289: they stay prospective), head-to-head induction/maintenance comparisons (no placebo/background-only arm),
flare-prevention or withdrawal designs, phase 1 studies, supplements/diet/hormones/statins/devices, and secondary or
pooled analyses.

## 3. Eligibility rules applied in run 3 (mechanical, before any score is computed)
A candidate is INCLUDED only if all hold; each exclusion is logged with its reason.
1. **Not already tabulated.** The registry ID given in the paper or CT.gov record is not in the tables, the
   `build_sle_table.py` maps / `EXCLUDED` set, `build_resp_table.py`, or the W1 frozen set (e.g. N07 SLEek and N19
   Lupuzor are dropped here if their NCT is NCT03978520 / NCT01135459).
2. **Add-on randomised design.** Investigational drug (+ background therapy) vs placebo or background therapy alone;
   for multi-period designs only the randomised placebo-controlled period counts.
3. **Open source.** CT.gov posted results, or a full text / supplement downloadable without an account (PMC, publisher
   open access). Paywalled with no open supplement -> excluded "not open". Unpublished -> excluded "no results".
4. **Index.** SLE: SRI-4 responders (primary analysis). If SRI-4 is not reported but BICLA is, the trial enters the
   BICLA secondary only. LN: complete renal response as the trial defines it; trials reporting only continuous
   proteinuria -> excluded.
5. **Visits** (identical to PREREG-RT1 section 4). SLE: early = latest reported visit with 12 <= week <= 26; target =
   latest with 48 <= week <= 52. LN: early = latest with week <= 26 (and >= 12); target = latest with 48 <= week <= 104.
   Months convert as week = month x 52/12, rounded to the nearest reported visit.
6. **Numbers.** By-arm responder counts and arm N for the randomised (ITT / full analysis set) population, from tables
   or text. Percentages with arm N are converted to counts (rounded). Values read off figures are NOT used unless the
   figure prints the number. Non-responder imputation as reported; no re-imputation.
7. **Arms.** Each active dose arm vs the shared control is a separate arm (as RT-1).
8. **Date for the time-forward order.** CT.gov `resultsFirstPostDate` if results are posted; otherwise the earliest
   PubMed publication date (electronic). Trials are inserted among the 60 tabulated trials by this date; only new
   trials are scored, each forecast from all trials (old and new) of its stratum dated strictly earlier.
9. **Minimum n:** >= 10 new trials with both visits in the primary (SRI-4 / CRR) set, else UNDERPOWERED.

## 4. Provenance plan
Run 3 downloads every source through `scripts/provenance/fetch.sh` into `newtrials/sources/` with a committed
manifest, so fetch_utc > PREREG commit time can be checked; the extraction table is committed before the scoring
commit. If this holds, the test is the first in the project that can be called verifiably pre-registered (rule 12b,
round 133), with the analyst-knowledge disclosure above.

## 5. Expected yield (stated before extraction, from registry designs and titles only)
Several candidates are short (12-24 wk: N09, N10, N11, N18, N19 likely), several CT.gov studies have no publication
(N29-N36), and some SLE trials report BICLA or BILAG rather than SRI-4. A realistic range is about 5-14 qualifying trials,
so **UNDERPOWERED is a likely outcome** and would be reported as such.
