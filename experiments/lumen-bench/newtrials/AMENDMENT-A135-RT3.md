# Amendment A135 to PREREG-RT3 (audit round 135) — labels and one reported sensitivity; models, metric, pass rules unchanged

Committed 2026-11-24 ~15:35 UTC by the round-135 audit, BEFORE any run-3 fetch, extraction or scoring. PREREG-RT3.md
(68900fb) is not edited. This file changes no model, metric, subset, minimum n or decision rule (sections 5-9 stand), so
it cannot move the decision; it fixes what the result may be CALLED and adds one reported sensitivity.

## Findings that motivate it (from committed files only; no new outcome data opened)
1. **N30 JASMINE (nipocalimab, NCT04882878) early AND late by-arm numbers are already in the repo.** Round 131
   (`prospective/q61/sweep_2026-11-23.md`, commit b748458, 2026-11-23 22:27 UTC) records SRI-4 wk 24 and wk 52 percentages
   by arm. PREREG-RT3 section 2 ("no early-visit by-arm counts") is therefore wrong for N30.
2. **Abstracts of all 28 PubMed-sourced candidates (N01-N28) were committed before the PREREG** inside
   `raw/pubmed_efetch.xml` (836f753, 13:38 UTC; disclosed in CATALOGUE section 1 as "downloaded but not parsed").
   A blind flag scan (flags only, no text printed; audit-r135/abstract_flags.tsv) finds numeric content in 27/28 abstracts
   and an early-visit (12-26 wk / 3-6 mo) sentence carrying a number in 5 (N01, N11, N12, N19, N24).
   Rule 12b reserves "verifiably pre-registered" for outcome data that ARRIVE after the PREREG; re-fetching these sources in
   run 3 would give later fetch_utc values (PREREG section 11) while the outcome text was already committed earlier.

## Amended label rule (replaces PREREG-RT3 section 11, second sentence)
- Trials N01-N28 and N30 are **12b-timed (outcome text committed before PREREG; not parsed; earlier looks not excludable)**.
- Only candidates whose outcome data first enter the repo via run-3 `fetch.sh` (N29, N31-N36 if included) can be called
  verifiably pre-registered. The run-3 result is called "verifiably pre-registered" only if the decision is reached on that
  clean subset alone with >= 10 trials; otherwise the whole result is "12b-timed".
- Supported/Emerging thresholds (>= 10 new trials) are unchanged.

## Added reported sensitivity (descriptive; does not change the decision)
- S7: primary subset without N30, and primary subset restricted to the clean candidates (N29, N31-N36), each with n
  trials/arms and all four model IV-MSEs. If the S7-without-N30 decision category differs from the primary one, the
  round doc says so in its first line.

Research and simulation only. Not clinical advice.
