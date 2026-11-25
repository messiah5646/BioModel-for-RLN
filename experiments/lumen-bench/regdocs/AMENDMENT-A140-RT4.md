# AMENDMENT to PREREG-RT4 from audit round 140 (2026-11-25), before any RT-4 run-3 fetch

Written and pushed BEFORE run 3: `experiments/lumen-bench/regdocs/sources/` does not exist at this commit.
PREREG-RT4.md (67f2d26) and AMENDMENT-R20261125-RT4.md (414a21e) are left unedited.

## Finding (audit r140)
1. **LOTUS (NCT03517722) early AND late by-arm SRI-4 values are already in a committed, frozen PREREG input.**
   `data/real/ctgov/lupus_with_results.json.gz` (md5 1b31cd52..., listed in PREREG-RT4 "Input md5") holds CT.gov
   results for "SRI-4 Composite Response at Week 24" and "... at Week 52" by arm. They were never tabulated
   (`build_sle_table.py` EXC regex drops 'composite response'), and round 138 read titles only; the audit read titles
   only. Still, PREREG-RT4 section 2 ("early values new") cannot be claimed for LOTUS.
   MUSE (NCT01438489) CT.gov results in the same file are OCS-taper-conditioned (Day 169 / Day 365) and fail the
   PREREG section 4 plain-SRI-4 rule; the same caveat applies to MUSE if run 3 includes it.
2. **Telitacicept ph2b source PMC10958275 full text was committed in round 136** (newtrials/sources, fetched
   2026-11-24 16:24Z, before PREREG-RT4). Masked-digit scan by the audit: its text holds no by-arm early SRI-4 value
   (early numbers are in Fig 2B / an unfetched supplement, as RT-3 exclusions.csv says). Disclosure only.
3. No other Table A trial has early-visit responder outcomes in the committed CT.gov file (titles checked: BLISS-LN
   wk104 only; BLISS-SC/BLISS-76/PLUTO/BLISS-NEA wk52 (+wk76, OLE); TULIP-2 BICLA wk52; ACCESS and telitacicept absent).

## Added (descriptive; models, metric, minimum n, decision rules and the shipped branch are UNCHANGED)
- **S9 clean-early subset.** Primary and S1 re-scored without LOTUS and MUSE (trials whose early value sat in a
  committed pre-PREREG file). Report n and which section-9 rule S9 would select; if it differs from the primary, the
  round doc says so in its first lines.
- **Label text** for rt4_result.json and the round doc: "12b-timed (late outcomes already seen; early values new for
  8 of 10 Table A trials; LOTUS/MUSE early values were in committed CT.gov data, unread)".
- **Population check for LOTUS** (PREREG section 4 rule): if run 3 extracts LOTUS from PMC9606504, the CT.gov counts may
  be used as the population cross-check (> 2 responders difference in any arm -> excluded from primary, kept in S3).
- `rt4_result.json` carries key `S9_clean_early`.

Research and simulation only. Not clinical or drug-development guidance.
