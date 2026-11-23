# RESULT P6 / Q56 / Q59 — round 118 (2026-11-23, scheduled run)

PREREG: `PREREG-P6.md` (commit 435dcef, pushed before any titre value was opened). Data record: commit 5fe2cfd.

## Decision this can change
Add a slow CD19-negative plasma-cell reservoir state to LUMEN (if two-pool beats one-pool on held-out curves) or not.

## Verdict: UNDERPOWERED -> decision B (state NOT added)
Pass rule 3 requires >= 8 CD19/CD20-only (C1) curves with baseline + >= 3 numeric post-treatment titres, last >= 90 d.
**Qualifying C1 curves found: 0.** T1/T2/T3 were not run; no model was fitted.

What was searched (open sources only):
| Source | What is numeric | Why it fails inclusion |
|---|---|---|
| Mackensen 2022 Nat Med (CD19 CAR-T, n 5), source data | anti-dsDNA baseline + month 3 per patient | 1 post-treatment point |
| Relma-cel phase 1, EClinicalMedicine 2025, PMC12083988 + supplements (n 8) | baseline anti-dsDNA only | time course in figures only |
| Erlangen long-term follow-up, Z Rheumatol 2025, PMC12488764 (n 11) | seroconversion yes/no | no titres |
| BCMA-CD19 cCAR extended follow-up, J Hematol Oncol 2026, PMC12857108 (C2) | SLEDAI | titres in figures only |
| Feng 2025 Nat Med (CD19+BCMA co-infusion, n 15, C2), supplementary tables | BCR, AEs, proteomics | no titre table |
| CT.gov results (lupus_with_results) | anti-dsDNA by visit exists only for non-depleting drugs (belimumab withdrawal, anti-IL-21, KPG-818); NOBILITY single visit | wrong class / 1 point |
| PMC open-access table scan, 558 articles (rituximab/obinutuzumab, CAR-T, daratumumab, bispecific + lupus + anti-dsDNA) | 9 numeric anti-dsDNA rows | none is a C1 time series (scan_rows.txt; scan_cols.txt empty) |

Per-patient titre time courses after CD19/CD20 depletion exist publicly only as figure panels. The PREREG allowed
digitised values only as a sensitivity set, and the n >= 8 gate is on the numeric set, so digitising could not change
the decision; it was not done.

## Descriptive, not pre-registered as a test (Emerging, n 5 patients, one source)
Mackensen month-3 / baseline anti-dsDNA ratios: 0.013, 0.028, 0.083, 0.66 (baseline 14.3, near normal: uninformative),
0.114; median 0.083. Pure IgG wash-out with production switched off on day 0 (half-life ~21 d) gives ~0.05 at day 90.
The residual at month 3 is an UPPER bound on the CD19-unreachable share (any fast-pool remainder adds to it), so in
4 informative patients that share is <= 1-11% of baseline anti-dsDNA production. This fits round 47's corrected
reading (low autoreactive long-lived-plasma-cell share) and says the round-47 sweep above g_l ~0.15 is not needed.

## Why the test is structurally weak even with more data (identifiability note)
In responders anti-dsDNA falls below the assay limit by month 2-3. With the LOD floor applied (PREREG rule), one-pool
and two-pool predictions coincide whenever the plateau (1-f) x baseline is below the LOD — true for every informative
patient above. Responder anti-dsDNA curves therefore cannot identify the reservoir. It is identifiable only from
(a) non-responders / relapsers after CD19-only depletion, or (b) paired protective (vaccine) IgG, which is a
long-lived-plasma-cell product by construction. Q59's relapse question: the Erlangen series has 1 non-seroconverter
of 11 (patient 9, SLEDAI 22 at last visit) — n 1, descriptive only.

## Decision taken
- LUMEN keeps ONE antibody pool. Slow reservoir state stays **Hypothetical**; preprint says it is not identifiable from
  open titre curves.
- Record the bound: if a future LUMEN decision needs the CD19-unreachable autoantibody share, use <= 0.11 (upper bound,
  n 4 informative patients), not round 47's 0.6.
- Re-open only if an open numeric source appears with >= 8 CD19/CD20-only patients x >= 3 post-treatment titres,
  ideally with vaccine IgG or relapse labels; then run PREREG-P6 unchanged (W2 weekly scan checks this).

## Tests counted
0 pre-registered tests run (underpowered by rule). 1 descriptive computation (ratios). 0 models fitted.
