# Q16 / round 106 result: damage floor vs slower rate

PREREG faef233, amendment 1 3327efa (both committed before any outcome extraction).

## Primary (scored) test: BLOCKED, not failed
No open source met the eligibility rules: no open paper, supplement or repository gives proteinuria at >=3 visits,
or time-to-response curves, or patient-level longitudinal data, split by baseline NIH chronicity index (CI).
Search: PubMed, Europe PMC (~15 field-restricted queries incl. figure captions/tables; ~900 PMC full texts scanned
via E-utilities for CI-group captions/columns), Dryad, Zenodo, figshare, GEO sample index (62,446 samples:
CI appears only in GSE99967 (blood, cross-sectional) and GSE42648 (high/low only, no proteinuria)).
Near misses: AURORA repeat-biopsy sub-study (PMC12479183; patient-level CI but a single month-18 response);
IFTA paper (PMC13130736; final response only); class IV+-V paper (PMC12119475; KM by CI>4 but outcome is kidney
survival); several HR/OR-only papers (PMC12914276, PMC9457419, PMC11429000). Paywalled candidates (Hill 2001 KI,
Rijnink 2017 CJASN) are outside the open-data policy.

## Exploratory (NOT scored): composite-risk KM, Saitama triple-therapy cohort
Arthritis Res Ther 2026, doi 10.1186/s13075-026-03876-w, PMC13523271, Fig 2: cumulative CR by composite risk group
(low n=67 / moderate n=32 / high n=71; CI>3 is one of the risk factors: 0% / 12.5% / 38% of each group).
Digitised at 3, 6, 12, 24, 36 months; PREREG amendment-1 models applied unchanged (exploratory_saitama_km.py).

| groups | pooled AICc | floor AICc | rate AICc | floor params | held-out MAE at 36 mo (fit on 3,6 mo): pooled / floor / rate / carry-forward |
|---|---|---|---|---|---|
| low vs high | -21.4 | **-59.2** | -38.1 | k 0.37/mo, c 0.93 vs 0.38 | 0.26 / **0.053** / 0.225 / 0.095 |
| all three | -41.7 | **-95.2** | -56.9 | k 0.35/mo, c 0.94 / 0.67 / 0.38 | 0.19 / **0.057** / 0.21 / 0.11 |

The floor (different eventual fraction, shared rate) beats the rate model by >20 AICc and beats both trivial
baselines on the held-out 36-month point. The rate model needs an implausible k of 0.03/mo for the high group.

## What it means (and doesn't)
- Shape supports a CEILING on response, reached at the same speed, not slower response, in high-risk patients.
- It does NOT isolate chronicity: risk groups also differ in eGFR, proteinuria and other factors; CI>3 is in only
  38% of the high-risk group.
- An aggregate "floor" is the same signature as a MIXTURE (a fraction who never respond), so this cannot tell a
  per-patient proteinuria floor from a non-responder subgroup. That is the identifiability limit stated in the PREREG.
- Evidence label: Emerging (one cohort, aggregate, composite grouping, exploratory).
