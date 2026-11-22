# RESULT round 97 (Q11): blood IFN / plasmablast scores add nothing to UPCR + SLEDAI for MMF response

PREREG frozen at 181690a. Code: run.py. Output: result.json, patients.csv (one baseline row per patient).
Cohort: GSE224705 = treated LN patients from the Hopkins Lupus Cohort (López-Domínguez 2024, doi 10.21203/rs.3.rs-3783877/v1).
Label = per-patient responder flag (constant within every patient: 44/44 MMF, 20/20 AZA).

| Model (leave-one-patient-out AUC) | MMF (n 43; 33 R / 10 NR) | AZA (n 18; 11 R / 7 NR) |
|---|---|---|
| Chance | 0.50 | 0.50 |
| SLEDAI alone | 0.65 | 0.82 |
| UPCR alone | **0.84** | 0.00 (LOO artefact; raw 0.70 inverted) |
| IFN6 + PB4 alone | 0.00 (LOO artefact; raw 0.50) | 0.61 |
| Base: UPCR + SLEDAI | **0.87** | **0.83** |
| Full: + IFN6 + PB4 | 0.84 | 0.79 |
| Delta full - base (95% bootstrap) | -0.025 (-0.08 to +0.02) | -0.04 (-0.20 to +0.13) |

- **H1 (MMF): FAIL.** Delta is negative; the CI excludes the +0.05 bar (upper bound +0.02), so this is "no added value", not just low power.
- **H2 (AZA): FAIL.** Negative delta; CI too wide to say much (n 18).
- **H3 trivial baselines:** UPCR alone carries almost all the MMF signal. Non-responders started with median UPCR 0.78 vs 0.11 g/g; raw IFN6 AUC 0.497.
- Sensitivity (UPCR >= 0.5 g/g, n 13, 4 R / 9 NR): too small; delta -0.31.
- Exploratory 5-fold x 50 (not pre-registered; LOO-AUC is biased to ~0 for uninformative features): same ranking.
  MMF base 0.86, full 0.84, IFN6 0.30. AZA base 0.83, full 0.81, IFN6 alone 0.81.
- Exploratory signal: in AZA, a higher baseline IFN6 went with NON-response (raw AUC 0.13, i.e. 0.87 inverted, n 18). But IFN6 tracks SLEDAI there (rho 0.56), so it adds nothing once SLEDAI is in the model. Emerging at most.

Deviations: the UPCR field is in g/g despite its "mg/g" name, so the sensitivity threshold was applied as 0.5 g/g. The TXNDC5 probe is annotated TXNDC5///BLOC1S5.
Limits: baseline = first sample in the cohort, which is often already on the drug (not pre-treatment). The response definition is the authors' own (field name sri4_response). Single centre.
