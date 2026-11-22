# Round 99 result (Q14/Q29) — run against PREREG.md (commit 0e99c70)

Positive control (QC, not a test): IFN6 separates SLE from healthy, AUC 0.93 (GSE65391), 0.91 (GSE72326), 0.80 (GSE99967).

| Held-out AUC for active nephritis | T1 serology | T0 erSLEDAI | BASE | FULL (+IFN6) | IFN6 alone | FULL-BASE [95% CI] |
|---|---|---|---|---|---|---|
| LOCO train A -> test B (n 60, 23 LN) | 0.40 | 0.76 | 0.66 | 0.64 | 0.53 | -0.02 [-0.08, +0.04] |
| LOCO train B -> test A (n 158, 76 LN) | 0.48 | 0.62 | 0.54 | 0.55 | 0.37 | +0.01 [-0.02, +0.04] |
| Within-A 5x20 CV | 0.53 | 0.61 | 0.63 | 0.64 | 0.62 | +0.007 |
| Within-B 5x20 CV | 0.60 | 0.73 | 0.72 | 0.70 | 0.37 | -0.021 |

- H1 FAIL: IFN6 adds nothing beyond extra-renal SLEDAI + serology (both LOCO deltas < 0.05, CIs span 0).
- H2 FAIL: IFN6 alone 0.53 and 0.37 in the LOCO test cohorts (rule: >= 0.60 and > serology).
- H3 FAIL, and inverted: GSE99967 biopsy-proven active LN vs active SLE without LN, IFN6 AUC 0.37 [0.21, 0.53] — LN patients
  have LOWER blood IFN than active non-renal SLE.
- Direction differs by cohort: pediatric A raw AUC 0.63 (IFN higher with active renal items, rho 0.44 with erSLEDAI, i.e. it
  tracks global activity); adult B 0.53; adult biopsy-proven C 0.37.
- Sensitivity (continuous anti-dsDNA in B): BASE 0.72, FULL 0.70 — same answer.
- Caveat: transported serology coefficients do worse than chance (0.40, 0.48): the binary low-complement/anti-dsDNA codings
  are not equivalent across cohorts; erSLEDAI is the only baseline that transports.
