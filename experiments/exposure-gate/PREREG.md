# PREREG Q15 (round 105): does early drug exposure / B-cell engagement add >= 0.05 AUC over month-3 UPCR for 12-month complete renal response?

Frozen 2026-11-22, before extracting any effect size beyond what round 93 already quoted
(round 93 quoted: MPA AUC 42.6 vs 33.1 in CR vs NR, n=100; LUNAR CR 47% with complete depletion vs 13% without, n=68;
month-3 urine CD206/CD163 0.89-0.91 vs month-3 proteinuria 0.80, AMP). No patient-level data has been opened.

## Part A — the real test (patient-level), as written in the queue
Data required (open only: GEO, ClinicalTrials.gov, Zenodo/figshare, paper supplements/tables): per-patient
(i) month-1 (<= 6 weeks) exposure (MPA AUC0-12, trough, or peripheral CD19 count after anti-CD20),
(ii) month-3 UPCR, (iii) 12-month complete renal response (any standard definition).
Search: PubMed/PMC full text of MPA-exposure and B-cell-depletion LN cohorts (open-access first), GEO sample index
(data/real/catalogues/geo_sample_index.tsv.gz), ClinicalTrials.gov posted results (aggregate only -> cannot serve).
- If >= 1 dataset with >= 20 patients and all three fields: model M0 = logistic(month-3 log UPCR), M1 = M0 + exposure;
  leave-one-out AUC; PASS if AUC(M1) - AUC(M0) >= 0.05 AND bootstrap 95% CI lower bound > 0. Trivial baselines also
  reported: majority class, baseline UPCR alone.
- If no dataset qualifies: Part A is recorded BLOCKED (not open), not failed.

## Part B — aggregate-level bound (fallback, explicitly weaker; label at most Emerging)
Binormal equal-variance model. A marker with standalone AUC A has separation d = sqrt(2) * Phi^-1(A).
Two markers with separations d1 (month-3 UPCR), d2 (exposure) and within-class correlation rho:
  d12^2 = (d1^2 + d2^2 - 2 rho d1 d2) / (1 - rho^2);  AUC12 = Phi(d12 / sqrt 2);  dAUC = AUC12 - Phi(d1/sqrt 2).
Inputs, extracted by a fixed rule: for each marker take the reported ROC AUC for a renal-response outcome from the
largest open-text LN cohort; if only a 2x2 table, AUC = (sens + spec)/2; if only medians/IQR, convert via
d = (median_R - median_NR) / pooled SD (SD = IQR/1.35) and record as lower-quality. Month-3 UPCR AUC: largest
published LN cohort reporting month-3 (or 12-week) proteinuria AUC for 12-month CR; if several, use the median.
Within-class rho between exposure and month-3 UPCR is not reported anywhere we know; scan rho in [-0.2, 0.6].
Decision (Part B):
- B-PASS ("plausible gate") if dAUC >= 0.05 for all rho in [0, 0.4] for at least one exposure marker.
- B-CONTINGENT if dAUC >= 0.05 only for some rho in [-0.2, 0.4].
- B-FAIL if dAUC < 0.05 for all rho in [-0.2, 0.6].
Caveat frozen now: the binormal model ignores mediation (exposure acting through month-3 UPCR makes conditional value
shrink, which is partly captured by rho > 0) and different outcomes/time points across studies; Part B cannot
pass the queue rule, only say whether Part A is worth pursuing.
