# PREREG Q16 (round 106): does baseline chronicity cap the proteinuria response (floor) or slow it (rate)?

Committed BEFORE any outcome numbers (proteinuria-by-chronicity time courses) are extracted. Written 2026-11-22.

## Question
Round 93 proposed a "damage floor" route to non-response: high baseline NIH chronicity index (CI) leaves a residual
proteinuria plateau that immunosuppression cannot remove. Alternative: chronicity slows the fall (rate) but the
eventual level is the same. Six-state-model reading: damage D sets a floor term in UPCR; activity sets the decaying term.

## Data eligibility (open only; Open-data policy 2026-11-22)
Any open source (PMC full text, open supplement, GEO characteristics, ClinicalTrials.gov) reporting proteinuria
(UPCR or 24h protein, mean or median, absolute or % of baseline) at >= 3 post-baseline visits spanning >= 12 months,
separately for a high- and a low-chronicity group (CI >= 3 or 4 vs below, or tertiles, or authors' split), in treated LN.
Patient-level open data with baseline CI and >= 3 proteinuria visits also qualify (preferred).
Minimum: 1 eligible dataset to run; otherwise BLOCKED (not failed).
Secondary evidence (not scored for the shape test): response rates (CR/PR at 2 time points) by CI group. With 2 time points
the shape test is impossible; reported descriptively only.

## Models (fit to R_g(t) = P_g(t)/P_g(0), group g in {low, high}; t in months; least squares on log R if medians, R if means)
- M0 pooled:   R(t) = f + (1-f) exp(-k t)                            (2 params; trivial "chronicity does not matter")
- MF floor:    R_g(t) = f_g + (1-f_g) exp(-k t)                      (3 params: k, f_low, f_high)
- MK rate:     R_g(t) = f + (1-f) exp(-k_g t)                        (3 params: f, k_low, k_high)
- MB both:     group-specific f and k                                 (4 params; reported only)
Bounds: 0 <= f < 1, 0.01 <= k <= 5 /month.

## Decision rules (per dataset; then majority across datasets)
1. Chronicity matters: best of MF/MK beats M0 by AICc >= 2. Else verdict "no chronicity effect on shape".
2. FLOOR supported if AICc(MF) <= AICc(MK) - 2 AND f_high > f_low.
   RATE supported if AICc(MK) <= AICc(MF) - 2 AND k_high < k_low.
   Otherwise "not separable" (equal parameter counts, so AICc difference = SSE ratio).
3. Held-out check (trivial baselines): refit MF and MK on visits <= 6 months (>= 2 post-baseline points required);
   predict the last reported visit per group. Baselines: B1 carry-forward of each group's last fitted R;
   B2 pooled M0 fit on the same early visits. The winning shape must beat both B1 and B2 in mean absolute error
   on held-out visits, else verdict downgraded to "fits but does not predict".
4. Evidence label: one dataset, aggregate = Emerging at most; >= 2 concordant datasets incl. patient-level = Supported.

## Expected answer (stated first)
FLOOR (H1): high-CI group plateaus higher, similar rate. Confidence low (~55%); aggregate medians from small cohorts
are noisy and selection (non-responders re-biopsied/switched) can flatten late curves.

## Known confounds to report
High-CI patients often have lower eGFR, longer disease duration, more prior therapy; treatment differs across cohorts.
Floor vs rate in aggregate curves is also produced by mixing (a fraction of non-responders), which this test cannot
separate from a per-patient floor. Stated, not fixed.
