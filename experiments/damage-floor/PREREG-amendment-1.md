# PREREG Q16 amendment 1 (committed before any outcome extraction; 2026-11-22)

Reason: a first PubMed pass (titles/abstract metadata only) found no open paper with multi-visit proteinuria by chronicity
group. Time-to-response curves (Kaplan-Meier / cumulative incidence of complete or partial renal response) stratified by
chronicity carry the same floor-vs-rate question and are added as a second eligible data type. The original R(t) rules
are unchanged.

## Added data type: cumulative response C_g(t) (fraction responded by t), >= 3 read time points spanning >= 12 months
(digitised from a published open-access KM figure or read from a table), groups high vs low chronicity.
- M0c pooled: C(t) = c (1 - exp(-k t))
- MFc floor (cure fraction differs): C_g(t) = c_g (1 - exp(-k t))
- MKc rate: C_g(t) = c (1 - exp(-k_g t))
Same AICc rules 1-2 (FLOOR if MFc wins by >=2 AND c_high < c_low; RATE if MKc wins by >=2 AND k_high < k_low)
and same held-out rule 3 (fit t <= 6 mo, predict last read time; baselines B1 carry-forward, B2 pooled).
Digitisation: read at fixed months {3,6,9,12,18,24,36,48,60} where curves exist; reading error assumed +/-0.02.
A HR-only result (Cox, no curve) does NOT qualify: a HR cannot separate floor from rate.
Two-time-point response rates (e.g. CR at 6 and 12 mo by CI group): descriptive only, as in the original PREREG.
