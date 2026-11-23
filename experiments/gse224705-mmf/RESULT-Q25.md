# RESULT round 121 (Q25): blood IFN6 does not fall within patients in MMF (or AZA) responders

PREREG: PREREG-Q25-trajectory.md (a120f82), download.sh (c85f59a), analysis q25_trajectory.py. Output q25_result.json.
Patients with >= 2 samples; first vs last sample; median interval 109 d (MMF), 104 d (AZA).

| | MMF (n 43; 34 R / 9 NR) | AZA (n 20; 11 R / 9 NR) |
|---|---|---|
| MSE C0 carry-forward | 0.099 | 0.081 |
| MSE C1 set-point regression (last = a + b*first) | **0.063** (b 0.82) | 0.102 (b 1.00) |
| MSE C2 + responder term | 0.067 (ratio 1.05; rule <= 0.90) | 0.107 (ratio 1.05) |
| responder coefficient c (full fit) | +0.07 (wrong sign; rule < 0) | +0.23 (wrong sign) |
| median delta IFN6, R vs NR | +0.03 vs +0.00, MW p 0.99 | +0.14 vs -0.01, MW p 0.29 |
| H1 / H2 | **FAIL** (all 3 conditions) | **FAIL** |

- H3: adding delta-SLEDAI leaves c unchanged (MMF 0.072 -> 0.071); there is no response coupling to explain away.
- Sensitivity dt >= 180 d: too few non-responders (MMF 2, AZA 3); per-year slope MW p 0.89 / 0.36.
- Exploratory (not pre-registered, q25_exploratory_ci.txt): MMF bootstrap 95% CI for c = -0.12 to +0.30, i.e. a
  responder-specific fall larger than ~0.11 between-patient SD is excluded; first-last correlation 0.97 (MMF), 0.96 (AZA).
- Limits: response is SRI-4 (global), not renal; first sample is often already on drug; short median interval
  (109 d); 9 MMF non-responders; single centre. Deviation stated in PREREG: z reference = all samples of the drug group.
