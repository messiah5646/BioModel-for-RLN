# PREREG Q25 (round 121): does blood IFN6 fall within patients only in MMF responders?

Frozen 2026-11-23 (scheduled run) BEFORE any follow-up-visit expression value was looked at.
Already seen (round 97, PREREG 181690a): baseline (first-sample) IFN6/PB4 per patient, sri4_response labels,
UPCR/SLEDAI per sample, design (130 MMF samples / 44 patients; 54 AZA / 20). NOT seen: IFN6 at any non-first visit.

## Decision this can change
LUMEN's IFN set-point layer (setpoint.py; Supported, n ~450, round 117) is response-independent: IFN_i(t) = own
baseline x therapy multiplier. KURDS-type claims say IFN falls only in responders.
- If H1 PASSES -> add a response-coupled term to the set-point layer (IFN falls with clinical response) and make
  on-treatment IFN an eligible early read-out for the forecaster's Task C' tests.
- If H1 FAILS -> keep the set point response-independent; record "IFN falls only in responders" as not supported
  in open data; IFN is not used as an on-treatment response read-out.

## Data and units
GSE224705 (GPL13158 whole blood), same parsing as run.py (round 97): highest-mean probe per gene,
IFN6 = mean per-gene z of IFI27, IFI44L, IFIT1, ISG15, RSAD2, SIGLEC1, z-scaled on ALL samples of the drug group.
(Deviation from r97, stated in advance: z reference = all samples of that drug group, because follow-up samples
are now used; the scaling is linear so it does not change within-patient differences up to a constant.)
Patient-level: first and LAST sample of each patient within the drug group; patients with 1 sample dropped.
Label: sri4_response as in round 97 (YES = responder). dt = days between first and last sample.

## Models (leave-one-patient-out, predict last-visit IFN6)
- C0 carry-forward: last = first.
- C1 set-point regression (trivial baseline to beat): last = a + b*first, fitted on the other patients.
- C2 response-coupled: last = a + b*first + c*R, fitted on the other patients.

## Hypotheses and pass rules
- **H1 (primary, MMF).** MSE(C2) <= 0.90 x MSE(C1) AND full-data c < 0 (responders end lower) AND
  Mann-Whitney on delta = last - first, R vs NR, two-sided p < 0.05. All three, else FAIL.
- **H2 (replication, AZA).** Same rule. Reported; n ~20, expected underpowered.
- **H3 (confound, reported, not a pass rule).** SRI-4 is defined by a SLEDAI fall. Add delta-SLEDAI to C2
  (C3: last = a + b*first + c*R + d*dSLEDAI) and report c with and without it; if c shrinks by > 50% the
  response coupling is "activity tracking", not response per se.
- Sensitivity (reported): restrict to dt >= 180 days; per-year slope instead of last - first.

## Expected
~10 MMF non-responders (r97), so wide CIs; the round-117 set point (ICC 0.87-0.93) predicts H1 FAIL.
Research only, not clinical guidance.
