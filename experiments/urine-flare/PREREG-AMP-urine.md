# Pre-registration: urine CD163 / MRC1 vs UPCR for 12-month response, AMP Phase II urine proteomics (frozen 2026-11-22)
Written before any access to syn53124971 (access not yet granted; nobody on this project has seen the patient-level file).

Data: ph2_urine_proteomics.txt + ph2_urine_proteomics_basicclinicodemographic.txt (Synapse syn53124971).
Patients: proliferative or mixed LN with a baseline (biopsy) and a month-3 urine sample and a 12-month response label.
Outcome: complete response (CR) at 12 months vs not (as labelled in the file; partial responders counted as not-CR).

Features (urine, as given, log2; normalised to urine creatinine if the file provides it, else as-is):
- M0: month-3 UPCR alone (trivial baseline 1).
- M1: month-3 UPCR + month-3 urine CD163.
- M2: month-3 UPCR + month-3 urine CD163 + MRC1 (CD206).
- Trivial baseline 2: percent fall in UPCR from baseline to month 3.
Model: logistic regression, no tuning; 5-fold stratified CV repeated 20x; report mean AUC and 95% bootstrap CI.

Pass rule (H1): M1 AUC exceeds max(M0, trivial 2) by >= 0.05 AND the bootstrap 95% CI of the paired AUC difference excludes 0.
H2: M2 exceeds M1 by >= 0.03 (otherwise MRC1 adds nothing).
H3 (flare, only if the file carries flare dates): month-of-flare uCD163 > preceding sample in >= 75% of flare episodes,
and uCD163 rises in the sample before the flare in >= 75% of episodes with a pre-flare sample.
If n < 30 patients with both samples, the result is labelled Emerging whatever the outcome.
