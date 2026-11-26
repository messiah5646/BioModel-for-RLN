# RESULT VC-1 (round 142, corrected round 143): the learned virtual-cell predictor FAILS; other patients' month-3 average beats own baseline (library-confounded)

PREREG-VC1 0ea016d (04:26:42Z) -> amendment A1 674d5fa (04:46:39Z, synthetic code checks only) -> real run started
05:27:33Z (60 min 51 s after the PREREG), finished 06:37:25Z. Code as run = run_vc1.py at 674d5fa. The only later change
adds the stamp.py call it had omitted; it was stamped after the run, which is disclosed in the JSON.
Label: **12b-timed (data on disk before the PREREG; earlier looks not excludable)**. Emerging, n 11 patients.

## Verdict: FAIL (rule 2b)
Leave-one-patient-out, 11 SLE patients, 15-23 scored cell types each (36 cell types; 2,000 Base-chosen HVGs).
Primary metric: MSE of predicted vs observed month-3 pseudobulk.

| Model | mean MSE | vs best trivial (B2) |
|---|---|---|
| B0 carry-forward (own baseline) | 0.0535 | 1.67x |
| B1 baseline + average shift | 0.0584 | 1.82x |
| **B2 population-average M3** | **0.0321** | 1.00 |
| B3 shrink (fitted a) | 0.0332 | 1.03x |
| J1 linear latent (residual, A1) | 0.0493 | 1.54x, wins 0/11 |
| **J2 LUMEN-JEPA MLP, augmented (A1)** | **0.0423** | **1.32x, wins 0/11 -> FAIL** |
| J2 without augmentation | 0.0443 | 1.38x |
| J2orig (PREREG-as-written, S10) | 0.0392 | 1.22x (also FAIL; better than A1 form on real data) |
| noise ceiling (half-split) | 0.0051 | |

- **S7:** patient-bootstrap 95% interval of J2/B2 is 1.20 to 1.44, so the loss is not noise.
- **S1 delta correlation** (median over units): B3 0.66 > J2orig 0.54 > J1 0.50 > J2 0.49 > B1 0.21. VCC-style score
  (0 = mean shift, 1 = noise ceiling): B3 0.47, J2 0.30, J1 0.17.
- **S4:** the ordering B3/B2 < J2 < B0 < B1 holds in every compartment (T, B, monocytes, myeloid, LDG).
- **Augmentation:** it buys about 5% (J2 0.0423 vs 0.0443 without), not enough to matter.
- **Amendment A1 disclosure:** the residual-over-mean-shift form passed the synthetic controls but scored *worse* on real
  data than the as-written form. The reason: the average shift (B1) is itself a bad predictor here. Neither form changes
  the verdict.
- **Not run:** S2 composition and S3 (D28 target) were not computed in this run. S6 STAR-T was skipped (reason in A1).
- **Positive control (S11):** on synthetic data with a planted baseline-dependent response, J2 beat the best trivial
  rule by 16%. On a planted shared shift it stayed within 7% of the mean shift. So the model can find a personal signal when one
  exists; at n 11 the real data contain none that beats the population average.

## Decision applied (PREREG section 2)
FAIL -> LUMEN's cell-state layer, if one is built, is the **population-average post-therapy state (B2)** or shrink (B3),
not a learned model. No deep model is shipped. Learned latent dynamics add nothing at n 11, the same lesson as VCC 2026 (r101).

## What the data say (descriptive, post hoc; corrected in round 143 under reviewer R-20261126-1)
Round 142 first read this as "month-3 state is set by the therapy more than by the patient" and "converges to a shared
post-therapy profile". The reviewer was right that the numbers did not support that reading as written. Round 143
(`vc1_descriptive.py`, `vc1_descriptive.json`, re-fetched h5ads with md5 matching the PREREG, B0-B3 reproduce exactly)
measured it:

| LOPO, n 11 patients (n 10 without P9) | B0 own Base | B2base others' mean Base | B2 others' mean target | B3 |
|---|---|---|---|---|
| target M3 | 0.0535 | 0.0356 | **0.0321** | 0.0332 |
| target M3, without P9 | 0.0554 | 0.0360 | **0.0318** | 0.0327 |
| target D28 (S3, pre-registered, not run in r142) | 0.0451 | 0.0385 | 0.0326 | **0.0285** |
| target D28, without P9 | 0.0452 | 0.0389 | 0.0327 | **0.0294** |

- **What was measured:** the other patients' month-3 average predicts a patient's month-3 pseudobulk better than the
  patient's own baseline does (library-confounded: every patient x timepoint is its own library).
- **Most of that win is pooling, not a post-therapy state.** Replacing the patient's own Base library with the other
  patients' *Base* average already cuts the error from 0.0535 to 0.0356. Using their *M3* average instead cuts it a
  further 10% (B2base/B2 = 1.11, patient-bootstrap 95% 1.04-1.18; 8 of 11 patients; without P9 1.13, 1.02-1.25).
- **A small visit-specific shift exists at both visits.** At D28, B2base/B2 = 1.18 (1.07-1.31). It cannot be separated
  from visit-correlated technical effects (libraries of the same visit may share processing).
- **Dispersion (committed code; replaces the uncommitted 0.0204 -> 0.0163 numbers).** Noise-corrected between-patient
  dispersion, same patients at both visits, median over cell types: Base 0.0251 -> M3 0.0182 (lower in 22 of 26; without
  P9 18 of 25); D28 0.0211 -> M3 0.0200 (20 of 26; without P9 21 of 26). The half-split correction removes cell-sampling
  noise, not library effects.
- **Contrast with round 96, same 11 patients:** for the 6-gene monocyte IFN score, the post-treatment population average
  was the *worst* predictor (0.914 vs carry-forward 0.410). So a pooled profile beating carry-forward on 2,000 HVGs does
  not carry over to the one programme LUMEN has tested.
- **P9** (`S14_M13` M3 library): excluding it changes no conclusion above.

**Reading:** by the reviewer's pre-stated rule (B2 beats B2base by >= 10% and the dispersion drop holds without P9 and at
D28 -> M3) "convergence" is kept only as **Hypothetical (n 11, library-confounded, post hoc)**; the margin is narrow
(1.11 vs the 1.10 bar). The main reason carry-forward loses is that patient-specific pseudobulk signal largely does not
persist across two libraries at n 11. One therapy only, so "therapy-specific" is untested.

## Labels
- **Emerging (n 11 patients, 12b-timed): FAIL** for the learned virtual-cell predictor (J2 1.32x worse than the population
  mean, 0/11 wins).
- **Hypothetical (n 11, post hoc, library-confounded):** a small visit-specific population shift after CD3xCD19 therapy
  (B2 beats the pooled baseline average by 11% at M3, 18% at D28); between-patient dispersion modestly lower at M3.
- **Hypothetical:** a learned predictor might help at n >= 50 patients, or with within-library multiplexed designs.
  Nothing here supports it yet.
- LUMEN fallback cell-state rule (if built): **population average at the target visit**; whether it is therapy-specific
  is not shown.

Research and simulation only. Not clinical or drug-development guidance.
