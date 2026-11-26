# RESULT VC-1 (round 142): the learned virtual-cell predictor FAILS; the population-average month-3 state beats everything

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

## What the data say (descriptive, post hoc reading of pre-registered baselines; Emerging, n 11)
- **Month-3 state is set by the therapy more than by the patient.** A patient's own baseline (B0) predicts their month-3
  blood transcriptome 1.7x worse than the average month-3 profile of the other 10 patients does. The average
  Base->M3 shift (B1) is worse still, so patients do not move in parallel. They move toward a common post-therapy state.
- **Between-patient spread shrinks modestly by month 3.** Noise-corrected dispersion has median 0.0204 at Base and 0.0163
  at M3, lower in 20 of 25 cell types. This matches round 119: blood IFN converges after BCMA/CD19 CAR-T rather than scaling.
- **Confound that cannot be removed here:** each patient x timepoint is its own library (`sample` S01_Base, S01_M3, ...).
  Per-library technical effects therefore count against carry-forward (two libraries) more than against the
  population mean (averaged over 10). Part of the convergence may be technical.
- **Data-quality flag:** patient P9's month-3 library is named `S14_M13` in `obs.sample`, while `timepoint` says M3.
  It may be a month-13 sample. P9 and P1 are the only two patients where B0 beat B2.

## Labels
- **Emerging (n 11 patients, 12b-timed): FAIL** for the learned virtual-cell predictor (J2 1.32x worse than the population
  mean, 0/11 wins).
- **Emerging (n 11, descriptive):** after CD3xCD19 therapy, month-3 blood cell states are closer to a shared post-therapy
  state than to each patient's own baseline. This is confounded by per-library effects.
- **Hypothetical:** a learned predictor might help at n >= 50 patients, or with within-library multiplexed designs.
  Nothing here supports it yet.

Research and simulation only. Not clinical or drug-development guidance.
