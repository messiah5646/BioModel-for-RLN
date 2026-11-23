# PREREG P5 (= Q32), round 117: is blood IFN a patient set point? ILLUMINATE-1/2 (tabalumab vs placebo)

Committed BEFORE the expression files are downloaded (rule 12a push 1).

**Decision this can change:** keep or drop LUMEN's IFN set-point layer (round 96: IFN_i(t) = own baseline x shared
therapy multiplier). If the set point holds in placebo and the set-point model beats the population (arm-mean) model
out of trial -> keep the layer. If persistence fails -> drop it and model blood IFN as a population/arm-level quantity.

## What was seen before this PREREG
GEO metadata only (GSE88885 ILLUMINATE-1: 908 samples, 309 subjects; GSE88886 ILLUMINATE-2: 418 samples, 147 subjects;
baseline / week 16 / week 52; placebo, tabalumab 120 mg Q2W, Q4W; baseline SLEDAI, anti-dsDNA, C3/C4; platform GPL17586).
No expression values. The analyst knows tabalumab's clinical results (ILLUMINATE-1 missed, ILLUMINATE-2 met SRI-5) but not
these IFN data.

## Why the round-96 multiplier (0.73) is not used as a number
It was measured on single-cell log1p(CP10K) in monocytes; GEO here is whole-blood microarray log2. A ratio of scores on
different log scales is not comparable. This test checks the set-point STRUCTURE; the multiplier is re-estimated.

## Score (fixed now)
IFN6 = per sample, mean over IFI27, IFI44L, IFIT1, ISG15, RSAD2, SIGLEC1 of (log2 expression - that gene's median log2
expression across all baseline samples of the same trial). Probesets mapped to genes by the GPL17586 annotation (or the
file's own gene column); several probesets per gene -> mean. If a gene is absent, the score uses the rest (reported).

## Primary analyses (patients with all needed visits)
**A. Persistence in placebo (set-point core).** Spearman(IFN6 baseline, IFN6 week 52), placebo arm, each trial.
- Passes if rho >= 0.5 in BOTH trials. Fails if rho < 0.5 in both. Mixed -> inconclusive.
- Reported with it: intraclass correlation (patient share of variance) over baseline/w16/w52 in placebo; round 89 said
  54% (n 6), Q31's band 39-69%.

**B. Out-of-trial prediction of week-52 IFN6** (train on one trial, test on the other, both directions; all arms).
- M0 arm mean: mean week-52 IFN6 of the same arm in the training trial (population model).
- M1 carry-forward: own baseline.
- M2 set point: own baseline + mean (w52 - baseline) change of the same arm in the training trial.
- M3 shrunk set point: OLS w52 ~ baseline + arm, fitted on the training trial.
- Score: MSE on the test trial. **Set-point model wins** if min(M2, M3) beats M0 by >= 10% in BOTH directions.
  **Loses** if M0 is better than both M2 and M3 in both directions.

**Decision rule:** keep the layer if A passes and B wins. Drop it if A fails. Anything else: keep as Emerging, no change.

## Secondary (reported, not decisive)
- C. Arm effect: tabalumab minus placebo change (w52 - baseline) in IFN6, each trial, with 95% CI. Direction predicted by
  the round-96 logic (partial B-lineage hit): tabalumab change < placebo change, but smaller than after B-cell depletion
  (not quantified across scales).
- D. Week-16 vs week-52: does the w16 value predict w52 better than baseline? (Early-read-out question for IFN.)
- E. Does baseline anti-dsDNA or C3 predict the IFN6 change (round-113 follow-up at population scale)? Spearman, each trial.
- QC: none available per sample beyond the preprocessed values; batch = trial (handled by within-trial gene medians).

## Labels
n ~ 450 patients: a pass or fail here can be labelled Supported (pre-registered, n >= 10, two independent trials).
Not blind to the trial's clinical results.

## Not allowed after the download
Changing genes, score, visit choice, thresholds or models.
