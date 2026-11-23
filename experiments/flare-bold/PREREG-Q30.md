# PREREG Q30 (round 122): does a baseline blood plasmablast score predict time to flare in BOLD (GSE92776) beyond baseline SLEDAI?

Written 2026-11-23 before any expression value of GSE92776 was downloaded or read.
Already seen (sample characteristics in data/real/catalogues/geo_sample_index.tsv.gz): group A = longitudinal BOLD patients
(~41 with a BASELINE sample; ~34 have a FLARE visit, ~6 do not), group B = 58 baseline-only patients, group C = 99 healthy.
Visit dates and SLEDAI per sample are in the characteristics. No expression, no score, no flare timing computed.

**Decision this can change:** whether LUMEN's planned flare/relapse layer (Q9 line; maintenance, withdrawal and refractory
switching simulations) takes a blood B-lineage (plasmablast) covariate or stays clinical-only.
If H1 passes -> plasmablast score becomes an Emerging flare-hazard covariate candidate, to be replicated in GSE121239 before it
is wired in. If H1 fails -> the flare layer stays clinical-only and the D1/D4 "plasmablasts = activity" claim is not used for
forecasting.

Rule 8 note: blood IFN as a flare marker is settled (rounds 89, 99, 104). IFN6 is computed here only as a reported comparator;
it carries no hypothesis and no decision.

## Data
GSE92776 series matrix (GEO, open). Patients: subject_type=Patient, group_id=A, with a BASELINE sample.
t0 = baseline visit date. Event = first FLARE visit (nicevisit=FLARE) date after t0. Non-flarers censored at their last visit date.
Patients with a flare but no baseline sample are excluded.

## Scores (frozen, same genes as round 97, experiments/gse224705-mmf/run.py)
- PB4 = mean z of JCHAIN (or IGJ), MZB1, TNFRSF17, TXNDC5 (whichever present; probes collapsed to gene by mean; log2 scale as deposited,
  log2(x+1) if max > 100). z computed across the analysed baseline samples.
- IFN6 = mean z of IFI27, IFI44L, IFIT1, ISG15, RSAD2, SIGLEC1 (comparator only).

## Models (Cox proportional hazards, lifelines, penalizer 0.01)
- M0 null (C = 0.5); M1 baseline SLEDAI; M2 SLEDAI + PB4; M3 SLEDAI + IFN6 (comparator).
- Score = leave-one-patient-out Harrell C-index of the out-of-fold linear predictor (refit without the patient, predict the patient).

## H1 (primary, decides)
PASS if LOO C(M2) - LOO C(M1) >= 0.05 AND the bootstrap (2,000 patient resamples, seed 122) 95% CI of the M2 PB4 log-HR excludes 0
with HR > 1 (higher plasmablast score -> earlier flare). Otherwise FAIL. If fewer than 20 events, report UNDERPOWERED instead.

## Robustness (reported, cannot rescue a FAIL)
R1 add baseline lymphocyte % (where recorded) to M1 and M2; R2 PB4 alone vs null; R3 Spearman of PB4 with days-to-flare among flarers.

## H2 (secondary, descriptive for the D1 claim)
Within patient, PB4 at the FLARE visit minus PB4 at baseline (paired Wilcoxon, two-sided) and FLARE minus IMPROVING. PB4 here is
z-scored across all group-A samples. Reported with n; no decision.

## Trivial baselines
Null (C 0.5) and SLEDAI-only, as above. Report every number including failures.
