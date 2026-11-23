# PREREG P6 / Q56 (re-scoped by Q59) — round 118
Written 2026-11-23 (scheduled run) BEFORE any autoantibody time-course values were opened. Only paper titles/metadata
have been looked at so far. Frozen at this commit; amendments go in separate files.

## Decision this can change
Whether LUMEN gets a slow CD19-negative plasma-cell reservoir state (one extra parameter: the fraction of autoantibody
production that CD19/CD20-only depletion cannot reach).
- If the two-pool model wins (rules below) -> add the state to LUMEN's patient-state layer with the fitted plateau
  fraction, and use it to represent refractoriness/relapse after CD19/CD20-only depletion.
- If it does not win -> do not add it; LUMEN keeps one antibody pool; round-110 'refractory reservoir' stays Hypothetical
  and the preprint says the reservoir is not identifiable from open titre curves.

## Data (open only)
Curves = one patient's (or one published arm's) autoantibody titre (anti-dsDNA first; anti-Sm/anti-C1q/ANA titres
only if numeric) versus time after a B-lineage depleting therapy, taken from open full text, open supplements or
source-data files. Inclusion: baseline value + >= 3 post-treatment values, last >= 90 days. Numbers from tables or
source data only for the PRIMARY set; values read off figures (digitised) are allowed only in a sensitivity set,
flagged `digitised=1`.
Therapy classes: C1 = CD19/CD20-only (CD19 CAR-T, rituximab, obinutuzumab, CD19xCD3 engagers);
C2 = plasma-cell-hitting (daratumumab, CD19+BCMA or BCMA CAR-T).
Values reported as negative / below the assay limit are set to the reported lower limit (LOD) and every model's
prediction is floored at the same LOD before scoring, so an assay floor cannot fake a plateau.

## Models (y = log10(titre / baseline titre), t in days after first dose)
- M0 carry-forward: y_hat = 0 (titre stays at baseline). [Trivial baseline]
- M0b population-mean: y_hat(t) = mean y of the training curves at the nearest visit (+-30 d), else 0.
- M1 one pool:  titre(t)/titre0 = exp(-k t)                          (1 parameter, k >= 0)
- M2 two pool:  titre(t)/titre0 = f exp(-k t) + (1 - f)               (2 parameters, 0<=f<=1, k>=0)
  (slow pool treated as flat over the observed window; IgG turnover sets the fast decay.)
Fits: least squares on y with the LOD floor applied inside the model.

## Tests
T1 (primary, held-out): leave-one-curve-out. Shared parameters fitted on all other curves of the same class
(if the class has < 4 curves, pooled across classes); predict every post-baseline point of the held-out curve from its
baseline only. Score = mean absolute error in log10 units over all held-out points, curves weighted equally.
T2 (within-curve): per-curve AICc of M2 vs M1 (n points per curve; Gaussian on y). Summed dAICc and the count of
curves where M2 is better by >= 2.
T3 (mechanistic, directional): fitted plateau (1 - f) for C1 > C2.

## Pass rule (all must hold, computed on C1 = the refractoriness-relevant class, primary numeric set)
1. T1: MAE(M2) < MAE(M1) by >= 10% relative, AND MAE(M2) < MAE(M0) and < MAE(M0b).
2. T2: summed dAICc (M1 - M2) >= 4 AND M2 better by >= 2 in at least half of C1 curves.
3. n >= 8 C1 curves. If fewer than 8 curves qualify, the test is reported as UNDERPOWERED and the state is NOT added
   (decision B), whatever the point estimates say.
T3 is reported but does not gate the decision (C2 is expected to be tiny).

## Q59 secondary (descriptive only, no pass rule)
Where the same source reports relapse/flare per patient after C1 therapy, compare the fitted per-curve plateau in
relapsing vs non-relapsing patients. No claim above Emerging.

## Labels
Supported only if pass and n >= 10 curves from >= 2 independent sources; otherwise Emerging at most.
