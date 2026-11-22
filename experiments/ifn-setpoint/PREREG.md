# Pre-registration, round 96 (queue Q10): blood IFN as a therapy-lowered patient set point, joint over CAR-T and CD3xCD19
Frozen 2026-11-22, after reading sample labels only (obs: sample, donor, timepoint, dose). No per-donor IFN values
are in the repo; they are recomputed here for the first time since round 92.

## What is already known (disclosed, so it cannot be "predicted")
Direction counts only: CAR-T monocyte score fell in 7/7 (all-cell 6/7); CD3xCD19 monocyte score fell at M3 in 8/11,
at D28 in 6/11; CAR-T Spearman(pre-score, drop) = 0.32. Round 89 (GSE254176): patient identity explains 54% of variance.

## Data and score
- CD3xCD19 (CELLxGENE 3979a248 Monocytes): 11 donors x Base, D28, M3 (S14's third sample is M13; it is used as its
  late sample and flagged). Monocytes = IFI+, Classical, CD16+ Non-Classical; granulocytes excluded. Score = mean over
  cells of the mean CELLxGENE log-normalised expression of IFI27, IFI44L, IFIT1, ISG15, RSAD2, SIGLEC1.
- CAR-T (GSE263931): 7 donors x pre, post (scRNA samples only). Raw counts, cells > 500 UMI, normalise to 1e4 + log1p,
  monocyte gate LYZ > 2 (log scale; same crude gate as the ifn-cart round), same 6-gene score.
- The two datasets are not on a common scale, so every model works on the ratio post/base within a trial.

## Model (patient-state layer)
IFN_i(t) = S_i * f(t), S_i = the donor's own baseline set point, f(t) = the therapy-lowered fraction (f(0) = 1).
"Joint" = one f for "late after deep B-cell depletion" shared by CAR-T post and CD3xCD19 M3.

## Predictions for a held-out donor (leave-one-donor-out within the pooled 18 donors; parameters from the others only)
- B0 carry-forward: post = base.
- B1 population level: post = mean post score of the training donors in the same trial (ignores the set point).
- M1 set point x trial-specific fall: post = base x median(post/base) of training donors, same trial and timepoint.
- M2 set point x joint fall: post = base x median(post/base) over training donors from BOTH trials (late targets only).
Error = |predicted - observed| / observed baseline of that donor (relative absolute error, scale-free).
Targets: "late" = CAR-T post (7) + CD3xCD19 M3/M13 (11) = 18; "early" = CD3xCD19 D28 (11), M1 only vs B0/B1.

## Hypotheses and pass rules
- H1 set point persists: on late targets, M1 mean error < B1 mean error by >= 10% (relative), AND Spearman(base, late)
  >= 0.5 within each trial.
- H2 therapy lowers it: on late targets, M1 mean error < B0 mean error by >= 10%.
- H3 one shared fall across the two mechanisms: M2 mean error <= 1.10 x M1 mean error on late targets.
- H4 delay (descriptive, not a test, because the D28 direction count is already known): median f at D28 vs M3.
Each passes or fails on its own. n = 18 donors in two small trials: any pass is Emerging.
If H1 and H2 pass, IFN enters LUMEN as a patient-state set point with a therapy multiplier; if H1 fails, the set point is
not portable and IFN should be modelled as a population level.
