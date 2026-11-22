# PREREG amendment (Q43, round 104): "copy the nearest measured context" baseline

Frozen 2026-11-22, before C5 (Q37) expression is loaded. Does NOT edit PREREG.md, PREREG-C1-round100.md or PREREG-C3C4-round102.md.
Motivation: Virtual Cell Challenge 2026 winners (round 101) mostly transfer a measured effect from the nearest context.

1. New baseline **B5 (nearest-context copy)** for every future classification contrast in this protocol (C5 onward, and any rerun
   of C1-C4 that is labelled as a rerun): the score of a test sample = mean training label of its k = 3 nearest training samples,
   Euclidean distance on standardised B2 (CLR composition) features. Standardisation fitted on training samples only.
   B5 joins B1-B4 in "best simple baseline" for decision rule 1.
2. cell-eval2 scaling (0 = context mean, 1 = replicate ceiling) applies to RESPONSE-PREDICTION tasks (predicting a post-treatment
   profile), not to label classification. It is adopted for any future response-prediction contrast in this protocol; for
   classification contrasts it is not applicable and is not reported. Declared now so it cannot be chosen after results.
