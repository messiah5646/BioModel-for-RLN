# PREREG amendment (Q45, round 104): decision rule 3 must match on ALL B0 metrics

Frozen 2026-11-22, before C5 (Q37) expression is loaded. Does NOT edit earlier frozen PREREGs; applies to C5 onward.
Reason: in C4 (round 102), matching on log counts, log genes and mito only left B0 at AUC 0.93 because ribosomal fraction
(not matched) itself separates SLE from controls.

Rule 3 (replaces the matching set for future contrasts): nearest-neighbour matching within cell type uses every per-cell B0
metric that is not the hypothesis: log counts, log genes, top-20 %, mito, ribo, stress, HB, platelet and IG fractions
(standardised; caliper 0.5 SD x sqrt(number of metrics)). If B0 is still > 0.6 AUC after matching, the result is reported as
"QC metrics carry biology (or batch) that matching cannot remove", rule 3 is recorded as FAILED, and no latent-space gain is
claimed for that contrast.
