# Pre-registration, round 89: what changes in blood before a lupus flare (GSE254176), 2026-11-21
Written after reading sample labels only (labels.csv); no expression values loaded.

Data: 19 PBMC scRNA samples, 6 SLE patients (LN cohort), labels Before flare / Flare / Remission.
Episodes Before->Flare (consecutive): P1 1->2, P2 1->2, P3 1->2, P3 4->5, P4 2->3, P4 5->6 (6 episodes).
Run-up pairs Remission->Before flare: P3 3->4, P4 1->2, P4 4->5 (3 pairs).

Scores (same as ifn-cart): monocyte 6-gene IFN-I score (IFI27, IFI44L, IFIT1, ISG15, RSAD2, SIGLEC1; LYZ-high gate);
plasmablast fraction = cells with JCHAIN and MZB1 both above zero and CD3E zero, over all QC cells (>500 UMI).

H1: IFN score Flare > Before-flare in >= 5 of 6 episodes.
H2: IFN score Before-flare > preceding Remission in 3 of 3 run-up pairs (early warning).
H3: plasmablast fraction Flare > Before-flare in >= 5 of 6.
H4: plasmablast fraction Before-flare > preceding Remission in 3 of 3.
Each passes or fails on its own. n is tiny: a pass is Emerging evidence, never more.
