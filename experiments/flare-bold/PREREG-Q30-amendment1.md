# PREREG Q30 amendment 1 (committed after download, before any score or outcome relation was computed)

Seen on download (platform + data_processing text only): GPL22843 is a ~300-gene TaqMan OpenArray qPCR panel. Values are
batch-corrected delta Ct (target Ct minus median of endogenous controls), so HIGHER value = LOWER expression.

Changes, fixed now:
1. Expression = -deltaCt for every gene (no log transform). Multiple assays for one gene are averaged.
2. PB4 genes present on the panel: TNFRSF17 (BCMA) only. JCHAIN/IGJ, MZB1, TXNDC5 are absent. Per the frozen rule ("whichever present"),
   the primary plasmablast score is z(-dCt TNFRSF17), a single-gene score. This is a weaker read-out than PB4; the round doc must say so.
   No substitute genes (CD38, CD27) are added to the primary score.
3. Missing values (undetermined assays): if TNFRSF17 is missing in >20% of analysed baseline samples, report UNDERPOWERED-DATA and stop
   H1; otherwise drop patients with a missing value and report n.
4. IFN6 comparator: all six genes are on the panel; unchanged.
Everything else in PREREG-Q30.md is unchanged.
