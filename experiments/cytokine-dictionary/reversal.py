"""Exploratory reversal screen (hypothesis-generating, not a test).
Flare signature per cell type = mean of the unit-normalised responses to IFN-beta, IFN-gamma and IL-15, the three
signatures the screen's authors report as enriched in lupus flares. For every other cytokine q and cell type c:
  cosine   direction of q's response against the flare signature (negative = opposing)
  cancel   -<R(q,c), F_c> / <F_c, F_c> : share of the flare signature that q's response would cancel if effects added
Healthy donor cells, 24 h, one dose. Opposition here is not evidence of benefit in patients."""
import sys, json, numpy as np, pandas as pd
raw = pd.read_csv(sys.argv[1], index_col=0); d = raw[(raw.adj_p_value < 0.05) & (raw.log_fc.abs() > 0.5)]
CELLS = ["CD14_Mono", "CD16_Mono", "cDC", "pDC", "Naive_B_cell", "Intermediate_B_cell", "Plasmablast", "CD4_Naive_T_cell", "CD4_Memory_T_cell", "CD8_Naive_T_cell", "CD8_Memory_T_cell", "Treg", "NK_CD56low", "NK_CD56hi"]
FLARE = ["IFN-beta", "IFN-gamma", "IL-15"]; rows = []
print("cell types present:", sorted(d.celltype.unique()))
for c in CELLS:
    sub = d[d.celltype == c]
    if sub.empty: continue
    P = sub.pivot_table(index="gene", columns="cytokine", values="log_fc", aggfunc="first").fillna(0.0)
    if not all(f in P.columns for f in FLARE): continue
    F = sum(P[f] / np.linalg.norm(P[f]) for f in FLARE) / 3; scale = np.mean([np.linalg.norm(P[f]) for f in FLARE]); Fs = F * scale
    core = F[F > 0].sort_values(ascending=False).index[:100]                 # 100 genes most strongly up in the flare signature
    for q in P.columns:
        y = P[q]; n = int((y != 0).sum())
        if n < 20: continue
        rows.append(dict(cytokine=q, cell=c, genes=n, cosine=float(y @ F / (np.linalg.norm(y) * np.linalg.norm(F))), cancel=float(-(y @ Fs) / (Fs @ Fs)),
                         core_down=int((y[core] < 0).sum()), core_up=int((y[core] > 0).sum())))
R = pd.DataFrame(rows); R.to_csv("reversal.csv", index=False)
X = R[~R.cytokine.str.startswith("IFN") & (R.cytokine != "IL-15")]
agg = X.groupby("cytokine").agg(cells=("cell", "size"), opposing_cells=("cosine", lambda v: int((v < -0.05).sum())), aligned_cells=("cosine", lambda v: int((v > 0.05).sum())),
                                mean_cos=("cosine", "mean"), min_cos=("cosine", "min"), mean_cancel=("cancel", "mean"), core_down=("core_down", "sum"), core_up=("core_up", "sum"))
pd.set_option("display.width", 220)
print("== strongest opposers (by mean cancel share)"); print(agg.sort_values("mean_cancel", ascending=False).head(14).round(3).to_string())
print("== strongest amplifiers"); print(agg.sort_values("mean_cancel").head(12).round(3).to_string())
piv = X.pivot(index="cytokine", columns="cell", values="cancel").round(2)
top = list(agg.sort_values("mean_cancel", ascending=False).head(8).index) + list(agg.sort_values("mean_cancel").head(8).index)
print("== cancel share by cell type (positive = opposes flare signature)"); print(piv.loc[top, [c for c in CELLS if c in piv.columns]].to_string())
