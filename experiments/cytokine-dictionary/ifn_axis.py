"""T8 interferon decomposition. For each cytokine q and cell type c, fit R(q,c) = a*R(IFN-beta,c) + b*R(IFN-gamma,c), a,b >= 0.
Criterion fixed before running: across non-interferon cytokines with >= 200 changed genes in CD14 monocytes, the IFN-gamma
coefficient b must track induction of the IFNG gene (largest significant log2FC in any cell type, 0 if none): Spearman > 0.5."""
import sys, json, numpy as np, pandas as pd
from scipy.optimize import nnls
from scipy.stats import spearmanr
raw = pd.read_csv(sys.argv[1], index_col=0); d = raw[(raw.adj_p_value < 0.05) & (raw.log_fc.abs() > 0.5)]
IFNS = [q for q in d.cytokine.unique() if q.startswith("IFN")]
def decomp(cell):
    sub = d[d.celltype == cell]; P = sub.pivot_table(index="gene", columns="cytokine", values="log_fc", aggfunc="first").fillna(0.0)
    A = P[["IFN-beta", "IFN-gamma"]].values; rows = []
    for q in P.columns:
        y = P[q].values; n = int((y != 0).sum())
        if n < 20: continue
        (a, b), _ = nnls(A, y); e = lambda v: float((v ** 2).sum())
        tot = e(y); fit = 1 - e(y - A @ [a, b]) / tot
        only_b = 1 - e(y - A[:, 1] * max(float(A[:, 1] @ y / (A[:, 1] @ A[:, 1])), 0)) / tot
        only_a = 1 - e(y - A[:, 0] * max(float(A[:, 0] @ y / (A[:, 0] @ A[:, 0])), 0)) / tot
        rows.append(dict(cytokine=q, cell=cell, genes=n, a_typeI=round(a, 3), b_gamma=round(b, 3), explained=round(fit, 3),
                         gain_from_gamma=round(fit - only_a, 3), gain_from_typeI=round(fit - only_b, 3)))
    return pd.DataFrame(rows)
ifng = d[d.gene == "IFNG"].groupby("cytokine").log_fc.max(); t1 = d[d.gene.str.match(r"^IFN(A\d+|B1|W1|E|K)$")].groupby("cytokine").log_fc.max()
out = {}
for cell in ("CD14_Mono", "CD16_Mono", "cDC", "Naive_B_cell", "CD4_Naive_T_cell", "NK_CD56low"):
    D = decomp(cell); D["IFNG_log2fc"] = D.cytokine.map(ifng).fillna(0).clip(lower=0).round(2); D["typeI_gene_log2fc"] = D.cytokine.map(t1).fillna(0).clip(lower=0).round(2)
    S = D[(D.genes >= 200) & ~D.cytokine.isin(IFNS)]
    rho = spearmanr(S.b_gamma, S.IFNG_log2fc)[0] if len(S) > 3 else np.nan
    rho2 = spearmanr(S.gain_from_gamma, S.IFNG_log2fc)[0] if len(S) > 3 else np.nan
    out[cell] = dict(n=len(S), spearman_b_vs_IFNG=None if np.isnan(rho) else round(float(rho), 2), spearman_gammagain_vs_IFNG=None if np.isnan(rho2) else round(float(rho2), 2))
    print("==", cell, out[cell])
    if cell in ("CD14_Mono", "Naive_B_cell"):
        print(S.sort_values("explained", ascending=False)[["cytokine", "genes", "a_typeI", "b_gamma", "explained", "gain_from_gamma", "gain_from_typeI", "IFNG_log2fc", "typeI_gene_log2fc"]].head(22).to_string(index=False))
    D.to_csv(f"ifn_axis_{cell}.csv", index=False)
json.dump(out, open("ifn_axis.json", "w"), indent=1)
