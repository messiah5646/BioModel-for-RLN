"""Criterion C1 on real data, as far as this table allows: does behavioural distance between cell types,
measured on one half of the cytokines, predict how differently they respond to the other half?"""
import sys, hashlib, numpy as np, pandas as pd
from scipy.stats import spearmanr
MYE = ["CD14_Mono", "CD16_Mono", "cDC"]
LYM = ["CD4_Naive_T_cell", "CD8_Naive_T_cell", "CD4_Memory_T_cell", "CD8_Memory_T_cell", "Treg", "MAIT", "NK_CD56low", "NK_CD56hi", "Naive_B_cell", "Intermediate_B_cell"]
SUB = dict(CD4_Naive_T_cell="T", CD8_Naive_T_cell="T", CD4_Memory_T_cell="T", CD8_Memory_T_cell="T", Treg="T", MAIT="T", NK_CD56low="NK", NK_CD56hi="NK",
           Naive_B_cell="B", Intermediate_B_cell="B", CD14_Mono="M", CD16_Mono="M", cDC="M")
CELLS = MYE + LYM
raw = pd.read_csv(sys.argv[1], index_col=0)
d = raw[(raw.adj_p_value < 0.05) & (raw.log_fc.abs() > 0.5) & raw.celltype.isin(CELLS)]
cyts = sorted(d.cytokine.unique()); genes = sorted(d.gene.unique())
qi = {q: i for i, q in enumerate(cyts)}; gi = {g: i for i, g in enumerate(genes)}; cx = {c: i for i, c in enumerate(CELLS)}
R = np.zeros((len(cyts), len(CELLS), len(genes)), np.float32)
R[d.cytokine.map(qi).values, d.celltype.map(cx).values, d.gene.map(gi).values] = d.log_fc.values
J = np.einsum("qcg,qdg->qcd", R, R, optimize=True).astype(float)
dev = np.array([int(hashlib.md5(q.encode()).hexdigest(), 16) % 2 == 0 for q in cyts])
def dist(mask):
    S = J[mask].sum(0); n = np.sqrt(np.diag(S)); return 1 - S / np.outer(n, n)
iu = np.triu_indices(len(CELLS), 1); a, b = dist(dev)[iu], dist(~dev)[iu]
same = np.array([SUB[CELLS[i]] == SUB[CELLS[j]] for i, j in zip(*iu)])
print("cytokines dev/test", int(dev.sum()), int((~dev).sum()), "| cell-type pairs", len(a))
print("Spearman, behavioural distance (dev cytokines) vs response divergence (test cytokines): all pairs %.2f" % spearmanr(a, b)[0])
print("  lineage label alone (same vs different) vs test divergence: %.2f" % spearmanr(~same, b)[0])
print("  within same-lineage pairs (n=%d): %.2f ; within different-lineage pairs (n=%d): %.2f" % (same.sum(), spearmanr(a[same], b[same])[0], (~same).sum(), spearmanr(a[~same], b[~same])[0]))
# strongest-cytokine dependence: drop the 10 largest responders from both halves
tot = J.trace(axis1=1, axis2=2); big = np.argsort(tot)[-10:]; m = np.ones(len(cyts), bool); m[big] = False
a2, b2 = dist(dev & m)[iu], dist(~dev & m)[iu]
print("  without the 10 strongest cytokines: all %.2f ; within-lineage %.2f" % (spearmanr(a2, b2)[0], spearmanr(a2[same], b2[same])[0]))
# confound check: sparse thresholded vectors make distance depend on statistical power (cell number)
ncell = raw[raw.celltype.isin(CELLS)].groupby("celltype").median_num_cells_pbs.median().reindex(CELLS).values
pw = np.abs(np.log(ncell)[:, None] - np.log(ncell)[None])[iu]
nd = np.log((R != 0).sum((0, 2))); pw2 = np.abs(nd[:, None] - nd[None])[iu]
from scipy.stats import rankdata
def partial(x, y, zs):
    Z = np.column_stack([np.ones(len(x))] + [rankdata(z) for z in zs]); rx, ry = rankdata(x), rankdata(y)
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]; ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1])
print("power confound: cell-number difference vs test divergence %.2f ; changed-gene-count difference vs test divergence %.2f" % (spearmanr(pw, b)[0], spearmanr(pw2, b)[0]))
print("partial Spearman of behavioural distance with test divergence, controlling lineage, cell-number and gene-count differences: %.2f" % partial(a, b, [(~same).astype(float), pw, pw2]))
print("same control applied to lineage label alone: %.2f" % partial((~same).astype(float), b, [pw, pw2]))
