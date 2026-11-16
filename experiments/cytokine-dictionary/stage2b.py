"""Route-informed prior: does knowing which ligands are induced help the similarity model?

Linear model as in stage2.py, but the ridge penalty on cytokine s is divided by (1 + gamma * a_s),
where a_s is the scaled induction of s's ligand under q in the observed block. gamma = 0 is the plain linear model.
gamma is chosen on a development half of the cytokines (split by name hash) and scored once on the other half.
Criterion fixed before running: skill gain on the test half > 0 with a bootstrap 95% interval excluding zero.
"""
import sys, json, hashlib
import numpy as np, pandas as pd
from scipy.optimize import nnls
MYE = ["CD14_Mono", "CD16_Mono", "cDC"]
LYM = ["CD4_Naive_T_cell", "CD8_Naive_T_cell", "CD4_Memory_T_cell", "CD8_Memory_T_cell", "Treg", "MAIT", "NK_CD56low", "NK_CD56hi", "Naive_B_cell", "Intermediate_B_cell"]
CELLS = MYE + LYM; RIDGE = 0.1; GAMMAS = [0, 1, 3, 10, 30]
raw = pd.read_csv(sys.argv[1], index_col=0); ci = pd.read_excel(sys.argv[2], sheet_name="all_cytokines").set_index("name")
d = raw[(raw.adj_p_value < 0.05) & (raw.log_fc.abs() > 0.5) & raw.celltype.isin(CELLS)]
cyts = sorted(d.cytokine.unique()); genes = sorted(d.gene.unique())
qi = {q: i for i, q in enumerate(cyts)}; gi = {g: i for i, g in enumerate(genes)}; cidx = {c: i for i, c in enumerate(CELLS)}
Q, C = len(cyts), len(CELLS)
R = np.zeros((Q, C, len(genes)), np.float32)
R[d.cytokine.map(qi).values, d.celltype.map(cidx).values, d.gene.map(gi).values] = d.log_fc.values
K = np.einsum("qcg,pcg->cqp", R, R, optimize=True).astype(np.float64); nde = (R != 0).sum(-1)
frac = raw[raw.celltype.isin(CELLS)].groupby("celltype").median_num_cells_pbs.median().reindex(CELLS); frac = (frac / frac.sum()).values
expr = raw[raw.celltype.isin(CELLS)].groupby(["gene", "celltype"]).logCPM.median()
lig = {s: [g.strip() for g in str(ci.loc[s, "gene"]).split(",")] for s in cyts}
def induction(cells, shuffle=None):
    a = np.zeros((Q, Q))
    for s in cyts:
        for gname in lig[s]:
            supply = {c: frac[cidx[c]] * 2 ** expr[(gname, c)] for c in CELLS if (gname, c) in expr.index}; tot = sum(supply.values())
            if tot == 0: continue
            for r in d[(d.gene == gname) & d.celltype.isin(cells)].itertuples():
                a[qi[r.cytokine], qi[s]] += supply[r.celltype] / tot * r.log_fc / len(lig[s])
    a = np.maximum(a, 0); np.fill_diagonal(a, 0)
    if shuffle is not None:                                       # negative control: same induction values, wrong ligands
        for q in range(Q): a[q] = a[q, shuffle.permutation(Q)]
        np.fill_diagonal(a, 0)
    return a / np.maximum(a.max(1, keepdims=True), 1e-9)
def run(shuffle=None):
    rows = []
    for task, hidden in (("myeloid", MYE), ("lymphoid", LYM)):
        hid = [cidx[c] for c in hidden]; obs = [i for i in range(C) if i not in hid]
        a = induction([CELLS[i] for i in obs], shuffle); Kobs = K[obs].sum(0)
        for q in range(Q):
            keep = np.arange(Q) != q; Gm = Kobs[np.ix_(keep, keep)]; b = Kobs[keep, q]; base = RIDGE * (np.trace(Gm) / len(b) + 1e-9)
            for gam in GAMMAS:
                L = np.linalg.cholesky(Gm + np.diag(base / (1 + gam * a[q, keep])) + 1e-9 * np.eye(len(b)))
                w = np.zeros(Q); w[keep] = nnls(L.T, np.linalg.solve(L, b))[0]
                for c in hid:
                    if nde[q, c] >= 20:
                        rows.append((task, cyts[q], CELLS[c], gam, K[c, q, q] - 2 * w @ K[c, :, q] + w @ K[c] @ w, K[c, q, q]))
    return pd.DataFrame(rows, columns=["task", "cytokine", "cell", "gamma", "sse", "sse0"])
def skill(df): return df.groupby(["gamma", "cytokine"]).apply(lambda g: 1 - g.sse.sum() / g.sse0.sum(), include_groups=False).unstack("gamma")
dev = [q for q in cyts if int(hashlib.md5(q.encode()).hexdigest(), 16) % 2 == 0]
out = {}
for name, shuf in (("real ligand routes", None), ("shuffled ligand routes (control)", np.random.default_rng(1))):
    S = skill(run(shuf)); Sd, St = S[S.index.isin(dev)], S[~S.index.isin(dev)]
    best = Sd.mean().idxmax(); diff = St[best] - St[0]; rng = np.random.default_rng(0)
    bs = [diff.sample(len(diff), replace=True, random_state=int(rng.integers(1e9))).mean() for _ in range(4000)]
    out[name] = dict(dev_skill_by_gamma={int(g): round(float(v), 3) for g, v in Sd.mean().items()}, chosen_gamma=int(best),
                     test_n_cytokines=len(St), test_skill_plain=round(float(St[0].mean()), 3), test_skill_route_prior=round(float(St[best].mean()), 3),
                     test_gain=round(float(diff.mean()), 3), ci95=[round(float(np.percentile(bs, 2.5)), 3), round(float(np.percentile(bs, 97.5)), 3)],
                     test_skill_by_gamma={int(g): round(float(v), 3) for g, v in St.mean().items()},
                     biggest_gains=diff.sort_values().tail(5).round(2).to_dict(), biggest_losses=diff.sort_values().head(3).round(2).to_dict())
    print(json.dumps(out[name]), flush=True)
json.dump(out, open(sys.argv[3], "w"), indent=1)
