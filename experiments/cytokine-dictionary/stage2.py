"""Real-data predictive test: can induced-cytokine routes predict a held-out cell lineage's response?

For each cytokine q, one lineage block (myeloid or lymphoid) is hidden. Models see q's response
in the other block and every other cytokine's full response dictionary.

  zero      predict no change
  ctmean    average response of that cell type to all other cytokines
  linear    response-similarity model: non-negative ridge fit of q's observed response on other
            cytokines' observed responses; same coefficients applied in the hidden block
  mediator  ligand-route model: coefficient for cytokine s = potency_s x induction of s's ligand
            gene under q in the OBSERVED block; potencies fitted on the other cytokines only
  combined  u*linear + v*mediator

Criteria fixed before the first run:
  ADDS VALUE   combined beats linear on mean skill over responding hidden pairs,
               bootstrap 95% interval over cytokines excluding zero
  SUFFICIENT   mediator within 0.05 of linear on the same measure
Input is the pooled significant-gene table (adj p < 0.05, |log2FC| > 0.5); non-significant genes count as zero.
"""
import sys, json
import numpy as np, pandas as pd
from scipy.optimize import nnls

MYE = ["CD14_Mono", "CD16_Mono", "cDC"]
LYM = ["CD4_Naive_T_cell", "CD8_Naive_T_cell", "CD4_Memory_T_cell", "CD8_Memory_T_cell", "Treg", "MAIT",
       "NK_CD56low", "NK_CD56hi", "Naive_B_cell", "Intermediate_B_cell"]
CELLS = MYE + LYM
MIN_GENES = 20
RIDGE = float(sys.argv[4]) if len(sys.argv) > 4 else 0.1

raw = pd.read_csv(sys.argv[1], index_col=0)
ci = pd.read_excel(sys.argv[2], sheet_name="all_cytokines").set_index("name")
d = raw[(raw.adj_p_value < 0.05) & (raw.log_fc.abs() > 0.5) & raw.celltype.isin(CELLS)]
cyts = sorted(d.cytokine.unique()); genes = sorted(d.gene.unique())
qi = {q: i for i, q in enumerate(cyts)}; gi = {g: i for i, g in enumerate(genes)}; cidx = {c: i for i, c in enumerate(CELLS)}
Q, C, G = len(cyts), len(CELLS), len(genes)
R = np.zeros((Q, C, G), np.float32)
R[d.cytokine.map(qi).values, d.celltype.map(cidx).values, d.gene.map(gi).values] = d.log_fc.values
K = np.einsum("qcg,pcg->cqp", R, R, optimize=True).astype(np.float64)       # K[c, q, p] = <R(q,c), R(p,c)>
nde = (R != 0).sum(-1)

# ---- ligand induction a[q, s, block]: expression-share-weighted log2FC of s's ligand gene(s) under q
frac = raw[raw.celltype.isin(CELLS)].groupby("celltype").median_num_cells_pbs.median().reindex(CELLS); frac = (frac / frac.sum()).values
expr = raw[raw.celltype.isin(CELLS)].groupby(["gene", "celltype"]).logCPM.median()
lig = {s: [g.strip() for g in str(ci.loc[s, "gene"]).split(",")] for s in cyts}
def induction(block):
    cells = [c for c in CELLS if c in block]
    a = np.zeros((Q, Q))
    for s in cyts:
        for gname in lig[s]:
            supply = {c: frac[cidx[c]] * 2 ** expr[(gname, c)] for c in CELLS if (gname, c) in expr.index}
            tot = sum(supply.values())
            if tot == 0: continue
            rows = d[(d.gene == gname) & d.celltype.isin(cells)]
            for r in rows.itertuples():
                a[qi[r.cytokine], qi[s]] += supply[r.celltype] / tot * r.log_fc / len(lig[s])
    a = np.maximum(a, 0); np.fill_diagonal(a, 0)
    return a

def nn_solve(Gm, b, lam):
    Gm = Gm + lam * np.eye(len(b)) * (np.trace(Gm) / len(b) + 1e-9)
    L = np.linalg.cholesky(Gm + 1e-9 * np.eye(len(b)))
    return nnls(L.T, np.linalg.solve(L, b))[0]

def sse(w, c, q):                                                           # ||R(q,c) - sum_s w_s R(s,c)||^2 via kernels
    return K[c, q, q] - 2 * w @ K[c, :, q] + w @ K[c] @ w
def cosw(w, c, q):
    n = w @ K[c] @ w
    return float(w @ K[c, :, q] / np.sqrt(n * K[c, q, q])) if n > 0 and K[c, q, q] > 0 else 0.0

rows = []
for task, hidden in (("predict myeloid from lymphoid", MYE), ("predict lymphoid from myeloid", LYM)):
    hid = [cidx[c] for c in hidden]; obs = [i for i in range(C) if i not in hid]
    a = induction([CELLS[i] for i in obs])
    Kobs, Khid = K[obs].sum(0), K[hid]
    for q in range(Q):
        keep = np.arange(Q) != q
        # linear: fit in the observed block
        beta = np.zeros(Q); beta[keep] = nn_solve(Kobs[np.ix_(keep, keep)], Kobs[keep, q], RIDGE)
        # mediator potencies: fitted on the other cytokines' hidden-block responses (never q's)
        A = a.copy(); A[q, :] = 0; A[:, q] = 0                               # q is neither a training target nor a usable feature
        Gm = np.zeros((Q, Q)); b = np.zeros(Q)
        for c in hid:
            Gm += (A.T @ A) * K[c]; b += (A * K[c].T).sum(0)                 # b_s = sum_q' a[q',s] K[c][s,q']
        kappa = np.zeros(Q); kappa[keep] = nn_solve(Gm[np.ix_(keep, keep)], b[keep], RIDGE)
        wmed = kappa * a[q]; wmed[q] = 0
        wmean = keep / keep.sum()
        # combined weights (u, v): fitted on training cytokines, using their own linear and mediator predictions
        num = np.zeros((2, 2)); rhs = np.zeros(2)
        for qp in np.where(keep)[0][::3]:                                    # every third training cytokine, for speed
            k2 = keep & (np.arange(Q) != qp)
            bl = np.zeros(Q); bl[k2] = nn_solve(Kobs[np.ix_(k2, k2)], Kobs[k2, qp], RIDGE)
            wm = kappa * a[qp]; wm[qp] = 0; wm[q] = 0
            for c in hid:
                P = np.array([[bl @ K[c] @ bl, bl @ K[c] @ wm], [bl @ K[c] @ wm, wm @ K[c] @ wm]])
                num += P; rhs += np.array([bl @ K[c, :, qp], wm @ K[c, :, qp]])
        uv = nn_solve(num, rhs, 1e-6)
        wcomb = uv[0] * beta + uv[1] * wmed
        for c in hid:
            if K[c, q, q] == 0 and nde[q, c] == 0: z = 1.0
            z = K[c, q, q]
            rec = dict(task=task, cytokine=cyts[q], cell=CELLS[c], genes=int(nde[q, c]), u=uv[0], v=uv[1],
                       n_routes=int((wmed > 0).sum()))
            for name, w in (("ctmean", wmean), ("linear", beta), ("mediator", wmed), ("combined", wcomb)):
                rec["sse_" + name] = sse(w, c, q); rec["cos_" + name] = cosw(w, c, q)
            rec["sse_zero"] = z
            rows.append(rec)
df = pd.DataFrame(rows); df.to_csv(sys.argv[3], index=False)

def summarise(sub, label):
    models = ["ctmean", "linear", "mediator", "combined"]
    percy = sub.groupby("cytokine").apply(lambda g: pd.Series(
        {m: 1 - g["sse_" + m].sum() / g.sse_zero.sum() for m in models} | {"cos_" + m: g["cos_" + m].mean() for m in models}))
    rng = np.random.default_rng(0)
    def ci(x):
        bs = [x.sample(len(x), replace=True, random_state=int(rng.integers(1e9))).mean() for _ in range(2000)]
        return [round(float(np.percentile(bs, 2.5)), 3), round(float(np.percentile(bs, 97.5)), 3)]
    out = dict(label=label, n_pairs=len(sub), n_cytokines=len(percy),
               skill={m: round(float(percy[m].mean()), 3) for m in models},
               cosine={m: round(float(percy["cos_" + m].mean()), 3) for m in models},
               combined_minus_linear=dict(mean=round(float((percy.combined - percy.linear).mean()), 3), ci95=ci(percy.combined - percy.linear)),
               mediator_minus_linear=dict(mean=round(float((percy.mediator - percy.linear).mean()), 3), ci95=ci(percy.mediator - percy.linear)),
               mediator_minus_ctmean=dict(mean=round(float((percy.mediator - percy.ctmean).mean()), 3), ci95=ci(percy.mediator - percy.ctmean)),
               cos_mediator_minus_linear=dict(mean=round(float((percy.cos_mediator - percy.cos_linear).mean()), 3), ci95=ci(percy.cos_mediator - percy.cos_linear)))
    print(json.dumps(out)); return out

resp = df[df.genes >= MIN_GENES]
res = [summarise(resp, "all responding hidden pairs")]
for t in df.task.unique(): res.append(summarise(resp[resp.task == t], t))
pd.set_option("display.width", 250)
show = resp[resp.cytokine.isin(["IL-12", "IL-2", "IL-15", "IL-7", "IL-21", "GM-CSF", "TSLP", "IFN-gamma", "IL-1-beta"]) & resp.cell.isin(["CD14_Mono", "CD4_Naive_T_cell", "NK_CD56low"])]
print(show.assign(**{m: 1 - show["sse_" + m] / show.sse_zero for m in ["ctmean", "linear", "mediator", "combined"]})[
    ["task", "cytokine", "cell", "genes", "n_routes", "ctmean", "linear", "mediator", "combined", "cos_linear", "cos_mediator"]].round(2).to_string())
json.dump(res, open(sys.argv[3].replace(".csv", ".json"), "w"), indent=1)
