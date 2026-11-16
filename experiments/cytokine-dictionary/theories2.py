"""Round 2 of theory tests on the real table. Criteria fixed before the first run:

T6 shared latent activation, cell-type-specific readout (low-rank tensor R[q,c,g] = sum_r U[q,r] V[c,r,g]):
   rank chosen inside the development cytokines; PASS if it beats the similarity model on test cytokines with a
   95% interval above zero; COMPARABLE if the interval includes zero; FAIL if below.
T7 opposing regulatory axes: six textbook antagonistic pairs; PASS if at least 4 of 6 have a negative response
   cosine in CD14 monocytes that is below the 5th percentile of that cytokine's other pairs.
Robustness: the same comparison repeated on genes that are changed against at least 4 of 6 individual PBS wells.

Regime for T6: for each task the hidden lineage block is hidden for ALL test cytokines at once. Every model is
scored under that same regime, so numbers differ from the leave-one-out numbers in theories.py.
"""
import sys, json, hashlib
import numpy as np, pandas as pd, torch
from scipy.optimize import nnls
torch.set_num_threads(2)
MYE = ["CD14_Mono", "CD16_Mono", "cDC"]
LYM = ["CD4_Naive_T_cell", "CD8_Naive_T_cell", "CD4_Memory_T_cell", "CD8_Memory_T_cell", "Treg", "MAIT", "NK_CD56low", "NK_CD56hi", "Naive_B_cell", "Intermediate_B_cell"]
CELLS = MYE + LYM; UV = (0.86, 0.48)            # stacking weights fitted on development cytokines in theories.py
raw = pd.read_csv(sys.argv[1], index_col=0)
rng = np.random.default_rng(0)
def boot(x, n=4000):
    x = np.asarray(x); bs = [x[rng.integers(0, len(x), len(x))].mean() for _ in range(n)]
    return [round(float(np.percentile(bs, 2.5)), 3), round(float(np.percentile(bs, 97.5)), 3)]
def nn_solve(Gm, b, lam=0.1):
    Gm = Gm + lam * np.eye(len(b)) * (np.trace(Gm) / len(b) + 1e-9)
    L = np.linalg.cholesky(Gm + 1e-9 * np.eye(len(b))); return nnls(L.T, np.linalg.solve(L, b))[0]

def build(min_wells):
    d = raw[(raw.adj_p_value < 0.05) & (raw.log_fc.abs() > 0.5) & raw.celltype.isin(CELLS) & (raw.num_DE_pbs_wells >= min_wells)]
    cyts = sorted(raw[raw.celltype.isin(CELLS)].cytokine.unique()); genes = sorted(d.gene.unique())
    qi = {q: i for i, q in enumerate(cyts)}; gi = {g: i for i, g in enumerate(genes)}; cx = {c: i for i, c in enumerate(CELLS)}
    R = np.zeros((len(cyts), len(CELLS), len(genes)), np.float32)
    R[d.cytokine.map(qi).values, d.celltype.map(cx).values, d.gene.map(gi).values] = d.log_fc.values
    return cyts, R

def tensor_fit(Z, mask, rank, seed, steps=400):
    """Z [Q, C, D] rotated responses; mask [Q, C] observed entries."""
    g = torch.Generator().manual_seed(seed); Q, C, D = Z.shape
    U = (torch.randn(Q, rank, generator=g) * 0.3).requires_grad_(); V = (torch.randn(C, rank, D, generator=g) * 0.02).requires_grad_()
    Zt = torch.tensor(Z); m = torch.tensor(mask, dtype=torch.float32)[:, :, None]
    opt = torch.optim.Adam([U, V], lr=0.03); sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, steps)
    for _ in range(steps):
        opt.zero_grad(); pred = torch.matmul(U, V).transpose(0, 1)
        loss = (((pred - Zt) ** 2) * m).sum() / m.sum() + 1e-3 * (U.pow(2).mean() + V.pow(2).mean() * D)
        loss.backward(); opt.step(); sch.step()
    with torch.no_grad(): return torch.matmul(U, V).transpose(0, 1).numpy()

def evaluate(tag, min_wells, do_tensor=True):
    cyts, R = build(min_wells); Q, C, G = R.shape
    X = R.reshape(Q * C, G).astype(np.float64)
    B = np.linalg.svd(X, full_matrices=False)[2]                      # exact rotation onto the row space
    Z = (X @ B.T).reshape(Q, C, -1).astype(np.float32); nde = (R != 0).sum(-1)
    M = np.einsum("qcd,sed->qcse", Z.astype(np.float64), Z.astype(np.float64), optimize=True)
    dev = np.array([int(hashlib.md5(q.encode()).hexdigest(), 16) % 2 == 0 for q in cyts]); devi, testi = np.where(dev)[0], np.where(~dev)[0]
    rows = []; rank_log = {}
    for task, hidden in (("myeloid", MYE), ("lymphoid", LYM)):
        hid = [CELLS.index(c) for c in hidden]; vis = [i for i in range(C) if i not in hid]
        # ---- tensor: choose rank inside development cytokines (hide the block for a random half of them)
        if do_tensor:
            dev2 = rng.permutation(devi)[: len(devi) // 2]; mask = np.ones((Q, C), bool); mask[np.ix_(testi, hid)] = False
            mt = mask.copy(); mt[np.ix_(dev2, hid)] = False; sc = {}
            for rk in (4, 8, 16, 32):
                P = tensor_fit(Z, mt, rk, 0); num = den = 0.0
                for q in dev2:
                    for c in hid:
                        if nde[q, c] >= 20: num += ((P[q, c] - Z[q, c]) ** 2).sum(); den += (Z[q, c] ** 2).sum()
                sc[rk] = round(float(1 - num / den), 3)
            best = max(sc, key=sc.get); rank_log[task] = dict(dev_skill_by_rank=sc, chosen=best)
            P = (tensor_fit(Z, mask, best, 1) + tensor_fit(Z, mask, best, 2)) / 2
        # ---- similarity and cell-axis models under the same regime: only development cytokines can be used as dictionary
        Kobs = sum(M[:, c, :, c] for c in vis)
        Avv = sum(M[s, vis][:, s][:, vis] for s in devi)                                    # sum_s <R(s,v), R(s,v')>
        gam_by_cell = {c: nn_solve(Avv, sum(M[s, c, s, vis] for s in devi)) for c in hid}
        for q in testi:
            beta = np.zeros(Q); beta[devi] = nn_solve(Kobs[np.ix_(devi, devi)], Kobs[devi, q])
            for c in hid:
                if nde[q, c] < 20: continue
                gam = gam_by_cell[c]
                tt = M[q, c, q, c]; pc = Z[:, c].astype(np.float64).T @ beta; pl = Z[q, vis].astype(np.float64).T @ gam; y = Z[q, c].astype(np.float64)
                rec = dict(task=task, cytokine=cyts[q], cell=CELLS[c], genes=int(nde[q, c]), tt=tt,
                           sse_similarity=((pc - y) ** 2).sum(), sse_cell=((pl - y) ** 2).sum(), sse_twoway=((UV[0] * pc + UV[1] * pl - y) ** 2).sum())
                if do_tensor:
                    rec["sse_tensor"] = ((P[q, c] - y) ** 2).sum()
                    rec["sse_tensor_plus_similarity"] = (((P[q, c] + pc) / 2 - y) ** 2).sum()
                rows.append(rec)
    df = pd.DataFrame(rows); models = [c[4:] for c in df.columns if c.startswith("sse_")]
    per = df.groupby("cytokine").apply(lambda g: pd.Series({m: 1 - g["sse_" + m].sum() / g.tt.sum() for m in models}), include_groups=False)
    res = dict(tag=tag, n_genes=G, n_test_cytokines=len(per), n_pairs=len(df), rank=rank_log, skill={m: round(float(per[m].mean()), 3) for m in models},
               twoway_minus_similarity=dict(mean=round(float((per.twoway - per.similarity).mean()), 3), ci95=boot(per.twoway - per.similarity)))
    if do_tensor:
        res["tensor_minus_similarity"] = dict(mean=round(float((per.tensor - per.similarity).mean()), 3), ci95=boot(per.tensor - per.similarity))
        res["tensor_minus_twoway"] = dict(mean=round(float((per.tensor - per.twoway).mean()), 3), ci95=boot(per.tensor - per.twoway))
        res["tensor_plus_similarity_minus_twoway"] = dict(mean=round(float((per.tensor_plus_similarity - per.twoway).mean()), 3), ci95=boot(per.tensor_plus_similarity - per.twoway))
        strong = df[df.genes >= 200]
        ps = strong.groupby("cytokine").apply(lambda g: pd.Series({m: 1 - g["sse_" + m].sum() / g.tt.sum() for m in models}), include_groups=False)
        res["strong_pairs_skill"] = dict(n_cytokines=len(ps), **{m: round(float(ps[m].mean()), 3) for m in models})
    print(json.dumps(res), flush=True); df.to_csv(f"theories2_{tag}.csv", index=False)
    return res, (cyts, R)

out = {}
out["T6_all_genes"], (cyts, R) = evaluate("all", 0)
out["T6_robust_genes"], _ = evaluate("robust", 4)

# ---------------- T7 opposing axes in monocytes
qi = {q: i for i, q in enumerate(cyts)}; nde = (R != 0).sum(-1); res = []
for a, b in (("IL-10", "IL-1-beta"), ("IL-10", "IFN-gamma"), ("IL-10", "TNF-alpha"), ("IL-4", "IFN-gamma"), ("IL-13", "IFN-gamma"), ("IL-4", "IL-1-beta")):
    for cell in ("CD14_Mono", "CD16_Mono", "cDC"):
        c = CELLS.index(cell); Y = R[:, c, :].astype(np.float64); n = np.linalg.norm(Y, axis=1)
        cosv = lambda s: float(Y[qi[a]] @ Y[s] / (n[qi[a]] * n[s]))
        null = np.array([cosv(s) for s in range(len(cyts)) if s not in (qi[a], qi[b]) and nde[s, c] >= 20])
        ya, yb = Y[qi[a]], Y[qi[b]]; sh = (ya != 0) & (yb != 0)
        res.append(dict(pair=f"{a} vs {b}", cell=cell, cosine=round(cosv(qi[b]), 3), null_p05=round(float(np.percentile(null, 5)), 3), null_median=round(float(np.median(null)), 3),
                        shared_genes=int(sh.sum()), opposite_sign_share=round(float((np.sign(ya[sh]) != np.sign(yb[sh])).mean()), 2) if sh.sum() else None))
out["T7"] = res
t7 = pd.DataFrame(res); print(t7.to_string())
json.dump(out, open(sys.argv[2], "w"), indent=1)
