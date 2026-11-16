"""T6, second attempt (declared after the first failed): the first tensor fit left latent coordinates that the
visible cells do not read at their random starting values. Fix: near-zero initialisation and a weight-decay
strength tuned, with rank, inside the development cytokines. Same regime, pairs and criterion as theories2.py."""
import sys, json, hashlib
import numpy as np, pandas as pd, torch
sys.argv = [sys.argv[0], sys.argv[1], "/dev/null"]
exec(open("theories2.py").read().split("out = {}")[0])
def fit(Z, mask, rank, wd, seed, steps=500):
    g = torch.Generator().manual_seed(seed); Q, C, D = Z.shape
    U = (torch.randn(Q, rank, generator=g) * 0.01).requires_grad_(); V = (torch.randn(C, rank, D, generator=g) * 0.01).requires_grad_()
    Zt = torch.tensor(Z); m = torch.tensor(mask, dtype=torch.float32)[:, :, None]
    opt = torch.optim.Adam([U, V], lr=0.03); sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, steps)
    scale = float((Zt ** 2 * m).sum() / m.sum())
    for _ in range(steps):
        opt.zero_grad(); pred = torch.matmul(U, V).transpose(0, 1)
        loss = (((pred - Zt) ** 2) * m).sum() / m.sum() / scale + wd * (U.pow(2).sum(1).mean() + V.pow(2).sum((1, 2)).mean())
        loss.backward(); opt.step(); sch.step()
    with torch.no_grad(): return torch.matmul(U, V).transpose(0, 1).numpy()
cyts, R = build(0); Q, C, G = R.shape; X = R.reshape(Q * C, G).astype(np.float64)
B = np.linalg.svd(X, full_matrices=False)[2]; Z = (X @ B.T).reshape(Q, C, -1).astype(np.float32); nde = (R != 0).sum(-1)
dev = np.array([int(hashlib.md5(q.encode()).hexdigest(), 16) % 2 == 0 for q in cyts]); devi, testi = np.where(dev)[0], np.where(~dev)[0]
rng = np.random.default_rng(0); rows = []; log = {}
def sk(P, qs, hid):
    num = den = 0.0
    for q in qs:
        for c in hid:
            if nde[q, c] >= 20: num += ((P[q, c] - Z[q, c]) ** 2).sum(); den += (Z[q, c] ** 2).sum()
    return float(1 - num / den)
for task, hidden in (("myeloid", MYE), ("lymphoid", LYM)):
    hid = [CELLS.index(c) for c in hidden]
    dev2 = rng.permutation(devi)[: len(devi) // 2]; mask = np.ones((Q, C), bool); mask[np.ix_(testi, hid)] = False
    mt = mask.copy(); mt[np.ix_(dev2, hid)] = False
    grid = {(rk, wd): round(sk(fit(Z, mt, rk, wd, 0), dev2, hid), 3) for rk in (4, 8, 16, 32) for wd in (1e-4, 1e-3, 1e-2, 1e-1)}
    best = max(grid, key=grid.get); log[task] = dict(dev_grid={f"r{k[0]}_wd{k[1]:g}": v for k, v in grid.items()}, chosen=list(best))
    P = np.mean([fit(Z, mask, best[0], best[1], s) for s in (1, 2, 3)], 0)
    for q in testi:
        for c in hid:
            if nde[q, c] >= 20: rows.append(dict(task=task, cytokine=cyts[q], cell=CELLS[c], sse_tensor2=float(((P[q, c] - Z[q, c]) ** 2).sum())))
t = pd.DataFrame(rows).merge(pd.read_csv("theories2_all.csv"), on=["task", "cytokine", "cell"])
models = ["similarity", "cell", "twoway", "tensor2"]
per = t.groupby("cytokine").apply(lambda g: pd.Series({m: 1 - g["sse_" + m].sum() / g.tt.sum() for m in models}), include_groups=False)
t["sse_tensor2_plus_similarity"] = np.nan
res = dict(tuning=log, n_pairs=len(t), skill={m: round(float(per[m].mean()), 3) for m in models},
           tensor2_minus_similarity=dict(mean=round(float((per.tensor2 - per.similarity).mean()), 3), ci95=boot(per.tensor2 - per.similarity)),
           tensor2_minus_twoway=dict(mean=round(float((per.tensor2 - per.twoway).mean()), 3), ci95=boot(per.tensor2 - per.twoway)))
print(json.dumps(res)); json.dump(res, open("theories2b.json", "w"), indent=1); t.to_csv("theories2b.csv", index=False)
