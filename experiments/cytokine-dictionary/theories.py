"""Five new theories tested on the real significant-gene table. Criteria were fixed before the first run:

T1 shared programs    cross-lineage transfer ratio of top-10 gene programs: pass >= 0.5, fail < 0.25
T2 receptor sharing   shared-chain pairs more similar than same-family pairs without a shared chain, permutation p < 0.01
T3 canonical axes     5 hub cytokines chosen on dev capture >= 50% of test response energy: pass; < 25%: fail
T4 cell axis          cell-axis model adds to the cytokine-similarity model on test cytokines, 95% interval excluding 0
T5 tonic signalling   antagonist vs agonist response cosine < 0 and below the 5th percentile of that antagonist's other pairs
"""
import sys, json, hashlib
import numpy as np, pandas as pd
from scipy.optimize import nnls

MYE = ["CD14_Mono", "CD16_Mono", "cDC"]
LYM = ["CD4_Naive_T_cell", "CD8_Naive_T_cell", "CD4_Memory_T_cell", "CD8_Memory_T_cell", "Treg", "MAIT", "NK_CD56low", "NK_CD56hi", "Naive_B_cell", "Intermediate_B_cell"]
SUB = dict(CD4_Naive_T_cell="T", CD8_Naive_T_cell="T", CD4_Memory_T_cell="T", CD8_Memory_T_cell="T", Treg="T", MAIT="T", NK_CD56low="NK", NK_CD56hi="NK",
           Naive_B_cell="B", Intermediate_B_cell="B", CD14_Mono="M", CD16_Mono="M", cDC="M")
CELLS = MYE + LYM
raw = pd.read_csv(sys.argv[1], index_col=0); ci = pd.read_excel(sys.argv[2], sheet_name="all_cytokines").set_index("name")
d = raw[(raw.adj_p_value < 0.05) & (raw.log_fc.abs() > 0.5) & raw.celltype.isin(CELLS)]
cyts = sorted(d.cytokine.unique()); genes = sorted(d.gene.unique())
qi = {q: i for i, q in enumerate(cyts)}; gi = {g: i for i, g in enumerate(genes)}; cx = {c: i for i, c in enumerate(CELLS)}
Q, C, G = len(cyts), len(CELLS), len(genes)
R = np.zeros((Q, C, G), np.float32)
R[d.cytokine.map(qi).values, d.celltype.map(cx).values, d.gene.map(gi).values] = d.log_fc.values
X = R.reshape(Q * C, G).astype(np.float64); M = (X @ X.T).reshape(Q, C, Q, C)          # M[q,c,s,c'] = <R(q,c), R(s,c')>
nde = (R != 0).sum(-1)
dev = np.array([int(hashlib.md5(q.encode()).hexdigest(), 16) % 2 == 0 for q in cyts]); test = ~dev
rng = np.random.default_rng(0); out = {}
def boot(x, n=4000):
    x = np.asarray(x); bs = [x[rng.integers(0, len(x), len(x))].mean() for _ in range(n)]
    return [round(float(np.percentile(bs, 2.5)), 3), round(float(np.percentile(bs, 97.5)), 3)]

# ------------------------------------------------------------------ T1 shared gene programs
k = 10; V = {}
for c in range(C):
    V[c] = np.linalg.svd(R[dev][:, c, :].astype(np.float64), full_matrices=False)[2][:k]
ev = np.zeros((C, C))
for a in range(C):
    for b in range(C):
        Y = R[test][:, b, :].astype(np.float64); ev[a, b] = ((Y @ V[a].T) ** 2).sum() / (Y ** 2).sum()
ratio = ev / np.diag(ev)[None, :]
lin = np.array([SUB[c] for c in CELLS]); mye = np.array([c in MYE for c in CELLS])
off = ~np.eye(C, dtype=bool)
same_sub = (lin[:, None] == lin[None]) & off; cross_block = mye[:, None] != mye[None]; cross_sub_lym = (lin[:, None] != lin[None]) & ~mye[:, None] & ~mye[None]
out["T1"] = dict(own_programs_explain_test=round(float(np.median(np.diag(ev))), 3),
                 transfer_ratio=dict(within_subtype=round(float(np.median(ratio[same_sub])), 3),
                                     across_lymphoid_subtypes=round(float(np.median(ratio[cross_sub_lym])), 3),
                                     myeloid_vs_lymphoid=round(float(np.median(ratio[cross_block])), 3),
                                     myeloid_vs_lymphoid_range=[round(float(ratio[cross_block].min()), 3), round(float(ratio[cross_block].max()), 3)]),
                 top1_program_share_of_own=round(float(np.median([((R[test][:, c, :].astype(np.float64) @ V[c][:1].T) ** 2).sum() / (R[test][:, c, :].astype(np.float64) ** 2).sum() for c in range(C)])), 3))
print("T1", json.dumps(out["T1"]), flush=True)

# ------------------------------------------------------------------ T2 receptor-chain sharing vs family vs cascade link
def chains(q):
    s = ci.loc[q, "receptor gene"]
    return set() if pd.isna(s) else {x.strip() for x in str(s).replace(";", ",").split(",") if x.strip()}
rec = {q: chains(q) for q in cyts}; fam = {q: ci.loc[q, "family_level_2"] for q in cyts}
strong = [i for i in range(Q) if (nde[i] >= 20).sum() >= 3]
Sq = np.einsum("qcsc->qs", M); nrm = np.sqrt(np.diag(Sq)); S = Sq / np.outer(nrm, nrm)
lig = {s: [g.strip() for g in str(ci.loc[s, "gene"]).split(",")] for s in cyts}
ind = np.zeros((Q, Q), bool)
for s in cyts:
    hit = d[d.gene.isin(lig[s]) & (d.log_fc > 1)]
    for q in hit.cytokine.unique(): ind[qi[q], qi[s]] = True
np.fill_diagonal(ind, False)
pairs = [(i, j) for a, i in enumerate(strong) for j in strong[a + 1:]]
sim = np.array([S[i, j] for i, j in pairs])
share = np.array([len(rec[cyts[i]] & rec[cyts[j]]) > 0 for i, j in pairs])
samef = np.array([fam[cyts[i]] == fam[cyts[j]] and not pd.isna(fam[cyts[i]]) for i, j in pairs])
casc = np.array([ind[i, j] or ind[j, i] for i, j in pairs])
obs = sim[share].mean() - sim[samef & ~share].mean()
perm = []
for _ in range(2000):
    p = rng.permutation(len(strong)); mp = {strong[a]: strong[p[a]] for a in range(len(strong))}
    sh = np.array([len(rec[cyts[mp[i]]] & rec[cyts[mp[j]]]) > 0 for i, j in pairs])
    sf = np.array([fam[cyts[mp[i]]] == fam[cyts[mp[j]]] and not pd.isna(fam[cyts[mp[i]]]) for i, j in pairs])
    if sh.sum() and (sf & ~sh).sum(): perm.append(sim[sh].mean() - sim[sf & ~sh].mean())
top = sim >= np.quantile(sim, 0.95)
def auc(lab):
    from scipy.stats import rankdata
    r = rankdata(sim); n1 = lab.sum(); n0 = (~lab).sum(); return round(float((r[lab].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)), 3)
out["T2"] = dict(n_cytokines=len(strong), n_pairs=len(pairs),
                 mean_similarity=dict(shared_chain=round(float(sim[share].mean()), 3), same_family_no_shared_chain=round(float(sim[samef & ~share].mean()), 3),
                                      cascade_linked_no_shared_chain=round(float(sim[casc & ~share].mean()), 3), unrelated=round(float(sim[~share & ~samef & ~casc].mean()), 3)),
                 n=dict(shared_chain=int(share.sum()), same_family_no_shared_chain=int((samef & ~share).sum()), cascade_linked_no_shared_chain=int((casc & ~share).sum())),
                 shared_minus_family=round(float(obs), 3), permutation_p=round(float((np.array(perm) >= obs).mean()), 4),
                 auc_predicting_similarity=dict(shared_chain=auc(share), same_family=auc(samef), cascade_link=auc(casc)),
                 top5pct_most_similar=dict(n=int(top.sum()), share_chain=round(float(share[top].mean()), 2), same_family=round(float(samef[top].mean()), 2),
                                           cascade_linked=round(float(casc[top].mean()), 2), none_of_these=round(float((~share & ~samef & ~casc)[top].mean()), 2)),
                 base_rates=dict(share_chain=round(float(share.mean()), 3), same_family=round(float(samef.mean()), 3), cascade_linked=round(float(casc.mean()), 3)))
print("T2", json.dumps(out["T2"]), flush=True)

# ------------------------------------------------------------------ T3 canonical axes (hub cytokines)
Mc = np.stack([M[:, c, :, c] for c in range(C)])                                           # [C, Q, Q]
def explained(basis, targets):
    num = 0.0; den = 0.0
    for c in range(C):
        Gb = Mc[c][np.ix_(basis, basis)] + 1e-8 * np.eye(len(basis)); B = Mc[c][np.ix_(basis, targets)]
        num += (B * np.linalg.solve(Gb, B)).sum(); den += np.diag(Mc[c])[targets].sum()
    return num / den
devi, testi = np.where(dev)[0], np.where(test)[0]; basis = []; path = []
for step in range(8):
    best = max((s for s in range(Q) if s not in basis), key=lambda s: explained(basis + [s], [t for t in devi if t not in basis + [s]]))
    basis.append(best)
    path.append(dict(k=step + 1, hub=cyts[best], dev=round(float(explained(basis, [t for t in devi if t not in basis])), 3),
                     test=round(float(explained(basis, [t for t in testi if t not in basis])), 3)))
b5 = basis[:5]; t5 = [t for t in testi if t not in b5]
rnd = [explained(list(rb), [t for t in testi if t not in rb]) for rb in (rng.choice(Q, 5, replace=False) for _ in range(300))]
svd_num = svd_den = 0.0
for c in range(C):
    Vt = np.linalg.svd(R[dev][:, c, :].astype(np.float64), full_matrices=False)[2][:5]; Y = R[t5][:, c, :].astype(np.float64)
    svd_num += ((Y @ Vt.T) ** 2).sum(); svd_den += (Y ** 2).sum()
out["T3"] = dict(path=path, test_energy_5_hubs=path[4]["test"], random_5_mean=round(float(np.mean(rnd)), 3), random_5_p95=round(float(np.percentile(rnd, 95)), 3),
                 best_5_dim_subspace_from_dev=round(float(svd_num / svd_den), 3))
print("T3", json.dumps(out["T3"]), flush=True)

# ------------------------------------------------------------------ T4 cell axis vs cytokine axis, hidden lineage block
def nn_solve(Gm, b, lam=0.1):
    Gm = Gm + lam * np.eye(len(b)) * (np.trace(Gm) / len(b) + 1e-9)
    L = np.linalg.cholesky(Gm + 1e-9 * np.eye(len(b))); return nnls(L.T, np.linalg.solve(L, b))[0]
rows = []
for task, hidden in (("myeloid", MYE), ("lymphoid", LYM)):
    hid = [cx[c] for c in hidden]; vis = [i for i in range(C) if i not in hid]
    Kobs = sum(Mc[c] for c in vis)
    Avv = np.einsum("qaqb->qab", M)[:, vis][:, :, vis]                                   # [Q, V, V]
    for q in range(Q):
        keep = np.arange(Q) != q
        beta = np.zeros(Q); beta[keep] = nn_solve(Kobs[np.ix_(keep, keep)], Kobs[keep, q])
        for c in hid:
            if nde[q, c] < 20: continue
            bvc = np.array([M[qq, c, qq, vis] for qq in range(Q)])                         # [Q, V] = <R(q',c), R(q',v)>
            gam = nn_solve(Avv[keep].sum(0), bvc[keep].sum(0))
            tt = M[q, c, q, c]
            pc_t = beta @ M[:, c, q, c]; pc_pc = beta @ Mc[c] @ beta                       # cytokine-axis prediction
            pl_t = gam @ M[q, vis, q, c]; pl_pl = gam @ Avv[q] @ gam                       # cell-axis prediction
            pc_pl = beta @ (M[:, c, q, :][:, vis] @ gam)
            rows.append(dict(task=task, q=q, cytokine=cyts[q], cell=CELLS[c], tt=tt, pc_t=pc_t, pc_pc=pc_pc, pl_t=pl_t, pl_pl=pl_pl, pc_pl=pc_pl))
df = pd.DataFrame(rows); isdev = df.q.map(lambda x: dev[x])
D = df[isdev]
uv = nn_solve(np.array([[D.pc_pc.sum(), D.pc_pl.sum()], [D.pc_pl.sum(), D.pl_pl.sum()]]), np.array([D.pc_t.sum(), D.pl_t.sum()]), 1e-9)
T = df[~isdev].copy()
T["sse_cyt"] = T.tt - 2 * T.pc_t + T.pc_pc; T["sse_cell"] = T.tt - 2 * T.pl_t + T.pl_pl
T["sse_both"] = T.tt - 2 * (uv[0] * T.pc_t + uv[1] * T.pl_t) + uv[0] ** 2 * T.pc_pc + uv[1] ** 2 * T.pl_pl + 2 * uv[0] * uv[1] * T.pc_pl
per = T.groupby("cytokine").apply(lambda g: pd.Series({m: 1 - g["sse_" + m].sum() / g.tt.sum() for m in ("cyt", "cell", "both")}), include_groups=False)
bytask = {t: {m: round(float(v), 3) for m, v in T[T.task == t].groupby("cytokine").apply(
    lambda g: pd.Series({m: 1 - g["sse_" + m].sum() / g.tt.sum() for m in ("cyt", "cell", "both")}), include_groups=False).mean().items()} for t in ("myeloid", "lymphoid")}
out["T4"] = dict(n_test_cytokines=len(per), n_test_pairs=len(T), stack_weights=dict(cytokine_axis=round(float(uv[0]), 2), cell_axis=round(float(uv[1]), 2)),
                 test_skill={m: round(float(per[m].mean()), 3) for m in ("cyt", "cell", "both")},
                 both_minus_cyt=dict(mean=round(float((per.both - per.cyt).mean()), 3), ci95=boot(per.both - per.cyt)),
                 cell_axis_ci95=boot(per.cell), by_hidden_block=bytask,
                 biggest_gains=(per.both - per.cyt).sort_values().tail(6).round(2).to_dict(), biggest_losses=(per.both - per.cyt).sort_values().head(4).round(2).to_dict())
print("T4", json.dumps(out["T4"]), flush=True)

# ------------------------------------------------------------------ T5 tonic signalling: antagonist vs agonist
res = []
for ant, ago in (("IL-1Ra", "IL-1-beta"), ("IL-1Ra", "IL-1-alpha"), ("IL-36Ra", "IL-36-alpha")):
    if ant not in qi or ago not in qi: continue
    a, g = qi[ant], qi[ago]
    for c in range(C):
        if nde[a, c] < 20 or nde[g, c] < 20: continue
        cosv = lambda s: Mc[c][a, s] / np.sqrt(Mc[c][a, a] * Mc[c][s, s])
        null = np.array([cosv(s) for s in range(Q) if s not in (a, g) and nde[s, c] >= 20])
        res.append(dict(antagonist=ant, agonist=ago, cell=CELLS[c], genes_antagonist=int(nde[a, c]), genes_agonist=int(nde[g, c]), cosine=round(float(cosv(g)), 3),
                        null_p05=round(float(np.percentile(null, 5)), 3), null_median=round(float(np.median(null)), 3), percentile_in_null=round(float((null < cosv(g)).mean()), 3)))
out["T5"] = res
print("T5", json.dumps(res), flush=True)
json.dump(out, open(sys.argv[3], "w"), indent=1)
