"""Round 100 C1 evaluation per PREREG-C1-round100.md. Usage: python3 -I c1_eval.py <work_dir> <out_json>"""
import sys, glob, json, itertools, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
W, OUT = sys.argv[1:3]
C = pd.read_csv(f'{W}/cells.csv.gz', index_col=0); S = pd.read_csv(f'{W}/sub_meta.csv.gz', index_col=0).reset_index(drop=True)  # full subsample; latents index into it
PB = pd.read_pickle(f'{W}/pseudobulk.pkl')
C['keep'] = ~(C.qc_fail | C.doublet)
smp = C.groupby('sample').agg(donor=('donor_id', 'first'), tp=('timepoint', 'first')).sort_index()
both = smp.groupby('donor').tp.nunique(); smp = smp[smp.donor.map(both) == 2]
smp['y'] = (smp.tp == 'M3').astype(int); y = smp.y.values; donors = smp.donor.values; samples = smp.index.values
C = C[C['sample'].isin(samples)]; K = C[C.keep]
# removal report
rem = C.groupby(['sample']).agg(n=('keep', 'size'), qc=('qc_fail', 'sum'), dbl=('doublet', 'sum')).join(smp.tp)
prot = {t: float(C[C.major == t].qc_fail.mean()) for t in C.major.unique()}
# ---- B0
def b0(Kc, Call):
    nonB = Kc[~Kc.major.isin(['B_plasma', 'platelet_other'])]
    f = Kc.groupby('sample')[['log_counts', 'log_genes', 'mito', 'ribo', 'doublet_score', 'stress']].median()
    f['ig'] = nonB.groupby('sample').ig.median(); f['hb'] = nonB.groupby('sample').hb.median()
    f['log_cells'] = np.log(Call.groupby('sample').size())
    return f.loc[samples].fillna(0)
B0 = b0(K, C)
B1 = K[K.major == 'monocyte'].groupby('sample').ifn6.mean().loc[samples].to_frame()
cnt = K.groupby(['sample', 'major']).size().unstack(fill_value=0).loc[samples] + 0.5
B2 = np.log(cnt).sub(np.log(cnt).mean(1), axis=0)
# B3 pseudobulk: types with >=20 cells in every sample
tc = K.groupby(['sample', 'major']).size().unstack(fill_value=0).loc[samples]
pb_types = [t for t in tc.columns if (tc[t] >= 20).all()]
def b3(train_mask):
    blocks = []
    for t in pb_types:
        P = PB.loc[[f'{s}|{t}' for s in samples]].values.astype(float)
        L = np.log1p(P / P.sum(1, keepdims=True) * 1e6)
        v = L[train_mask].var(0); hv = np.argsort(-v)[:2000]; blocks.append(L[:, hv])
    return np.hstack(blocks)
# latents
types_all = [t for t in S.major.unique() if (S.groupby(['sample', 'major']).size().unstack(fill_value=0).reindex(samples).fillna(0)[t] >= 10).all()]
lat = {}
for f in glob.glob(f'{W}/latent/*.npz'):
    z = np.load(f); d = f.split('/')[-1][:-4]; lat[d] = z
def sample_means(Z, idx, cell_mask=None):
    m = S.iloc[idx].copy(); m['row'] = np.arange(len(idx))
    if cell_mask is not None: m = m[cell_mask[idx]]
    out = []
    for s in samples:
        row = []
        for t in types_all:
            r = m[(m['sample'] == s) & (m.major == t)].row.values
            row.append(Z[r].mean(0) if len(r) else np.full(Z.shape[1], np.nan))
        out.append(np.concatenate(row))
    return np.array(out)
def fold_latent(d, key, cell_mask=None):
    z = lat[d]; Z = np.zeros((len(S), z[f'{key}_tr'].shape[1])); Z[z['tr_idx']] = z[f'{key}_tr']; Z[z['te_idx']] = z[f'{key}_te']
    return sample_means(Z, np.arange(len(S)), cell_mask)
def lodo(feat_fn, hi_dim):
    p = np.zeros(len(samples))
    for d in np.unique(donors):
        te = donors == d; tr = ~te
        Xf = feat_fn(d, tr)
        Xf = np.nan_to_num(Xf - np.nanmean(Xf[tr], 0))
        if hi_dim:
            pc = PCA(n_components=min(10, tr.sum() - 1), random_state=0).fit(Xf[tr]); Xf = pc.transform(Xf)
        sc = StandardScaler().fit(Xf[tr]); clf = LogisticRegression(C=1.0, max_iter=5000).fit(sc.transform(Xf[tr]), y[tr])
        p[te] = clf.predict_proba(sc.transform(Xf[te]))[:, 1]
    return p
def metrics(p):
    d = np.array([p[(donors == u) & (y == 1)][0] - p[(donors == u) & (y == 0)][0] for u in np.unique(donors)])
    obs = d.mean(); null = [np.mean(d * np.array(s)) for s in itertools.product([1, -1], repeat=len(d))]
    return dict(auc=round(roc_auc_score(y, p), 3), paired_acc=round(float((d > 0).mean()), 3), signflip_p=round(float(np.mean(np.array(null) >= obs - 1e-12)), 4))
def run(cell_mask=None, B0x=None):
    B0u = B0 if B0x is None else B0x
    P = {'B0': lodo(lambda d, tr: B0u.values, False), 'B1': lodo(lambda d, tr: B1.values, False),
         'B2': lodo(lambda d, tr: B2.values, False), 'B3': lodo(lambda d, tr: b3(tr), True),
         'B4': lodo(lambda d, tr: fold_latent(d, 'pca', cell_mask), True), 'D1': lodo(lambda d, tr: fold_latent(d, 'scvi', cell_mask), True)}
    P['D1+B0'] = lodo(lambda d, tr: np.hstack([PCA(10, random_state=0).fit(np.nan_to_num(fold_latent(d, 'scvi', cell_mask)[tr])).transform(np.nan_to_num(fold_latent(d, 'scvi', cell_mask))), B0u.values]), False)
    return P
rng = np.random.default_rng(0); ud = np.unique(donors)
def boot(pa, pb):
    out = []
    for _ in range(2000):
        bd = rng.choice(ud, len(ud)); ii = np.concatenate([np.flatnonzero(donors == u) for u in bd])
        out.append(roc_auc_score(y[ii], pa[ii]) - roc_auc_score(y[ii], pb[ii]))
    return [round(float(np.percentile(out, 2.5)), 3), round(float(np.percentile(out, 97.5)), 3)]
P = run(); R = {k: metrics(v) for k, v in P.items()}
best = max(['B1', 'B2', 'B3', 'B4'], key=lambda k: R[k]['auc'])
r1 = dict(best=best, delta=round(R['D1']['auc'] - R[best]['auc'], 3), ci=boot(P['D1'], P[best]))
r1['pass'] = bool(r1['delta'] >= 0.05 and r1['ci'][0] > 0)
r2 = dict(delta=round(R['D1']['auc'] - R['B0']['auc'], 3), d1b0_gt_b0=bool(R['D1+B0']['auc'] > R['B0']['auc']))
r2['pass'] = bool(r2['delta'] >= 0.05 and r2['d1b0_gt_b0'])
# rule 3: QC matching within donor x major type on subsampled cells
mask = np.zeros(len(S), bool); cols = ['log_counts', 'log_genes', 'mito']; sd = S[cols].std().values
for (d, t), g in S.groupby(['donor_id', 'major']):
    a = g[g.timepoint == 'M3']; b = g[g.timepoint == 'Base']
    if len(a) == 0 or len(b) == 0: continue
    used = set(); A_ = a[cols].values / sd; B_ = b[cols].values / sd
    for i in range(len(a)):
        dist = np.sqrt(((B_ - A_[i]) ** 2).sum(1)); order = np.argsort(dist)
        for j in order:
            if dist[j] > 0.5 * np.sqrt(len(cols)): break
            if j not in used: used.add(j); mask[a.index[i]] = True; mask[b.index[j]] = True; break
Sm = S[mask].copy(); Sm['keep'] = True
B0m = Sm.groupby('sample')[['log_counts', 'log_genes', 'mito', 'ribo', 'doublet_score', 'stress']].median().reindex(samples)
nonBm = Sm[~Sm.major.isin(['B_plasma', 'platelet_other'])]
B0m['ig'] = nonBm.groupby('sample').ig.median(); B0m['hb'] = nonBm.groupby('sample').hb.median(); B0m['log_cells'] = np.log(Sm.groupby('sample').size().reindex(samples))
B0m = B0m.fillna(0)
Pm = run(mask, B0m); Rm = {k: metrics(v) for k, v in Pm.items() if k in ('B0', 'B4', 'D1')}
r3 = dict(matched_cells=int(mask.sum()), B0_matched=Rm['B0']['auc'], B4_matched=Rm['B4']['auc'], D1_matched=Rm['D1']['auc'],
          b0_check=bool(Rm['B0']['auc'] <= 0.6), d1_gain_survives=bool(Rm['D1']['auc'] - max(Rm['B4']['auc'], R[best]['auc'] if best != 'B4' else 0) >= 0.05))
r3['pass'] = bool(r3['b0_check'] and r3['d1_gain_survives'])
res = dict(results=R, rule1=r1, rule2=r2, rule3=r3, rule4='not applicable (library == donor x timepoint)', H1=bool(r1['pass'] and r2['pass'] and r3['pass']),
           n_samples=len(samples), pb_types=pb_types, latent_types=types_all, folds=len(lat),
           scvi_epochs={d: int(z['epochs']) for d, z in lat.items()}, qc_fail_rate_by_type={k: round(v, 3) for k, v in prot.items()},
           removal=rem.reset_index().to_dict('records'))
json.dump(res, open(OUT, 'w'), indent=1, default=str)
print(json.dumps({k: res[k] for k in ['results', 'rule1', 'rule2', 'rule3', 'H1']}, indent=1))
