"""Round 102 (Q36) evaluation of C3/C4 per PREREG-C3C4-round102.md. Usage: python3 -I p_eval.py <work_dir> <C3|C4> <out_json>"""
import sys, glob, json, itertools, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.neighbors import KNeighborsClassifier
W, C, OUT = sys.argv[1:4]
Call = pd.read_pickle(f'{W}/cells.pkl'); Call = Call[Call.contrast == C]
S = pd.read_pickle(f'{W}/sub_meta.pkl'); S = S[S.contrast == C].reset_index()
PB = pd.read_pickle(f'{W}/pseudobulk.pkl')
smp = Call.groupby('sample').agg(donor=('donor_id', 'first'), y=('y', 'first')).sort_index()
samples = smp.index.values; y = smp.y.values; donors = smp.donor.astype(str).values
fo = pd.read_csv(f'{W}/latent_{C}/folds.csv', index_col=0).iloc[:, 0]
fold0 = (pd.Series(donors).map(fo) if C == 'C3' else pd.Series(samples).map(fo)).values
def b0(K, Ca):
    nonB = K[~K.author_cell_type.isin(['B', 'PB'])]
    f = K.groupby('sample')[['log_counts', 'log_genes', 'mito', 'ribo', 'stress']].median()
    f['ig'] = nonB.groupby('sample').ig.median(); f['hb'] = nonB.groupby('sample').hb.median()
    f['log_cells'] = np.log(Ca.groupby('sample').size()); return f.reindex(samples).fillna(0)
def basic(keepcol):
    K = Call[Call[keepcol]]
    B1 = K[K.author_cell_type.isin(['cM', 'ncM'])].groupby('sample').ifn6.mean().reindex(samples).to_frame()
    cnt = K.groupby(['sample', 'author_cell_type'], observed=True).size().unstack(fill_value=0).reindex(samples).fillna(0) + 0.5
    B2 = np.log(cnt).sub(np.log(cnt).mean(1), axis=0); return b0(K, Call), B1, B2, K
B0, B1, B2, K = basic('keep')
tc = K.groupby(['sample', 'author_cell_type'], observed=True).size().unstack(fill_value=0).reindex(samples).fillna(0)
pb_types = [t for t in tc.columns if (tc[t] >= 20).all()]
PBL = {t: (lambda P: np.log1p(P / P.sum(1, keepdims=True) * 1e6))(PB.loc[[f'{C}#{s}#{t}' for s in samples]].values.astype(float)) for t in pb_types}
def b3(tr): return np.hstack([L[:, np.argsort(-L[tr].var(0))[:2000]] for L in PBL.values()])
sc_ = S.groupby(['sample', 'author_cell_type'], observed=True).size().unstack(fill_value=0).reindex(samples).fillna(0)
types_all = [t for t in sc_.columns if (sc_[t] >= 10).all()]
lat = {f.split('/')[-1][:-4]: np.load(f) for f in glob.glob(f'{W}/latent_{C}/*.npz')}
sidx = {s: {t: np.flatnonzero((S['sample'] == s).values & (S.author_cell_type == t).values) for t in types_all} for s in samples}
def fold_latent(d, key, cell_mask=None):
    z = lat[d]; Z = np.zeros((len(S), z[f'{key}_tr'].shape[1])); Z[z['tr_idx']] = z[f'{key}_tr']; Z[z['te_idx']] = z[f'{key}_te']
    out = []
    for s in samples:
        row = []
        for t in types_all:
            r = sidx[s][t]; r = r[cell_mask[r]] if cell_mask is not None else r
            row.append(Z[r].mean(0) if len(r) else np.full(Z.shape[1], np.nan))
        out.append(np.concatenate(row))
    return np.array(out)
cache = {}
def feats(name, d, tr, cell_mask=None, B0u=None):
    if name == 'B0': return (B0 if B0u is None else B0u).values
    if name == 'B1': return B1.values
    if name == 'B2': return B2.values
    if name == 'B3': return b3(tr)
    k = (name, d, None if cell_mask is None else 'm')
    if k not in cache: cache[k] = fold_latent(d, {'B4': 'pca', 'D1': 'scvi'}[name], cell_mask)
    return cache[k]
HI = {'B3', 'B4', 'D1'}
def prep(Xf, tr, hi):
    Xf = np.nan_to_num(Xf - np.nanmean(Xf[tr], 0))
    if hi: Xf = PCA(n_components=min(10, tr.sum() - 1), random_state=0).fit(Xf[tr]).transform(Xf)
    return Xf
def cv(name, folds, yy=None, cell_mask=None, B0u=None, d1b0=False):
    yy = y if yy is None else yy; p = np.zeros(len(samples))
    for d in np.unique(folds):
        te = folds == d; tr = ~te
        if d1b0: Xf = np.hstack([prep(feats('D1', d, tr, cell_mask), tr, True), (B0 if B0u is None else B0u).values])
        else: Xf = prep(feats(name, d, tr, cell_mask, B0u), tr, name in HI)
        s = StandardScaler().fit(Xf[tr]); p[te] = LogisticRegression(C=1.0, max_iter=5000).fit(s.transform(Xf[tr]), yy[tr]).predict_proba(s.transform(Xf[te]))[:, 1]
    return p
def paired(p):
    d = np.array([p[(donors == u) & (y == 1)][0] - p[(donors == u) & (y == 0)][0] for u in np.unique(donors)])
    null = [np.mean(d * np.array(sg)) for sg in itertools.product([1, -1], repeat=len(d))]
    return dict(paired_acc=round(float((d > 0).mean()), 3), signflip_p=round(float(np.mean(np.array(null) >= d.mean() - 1e-12)), 4))
NAMES = ['B0', 'B1', 'B2', 'B3', 'B4', 'D1']
P = {n: cv(n, fold0) for n in NAMES}; P['D1+B0'] = cv(None, fold0, d1b0=True)
R = {n: dict(auc=round(roc_auc_score(y, p), 3)) for n, p in P.items()}
if C == 'C3':
    for n in R: R[n].update(paired(P[n]))
else:
    # 20 repeats for B0-B3; permutation null (classifier refit only, latents fixed and label-blind) on seed-0 split
    for n in ['B0', 'B1', 'B2', 'B3']:
        a = []
        for sd in range(20):
            f = np.empty(len(samples), object)
            for i, (_, te) in enumerate(StratifiedKFold(5, shuffle=True, random_state=sd).split(samples, y)): f[te] = f'f{i}'
            a.append(roc_auc_score(y, cv(n, f)))
        R[n]['auc_20rep_mean'] = round(float(np.mean(a)), 3); R[n]['auc_20rep_sd'] = round(float(np.std(a)), 3)
    rng = np.random.default_rng(1)
    for n in NAMES:
        obs = R[n]['auc']; null = [roc_auc_score(yp, cv(n, fold0, yy=yp)) for yp in (rng.permutation(y) for _ in range(1000 if n not in HI else 200))]
        R[n]['perm_p'] = round(float((1 + np.sum(np.array(null) >= obs)) / (1 + len(null))), 4); R[n]['perm_n'] = len(null)
rng = np.random.default_rng(0); ud = np.unique(donors)
def boot(pa, pb):
    out = []
    for _ in range(2000):
        bd = rng.choice(ud, len(ud)); ii = np.concatenate([np.flatnonzero(donors == u) for u in bd])
        if len(np.unique(y[ii])) < 2: continue
        out.append(roc_auc_score(y[ii], pa[ii]) - roc_auc_score(y[ii], pb[ii]))
    return [round(float(np.percentile(out, 2.5)), 3), round(float(np.percentile(out, 97.5)), 3)]
best = max(['B1', 'B2', 'B3', 'B4'], key=lambda k: R[k]['auc'])
r1 = dict(best=best, delta=round(R['D1']['auc'] - R[best]['auc'], 3), ci=boot(P['D1'], P[best])); r1['pass'] = bool(r1['delta'] >= 0.05 and r1['ci'][0] > 0)
r2 = dict(delta=round(R['D1']['auc'] - R['B0']['auc'], 3), d1b0=R['D1+B0']['auc'], d1b0_gt_b0=bool(R['D1+B0']['auc'] > R['B0']['auc'])); r2['pass'] = bool(r2['delta'] >= 0.05 and r2['d1b0_gt_b0'])
# rule 3: QC matching
mask = np.zeros(len(S), bool); cols = ['log_counts', 'log_genes', 'mito']; sdv = S[cols].std().values
grp = ['donor_id', 'author_cell_type'] if C == 'C3' else ['library_uuid', 'author_cell_type']
for _, g in S.groupby(grp, observed=True):
    a = g[g.y == 1]; b = g[g.y == 0]
    if len(a) == 0 or len(b) == 0: continue
    used = set(); A_ = a[cols].values / sdv; B_ = b[cols].values / sdv
    for i in range(len(a)):
        dist = np.sqrt(((B_ - A_[i]) ** 2).sum(1))
        for j in np.argsort(dist):
            if dist[j] > 0.5 * np.sqrt(3): break
            if j not in used: used.add(j); mask[a.index[i]] = True; mask[b.index[j]] = True; break
Sm = S[mask]; nonBm = Sm[~Sm.author_cell_type.isin(['B', 'PB'])]
B0m = Sm.groupby('sample')[['log_counts', 'log_genes', 'mito', 'ribo', 'stress']].median().reindex(samples)
B0m['ig'] = nonBm.groupby('sample').ig.median(); B0m['hb'] = nonBm.groupby('sample').hb.median(); B0m['log_cells'] = np.log(Sm.groupby('sample').size().reindex(samples)); B0m = B0m.fillna(0)
Rm = {n: round(roc_auc_score(y, cv(n, fold0, cell_mask=mask, B0u=B0m)), 3) for n in ['B0', 'B4', 'D1']}
nonlat_best = max(R[k]['auc'] for k in ['B1', 'B2', 'B3'])
r3 = dict(matched_cells=int(mask.sum()), **{f'{k}_matched': v for k, v in Rm.items()}, b0_check=bool(Rm['B0'] <= 0.6),
          d1_gain_survives=bool(Rm['D1'] - max(Rm['B4'], nonlat_best) >= 0.05)); r3['pass'] = bool(r3['b0_check'] and r3['d1_gain_survives'])
# rule 4: batch probe on training cells of first fold
d = sorted(lat)[0]; z = lat[d]; Ztr = z['scvi_tr']; ob = S.iloc[z['tr_idx']]
def skill(lab):
    lab = pd.factorize(lab)[0]; ch = 1 / len(np.unique(lab))
    pr = cross_val_predict(KNeighborsClassifier(15), Ztr, lab, cv=StratifiedKFold(5, shuffle=True, random_state=0))
    return round(float((balanced_accuracy_score(lab, pr) - ch) / (1 - ch)), 3)
r4 = dict(fold=d, label_skill=skill(ob.y.values), library_skill=skill(ob.library_uuid.astype(str).values)); r4['pass'] = bool(r4['label_skill'] >= r4['library_skill'])
# strict QC (B0-B2)
B0s, B1s, B2s, _ = basic('keep_strict'); strict = {}
for n, Fm in [('B0', B0s), ('B1', B1s), ('B2', B2s)]:
    pp = np.zeros(len(samples))
    for dd in np.unique(fold0):
        te = fold0 == dd; tr = ~te; Xf = prep(Fm.values, tr, False); s = StandardScaler().fit(Xf[tr])
        pp[te] = LogisticRegression(C=1.0, max_iter=5000).fit(s.transform(Xf[tr]), y[tr]).predict_proba(s.transform(Xf[te]))[:, 1]
    strict[n] = round(roc_auc_score(y, pp), 3)
# descriptive: IFN6 direction
ifn = dict(pos_mean=round(float(B1.values[y == 1].mean()), 3), neg_mean=round(float(B1.values[y == 0].mean()), 3))
res = dict(contrast=C, n_samples=len(samples), n_pos=int(y.sum()), results=R, rule1=r1, rule2=r2, rule3=r3, rule4=r4,
           H1=bool(r1['pass'] and r2['pass'] and r3['pass'] and r4['pass']), strict_qc=strict, ifn6_means=ifn, pb_types=pb_types, latent_types=types_all,
           scvi_epochs={k: int(v['epochs']) for k, v in lat.items()}, cells=int(len(Call)), removed=int((~Call.keep).sum()), removed_strict=int((~Call.keep_strict).sum()),
           protected_kept=int((Call.protected & Call.qc_fail).sum()), sub_cells=int(len(S)))
json.dump(res, open(OUT, 'w'), indent=1, default=str); print(json.dumps(res, indent=1, default=str))
