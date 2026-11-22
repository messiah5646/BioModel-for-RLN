"""Round 104 (Q37/C5) evaluation per PREREG-C5-round104.md: frozen Perez-C3 classifiers, single pass on GSE254176; bridging-donor check.
Usage: python3 -I c5_eval.py <work_dir> <labels.tsv> <out_json>"""
import sys, json, numpy as np, pandas as pd, scipy.sparse as sp
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score, balanced_accuracy_score
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
from sklearn.model_selection import StratifiedKFold, cross_val_predict
W, LAB, OUT = sys.argv[1:4]; R = {}
genes = pd.read_csv(f'{W}/genes.csv').iloc[:, 0].astype(str).values
Pc = pd.read_pickle(f'{W}/perez_cells.pkl'); Pk = Pc[Pc.is_c3 & Pc.keep].rename(columns={'author_cell_type': 'ct'}); Pk['ct'] = Pk.ct.astype(str)
Ps = pd.read_pickle(f'{W}/perez_sub_meta.pkl'); P3 = Ps[Ps.sub3.values].reset_index(drop=True); P3['ct'] = P3.author_cell_type.astype(str)
PB = pd.read_pickle(f'{W}/perez_pb.pkl')
C = pd.read_pickle(f'{W}/gse_cells_annot.pkl').rename(columns={'ctype': 'ct'})
lab = pd.read_csv(LAB, sep='\t', header=None, names=['gsm', 'title', 'state']); lab['patient'] = lab.title.str.extract(r'(\d+), timepoint')[0].astype(int)
lab['tp'] = lab.title.str.extract(r'timepoint (\d+)')[0].astype(int); lab['y'] = (lab.state == 'Flare').astype(int); lab = lab.set_index('gsm')
L = np.load(f'{W}/c5_latent.npz'); SV = np.load(f'{W}/c5_sysvi.npz')
ps = np.array(sorted(Pk['sample'].unique())); yp = np.array([int(s.endswith('|flare')) for s in ps]); gs = lab.index.values; yg = lab.y.values
Ck = C[C.keep]
TYPES = sorted(Pk.ct.unique())
R['gse_type_counts'] = Ck.groupby(['gsm', 'ct']).size().unstack(fill_value=0).to_dict(orient='index')
R['gse_ctype_maxp_median'] = float(C.ctype_maxp.median()); R['gse_cells'] = dict(total=len(C), kept=int(C.keep.sum()), sub=int(C['sub'].sum()))
Cs = C[C['sub']]
sp_cnt = P3.groupby(['sample', 'ct']).size().unstack(fill_value=0).reindex(ps).fillna(0); sg_cnt = Cs.groupby(['gsm', 'ct']).size().unstack(fill_value=0).reindex(index=gs, columns=TYPES).fillna(0)
types_lat = [t for t in TYPES if (sp_cnt.get(t, 0) >= 10).all() and (sg_cnt[t] >= 10).all()]
kp_cnt = Pk.groupby(['sample', 'ct']).size().unstack(fill_value=0).reindex(ps).fillna(0); kg_cnt = Ck.groupby(['gsm', 'ct']).size().unstack(fill_value=0).reindex(index=gs, columns=TYPES).fillna(0)
types_pb = [t for t in TYPES if (kp_cnt.get(t, 0) >= 20).all() and (kg_cnt[t] >= 20).all()]
R['types_lat'] = types_lat; R['types_pb'] = types_pb
def b0(K, key, smp):
    nonB = K[~K.ct.isin(['B', 'PB'])]; f = K.groupby(key)[['log_counts', 'log_genes', 'mito', 'ribo', 'stress']].median()
    f['ig'] = nonB.groupby(key).ig.median(); f['hb'] = nonB.groupby(key).hb.median(); f['log_cells'] = np.log(K.groupby(key).size()); return f.reindex(smp).values
def b1(K, key, smp): return K[K.ct.isin(['cM', 'ncM'])].groupby(key).ifn6.mean().reindex(smp).values[:, None]
def b2(K, key, smp):
    c = K.groupby([key, 'ct']).size().unstack(fill_value=0).reindex(index=smp, columns=TYPES).fillna(0) + 0.5; l = np.log(c.values); return l - l.mean(1, keepdims=True)
# B3 pseudobulk
gpb = {}
for g in gs:
    X = sp.load_npz(f'{W}/gse/{g}.npz').tocsr(); cc = C[C.gsm == g]
    for t in types_pb: gpb[(g, t)] = np.asarray(X[cc.cell.values[(cc.keep & (cc.ct == t)).values]].sum(0)).ravel()
def lcpm(M): return np.log1p(M / M.sum(1, keepdims=True) * 1e6)
B3p, B3g = [], []
for t in types_pb:
    Lp = lcpm(PB.loc[[f'{s}#{t}' for s in ps]].values.astype(float)); Lg = lcpm(np.array([gpb[(g, t)] for g in gs], float))
    top = np.argsort(-Lp.var(0))[:2000]; B3p.append(Lp[:, top]); B3g.append(Lg[:, top])
B3p, B3g = np.hstack(B3p), np.hstack(B3g)
def lat_feat(Z, meta, key, smp, mask=None):
    out = []
    for s in smp:
        row = []
        for t in types_lat:
            r = np.flatnonzero((meta[key] == s).values & (meta.ct == t).values)
            if mask is not None: r = r[mask[r]]
            row.append(Z[r].mean(0) if len(r) else np.full(Z.shape[1], np.nan))
        out.append(np.concatenate(row))
    return np.array(out)
Cfull = C.reset_index(drop=True); subidx = np.flatnonzero(Cfull['sub'].values); Csub = Cfull.iloc[subidx].reset_index(drop=True)
zs_g = np.full((len(Cfull), SV['g'].shape[1]), np.nan); zs_g[SV['gidx']] = SV['g']
F = {'B0': (b0(Pk, 'sample', ps), b0(Ck, 'gsm', gs)), 'B1': (b1(Pk, 'sample', ps), b1(Ck, 'gsm', gs)), 'B2': (b2(Pk, 'sample', ps), b2(Ck, 'gsm', gs)),
     'B3': (B3p, B3g), 'B4': (lat_feat(L['pca_p'], P3, 'sample', ps), lat_feat(L['pca_g'][subidx], Csub, 'gsm', gs)),
     'D1': (lat_feat(L['scvi_p'], P3, 'sample', ps), lat_feat(L['scvi_g'][subidx], Csub, 'gsm', gs)),
     'D2': (lat_feat(SV['p'], P3, 'sample', ps), lat_feat(zs_g[subidx], Csub, 'gsm', gs))}
HI = {'B3', 'B4', 'D1', 'D2'}
def fit_apply(Xp, Xg, hi, domain_std=False):
    if domain_std:
        Xp = (Xp - np.nanmean(Xp, 0)) / (np.nanstd(Xp, 0) + 1e-9); Xg = (Xg - np.nanmean(Xg, 0)) / (np.nanstd(Xg, 0) + 1e-9)
    mu = np.nanmean(Xp, 0); Xp = np.nan_to_num(Xp - mu); Xg = np.nan_to_num(Xg - mu)
    if hi: pc = PCA(10, random_state=0).fit(Xp); Xp, Xg = pc.transform(Xp), pc.transform(Xg)
    s = StandardScaler().fit(Xp); m = LogisticRegression(C=1.0, max_iter=5000).fit(s.transform(Xp), yp)
    return m.predict_proba(s.transform(Xg))[:, 1], m.predict_proba(s.transform(Xp))[:, 1]
def b5(Xp, Xg):
    s = StandardScaler().fit(Xp); nn = NearestNeighbors(n_neighbors=3).fit(s.transform(Xp)); return yp[nn.kneighbors(s.transform(Xg))[1]].mean(1)
EP = [(1, 1, 2), (2, 1, 2), (3, 1, 2), (3, 4, 5), (4, 2, 3), (4, 5, 6)]
gidx = {(r.patient, r.tp): i for i, r in enumerate(lab.itertuples())}; pat = lab.patient.values
def metrics(p):
    out = dict(auc=round(roc_auc_score(yg, p), 3))
    d = [p[gidx[(a, f)]] - p[gidx[(a, b)]] for a, b, f in EP]; out['paired_k_of_6'] = int(np.sum(np.array(d) > 0))
    m2 = ~np.isin(pat, [5, 6]); pc = p - pd.Series(p).groupby(pat).transform('mean').values
    out['within_patient_auc'] = round(roc_auc_score(yg[m2], pc[m2]), 3); return out
rng = np.random.default_rng(0); up = np.unique(pat); BOOT = [rng.choice(up, len(up)) for _ in range(2000)]
def boot_auc(p, q=None):
    v = []
    for bs in BOOT:
        ix = np.concatenate([np.flatnonzero(pat == u) for u in bs]); yy = yg[ix]
        if yy.min() == yy.max(): continue
        a = roc_auc_score(yy, p[ix]); v.append(a - (roc_auc_score(yy, q[ix]) if q is not None else 0))
    return [round(float(np.percentile(v, 2.5)), 3), round(float(np.percentile(v, 97.5)), 3)]
Pg, Ptr = {}, {}
for n, (xp, xg) in F.items(): Pg[n], Ptr[n] = fit_apply(xp, xg, n in HI)
Pg['B5'] = b5(F['B2'][0], F['B2'][1])
mu1 = np.nanmean(F['D1'][0], 0); pc1 = PCA(10, random_state=0).fit(np.nan_to_num(F['D1'][0] - mu1))
XD = lambda X: pc1.transform(np.nan_to_num(X - mu1))
Pg['D1+B0'] = fit_apply(np.hstack([XD(F['D1'][0]), F['B0'][0]]), np.hstack([XD(F['D1'][1]), F['B0'][1]]), False)[0]
R['train_in_sample_auc'] = {n: round(roc_auc_score(yp, v), 3) for n, v in Ptr.items()}
R['test'] = {n: dict(**metrics(p), ci=boot_auc(p)) for n, p in Pg.items()}
R['test_domain_std'] = {n: metrics(fit_apply(xp_, xg_, n in HI, True)[0]) for n, (xp_, xg_) in F.items()}
base = ['B1', 'B2', 'B3', 'B4', 'B5']; best = max(base, key=lambda n: R['test'][n]['auc'])
rule = {}
for D in ['D1', 'D2']:
    gain = R['test'][D]['auc'] - R['test'][best]['auc']; ci = boot_auc(Pg[D], Pg[best])
    rule[D] = dict(best_baseline=best, gain=round(gain, 3), ci=ci, rule1=bool(gain >= 0.05 and ci[0] > 0),
                   rule2=bool(R['test'][D]['auc'] - R['test']['B0']['auc'] >= 0.05 and R['test']['D1+B0']['auc'] > R['test']['B0']['auc']))
# rule 3: QC matching within patient x type, Flare cells 1:1 to non-Flare cells of same patient on 9 metrics
met = ['log_counts', 'log_genes', 'top20', 'mito', 'ribo', 'stress', 'hb', 'plt', 'ig']
K = Cfull[Cfull.keep].copy(); K['patient'] = K.gsm.map(lab.patient); K['flare'] = K.gsm.map(lab.y)
Z = (K[met] - K[met].mean()) / (K[met].std() + 1e-12); cal = 0.5 * np.sqrt(len(met)); matched = pd.Series(False, index=K.index)
for (pt, t), g in K.groupby(['patient', 'ct']):
    if pt in (5, 6): matched[g.index] = True; continue
    f = g.index[g.flare == 1]; nf = g.index[g.flare == 0]
    if len(f) == 0 or len(nf) == 0: continue
    nn = NearestNeighbors(n_neighbors=min(len(nf), 30)).fit(Z.loc[nf].values); dist, ind = nn.kneighbors(Z.loc[f].values); used = set()
    for i in np.argsort(dist[:, 0]):
        for dd, j in zip(dist[i], ind[i]):
            if dd > cal: break
            if j not in used: used.add(j); matched[f[i]] = True; matched[nf[j]] = True; break
mm = np.zeros(len(Cfull), bool); mm[matched.index[matched.values]] = True
Km = Cfull[mm]; sm = mm[subidx]
R['matching'] = dict(kept=len(K), matched=int(mm.sum()), per_sample=Km.groupby('gsm').size().to_dict())
Fm = {'B0': b0(Km, 'gsm', gs), 'B4': lat_feat(L['pca_g'][subidx], Csub, 'gsm', gs, sm), 'D1': lat_feat(L['scvi_g'][subidx], Csub, 'gsm', gs, sm),
      'B1': b1(Km, 'gsm', gs), 'B2': b2(Km, 'gsm', gs)}
Pm = {n: fit_apply(F[n][0], Fm[n], n in HI)[0] for n in Fm}; Pm['B5'] = b5(F['B2'][0], Fm['B2'])
R['test_matched'] = {n: metrics(p) for n, p in Pm.items()}
bm = max(['B1', 'B2', 'B4', 'B5'], key=lambda n: R['test_matched'][n]['auc'])
rule['D1']['rule3'] = bool(R['test_matched']['B0']['auc'] <= 0.6 and R['test_matched']['D1']['auc'] - R['test_matched'][bm]['auc'] >= 0.05)
rule['D1']['rule3_detail'] = dict(B0_matched=R['test_matched']['B0']['auc'], best_matched=bm, gain_matched=round(R['test_matched']['D1']['auc'] - R['test_matched'][bm]['auc'], 3))
# rule 4: batch probe on GSE sub cells in D1 latent
zsub = L['scvi_g'][subidx]; pt_lab = Csub.gsm.map(lab.patient).values; fl = Csub.gsm.map(lab.y).values
def skill(yv):
    pr = cross_val_predict(KNeighborsClassifier(15), zsub, yv, cv=StratifiedKFold(5, shuffle=True, random_state=0)); k = len(np.unique(yv))
    return round(float((balanced_accuracy_score(yv, pr) - 1 / k) / (1 - 1 / k)), 3)
sp_, sl_ = skill(pt_lab), skill(fl); rule['D1']['rule4'] = bool(sl_ >= sp_); rule['D1']['rule4_detail'] = dict(patient_skill=sp_, label_skill=sl_)
rule['H1'] = bool(all(rule['D1'][k] for k in ['rule1', 'rule2', 'rule3', 'rule4']))
R['rules'] = rule
R['transfers'] = {n: bool(v['auc'] >= 0.70 and v['paired_k_of_6'] >= 5) for n, v in R['test'].items()}
R['test_scores'] = {n: dict(zip(gs, np.round(p, 3).tolist())) for n, p in Pg.items()}
# bridging-donor check (descriptive)
try:
    BL = np.load(f'{W}/bridge_latent.npz'); Bm = Ps[Ps.subb.values].reset_index(drop=True); Bm['ct'] = Bm.author_cell_type.astype(str)
    bsm = np.array(sorted(Bm.bsample.unique())); cnt = Bm.groupby(['bsample', 'ct']).size().unstack(fill_value=0).reindex(bsm).fillna(0)
    tb = [t for t in cnt.columns if (cnt[t] >= 5).all()]
    don = np.array([s.split('|c')[0] for s in bsm]); coh = np.array([s.split('|c')[1] for s in bsm])
    def vec(Z):
        return np.array([np.concatenate([Z[(Bm.bsample == s).values & (Bm.ct == t).values].mean(0) for t in tb]) for s in bsm])
    Bk = Pc[Pc.is_br & Pc.keep].rename(columns={'author_cell_type': 'ct'}); Bk['ct'] = Bk.ct.astype(str)
    reps = {k: vec(BL[k]) for k in ['pca', 'scvi', 'sysvi']}; reps['B2'] = b2(Bk, 'bsample', bsm)
    out = {}
    for k, V in reps.items():
        V = StandardScaler().fit_transform(V); Vn = V / np.linalg.norm(V, axis=1, keepdims=True); S_ = Vn @ Vn.T; ranks = []
        for i in range(len(bsm)):
            for c2 in set(coh[don == don[i]]) - {coh[i]}:
                cand = np.flatnonzero(coh == c2); tgt = np.flatnonzero((coh == c2) & (don == don[i]))[0]
                ranks.append(int((S_[i, cand] > S_[i, tgt]).sum()) + 1)
        ranks = np.array(ranks); pr = cross_val_predict(KNeighborsClassifier(5), V, coh, cv=StratifiedKFold(5, shuffle=True, random_state=0)); kk = len(np.unique(coh))
        out[k] = dict(top1=round(float((ranks == 1).mean()), 3), median_rank=float(np.median(ranks)), n_queries=len(ranks),
                      cohort_probe_skill=round(float((balanced_accuracy_score(coh, pr) - 1 / kk) / (1 - 1 / kk)), 3))
    R['bridging'] = dict(types=tb, n_bsamples=len(bsm), cohort_sizes=pd.Series(coh).value_counts().to_dict(), result=out)
except FileNotFoundError: R['bridging'] = 'not run'
json.dump(R, open(OUT, 'w'), indent=1, default=str); print(json.dumps({k: R[k] for k in ['types_lat', 'types_pb', 'test', 'test_matched', 'rules', 'transfers', 'bridging', 'train_in_sample_auc']}, indent=1, default=str))
