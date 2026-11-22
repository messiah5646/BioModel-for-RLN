"""Round 98 / Q13: QC-adjusted LODO AUC (see PREREG-Q13-qc-adjusted.md).
Usage: python -I evaluate_q13.py <mono.h5ad> <q13_dir> <outdir>"""
import sys, json, numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression, LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
h5, d, out = sys.argv[1:4]
meta = pd.read_csv(f'{d}/meta.csv', index_col=0); y = meta['y'].to_numpy(); don = meta['donor_id'].to_numpy()
a = ad.read_h5ad(h5, backed='r')
q = a.obs.loc[meta.index, ['nCount_RNA', 'nFeature_RNA', 'percent.mito']].astype(float)
Q = np.column_stack([np.log1p(q.nCount_RNA), np.log1p(q.nFeature_RNA), q['percent.mito']])
X = sp.load_npz(f'{d}/X_lognorm.npz').toarray().astype(np.float32)
E = np.load(f'{d}/emb.npy'); assert not np.isnan(E).any()
F1 = meta[['ifn']].to_numpy()
clf = lambda: make_pipeline(StandardScaler(), LogisticRegression(C=1.0, class_weight='balanced', max_iter=5000))
def auc(F, tr, te): return roc_auc_score(y[te], clf().fit(F[tr], y[tr]).predict_proba(F[te])[:, 1])
def resid(F, tr):
    lr = LinearRegression().fit(Q[tr], F[tr]); return F - lr.predict(Q)
donors = sorted(set(don), key=lambda s: int(s[1:])); rows = []
for dn in donors:
    te = don == dn; tr = ~te; r = {'donor': dn}
    hv = np.argsort(-X[tr].var(0))[:2000]; F2 = PCA(50, random_state=0).fit(X[tr][:, hv]).transform(X[:, hv])
    feats = {'F1_ifn': F1, 'F2_pca50': F2, 'F3_gf': E}
    r['F4_qc'] = auc(Q, tr, te)
    for k, F in feats.items():
        r[k + '_raw'] = auc(F, tr, te)
        r[k + '_resid'] = auc(resid(F, tr), tr, te)
        r[k + '_plusQC'] = auc(np.column_stack([Q, F]), tr, te)
    rows.append(r)
R = pd.DataFrame(rows).set_index('donor'); R.to_csv(f'{out}/q13_lodo_auc.csv')
m = R.mean(); print(R.round(3).to_string()); print(m.round(3).to_string())
# A3 QC-matched subsample (secondary)
Z = (Q - Q.mean(0)) / Q.std(0); keep = np.zeros(len(y), bool); npairs = {}
for dn in donors:
    ii = np.flatnonzero(don == dn); m3 = ii[y[ii] == 1]; b = list(ii[y[ii] == 0]); n = 0
    for i in m3:
        if not b: break
        dist = np.linalg.norm(Z[b] - Z[i], axis=1); j = int(np.argmin(dist))
        if dist[j] <= 0.25: keep[i] = keep[b[j]] = True; b.pop(j); n += 1
    npairs[dn] = n
rows3 = []
for dn in donors:
    if npairs[dn] < 20: continue
    te = (don == dn) & keep; tr = (don != dn) & keep; r = {'donor': dn, 'pairs': npairs[dn]}
    hv = np.argsort(-X[tr].var(0))[:2000]; F2 = PCA(50, random_state=0).fit(X[tr][:, hv]).transform(X[:, hv])
    for k, F in {'F4_qc': Q, 'F1_ifn': F1, 'F2_pca50': F2, 'F3_gf': E}.items(): r[k] = auc(F, tr, te)
    rows3.append(r)
R3 = pd.DataFrame(rows3); R3.to_csv(f'{out}/q13_matched_auc.csv', index=False)
print('matched pairs', npairs); print(R3.round(3).to_string(index=False)); print(R3.drop(columns=['donor']).mean().round(3))
# verdicts
h1a = m['F1_ifn_resid'] >= 0.60; inc = R['F1_ifn_plusQC'] - R['F4_qc']
h1b = inc.mean() >= 0.03 and (inc > 0).sum() >= 7
c1 = m['F3_gf_resid'] - m['F1_ifn_resid'] >= 0.05; c2 = int((R['F3_gf_resid'] > R['F1_ifn_resid']).sum())
c3 = m['F3_gf_resid'] >= m['F2_pca50_resid'] - 0.02
edge = m['F2_pca50_resid'] - m['F1_ifn_resid']
v = {'mean': m.round(4).to_dict(),
     'H1': {'resid_ifn': round(m['F1_ifn_resid'], 4), 'a': bool(h1a), 'inc_mean': round(inc.mean(), 4),
            'inc_donors_up': int((inc > 0).sum()), 'b': bool(h1b), 'PASS': bool(h1a and h1b)},
     'H2': {'delta': round(m['F3_gf_resid'] - m['F1_ifn_resid'], 4), 'c1': bool(c1), 'donors': c2, 'c2': c2 >= 8,
            'c3': bool(c3), 'KEPT': bool(c1 and c2 >= 8 and c3)},
     'H3': {'edge_resid': round(edge, 4), 'edge_raw': round(m['F2_pca50_raw'] - m['F1_ifn_raw'], 4),
            'largely_technical': bool(edge <= 0.069)},
     'A3_matched_mean': R3.drop(columns=['donor']).mean().round(4).to_dict(), 'A3_pairs': npairs}
json.dump(v, open(f'{out}/q13_verdict.json', 'w'), indent=1, default=int); print(json.dumps(v, indent=1, default=int))
