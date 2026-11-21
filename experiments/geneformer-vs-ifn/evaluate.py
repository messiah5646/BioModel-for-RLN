"""Leave-one-donor-out AUC for F1 (IFN score), F2 (PCA-50 of HVG), F3 (Geneformer). Usage: python -I evaluate.py <q3_dir> <outdir>"""
import sys, os, json, numpy as np, pandas as pd, scipy.sparse as sp
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score
d, out = sys.argv[1:3]
meta = pd.read_csv(f'{d}/meta.csv', index_col=0); y = meta['y'].to_numpy(); don = meta['donor_id'].to_numpy()
X = sp.load_npz(f'{d}/X_lognorm.npz').toarray().astype(np.float32)
E = np.load(f'{d}/emb.npy') if os.path.exists(f'{d}/emb.npy') else None
clf = lambda: make_pipeline(StandardScaler(), LogisticRegression(C=1.0, class_weight='balanced', max_iter=5000))
rows = []
for dn in sorted(set(don), key=lambda s: int(s[1:])):
    te = don == dn; tr = ~te; r = {'donor': dn, 'n_test': int(te.sum())}
    f1 = meta[['ifn']].to_numpy()
    r['F1_ifn'] = roc_auc_score(y[te], clf().fit(f1[tr], y[tr]).predict_proba(f1[te])[:, 1])
    hv = np.argsort(-X[tr].var(0))[:2000]; p = PCA(50, random_state=0).fit(X[tr][:, hv])
    Z = p.transform(X[:, hv])
    r['F2_pca50'] = roc_auc_score(y[te], clf().fit(Z[tr], y[tr]).predict_proba(Z[te])[:, 1])
    if E is not None:
        r['F3_geneformer'] = roc_auc_score(y[te], clf().fit(E[tr], y[tr]).predict_proba(E[te])[:, 1])
    rows.append(r)
R = pd.DataFrame(rows); R.to_csv(f'{out}/lodo_auc.csv', index=False)
print(R.round(3).to_string(index=False)); print(R.drop(columns=['donor', 'n_test']).mean().round(3))
if E is not None:
    m = R.drop(columns=['donor', 'n_test']).mean()
    c1 = m['F3_geneformer'] - m['F1_ifn'] >= 0.05
    c2 = int((R['F3_geneformer'] > R['F1_ifn']).sum())
    c3 = m['F3_geneformer'] >= m['F2_pca50'] - 0.02
    v = {'mean': m.round(4).to_dict(), 'c1_delta': round(m['F3_geneformer'] - m['F1_ifn'], 4), 'c1': bool(c1),
         'c2_donors_F3_gt_F1': c2, 'c2': c2 >= 8, 'c3': bool(c3), 'KEPT': bool(c1 and c2 >= 8 and c3)}
    json.dump(v, open(f'{out}/verdict.json', 'w'), indent=1); print(v)
