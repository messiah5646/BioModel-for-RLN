"""EXPLORATORY (not pre-registered): technical-only baseline from QC covariates. Usage: python -I exploratory_qc.py <mono.h5ad> <q3_dir> <outdir>"""
import sys, numpy as np, pandas as pd, anndata as ad
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
h5, d, out = sys.argv[1:4]
a = ad.read_h5ad(h5, backed='r'); meta = pd.read_csv(f'{d}/meta.csv', index_col=0)
q = a.obs.loc[meta.index, ['nCount_RNA', 'nFeature_RNA', 'percent.mito']].astype(float)
q['nCount_RNA'] = np.log1p(q['nCount_RNA']); q['nFeature_RNA'] = np.log1p(q['nFeature_RNA'])
Q = q.to_numpy(); y = meta['y'].to_numpy(); don = meta['donor_id'].to_numpy(); rows = []
for dn in sorted(set(don), key=lambda s: int(s[1:])):
    te = don == dn; tr = ~te
    m = make_pipeline(StandardScaler(), LogisticRegression(class_weight='balanced', max_iter=5000)).fit(Q[tr], y[tr])
    rows.append({'donor': dn, 'F4_qc_only': roc_auc_score(y[te], m.predict_proba(Q[te])[:, 1])})
R = pd.DataFrame(rows); R.to_csv(f'{out}/exploratory_qc_auc.csv', index=False); print(R.round(3).to_string(index=False)); print('mean', R.F4_qc_only.mean().round(3))
print(q.assign(tp=meta.timepoint).groupby('tp').median().round(3))
