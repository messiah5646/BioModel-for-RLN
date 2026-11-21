"""Select cells, compute F1 inputs and Geneformer tokens. Usage: python -I prep.py <mono.h5ad> <gf_dir> <outdir>"""
import sys, pickle, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
h5, gfd, out = sys.argv[1:4]
a = ad.read_h5ad(h5, backed='r')
obs = a.obs
keep = obs['timepoint'].isin(['Base', 'M3']) & (obs['anno.fig.ed8'] != 'Granulocyte')
rng = np.random.default_rng(0)
idx = []
for s, g in obs[keep].groupby('sample', observed=True):
    ii = np.flatnonzero(obs.index.isin(g.index))
    idx += sorted(rng.choice(ii, size=min(100, len(ii)), replace=False).tolist())
idx = np.array(sorted(idx))
X = a.X[idx]; X = sp.csr_matrix(X)
R = sp.csr_matrix(a.raw.X[idx])
meta = obs.iloc[idx][['sample', 'donor_id', 'timepoint', 'dose', 'anno.fig.ed8']].copy()
meta['y'] = (meta['timepoint'] == 'M3').astype(int)
print('cells', len(idx), 'raw integer:', np.allclose(R.data[:1000], np.round(R.data[:1000])))
var = a.var.index.to_numpy(); rawvar = a.raw.var.index.to_numpy()
sym = a.var['feature_name'].astype(str).to_numpy()
ifn = ['IFI27', 'IFI44L', 'IFIT1', 'ISG15', 'RSAD2', 'SIGLEC1']
cols = [np.flatnonzero(sym == g)[0] for g in ifn]
meta['ifn'] = np.asarray(X[:, cols].mean(axis=1)).ravel()
meta.to_csv(f'{out}/meta.csv')
sp.save_npz(f'{out}/X_lognorm.npz', X)
tok = pickle.load(open(f'{gfd}/geneformer/token_dictionary_gc104M.pkl', 'rb'))
med = pickle.load(open(f'{gfd}/geneformer/gene_median_dictionary_gc104M.pkl', 'rb'))
ok = np.array([g in tok and g in med for g in rawvar])
print('raw genes', len(rawvar), 'in Geneformer vocab', ok.sum())
gid = rawvar[ok]; mvec = np.array([med[g] for g in gid]); tvec = np.array([tok[g] for g in gid])
Rk = R[:, np.flatnonzero(ok)].tocsr()
tot = np.asarray(R.sum(axis=1)).ravel()  # total counts over all genes, as in the Geneformer tokenizer
cls, eos = tok['<cls>'], tok['<eos>']
seqs = []
for i in range(Rk.shape[0]):
    s, e = Rk.indptr[i], Rk.indptr[i + 1]
    j = Rk.indices[s:e]; v = Rk.data[s:e] / tot[i] * 1e4 / mvec[j]
    order = j[np.argsort(-v, kind='stable')]
    seqs.append([cls] + tvec[order][:2046].tolist() + [eos])
L = [len(x) for x in seqs]; print('token length median', np.median(L), 'max', max(L), 'truncated', sum(l >= 2048 for l in L))
pickle.dump(seqs, open(f'{out}/tokens.pkl', 'wb'))
