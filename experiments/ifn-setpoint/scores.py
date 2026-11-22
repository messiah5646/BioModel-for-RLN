"""Per-sample monocyte 6-gene IFN score for CD3xCD19 (CELLxGENE) and CAR-T (GSE263931). Writes scores.csv."""
import sys, glob, os, numpy as np, pandas as pd, scipy.io, scipy.sparse as sp, anndata as ad
G = ['IFI27','IFI44L','IFIT1','ISG15','RSAD2','SIGLEC1']
mono_path, cart_dir, out = sys.argv[1], sys.argv[2], sys.argv[3]
rows = []
a = ad.read_h5ad(mono_path, backed='r')
sym = a.var['feature_name'].astype(str).values if 'feature_name' in a.var else a.var_names.values
idx = [int(np.where(sym == g)[0][0]) for g in G]
parts = []
for i in range(0, a.n_obs, 5000):
    Xc = a.X[i:i+5000]
    Xc = Xc[:, idx]
    parts.append(np.asarray(Xc.toarray() if sp.issparse(Xc) else Xc).mean(1))
cs = np.concatenate(parts)
keepm = (a.obs['anno.fig.ed8'] != 'Granulocyte').values
df = pd.DataFrame({'sample': a.obs['sample'].astype(str).values, 'donor': a.obs['donor_id'].astype(str).values,
                   'tp': a.obs['timepoint'].astype(str).values, 'dose': a.obs['dose'].astype(str).values, 's': cs})[keepm]
for (smp, d, tp, dose), g in df.groupby(['sample','donor','tp','dose']):
    rows.append(dict(trial='CD3xCD19', donor=d, sample=smp, tp=tp, dose=dose, n=len(g), score=g.s.mean()))
for m in sorted(glob.glob(os.path.join(cart_dir, '*_matrix.mtx.gz'))):
    base = m[:-len('_matrix.mtx.gz')]; name = os.path.basename(base).split('_', 1)[1]  # pat1_pre
    pat, tp = name.split('_')
    M = scipy.io.mmread(m).tocsc()  # genes x cells
    feat = pd.read_csv(base + '_features.tsv.gz', sep='\t', header=None)
    s = feat[1].astype(str).values
    tot = np.asarray(M.sum(0)).ravel(); keep = tot > 500
    M = M[:, keep]; tot = tot[keep]
    def lg(g):
        r = M[np.where(s == g)[0][0], :].toarray().ravel()
        return np.log1p(r / tot * 1e4)
    mono = lg('LYZ') > 2
    sc = np.mean([lg(g) for g in G], axis=0)
    rows.append(dict(trial='CART', donor=pat, sample=name, tp=tp, dose='', n=int(mono.sum()), score=sc[mono].mean()))
pd.DataFrame(rows).to_csv(out, index=False)
print(pd.DataFrame(rows).to_string())
