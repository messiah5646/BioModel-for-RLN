"""Merge per-sample outputs of c1_prep.py. Usage: python3 -I c1_aggregate.py <a319_dir> <out_dir>"""
import sys, glob, pickle, numpy as np, pandas as pd, scipy.sparse as sp, h5py
D, O = sys.argv[1:3]
with h5py.File(f'{D}/0f689b40-f76e-40bb-93d7-0cf1ec30d144.h5ad', 'r') as f:
    g = f['raw/var']; k = 'feature_name'
    sym = g[k]['categories'][:][g[k]['codes'][:]].astype(str) if isinstance(g[k], h5py.Group) else g[k][:].astype(str)
Ms, subs, subX, pbs = [], [], [], {}
for fn in sorted(glob.glob(f'{O}/per_sample/*.pkl')):
    p = pickle.load(open(fn, 'rb')); Ms.append(p['M']); subs += p['subs']; subX += p['subX']; pbs.update(p['pbs']); del p
    print(fn.split('/')[-1], flush=True)
pd.concat(Ms).to_csv(f'{O}/cells.csv.gz')
pd.DataFrame(np.vstack(list(pbs.values())), index=[f'{s}|{t}' for s, t in pbs], columns=sym).astype(np.int64).to_pickle(f'{O}/pseudobulk.pkl')
S = pd.concat(subs); S.to_csv(f'{O}/sub_meta.csv.gz'); sp.save_npz(f'{O}/sub_X.npz', sp.vstack(subX).tocsr())
pd.Series(sym).to_csv(f'{O}/genes.csv', index=False); print('cells', sum(len(m) for m in Ms), 'sub', len(S))
