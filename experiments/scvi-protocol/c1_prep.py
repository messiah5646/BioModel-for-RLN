"""Round 100 C1 stage 1: per-sample QC, doublets, ambient/stress scores, pseudobulk, subsample. Follows PREREG-C1-round100.md.
Usage: python3 -I c1_prep.py <a319_dir> <out_dir>"""
import sys, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, scanpy as sc
D, O = sys.argv[1:3]; sc.settings.n_jobs = 2
F = {'T': '00260ebe-ef0d-400f-b937-8869459e2940', 'B': '65c2c3c1-150c-49c6-87a1-76ef4612a2b4', 'MY': '0f689b40-f76e-40bb-93d7-0cf1ec30d144',
     'MONO': 'd0bce338-4177-4e01-b5c4-62bcb1f4457c', 'LDG': '729e31d1-3d93-4a85-9602-df56638f54b2'}
A = {k: ad.read_h5ad(f'{D}/{v}.h5ad', backed='r') for k, v in F.items()}
my = set(A['MY'].obs_names)
ov = {k: float(np.mean([x in my for x in A[k].obs_names])) for k in ['MONO', 'LDG']}
print('fraction of barcodes inside Myeloid:', ov)
use = ['T', 'B', 'MY'] + [k for k in ['MONO', 'LDG'] if ov[k] < 0.5]
sym = A['MY'].raw.var['feature_name'].astype(str).to_numpy() if 'feature_name' in A['MY'].raw.var else A['MY'].var['feature_name'].astype(str).to_numpy()
assert all((A[k].raw.var.index == A['MY'].raw.var.index).all() for k in use)
def gi(gs): return np.flatnonzero(np.isin(sym, gs))
mt = np.flatnonzero(pd.Series(sym).str.startswith('MT-').to_numpy()); ribo = np.flatnonzero(pd.Series(sym).str.match(r'^RP[SL]\d').to_numpy())
HB = gi(['HBA1', 'HBA2', 'HBB']); PL = gi(['PPBP', 'PF4'])
IG = np.flatnonzero(pd.Series(sym).str.match(r'^(IGHG\d|IGKC|IGLC\d|JCHAIN|IGHA\d|IGHM)$').to_numpy())
ST = gi(['FOS', 'FOSB', 'JUN', 'JUNB', 'HSPA1A', 'HSPA1B', 'DUSP1', 'IER2'])
IFN = ['IFI27', 'IFI44L', 'IFIT1', 'ISG15', 'RSAD2', 'SIGLEC1']; IFNi = [np.flatnonzero(sym == g)[0] for g in IFN]
def major(ct):
    c = ct.lower()
    if 'b cell' in c or 'plasma' in c: return 'B_plasma'
    if 'natural killer' in c and 'nk t' not in c: return 'NK'
    if 'cd4' in c or 't-helper' in c or 'regulatory t' in c: return 'CD4_T'
    if 'cd8' in c: return 'CD8_T'
    if ' t cell' in c or c.startswith('t cell') or 'nk t' in c: return 'other_T'
    if 'monocyte' in c: return 'monocyte'
    if 'neutrophil' in c or 'granulocyte' in c: return 'neutrophil'
    if 'dendritic' in c: return 'DC'
    return 'platelet_other'

samples = sorted({s for k in use for s in A[k].obs.loc[A[k].obs.timepoint.isin(['Base', 'M3']), 'sample'].unique()})
import os, pickle
os.makedirs(f'{O}/per_sample', exist_ok=True)
if len(sys.argv) > 3 and sys.argv[3] == 'list': print(' '.join(samples)); sys.exit()
if len(sys.argv) > 3: samples = sys.argv[3:]
for s in samples:
    rows, pbs, subs, subX = [], {}, [], []
    if os.path.exists(f'{O}/per_sample/{s}.pkl'): continue
    Xs, Ms = [], []; rng = np.random.default_rng(int(''.join(c for c in s if c.isdigit())))  # fixed per-sample seed
    for k in use:
        o = A[k].obs; idx = np.flatnonzero((o['sample'] == s).to_numpy() & o.timepoint.isin(['Base', 'M3']).to_numpy())
        if len(idx) == 0: continue
        X = sp.csr_matrix(A[k].raw.X[idx]); Xs.append(X)
        m = o.iloc[idx][['sample', 'donor_id', 'timepoint', 'dose', 'cell_type']].copy(); m['file'] = k; Ms.append(m)
    X = sp.vstack(Xs).tocsr(); M = pd.concat(Ms)
    M['major'] = M.cell_type.astype(str).map(major)
    tot = np.asarray(X.sum(1)).ravel(); ng = np.diff(X.indptr)
    top20 = np.array([np.sort(X.data[X.indptr[i]:X.indptr[i+1]])[-20:].sum() for i in range(X.shape[0])]) / tot
    def frac(ix): return np.asarray(X[:, ix].sum(1)).ravel() / tot
    M['log_counts'] = np.log1p(tot); M['log_genes'] = np.log1p(ng); M['top20'] = top20
    M['mito'] = frac(mt) * 100; M['ribo'] = frac(ribo) * 100; M['hb'] = frac(HB); M['ig'] = frac(IG); M['plt'] = frac(PL)
    L = X.multiply(1e4 / tot[:, None]).tocsr(); L.data = np.log1p(L.data)
    M['stress'] = np.asarray(L[:, ST].mean(1)).ravel(); M['ifn6'] = np.asarray(L[:, IFNi].mean(1)).ravel()
    def mad_out(x, n, side='both'):
        med = np.median(x); mad = np.median(np.abs(x - med)) * 1.4826 + 1e-9
        return (x > med + n * mad) | ((x < med - n * mad) if side == 'both' else False)
    lm = np.log((M.mito / 100 + 1e-4) / (1 - M.mito / 100 + 1e-4))
    M['qc_fail'] = mad_out(M.log_counts.values, 5) | mad_out(M.log_genes.values, 5) | mad_out(M.top20.values, 5) | (mad_out(lm.values, 3, 'upper') & (M.mito.values > 10))
    try:
        sa = sc.AnnData(X.astype(np.float32)); sc.pp.scrublet(sa, random_state=0)
        M['doublet_score'] = sa.obs.doublet_score.values; M['doublet'] = sa.obs.predicted_doublet.fillna(False).values.astype(bool)
    except Exception as e:
        print(s, 'scrublet failed', e); M['doublet_score'] = np.nan; M['doublet'] = False
    keep = ~(M.qc_fail.values | M.doublet.values)
    for t, g in M[keep].groupby('major'):
        ii = np.flatnonzero(keep & (M.major.values == t))
        pbs[(s, t)] = np.asarray(X[ii].sum(0)).ravel()
        pick = np.sort(rng.choice(ii, size=min(150, len(ii)), replace=False))
        subs.append(M.iloc[pick].assign(cell_id=M.index[pick])); subX.append(X[pick])
    rows.append(M)
    pickle.dump(dict(M=M, pbs=pbs, subs=subs, subX=subX), open(f'{O}/per_sample/{s}.pkl', 'wb'))
    print(s, M.timepoint.iloc[0], len(M), 'qc_fail', int(M.qc_fail.sum()), 'doublet', int(M.doublet.sum()), flush=True)
if len(sys.argv) > 3: sys.exit()
P = [pickle.load(open(f'{O}/per_sample/{s}.pkl', 'rb')) for s in samples]
C = pd.concat([p['M'] for p in P]); pbs = {k: v for p in P for k, v in p['pbs'].items()}
subs = [x for p in P for x in p['subs']]; subX = [x for p in P for x in p['subX']]
C = pd.concat(rows); C.to_csv(f'{O}/cells.csv.gz')
pd.DataFrame({f'{s}|{t}': v for (s, t), v in pbs.items()}, index=sym).T.astype(np.int64).to_pickle(f'{O}/pseudobulk.pkl')
S = pd.concat(subs); S.to_csv(f'{O}/sub_meta.csv.gz'); sp.save_npz(f'{O}/sub_X.npz', sp.vstack(subX).tocsr())
pd.Series(sym).to_csv(f'{O}/genes.csv', index=False)
print('done', len(C), 'cells;', len(S), 'subsampled')
