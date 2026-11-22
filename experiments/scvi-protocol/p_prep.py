"""Round 102 (Q36) stage 1 for Perez C3/C4: stream raw.X, per-cell QC metrics, MAD flags per library x donor, protected-population rule,
pseudobulk of kept cells, label-blind subsample. Follows PREREG-C3C4-round102.md.
Usage: python3 -I p_prep.py <h5ad> <work_dir>   (work_dir holds obs.pkl from p_design.py)"""
import sys, h5py, numpy as np, pandas as pd, scipy.sparse as sp, pickle, time
H, W = sys.argv[1:3]
f = h5py.File(H, 'r'); X = f['raw/X']; D, I, P = X['data'], X['indices'], X['indptr'][:]
v = f['raw/var']['feature_name']; sym = v['categories'][:].astype(str)[v['codes'][:]]; G = len(sym)
M = pd.read_pickle(f'{W}/obs.pkl'); M['row'] = np.arange(len(M))
c3d = M[M.disease_state.isin(['flare', 'treated'])].groupby('donor_id').disease_state.nunique(); c3d = c3d[c3d == 2].index
sel3 = (M.Processing_Cohort == '3.0') & M.donor_id.isin(c3d) & M.disease_state.isin(['flare', 'treated'])
sel4 = M.Processing_Cohort == '4.0'
M = M[sel3 | sel4].copy(); M['contrast'] = np.where(sel3[M.index], 'C3', 'C4')
M['sample'] = np.where(M.contrast == 'C3', M.donor_id.astype(str) + '|' + M.disease_state.astype(str), M.donor_id.astype(str))
M['y'] = np.where(M.contrast == 'C3', (M.disease_state == 'flare').astype(int), (M.disease == 'systemic lupus erythematosus').astype(int))
print(M.groupby('contrast').agg(cells=('row', 'size'), samples=('sample', 'nunique'), pos=('y', 'mean')))
def gi(gs): return np.flatnonzero(np.isin(sym, gs))
s = pd.Series(sym)
mt = np.flatnonzero(s.str.startswith('MT-').values); ribo = np.flatnonzero(s.str.match(r'^RP[SL]\d').values)
HB = gi(['HBA1', 'HBA2', 'HBB']); PL = gi(['PPBP', 'PF4'])
IG = np.flatnonzero(s.str.match(r'^(IGHG\d|IGKC|IGLC\d|JCHAIN|IGHA\d|IGHM)$').values)
ST = gi(['FOS', 'FOSB', 'JUN', 'JUNB', 'HSPA1A', 'HSPA1B', 'DUSP1', 'IER2'])
IFN = gi(['IFI27', 'IFI44L', 'IFIT1', 'ISG15', 'RSAD2', 'SIGLEC1']); assert len(IFN) == 6, len(IFN)
def ind(ix): z = np.zeros(G, np.float32); z[ix] = 1; return z
masks = np.stack([ind(mt), ind(ribo), ind(HB), ind(PL), ind(IG)], 1)
rows = np.sort(M.row.values); B = 20000
def blocks():
    for a in range(0, len(P) - 1, B):
        b = min(a + B, len(P) - 1); r = rows[(rows >= a) & (rows < b)]
        if len(r) == 0: continue
        lo, hi = P[a], P[b]
        Xb = sp.csr_matrix((D[lo:hi], I[lo:hi], P[a:b + 1] - lo), shape=(b - a, G))
        yield r, Xb[r - a]
t0 = time.time(); feats = []
for r, Xr in blocks():
    tot = np.asarray(Xr.sum(1)).ravel(); ng = np.diff(Xr.indptr)
    top20 = np.array([np.sort(Xr.data[Xr.indptr[i]:Xr.indptr[i + 1]])[-20:].sum() for i in range(Xr.shape[0])]) / tot
    fr = np.asarray(Xr @ masks) / tot[:, None]
    L = Xr.multiply(1e4 / tot[:, None]).tocsr(); L.data = np.log1p(L.data)
    feats.append(pd.DataFrame(dict(row=r, log_counts=np.log1p(tot), log_genes=np.log1p(ng), top20=top20, mito=fr[:, 0] * 100, ribo=fr[:, 1] * 100,
                                   hb=fr[:, 2], plt=fr[:, 3], ig=fr[:, 4], stress=np.asarray(L[:, ST].mean(1)).ravel(), ifn6=np.asarray(L[:, IFN].mean(1)).ravel())))
    print('pass1', r[-1], round(time.time() - t0), flush=True)
F = pd.concat(feats).set_index('row'); M = M.join(F, on='row')
def mad_out(x, n, side='both'):
    med = np.median(x); mad = np.median(np.abs(x - med)) * 1.4826 + 1e-9
    return (x > med + n * mad) | ((x < med - n * mad) if side == 'both' else False)
lm = np.log((M.mito / 100 + 1e-4) / (1 - M.mito / 100 + 1e-4)); M['lm'] = lm
def flags(g, n):
    return mad_out(g.log_counts.values, n) | mad_out(g.log_genes.values, n) | mad_out(g.top20.values, n) | (mad_out(g.lm.values, 3, 'upper') & (g.mito.values > 15))
for col, n in [('qc_fail', 5), ('qc_fail_strict', 3)]:
    M[col] = False
    for _, g in M.groupby(['library_uuid', 'donor_id'], observed=True): M.loc[g.index, col] = flags(g, n)
M['protected'] = False
for c, g in M.groupby('contrast'):
    base = g.qc_fail.mean(); rt = g.groupby('author_cell_type', observed=True).qc_fail.mean()
    prot = rt[rt > 3 * base].index.tolist(); print(c, 'overall flag', round(base, 4), 'by type', rt.round(4).to_dict(), 'protected', prot)
    M.loc[g.index[g.author_cell_type.isin(prot).values], 'protected'] = True
M['keep'] = ~M.qc_fail | M.protected; M['keep_strict'] = ~M.qc_fail_strict | M.protected
# label-blind subsample of kept cells
rng = np.random.default_rng(0); cap = {'C3': 100, 'C4': 30}; pick = []
for (c, smp, t), g in M[M.keep].groupby(['contrast', 'sample', 'author_cell_type'], observed=True):
    pick.append(rng.choice(g.row.values, size=min(cap[c], len(g)), replace=False))
pick = np.sort(np.concatenate(pick)); M['sub'] = M.row.isin(pick)
# pass 2: pseudobulk (kept cells) + subsample counts
key = M.contrast + '#' + M['sample'] + '#' + M.author_cell_type.astype(str); codes, uk = pd.factorize(key[M.keep]); keeprows = M.row[M.keep].values
cmap = pd.Series(codes, index=keeprows); PBm = np.zeros((len(uk), G), np.float64); subX = []; subrows = []
rows = np.sort(M.row.values); t0 = time.time()
for r, Xr in blocks():
    k = np.isin(r, keeprows); 
    if k.any():
        cc = cmap.loc[r[k]].values; A = sp.csr_matrix((np.ones(len(cc)), (cc, np.arange(len(cc)))), shape=(len(uk), len(cc)))
        PBm += (A @ Xr[k]).toarray()
    s_ = np.isin(r, pick)
    if s_.any(): subX.append(Xr[s_]); subrows.append(r[s_])
    print('pass2', r[-1], round(time.time() - t0), flush=True)
M.drop(columns=['lm']).to_pickle(f'{W}/cells.pkl')
pd.DataFrame(PBm.astype(np.int64), index=uk, columns=sym).to_pickle(f'{W}/pseudobulk.pkl')
sr = np.concatenate(subrows); assert (np.diff(sr) > 0).all()
sp.save_npz(f'{W}/sub_X.npz', sp.vstack(subX).tocsr()); M.set_index('row').loc[sr].to_pickle(f'{W}/sub_meta.pkl')
pd.Series(sym).to_csv(f'{W}/genes.csv', index=False)
print('done', len(M), 'cells; kept', int(M.keep.sum()), '; sub', len(sr))
