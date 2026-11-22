"""Round 104 (Q37/C5) Perez side: C3 training cells (cohort 3, flare/treated pairs) + bridging-donor cells (donors in >1 Processing_Cohort).
Stage-1 QC as p_prep.py (5 MAD, mito 3 MAD & >15%, per library x donor, protected-population rule), label-blind subsample
(C3 <= 100 / sample x type; bridging <= 30 / sample x type), pseudobulk of kept C3 cells. Genes restricted to the GSE254176 symbol set.
Usage: python3 -I c5_prep_perez.py <h5ad> <gse_features.tsv.gz> <work_dir>"""
import sys, os, h5py, numpy as np, pandas as pd, scipy.sparse as sp, time
H, FEAT, W = sys.argv[1:4]; os.makedirs(W, exist_ok=True)
f = h5py.File(H, 'r'); X = f['raw/X']; D, I, P = X['data'], X['indices'], X['indptr'][:]
v = f['raw/var']['feature_name']; sym = v['categories'][:].astype(str)[v['codes'][:]]; G = len(sym)
g = f['obs']
def col(k):
    x = g[k]
    if isinstance(x, h5py.Group): return pd.Categorical.from_codes(x['codes'][:], x['categories'][:].astype(str))
    return x[:].astype(str) if x.dtype.kind in 'SO' else x[:]
M = pd.DataFrame({k: col(k) for k in ['Processing_Cohort', 'disease', 'disease_state', 'donor_id', 'library_uuid', 'author_cell_type']})
M['row'] = np.arange(len(M))
c3d = M[M.disease_state.isin(['flare', 'treated'])].groupby('donor_id', observed=True).disease_state.nunique(); c3d = c3d[c3d == 2].index
sel3 = (M.Processing_Cohort == '3.0') & M.donor_id.isin(c3d) & M.disease_state.isin(['flare', 'treated'])
br = M.groupby('donor_id', observed=True).Processing_Cohort.nunique(); brd = br[br > 1].index
selb = M.donor_id.isin(brd)
M = M[sel3 | selb].copy(); M['is_c3'] = sel3[M.index].values; M['is_br'] = selb[M.index].values
M['sample'] = np.where(M.is_c3, M.donor_id.astype(str) + '|' + M.disease_state.astype(str), '')
M['bsample'] = M.donor_id.astype(str) + '|c' + M.Processing_Cohort.astype(str)
M['y'] = (M.disease_state == 'flare').astype(int)
print('C3 cells', int(M.is_c3.sum()), 'samples', M[M.is_c3]['sample'].nunique(), '| bridging donors', len(brd), 'cells', int(M.is_br.sum()),
      'bsamples', M[M.is_br].bsample.nunique(), flush=True)
gs = pd.read_csv(FEAT, sep='\t', header=None)[1].astype(str)
s = pd.Series(sym); keepg = np.flatnonzero(s.isin(set(gs)).values & ~s.duplicated().values); print('genes shared', len(keepg), 'of', G)
def gi(x): return np.flatnonzero(np.isin(sym, x))
mt = np.flatnonzero(s.str.startswith('MT-').values); ribo = np.flatnonzero(s.str.match(r'^RP[SL]\d').values)
HB = gi(['HBA1', 'HBA2', 'HBB']); PL = gi(['PPBP', 'PF4'])
IG = np.flatnonzero(s.str.match(r'^(IGHG\d|IGKC|IGLC\d|JCHAIN|IGHA\d|IGHM)$').values)
ST = gi(['FOS', 'FOSB', 'JUN', 'JUNB', 'HSPA1A', 'HSPA1B', 'DUSP1', 'IER2'])
IFN = gi(['IFI27', 'IFI44L', 'IFIT1', 'ISG15', 'RSAD2', 'SIGLEC1']); assert len(IFN) == 6
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
    if len(feats) % 10 == 0: print('pass1', r[-1], round(time.time() - t0), flush=True)
F = pd.concat(feats).set_index('row'); M = M.join(F, on='row')
def mad_out(x, n, side='both'):
    med = np.median(x); mad = np.median(np.abs(x - med)) * 1.4826 + 1e-9
    return (x > med + n * mad) | ((x < med - n * mad) if side == 'both' else False)
M['lm'] = np.log((M.mito / 100 + 1e-4) / (1 - M.mito / 100 + 1e-4))
M['qc_fail'] = False
for _, gg in M.groupby(['library_uuid', 'donor_id'], observed=True):
    M.loc[gg.index, 'qc_fail'] = mad_out(gg.log_counts.values, 5) | mad_out(gg.log_genes.values, 5) | mad_out(gg.top20.values, 5) | (mad_out(gg.lm.values, 3, 'upper') & (gg.mito.values > 15))
base = M.qc_fail.mean(); rt = M.groupby('author_cell_type', observed=True).qc_fail.mean(); prot = rt[rt > 3 * base].index.tolist()
print('flag rate', round(base, 4), 'protected', prot, flush=True)
M['keep'] = ~M.qc_fail | M.author_cell_type.isin(prot)
rng = np.random.default_rng(0); pick3, pickb = [], []
for (_, _t), gg in M[M.keep & M.is_c3].groupby(['sample', 'author_cell_type'], observed=True): pick3.append(rng.choice(gg.row.values, min(100, len(gg)), replace=False))
for (_, _t), gg in M[M.keep & M.is_br].groupby(['bsample', 'author_cell_type'], observed=True): pickb.append(rng.choice(gg.row.values, min(30, len(gg)), replace=False))
pick3 = np.concatenate(pick3); pickb = np.concatenate(pickb); pick = np.union1d(pick3, pickb)
M['sub3'] = M.row.isin(pick3); M['subb'] = M.row.isin(pickb)
c3k = M[M.keep & M.is_c3]; key = c3k['sample'] + '#' + c3k.author_cell_type.astype(str); codes, uk = pd.factorize(key)
cmap = pd.Series(codes, index=c3k.row.values); PBm = np.zeros((len(uk), len(keepg))); subX, subrows = [], []
t0 = time.time()
for r, Xr in blocks():
    Xr = Xr[:, keepg]; k = np.isin(r, c3k.row.values)
    if k.any():
        cc = cmap.loc[r[k]].values; A = sp.csr_matrix((np.ones(len(cc)), (cc, np.arange(len(cc)))), shape=(len(uk), len(cc))); PBm += (A @ Xr[k]).toarray()
    s_ = np.isin(r, pick)
    if s_.any(): subX.append(Xr[s_]); subrows.append(r[s_])
M.drop(columns=['lm']).to_pickle(f'{W}/perez_cells.pkl')
pd.DataFrame(PBm.astype(np.int64), index=uk, columns=sym[keepg]).to_pickle(f'{W}/perez_pb.pkl')
sr = np.concatenate(subrows); assert (np.diff(sr) > 0).all()
sp.save_npz(f'{W}/perez_sub_X.npz', sp.vstack(subX).tocsr()); M.set_index('row').loc[sr].to_pickle(f'{W}/perez_sub_meta.pkl')
pd.Series(sym[keepg]).to_csv(f'{W}/genes.csv', index=False)
print('done; C3 sub', len(pick3), 'bridging sub', len(pickb), flush=True)
