"""Round 103 K1 stage 1: KPMP reference subset + E-MTAB-13596 query QC. Follows PREREG-K1-round103.md.
Usage: python -I k1_prep.py <kpmp.h5ad> <emtab_dir> <outdir>"""
import sys, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio, anndata as ad, scanpy as sc
kp, qd, out = sys.argv[1:4]
rng = np.random.default_rng(0)
sc.settings.n_jobs = 2  # scrublet fails with n_jobs=-1 on this 2-CPU box

# ---------- reference ----------
import os
SKIP_REF = os.path.exists(f'{out}/k1_ref.h5ad')
if not SKIP_REF:
    a = ad.read_h5ad(kp, backed='r')
    o = a.obs[['author_cell_type', 'subclass.l1', 'donor_id', 'LibraryID', 'disease']].copy()
    imm = o['subclass.l1'].isin(['T', 'B', 'MYL']).values
    lab = np.where(imm, o['author_cell_type'].astype(str), o['subclass.l1'].astype(str))
    lab = pd.Series(lab, index=o.index)
    vc = lab[imm].value_counts()
    small = vc[vc < 50].index
    lab[imm & lab.isin(small).values] = o['subclass.l1'].astype(str)[imm & lab.isin(small).values]
    o['label'] = lab.values; o['immune'] = imm
    keep = []
    pos = pd.Series(np.arange(len(o)), index=o.index)
    oi = o[imm]
    for _, g in oi.groupby(['donor_id', 'label'], observed=True):
        p = pos[g.index].values
        keep.extend(rng.choice(p, min(300, len(p)), replace=False))
    on = o[~imm]
    for _, g in on.groupby('label', observed=True):
        p = pos[g.index].values
        keep.extend(rng.choice(p, min(1000, len(p)), replace=False))
    keep = np.sort(np.array(keep))
    print('reference cells kept', len(keep), 'immune', int(imm[keep].sum()), 'merged small labels', list(small))
    X = a.raw.X[keep] if a.raw is not None else a.X[keep]
    X = sp.csr_matrix(X)
    ref = ad.AnnData(X=X, obs=o.iloc[keep].copy(), var=pd.DataFrame(index=a.raw.var_names if a.raw is not None else a.var_names))
    ref.var['symbol'] = a.var['feature_name'].astype(str).reindex(ref.var_names).values
    ref.obs['batch'] = ref.obs['LibraryID'].astype(str)
    ref.obs['dataset'] = 'KPMP'
    assert np.allclose(ref.X[:50].data, np.round(ref.X[:50].data)), 'reference not raw counts'
    ref.write_h5ad(f'{out}/k1_ref.h5ad')
    print(ref.obs.label.value_counts().to_string())

# ---------- query QC ----------
ST = ['FOS', 'FOSB', 'JUN', 'JUNB', 'HSPA1A', 'HSPA1B', 'DUSP1', 'IER2']
rows, qs = [], []
for s in ['AGN1', 'AGN2', 'AGN3', 'AGN4', 'AGN5', 'LN', 'NC']:
    M = sio.mmread(f'{qd}/{s}_matrix.mtx').T.tocsr()
    f = pd.read_csv(f'{qd}/{s}_features.tsv', sep='\t', header=None)
    b = pd.read_csv(f'{qd}/{s}_barcodes.tsv', sep='\t', header=None)[0]
    g = f[2] == 'Gene Expression'
    M = M[:, np.flatnonzero(g.values)]; f = f[g.values]
    q = ad.AnnData(M.astype(np.float32), obs=pd.DataFrame(index=s + '_' + b.values), var=pd.DataFrame({'symbol': f[1].values}, index=f[0].values))
    q.var_names_make_unique()
    sym = np.asarray(q.var.symbol.values, dtype=str)
    tot = np.asarray(M.sum(1)).ravel(); ng = np.asarray((M > 0).sum(1)).ravel()
    mt = np.char.startswith(sym, 'MT-')
    mito = np.asarray(M[:, mt].sum(1)).ravel() / np.maximum(tot, 1)
    top20 = np.array([np.sort(M.data[M.indptr[i]:M.indptr[i + 1]])[-20:].sum() for i in range(M.shape[0])]) / np.maximum(tot, 1)
    def mad_out(x, n, side='both'):
        med = np.median(x); mad = np.median(np.abs(x - med)) * 1.4826 + 1e-9
        return (x > med + n * mad) | ((x < med - n * mad) if side == 'both' else False)
    lm = np.log((mito + 1e-4) / (1 - mito + 1e-4))
    fail = mad_out(np.log1p(tot), 5) | mad_out(np.log1p(ng), 5) | mad_out(top20, 5) | (mad_out(lm, 3, 'upper') & (mito > 0.5))
    try:
        t = q.copy(); sc.pp.scrublet(t, random_state=0); dbl = t.obs.predicted_doublet.fillna(False).values.astype(bool); ds = t.obs.doublet_score.values
    except Exception as e:
        print(s, 'scrublet failed', e); dbl = np.zeros(len(q), bool); ds = np.full(len(q), np.nan)
    L = q.copy(); sc.pp.normalize_total(L, target_sum=1e4); sc.pp.log1p(L)
    sti = np.isin(sym, ST)
    q.obs['sample'] = s; q.obs['log_counts'] = np.log1p(tot); q.obs['log_genes'] = np.log1p(ng); q.obs['mito'] = mito
    q.obs['top20'] = top20; q.obs['doublet_score'] = ds; q.obs['doublet'] = dbl; q.obs['qc_fail'] = fail
    q.obs['stress'] = np.asarray(L.X[:, sti].mean(1)).ravel()
    rows.append(dict(sample=s, cells=len(q), qc_fail=int(fail.sum()), doublet=int(dbl.sum()), removed=int((fail | dbl).sum()),
                     median_log_counts=float(np.median(np.log1p(tot))), median_mito=float(np.median(mito))))
    print(rows[-1], flush=True)
    qs.append(q)
Q = ad.concat(qs, join='outer', merge='first')
Q.var['symbol'] = pd.concat([x.var.symbol for x in qs]).groupby(level=0).first().reindex(Q.var_names).values
Q.obs['batch'] = Q.obs['sample']; Q.obs['dataset'] = 'E-MTAB-13596'
Q.write_h5ad(f'{out}/k1_query_all.h5ad')
pd.DataFrame(rows).to_csv(f'{out}/k1_query_qc.csv', index=False)
