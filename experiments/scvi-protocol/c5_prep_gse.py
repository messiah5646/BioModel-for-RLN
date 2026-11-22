"""Round 104 (Q37/C5) GSE254176 side: per-sample Stage-1 QC (PREREG-C5-round104.md), genes aligned to the Perez/GSE shared list.
Writes per-sample counts (aligned genes) and a per-cell metrics table. No labels used here.
Usage: python3 -I c5_prep_gse.py <raw_dir> <labels.tsv> <work_dir>"""
import sys, os, glob, numpy as np, pandas as pd, scipy.sparse as sp, scipy.io as sio
R, LAB, W = sys.argv[1:4]; os.makedirs(f'{W}/gse', exist_ok=True)
genes = pd.read_csv(f'{W}/genes.csv').iloc[:, 0].astype(str).values; gpos = {g: i for i, g in enumerate(genes)}
lab = pd.read_csv(LAB, sep='\t', header=None, names=['gsm', 'title', 'state'])
s = pd.Series(genes)
def idx(mask): return np.flatnonzero(mask)
mt = idx(s.str.startswith('MT-').values); ribo = idx(s.str.match(r'^RP[SL]\d').values)
HB = idx(s.isin(['HBA1', 'HBA2', 'HBB']).values); PL = idx(s.isin(['PPBP', 'PF4']).values)
IG = idx(s.str.match(r'^(IGHG\d|IGKC|IGLC\d|JCHAIN|IGHA\d|IGHM)$').values)
ST = idx(s.isin(['FOS', 'FOSB', 'JUN', 'JUNB', 'HSPA1A', 'HSPA1B', 'DUSP1', 'IER2']).values)
IFN = idx(s.isin(['IFI27', 'IFI44L', 'IFIT1', 'ISG15', 'RSAD2', 'SIGLEC1']).values); assert len(IFN) == 6
def mad_out(x, n, side='both'):
    med = np.median(x); mad = np.median(np.abs(x - med)) * 1.4826 + 1e-9
    return (x > med + n * mad) | ((x < med - n * mad) if side == 'both' else False)
out = []
for gsm in lab.gsm:
    pre = glob.glob(f'{R}/{gsm}_*matrix.mtx.gz')[0][:-len('matrix.mtx.gz')]
    Xg = sio.mmread(pre + 'matrix.mtx.gz').T.tocsr().astype(np.float32)  # cells x genes (full 36,601, incl. symbols)
    fsym = pd.read_csv(pre + 'features.tsv.gz', sep='\t', header=None)[1].astype(str).values
    tot_all = np.asarray(Xg.sum(1)).ravel(); ng_all = np.diff(Xg.indptr)
    first = pd.Series(fsym).drop_duplicates().index.values; colmap = [(j, gpos[fsym[j]]) for j in first if fsym[j] in gpos]
    src, dst = map(np.array, zip(*colmap))
    T = sp.csr_matrix((np.ones(len(src), np.float32), (src, dst)), shape=(Xg.shape[1], len(genes)))
    Xa = (Xg @ T).tocsr()
    k200 = ng_all >= 200; Xa = Xa[k200]; tot = tot_all[k200]; ng = ng_all[k200]
    full = Xg[k200]; top20 = np.array([np.sort(full.data[full.indptr[i]:full.indptr[i + 1]])[-20:].sum() for i in range(full.shape[0])]) / tot
    def frac(ix): return np.asarray(Xa[:, ix].sum(1)).ravel() / tot
    L = Xa.multiply(1e4 / tot[:, None]).tocsr(); L.data = np.log1p(L.data)
    m = pd.DataFrame(dict(gsm=gsm, cell=np.arange(Xa.shape[0]), log_counts=np.log1p(tot), log_genes=np.log1p(ng), top20=top20,
                          mito=frac(mt) * 100, ribo=frac(ribo) * 100, hb=frac(HB), plt=frac(PL), ig=frac(IG),
                          stress=np.asarray(L[:, ST].mean(1)).ravel(), ifn6=np.asarray(L[:, IFN].mean(1)).ravel()))
    lm = np.log((m.mito / 100 + 1e-4) / (1 - m.mito / 100 + 1e-4)).values
    m['qc_fail'] = mad_out(m.log_counts.values, 5) | mad_out(m.log_genes.values, 5) | mad_out(m.top20.values, 5) | (mad_out(lm, 3, 'upper') & (m.mito.values > 15))
    sp.save_npz(f'{W}/gse/{gsm}.npz', Xa); out.append(m)
    print(gsm, Xg.shape[0], 'barcodes; >=200 genes', int(k200.sum()), '; flagged', int(m.qc_fail.sum()), flush=True)
pd.concat(out, ignore_index=True).to_pickle(f'{W}/gse_cells.pkl')
