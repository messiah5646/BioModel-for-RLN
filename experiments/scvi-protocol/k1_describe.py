"""Round 103 K1 descriptive report (PREREG-K1-round103.md). Usage: python -I k1_describe.py <workdir>"""
import sys, json, numpy as np, pandas as pd, anndata as ad, scanpy as sc
wd = sys.argv[1]
Q = ad.read_h5ad(f'{wd}/k1_query_all.h5ad')
o = pd.read_csv(f'{wd}/k1_query_obs.csv', index_col=0).reindex(Q.obs_names)
for c in ['pred_scanvi', 'pred_knn_pca', 'maxprob', 'knn_unc']: Q.obs[c] = o[c].values
NONIMM = ['PT', 'TAL', 'PC', 'IC', 'EC', 'DTL/ATL', 'DCT/CNT', 'VSMC/MC/FIB', 'POD/PEC']
MYE = ['MON', 'ncMON', 'MAC-M2', 'MDC', 'cDC', 'cycMNP', 'pDC']
# Stage 1 removal with the protected-population rule (per sample: a label whose flag rate > 3x the sample rate keeps its flagged cells)
fl = (Q.obs.qc_fail | Q.obs.doublet).values
keep = ~fl; prot = []
for s, g in Q.obs.groupby('sample'):
    r = (g.qc_fail | g.doublet).mean()
    for l, gl in g.groupby('pred_scanvi'):
        rl = (gl.qc_fail | gl.doublet).mean()
        if len(gl) >= 20 and r > 0 and rl > 3 * r:
            keep[Q.obs_names.get_indexer(gl.index)] = True; prot.append((s, l, len(gl), round(rl, 3), round(r, 3)))
print('protected (sample, label, n, label flag rate, sample flag rate):', prot)
Q = Q[keep].copy()
L = Q.copy(); sc.pp.normalize_total(L, target_sum=1e4); sc.pp.log1p(L)
L.var_names = pd.Index(np.asarray(L.var.symbol, dtype=str)); L.var_names_make_unique()
out = {'protected': prot}
ob = Q.obs
ob['unresolved'] = ob.knn_unc > 0.5
ob['nonimmune'] = ob.pred_scanvi.isin(NONIMM)
summ = ob.groupby('sample').agg(cells=('pred_scanvi', 'size'), nonimmune_frac=('nonimmune', 'mean'), unresolved_frac=('unresolved', 'mean'),
                                median_maxprob=('maxprob', 'median'), agree_knn_pca=('pred_knn_pca', lambda x: np.nan))
summ['agree_knn_pca'] = ob.assign(a=ob.pred_scanvi == ob.pred_knn_pca).groupby('sample').a.mean()
print(summ.round(3).to_string()); out['summary'] = summ.round(4).reset_index().to_dict('records')
im = ob[~ob.nonimmune]
prop = pd.crosstab(im.pred_scanvi, im['sample'], normalize='columns').round(3)
print(prop.to_string()); out['immune_props'] = prop.to_dict()
unr = im.groupby(['sample', 'pred_scanvi']).unresolved.mean().unstack(0).round(2)
print('unresolved fraction by label x sample\n', unr.to_string()); out['unresolved_by_label'] = unr.to_dict()
# marker cross-check
MK = {'T': 'CD3E', 'T-CYT': 'CD8A', 'NKT': 'NKG7', 'NK1': 'NKG7', 'B': 'MS4A1', 'PL': 'MZB1', 'MON': 'CD14', 'ncMON': 'FCGR3A',
      'MAC-M2': 'C1QA', 'cDC': 'CD1C', 'MDC': 'C1QA', 'pDC': 'LILRA4', 'MAST': 'TPSAB1', 'cycT': 'MKI67', 'cycMNP': 'MKI67'}
genes = sorted(set(MK.values()) | {'PTPRC', 'CLEC9A', 'FCER1A', 'CD163', 'MRC1', 'LYVE1', 'APOE', 'SPP1', 'FCN1', 'S100A8', 'IGJ', 'JCHAIN', 'CSF3R', 'FCGR3B', 'CXCR2', 'VCAN', 'CX3CR1'})
genes = [g for g in genes if g in L.var_names]
mm = pd.DataFrame(L[:, genes].X.toarray(), columns=genes, index=L.obs_names).groupby(L.obs.pred_scanvi.values).mean()
print(mm.round(2).to_string())
flags = {l: (mm.loc[l, MK[l]] < mm[MK[l]].quantile(0.75)) for l in mm.index if l in MK and MK[l] in mm.columns}
out['marker_flags_own_marker_below_q75'] = {k: bool(v) for k, v in flags.items()}
out['marker_means'] = mm.round(3).to_dict()
# programme scores (descriptive)
R91 = ['CD99', 'MAP4K4', 'CSF1R', 'PTPN6', 'TGFB1', 'HLA-DRA', 'IFI16', 'LILRB1', 'LILRB2', 'C1S']
IFN6 = ['IFI27', 'IFI44L', 'IFIT1', 'ISG15', 'RSAD2', 'SIGLEC1']
sc.tl.score_genes(L, [g for g in R91 if g in L.var_names], score_name='r91', random_state=0)
sc.tl.score_genes(L, [g for g in IFN6 if g in L.var_names], score_name='ifn6', random_state=0)
L.obs['group'] = L.obs['sample'].map(lambda s: 'LN' if s == 'LN' else ('NC' if s == 'NC' else 'AGN'))
L.obs['comp'] = np.where(L.obs.pred_scanvi.isin(MYE), 'myeloid', np.where(L.obs.pred_scanvi.isin(['T', 'T-CYT', 'NKT', 'NK1', 'cycT']), 'T/NK', 'other'))
sco = L.obs[L.obs.comp != 'other'].groupby(['comp', 'sample'])[['r91', 'ifn6']].median().round(3)
print(sco.to_string()); out['scores_median_by_comp_sample'] = {f'{a}|{b}': v for (a, b), v in sco.to_dict('index').items()}
# macrophage state split inside LN vs NC myeloid (resident-like vs infiltrating-like), descriptive
my = L[L.obs.comp == 'myeloid']
ms = pd.DataFrame(my[:, [g for g in ['CD14', 'FCN1', 'S100A8', 'C1QA', 'CD163', 'MRC1', 'LYVE1', 'APOE', 'SPP1'] if g in my.var_names]].X.toarray(),
                  columns=[g for g in ['CD14', 'FCN1', 'S100A8', 'C1QA', 'CD163', 'MRC1', 'LYVE1', 'APOE', 'SPP1'] if g in my.var_names], index=my.obs_names)
ms = ms.groupby(my.obs['sample'].values).mean().round(2); print(ms.to_string()); out['myeloid_markers_by_sample'] = ms.to_dict()
json.dump(out, open(f'{wd}/k1_describe.json', 'w'), indent=1, default=str)
