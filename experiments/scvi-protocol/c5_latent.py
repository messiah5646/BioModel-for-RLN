"""Round 104 (Q37/C5) latents. Step 'c5': HVGs + PCA-50 + scVI on Perez C3 cells; scArches map of GSE254176 (fine-tune on a label-blind
<= 1,500 cells/sample subsample, then encode all cells); kNN (k 15) cell-type transfer; GSE subsample; sysVI (transductive, label-free)
on Perez C3 + GSE subsample. Step 'bridge': PCA-50, scVI, sysVI on bridging-donor cells.
Usage: python3 -I c5_latent.py <work_dir> c5|bridge"""
import sys, os, numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad, torch, scvi
from sklearn.neighbors import KNeighborsClassifier
W, STEP = sys.argv[1:3]; torch.set_num_threads(2); scvi.settings.seed = 0
genes = pd.read_csv(f'{W}/genes.csv').iloc[:, 0].astype(str).values
Mp = pd.read_pickle(f'{W}/perez_sub_meta.pkl'); Xp = sp.load_npz(f'{W}/perez_sub_X.npz').tocsr()
def mk(X, obs):
    A = ad.AnnData(X.astype(np.float32), obs=obs.reset_index(drop=True)); A.var_names = pd.Index(genes); A.obs_names = [f'c{i}' for i in range(A.n_obs)]; return A
def hvg(A, key):
    try: sc.pp.highly_variable_genes(A, flavor='seurat_v3', n_top_genes=2000, batch_key=key)
    except ValueError: sc.pp.highly_variable_genes(A, flavor='seurat_v3', n_top_genes=2000, batch_key=key, span=1.0); print('HVG span 1.0 fallback', flush=True)
    return A.var_names[A.var.highly_variable]
def lognorm(X):
    X = X.tocsr().astype(np.float32); t = np.asarray(X.sum(1)).ravel(); L = sp.diags(1e4 / np.maximum(t, 1)) @ X; L = L.tocsr(); L.data = np.log1p(L.data); return L
def pca_fit(L):
    mu = np.asarray(L.mean(0)).ravel(); sd = np.sqrt(np.asarray(L.power(2).mean(0)).ravel() - mu ** 2) + 1e-6
    Z = (L.toarray() - mu) / sd; V = np.linalg.svd(Z, full_matrices=False)[2][:50].T; return mu, sd, V, Z @ V
def scvi_fit(A, key):
    scvi.model.SCVI.setup_anndata(A, batch_key=key)
    m = scvi.model.SCVI(A, n_latent=20, n_layers=2, gene_likelihood='nb', use_layer_norm='both', use_batch_norm='none', encode_covariates=True, dropout_rate=0.2)
    m.train(max_epochs=200, early_stopping=True, accelerator='cpu', enable_progress_bar=False); return m
def sysvi_fit(A, system, batch):
    from scvi.external import SysVI
    A = A.copy(); A.X = lognorm(A.X)  # SysVI models log-normalised expression (Gaussian), per scvi-tools docs
    SysVI.setup_anndata(A, batch_key=system, categorical_covariate_keys=[batch])
    m = SysVI(A); m.train(max_epochs=200, early_stopping=True, accelerator='cpu', enable_progress_bar=False); return m.get_latent_representation(A)
if STEP == 'c5':
    k3 = Mp.sub3.values; P = mk(Xp[k3], Mp[k3][['sample', 'donor_id', 'library_uuid', 'author_cell_type', 'y']].astype(str))
    hv = hvg(P, 'library_uuid'); pd.Series(hv).to_csv(f'{W}/c5_hvg.csv', index=False)
    mu, sd, V, pca_p = pca_fit(lognorm(P[:, hv].X)); Ps = P[:, hv].copy(); m = scvi_fit(Ps, 'library_uuid'); zp = m.get_latent_representation()
    print('scVI perez epochs', len(m.history['elbo_train']), flush=True)
    C = pd.read_pickle(f'{W}/gse_cells.pkl'); hvi = pd.Index(genes).get_indexer(hv)
    Xg = sp.vstack([sp.load_npz(f'{W}/gse/{g}.npz').tocsr()[:, hvi] for g in C.gsm.unique()]).tocsr()
    G = ad.AnnData(Xg.astype(np.float32), obs=pd.DataFrame({'library_uuid': C.gsm.values})); G.var_names = pd.Index(hv); G.obs_names = [f'g{i}' for i in range(G.n_obs)]
    rng = np.random.default_rng(0); ft = np.sort(np.concatenate([rng.choice(np.flatnonzero((C.gsm == g).values), min(1500, int((C.gsm == g).sum())), replace=False) for g in C.gsm.unique()]))
    Gf = G[ft].copy(); scvi.model.SCVI.prepare_query_anndata(Gf, m); q = scvi.model.SCVI.load_query_data(Gf, m)
    q.train(max_epochs=100, plan_kwargs=dict(weight_decay=0.0), accelerator='cpu', enable_progress_bar=False)
    scvi.model.SCVI.prepare_query_anndata(G, q); zg = q.get_latent_representation(G)
    knn = KNeighborsClassifier(15).fit(zp, P.obs.author_cell_type.values); C['ctype'] = knn.predict(zg); C['ctype_maxp'] = knn.predict_proba(zg).max(1)
    pca_g = ((lognorm(Xg).toarray() - mu) / sd) @ V
    # protected-population rule with transferred labels, then label-blind GSE subsample (<= 100 / sample x type)
    base = C.qc_fail.mean(); rt = C.groupby('ctype').qc_fail.mean(); prot = rt[rt > 3 * base].index.tolist(); print('GSE flag', round(base, 4), 'by type', rt.round(3).to_dict(), 'protected', prot, flush=True)
    C['keep'] = ~C.qc_fail | C.ctype.isin(prot)
    C['sub'] = False
    for _, g in C[C.keep].groupby(['gsm', 'ctype']): C.loc[rng.choice(g.index.values, min(100, len(g)), replace=False), 'sub'] = True
    C.to_pickle(f'{W}/gse_cells_annot.pkl')
    np.savez(f'{W}/c5_latent.npz', pca_p=pca_p, scvi_p=zp, pca_g=pca_g, scvi_g=zg, epochs=len(m.history['elbo_train']))
if STEP in ('c5', 'c5sysvi'):  # 'c5sysvi' resumes here from saved files (bug fix: C.sub -> C['sub'])
    k3 = Mp.sub3.values; P = mk(Xp[k3], Mp[k3][['sample', 'donor_id', 'library_uuid', 'author_cell_type', 'y']].astype(str))
    hv = pd.read_csv(f'{W}/c5_hvg.csv').iloc[:, 0].astype(str).values; Ps = P[:, hv].copy(); C = pd.read_pickle(f'{W}/gse_cells_annot.pkl')
    hvi = pd.Index(genes).get_indexer(hv); Xg = sp.vstack([sp.load_npz(f'{W}/gse/{g}.npz').tocsr()[:, hvi] for g in C.gsm.unique()]).tocsr()
    # sysVI: transductive, label-free; Perez C3 sub + GSE sub, system = dataset, covariate = library/GSM
    gs = C['sub'].values; S = ad.concat([Ps, ad.AnnData(Xg[gs].astype(np.float32), obs=pd.DataFrame({'library_uuid': C.gsm.values[gs]}), var=pd.DataFrame(index=hv))], join='inner')
    S.obs['system'] = ['perez'] * Ps.n_obs + ['gse'] * int(gs.sum()); S.obs['batch'] = S.obs.library_uuid.astype(str); S.obs_names = [f's{i}' for i in range(S.n_obs)]
    zs = sysvi_fit(S, 'system', 'batch'); np.savez(f'{W}/c5_sysvi.npz', p=zs[:Ps.n_obs], g=zs[Ps.n_obs:], gidx=np.flatnonzero(gs)); print('c5 done', flush=True)
if STEP == 'bridge':
    kb = Mp.subb.values; B = mk(Xp[kb], Mp[kb][['bsample', 'donor_id', 'library_uuid', 'Processing_Cohort', 'author_cell_type']].astype(str))
    hv = hvg(B, 'library_uuid'); Bs = B[:, hv].copy(); pca_b = pca_fit(lognorm(Bs.X))[3]
    m = scvi_fit(Bs.copy(), 'library_uuid'); zb = m.get_latent_representation()
    S = Bs.copy(); zs = sysvi_fit(S, 'Processing_Cohort', 'library_uuid')
    np.savez(f'{W}/bridge_latent.npz', pca=pca_b, scvi=zb, sysvi=zs); print('bridge done', flush=True)
