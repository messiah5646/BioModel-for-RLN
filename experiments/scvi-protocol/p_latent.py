"""Round 102 (Q36): per-fold PCA-50 (B4) and scVI (D1, batch_key library_uuid) latents for the subsampled cells.
Usage: python3 -I p_latent.py <work_dir> <C3|C4> [fold ...]. Folds: C3 = held-out donor; C4 = 5 stratified donor folds, seed 0."""
import sys, os, numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad, torch, scvi
from sklearn.model_selection import StratifiedKFold
W, C = sys.argv[1:3]; torch.set_num_threads(2); scvi.settings.seed = 0
M = pd.read_pickle(f'{W}/sub_meta.pkl'); X = sp.load_npz(f'{W}/sub_X.npz').tocsr()
genes = pd.read_csv(f'{W}/genes.csv').iloc[:, 0].astype(str).values
k = (M.contrast == C).values; M = M[k]; X = X[k]
A = ad.AnnData(X.astype(np.float32), obs=M.reset_index(drop=True)[['sample', 'donor_id', 'library_uuid', 'author_cell_type', 'y']].astype({'sample': str, 'donor_id': str, 'library_uuid': str}))
A.var_names = pd.Index(genes); A.var_names_make_unique(); A.obs_names = [f'c{i}' for i in range(A.n_obs)]
def fold_of():
    if C == 'C3': return {d: d for d in A.obs.donor_id.unique()}
    sm = A.obs.groupby('sample').y.first(); out = {}
    for i, (_, te) in enumerate(StratifiedKFold(5, shuffle=True, random_state=0).split(sm.index, sm.values)):
        for s in sm.index[te]: out[s] = f'f{i}'
    return out
fo = fold_of(); key = A.obs.donor_id if C == 'C3' else A.obs['sample']; A.obs['fold'] = key.map(fo).values
folds = sys.argv[3:] or sorted(A.obs.fold.unique()); os.makedirs(f'{W}/latent_{C}', exist_ok=True)
pd.Series(fo).to_csv(f'{W}/latent_{C}/folds.csv')
for d in folds:
    if os.path.exists(f'{W}/latent_{C}/{d}.npz'): continue
    tr = A[A.obs.fold != d].copy(); te = A[A.obs.fold == d].copy()
    try: sc.pp.highly_variable_genes(tr, flavor='seurat_v3', n_top_genes=2000, batch_key='library_uuid')
    except ValueError:  # deviation logged in PREREG-C3C4-round102.md: loess singular in small libraries -> span 1.0
        sc.pp.highly_variable_genes(tr, flavor='seurat_v3', n_top_genes=2000, batch_key='library_uuid', span=1.0); print('HVG span 1.0 fallback', d, flush=True)
    hv = tr.var_names[tr.var.highly_variable]
    def ln(a):
        b = a[:, hv].copy(); sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b); return b
    trl, tel = ln(tr), ln(te)
    mu = np.asarray(trl.X.mean(0)).ravel(); sd = np.sqrt(np.asarray(trl.X.power(2).mean(0)).ravel() - mu ** 2) + 1e-6
    Z = (trl.X.toarray() - mu) / sd; U, S, Vt = np.linalg.svd(Z, full_matrices=False); V = Vt[:50].T
    pca_tr = Z @ V; pca_te = ((tel.X.toarray() - mu) / sd) @ V; del Z, U
    trs = tr[:, hv].copy(); tes = te[:, hv].copy()
    scvi.model.SCVI.setup_anndata(trs, batch_key='library_uuid')
    m = scvi.model.SCVI(trs, n_latent=20, n_layers=2, gene_likelihood='nb', use_layer_norm='both', use_batch_norm='none',
                        encode_covariates=True, dropout_rate=0.2)
    m.train(max_epochs=200, early_stopping=True, accelerator='cpu', enable_progress_bar=False)
    scvi.model.SCVI.prepare_query_anndata(tes, m); q = scvi.model.SCVI.load_query_data(tes, m)
    q.train(max_epochs=100, plan_kwargs=dict(weight_decay=0.0), accelerator='cpu', enable_progress_bar=False)
    np.savez(f'{W}/latent_{C}/{d}.npz', tr_idx=np.flatnonzero((A.obs.fold != d).values), te_idx=np.flatnonzero((A.obs.fold == d).values),
             pca_tr=pca_tr, pca_te=pca_te, scvi_tr=m.get_latent_representation(), scvi_te=q.get_latent_representation(),
             epochs=len(m.history['elbo_train']))
    print(C, 'fold', d, 'done; epochs', len(m.history['elbo_train']), flush=True)
