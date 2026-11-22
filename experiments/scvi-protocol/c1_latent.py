"""Round 100 C1: per-fold PCA-50 (B4) and scVI (D1) latents for the subsampled cells. Usage: python3 -I c1_latent.py <work_dir> [fold ...]"""
import sys, os, numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad, torch, scvi
W = sys.argv[1]; torch.set_num_threads(2); scvi.settings.seed = 0; sc.settings.n_jobs = 2
M = pd.read_csv(f'{W}/sub_meta.csv.gz', index_col=0); X = sp.load_npz(f'{W}/sub_X.npz').tocsr()
genes = pd.read_csv(f'{W}/genes.csv').iloc[:, 0].astype(str).values
A = ad.AnnData(X.astype(np.float32), obs=M.reset_index(drop=True).astype({'donor_id': str, 'sample': str}))
A.var_names = pd.Index(genes).astype(str); A.var_names_make_unique(); A.obs_names = [f'c{i}' for i in range(A.n_obs)]
donors = sorted(A.obs.donor_id.unique()); folds = sys.argv[2:] or donors
os.makedirs(f'{W}/latent', exist_ok=True)
for d in folds:
    if os.path.exists(f'{W}/latent/{d}.npz'): continue
    tr = A[A.obs.donor_id != d].copy(); te = A[A.obs.donor_id == d].copy()
    sc.pp.highly_variable_genes(tr, flavor='seurat_v3', n_top_genes=2000, batch_key='donor_id')
    hv = tr.var_names[tr.var.highly_variable]
    # B4: PCA-50 on log-normalised HVGs, fit on training cells
    def ln(a):
        b = a[:, hv].copy(); sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b); return b
    trl, tel = ln(tr), ln(te)
    mu = np.asarray(trl.X.mean(0)).ravel(); sd = np.sqrt(np.asarray(trl.X.power(2).mean(0)).ravel() - mu ** 2) + 1e-6
    Z = (trl.X.toarray() - mu) / sd; U, S, Vt = np.linalg.svd(Z, full_matrices=False); V = Vt[:50].T
    pca_tr = Z @ V; pca_te = ((tel.X.toarray() - mu) / sd) @ V
    # D1: scVI reference on training donors, scArches query for held-out donor
    trs = tr[:, hv].copy(); tes = te[:, hv].copy()
    scvi.model.SCVI.setup_anndata(trs, batch_key='donor_id')
    m = scvi.model.SCVI(trs, n_latent=20, n_layers=2, gene_likelihood='nb', use_layer_norm='both', use_batch_norm='none',
                        encode_covariates=True, dropout_rate=0.2)
    m.train(max_epochs=200, early_stopping=True, accelerator='cpu', enable_progress_bar=False)
    scvi.model.SCVI.prepare_query_anndata(tes, m)
    q = scvi.model.SCVI.load_query_data(tes, m)
    q.train(max_epochs=100, plan_kwargs=dict(weight_decay=0.0), accelerator='cpu', enable_progress_bar=False)
    np.savez(f'{W}/latent/{d}.npz', tr_idx=np.flatnonzero((A.obs.donor_id != d).values), te_idx=np.flatnonzero((A.obs.donor_id == d).values),
             pca_tr=pca_tr, pca_te=pca_te, scvi_tr=m.get_latent_representation(), scvi_te=q.get_latent_representation(),
             epochs=len(m.history['elbo_train']))
    print('fold', d, 'done; epochs', len(m.history['elbo_train']), flush=True)
