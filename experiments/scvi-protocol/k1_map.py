"""Round 103 K1: held-out sanity gate (scANVI vs kNN-PCA vs majority) and KPMP -> E-MTAB-13596 mapping.
Follows PREREG-K1-round103.md. Usage: python -I k1_map.py <workdir> <stage: val|full>"""
import sys, json, numpy as np, pandas as pd, anndata as ad, scanpy as sc, scvi, torch
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
from sklearn.metrics import accuracy_score, f1_score
torch.set_num_threads(2); scvi.settings.seed = 0; sc.settings.n_jobs = 2
wd, stage = sys.argv[1:3]
IMM_L1 = ['T', 'B', 'MYL']
ref = ad.read_h5ad(f'{wd}/k1_ref.h5ad')
vc = ref.obs.label.value_counts()
ref.obs['label_train'] = ref.obs.label.astype(str).where(ref.obs.label.map(vc) >= 10, 'Unknown')
print('labels dropped (<10 cells):', list(vc[vc < 10].index))

SCVI = dict(n_latent=20, n_layers=2, dropout_rate=0.2, gene_likelihood='nb', encode_covariates=True,
            deeply_inject_covariates=False, use_layer_norm='both', use_batch_norm='none')

def hvg(a):
    t = a.copy(); sc.pp.highly_variable_genes(t, n_top_genes=2000, flavor='seurat_v3', batch_key='batch', span=1.0)
    return t.var_names[t.var.highly_variable].values

def fit_ref(a, genes):
    a = a[:, genes].copy()
    scvi.model.SCVI.setup_anndata(a, batch_key='batch')
    m = scvi.model.SCVI(a, **SCVI)
    m.train(max_epochs=200, early_stopping=True, accelerator='cpu', batch_size=512)
    print('scVI epochs', len(m.history['elbo_train']), flush=True)
    s = scvi.model.SCANVI.from_scvi_model(m, unlabeled_category='Unknown', labels_key='label_train', adata=a)
    s.train(max_epochs=20, n_samples_per_label=100, accelerator='cpu', batch_size=512)
    return a, s

def map_query(s, q, genes):
    q = q.copy()
    miss = [g for g in genes if g not in q.var_names]
    import scipy.sparse as sp
    have = [g for g in genes if g in q.var_names]
    X = sp.hstack([sp.csr_matrix(q[:, have].X), sp.csr_matrix((q.n_obs, len(miss)), dtype=np.float32)]).tocsr()
    q = ad.AnnData(X, obs=q.obs.copy(), var=pd.DataFrame(index=have + miss))[:, list(genes)].copy()
    q.obs['label_train'] = 'Unknown'
    scvi.model.SCANVI.prepare_query_anndata(q, s)
    sq = scvi.model.SCANVI.load_query_data(q, s)
    sq.train(max_epochs=100, plan_kwargs=dict(weight_decay=0.0), accelerator='cpu', batch_size=512)
    return q, sq, len(miss)

def lognorm(a):
    t = a.copy(); sc.pp.normalize_total(t, target_sum=1e4); sc.pp.log1p(t)
    X = t.X; return X.toarray() if hasattr(X, 'toarray') else np.asarray(X)

if stage == 'val':
    donors = np.array(sorted(ref.obs.donor_id.astype(str).unique()))
    rng = np.random.default_rng(0); test_d = rng.choice(donors, int(round(0.2 * len(donors))), replace=False)
    te = ref.obs.donor_id.astype(str).isin(test_d).values
    tr_a, te_a = ref[~te].copy(), ref[te].copy()
    genes = hvg(tr_a)
    a, s = fit_ref(tr_a, genes)
    q, sq, nm = map_query(s, te_a, genes)
    pred_s = sq.predict(q)
    # kNN-PCA baseline
    Xtr, Xte = lognorm(tr_a[:, genes]), lognorm(te_a[:, genes])
    p = PCA(50, random_state=0).fit(Xtr)
    knn = KNeighborsClassifier(15).fit(p.transform(Xtr), tr_a.obs.label_train.values)
    pred_k = knn.predict(p.transform(Xte))
    maj = tr_a.obs.label_train[tr_a.obs['subclass.l1'].isin(IMM_L1)].value_counts().index[0]
    y = te_a.obs.label.astype(str).values
    imm = te_a.obs['subclass.l1'].isin(IMM_L1).values & (te_a.obs.label_train != 'Unknown').values
    labs = sorted(set(y[imm]))
    res = {'test_donors': list(test_d), 'n_test_cells': int(te.sum()), 'n_test_immune': int(imm.sum())}
    for k, pr in [('majority', np.full(len(y), maj)), ('knn_pca', pred_k), ('scanvi', pred_s)]:
        res[k] = dict(acc=float(accuracy_score(y[imm], pr[imm])), macro_f1=float(f1_score(y[imm], pr[imm], labels=labs, average='macro', zero_division=0)),
                      acc_all=float(accuracy_score(y, pr)))
    res['per_label_f1_scanvi'] = dict(zip(labs, map(float, f1_score(y[imm], pred_s[imm], labels=labs, average=None, zero_division=0))))
    res['per_label_f1_knn'] = dict(zip(labs, map(float, f1_score(y[imm], pred_k[imm], labels=labs, average=None, zero_division=0))))
    g = res['scanvi']['macro_f1'] >= res['knn_pca']['macro_f1'] - 0.05 and res['scanvi']['macro_f1'] > res['majority']['macro_f1']
    res['gate_scanvi_primary'] = bool(g)
    print(json.dumps(res, indent=1)); json.dump(res, open(f'{wd}/k1_val.json', 'w'), indent=1)

if stage == 'full':
    Q = ad.read_h5ad(f'{wd}/k1_query_all.h5ad')
    genes = hvg(ref)
    a, s = fit_ref(ref, genes)
    q, sq, nm = map_query(s, Q, genes)
    print('query genes zero-filled', nm)
    Q.obs['pred_scanvi'] = sq.predict(q)
    pr = sq.predict(q, soft=True); Q.obs['maxprob'] = pr.max(1).values
    zr = s.get_latent_representation(a); zq = sq.get_latent_representation(q)
    nn = NearestNeighbors(n_neighbors=15).fit(zr)
    _, idx = nn.kneighbors(zq)
    rl = a.obs.label_train.values[idx]
    Q.obs['knn_unc'] = (rl != Q.obs.pred_scanvi.values[:, None]).mean(1)
    # kNN-PCA secondary labels
    Xr, Xq = lognorm(a), lognorm(q)
    p = PCA(50, random_state=0).fit(Xr)
    Q.obs['pred_knn_pca'] = KNeighborsClassifier(15).fit(p.transform(Xr), a.obs.label_train.values).predict(p.transform(Xq))
    np.save(f'{wd}/k1_zq.npy', zq)
    Q.obs.to_csv(f'{wd}/k1_query_obs.csv')
    print(Q.obs.groupby('sample').pred_scanvi.value_counts().unstack(fill_value=0).T.to_string())
