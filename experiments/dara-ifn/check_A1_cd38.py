"""Round 114 self-audit: verify assumption A1 (GEO 'Patient N' == paper 'SLE #N') with independent data.
Paper Fig 3b (flow): % CD38+ among CD4 T cells at week 0, per SLE #. Sequencing: CITE-seq CD38 ADT on TCR+ CD4+ cells
in week-0 single-patient libraries. Usage: python check_A1_cd38.py DATADIR"""
import sys, glob, numpy as np, pandas as pd, h5py, scipy.sparse as sp, itertools
from scipy.stats import spearmanr
D = sys.argv[1]
FLOW_W0 = {1: 52.7, 2: 34.6, 3: 5.56, 4: 20.5, 5: 39.2, 6: 15.9, 7: 30, 8: 13.8, 9: 69.5, 10: 27.1}
W0_SINGLE = {'GSM8919721': 1, 'GSM8919724': 2, 'GSM8919733': 3, 'GSM8919739': 4, 'GSM8919757': 9}
TCR = {'GSM8919721': 'GSM8919723', 'GSM8919724': 'GSM8919726', 'GSM8919733': 'GSM8919735', 'GSM8919739': 'GSM8919741', 'GSM8919757': 'GSM8919759'}
rows = []
for g, p in W0_SINGLE.items():
    m = h5py.File(glob.glob(f'{D}/{g}_*.h5')[0], 'r')['matrix']
    X = sp.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=m['shape'][:]).tocsr()
    n = np.array(m['features']['name']).astype(str); ft = np.array(m['features']['feature_type']).astype(str)
    bc = np.array(m['barcodes']).astype(str)
    Xg = X[ft == 'Gene Expression']; gn = n[ft == 'Gene Expression']
    def gene(g): return np.asarray(Xg[np.where(gn == g)[0][0]].todense()).ravel()
    t = pd.read_csv(glob.glob(f'{D}/{TCR[g]}_*.csv.gz')[0]); tb = set(t[t.productive.astype(str).str.lower() == 'true'].barcode)
    tot = np.asarray(Xg.sum(0)).ravel()
    isT = np.isin(bc, list(tb)) & (tot > 500)
    sel = isT & (gene('CD8A') + gene('CD8B') == 0)  # CD4-lineage proxy: TCR+, no CD8 mRNA
    cd38 = gene('CD38')
    rows.append(dict(patient=p, n_cd4=int(sel.sum()), cd38_mrna_pos_frac=float((cd38[sel] > 0).mean()),
                     cd38_cp10k=float((cd38[sel] / tot[sel] * 1e4).mean()), flow_pct=FLOW_W0[p]))
r = pd.DataFrame(rows); print(r.round(3).to_string())
for c in ['cd38_mrna_pos_frac', 'cd38_cp10k']:
    rho = spearmanr(r[c], r.flow_pct)[0]
    perms = [spearmanr(r[c], r.flow_pct.values[list(pp)])[0] for pp in itertools.permutations(range(len(r)))]
    print(f'{c}: Spearman vs paper flow = {rho:.2f}; exact p(rho >= obs) = {np.mean(np.array(perms) >= rho - 1e-12):.3f}')
