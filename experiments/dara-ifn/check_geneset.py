"""Round 114 self-audit: is the P2b result an artefact of the IFN6 gene set (IFI27, SIGLEC1 are myeloid genes)?
Rescore T cells with the paper's own T-cell IFN genes (MX1, IFI44L, XAF1, ISG15, IFIT1; Ostendorf 2026 Fig 3I text).
Only clonotype-confirmed samples. NOT pre-registered; a robustness check. Usage: python check_geneset.py DATADIR"""
import sys, glob, numpy as np, pandas as pd, h5py, scipy.sparse as sp
from scipy.stats import spearmanr
D = sys.argv[1]
SETS = {'IFN6': ['IFI27', 'IFI44L', 'IFIT1', 'ISG15', 'RSAD2', 'SIGLEC1'], 'paper5': ['MX1', 'IFI44L', 'XAF1', 'ISG15', 'IFIT1']}
# (expr GSM, TCR GSM, patient, week, hashtag filter) ; P1/P2 w36 share GSM8919748 (hashtag 1 = P1, 2 = P2; demux.csv)
S = [('GSM8919721', 'GSM8919723', 1, 0, None), ('GSM8919727', 'GSM8919729', 1, 9, None), ('GSM8919748', 'GSM8919750', 1, 36, 1),
     ('GSM8919724', 'GSM8919726', 2, 0, None), ('GSM8919730', 'GSM8919732', 2, 9, None), ('GSM8919748', 'GSM8919750', 2, 36, 2),
     ('GSM8919733', 'GSM8919735', 3, 0, None), ('GSM8919736', 'GSM8919738', 3, 9, None), ('GSM8919754', 'GSM8919756', 3, 36, None),
     ('GSM8919739', 'GSM8919741', 4, 0, None), ('GSM8919745', 'GSM8919747', 4, 9, None),
     ('GSM8919757', 'GSM8919759', 9, 0, None), ('GSM8919772', 'GSM8919774', 9, 36, None)]
DS = {1: (331, 73.4, 94.6), 2: (1769, 201, 127.2), 3: (64.9, 27.8, 31.3), 4: (167.1, 61.7, 66.2), 9: (144.7, 114, 197.4)}
rows = []
for ge, tc, p, w, ht in S:
    m = h5py.File(glob.glob(f'{D}/{ge}_*.h5')[0], 'r')['matrix']
    X = sp.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=m['shape'][:]).tocsr()
    n = np.array(m['features']['name']).astype(str); ft = np.array(m['features']['feature_type']).astype(str)
    bc = np.array(m['barcodes']).astype(str)
    Xg = X[ft == 'Gene Expression']; gn = n[ft == 'Gene Expression']; tot = np.asarray(Xg.sum(0)).ravel()
    t = pd.read_csv(glob.glob(f'{D}/{tc}_*.csv.gz')[0]); tb = set(t[t.productive.astype(str).str.lower() == 'true'].barcode)
    sel = np.isin(bc, list(tb)) & (tot > 500)
    if ht:
        A = X[ft == 'Antibody Capture']; an = list(n[ft == 'Antibody Capture'])
        h1 = np.asarray(A[an.index('Hashtag1')].todense()).ravel(); h2 = np.asarray(A[an.index('Hashtag2')].todense()).ravel()
        sel &= (h1 > h2) if ht == 1 else (h2 > h1)
    r = dict(patient=p, week=w, nT=int(sel.sum()))
    for k, G in SETS.items():
        r[k] = float(np.mean([np.log1p(np.asarray(Xg[np.where(gn == g)[0][0]].todense()).ravel()[sel] / tot[sel] * 1e4) for g in G], axis=0).mean())
        for g in G: r[f'{k}:{g}_pct_pos'] = None
    rows.append(r)
d = pd.DataFrame(rows)[['patient', 'week', 'nT', 'IFN6', 'paper5']]
print(d.round(3).to_string())
for k in SETS:
    w = d.pivot(index='patient', columns='week', values=k)
    m9 = (w[9] / w[0]).dropna(); m36 = (w[36] / w[0]).dropna()
    rho = spearmanr(np.log([DS[p][2] / DS[p][0] for p in m36.index]), np.log(m36.values))[0]
    print(f'{k}: m9 {m9.round(2).to_dict()} median {m9.median():.2f} | m36 {m36.round(2).to_dict()} median {m36.median():.2f} | rho36 vs anti-dsDNA {rho:.2f} (n {len(m36)})')
