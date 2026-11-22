# Round 99 (Q14/Q29): blood IFN6 vs active nephritis beyond extra-renal SLEDAI + serology. Follows PREREG.md exactly.
# usage: python3 -I run.py <dir with GSE65391/GSE72326/GSE99967 series matrices + GPL10558.annot.gz> <out dir>
import sys, gzip, json, re, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
D, O = sys.argv[1], sys.argv[2]
IFN = ['IFI27', 'IFI44L', 'IFIT1', 'ISG15', 'RSAD2', 'SIGLEC1']
ENTREZ = {'IFI27': '3429', 'IFI44L': '10964', 'IFIT1': '3434', 'ISG15': '9636', 'RSAD2': '91543', 'SIGLEC1': '6614'}

def read_matrix(gse):
    meta, rows = {}, []
    with gzip.open(f'{D}/{gse}_series_matrix.txt.gz', 'rt') as f:
        chars = []
        for line in f:
            if line.startswith('!Sample_geo_accession'): meta['gsm'] = [x.strip('"') for x in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!Sample_title'): meta['title'] = [x.strip('"') for x in line.rstrip('\n').split('\t')[1:]]
            elif line.startswith('!Sample_characteristics_ch1'): chars.append([x.strip('"') for x in line.rstrip('\n').split('\t')[1:]])
            elif line.startswith('"ID_REF"'): hdr = [x.strip('"') for x in line.rstrip('\n').split('\t')]; break
        X = pd.read_csv(f, sep='\t', header=None, names=hdr, comment='!', index_col=0, dtype={0: str})
    X.index = X.index.astype(str).str.strip('"')
    recs = []
    for j in range(len(meta['gsm'])):
        r = {'gsm': meta['gsm'][j], 'title': meta['title'][j]}
        for c in chars:
            if ': ' in c[j]: k, v = c[j].split(': ', 1); r[k.strip().lower()] = v.strip()
        recs.append(r)
    return X, pd.DataFrame(recs).set_index('gsm')

def ifn_illumina(X):
    ann = pd.read_csv(f'{D}/GPL10558.annot.gz', sep='\t', skiprows=28, usecols=['ID', 'Gene symbol'], dtype=str).dropna()
    m = {p: s for p, s in zip(ann.ID, ann['Gene symbol']) if s in IFN}
    E = X.loc[X.index.intersection(list(m))].apply(pd.to_numeric, errors='coerce')
    E['gene'] = [m[p] for p in E.index]
    E = E.assign(mm=E.drop(columns='gene').mean(axis=1)).sort_values('mm').groupby('gene').tail(1).set_index('gene').drop(columns='mm')
    return E.T

def zscore_ifn(E, idx):
    E = E.loc[idx, IFN].astype(float)
    if np.nanmedian(E.values) > 50: E = np.log2(E.clip(lower=1))  # Illumina non-log safeguard
    z = (E - E.mean()) / E.std()
    return z.mean(axis=1)

num = lambda s: pd.to_numeric(s, errors='coerce')
out = {}

# ---- cohort A: GSE65391
X, M = read_matrix('GSE65391'); E = ifn_illumina(X)
M = M[M['disease state'] == 'SLE'].copy(); M['visit'] = num(M.visit)
A = M.sort_values('visit').groupby('subject').head(1).copy()
items = ['proteinuria', 'hematuria', 'pyuria', 'urinary_casts']
A['y'] = (num(A.renal) == 1).astype(float); A.loc[num(A.renal).isna(), 'y'] = np.nan
A['erSLEDAI'] = num(A.sledai) - sum(num(A[i]) for i in items)
A['lowC'] = (num(A.low_complement) > 0).astype(float); A.loc[num(A.low_complement).isna(), 'lowC'] = np.nan
A['dsDNA'] = (num(A.increased_dna_binding) > 0).astype(float); A.loc[num(A.increased_dna_binding).isna(), 'dsDNA'] = np.nan
A['dsDNA_c'] = A['dsDNA']
A = A.dropna(subset=['y', 'erSLEDAI', 'lowC', 'dsDNA'])
A['IFN6'] = zscore_ifn(E, A.index)
A['cohort'] = 'A_GSE65391'

# ---- cohort B: GSE72326
X, M = read_matrix('GSE72326'); E = ifn_illumina(X)
M = M[M['group'] == 'SLE'].copy()
M['pt'] = M.title.str.replace(r'V\d+$', '', regex=True); M['v'] = num(M.title.str.extract(r'V(\d+)$')[0])
B = M.sort_values('v').groupby('pt').head(1).copy()
yn = lambda s: s.map({'Y': 1.0, 'N': 0.0})
B['y'] = yn(B['renal'])
renal_items = yn(B['haematuria']) + yn(B['leukocyturia']) + yn(B['uriry casts']) + (num(B['24h proteinuria']) > 0.5).astype(float).where(num(B['24h proteinuria']).notna())
B['erSLEDAI'] = (num(B.sledai) - 4 * renal_items).clip(lower=0)
B['lowC'] = ((yn(B['low c3']) == 1) | (yn(B['low c4']) == 1)).astype(float).where(yn(B['low c3']).notna() | yn(B['low c4']).notna())
B['dsDNA_c'] = np.log1p(num(B['a-dsd']))
B['dsDNA'] = (num(B['a-dsd']) >= 10).astype(float).where(num(B['a-dsd']).notna())
B = B.dropna(subset=['y', 'erSLEDAI', 'lowC', 'dsDNA'])
B['IFN6'] = zscore_ifn(E, B.index)
B['cohort'] = 'B_GSE72326'

MODELS = {'T1 serology': ['lowC', 'dsDNA'], 'T0 erSLEDAI': ['erSLEDAI'], 'BASE': ['erSLEDAI', 'lowC', 'dsDNA'],
          'FULL': ['erSLEDAI', 'lowC', 'dsDNA', 'IFN6'], 'IFN-alone': ['IFN6']}
mk = lambda: make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=1000))
rng = np.random.default_rng(0)

def boot_delta(y, a, b, n=2000):
    ds = []; idx = np.arange(len(y))
    for _ in range(n):
        i = rng.choice(idx, len(idx))
        if y[i].min() == y[i].max(): continue
        ds.append(roc_auc_score(y[i], a[i]) - roc_auc_score(y[i], b[i]))
    return [round(float(np.percentile(ds, 2.5)), 3), round(float(np.percentile(ds, 97.5)), 3)]

def describe(d):
    return {'n': int(len(d)), 'n_active_LN': int(d.y.sum()), 'raw_auc_erSLEDAI': round(roc_auc_score(d.y, d.erSLEDAI), 3),
            'raw_auc_IFN6': round(roc_auc_score(d.y, d.IFN6), 3),
            'spearman_IFN6_erSLEDAI': round(float(spearmanr(d.IFN6, d.erSLEDAI)[0]), 3),
            'mean_IFN6_LN': round(float(d.IFN6[d.y == 1].mean()), 3), 'mean_IFN6_noLN': round(float(d.IFN6[d.y == 0].mean()), 3)}

out['describe'] = {'A': describe(A), 'B': describe(B)}
# LOCO
loco = {}
for tr, te, nm in [(A, B, 'train A -> test B'), (B, A, 'train B -> test A')]:
    P = {k: mk().fit(tr[f], tr.y).predict_proba(te[f])[:, 1] for k, f in MODELS.items()}
    y = te.y.values; r = {k: round(roc_auc_score(y, p), 3) for k, p in P.items()}
    r['FULL-BASE'] = round(r['FULL'] - r['BASE'], 3); r['FULL-BASE_CI'] = boot_delta(y, P['FULL'], P['BASE'])
    r['IFN-T1_CI'] = boot_delta(y, P['IFN-alone'], P['T1 serology'])
    loco[nm] = r
out['LOCO'] = loco
# within-cohort repeated CV
cv = {}
for d, nm in [(A, 'A'), (B, 'B')]:
    res = {k: [] for k in MODELS}
    for s in range(20):
        oof = {k: np.zeros(len(d)) for k in MODELS}
        for tri, tei in StratifiedKFold(5, shuffle=True, random_state=s).split(d, d.y):
            for k, f in MODELS.items(): oof[k][tei] = mk().fit(d.iloc[tri][f], d.y.iloc[tri]).predict_proba(d.iloc[tei][f])[:, 1]
        for k in MODELS: res[k].append(roc_auc_score(d.y, oof[k]))
    cv[nm] = {k: round(float(np.mean(v)), 3) for k, v in res.items()}
    cv[nm]['FULL-BASE'] = round(cv[nm]['FULL'] - cv[nm]['BASE'], 3)
out['withinCV'] = cv
# sensitivity: continuous dsDNA in B (not harmonised with A, so within-B CV only)
Bs = B.dropna(subset=['dsDNA_c']); res = {'BASEc': [], 'FULLc': []}
for s in range(20):
    oof = {k: np.zeros(len(Bs)) for k in res}
    for tri, tei in StratifiedKFold(5, shuffle=True, random_state=s).split(Bs, Bs.y):
        for k, f in [('BASEc', ['erSLEDAI', 'lowC', 'dsDNA_c']), ('FULLc', ['erSLEDAI', 'lowC', 'dsDNA_c', 'IFN6'])]:
            oof[k][tei] = mk().fit(Bs.iloc[tri][f], Bs.y.iloc[tri]).predict_proba(Bs.iloc[tei][f])[:, 1]
    for k in res: res[k].append(roc_auc_score(Bs.y, oof[k]))
out['sens_B_continuous_dsDNA'] = {k: round(float(np.mean(v)), 3) for k, v in res.items()} | {'n': int(len(Bs))}

# ---- cohort C: GSE99967 (secondary)
X, M = read_matrix('GSE99967')
M = M[M['disease state'] == 'SLE'].copy(); M['y'] = M['nephritis state'].map({'active LN (ALN)': 1.0, 'without LN (ANLN)': 0.0})
M = M.dropna(subset=['y'])
E = X.loc[[ENTREZ[g] + '_at' for g in IFN], M.index].astype(float).T; E.columns = IFN
M['IFN6'] = ((E - E.mean()) / E.std()).mean(axis=1)
yC, sC = M.y.values, M.IFN6.values; bs = []
for _ in range(2000):
    i = rng.choice(len(yC), len(yC))
    if yC[i].min() != yC[i].max(): bs.append(roc_auc_score(yC[i], sC[i]))
out['C'] = {'n': int(len(M)), 'n_ALN': int(yC.sum()), 'auc_IFN6': round(roc_auc_score(yC, sC), 3),
            'auc_CI': [round(float(np.percentile(bs, 2.5)), 3), round(float(np.percentile(bs, 97.5)), 3)],
            'raw_auc_total_SLEDAI_(includes renal items)': round(roc_auc_score(yC, num(M.sledai)), 3)}

# ---- verdicts
L = out['LOCO'].values()
out['H1_PASS'] = all(r['FULL-BASE'] >= 0.05 and r['FULL-BASE_CI'][0] > 0 for r in L)
out['H2_PASS'] = all(r['IFN-alone'] > r['T1 serology'] and r['IFN-alone'] >= 0.60 for r in L)
out['H3_PASS'] = out['C']['auc_IFN6'] >= 0.65 and out['C']['auc_CI'][0] > 0.5
pd.concat([A, B])[['cohort', 'y', 'erSLEDAI', 'lowC', 'dsDNA', 'IFN6']].to_csv(f'{O}/patients.csv')
json.dump(out, open(f'{O}/result.json', 'w'), indent=1); print(json.dumps(out, indent=1))
