"""Q57: regenerate round-91 numbers (result.txt) from open GEO files for GSE200306.
Usage: python regenerate_r91.py DATADIR   (DATADIR holds matrix.txt.gz and the two normalized *_expression_*.txt.gz)
Round 91 committed no script; this is a reconstruction from its PREREG + RESULT text. Choices not stated there are
marked RECONSTRUCTED below and reported, so a mismatch can be traced."""
import sys, os, re, gzip, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

D = sys.argv[1]
# ---- sample metadata: patient, biopsy (A=1st, B=2nd), compartment, response, proteinuria, creatinine ----
rows = [l.rstrip('\n').split('\t') for l in gzip.open(os.path.join(D, 'matrix.txt.gz'), 'rt') if l.startswith('!Sample_')]
meta = {}
for r in rows:
    v = [x.strip('"') for x in r[1:]]
    if r[0] == '!Sample_characteristics_ch1' and ': ' in v[0]:
        meta[v[0].split(': ')[0]] = [x.split(': ', 1)[1] if ': ' in x else np.nan for x in v]
    elif r[0] == '!Sample_title':
        meta['title'] = v
    elif r[0] == '!Sample_geo_accession':
        meta['gsm'] = v
meta = pd.DataFrame(meta)
def parse_title(t):
    m = re.match(r'\s*patient (\d+), (first|second) biopsy, (glomeruli|tubulointerstitium)', t, re.I)
    return (int(m.group(1)), 'A' if m.group(2).lower() == 'first' else 'B', 'G' if m.group(3).lower().startswith('g') else 'TI') if m else (None, None, None)
meta[['patient', 'bx', 'comp']] = meta.title.apply(lambda t: pd.Series(parse_title(t)))
meta = meta.dropna(subset=['patient'])
meta['resp'] = meta['clinical response group'].str.lower().map({'complete response': 'CR', 'partial response': 'PR', 'no response': 'NR'})
meta['prot'] = pd.to_numeric(meta['proteinuria (g/24hr)'], errors='coerce')
meta['creat'] = pd.to_numeric(meta['serum creatinine (mg/dl)'], errors='coerce')
resp = meta.groupby('patient').resp.first()
print('matched keys:', meta[['patient', 'bx']].drop_duplicates().shape[0], dict(meta.drop_duplicates(['patient', 'bx']).resp.value_counts()))

# ---- expression: column name -> (patient, biopsy). RECONSTRUCTED parse: X|Avg_ + digits + A/a/B/b ----
def load(fname):
    e = pd.read_csv(os.path.join(D, fname), sep='\t').set_index('Name')
    e = e[e['Code.Class'] == 'Endogenous']
    keep = {}
    for c in e.columns:
        m = re.match(r'^(?:X|Avg_)0*(\d+)([AaBb])', c)
        if m and not re.search(r'XN|Y', c): keep[c] = (int(m.group(1)), m.group(2).upper())
    x = e[list(keep)].astype(float).T
    x.index = pd.MultiIndex.from_tuples([keep[c] for c in x.index], names=['patient', 'bx'])
    return x.groupby(level=[0, 1]).mean()
SOURCE = sys.argv[2] if len(sys.argv) > 2 else 'matrix'
def load_matrix(comp):
    lines = [l.rstrip('\n').split('\t') for l in gzip.open(os.path.join(D, 'matrix.txt.gz'), 'rt')]
    i0 = next(i for i, l in enumerate(lines) if l[0] == '!series_matrix_table_begin')
    i1 = next(i for i, l in enumerate(lines) if l[0] == '!series_matrix_table_end')
    hdr = [x.strip('"') for x in lines[i0 + 1]]
    tab = pd.DataFrame([[x.strip('"') for x in l] for l in lines[i0 + 2:i1]], columns=hdr).set_index('ID_REF').replace('', np.nan).astype(float)
    mm = meta[meta.comp == comp].set_index('gsm')
    x = tab[[g for g in mm.index if g in tab.columns]].T
    x.index = pd.MultiIndex.from_arrays([mm.loc[x.index, 'patient'].astype(int), mm.loc[x.index, 'bx']], names=['patient', 'bx'])
    return x.groupby(level=[0, 1]).mean()
if SOURCE == 'matrix':
    G, T = load_matrix('G'), load_matrix('TI')
else:
    G = load([f for f in os.listdir(D) if 'Glom' in f][0]); T = load([f for f in os.listdir(D) if '_TI_' in f][0])
G, T = G.dropna(axis=1), T.dropna(axis=1)  # drop genes with any missing value
print('source:', SOURCE, 'glom', G.shape, 'TI', T.shape)

MAC = ['CD163', 'CD14', 'CSF1R', 'C1QA', 'C1QB']
def score(x): return x[[g for g in MAC if g in x.columns]].mean(1)
print('H1 as pre-registered (CD68):', 'NOT TESTABLE - CD68 not on the panel.' if 'CD68' not in G.columns else 'testable')
for name, x in [('glomerular', G), ('TI', T)]:
    s = score(x).xs('A', level='bx'); r = resp.reindex(s.index)
    cr, nr = s[r == 'CR'], s[r == 'NR']
    print(f'EXPLORATORY {name} macrophage score Bx1: CR n={len(cr)} median {cr.median():.2f} | NR n={len(nr)} median {nr.median():.2f} | '
          f'p(CR>NR)={mannwhitneyu(cr, nr, alternative="greater").pvalue:.3f} p(NR>CR)={mannwhitneyu(nr, cr, alternative="greater").pvalue:.3f}')

# ---- H2: CR vs not-CR at Bx1; clinical from the Bx1 metadata row of the same compartment (RECONSTRUCTED) ----
def h2(x, comp, label):
    X1 = x.xs('A', level='bx')
    cl = meta[(meta.bx == 'A') & (meta.comp == comp)].drop_duplicates('patient').set_index('patient')[['prot', 'creat']]
    idx = [p for p in X1.index if p in cl.index and cl.loc[p].notna().all() and pd.notna(resp.get(p))]
    y = np.array([resp[p] == 'CR' for p in idx]).astype(int); C = cl.loc[idx].values; E = X1.loc[idx].values
    p0, p1 = np.zeros(len(y)), np.zeros(len(y))
    for i in range(len(y)):
        tr = np.arange(len(y)) != i
        sc = StandardScaler().fit(C[tr]); m0 = LogisticRegression(max_iter=1000).fit(sc.transform(C[tr]), y[tr])
        p0[i] = m0.predict_proba(sc.transform(C[i:i + 1]))[0, 1]
        se = StandardScaler().fit(E[tr]); pca = PCA(3).fit(se.transform(E[tr]))
        Z = np.hstack([C, np.zeros((len(y), 3))]); Ztr = np.hstack([C[tr], pca.transform(se.transform(E[tr]))])
        Zi = np.hstack([C[i:i + 1], pca.transform(se.transform(E[i:i + 1]))])
        sz = StandardScaler().fit(Ztr); m1 = LogisticRegression(max_iter=1000).fit(sz.transform(Ztr), y[tr])
        p1[i] = m1.predict_proba(sz.transform(Zi))[0, 1]
    a0, a1 = roc_auc_score(y, p0), roc_auc_score(y, p1)
    print(f'H2 {label}: n={len(y)} (CR={y.sum()}) | clinical AUC {a0:.3f} | +3 PCs AUC {a1:.3f} | gain {a1 - a0:+.3f} -> {"PASS" if a1 - a0 >= 0.05 else "FAIL"}')
h2(T, 'TI', 'TI'); h2(G, 'G', 'glom')

# ---- H3: TI at Bx2, NR vs CR, effect size = standardised mean difference (RECONSTRUCTED) ----
X2 = T.xs('B', level='bx'); r2 = resp.reindex(X2.index)
nr, cr = X2[r2 == 'NR'], X2[r2 == 'CR']
sd = np.sqrt(((len(nr) - 1) * nr.var() + (len(cr) - 1) * cr.var()) / (len(nr) + len(cr) - 2))
d = ((nr.mean() - cr.mean()) / sd).dropna().sort_values()
print(f'H3 TI at Bx2, NR n={len(nr)} vs CR n={len(cr)}. Higher in NR: {d[::-1][:10].round(2).to_dict()}')
print(f'   Higher in CR: {d[:5].round(2).to_dict()}')
rk = d[::-1].reset_index(); rk.columns = ['gene', 'd']
for g in ['TGFB1', 'CSF1R', 'CD99', 'MAP4K4', 'PTPN6', 'HLA-DRA', 'IFI16', 'LILRB1', 'LILRB2', 'C1S']:
    r = rk.index[rk.gene == g]
    print(f'   rank of {g} (higher in NR): {int(r[0]) + 1 if len(r) else "absent"} of {len(rk)}, d {float(rk.d[r[0]]) if len(r) else float("nan"):.2f}')
