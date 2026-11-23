"""P5 analysis (push 3) per PREREG-P5.md. Usage: python3 -I analyse_p5.py DATADIR OUTDIR
DATADIR: GSE88885/6 series matrices, preprocessed expression, gpl_data.txt (GPL17586 'view=data' table)."""
import sys, os, gzip, re, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr

D, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
G = ['IFI27', 'IFI44L', 'IFIT1', 'ISG15', 'RSAD2', 'SIGLEC1']

# ---- probeset -> gene (GPL17586 gene_assignment: 'acc // SYMBOL // ...') ----
_start = next(i for i, l in enumerate(open(os.path.join(D, 'gpl_data.txt'))) if l.startswith('!platform_table_begin'))
gpl = pd.read_csv(os.path.join(D, 'gpl_data.txt'), sep='\t', skiprows=_start + 1, low_memory=False,
                  usecols=['ID', 'gene_assignment'], on_bad_lines='skip')
gpl = gpl.dropna(subset=['gene_assignment'])
gpl['sym'] = gpl.gene_assignment.astype(str).str.split(' // ').str[1]
ids = {g: set(gpl.ID[gpl.sym == g]) for g in G}

# ---- sample metadata ----
def meta(g):
    rows = [l.rstrip('\n').split('\t') for l in gzip.open(os.path.join(D, f'{g}.matrix.gz'), 'rt') if l.startswith('!Sample_')]
    m = {}
    for r in rows:
        v = [x.strip('"') for x in r[1:]]
        if r[0] == '!Sample_characteristics_ch1': m[v[0].split(':')[0]] = [x.split(': ', 1)[1] if ': ' in x else '' for x in v]
        elif r[0] in ('!Sample_geo_accession', '!Sample_supplementary_file'): m[r[0]] = v
    m = pd.DataFrame(m)
    m['cel'] = m['!Sample_supplementary_file'].map(lambda s: re.sub(r'\.gz$', '', s.split('/')[-1]).split('_', 1)[-1])
    m['trial'] = g
    m['arm'] = m[[c for c in m.columns if c.startswith('treatment')][0]]
    return m.rename(columns={'subject_id': 'subj', 'time': 'visit', 'antidsdna_at_baseline (iu)': 'dsdna', 'c3_at_baseline (g/l)': 'c3'})

def scores(g, fname):
    m = meta(g)
    want = set().union(*ids.values())
    rows, hdr = {}, None
    with gzip.open(os.path.join(D, fname), 'rt') as f:
        hdr = [h.strip('"') for h in f.readline().rstrip('\n').split('\t')]
        for line in f:
            pid = line.split('\t', 1)[0].strip('"').split(':')[0]
            if pid in want: rows[pid] = np.array(line.rstrip('\n').split('\t')[1:], dtype=float)
    X = pd.DataFrame(rows, index=hdr[1:]).T  # probesets x CEL
    genes = {}
    for gname, ps in ids.items():
        p = [x for x in ps if x in X.index]
        if p: genes[gname] = X.loc[p].mean(0)
    E = pd.DataFrame(genes)  # CEL x gene
    E.index = [c.split('_', 1)[-1] if '_' in c else c for c in E.index]
    m = m.set_index('cel'); common = [c for c in m.index if c in E.index]
    E = E.loc[common]; m = m.loc[common]
    base = m.visit == 'baseline'
    Z = E - E[base.values].median(0)  # gene minus its baseline median within trial
    m['ifn6'] = Z.mean(1).values
    return m.reset_index(), list(genes), len(common)

out = {}
data = []
for g, fn in [('GSE88885', 'GSE88885_ILLUMINATE1_pd_preprocessed.txt.gz'), ('GSE88886', 'GSE88886_ILLUMINATE2_pd_preprocessed.txt.gz')]:
    m, genes_used, n_matched = scores(g, fn)
    out[f'{g}_genes_used'] = genes_used; out[f'{g}_samples_matched'] = n_matched
    data.append(m)
M = pd.concat(data, ignore_index=True)
W = M.pivot_table(index=['trial', 'subj'], columns='visit', values='ifn6').reset_index()
# arm = the (single) treatment recorded on post-baseline samples
arm = M[M.visit != 'baseline'].groupby(['trial', 'subj']).arm.agg(lambda s: s.iloc[0])
W = W.merge(arm.rename('arm').reset_index(), on=['trial', 'subj'], how='left')
cov = M[M.visit == 'baseline'][['trial', 'subj', 'dsdna', 'c3']]
W = W.merge(cov, on=['trial', 'subj'], how='left')
W.to_csv(os.path.join(OUT, 'per_patient.csv'), index=False)

# ---- A: persistence in placebo ----
A = {}
for g in ['GSE88885', 'GSE88886']:
    p = W[(W.trial == g) & (W.arm == 'Placebo')].dropna(subset=['baseline', 'week52'])
    rho = float(spearmanr(p.baseline, p.week52)[0])
    p3 = W[(W.trial == g) & (W.arm == 'Placebo')].dropna(subset=['baseline', 'week16', 'week52'])
    L = p3[['baseline', 'week16', 'week52']].values; k = 3; n = len(L)
    msb = k * ((L.mean(1) - L.mean()) ** 2).sum() / (n - 1); msw = ((L - L.mean(1, keepdims=True)) ** 2).sum() / (n * (k - 1))
    icc = float((msb - msw) / (msb + (k - 1) * msw))
    A[g] = dict(n=len(p), rho_base_w52=round(rho, 3), icc_3visits=round(icc, 3), n_icc=n)
passA = [v['rho_base_w52'] >= 0.5 for v in A.values()]
A['verdict'] = 'PASS' if all(passA) else ('FAIL' if not any(passA) else 'inconclusive')
out['A_persistence_placebo'] = A

# ---- B: out-of-trial prediction of week 52 ----
B = {}
for tr, te in [('GSE88885', 'GSE88886'), ('GSE88886', 'GSE88885')]:
    T = W[(W.trial == tr)].dropna(subset=['baseline', 'week52', 'arm']); S = W[(W.trial == te)].dropna(subset=['baseline', 'week52', 'arm'])
    am = T.groupby('arm').week52.mean(); ch = (T.week52 - T.baseline).groupby(T.arm).mean()
    Xtr = np.column_stack([np.ones(len(T)), T.baseline] + [(T.arm == a).astype(float) for a in sorted(T.arm.unique())[1:]])
    beta, *_ = np.linalg.lstsq(Xtr, T.week52.values, rcond=None)
    Xte = np.column_stack([np.ones(len(S)), S.baseline] + [(S.arm == a).astype(float) for a in sorted(T.arm.unique())[1:]])
    pred = dict(M0=S.arm.map(am).values, M1=S.baseline.values, M2=(S.baseline + S.arm.map(ch)).values, M3=Xte @ beta)
    mse = {k: round(float(np.mean((S.week52.values - v) ** 2)), 4) for k, v in pred.items()}
    best = min(mse['M2'], mse['M3'])
    B[f'train_{tr}_test_{te}'] = dict(n_test=len(S), mse=mse, M3_slope_on_baseline=round(float(beta[1]), 3),
                                      setpoint_gain_vs_M0=round(1 - best / mse['M0'], 3))
gains = [v['setpoint_gain_vs_M0'] for v in B.values()]
B['verdict'] = ('set-point model WINS' if all(x >= 0.10 for x in gains) else
                'set-point model LOSES' if all(min(v['mse']['M2'], v['mse']['M3']) > v['mse']['M0'] for k, v in B.items() if k != 'verdict') else 'inconclusive')
out['B_out_of_trial_w52'] = B
out['DECISION'] = ('keep set-point layer' if A['verdict'] == 'PASS' and B['verdict'].endswith('WINS') else
                   'drop set-point layer' if A['verdict'] == 'FAIL' else 'keep as Emerging, no change')

# ---- secondary C, D, E ----
C = {}
for g in ['GSE88885', 'GSE88886']:
    d = W[W.trial == g].dropna(subset=['baseline', 'week52', 'arm']).copy(); d['chg'] = d.week52 - d.baseline
    pl = d[d.arm == 'Placebo'].chg; tb = d[d.arm != 'Placebo'].chg
    diff = tb.mean() - pl.mean(); se = np.sqrt(tb.var() / len(tb) + pl.var() / len(pl))
    C[g] = dict(placebo_change=round(float(pl.mean()), 3), tabalumab_change=round(float(tb.mean()), 3),
                tab_minus_placebo=round(float(diff), 3), ci95=[round(float(diff - 1.96 * se), 3), round(float(diff + 1.96 * se), 3)],
                n_pl=len(pl), n_tab=len(tb))
out['C_arm_effect'] = C
Dd = {}
for g in ['GSE88885', 'GSE88886']:
    d = W[W.trial == g].dropna(subset=['baseline', 'week16', 'week52'])
    Dd[g] = dict(n=len(d), rho_w16_w52=round(float(spearmanr(d.week16, d.week52)[0]), 3), rho_base_w52=round(float(spearmanr(d.baseline, d.week52)[0]), 3),
                 mse_carry_w16=round(float(np.mean((d.week52 - d.week16) ** 2)), 4), mse_carry_base=round(float(np.mean((d.week52 - d.baseline) ** 2)), 4))
out['D_w16_vs_baseline'] = Dd
E = {}
for g in ['GSE88885', 'GSE88886']:
    d = W[W.trial == g].dropna(subset=['baseline', 'week52']).copy(); d['chg'] = d.week52 - d.baseline
    for c in ['dsdna', 'c3']:
        x = pd.to_numeric(d[c], errors='coerce'); ok = x.notna()
        E[f'{g}_{c}_vs_change'] = round(float(spearmanr(x[ok], d.chg[ok])[0]), 3)
out['E_serology_vs_IFN_change'] = E
json.dump(out, open(os.path.join(OUT, 'result.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
