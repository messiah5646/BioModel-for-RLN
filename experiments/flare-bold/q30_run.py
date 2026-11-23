"""Round 122 (Q30): baseline TNFRSF17 (plasmablast) score vs time to flare in BOLD (GSE92776), per PREREG-Q30.md + amendment 1.
Usage: python q30_run.py /path/to/GSE92776_series_matrix.txt.gz /path/to/GPL22843.txt > q30_result.txt"""
import sys, gzip, json, numpy as np, pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy.stats import wilcoxon, spearmanr

mx, gpl = sys.argv[1], sys.argv[2]
meta, rows, inside = {}, [], False
for line in gzip.open(mx, 'rt'):
    line = line.rstrip('\n')
    if line.startswith('!series_matrix_table_begin'): inside = True; continue
    if line.startswith('!series_matrix_table_end'): break
    if inside: rows.append([c.strip('"') for c in line.split('\t')]); continue
    if line.startswith('!Sample_geo_accession'): gsm = [c.strip('"') for c in line.split('\t')[1:]]
    if line.startswith('!Sample_characteristics_ch1'):
        for j, c in enumerate(line.split('\t')[1:]):
            c = c.strip('"')
            if ': ' in c:
                k, v = c.split(': ', 1); meta.setdefault(k, [None] * len(gsm))[j] = v
M = pd.DataFrame(meta, index=gsm)
X = pd.DataFrame([r[1:] for r in rows[1:]], index=[r[0] for r in rows[1:]], columns=rows[0][1:]).replace('', np.nan).astype(float)
al = [l.rstrip('\n').split('\t') for l in open(gpl) if '\t' in l and not l.startswith(('#', '!', '^', 'ID\t'))]
ann = pd.Series({a[0]: a[1] for a in al if len(a) >= 2})
E = (-X).groupby(ann.reindex(X.index).values).mean().T          # samples x genes, expression = -dCt
IFN = ['IFI27', 'IFI44L', 'IFIT1', 'ISG15', 'RSAD2', 'SIGLEC1']

M['date'] = pd.to_datetime(M['date_of_visit_sample_drawn_date']); M['sledai'] = pd.to_numeric(M.get('sledai_score'), errors='coerce')
M['lym'] = pd.to_numeric(M.get('lymphocytes'), errors='coerce')
A = M[(M.subject_type == 'Patient') & (M.group_id == 'A')].copy()
out = {'n_groupA_samples': len(A)}
pts = []
for d, g in A.groupby('donor_id'):
    b = g[g.nicevisit == 'BASELINE']
    if len(b) == 0: continue
    b = b.sort_values('date').iloc[0]; f = g[(g.nicevisit == 'FLARE') & (g.date > b.date)].sort_values('date')
    ev = int(len(f) > 0); tend = f.date.iloc[0] if ev else g.date.max()
    pts.append(dict(donor=d, gsm=b.name, event=ev, T=(tend - b.date).days, sledai=b.sledai, lym=b.lym))
P = pd.DataFrame(pts).set_index('gsm')
P['tnf17_raw'] = E.loc[P.index, 'TNFRSF17']
out['n_baseline_pts'] = len(P); out['tnfrsf17_missing_frac'] = float(P.tnf17_raw.isna().mean())
Zi = E.loc[P.index, IFN]; Zi = (Zi - Zi.mean()) / Zi.std(); P['IFN6'] = Zi.mean(axis=1)
P = P[P['T'] > 0]
if out['tnfrsf17_missing_frac'] > 0.2:
    out['H1'] = 'UNDERPOWERED-DATA'; print(json.dumps(out, indent=1)); sys.exit()
P = P.dropna(subset=['tnf17_raw', 'sledai']); P['PB'] = (P.tnf17_raw - P.tnf17_raw.mean()) / P.tnf17_raw.std()
out.update(n_analysed=len(P), events=int(P.event.sum()), median_T_days=float(P[P.event == 1]['T'].median()))
if P.event.sum() < 20: out['H1'] = 'UNDERPOWERED'

def loo_c(cols, D):
    lp = []
    for i in D.index:
        tr = D.drop(i); cph = CoxPHFitter(penalizer=0.01).fit(tr[cols + ['T', 'event']], 'T', 'event')
        lp.append(float(cph.predict_partial_hazard(D.loc[[i], cols]).iloc[0]))
    return concordance_index(D['T'], -np.array(lp), D['event'])

C = {k: loo_c(v, P) for k, v in {'M1_sledai': ['sledai'], 'M2_sledai_PB': ['sledai', 'PB'], 'M3_sledai_IFN6': ['sledai', 'IFN6'],
                                 'R2_PB_alone': ['PB'], 'IFN6_alone': ['IFN6']}.items()}
C['M0_null'] = 0.5; out['LOO_C'] = C
full = CoxPHFitter(penalizer=0.01).fit(P[['sledai', 'PB', 'T', 'event']], 'T', 'event'); out['M2_PB_logHR_full'] = float(full.params_['PB'])
rng = np.random.default_rng(122); bs = []
for _ in range(2000):
    S = P.sample(len(P), replace=True, random_state=rng.integers(1 << 31)).reset_index(drop=True)
    if S.event.sum() < 5: continue
    try: bs.append(CoxPHFitter(penalizer=0.01).fit(S[['sledai', 'PB', 'T', 'event']], 'T', 'event').params_['PB'])
    except Exception: pass
out['M2_PB_logHR_CI95'] = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]; out['n_boot_ok'] = len(bs)
dC = C['M2_sledai_PB'] - C['M1_sledai']; out['delta_C_M2_M1'] = dC
if 'H1' not in out:
    out['H1'] = 'PASS' if (dC >= 0.05 and out['M2_PB_logHR_CI95'][0] > 0) else 'FAIL'
# R1 lymphocyte-adjusted
L = P.dropna(subset=['lym'])
if len(L) >= 15 and L.event.sum() >= 10:
    out['R1_n'] = len(L); out['R1_LOO_C'] = {'sledai_lym': loo_c(['sledai', 'lym'], L), 'sledai_lym_PB': loo_c(['sledai', 'lym', 'PB'], L)}
# R3 among flarers
F = P[P.event == 1]; r, p = spearmanr(F.PB, F['T']); out['R3_spearman_PB_vs_days_flarers'] = [float(r), float(p), len(F)]
r, p = spearmanr(F.IFN6, F['T']); out['R3_spearman_IFN6_vs_days_flarers'] = [float(r), float(p)]
# H2 within patient (z across all group-A samples)
ga = E.loc[A.index, 'TNFRSF17']; A['PBz'] = (ga - ga.mean()) / ga.std()
h2 = {}
for other in ['BASELINE', 'IMPROVING']:
    diffs = []
    for d, g in A.groupby('donor_id'):
        f = g[g.nicevisit == 'FLARE'].PBz.dropna(); o = g[g.nicevisit == other].PBz.dropna()
        if len(f) and len(o): diffs.append(f.iloc[0] - o.mean())
    diffs = np.array(diffs)
    h2[f'FLARE_minus_{other}'] = dict(n=len(diffs), median=float(np.median(diffs)), wilcoxon_p=float(wilcoxon(diffs).pvalue) if len(diffs) > 5 else None)
out['H2'] = h2
out['sledai_baseline_median'] = float(P.sledai.median())
P.to_csv('q30_patients.csv')
print(json.dumps(out, indent=1))
