"""Round 121 (Q25). Usage: python3 -I q25_trajectory.py <data_dir> <out_dir>. Follows PREREG-Q25-trajectory.md."""
import sys, gzip, json, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
D, O = sys.argv[1], sys.argv[2]
IFN = ['IFI27','IFI44L','IFIT1','ISG15','RSAD2','SIGLEC1']
ann = pd.read_csv(f'{D}/gpl.txt.gz', sep='\t', skiprows=27, usecols=['ID','Gene symbol'], dtype=str).dropna()
want = {p: s.split('///')[0] for p, s in zip(ann.ID, ann['Gene symbol']) if s.split('///')[0] in IFN}
rows = {}
with gzip.open(f'{D}/sm.txt.gz', 'rt') as f:
    for l in f:
        if l.startswith('"ID_REF"'): cols = [c.strip('"') for c in l.rstrip('\n').split('\t')[1:]]
        elif l[0] == '"':
            p, *v = l.rstrip('\n').split('\t'); p = p.strip('"')
            if p in want: rows[p] = np.array(v, float)
X = pd.DataFrame(rows, index=cols).T; X['gene'] = X.index.map(want)
X = X.assign(m=X.drop(columns='gene').mean(axis=1)).sort_values('m').groupby('gene').tail(1).set_index('gene').drop(columns='m').T
c = pd.read_csv(f'{D}/chars.csv'); c['days'] = pd.to_numeric(c.days_from_first_visit)
c['sledai'] = pd.to_numeric(c.sledai, errors='coerce')
def fit(A, y): return np.linalg.lstsq(A, y, rcond=None)[0]
res = {}
for drug in ['MMF', 'AZA']:
    s = c[c.drug_group == drug].sort_values('days').copy()
    Z = X.loc[s.gsm]; Z = (Z - Z.mean()) / Z.std(); s['IFN6'] = Z[IFN].mean(axis=1).values
    g = s.groupby('patient_id'); n_s = g.size()
    keep = n_s[n_s >= 2].index; s = s[s.patient_id.isin(keep)]; g = s.groupby('patient_id')
    P = pd.DataFrame({'first': g.IFN6.first(), 'last': g.IFN6.last(), 'dt': g.days.last() - g.days.first(),
                      'sl0': g.sledai.first(), 'sl1': g.sledai.last(), 'nvis': g.size(),
                      'R': (g.sri4_response.last() == 'YES').astype(int)})
    P['delta'] = P['last'] - P['first']; P['dsl'] = P.sl1 - P.sl0
    P.to_csv(f'{O}/q25_{drug}_patients.csv')
    def lopo(d, feats):
        e = []
        for i in d.index:
            tr = d.drop(i); A = np.column_stack([np.ones(len(tr))] + [tr[f] for f in feats])
            b = fit(A, tr['last'].values); a = np.r_[1, [d.loc[i, f] for f in feats]]
            e.append((d.loc[i, 'last'] - a @ b) ** 2)
        return float(np.mean(e))
    def block(d, tag):
        r = {'n': len(d), 'n_R': int(d.R.sum()), 'n_NR': int((1 - d.R).sum()), 'median_dt_days': float(d.dt.median())}
        if min(r['n_R'], r['n_NR']) < 3: r['note'] = 'too few in a class'; return r
        r['MSE_C0_carry'] = float(((d['last'] - d['first']) ** 2).mean())
        r['MSE_C1_setpoint'] = lopo(d, ['first']); r['MSE_C2_resp'] = lopo(d, ['first', 'R'])
        A = np.column_stack([np.ones(len(d)), d['first'], d.R]); r['c_full'] = float(fit(A, d['last'].values)[2])
        r['b_full_C1'] = float(fit(np.column_stack([np.ones(len(d)), d['first']]), d['last'].values)[1])
        r['median_delta_R'] = float(d[d.R == 1].delta.median()); r['median_delta_NR'] = float(d[d.R == 0].delta.median())
        r['MW_p'] = float(mannwhitneyu(d[d.R == 1].delta, d[d.R == 0].delta, alternative='two-sided').pvalue)
        r['ratio_C2_C1'] = r['MSE_C2_resp'] / r['MSE_C1_setpoint']
        r['PASS'] = bool(r['ratio_C2_C1'] <= 0.90 and r['c_full'] < 0 and r['MW_p'] < 0.05)
        dd = d.dropna(subset=['dsl'])
        A3 = np.column_stack([np.ones(len(dd)), dd['first'], dd.R, dd.dsl]); b3 = fit(A3, dd['last'].values)
        A2 = np.column_stack([np.ones(len(dd)), dd['first'], dd.R])
        r['H3'] = {'n': len(dd), 'c_without_dSLEDAI': float(fit(A2, dd['last'].values)[2]), 'c_with_dSLEDAI': float(b3[2]), 'd_dSLEDAI': float(b3[3])}
        return r
    res[drug] = {'n_patients_ge2_samples': len(P), 'n_patients_dropped_1_sample': int((n_s < 2).sum()), 'primary': block(P, drug)}
    Ps = P[P.dt >= 180]; res[drug]['sens_dt_ge180'] = block(Ps, drug + ' dt>=180')
    P2 = P[P.dt > 0].copy(); P2['delta'] = P2.delta / (P2.dt / 365.25)
    res[drug]['sens_per_year_slope_MW_p'] = float(mannwhitneyu(P2[P2.R == 1].delta, P2[P2.R == 0].delta).pvalue) if min(P2.R.sum(), (1-P2.R).sum()) >= 3 else None
    res[drug]['sens_per_year_median_R_NR'] = [float(P2[P2.R == 1].delta.median()), float(P2[P2.R == 0].delta.median())]
json.dump(res, open(f'{O}/q25_result.json', 'w'), indent=1); print(json.dumps(res, indent=1))
