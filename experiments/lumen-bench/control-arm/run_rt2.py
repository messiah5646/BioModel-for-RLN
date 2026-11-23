"""RT-20261123-2 (round 128): control-arm response model for LN and SLE trials.

Implements experiments/lumen-bench/control-arm/PREREG-RT2.md (commit d5b5f84) exactly; interpretations of
under-specified points are listed in INTERPRETATIONS and copied into the result JSON.
Run from the repo root:  python3 -I experiments/lumen-bench/control-arm/run_rt2.py
Research and simulation only. Not clinical guidance.
"""
import sys, os, gzip, json, re, hashlib
import numpy as np, pandas as pd
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'provenance'))
from stamp import stamp

OUT = os.path.join(ROOT, 'experiments', 'lumen-bench', 'control-arm')
PREREG = 'experiments/lumen-bench/control-arm/PREREG-RT2.md'
MD5 = {'data/real/ctgov/lupus_with_results.json.gz': '1b31cd52ccc62dc6124cd611342ba5d0',
       'experiments/lumen-bench/resp_table.csv': '3ce67bdb51c6597bce6666d914618778',
       'experiments/lumen-bench/sle/sri4_table.csv': 'b37dd7cd72c9e912ba079088532c3765'}
for p, h in MD5.items():
    assert hashlib.md5(open(os.path.join(ROOT, p), 'rb').read()).hexdigest() == h, f'md5 mismatch {p}'

INTERPRETATIONS = [
    'Measure location: measure_ctrl in the committed tables is the outcome index chosen by build_resp_table.py / '
    'build_sle_table.py; both builders were re-run this run and reproduce the committed CSVs byte-identically '
    '(md5 3ce67bdb.., b37dd7cd..), so the index IS the re-run selection.',
    'S3 visit-matched: unit re-selected as the latest control row with 48<=week<=56 (LN) or week==52 (SLE); '
    'trials without such a row are dropped from S3 only.',
    'S4 M1 interval: pred +/- 1.96 sqrt(tau2 + x\'(X\'WX)^-1 x + se_obs^2), tau2 = the DL tau2 used for the M1 weights.',
    'S5 time-forward: train on trials of the same stratum with strictly earlier resultsFirstPostDate; score only '
    'trials with >= 5 such earlier trials; same pass/FAIL thresholds, reported as secondary.',
    'S6: early = latest control row with week <= 26 and week < target week; error vs observed target logit; '
    'compared with the LOTO B1 prediction on the same subset (descriptive).',
]

ct = json.load(gzip.open(os.path.join(ROOT, 'data/real/ctgov/lupus_with_results.json.gz')))
D = {s['protocolSection']['identificationModule']['nctId']: s for s in ct}

def dec_year(d):
    y, m = int(d[:4]), int(d[5:7]); dd = int(d[8:10]) if len(d) >= 10 else 15
    return y + (m - 1) / 12 + (dd - 1) / 365.25

STR_RE = re.compile(r'predniso|corticosteroid|glucocorticoid|steroid|rescue (medication|therap)', re.I)
CAP = {'LN': 104, 'SLE': 52}
TABLE = {'LN': 'experiments/lumen-bench/resp_table.csv', 'SLE': 'experiments/lumen-bench/sle/sri4_table.csv'}

def pick(g, lo, hi, exact=None):
    g = g[(g.week >= lo) & (g.week <= hi)] if exact is None else g[g.week == exact]
    if not len(g): return None
    g = g[g.week == g.week.max()]
    return g.sort_values('n_ctrl', ascending=False, kind='stable').iloc[0]

def unit_table(sel):
    rows = []
    for st in ['LN', 'SLE']:
        T = pd.read_csv(os.path.join(ROOT, TABLE[st]))
        for (trial, nct), g in T.groupby(['trial', 'nct'], sort=True):
            r = sel(st, g)
            if r is None: continue
            om = D[nct]['resultsSection']['outcomeMeasuresModule']['outcomeMeasures'][int(r.measure_ctrl)]
            title = om.get('title', ''); desc = om.get('description', '') or ''
            date = D[nct]['protocolSection']['statusModule']['resultsFirstPostDateStruct']['date']
            rc, nc = float(r.r_ctrl), float(r.n_ctrl)
            rows.append(dict(stratum=st, trial=trial, nct=nct, week=int(r.week), r_ctrl=rc, n_ctrl=nc,
                             y=np.log((rc + .5) / (nc - rc + .5)), se2=1 / (rc + .5) + 1 / (nc - rc + .5),
                             year=dec_year(date), post_date=date, title=title, description=desc[:300],
                             STRINGENT=int(bool(STR_RE.search(title + ' ' + desc)))))
    return pd.DataFrame(rows)

U = unit_table(lambda st, g: pick(g, 0, CAP[st]))
# coding table written BEFORE any model is fitted (PREREG section 4)
U[['stratum', 'trial', 'nct', 'week', 'title', 'description', 'STRINGENT']].to_csv(
    os.path.join(OUT, 'rt2_coding.csv'), index=False)

def dl_tau2(y, se2):
    w = 1 / se2; mu = np.sum(w * y) / np.sum(w); Q = np.sum(w * (y - mu) ** 2); k = len(y)
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    return max(0.0, (Q - (k - 1)) / c) if k > 1 and c > 0 else 0.0

def fit_predict(tr, te, covs):
    """returns dict model -> (pred, var_for_interval), flags"""
    y, se2 = tr.y.values, tr.se2.values
    t2 = dl_tau2(y, se2); w = 1 / se2
    mu = np.sum(w * y) / np.sum(w); semu2 = 1 / np.sum(w)
    out = {'B1': (mu, t2 + semu2), 'B2': (float(np.mean(y)), np.nan)}
    use = [c for c in covs if not (c == 'STRINGENT' and tr.STRINGENT.nunique() < 2)]
    flag = len(use) < len(covs)
    ym = tr.year.mean()
    def X(d):
        cols = [np.ones(len(d))]
        for c in use: cols.append(d.year.values - ym if c == 'YEAR' else d.STRINGENT.values.astype(float))
        return np.column_stack(cols)
    Xt = X(tr); W = 1 / (se2 + t2)
    A = Xt.T @ (W[:, None] * Xt); beta = np.linalg.solve(A, Xt.T @ (W * y)); Ainv = np.linalg.inv(A)
    xe = X(te)[0]
    out['M'] = (float(xe @ beta), t2 + float(xe @ Ainv @ xe))
    return out, flag, use

def ivmse(d, col):
    w = 1 / d.se2
    return float(np.sum(w * (d[col] - d.y) ** 2) / np.sum(w))

def verdict(m, base):
    if m <= 0.90 * base: return 'PASS'
    if m >= 1.10 * base: return 'FAIL'
    return 'INCONCLUSIVE'

def loto(U, covs, tag):
    rows = []
    for st in ['LN', 'SLE']:
        S = U[U.stratum == st].reset_index(drop=True)
        for i in range(len(S)):
            te = S.iloc[[i]]; tr = S.drop(i)
            o, flag, use = fit_predict(tr, te, covs)
            rows.append(dict(stratum=st, trial=te.trial.iloc[0], y=te.y.iloc[0], se2=te.se2.iloc[0],
                             B1=o['B1'][0], B1_var=o['B1'][1], B2=o['B2'][0], M1=o['M'][0], M1_var=o['M'][1],
                             str_dropped=flag, spec=tag))
    return pd.DataFrame(rows)

def summarise(P):
    res = {}
    for name, d in [('combined', P), ('LN', P[P.stratum == 'LN']), ('SLE', P[P.stratum == 'SLE'])]:
        if not len(d): continue
        m = {k: ivmse(d, k) for k in ['B1', 'B2', 'M1']}
        base = min(m['B1'], m['B2'])
        res[name] = dict(n_trials=int(len(d)), ivmse=m, strongest_baseline='B1' if m['B1'] <= m['B2'] else 'B2',
                         ratio_M1_to_best=m['M1'] / base, verdict=verdict(m['M1'], base),
                         folds_str_dropped=int(d.str_dropped.sum()))
    return res

def coverage(P):
    out = {}
    for name, d in [('combined', P), ('LN', P[P.stratum == 'LN']), ('SLE', P[P.stratum == 'SLE'])]:
        c = {}
        for k in ['B1', 'M1']:
            half = 1.96 * np.sqrt(d[f'{k}_var'] + d.se2)
            c[k] = float(np.mean(np.abs(d.y - d[k]) <= half))
        out[name] = c
    return out

result = {'prereg': PREREG, 'interpretations': INTERPRETATIONS}
n = U.groupby('stratum').size().to_dict()
result['n_units'] = n
result['stringent_counts'] = U.groupby('stratum').STRINGENT.sum().astype(int).to_dict()
result['minimum_n_met'] = bool(len(U) >= 40 and n.get('LN', 0) >= 15 and n.get('SLE', 0) >= 25)

P = loto(U, ['YEAR', 'STRINGENT'], 'primary')
prim = summarise(P)
result['primary'] = prim
result['primary_verdict'] = prim['combined']['verdict'] if result['minimum_n_met'] else 'UNDERPOWERED'
result['S2_LN_guard'] = prim['LN']['verdict']
result['S1'] = {'YEAR_only': summarise(loto(U, ['YEAR'], 'year')),
                'STRINGENT_only': summarise(loto(U, ['STRINGENT'], 'str'))}
U3 = unit_table(lambda st, g: pick(g, 48, 56) if st == 'LN' else pick(g, 0, 0, exact=52))
result['S3_visit_matched'] = summarise(loto(U3, ['YEAR', 'STRINGENT'], 'S3'))
result['S4_coverage'] = coverage(P)

# S5 time-forward
rows = []
for st in ['LN', 'SLE']:
    S = U[U.stratum == st].sort_values('post_date').reset_index(drop=True)
    for i in range(len(S)):
        tr = S[S.post_date < S.post_date.iloc[i]]
        if len(tr) < 5: continue
        te = S.iloc[[i]]; o, flag, _ = fit_predict(tr, te, ['YEAR', 'STRINGENT'])
        rows.append(dict(stratum=st, trial=te.trial.iloc[0], y=te.y.iloc[0], se2=te.se2.iloc[0], B1=o['B1'][0],
                         B1_var=o['B1'][1], B2=o['B2'][0], M1=o['M'][0], M1_var=o['M'][1], str_dropped=flag, spec='S5'))
P5 = pd.DataFrame(rows)
result['S5_time_forward'] = summarise(P5)

# S6 drift of control rate with visit length
rows = []
for st in ['LN', 'SLE']:
    T = pd.read_csv(os.path.join(ROOT, TABLE[st]))
    for (trial, nct), g in T.groupby(['trial', 'nct']):
        u = U[(U.stratum == st) & (U.trial == trial)].iloc[0]
        e = pick(g[g.week < u.week], 0, 26)
        if e is None: continue
        ye = np.log((e.r_ctrl + .5) / (e.n_ctrl - e.r_ctrl + .5))
        pb = P[(P.stratum == st) & (P.trial == trial)].B1.iloc[0]
        rows.append(dict(stratum=st, trial=trial, early_week=int(e.week), target_week=int(u.week), y=u.y, se2=u.se2,
                         CF=ye, B1=pb))
S6 = pd.DataFrame(rows)
result['S6_drift'] = {'n_trials': int(len(S6)),
                      'ivmse_carry_forward': ivmse(S6, 'CF') if len(S6) else None,
                      'ivmse_B1_same_subset': ivmse(S6, 'B1') if len(S6) else None,
                      'mean_logit_change_early_to_target': float((S6.y - S6.CF).mean()) if len(S6) else None,
                      'trials': S6[['stratum', 'trial', 'early_week', 'target_week']].to_dict('records')}

# descriptive: raw control rates by stratum
result['control_rate_summary'] = {st: dict(median=float((d.r_ctrl / d.n_ctrl).median()),
                                           min=float((d.r_ctrl / d.n_ctrl).min()), max=float((d.r_ctrl / d.n_ctrl).max()))
                                  for st, d in U.groupby('stratum')}

pd.concat([P, P5]).to_csv(os.path.join(OUT, 'rt2_predictions.csv'), index=False)
json.dump(stamp(result, prereg=PREREG), open(os.path.join(OUT, 'rt2_result.json'), 'w'), indent=1, default=float)
print(json.dumps({k: result[k] for k in ['n_units', 'stringent_counts', 'minimum_n_met', 'primary_verdict',
                                          'S2_LN_guard', 'primary']}, indent=1, default=float))
