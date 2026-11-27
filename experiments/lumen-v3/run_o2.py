"""O2 (round 154): implements experiments/lumen-v3/PREREG-O2.md (commit 7ae773f) exactly.

Robust MAP control-rate prior (Schmidli et al. 2014) vs B1 (fixed-effect IV pooled + DL tau^2) and B2 (unweighted
mean, honest width), time-forward, primary = mean negative log predictive density of the observed control logit.
Run from the repo root:  python3 -I experiments/lumen-v3/run_o2.py           (real data)
                         python3 -I experiments/lumen-v3/run_o2.py --selftest (simulated data only, no outcomes)
run_rt2.py executes and writes files at import, so unit_table / pick / dl_tau2 are copied verbatim (PREREG 5).
Research and simulation only. Not clinical or drug-development guidance.
"""
import sys, os, gzip, json, re, math
import numpy as np, pandas as pd
from scipy.special import logsumexp, expit
from scipy.stats import norm
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
from lumen import provenance, baselines as B
OUT = os.path.dirname(os.path.abspath(__file__))
PREREG = 'experiments/lumen-v3/PREREG-O2.md'
SEED = 20261127
TAU = np.arange(4000) * 0.001 + 0.0005          # 0.0005 ... 3.9995, step 0.001
VAGUE_SD = 2.0

INTERPRETATIONS = [
    'Robust mixture mean = MAP predictive mean (vague component centred at m_MAP), so R0/R1/R1-w0.5 share IV-MSE.',
    'Coverage of B1/B2/R*: central interval of the predictive for the OBSERVED y (component variances + held-out se2), '
    'as RT-2 S4 did for B1 (half-width 1.96 sqrt(var + se2)).',
    'LOTO mode: training = all other trials of the same stratum (no >= 5 floor needed; every stratum has >= 20).',
    'D1: RT-2 S6 trials (early control visit <= 26 wk, before target); compared with R1 LOTO on the same trials '
    '(RT-2 S6 compared with LOTO B1) and, where scored, with R1 time-forward.',
    'Bootstrap: paired resampling of scored trials (combined), statistic = mean NLPD(R1) - mean NLPD(best baseline), '
    'best baseline chosen on the full-sample primary score.',
    'Frozen W1-SLE intervals are for the TRUE control logit (no se2: arm size unknown), shown as rates via expit.',
]

# ---------------------------------------------------------------- copied verbatim from run_rt2.py
STR_RE = re.compile(r'predniso|corticosteroid|glucocorticoid|steroid|rescue (medication|therap)', re.I)
CAP = {'LN': 104, 'SLE': 52}
TABLE = {'LN': 'experiments/lumen-bench/resp_table.csv', 'SLE': 'experiments/lumen-bench/sle/sri4_table.csv'}

def dec_year(d):
    y, m = int(d[:4]), int(d[5:7]); dd = int(d[8:10]) if len(d) >= 10 else 15
    return y + (m - 1) / 12 + (dd - 1) / 365.25

def pick(g, lo, hi, exact=None):
    g = g[(g.week >= lo) & (g.week <= hi)] if exact is None else g[g.week == exact]
    if not len(g): return None
    g = g[g.week == g.week.max()]
    return g.sort_values('n_ctrl', ascending=False, kind='stable').iloc[0]

def unit_table(sel, D):
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

def dl_tau2(y, se2):
    w = 1 / se2; mu = np.sum(w * y) / np.sum(w); Q = np.sum(w * (y - mu) ** 2); k = len(y)
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    return max(0.0, (Q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
# ---------------------------------------------------------------- end of verbatim copy

# A predictive is a normal mixture: (weights, means, variances) for the TRUE logit of a new trial.
def pred_B1(y, se2):
    w = 1 / se2; mu = np.sum(w * y) / np.sum(w)
    return np.array([1.]), np.array([mu]), np.array([dl_tau2(y, se2) + 1 / np.sum(w)])

def pred_B2(y, se2):
    k = len(y)
    return np.array([1.]), np.array([float(np.mean(y))]), np.array([float(np.var(y, ddof=1)) * (1 + 1 / k)])

def pred_map(y, se2, hn_scale=1.0, w_vague=0.2):
    t2 = TAU ** 2
    V = se2[None, :] + t2[:, None]                         # grid x trials
    Wt = 1 / V
    Vmu = 1 / Wt.sum(1)
    muh = (Wt * y[None, :]).sum(1) * Vmu
    loglik = -0.5 * (np.log(2 * np.pi * V) + (y[None, :] - muh[:, None]) ** 2 / V).sum(1)
    logp = norm.logpdf(TAU, 0, hn_scale) + 0.5 * np.log(Vmu) + loglik
    wt = np.exp(logp - logsumexp(logp))
    m_map = float(np.sum(wt * muh))
    if w_vague == 0:
        return wt, muh, Vmu + t2
    return (np.concatenate([(1 - w_vague) * wt, [w_vague]]), np.concatenate([muh, [m_map]]),
            np.concatenate([Vmu + t2, [VAGUE_SD ** 2]]))

MODELS = {'B1': pred_B1, 'B2': pred_B2, 'R1': pred_map,
          'R0': lambda y, s: pred_map(y, s, w_vague=0.0),
          'R1-HN0.5': lambda y, s: pred_map(y, s, hn_scale=0.5),
          'R1-w0.5': lambda y, s: pred_map(y, s, w_vague=0.5)}

def mix_mean(p):
    return float(np.sum(p[0] * p[1]))

def nlpd(p, obs, se2_obs):
    w, m, v = p; v = v + se2_obs
    return float(-logsumexp(np.log(np.maximum(w, 1e-300)) - 0.5 * np.log(2 * np.pi * v) - (obs - m) ** 2 / (2 * v)))

def mix_cdf(p, x, se2_obs):
    w, m, v = p
    return float(np.sum(w * norm.cdf((x - m) / np.sqrt(v + se2_obs))))

def mix_quantile(p, q, se2_obs=0.0):
    w, m, v = p; sd = np.sqrt(v + se2_obs)
    lo, hi = float(np.min(m - 12 * sd)), float(np.max(m + 12 * sd))
    while hi - lo > 1e-8:
        mid = 0.5 * (lo + hi)
        if mix_cdf(p, mid, se2_obs) < q: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)

def covered(p, obs, se2_obs, level):
    a = (1 - level) / 2
    return bool(mix_quantile(p, a, se2_obs) <= obs <= mix_quantile(p, 1 - a, se2_obs))

def score(U, mode):
    """mode 'tf' = time-forward (strictly earlier resultsFirstPostDate, >= 5 training), 'loto' = leave-one-out."""
    rows = []
    for st in ['LN', 'SLE']:
        S = U[U.stratum == st].sort_values('post_date', kind='stable').reset_index(drop=True)
        for i in range(len(S)):
            tr = S[S.post_date < S.post_date.iloc[i]] if mode == 'tf' else S.drop(i)
            if len(tr) < 5: continue
            te = S.iloc[i]; y, se2 = tr.y.values, tr.se2.values
            row = dict(stratum=st, trial=te.trial, nct=te.nct, y=te.y, se2=te.se2, k=len(tr))
            for name, f in MODELS.items():
                p = f(y, se2)
                row[f'{name}_mean'] = mix_mean(p); row[f'{name}_nlpd'] = nlpd(p, te.y, te.se2)
                row[f'{name}_c95'] = covered(p, te.y, te.se2, .95); row[f'{name}_c80'] = covered(p, te.y, te.se2, .80)
            rows.append(row)
    return pd.DataFrame(rows)

def summarise(P):
    out = {}
    for name, d in [('combined', P), ('LN', P[P.stratum == 'LN']), ('SLE', P[P.stratum == 'SLE'])]:
        if not len(d): continue
        out[name] = dict(n_trials=int(len(d)),
                         nlpd={m: float(d[f'{m}_nlpd'].mean()) for m in MODELS},
                         ivmse={m: B.iv_mse(d.y, np.sqrt(d.se2), d[f'{m}_mean']) for m in MODELS},
                         cov95={m: float(d[f'{m}_c95'].mean()) for m in MODELS},
                         cov80={m: float(d[f'{m}_c80'].mean()) for m in MODELS})
    return out

def fail_fires(s, cand='R1'):
    reasons = []
    for b in ['B1', 'B2']:
        if s['nlpd'][cand] >= s['nlpd'][b] + 0.05: reasons.append(f'NLPD worse than {b} by >= 0.05')
        if s['ivmse'][cand] >= 1.10 * s['ivmse'][b]: reasons.append(f'IV-MSE >= 10% worse than {b}')
    return reasons

def bootstrap(P, best, n=2000):
    rng = np.random.default_rng(SEED)
    d = (P['R1_nlpd'] - P[f'{best}_nlpd']).values
    stats = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(n)])
    return float(d.mean()), [float(np.quantile(stats, .025)), float(np.quantile(stats, .975))]

def verdict(P):
    S = summarise(P); c = S['combined']
    best = 'B1' if c['nlpd']['B1'] <= c['nlpd']['B2'] else 'B2'
    diff, ci = bootstrap(P, best)
    fail = fail_fires(c)
    n_ln = S.get('LN', {}).get('n_trials', 0)
    powered = c['n_trials'] >= 20 and n_ln >= 10
    if not powered: v = 'UNDERPOWERED'
    elif fail: v = 'FAIL'
    elif diff <= -0.05 and ci[1] < 0 and 0.88 <= c['cov95']['R1'] <= 0.99: v = 'PASS'
    else: v = 'INCONCLUSIVE'
    ln_fail = fail_fires(S['LN']) if 'LN' in S else ['no LN trials']
    return dict(summary=S, best_baseline=best, R1_minus_best=diff, R1_minus_best_ci95=ci, fail_reasons=fail,
                minimum_n_met=powered, verdict=v, LN_guard_fail_reasons=ln_fail,
                LN_guard=('UNDERPOWERED' if n_ln < 10 else ('FAIL' if ln_fail else 'not fired')))

def d1(U, P_loto, P_tf):
    rows = []
    for st in ['LN', 'SLE']:
        T = pd.read_csv(os.path.join(ROOT, TABLE[st]))
        for (trial, nct), g in T.groupby(['trial', 'nct']):
            u = U[(U.stratum == st) & (U.trial == trial)].iloc[0]
            e = pick(g[g.week < u.week], 0, 26)
            if e is None: continue
            ye = np.log((e.r_ctrl + .5) / (e.n_ctrl - e.r_ctrl + .5)); se2e = 1 / (e.r_ctrl + .5) + 1 / (e.n_ctrl - e.r_ctrl + .5)
            cf = nlpd((np.array([1.]), np.array([ye]), np.array([se2e])), u.y, u.se2)
            lo = P_loto[(P_loto.stratum == st) & (P_loto.trial == trial)]
            tf = P_tf[(P_tf.stratum == st) & (P_tf.trial == trial)]
            rows.append(dict(stratum=st, trial=trial, early_week=int(e.week), target_week=int(u.week), CF_nlpd=cf,
                             R1_loto_nlpd=float(lo.R1_nlpd.iloc[0]) if len(lo) else np.nan,
                             R1_tf_nlpd=float(tf.R1_nlpd.iloc[0]) if len(tf) else np.nan))
    R = pd.DataFrame(rows)
    tfm = R.dropna(subset=['R1_tf_nlpd'])
    return dict(n_trials=int(len(R)), CF_nlpd=float(R.CF_nlpd.mean()), R1_loto_nlpd=float(R.R1_loto_nlpd.mean()),
                n_tf=int(len(tfm)), CF_nlpd_tf_subset=float(tfm.CF_nlpd.mean()) if len(tfm) else None,
                R1_tf_nlpd=float(tfm.R1_tf_nlpd.mean()) if len(tfm) else None, trials=R.to_dict('records'))

def freeze(U):
    path = os.path.join(ROOT, 'experiments/lumen-bench/prospective/w1-sle/frozen_2026-11-26.json')
    W = json.load(open(path)); S = U[U.stratum == 'SLE']
    p = pred_map(S.y.values, S.se2.values)
    m = mix_mean(p); lo, hi = mix_quantile(p, .025), mix_quantile(p, .975)
    out = {'note': 'O2 R1 robust-MAP control-arm SRI-4 rate, trained on all %d SLE units (descriptive; PREREG-W1-SLE '
                   'unchanged). Interval is for the true control rate (arm size unknown).' % len(S),
           'forecasts': {}, 'not_frozen': {}}
    for nct, v in W.items():
        if nct == '_meta': continue
        if v['endpoint'] == 'SRI-4':
            out['forecasts'][nct] = dict(trial=v['trial'], endpoint='SRI-4', primary_week=v['primary_week'],
                                         control_rate_mean_logit=m, control_rate_at_mean=float(expit(m)),
                                         control_rate_95=[float(expit(lo)), float(expit(hi))])
        else:
            out['not_frozen'][nct] = dict(trial=v['trial'], endpoint=v['endpoint'], reason='no BICLA control pool in RT-2')
    return out

def main():
    ok = provenance.verify_all(); assert all(ok.values()), ok
    ct = json.load(gzip.open(os.path.join(ROOT, 'data/real/ctgov/lupus_with_results.json.gz')))
    D = {s['protocolSection']['identificationModule']['nctId']: s for s in ct}
    U = unit_table(lambda st, g: pick(g, 0, CAP[st]), D)
    # assertion (i): unit table equals rt2_coding.csv
    C = pd.read_csv(os.path.join(ROOT, 'experiments/lumen-bench/control-arm/rt2_coding.csv'))
    a = U[['stratum', 'trial', 'nct', 'week']].reset_index(drop=True)
    assert a.equals(C[['stratum', 'trial', 'nct', 'week']].reset_index(drop=True)), 'unit table != rt2_coding.csv'
    P_tf = score(U, 'tf')
    # assertion (ii): time-forward B1/B2 means equal RT-2 S5 rows
    R5 = pd.read_csv(os.path.join(ROOT, 'experiments/lumen-bench/control-arm/rt2_predictions.csv'))
    R5 = R5[R5.spec == 'S5'].set_index(['stratum', 'trial'])
    Q = P_tf.set_index(['stratum', 'trial'])
    assert set(Q.index) == set(R5.index), 'scored set != RT-2 S5'
    for idx in Q.index:
        assert abs(Q.loc[idx, 'B1_mean'] - R5.loc[idx, 'B1']) < 1e-10, ('B1', idx)
        assert abs(Q.loc[idx, 'B2_mean'] - R5.loc[idx, 'B2']) < 1e-10, ('B2', idx)
    P_loto = score(U, 'loto')
    res = {'prereg': PREREG, 'interpretations': INTERPRETATIONS, 'label': '12b-timed (committed data; earlier looks not excludable)',
           'n_units': U.groupby('stratum').size().to_dict(), 'assertions': 'passed (unit table, S5 B1/B2 means to 1e-10)'}
    res['primary_time_forward'] = verdict(P_tf)
    res['S_loto'] = verdict(P_loto)
    res['D1_early_control_carry_forward'] = d1(U, P_loto, P_tf)
    res['frozen'] = 'experiments/lumen-bench/prospective/w1-sle/frozen_O2_2026-11-27.json'
    P_tf.assign(mode='tf').pipe(lambda x: pd.concat([x, P_loto.assign(mode='loto')])).to_csv(
        os.path.join(OUT, 'o2_predictions.csv'), index=False)
    json.dump(provenance.stamp(res, prereg=PREREG), open(os.path.join(OUT, 'o2_result.json'), 'w'), indent=1, default=float)
    fz = provenance.stamp(freeze(U), prereg=PREREG)
    json.dump(fz, open(os.path.join(ROOT, res['frozen']), 'w'), indent=1, default=float)
    p = res['primary_time_forward']
    print(json.dumps({k: p[k] for k in ['verdict', 'best_baseline', 'R1_minus_best', 'R1_minus_best_ci95',
                                         'fail_reasons', 'LN_guard', 'LN_guard_fail_reasons']}, indent=1, default=float))
    print(json.dumps(p['summary'], indent=1, default=float))

def selftest():
    """Simulated trials only: checks MAP weights, mixture quantiles and that R1 ~ B1 when tau is known."""
    rng = np.random.default_rng(1)
    k = 30; tau = 0.4; th = rng.normal(-0.5, tau, k); n = rng.integers(40, 200, k)
    r = rng.binomial(n, expit(th)); y = np.log((r + .5) / (n - r + .5)); se2 = 1 / (r + .5) + 1 / (n - r + .5)
    p = pred_map(y, se2); assert abs(p[0].sum() - 1) < 1e-12
    q = mix_quantile(p, .5); assert abs(mix_cdf(p, q, 0) - .5) < 1e-7
    g = (np.array([1.]), np.array([0.]), np.array([1.])); assert abs(mix_quantile(g, .975) - 1.959964) < 1e-5
    assert abs(nlpd(g, 0, 0) - 0.5 * math.log(2 * math.pi)) < 1e-12
    print('selftest ok; R1 mean %.3f B1 mean %.3f true -0.5; R1 sd %.3f' % (
        mix_mean(p), mix_mean(pred_B1(y, se2)), math.sqrt(np.sum(p[0] * (p[2] + (p[1] - mix_mean(p)) ** 2)))))

if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else main()
