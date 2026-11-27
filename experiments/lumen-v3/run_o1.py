"""O1 (round 151): implements experiments/lumen-v3/PREREG-O1.md (commit 4e8975a) exactly.

Run from the repo root:  python3 -I experiments/lumen-v3/run_o1.py
Imports run_rt1.py unchanged (tables, DATE, CFG, late_arms, score_stratum). Research and simulation only.
"""
import sys, os, json, math
import numpy as np, pandas as pd, scipy.stats as S
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'experiments', 'lumen-bench', 'timeforward'))
from lumen import provenance, candidates as C, baselines as B
ok = provenance.verify_all(); assert all(ok.values()), ok
import run_rt1 as R1   # unchanged
OUT = os.path.dirname(os.path.abspath(__file__))
SEED = 20261127

def train_sets(R, st, mode):
    """Same training rule as run_rt1.score_stratum: yields (trial, training late arms) for scored trials."""
    trials = R[['trial', 'nct']].drop_duplicates()
    trials = trials.assign(date=trials.nct.map(R1.DATE))
    for _, tr in trials.iterrows():
        tt = trials[trials.date < tr.date].trial if mode == 'tf' else trials[trials.trial != tr.trial].trial
        train = R[R.trial.isin(tt)]
        if train.trial.nunique() < 5: continue
        yield tr.trial, R1.late_arms(train, st)

def candidate_preds(P):
    y, se = P.lor.values, P.se.values
    out = {}
    m, s, df = C.predictive_frequentist(y, se); out['O1a'] = dict(kind='t', mean=m, sd=s, df=df)
    m, s = C.predictive_bayes(y, se, half_normal_tau=0.5); out['O1b'] = dict(kind='n', mean=m, sd=s)
    m, s = C.predictive_bayes(y, se); out['O1c'] = dict(kind='n', mean=m, sd=s)
    return out

def nlpd(kind, obs, se_obs, p):
    if kind == 'n':
        v = p['sd'] ** 2 + se_obs ** 2
        return 0.5 * math.log(2 * math.pi * v) + (obs - p['mean']) ** 2 / (2 * v)
    if kind == 't':
        sc = math.sqrt(p['sd'] ** 2 + se_obs ** 2)
        return -S.t.logpdf(obs, p['df'], loc=p['mean'], scale=sc)
    if kind == 'mix':
        dens = sum(w * S.norm.pdf(obs, m, math.sqrt(s ** 2 + se_obs ** 2)) for w, m, s in p['components'])
        return -math.log(dens)

def covered(kind, obs, se_obs, p, level):
    if kind == 't':
        q = S.t.ppf(0.5 + level / 2, p['df']); return abs(obs - p['mean']) <= q * math.sqrt(p['sd'] ** 2 + se_obs ** 2)
    if kind == 'n':
        q = S.norm.ppf(0.5 + level / 2); return abs(obs - p['mean']) <= q * math.sqrt(p['sd'] ** 2 + se_obs ** 2)
    if kind == 'mix':   # central interval of the mixture by bisection on its CDF
        sd = lambda s: math.sqrt(s ** 2 + se_obs ** 2)
        F = lambda x: sum(w * S.norm.cdf(x, m, sd(s)) for w, m, s in p['components'])
        def q(a):
            lo, hi = -20.0, 20.0
            for _ in range(80):
                mid = (lo + hi) / 2; lo, hi = (mid, hi) if F(mid) < a else (lo, mid)
            return (lo + hi) / 2
        return q(0.5 - level / 2) <= obs <= q(0.5 + level / 2)

MODELS = ['M0', 'M1', 'O1a', 'O1b', 'O1c', 'O1d']

def build(mode):
    rows = []
    for st, R in [('LN', R1.LN), ('SLE', R1.SLE)]:
        sc = pd.DataFrame(R1.score_stratum(R, st, mode)); sc = sc[sc.subset == 'all']
        cand = {t: candidate_preds(P) for t, P in train_sets(R, st, mode)}
        for _, r in sc.iterrows():
            sd1 = math.sqrt(r.v_M1)
            preds = dict(M0=dict(kind='n', mean=0.0, sd=sd1), M1=dict(kind='n', mean=r.M1, sd=sd1), **cand[r.trial])
            mx = C.robust_mixture(r.M1, sd1, 0.2, 1.0); preds['O1d'] = dict(kind='mix', mean=mx['mean'], sd=mx['sd'], components=mx['components'])
            row = dict(mode=mode, stratum=st, trial=r.trial, arm=r.arm, obs=r.obs, se_obs=r.se_obs)
            for m in MODELS:
                p = preds[m]
                row[m + '_mean'] = p['mean']; row[m + '_sd'] = p['sd']; row[m + '_nlpd'] = nlpd(p['kind'], r.obs, r.se_obs, p)
                row[m + '_cov95'] = covered(p['kind'], r.obs, r.se_obs, p, 0.95)
                row[m + '_cov80'] = covered(p['kind'], r.obs, r.se_obs, p, 0.80)
            row['O1a_df'] = preds['O1a']['df']
            rows.append(row)
    return pd.DataFrame(rows)

def summ(D, label):
    T = D.groupby(['stratum', 'trial'])[[m + '_nlpd' for m in MODELS]].mean()   # arms averaged within trial
    out = dict(label=label, n_trials=int(len(T)), n_arms=int(len(D)))
    for m in MODELS:
        out[m] = dict(logscore=round(float(T[m + '_nlpd'].mean()), 4),
                      ivmse=round(B.iv_mse(D.obs, D.se_obs, D[m + '_mean']), 4),
                      cov95=round(float(D[m + '_cov95'].mean()), 3), cov80=round(float(D[m + '_cov80'].mean()), 3),
                      mean_sd=round(float(D[m + '_sd'].mean()), 3))
    rng = np.random.default_rng(SEED); d = (T['O1a_nlpd'] - T['M1_nlpd']).values
    bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(2000)]
    out['O1a_minus_M1'] = dict(diff=round(float(d.mean()), 4), ci95=[round(float(np.percentile(bs, 2.5)), 4), round(float(np.percentile(bs, 97.5)), 4)])
    out['O1a_minus_M0'] = round(float((T['O1a_nlpd'] - T['M0_nlpd']).mean()), 4)
    return out

def verdict(s):
    a, m1, m0 = s['O1a'], s['M1'], s['M0']
    if s['n_trials'] < 20: return 'UNDERPOWERED'
    if a['logscore'] - m0['logscore'] >= 0.05 or a['logscore'] - m1['logscore'] >= 0.05 \
            or a['ivmse'] >= 1.10 * m0['ivmse'] or a['ivmse'] >= 1.10 * m1['ivmse']:
        return 'FAIL'
    if m1['logscore'] - a['logscore'] >= 0.05 and s['O1a_minus_M1']['ci95'][1] < 0 and 0.88 <= a['cov95'] <= 0.99:
        return 'PASS'
    return 'INCONCLUSIVE'

def freeze_w1():
    """Descriptive prospective record: O1a/O1b for the 8 frozen W1-SLE trials, same training rows as lumen_forecast.forecast."""
    sys.path.insert(0, os.path.join(ROOT, 'experiments', 'lumen-bench'))
    import lumen_forecast as L
    F = json.load(open(os.path.join(ROOT, 'experiments/lumen-bench/prospective/w1-sle/frozen_2026-11-26.json')))
    out = {}
    for nct, f in F.items():
        if nct == '_meta': continue
        R = L.load(os.path.join(ROOT, 'experiments/lumen-bench', f['table']))
        P = R[R.week <= 104].sort_values('week').groupby(['trial', 'arm']).tail(1)
        m, s, df = C.predictive_frequentist(P.lor.values, P.se.values); q = S.t.ppf(0.975, df)
        mb, sb = C.predictive_bayes(P.lor.values, P.se.values, half_normal_tau=0.5)
        out[nct] = dict(trial=f['trial'], table=f['table'], M1_frozen=[f['M1_pooled']['prior_logOR'], f['M1_pooled']['prior_sd']],
                        O1a=dict(mean=round(m, 3), sd=round(s, 3), df=df, logOR_95=[round(m - q * s, 3), round(m + q * s, 3)]),
                        O1b=dict(mean=round(mb, 3), sd=round(sb, 3), logOR_95=[round(mb - 1.96 * sb, 3), round(mb + 1.96 * sb, 3)]))
    return out

def main():
    tf, lo = build('tf'), build('loto')
    res = dict(prereg='experiments/lumen-v3/PREREG-O1.md', md5_verify=ok)
    res['primary'] = summ(tf, 'time-forward, all subset, LN+SLE')
    res['verdict'] = verdict(res['primary'])
    res['S_loto'] = summ(lo, 'LOTO, all subset, LN+SLE'); res['S_loto']['verdict_by_rule_descriptive'] = verdict(res['S_loto'])
    res['S_LN'] = summ(tf[tf.stratum == 'LN'], 'time-forward LN only')
    res['S_SLE'] = summ(tf[tf.stratum == 'SLE'], 'time-forward SLE only')
    res['O1a_df_range'] = [int(tf.O1a_df.min()), int(tf.O1a_df.max())]
    provenance.write_once_text(os.path.join(OUT, 'o1_predictions.csv'), pd.concat([tf, lo]).to_csv(index=False))
    res = provenance.stamp(res, prereg='experiments/lumen-v3/PREREG-O1.md')
    provenance.write_once(os.path.join(OUT, 'o1_result.json'), res, indent=1, default=str)
    fz = provenance.stamp(dict(note='Descriptive prospective record per PREREG-O1 s10; no outcomes read; PREREG-W1-SLE unchanged',
                                forecasts=freeze_w1()), prereg='experiments/lumen-v3/PREREG-O1.md')
    provenance.write_frozen(os.path.join(ROOT, 'experiments/lumen-bench/prospective/w1-sle/frozen_O1_2026-11-27.json'), fz, indent=1)
    print(json.dumps(res, indent=1, default=str))

if __name__ == '__main__':
    main()
