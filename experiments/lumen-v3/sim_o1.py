"""Synthetic coverage check of predictive intervals for a NEW trial (no real outcome data)."""
import sys, numpy as np, math, scipy.stats as S
import os; sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from lumen import candidates as c
rng = np.random.default_rng(0)
def run(k, tau, reps):
    hits = {n: [] for n in ['M1_shipped_DL_fixed_z', 'REML_RE_t', 'Bayes_HN0.5', 'Bayes_Turner_general']}
    for _ in range(reps):
        se = rng.uniform(0.15, 0.6, k + 1); th = 0.3 + rng.normal(0, tau, k + 1); y = th + rng.normal(0, se)
        yt, st, new, sn = y[:k], se[:k], y[k], se[k]
        t2 = c.tau2_dl(yt, st); w = 1/st**2; mu = (w*yt).sum()/w.sum()
        pr = {'M1_shipped_DL_fixed_z': (mu, math.sqrt(t2 + 1/w.sum()), 1.96)}
        m, s, df = c.predictive_frequentist(yt, st); pr['REML_RE_t'] = (m, s, S.t.ppf(.975, df))
        pr['Bayes_HN0.5'] = (*c.predictive_bayes(yt, st, half_normal_tau=0.5, grid=100), 1.96)
        pr['Bayes_Turner_general'] = (*c.predictive_bayes(yt, st, grid=100), 1.96)
        for n, (m, s, z) in pr.items(): hits[n].append(abs(new - m) <= z*math.sqrt(s**2 + sn**2))
    return {n: round(float(np.mean(v)), 3) for n, v in hits.items()}
for k in (5, 10, 20):
    for tau in (0.1, 0.3, 0.6):
        print(k, tau, run(k, tau, 300), flush=True)
