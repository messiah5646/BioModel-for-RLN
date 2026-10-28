"""
Round 43/44: the placebo falsifier and the variability-power squeeze.

Round 43 found that the model reproduces a realistic control-arm response rate only
at between-patient CV ~0.70, above the 0.50 pre-declared plausible -- so it survives
conditionally. It also found, unpredicted, that the treated-control gap SHRINKS as
variability rises.

Round 44 attacks that. The alternative explanation is that CRR is a threshold on a
continuous quantity, so any two arms must converge as spread grows -- arithmetic, not
biology. The discriminator is whether a CONTINUOUS effect size (Cohen's d on
log UPCR_52) erodes too. It does, more significantly than the binary one, so the
squeeze is not a threshold artifact.

Population runs use rtol=1e-7 (verified: 2e-9 relative vs the tight solve) and
checkpoint per cell -- round 43 v1 lost an unsaved run to its own timeout.
"""
import json
import os
import numpy as np
from scipy.stats import norm
from model import make_params, simulate, observe, treated

JIT = ['p0', 'd_Ab', 'k_clear', 'k_a', 'g_l', 'd_Pl']   # top 6 BY SENSITIVITY
TG = np.array([0., 364.])
RESULTS = 'results/r44_squeeze.json'


def base_params():
    p = make_params(phi=1.0)
    p['k_dam'] *= 10
    p['k_tox'] *= 10
    return p


def one_patient(p):
    """AURORA-1-like CRR at week 52, plus the continuous endpoint."""
    try:
        Y = simulate(p, TG, rtol=1e-7, atol=1e-9)
    except RuntimeError:
        return False, np.nan, False
    O = np.array([observe(y, p) for y in Y])
    u, g, gb = O[-1, 0], O[-1, 1], O[0, 1]
    if g <= 0.05 * p['G0'] + 1e-9:          # eGFR clamp -> non-physiological
        return False, np.nan, False
    crr = bool(u <= 0.5 and (g >= 60.0 or g >= 0.8 * gb))
    return crr, float(np.log(max(u, 1e-12))), True


def arm(cv, n, treat, rng, base):
    sd = np.sqrt(np.log(1 + cv ** 2))
    hits, logu = [], []
    for _ in range(n):
        p = dict(base)
        for nm in JIT:
            p[nm] = base[nm] * np.exp(rng.normal(0, sd) - 0.5 * sd * sd)
        if treat:
            p = treated(p)
        c, lu, ok = one_patient(p)
        if ok:
            hits.append(c); logu.append(lu)
    return np.array(hits, dtype=float), np.array(logu)


def required_n(p1, p2, alpha=0.05, power=0.8):
    if abs(p1 - p2) < 1e-9:
        return np.inf
    z = norm.ppf(1 - alpha / 2) + norm.ppf(power)
    return z ** 2 * (p1 * (1 - p1) + p2 * (1 - p2)) / (p1 - p2) ** 2


def main(cvs=(0.40, 0.50, 0.60, 0.70, 0.80, 1.00), seeds=range(8), n=600):
    base = base_params()
    res = json.load(open(RESULTS)) if os.path.exists(RESULTS) else []
    done = {(r['cv'], r['seed']) for r in res}
    for cv in cvs:
        for seed in seeds:
            if (cv, seed) in done:
                continue
            rng = np.random.default_rng(1000 + seed)    # paired: same seed both arms
            hc, lc = arm(cv, n, False, rng, base)
            ht, lt = arm(cv, n, True, rng, base)
            pc, pt = hc.mean(), ht.mean()
            sp = np.sqrt((lc.var(ddof=1) + lt.var(ddof=1)) / 2)
            d = (lc.mean() - lt.mean()) / sp if sp > 0 else 0.0
            res.append(dict(cv=cv, seed=seed, control=float(pc), treated=float(pt),
                            diff_pp=float((pt - pc) * 100), cohen_d=float(d),
                            req_n=float(required_n(pc, pt)), nc=len(hc), nt=len(ht)))
            json.dump(res, open(RESULTS, 'w'), indent=1)    # checkpoint per cell
            print(f"cv={cv:.2f} s={seed} ctrl={pc*100:5.1f}% trt={pt*100:5.1f}% "
                  f"d={d:5.3f}", flush=True)


if __name__ == '__main__':
    main()
