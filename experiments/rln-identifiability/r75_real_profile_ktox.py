"""Round 75 -- prereg: project doc claude/prereg-round75.md (written before any k_tox-free fit).
Profile likelihood of k_tox on real AURORA 1 arm means: placebo only (P), treated only (T), both (PT)."""
import csv, json, sys, numpy as np
from scipy.optimize import least_squares
from model import make_params, simulate, observe
def _load(path):
    with open(path) as f: return list(csv.DictReader([l for l in f if not l.lstrip().startswith('#')]))
D = _load('../../data/real/aurora1_per_visit.csv')
BASE = {'treated': (4.14, 78.3, 179), 'placebo': (3.87, 77.4, 178)}
W = [2, 4, 8, 12, 16, 20, 24, 30, 36, 42, 48, 52]; TG = np.array([0.0] + [w * 7.0 for w in W])
OBS = {a: {int(r['week']): {k: float(v) for k, v in r.items() if k != 'arm'} for r in D if r['arm'] == a} for a in BASE}
def arm_curve(arm, s_inf, kdam, ktox, cni=0.0, leak=0.0, h=0.0):
    u0, g0, _ = BASE[arm]; p = make_params(0.0); p['p0'] = u0; p['G0'] = g0 / (1 - p['c_inf']); p['s_inf'] = s_inf * (1 - cni)
    p['k_dam'], p['k_tox'], p['k_prog'] = kdam, ktox, 0.0; p['h_drug'] = h; p['leak_red'] = leak; p['drug_on'], p['drug_off'] = 1.0, 1e9
    Y = simulate(p, TG, rtol=1e-8, atol=1e-10); O = np.array([observe(y, p, t) for t, y in zip(TG, Y)])
    return O[1:, :2] - O[0, :2]
sg = lambda z: 1 / (1 + np.exp(-z))
def resid(x, lk, arms):
    s_inf, kdam, ktox = np.exp(min(x[0], 0.0)), np.exp(x[1]), 10.0 ** lk; r = []
    try:
        for arm in arms:
            C = arm_curve(arm, s_inf, kdam, ktox, *((sg(x[2]), sg(x[3]), 0.15 * sg(x[4])) if arm == 'treated' else ()))
            n = BASE[arm][2]
            for i, w in enumerate(W):
                o = OBS[arm][w]; r += [(C[i, 0] - o['upcr_change']) / (o['upcr_sd'] / np.sqrt(n)), (C[i, 1] - o['egfr_change']) / (o['egfr_sd'] / np.sqrt(n))]
    except Exception:
        return np.full(len(arms) * len(W) * 2, 1e3)
    return np.array(r)
GRID = np.linspace(-7, -3, 9) if len(sys.argv) > 1 else np.linspace(-7, -1, 13); out = {}  # truncated grid: fits at >= -2 hang (deviation, declared in round doc)
for tag, arms in ((('PT', ['placebo', 'treated']), ('T', ['treated']), ('P', ['placebo'])) if len(sys.argv) > 1 else (('P', ['placebo']), ('T', ['treated']), ('PT', ['placebo', 'treated']))):
    npar = 2 + (3 if 'treated' in arms else 0); prof = []; warm = None; rng = np.random.default_rng(0)
    for lk in GRID:
        best = None; starts = ([warm] if warm is not None else []) + [np.r_[np.log(rng.uniform(0.1, 0.9)), np.log(10 ** rng.uniform(-6, -3)), rng.normal(0, 1, 3)][:npar] for _ in range(3)]
        for x0 in starts:
            try: r = least_squares(resid, x0, args=(lk, arms), method='lm', max_nfev=300)
            except Exception: continue
            if best is None or r.cost < best.cost: best = r
        warm = best.x; chi2 = 2 * best.cost
        prof.append(dict(log10_ktox=float(lk), chi2=float(chi2), s_inf=float(np.exp(min(best.x[0], 0))), k_dam=float(np.exp(best.x[1])),
                         **(dict(cni=float(sg(best.x[2])), leak=float(sg(best.x[3])), h_drug=float(0.15 * sg(best.x[4]))) if npar == 5 else {})))
        print(tag, prof[-1], flush=True)
    c = np.array([q['chi2'] for q in prof]); dof = len(arms) * len(W) * 2 - npar - 1; red = float(c.min() / dof)
    def interval(scale):
        ok = (c - c.min()) / scale <= 3.84; g = GRID[ok]
        return dict(lo=float(g.min()), hi=float(g.max()), width=float(g.max() - g.min()), bounded_below=bool(not ok[0]), bounded_above=bool(not ok[-1]), contiguous=bool(ok[np.argmax(ok):len(ok) - np.argmax(ok[::-1])].all()))
    out[tag] = dict(profile=prof, chi2_min=float(c.min()), dof=dof, reduced_chi2=red, best_log10_ktox=float(GRID[c.argmin()]), raw=interval(1.0), scaled=interval(max(red, 1.0)))
    print(tag, {k: v for k, v in out[tag].items() if k != 'profile'}, flush=True)
json.dump(out, open('results/r75_real_profile_ktox.json' if len(sys.argv) < 2 else 'results/r75_real_profile_ktox_truncated.json', 'w'), indent=1)
