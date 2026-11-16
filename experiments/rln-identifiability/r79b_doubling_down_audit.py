"""Round 79b -- DOUBLING-DOWN AUDIT requested by the user. Exploratory, not pre-registered.
Question: once the day-0 step is in, does the six-state MECHANISM add anything, or would a purely descriptive curve with
the same number of parameters do as well? If descriptive models match it, the 'progress' of rounds 77-79 belongs to the
step alone and continuing to elaborate the mechanism would be doubling down.
Descriptive rivals (no biology, no damage state), fitted to the same AURORA 1 arm means with the same weights:
  D8  shared step + per-arm treated offset + per-arm linear trend, per observable           (8 parameters, = mechanistic)
  D10 round 38's comparator: per arm, UPCR exponential-to-asymptote (3) + eGFR straight line (2)   (10 parameters)
  D6  shared step only + per-arm treated offset + ONE shared linear trend per observable     (6 parameters)"""
import csv, json, numpy as np
from scipy.optimize import least_squares
A1 = list(csv.DictReader([l for l in open('../../data/real/aurora1_per_visit.csv') if not l.lstrip().startswith('#')]))
N = {'placebo': 178, 'treated': 179}; ARMS = ['placebo', 'treated']; W = np.array([2, 4, 8, 12, 16, 20, 24, 30, 36, 42, 48, 52], float)
O = {a: np.array([[float(r['upcr_change']), float(r['upcr_sd']) / np.sqrt(N[a]), float(r['egfr_change']), float(r['egfr_sd']) / np.sqrt(N[a])] for r in A1 if r['arm'] == a]) for a in ARMS}
def D8(x, w):  return {'placebo': (x[0] + x[4] * w, x[1] + x[5] * w), 'treated': (x[0] + x[2] + x[6] * w, x[1] + x[3] + x[7] * w)}
def D6(x, w):  return {'placebo': (x[0] + x[4] * w, x[1] + x[5] * w), 'treated': (x[0] + x[2] + x[4] * w, x[1] + x[3] + x[5] * w)}
def D10(x, w): return {a: (x[5 * i] + x[5 * i + 1] * (1 - np.exp(-w / np.exp(x[5 * i + 2]))), x[5 * i + 3] + x[5 * i + 4] * w) for i, a in enumerate(ARMS)}
MODELS = {'D6': (D6, 6, lambda g: g.normal(0, 1, 6)), 'D8': (D8, 8, lambda g: g.normal(0, 1, 8)), 'D10': (D10, 10, lambda g: np.r_[g.normal(0, 1, 2), g.normal(2, 1), g.normal(0, 1, 2), g.normal(0, 1, 2), g.normal(2, 1), g.normal(0, 1, 2)])}
def res(x, f, idx):
    p = f(x, W[idx]); r = []
    for a in ARMS: r += list((p[a][0] - O[a][idx, 0]) / O[a][idx, 1]) + list((p[a][1] - O[a][idx, 2]) / O[a][idx, 3])
    return np.array(r)
def fit(f, mk, idx):
    best = None
    for s in range(30):
        r = least_squares(res, mk(np.random.default_rng(s)), args=(f, idx))
        if best is None or r.cost < best.cost: best = r
    return best
def rmse(f, x, out):
    p = f(x, W[out]); e = []
    for a in ARMS: e += list(p[a][0] - O[a][out, 0]) + list(p[a][1] - O[a][out, 2])
    return float(np.sqrt(np.mean(np.square(e))))
ALL = np.arange(12); S1 = (np.arange(7), np.array([7, 8, 9, 10])); S2 = (np.arange(9), np.array([9, 10, 11])); out = {}
for nm, (f, k, mk) in MODELS.items():
    b = fit(f, mk, ALL); chi = float(2 * b.cost)
    out[nm] = dict(k=k, chi2=chi, aic=chi + 2 * k, split1=rmse(f, fit(f, mk, S1[0]).x, S1[1]), split2=rmse(f, fit(f, mk, S2[0]).x, S2[1]))
m = json.load(open('results/r78_shared_step.json'))['fits']; rb = json.load(open('results/r78b_recheck.json'))
out['mech_step'] = dict(k=8, chi2=m['steps']['chi2'], aic=m['steps']['chi2'] + 16, split1=m['steps']['holdout_rmse'], split2=rb['5_second_split_fit_2to36_predict_42to52']['steps'])
out['mech_nostep'] = dict(k=6, chi2=m['none']['chi2'], aic=m['none']['chi2'] + 12, split1=m['none']['holdout_rmse'], split2=rb['5_second_split_fit_2to36_predict_42to52']['none'])
out['trivial'] = dict(split1=dict(persistence=1.3074, linear=1.1371), split2=dict(persistence=rb['5_second_split_fit_2to36_predict_42to52']['persistence'], linear=rb['5_second_split_fit_2to36_predict_42to52']['linear']))
print(f"{'model':12s} {'k':>3s} {'chi2':>8s} {'AIC':>8s} {'split1':>7s} {'split2':>7s}")
for nm in ('mech_nostep', 'mech_step', 'D6', 'D8', 'D10'): o = out[nm]; print(f"{nm:12s} {o['k']:3d} {o['chi2']:8.1f} {o['aic']:8.1f} {o['split1']:7.3f} {o['split2']:7.3f}")
print('trivial', out['trivial']); json.dump(out, open('results/r79b_doubling_down_audit.json', 'w'), indent=1)
