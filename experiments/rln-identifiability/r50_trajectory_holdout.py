"""
R50: trajectory holdout. Prereg: docs/rounds/prereg_round50.md.

Fit weeks 2-24 (28 points), predict weeks 30-48 (16 points never used in any round).
Week 52 excluded from scoring because rounds 48-49 fitted to it.

V1 activity-driven damage | V2 activity-independent progression |
V4 = V2 + the route-3 repair (drug protein-leak reduction enters the DAMAGE DRIVER).
All three carry 5 free parameters, so the comparison is parameter-count fair.
"""
import json
import numpy as np
from scipy.optimize import least_squares
from model import make_params, simulate, observe

import csv
def _load(path):
    with open(path) as f:
        lines = [l for l in f if not l.lstrip().startswith('#')]
    return list(csv.DictReader(lines))
D = _load('../../data/real/aurora1_per_visit.csv')
BASE = {'treated': (4.14, 78.3, 179), 'placebo': (3.87, 77.4, 178)}
FIT_W = [2, 4, 8, 12, 16, 20, 24]
OUT_W = [30, 36, 42, 48]
ALL_W = sorted(set(FIT_W + OUT_W + [52]))
TG = np.array([0.0] + [w * 7.0 for w in ALL_W])


def rows(arm):
    out = {}
    for r in D:
        if r['arm'] == arm:
            out[int(r['week'])] = {k: (float(v) if k != 'arm' else v) for k, v in r.items()}
    return out


OBS = {a: rows(a) for a in BASE}


def arm_params(variant, arm, s_inf, scale, cni=0.0, leak=0.0, h_drug=0.0):
    u0, g0, _ = BASE[arm]
    p = make_params(phi=0.0)
    p['p0'] = u0
    p['G0'] = g0 / (1.0 - p['c_inf'])
    p['s_inf'] = s_inf * (1.0 - cni)
    if variant == 'V1':
        p['k_dam'], p['k_prog'] = scale, 0.0
    else:
        p['k_dam'], p['k_prog'] = 0.0, scale
    p['k_tox'] = 0.2 * scale
    p['h_drug'] = h_drug
    p['drug_on'], p['drug_off'] = 1.0, 1e9
    if variant == 'V4':
        p['leak_red'] = leak           # enters the ODE damage driver AND the observable
        p['_hu'] = 0.0
    else:
        p['leak_red'] = 0.0
        p['_hu'] = leak                # observation-only, as in rounds 48-49
    return p


def curves(p):
    Y = simulate(p, TG, rtol=1e-9, atol=1e-11)
    O = np.array([observe(y, p, t) for t, y in zip(TG, Y)])
    if p.get('_hu', 0.0) > 0.0:
        O[1:, 0] *= (1.0 - p['_hu'])
    return O


def predict(variant, x):
    """x = [log s_inf, log scale, logit cni, logit leak, logit h_drug*]"""
    s_inf, scale = np.exp(x[0]), np.exp(x[1])
    cni, leak = 1 / (1 + np.exp(-x[2])), 1 / (1 + np.exp(-x[3]))
    h_drug = 0.15 / (1 + np.exp(-x[4]))        # bounded to a physiological range
    out = {}
    for arm in BASE:
        if arm == 'placebo':
            p = arm_params(variant, arm, s_inf, scale)
        else:
            p = arm_params(variant, arm, s_inf, scale, cni, leak, h_drug)
        O = curves(p)
        out[arm] = {w: (O[i + 1, 0] - O[0, 0], O[i + 1, 1] - O[0, 1])
                    for i, w in enumerate(ALL_W)}
    return out


def residuals(x, variant, weeks):
    r = []
    try:
        pr = predict(variant, x)
    except Exception:
        return [1e3] * (len(weeks) * 4)
    for arm in BASE:
        n = BASE[arm][2]
        for w in weeks:
            o = OBS[arm][w]
            se_u = o['upcr_sd'] / np.sqrt(n)
            se_g = o['egfr_sd'] / np.sqrt(n)
            r.append((pr[arm][w][0] - o['upcr_change']) / se_u)
            r.append((pr[arm][w][1] - o['egfr_change']) / se_g)
    return r


def fit(variant):
    best = None
    for seed in range(12):
        rng = np.random.default_rng(seed)
        x0 = np.array([np.log(rng.uniform(0.1, 0.9)), np.log(10 ** rng.uniform(-6, -3)),
                       rng.normal(0, 1), rng.normal(0, 1), rng.normal(0, 1)])
        try:
            r = least_squares(residuals, x0, args=(variant, FIT_W), method='lm', max_nfev=8000)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    return best


def heldout_rmse(pred_fn):
    e = []
    for arm in BASE:
        for w in OUT_W:
            o = OBS[arm][w]
            pu, pg = pred_fn(arm, w)
            e += [pu - o['upcr_change'], pg - o['egfr_change']]
    return float(np.sqrt(np.mean(np.square(e)))), e


print("R50 -- fit weeks 2-24 (28 pts), predict weeks 30-48 (16 pts, never used)")
print(f"{'model':>6} | {'fit cost':>10} | {'held-out RMSE':>14} | {'UPCR rmse':>10} | {'eGFR rmse':>10}")
print("-" * 64)
res = {}
for v in ['V1', 'V2', 'V4']:
    b = fit(v)
    pr = predict(v, b.x)
    rmse, errs = heldout_rmse(lambda a, w, pr=pr: pr[a][w])
    eu = np.sqrt(np.mean([errs[i] ** 2 for i in range(0, len(errs), 2)]))
    eg = np.sqrt(np.mean([errs[i] ** 2 for i in range(1, len(errs), 2)]))
    res[v] = dict(cost=float(b.cost), rmse=rmse, rmse_upcr=float(eu), rmse_egfr=float(eg),
                  x=[float(z) for z in b.x])
    print(f"{v:>6} | {b.cost:10.2f} | {rmse:14.4f} | {eu:10.4f} | {eg:10.4f}")

# ---- pre-declared trivial comparators ----
def persistence(arm, w):
    o24 = OBS[arm][24]
    return o24['upcr_change'], o24['egfr_change']

def linear(arm, w):
    a, b = OBS[arm][16], OBS[arm][24]
    du = (b['upcr_change'] - a['upcr_change']) / 8.0
    dg = (b['egfr_change'] - a['egfr_change']) / 8.0
    return b['upcr_change'] + du * (w - 24), b['egfr_change'] + dg * (w - 24)

for nm, fn in [('persistence', persistence), ('linear extrap', linear)]:
    rmse, errs = heldout_rmse(fn)
    eu = np.sqrt(np.mean([errs[i] ** 2 for i in range(0, len(errs), 2)]))
    eg = np.sqrt(np.mean([errs[i] ** 2 for i in range(1, len(errs), 2)]))
    res[nm] = dict(rmse=rmse, rmse_upcr=float(eu), rmse_egfr=float(eg))
    print(f"{nm:>13} | {'':>10} | {rmse:14.4f} | {eu:10.4f} | {eg:10.4f}")

json.dump(res, open('results/r50_holdout.json', 'w'), indent=1)
best_model = min(['V1', 'V2', 'V4'], key=lambda k: res[k]['rmse'])
best_triv = min(res['persistence']['rmse'], res['linear extrap']['rmse'])
print(f"\nH1 (V4 beats V1 and V2): "
      f"{'SUPPORTED' if res['V4']['rmse'] < min(res['V1']['rmse'], res['V2']['rmse']) else 'FALSIFIED'}")
print(f"H2 (V4 beats both trivial): "
      f"{'SUPPORTED' if res['V4']['rmse'] < best_triv else 'FALSIFIED'}"
      f"   (best trivial {best_triv:.4f})")
print(f"best overall: {best_model if res[best_model]['rmse'] < best_triv else 'A TRIVIAL RULE'}")


def h3_and_horizon():
    """H3 signs, plus the reconciliation with round 48's opposite verdict."""
    r = json.load(open('results/r50_holdout.json'))
    print("\nH3 -- sign of the week 24->48 eGFR trend")
    obs = {a: OBS[a][48]['egfr_change'] - OBS[a][24]['egfr_change'] for a in BASE}
    print(f"  observed: treated {obs['treated']:+.2f}, placebo {obs['placebo']:+.2f} "
          f"(both falling, treated faster)")
    allok = True
    for v in ['V1', 'V2', 'V4']:
        pr = predict(v, np.array(r[v]['x']))
        d = {a: pr[a][48][1] - pr[a][24][1] for a in BASE}
        ok = (d['treated'] < 0) and (d['placebo'] < 0) and (d['treated'] < d['placebo'])
        allok &= ok
        print(f"  {v}: treated {d['treated']:+.2f}, placebo {d['placebo']:+.2f}  "
              f"-> {'signs+ordering OK' if ok else 'WRONG'}")
    print(f"  H3: {'SUPPORTED' if allok else 'FALSIFIED for at least one variant'}")

    print("\nRECONCILIATION with round 48 -- the verdict depends on the HORIZON")
    print("  (RMSEs are on different scales between the two tests; only the ORDERING")
    print("   within each test is comparable)")
    print(f"  {'test':>34} | {'horizon':>10} | winner")
    print("  " + "-" * 66)
    print(f"  {'R50 weeks 30-48 (16 pts)':>34} | {'<=6 months':>10} | "
          f"TRIVIAL (linear {r['linear extrap']['rmse']:.2f} < best model {min(r[v]['rmse'] for v in ['V1','V2','V4']):.2f})")
    print(f"  {'R48 month 12-36 slope (2 pts)':>34} | {'2 years':>10} | "
          f"MODEL (1.18-1.32 < persistence 2.22, linear 3.26)")
    print("  -> the ordering REVERSES. Trivial rules win at short horizons, where the")
    print("     trajectory is near-flat; the model wins at long horizons, where linear")
    print("     extrapolation diverges. Round 48's 'beats the baselines' was horizon-")
    print("     specific and should not have been stated unconditionally.")


if __name__ == '__main__' and '--extra' in __import__('sys').argv:
    h3_and_horizon()
