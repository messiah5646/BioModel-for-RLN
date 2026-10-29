"""
R49: fix round 48's placebo-arm falsification. Prereg: docs/rounds/prereg_round49.md.

Three model variants, EACH exactly determined (2 free placebo parameters on the same 2
year-1 targets), so each gives a point prediction of the held-out month 12->36 slope and
no variant gets extra slack.
"""
import json
import numpy as np
from scipy.optimize import least_squares
from model import make_params, simulate, observe

MO12, MO36 = 365.0, 1095.0
TG = np.array([0.0, 364.0, MO12, 730.0, MO36])
PLAC = dict(u0=3.87, g0=77.4, dU=-1.88, dG=+1.5, slope_obs=-3.05)
TREA = dict(u0=4.14, g0=78.3, dU=-2.65, dG=-1.5, slope_obs=-0.72)


def arm(variant, u0, g0, s_inf, scale, cni=0.0, haemo=0.0, h_drug=0.0):
    p = make_params(phi=0.0)
    p['p0'] = u0
    p['G0'] = g0 / (1.0 - p['c_inf'])
    p['s_inf'] = s_inf * (1.0 - cni)
    if variant == 'V1':       # activity-driven damage only
        p['k_dam'], p['k_prog'] = scale, 0.0
    elif variant == 'V2':     # activity-independent progression only
        p['k_dam'], p['k_prog'] = 0.0, scale
    else:                     # V3: both, fixed ratio
        p['k_dam'], p['k_prog'] = scale, 0.3 * scale
    p['k_tox'] = 0.2 * scale
    p['h_drug'] = h_drug
    p['drug_on'], p['drug_off'] = 1.0, MO36 + 1
    p['_hu'] = haemo
    return p


def curves(p):
    Y = simulate(p, TG, rtol=1e-9, atol=1e-11)
    O = np.array([observe(y, p, t) for t, y in zip(TG, Y)])
    if p.get('_hu', 0.0) > 0.0:
        O[1:, 0] *= (1.0 - p['_hu'])
    return O


def yr1(O, d):
    return [(O[1, 0] - O[0, 0]) - d['dU'], (O[1, 1] - O[0, 1]) - d['dG']]


def slope(O):
    return (O[4, 1] - O[2, 1]) / ((MO36 - MO12) / 365.25)


def fit_placebo(variant):
    def res(x):
        s_inf, sc = np.exp(x)
        if not (0.01 < s_inf < 1.0):
            return [1e3, 1e3]
        return yr1(curves(arm(variant, PLAC['u0'], PLAC['g0'], s_inf, sc)), PLAC)
    best = None
    for seed in range(12):
        rng = np.random.default_rng(seed)
        x0 = np.log([rng.uniform(0.05, 0.95), 10 ** rng.uniform(-6.0, -2.5)])
        try:
            r = least_squares(res, x0, method='lm', max_nfev=5000)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    s_inf, sc = np.exp(best.x)
    return s_inf, sc, best.cost


print("R49 -- held-out model selection on the placebo slope (observed -3.05)")
print(f"{'variant':>8} {'s_inf':>8} {'scale':>11} | {'fit cost':>10} | "
      f"{'pred slope':>11} | {'held-out err':>12}")
print("-" * 72)
out = {}
for v in ['V1', 'V2', 'V3']:
    s_inf, sc, cost = fit_placebo(v)
    if cost > 1e-10:
        print(f"{v:>8} {'':>8} {'':>11} | {cost:10.2e} |  FIT FAILED (stopping rule)")
        out[v] = dict(failed=True, cost=float(cost)); continue
    O = curves(arm(v, PLAC['u0'], PLAC['g0'], s_inf, sc))
    sl = slope(O); err = abs(sl - PLAC['slope_obs'])
    out[v] = dict(s_inf=float(s_inf), scale=float(sc), cost=float(cost),
                  slope=float(sl), err=float(err))
    print(f"{v:>8} {s_inf:8.4f} {sc:11.3e} | {cost:10.2e} | {sl:+11.3f} | {err:12.3f}")

ok = {k: val for k, val in out.items() if not val.get('failed')}
bestv = min(ok, key=lambda k: ok[k]['err'])
print(f"\nH1 (V2 held-out error < 0.75):",
      f"{'SUPPORTED' if ok.get('V2', {}).get('err', 9e9) < 0.75 else 'FALSIFIED'}"
      f"  (V2 err {ok.get('V2',{}).get('err',float('nan')):.3f})")
print(f"H2 (V2 beats V1):",
      "SUPPORTED" if ok['V2']['err'] < ok['V1']['err'] else "FALSIFIED -- R48 diagnosis wrong")
print(f"best variant on held-out error: {bestv} (err {ok[bestv]['err']:.3f})")

# ---- H3: does the better structure buy identifiability on the treated arm? ----
print(f"\nH3 -- treated-arm admissible route-3 interval under {bestv}")
s_inf_b, sc_b = ok[bestv]['s_inf'], ok[bestv]['scale']
adm = []
for h in np.arange(0.0, 0.1201, 0.0025):
    def res_t(x, hh=h):
        cni, haemo = 1 / (1 + np.exp(-x))
        return yr1(curves(arm(bestv, TREA['u0'], TREA['g0'], s_inf_b, sc_b, cni, haemo, hh)), TREA)
    bt = None
    for seed in range(12):
        rng = np.random.default_rng(400 + seed)
        try:
            r = least_squares(res_t, rng.normal(0, 1.5, 2), method='lm', max_nfev=5000)
        except Exception:
            continue
        if bt is None or r.cost < bt.cost:
            bt = r
    if bt.cost < 1e-10:
        cni, haemo = 1 / (1 + np.exp(-bt.x))
        adm.append((float(h), float(slope(curves(arm(bestv, TREA['u0'], TREA['g0'],
                                                     s_inf_b, sc_b, cni, haemo, h))))))
if adm:
    hs = [a[0] for a in adm]; sls = [a[1] for a in adm]
    span = max(sls) - min(sls)
    print(f"  admissible h_drug [{min(hs):.4f}, {max(hs):.4f}]")
    print(f"  predicted slope   [{min(sls):+.3f}, {max(sls):+.3f}]  span {span:.3f} mL/min/yr")
    print(f"  contains observed -0.72: {'YES' if min(sls) <= -0.72 <= max(sls) else 'NO'}")
    print(f"  round 48 span was 1.474 -> H3 {'SUPPORTED' if span < 1.474 else 'FALSIFIED'}")
    out['treated_span'] = span; out['treated_interval'] = [min(sls), max(sls)]
else:
    print("  no admissible values found")
json.dump(out, open('results/r49_variants.json', 'w'), indent=1)
print("\nsaved -> results/r49_variants.json")
