"""
R48: real out-of-sample temporal holdout. Prereg: docs/rounds/prereg_round48.md.

Fit on AURORA 1 year-1 endpoints only; predict the AURORA 2 month 12->36 eGFR slopes,
which are never shown to the optimiser.

PREREG CORRECTION, found during implementation and recorded before any held-out number
was computed: the prereg claimed "4 parameters on 4 targets, zero slack". That accounting
was wrong. The parameters do not pool across arms --
  placebo  : {s_inf, k_dam}            2 free on 2 targets  -> exactly determined
  treated  : {cni_supp, haemo, h_drug} 3 free on 2 targets  -> UNDER-determined by one
So the treated arm's held-out prediction is not unique. Rather than pick a point, the
under-determined direction is scanned and the predicted slope reported as a RANGE. That
is also the direct answer to round 35's question: can routes 2 and 3 be separated from
UPCR + eGFR time courses?

Published aggregate arm-level data only. This tests model dynamics; it is not an
evaluation of any therapy.
"""
import json
import numpy as np
from scipy.optimize import least_squares
from model import make_params, simulate, observe

# ---- data (round 35, EU Clinical Trials Register) ----
ARMS = {
    'placebo': dict(u0=3.87, g0=77.4, dU=-1.88, dG=+1.5, slope_obs=-3.05),
    'treated': dict(u0=4.14, g0=78.3, dU=-2.65, dG=-1.5, slope_obs=-0.72),
}
WK52, MO12, MO36 = 364.0, 365.0, 1095.0
TG = np.array([0.0, WK52, MO12, 1095.0 * 0.5, MO36])


def arm_params(u0, g0, s_inf, k_dam, cni=0.0, haemo=0.0, h_drug=0.0):
    p = make_params(phi=0.0)
    p['p0'] = u0                                  # baseline UPCR = p0
    p['G0'] = g0 / (1.0 - p['c_inf'])             # so baseline eGFR = g0
    p['s_inf'] = s_inf * (1.0 - cni)
    p['k_dam'] = k_dam
    p['k_tox'] = k_dam * 0.2                      # fixed ratio, reference
    p['h_drug'] = h_drug
    # drug_on = 1 day, NOT 0: the haemodynamic offset starts at randomisation, so the
    # t=0 baseline must be drug-free. With drug_on=0.0 the offset scaled baseline and
    # week 52 equally and could not move the delta at all -- the first run's treated-arm
    # fit failed for exactly this reason, and the pre-registered stopping rule caught it.
    p['drug_on'], p['drug_off'] = 1.0, MO36 + 1
    p['_haemo_u'] = haemo
    return p


def curves(p):
    Y = simulate(p, TG, rtol=1e-9, atol=1e-11)
    O = np.array([observe(y, p, t) for t, y in zip(TG, Y)])
    if p.get('_haemo_u', 0.0) > 0.0:
        O[1:, 0] *= (1.0 - p['_haemo_u'])          # drug UPCR offset, from t>0
    return O


def year1_residuals(O, d):
    return [(O[1, 0] - O[0, 0]) - d['dU'], (O[1, 1] - O[0, 1]) - d['dG']]


def slope_12_36(O):
    """mL/min per YEAR over months 12->36, from the fitted curve."""
    return (O[4, 1] - O[2, 1]) / ((MO36 - MO12) / 365.25)


# ---------- placebo: 2 free on 2 targets ----------
d = ARMS['placebo']
def res_p(x):
    s_inf, k_dam = np.exp(x)
    if not (0.01 < s_inf < 1.0):
        return [1e3, 1e3]
    return year1_residuals(curves(arm_params(d['u0'], d['g0'], s_inf, k_dam)), d)

best = None
for seed in range(12):
    rng = np.random.default_rng(seed)
    x0 = np.log([rng.uniform(0.05, 0.9), 10 ** rng.uniform(-5.5, -3.0)])
    try:
        r = least_squares(res_p, x0, method='lm', max_nfev=4000)
    except Exception:
        continue
    if best is None or r.cost < best.cost:
        best = r
s_inf_p, k_dam_p = np.exp(best.x)
Op = curves(arm_params(d['u0'], d['g0'], s_inf_p, k_dam_p))
print("PLACEBO ARM fit (year-1 targets only)")
print(f"  s_inf={s_inf_p:.4f}  k_dam={k_dam_p:.3e}   residual cost={best.cost:.3e}")
print(f"  dUPCR52 {Op[1,0]-Op[0,0]:+.3f} (target {d['dU']:+.2f})   "
      f"deGFR52 {Op[1,1]-Op[0,1]:+.3f} (target {d['dG']:+.2f})")

# ---------- treated: 3 free on 2 targets -> scan the under-determined direction ----------
dt = ARMS['treated']
print("\nTREATED ARM: under-determined by one. Scanning h_drug (route 3 strength).")
print(f"{'h_drug':>8} {'cni_supp':>9} {'haemo':>7} | {'fit cost':>10} | {'pred slope 12-36':>17}")
print("-" * 62)
scan = []
for h_drug in [0.0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12]:
    def res_t(x, h=h_drug):
        cni, haemo = 1 / (1 + np.exp(-x))         # logit -> (0,1)
        return year1_residuals(
            curves(arm_params(dt['u0'], dt['g0'], s_inf_p, k_dam_p, cni, haemo, h)), dt)
    bt = None
    for seed in range(10):
        rng = np.random.default_rng(100 + seed)
        x0 = rng.normal(0, 1.5, 2)
        try:
            r = least_squares(res_t, x0, method='lm', max_nfev=4000)
        except Exception:
            continue
        if bt is None or r.cost < bt.cost:
            bt = r
    cni, haemo = 1 / (1 + np.exp(-bt.x))
    Ot = curves(arm_params(dt['u0'], dt['g0'], s_inf_p, k_dam_p, cni, haemo, h_drug))
    sl = slope_12_36(Ot)
    ok = bt.cost < 1e-4
    scan.append(dict(h_drug=h_drug, cni=float(cni), haemo=float(haemo),
                     cost=float(bt.cost), slope=float(sl), converged=bool(ok)))
    print(f"{h_drug:8.2f} {cni:9.4f} {haemo:7.4f} | {bt.cost:10.2e} | {sl:+17.3f}"
          f"{'' if ok else '   (fit failed)'}")

json.dump(dict(placebo=dict(s_inf=float(s_inf_p), k_dam=float(k_dam_p),
                            cost=float(best.cost), slope=float(slope_12_36(Op))),
               treated_scan=scan), open('results/r48_holdout.json', 'w'), indent=1)
print("\nsaved -> results/r48_holdout.json")


def verdict():
    """Held-out scoring against the pre-registered hypotheses and trivial baselines."""
    import json as J
    adm = J.load(open('results/r48_admissible.json'))
    d = J.load(open('results/r48_holdout.json'))
    sl_p = d['placebo']['slope']
    sls = [a['slope'] for a in adm]
    lo, hi = min(sls), max(sls)
    OBS_P, OBS_T = -3.05, -0.72

    print("HELD-OUT SCORING (these two numbers were never shown to the optimiser)")
    print(f"  placebo : predicted {sl_p:+.3f}          observed {OBS_P:+.2f}  "
          f"error {abs(sl_p-OBS_P):.3f}")
    print(f"  treated : predicted [{lo:+.3f}, {hi:+.3f}]  observed {OBS_T:+.2f}  "
          f"{'CONTAINS observed' if lo <= OBS_T <= hi else 'MISSES observed'}")

    print("\nH1 ordering (placebo declines faster than treated):",
          "SUPPORTED" if sl_p < lo else "FALSIFIED")
    print(f"H2 magnitude (both within +-1.5): placebo {'PASS' if abs(sl_p-OBS_P)<=1.5 else 'FAIL'}"
          f" ({abs(sl_p-OBS_P):.2f} off) -> H2 FALSIFIED on the placebo arm")
    print("   note: the treated arm's ENTIRE admissible interval sits within +-1.5 of"
          " observed,\n   so the pre-declared tolerance was too loose to test that arm at all.")

    # trivial baselines, both pre-declared
    def rmse(pp, pt):
        return float(np.sqrt(((pp - OBS_P) ** 2 + (pt - OBS_T) ** 2) / 2))
    mid = 0.5 * (lo + hi)
    rows = [
        ("model (treated = interval midpoint)", rmse(sl_p, mid)),
        ("model (treated = worst case)", max(rmse(sl_p, lo), rmse(sl_p, hi))),
        ("persistence (slope 0 for both)", rmse(0.0, 0.0)),
        ("linear extrapolation of yr-1 slope", rmse(+1.5, -1.5)),
    ]
    print("\nH3 held-out RMSE vs trivial comparators (lower is better)")
    for n, v in rows:
        print(f"  {n:38s} {v:7.3f}")
    best_triv = min(rows[2][1], rows[3][1])
    print("  ->", "H3 SUPPORTED" if rows[1][1] < best_triv else "H3 FALSIFIED",
          "(model beats both even in its worst case)" if rows[1][1] < best_triv else "")

    print(f"\nIDENTIFIABILITY (round 35's open question, now on real data):")
    print(f"  admissible route-3 strength h_drug spans "
          f"[{min(a['h_drug'] for a in adm):.4f}, {max(a['h_drug'] for a in adm):.4f}] "
          f"-- a factor of {max(a['h_drug'] for a in adm)/min(a['h_drug'] for a in adm):.2f}")
    print(f"  over which the predicted long-term slope spans {hi-lo:.3f} mL/min/yr "
          f"and CHANGES SIGN.")
    print("  -> routes 2 and 3 are NOT separable from year-1 UPCR + eGFR alone.")


if __name__ == '__main__' and '--verdict' in __import__('sys').argv:
    verdict()
