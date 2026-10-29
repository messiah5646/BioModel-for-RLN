"""
R52: locate the crossover horizon. Round 50 section 6's top item.

Round 50 found trivial rules beat the model over weeks 30-48 pooled, while round 48 found
the model beat them on a 2-year slope. That predicts a CROSSOVER. This scores both
per-horizon on the same fit, so the claim becomes a curve rather than two verdicts.

Reuses round 50's fits on weeks 2-24 -- nothing is re-fitted, so no holdout is re-burned.
"""
import json
import numpy as np
import r50_trajectory_holdout as R

HOR = [30, 36, 42, 48, 52]          # weeks; 52 was used by rounds 48-49, flagged in output
r50 = json.load(open('results/r50_holdout.json'))


def err_at(pred_fn, w):
    """RMSE across both arms and both endpoints at a single horizon, standardised by
    each endpoint's own scale so UPCR and eGFR are commensurable."""
    e = []
    for arm in R.BASE:
        o = R.OBS[arm][w]
        pu, pg = pred_fn(arm, w)
        # standardise by the endpoint's observed SD at that visit -> dimensionless
        e.append((pu - o['upcr_change']) / o['upcr_sd'])
        e.append((pg - o['egfr_change']) / o['egfr_sd'])
    return float(np.sqrt(np.mean(np.square(e))))


preds = {}
for v in ['V1', 'V2', 'V4']:
    pr = R.predict(v, np.array(r50[v]['x']))
    preds[v] = (lambda a, w, pr=pr: pr[a][w])
preds['persistence'] = R.persistence
preds['linear extrap'] = R.linear

print("R52 -- standardised error vs prediction horizon (fit = weeks 2-24)")
print("(error in units of the endpoint's own SD, so UPCR and eGFR are comparable)\n")
names = ['V2', 'V1', 'V4', 'persistence', 'linear extrap']
print(f"{'horizon':>9} {'ahead':>7} | " + " ".join(f"{n:>13}" for n in names) + " | winner")
print("-" * 98)
curve = {}
for w in HOR:
    row = {n: err_at(preds[n], w) for n in names}
    curve[w] = row
    best = min(row, key=row.get)
    tag = "  (wk52 used by R48/49)" if w == 52 else ""
    print(f"{'week '+str(w):>9} {str(w-24)+'wk':>7} | "
          + " ".join(f"{row[n]:13.4f}" for n in names)
          + f" | {best}{tag}")

print("\nMODEL vs BEST TRIVIAL, ratio (>1 means trivial wins)")
print(f"{'horizon':>9} | {'best model':>11} | {'best trivial':>12} | {'ratio':>7} | trend")
print("-" * 60)
prev = None
for w in HOR:
    bm = min(curve[w][n] for n in ['V1', 'V2', 'V4'])
    bt = min(curve[w]['persistence'], curve[w]['linear extrap'])
    ratio = bm / bt
    arrow = '' if prev is None else ('closing' if ratio < prev else 'widening')
    print(f"{'week '+str(w):>9} | {bm:11.4f} | {bt:12.4f} | {ratio:7.3f} | {arrow}")
    prev = ratio

r0 = min(curve[30][n] for n in ['V1','V2','V4']) / min(curve[30]['persistence'], curve[30]['linear extrap'])
r1 = min(curve[48][n] for n in ['V1','V2','V4']) / min(curve[48]['persistence'], curve[48]['linear extrap'])
print(f"\nratio at 6wk ahead {r0:.3f} -> at 24wk ahead {r1:.3f}")
if r1 < r0:
    slope = (r1 - r0) / 18.0
    need = (1.0 - r0) / slope if slope < 0 else float('inf')
    print(f"  The gap IS closing with horizon ({(r0-r1)/r0*100:.1f}% over 18 weeks).")
    print(f"  Linear extrapolation of the ratio reaches parity at ~{24 - 24 + need:.0f} weeks"
          f" ahead (~{need/4.3:.1f} months).")
    print("  CAUTION: extrapolating a ratio from 5 points is exactly the move round 50")
    print("  showed fails. This locates a DIRECTION, not a crossover point.")
else:
    print("  The gap is NOT closing within this window -> round 50's horizon account")
    print("  is not supported by the short-horizon data, and round 48's long-horizon")
    print("  win needs a different explanation.")
json.dump(curve, open('results/r52_crossover.json', 'w'), indent=1)
