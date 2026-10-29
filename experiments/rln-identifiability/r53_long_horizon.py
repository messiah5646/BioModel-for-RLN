"""
R53: the crossover test at LONG horizon, using data round 52 could not reach.

Round 52 could only see 28 weeks ahead and found the model losing at every horizon.
AURORA 2 extends to month 36. Months 18, 24 and 30 were never used in any round
(rounds 48/49 used only the month 12->36 slope, i.e. the endpoints).

Fit: AURORA 1 weeks 2-24 (reuses round 50's fits -- nothing re-fitted).
Held out: months 18, 24, 30 = 12 points at 15-27 MONTHS ahead.

CONFOUND, declared: AURORA 2 is a continuer subset of AURORA 1, so its level differs
(month-12 UPCR change -3.17 vs AURORA 1 week-52 -2.65). Both raw and offset-corrected
scores are reported. The offset-corrected version tests trajectory SHAPE, which is what
the crossover question is actually about, and it anchors the model at month 12 exactly
as the trivial rules are anchored.
"""
import csv, json
import numpy as np
import r50_trajectory_holdout as R

def load(p):
    with open(p) as f:
        return [r for r in csv.DictReader([l for l in f if not l.lstrip().startswith('#')])]

A2 = {}
for r in load('../../data/real/aurora2_per_visit.csv'):
    A2.setdefault(r['arm'], {})[int(r['month'])] = {
        'du': float(r['upcr_change']), 'dg': float(r['egfr_change']),
        'su': float(r['upcr_se']) * np.sqrt(float(r['upcr_n'])),
        'sg': float(r['egfr_se']) * np.sqrt(float(r['egfr_n']))}

HELD = [18, 24, 30]                 # never used in any round
r50 = json.load(open('results/r50_holdout.json'))
MO_W = {12: 52, 18: 78, 24: 104, 30: 130, 36: 156}   # month -> week


def model_pred(variant, arm, month):
    """Model fitted on AURORA 1 weeks 2-24, evaluated at the AURORA 2 visit."""
    x = np.array(r50[variant]['x'])
    s_inf, scale = np.exp(x[0]), np.exp(x[1])
    cni, leak = 1/(1+np.exp(-x[2])), 1/(1+np.exp(-x[3]))
    h = 0.15/(1+np.exp(-x[4]))
    w = MO_W[month]
    p = (R.arm_params(variant, arm, s_inf, scale) if arm == 'placebo'
         else R.arm_params(variant, arm, s_inf, scale, cni, leak, h))
    import model as M
    tg = np.array([0.0, w*7.0])
    Y = M.simulate(p, tg, rtol=1e-9, atol=1e-11)
    O = np.array([M.observe(y, p, t) for t, y in zip(tg, Y)])
    du = O[1,0]-O[0,0]
    if p.get('_hu',0.0) > 0: du = O[1,0]*(1-p['_hu'])-O[0,0]
    return du, O[1,1]-O[0,1]


def score(fn, offset=False):
    e = []
    for arm in A2:
        ou12, og12 = A2[arm][12]['du'], A2[arm][12]['dg']
        pu12, pg12 = fn(arm, 12)
        cu = (ou12 - pu12) if offset else 0.0
        cg = (og12 - pg12) if offset else 0.0
        for m in HELD:
            o = A2[arm][m]
            pu, pg = fn(arm, m)
            e += [(pu+cu-o['du'])/o['su'], (pg+cg-o['dg'])/o['sg']]
    return float(np.sqrt(np.mean(np.square(e))))


# trivial comparators, anchored at month 12 (last observation before the holdout)
def persistence(arm, m):
    return A2[arm][12]['du'], A2[arm][12]['dg']

def linear(arm, m):
    """month-12 anchor + the most recent AURORA 1 slope available (weeks 24->52)."""
    a, b = R.OBS[arm][24], R.OBS[arm][52]
    du = (b['upcr_change']-a['upcr_change'])/28.0
    dg = (b['egfr_change']-a['egfr_change'])/28.0
    dw = MO_W[m]-52
    return A2[arm][12]['du']+du*dw, A2[arm][12]['dg']+dg*dw

print("R53 -- held out: months 18/24/30 (15-27 months ahead), 12 points, never used")
print(f"{'method':>16} | {'raw RMSE':>10} | {'offset-corrected':>17}")
print("-"*50)
rows = {}
for v in ['V1','V2','V4']:
    fn = (lambda a,m,v=v: model_pred(v,a,m))
    rows[v] = (score(fn), score(fn, True))
    print(f"{v:>16} | {rows[v][0]:10.3f} | {rows[v][1]:17.3f}")
for nm, fn in [('persistence',persistence), ('linear extrap',linear)]:
    rows[nm] = (score(fn), score(fn, True))
    print(f"{nm:>16} | {rows[nm][0]:10.3f} | {rows[nm][1]:17.3f}")

bm = min(rows[v][1] for v in ['V1','V2','V4'])
bt = min(rows['persistence'][1], rows['linear extrap'][1])
print(f"\nOFFSET-CORRECTED (trajectory shape, the fair comparison):")
print(f"  best model {bm:.3f} | best trivial {bt:.3f} | ratio {bm/bt:.3f}")
print(f"  -> {'MODEL WINS at long horizon' if bm < bt else 'trivial still wins'}")
print(f"\nCompare round 52 (<=28 weeks ahead): ratio 1.87-3.00, trivial won at every horizon.")
json.dump({k:{'raw':v[0],'offset':v[1]} for k,v in rows.items()},
          open('results/r53_long_horizon.json','w'), indent=1)
