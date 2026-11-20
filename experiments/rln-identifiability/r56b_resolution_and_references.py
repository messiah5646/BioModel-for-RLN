"""Round 56 supplement (same day). (a) Interpolated interval ends: the pre-registered 41-point grid (0.2 ln-unit step)
is coarser than the intervals, so ends are found by linear interpolation of delta-cost between grid points.
(b) Pre-registered section-5 reference points not computed in the first script: s_inf from the placebo arm alone,
and the leak parameter (route-3 reference) under V2. (c) Same intervals with delta-cost scaled by reduced chi2 (reported, not used for verdicts)."""
import json, signal, numpy as np
from scipy.optimize import least_squares
src = open('r50_trajectory_holdout.py').read().split('print("R50')[0]; ns = {}; exec(src, ns)
residuals, ALL_W, BASE, OBS, arm_params, curves = ns['residuals'], ns['ALL_W'], ns['BASE'], ns['OBS'], ns['arm_params'], ns['curves']
d = json.load(open('results/r56_profiles.json')); TH = 3.84 / 2
def interval(pts, mapv, mapc, scale=1.0):
    x = np.array([p['ln_value'] for p in pts if p['cost'] is not None]); c = (np.array([p['cost'] for p in pts if p['cost'] is not None]) - mapc) / scale
    m = np.log(mapv); x = np.r_[x, m]; c = np.r_[c, 0.0]; o = np.argsort(x); x, c = x[o], c[o]; k = int(np.argmin(np.abs(x - m)))
    def cross(direction):
        j = k
        while 0 <= j + direction < len(x):
            if c[j + direction] > TH: f = (TH - c[j]) / (c[j + direction] - c[j]); return x[j] + f * (x[j + direction] - x[j]), False
            j += direction
        return x[j], True
    lo, ul = cross(-1); hi, uh = cross(+1)
    return dict(lo=float(np.exp(lo)), hi=float(np.exp(hi)), frac_width=float((np.exp(hi) - np.exp(lo)) / mapv), ln_width=float(hi - lo), lower_unbounded=ul, upper_unbounded=uh)
res = {}
for t, o in d.items():
    if t == 'stop_rule': continue
    red = 2 * o['map_cost'] / (4 * len(o['weeks']) - 5)
    res[t] = dict(raw=interval(o['points'], o['map'], o['map_cost']), scaled=interval(o['points'], o['map'], o['map_cost'], max(red, 1.0)), reduced_chi2=red)
    print(t, 'red chi2 %.2f' % red, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in res[t]['raw'].items()})
# ---- references
class TO(Exception): pass
signal.signal(signal.SIGALRM, lambda *a: (_ for _ in ()).throw(TO()))
def r_plac(x, v, weeks):
    signal.setitimer(signal.ITIMER_REAL, 10.0)
    try:
        p = arm_params(v, 'placebo', np.exp(x[0]), np.exp(x[1])); O = curves(p); out = []
        for i, w in enumerate(ALL_W):
            if w in weeks: o = OBS['placebo'][w]; n = BASE['placebo'][2]; out += [(O[i + 1, 0] - O[0, 0] - o['upcr_change']) / (o['upcr_sd'] / np.sqrt(n)), (O[i + 1, 1] - O[0, 1] - o['egfr_change']) / (o['egfr_sd'] / np.sqrt(n))]
    except TO: out = [1e3] * (2 * len(weeks))
    signal.setitimer(signal.ITIMER_REAL, 0); return np.array(out)
def prof(fun, nfree, idx, x_map, nat_map, logit=False):
    pts = []
    for dl in np.linspace(-4, 4, 41):
        nat = nat_map * np.exp(dl)
        if logit and nat >= 1: continue
        fx = np.log(nat / (1 - nat)) if logit else np.log(nat); free = [j for j in range(nfree) if j != idx]
        def rr(y):
            x = np.empty(nfree); x[free] = y; x[idx] = fx; return fun(x)
        best = min((least_squares(rr, y0, method='lm', max_nfev=3000) for y0 in [x_map[free]] + [np.random.default_rng(s).normal(x_map[free], 1) for s in range(2)]), key=lambda r: r.cost)
        pts.append(dict(ln_value=float(np.log(nat)), value=float(nat), cost=float(best.cost)))
    return pts
refs = {}
for v in ('V1', 'V2'):
    b = min((least_squares(r_plac, np.array([np.log(g.uniform(0.1, 0.9)), np.log(10 ** g.uniform(-6, -3))]), args=(v, ALL_W), method='lm', max_nfev=8000) for g in [np.random.default_rng(s) for s in range(12)]), key=lambda r: r.cost)
    pts = prof(lambda x: r_plac(x, v, ALL_W), 2, 0, b.x, float(np.exp(b.x[0])))
    refs[f'{v}_s_inf_placebo_only'] = dict(map=float(np.exp(b.x[0])), map_cost=float(b.cost), raw=interval(pts, float(np.exp(b.x[0])), float(b.cost)))
    print('REF', v, 's_inf placebo only', refs[f'{v}_s_inf_placebo_only'], flush=True)
bfull = min((least_squares(lambda x: np.asarray(residuals(x, 'V2', ALL_W)), np.array([np.log(g.uniform(0.1, 0.9)), np.log(10 ** g.uniform(-6, -3)), g.normal(0, 1), g.normal(0, 1), g.normal(0, 1)]), method='lm', max_nfev=8000) for g in [np.random.default_rng(s) for s in range(12)]), key=lambda r: r.cost)
lm = float(1 / (1 + np.exp(-bfull.x[3]))); pts = prof(lambda x: np.asarray(residuals(x, 'V2', ALL_W)), 5, 3, bfull.x, lm, logit=True)
refs['V2_leak_w52'] = dict(map=lm, map_cost=float(bfull.cost), raw=interval(pts, lm, float(bfull.cost))); print('REF V2 leak', refs['V2_leak_w52'], flush=True)
json.dump(dict(intervals=res, references=refs), open('results/r56b_resolution_and_references.json', 'w'), indent=1)
