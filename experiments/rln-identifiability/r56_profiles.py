"""Round 56 -- prereg: project doc claude/prereg-round56.md (written 2026-11-02, before any profile).
Run 2026-11-20. Declared deviations (written into the round doc): per-point restarts reduced from 12 to warm-start + 2
random (cost); MAP for the scan = weeks 2-52 fit; profiled parameters are scanned in ln(natural value)."""
import json, os, sys, time, signal, numpy as np
from scipy.optimize import least_squares
src = open('r50_trajectory_holdout.py').read().split('print("R50')[0]; ns = {}; exec(src, ns)
residuals, FIT_W, ALL_W = ns['residuals'], ns['FIT_W'], ns['ALL_W']
W52 = [w for w in ALL_W]                       # weeks 2..52, 12 visits, 48 points
OUTF = 'results/r56_profiles.json'; out = json.load(open(OUTF)) if os.path.exists(OUTF) else {}
class TO(Exception): pass
def _al(*a): raise TO()
signal.signal(signal.SIGALRM, _al)
def R(x, v, weeks):
    signal.setitimer(signal.ITIMER_REAL, 10.0)
    try: r = np.asarray(residuals(x, v, weeks), float)
    except TO: r = np.full(len(weeks) * 4, 1e3)
    signal.setitimer(signal.ITIMER_REAL, 0); return r
def fit_full(v, weeks, seeds=range(12)):
    best = None
    for s in seeds:
        g = np.random.default_rng(s); x0 = np.array([np.log(g.uniform(0.1, 0.9)), np.log(10 ** g.uniform(-6, -3)), g.normal(0, 1), g.normal(0, 1), g.normal(0, 1)])
        r = least_squares(R, x0, args=(v, weeks), method='lm', max_nfev=8000)
        if best is None or r.cost < best.cost: best = r
    return best
# ---- stopping rule: reproduce committed round-50 costs (weeks 2-24)
if 'stop_rule' not in out:
    ref = dict(V1=82.45358848170758, V2=80.88179718190125, V4=82.66574269799922); got = {}
    for v in ref:
        t = time.time(); got[v] = float(fit_full(v, FIT_W).cost); print('stop-rule', v, got[v], 'ref', ref[v], 'diff', got[v] - ref[v], f'{time.time()-t:.0f}s', flush=True)
    out['stop_rule'] = dict(ref=ref, got=got, passed=bool(all(abs(got[v] - ref[v]) <= 1e-6 for v in ref)))
    json.dump(out, open(OUTF, 'w'), indent=1)
if not out['stop_rule']['passed']:
    print('STOP RULE FAILED -- no profiles computed', out['stop_rule']); sys.exit(0)
# ---- profiles. index: 0 log s_inf, 1 log scale, 2 logit cni
def to_nat(i, z): return np.exp(z) if i < 2 else 1 / (1 + np.exp(-z))
def from_nat(i, n): return np.log(n) if i < 2 else np.log(n / (1 - n))
def profile(v, weeks, i, tag):
    if tag in out: return
    b = fit_full(v, weeks); mapn = float(to_nat(i, b.x[i])); grid = np.log(mapn) + np.linspace(-4, 4, 41)
    pts = []; order = np.argsort(np.abs(grid - np.log(mapn)))
    warm = {}
    for k in order:
        lnv = grid[k]; nat = float(np.exp(lnv))
        if i == 2 and nat >= 1: pts.append(dict(ln_value=float(lnv), value=nat, cost=None, note='cni>=1 outside domain')); continue
        fixed = from_nat(i, nat); free = [j for j in range(5) if j != i]
        def rr(y):
            x = np.empty(5); x[free] = y; x[i] = fixed; return R(x, v, weeks)
        nb = min(warm, key=lambda q: abs(q - lnv)) if warm else None
        starts = ([warm[nb]] if nb is not None else [b.x[free]]) + [np.random.default_rng(1000 + k * 7 + s).normal(b.x[free], 1.0) for s in range(2)]
        best = None
        for y0 in starts:
            r = least_squares(rr, y0, method='lm', max_nfev=3000)
            if best is None or r.cost < best.cost: best = r
        warm[lnv] = best.x; pts.append(dict(ln_value=float(lnv), value=nat, cost=float(best.cost)))
    pts.sort(key=lambda d: d['ln_value'])
    c = np.array([p['cost'] if p['cost'] is not None else np.inf for p in pts]); cmin = min(c.min(), b.cost); ok = (c - cmin) <= 3.84 / 2
    vals = np.array([p['value'] for p in pts]); inside = vals[ok]
    res = dict(variant=v, weeks=[int(w) for w in weeks], param=['s_inf', 'scale', 'cni'][i], map=mapn, map_cost=float(b.cost), points=pts,
               lo=float(inside.min()), hi=float(inside.max()), lower_unbounded=bool(ok[0]), upper_unbounded=bool(ok[-1]),
               contiguous=bool(ok[np.argmax(ok):len(ok) - np.argmax(ok[::-1])].all()), frac_width=float((inside.max() - inside.min()) / mapn))
    out[tag] = res; json.dump(out, open(OUTF, 'w'), indent=1)
    print('PROFILE', tag, {k: res[k] for k in ('map', 'map_cost', 'lo', 'hi', 'lower_unbounded', 'upper_unbounded', 'contiguous', 'frac_width')}, flush=True)
for v in ('V1', 'V2'):
    profile(v, W52, 1, f'{v}_scale_w52')
for v in ('V1', 'V2'):
    profile(v, FIT_W, 1, f'{v}_scale_w24')
for v in ('V1', 'V2'):
    for i, nm in ((0, 's_inf'), (2, 'cni')): profile(v, W52, i, f'{v}_{nm}_w52')
print('DONE', flush=True)
