"""Round 76 -- prereg: project doc claude/prereg-round76.md (written before any fit including AURORA 2).
Profile likelihood of k_tox on AURORA 1 (weeks 2-52) + AURORA 2 (months 12-36, shape only via per-series offsets)."""
import csv, json, sys, numpy as np
from scipy.optimize import least_squares
from model import make_params, simulate, observe
def _load(path):
    with open(path) as f: return list(csv.DictReader([l for l in f if not l.lstrip().startswith('#')]))
A1 = _load('../../data/real/aurora1_per_visit.csv'); A2 = _load('../../data/real/aurora2_per_visit.csv')
BASE = {'treated': (4.14, 78.3, 179), 'placebo': (3.87, 77.4, 178)}
W = [2, 4, 8, 12, 16, 20, 24, 30, 36, 42, 48, 52]; M = [12, 18, 24, 30, 36]
TG = np.array([0.0] + [w * 7.0 for w in W] + [m * 30.44 for m in M])
O1 = {a: np.array([[float(r['upcr_change']), float(r['upcr_sd']) / np.sqrt(BASE[a][2]), float(r['egfr_change']), float(r['egfr_sd']) / np.sqrt(BASE[a][2])] for r in A1 if r['arm'] == a]) for a in BASE}
O2 = {a: np.array([[float(r['upcr_change']), float(r['upcr_se']), float(r['egfr_change']), float(r['egfr_se'])] for r in A2 if r['arm'] == a]) for a in BASE}
sg = lambda z: 1 / (1 + np.exp(-np.clip(z, -50, 50)))
def arm_curve(arm, s_inf, kdam, ktox, kprog, cni=0.0, leak=0.0, h=0.0):
    u0, g0, _ = BASE[arm]; p = make_params(0.0); p['p0'] = u0; p['G0'] = g0 / (1 - p['c_inf']); p['s_inf'] = s_inf * (1 - cni)
    p['k_dam'], p['k_tox'], p['k_prog'] = kdam, ktox, kprog; p['h_drug'] = h; p['leak_red'] = leak; p['drug_on'], p['drug_off'] = 1.0, 1e9
    Y = simulate(p, TG, rtol=1e-8, atol=1e-10); O = np.array([observe(y, p, t) for t, y in zip(TG, Y)]); return O[1:, :2] - O[0, :2], float(Y[-1, 5])
def unpack(x, arms, kp):
    i = 2; kprog = 0.0
    if kp: kprog = np.exp(x[i]); i += 1
    drug = (sg(x[i]), sg(x[i + 1]), 0.15 * sg(x[i + 2])) if 'treated' in arms else (); i += 3 if 'treated' in arms else 0
    return np.exp(min(x[0], 0.0)), np.exp(x[1]), kprog, drug, x[i:]
def resid(x, lk, arms, kp):
    s_inf, kdam, kprog, drug, off = unpack(x, arms, kp); r = []
    try:
        for j, arm in enumerate(arms):
            C, _ = arm_curve(arm, s_inf, kdam, 10.0 ** lk, kprog, *(drug if arm == 'treated' else ()))
            a, b = C[:len(W)], C[len(W):]; o1, o2 = O1[arm], O2[arm]
            r += list((a[:, 0] - o1[:, 0]) / o1[:, 1]) + list((a[:, 1] - o1[:, 2]) / o1[:, 3])
            r += list((b[:, 0] + off[2 * j] - o2[:, 0]) / o2[:, 1]) + list((b[:, 1] + off[2 * j + 1] - o2[:, 2]) / o2[:, 3])
    except Exception:
        return np.full(len(arms) * 34, 1e3)
    return np.array(r)
GRID = np.linspace(-7, -3, 9); out = {}
RUNS = (('PT', ['placebo', 'treated'], False), ('PT_kprog', ['placebo', 'treated'], True), ('P', ['placebo'], False), ('T', ['treated'], False))
for tag, arms, kp in RUNS:
    npar = 2 + int(kp) + (3 if 'treated' in arms else 0) + 2 * len(arms); prof = []; warm = None; rng = np.random.default_rng(0)
    def x0():
        v = [np.log(rng.uniform(0.2, 0.9)), np.log(10 ** rng.uniform(-5, -3))] + ([np.log(10 ** rng.uniform(-7, -4))] if kp else []) + (list(rng.normal(0, 1, 3)) if 'treated' in arms else []) + [0.0] * (2 * len(arms))
        return np.array(v)
    for lk in GRID:
        best = None
        for s in ([warm] if warm is not None else []) + [x0(), x0()]:
            try: r = least_squares(resid, s, args=(lk, arms, kp), method='lm', max_nfev=400)
            except Exception: continue
            if best is None or r.cost < best.cost: best = r
        warm = best.x; s_inf, kdam, kprog, drug, off = unpack(best.x, arms, kp)
        Dfin = {a: arm_curve(a, s_inf, kdam, 10.0 ** lk, kprog, *(drug if a == 'treated' else ()))[1] for a in arms}
        q = dict(log10_ktox=float(lk), chi2=float(2 * best.cost), s_inf=float(s_inf), k_dam=float(kdam), k_prog=float(kprog), drug=[float(z) for z in drug], offsets=[float(z) for z in off], D36=Dfin)
        prof.append(q); print(tag, json.dumps(q), flush=True)
    c = np.array([q['chi2'] for q in prof]); dof = len(arms) * 34 - npar - 1; red = float(c.min() / dof)
    def interval(scale):
        ok = (c - c.min()) / scale <= 3.84; g = GRID[ok]; i0, i1 = int(np.argmax(ok)), len(ok) - int(np.argmax(ok[::-1]))
        return dict(lo=float(g.min()), hi=float(g.max()), bounded_below=bool(not ok[0]), bounded_above=bool(not ok[-1]), contiguous=bool(ok[i0:i1].all()))
    out[tag] = dict(profile=prof, chi2_min=float(c.min()), dof=dof, reduced_chi2=red, best_log10_ktox=float(GRID[c.argmin()]), raw=interval(1.0), scaled=interval(max(red, 1.0)))
    print('SUMMARY', tag, json.dumps({k: v for k, v in out[tag].items() if k != 'profile'}), flush=True)
    json.dump(out, open('results/r76_aurora12_profile_ktox.json', 'w'), indent=1)
