"""Round 78 -- prereg: project doc claude/prereg-round78.md (written before any two-arm fit with steps)."""
import csv, json, signal, numpy as np
from scipy.optimize import least_squares
from model import make_params, simulate, observe
ld = lambda p: list(csv.DictReader([l for l in open(p) if not l.lstrip().startswith('#')]))
A1, A2 = ld('../../data/real/aurora1_per_visit.csv'), ld('../../data/real/aurora2_per_visit.csv')
BASE = {'placebo': (3.87, 77.4, 178), 'treated': (4.14, 78.3, 179)}; ARMS = ['placebo', 'treated']
W = [2, 4, 8, 12, 16, 20, 24, 30, 36, 42, 48, 52]; M = [12, 18, 24, 30, 36]; TG = np.array([0.0] + [w * 7.0 for w in W] + [m * 30.44 for m in M])
O1 = {a: np.array([[float(r['upcr_change']), float(r['upcr_sd']) / np.sqrt(BASE[a][2]), float(r['egfr_change']), float(r['egfr_sd']) / np.sqrt(BASE[a][2])] for r in A1 if r['arm'] == a]) for a in ARMS}
O2 = {a: np.array([[float(r['egfr_change']), float(r['egfr_se'])] for r in A2 if r['arm'] == a]) for a in ARMS}
sg = lambda z: 1 / (1 + np.exp(-np.clip(z, -50, 50)))
class TO(Exception): pass
def _al(*a): raise TO()
signal.signal(signal.SIGALRM, _al)
def curve(arm, s_inf, kdam, ktox, cni=0.0, leak=0.0, h=0.0):
    u0, g0, _ = BASE[arm]; p = make_params(0.0); p['p0'] = u0; p['G0'] = g0 / (1 - p['c_inf']); p['s_inf'] = s_inf * (1 - cni)
    p['k_dam'], p['k_tox'] = kdam, ktox; p['h_drug'] = h; p['leak_red'] = leak; p['drug_on'], p['drug_off'] = 1.0, 1e9
    Y = simulate(p, TG, rtol=1e-8, atol=1e-10); O = np.array([observe(y, p, t) for t, y in zip(TG, Y)]); return O[1:, :2] - O[0, :2]
def unpack(x, steps):
    d = dict(s_inf=np.exp(min(x[0], 0)), kdam=np.exp(x[1]), ktox=np.exp(x[2]), cni=sg(x[3]), leak=sg(x[4]), h=0.15 * sg(x[5])); d['su'], d['sgf'] = (x[6], x[7]) if steps else (0.0, 0.0); return d
def model(x, steps):
    d = unpack(x, steps); out = {}
    for a in ARMS:
        C = curve(a, d['s_inf'], d['kdam'], d['ktox'], *((d['cni'], d['leak'], d['h']) if a == 'treated' else ()))
        C[:, 0] += d['su']; C[:, 1] += d['sgf']; out[a] = C
    return out
def resid(x, steps, idx):
    signal.setitimer(signal.ITIMER_REAL, 8.0)
    try: mo = model(x, steps)
    except Exception: signal.setitimer(signal.ITIMER_REAL, 0); return np.full(4 * len(idx), 1e3)
    signal.setitimer(signal.ITIMER_REAL, 0); r = []
    for a in ARMS:
        C, o = mo[a][idx], O1[a][idx]; r += list((C[:, 0] - o[:, 0]) / o[:, 1]) + list((C[:, 1] - o[:, 2]) / o[:, 3])
    return np.array(r)
def fit(steps, idx, nstart=8):
    best = None
    for s in range(nstart):
        g = np.random.default_rng(s); x0 = [np.log(g.uniform(0.2, 0.95)), np.log(10 ** g.uniform(-5, -3)), np.log(10 ** g.uniform(-8, -5)), g.normal(-1, 1), g.normal(0, 1), g.normal(0, 1)] + ([g.normal(-0.7, 0.4), g.normal(3, 1.5)] if steps else [])
        try: r = least_squares(resid, x0, args=(steps, idx), method='lm', max_nfev=500)
        except Exception: continue
        if best is None or r.cost < best.cost: best = r
    return best
def wls_slope(t, y, se):
    w = 1 / se ** 2; X = np.c_[np.ones_like(t), t]; C = np.linalg.inv(X.T @ (w[:, None] * X)); b = C @ X.T @ (w * y); return float(b[1]), float(np.sqrt(C[1, 1]))
ALL = np.arange(12); FITW = np.arange(7); OUTW = np.array([7, 8, 9, 10]); yrs = np.array(M) / 12.0; res = {}
obs_slope = {a: wls_slope(yrs, O2[a][:, 0], O2[a][:, 1]) for a in ARMS}; print('observed AURORA 2 eGFR slopes (mL/min/yr, se):', obs_slope)
for steps in (False, True):
    b = fit(steps, ALL); d = unpack(b.x, steps); mo = model(b.x, steps); rr = resid(b.x, steps, ALL); per = {a: float((rr[i * 24:(i + 1) * 24] ** 2).sum()) for i, a in enumerate(ARMS)}
    sl = {a: float(np.polyfit(yrs, mo[a][12:, 1], 1)[0]) for a in ARMS}
    z = {a: (sl[a] - obs_slope[a][0]) / obs_slope[a][1] for a in ARMS}
    yr1 = {a: float(mo[a][11, 1] - mo[a][0, 1]) for a in ARMS}
    # R50 split
    bh = fit(steps, FITW, nstart=6); mh = model(bh.x, steps); e = []
    for a in ARMS: e += list(mh[a][OUTW, 0] - O1[a][OUTW, 0]) + list(mh[a][OUTW, 1] - O1[a][OUTW, 2])
    eu = [mh[a][OUTW, 0] - O1[a][OUTW, 0] for a in ARMS]; eg = [mh[a][OUTW, 1] - O1[a][OUTW, 2] for a in ARMS]
    tag = 'steps' if steps else 'none'
    res[tag] = dict(chi2=float(2 * b.cost), chi2_by_arm=per, params={k: float(v) for k, v in d.items()}, a2_slope_model=sl, a2_slope_z=z, egfr_wk2_to_wk52_model=yr1,
                    holdout_rmse=float(np.sqrt(np.mean(np.square(e)))), holdout_rmse_upcr=float(np.sqrt(np.mean(np.square(eu)))), holdout_rmse_egfr=float(np.sqrt(np.mean(np.square(eg)))), holdout_fit_params={k: float(v) for k, v in unpack(bh.x, steps).items()})
    print(tag, json.dumps(res[tag]), flush=True)
n, s = res['none'], res['steps']
v = dict(H1=s['chi2'] < 0.5 * n['chi2'], H2=s['chi2_by_arm']['treated'] <= 1.1 * n['chi2_by_arm']['treated'], H3a=abs(s['a2_slope_z']['placebo']) <= 2, H3b=s['a2_slope_z']['treated'] < -2, H4=s['holdout_rmse'] > 1.137)
v = {k: bool(x) for k, x in v.items()}; print(v); json.dump(dict(observed_a2_slope=obs_slope, fits=res, verdict=v, trivial=dict(persistence=1.307, linear=1.137)), open('results/r78_shared_step.json', 'w'), indent=1)
