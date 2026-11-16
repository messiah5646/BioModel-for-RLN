"""Round 79 -- prereg: project doc claude/prereg-round79.md (written before any joint fit with the step)."""
import csv, json, signal, os, numpy as np
from scipy.optimize import least_squares
from model import make_params, simulate, observe
ld = lambda p: list(csv.DictReader([l for l in open(p) if not l.lstrip().startswith('#')]))
A1, A2 = ld('../../data/real/aurora1_per_visit.csv'), ld('../../data/real/aurora2_per_visit.csv')
BASE = {'placebo': (3.87, 77.4, 178), 'treated': (4.14, 78.3, 179)}; ARMS = ['placebo', 'treated']
W = [2, 4, 8, 12, 16, 20, 24, 30, 36, 42, 48, 52]; M = [12, 18, 24, 30, 36]; TG = np.array([0.0] + [w * 7.0 for w in W] + [m * 30.44 for m in M])
O1 = {a: np.array([[float(r['upcr_change']), float(r['upcr_sd']) / np.sqrt(BASE[a][2]), float(r['egfr_change']), float(r['egfr_sd']) / np.sqrt(BASE[a][2])] for r in A1 if r['arm'] == a]) for a in ARMS}
O2 = {a: np.array([[float(r['upcr_change']), float(r['upcr_se']), float(r['egfr_change']), float(r['egfr_se'])] for r in A2 if r['arm'] == a]) for a in ARMS}
sg = lambda z: 1 / (1 + np.exp(-np.clip(z, -50, 50)))
class TO(Exception): pass
def _al(*a): raise TO()
signal.signal(signal.SIGALRM, _al)
def curve(arm, s_inf, kdam, ktox, cni=0.0, leak=0.0, h=0.0):
    u0, g0, _ = BASE[arm]; p = make_params(0.0); p['p0'] = u0; p['G0'] = g0 / (1 - p['c_inf']); p['s_inf'] = s_inf * (1 - cni)
    p['k_dam'], p['k_tox'] = kdam, ktox; p['h_drug'] = h; p['leak_red'] = leak; p['drug_on'], p['drug_off'] = 1.0, 1e9
    Y = simulate(p, TG, rtol=1e-8, atol=1e-10); O = np.array([observe(y, p, t) for t, y in zip(TG, Y)]); return O[1:, :2] - O[0, :2], float(Y[-1, 5])
# x layout: [log s_inf, log kdam, leak, h, su, sg, off*4] + [cni] if free + [log ktox] if free
def unpack(x, cni_free, kt):   # kt: 'free' | float (fixed value)
    d = dict(s_inf=np.exp(min(x[0], 0)), kdam=np.exp(x[1]), leak=sg(x[2]), h=0.15 * sg(x[3]), su=x[4], sgf=x[5], off=list(x[6:10])); i = 10
    d['cni'] = sg(x[i]) if cni_free else 0.0; i += int(cni_free); d['ktox'] = np.exp(x[i]) if kt == 'free' else kt; return d
def predict(d):
    out = {}
    for a in ARMS:
        C, D = curve(a, d['s_inf'], d['kdam'], d['ktox'], *((d['cni'], d['leak'], d['h']) if a == 'treated' else ())); C[:, 0] += d['su']; C[:, 1] += d['sgf']; out[a] = (C, D)
    return out
def resid(x, cni_free, kt):
    d = unpack(x, cni_free, kt); signal.setitimer(signal.ITIMER_REAL, 8.0)
    try: mo = predict(d)
    except Exception: signal.setitimer(signal.ITIMER_REAL, 0); return np.full(68, 1e3)
    signal.setitimer(signal.ITIMER_REAL, 0); r = []
    for j, a in enumerate(ARMS):
        C = mo[a][0]; a1, a2 = C[:12], C[12:]; o1, o2 = O1[a], O2[a]
        r += list((a1[:, 0] - o1[:, 0]) / o1[:, 1]) + list((a1[:, 1] - o1[:, 2]) / o1[:, 3]) + list((a2[:, 0] + d['off'][2 * j] - o2[:, 0]) / o2[:, 1]) + list((a2[:, 1] + d['off'][2 * j + 1] - o2[:, 2]) / o2[:, 3])
    return np.array(r)
def start(g, cni_free, kt):
    return np.array([np.log(g.uniform(0.3, 0.95)), np.log(10 ** g.uniform(-4.5, -3)), g.normal(-1, 1), g.normal(0, 1), g.normal(-1, 0.3), g.normal(4, 1), 0, 0, 0, 0] + ([g.normal(-1, 1.5)] if cni_free else []) + ([np.log(10 ** g.uniform(-7, -4))] if kt == 'free' else []))
def fit(cni_free, kt, seeds, warm=None):
    best = None; costs = []
    for x0 in ([warm] if warm is not None else []) + [start(np.random.default_rng(s), cni_free, kt) for s in seeds]:
        try: r = least_squares(resid, x0, args=(cni_free, kt), method='lm', max_nfev=500)
        except Exception: continue
        costs.append(round(float(2 * r.cost), 3))
        if best is None or r.cost < best.cost: best = r
    return best, sorted(costs)
def summarise(b, cni_free, kt, costs):
    d = unpack(b.x, cni_free, kt); mo = predict(d); rr = resid(b.x, cni_free, kt); yrs = np.array(M) / 12.0; parts = {}
    for j, a in enumerate(ARMS):
        s = rr[j * 34:(j + 1) * 34]; parts[a] = dict(a1_upcr=float((s[:12] ** 2).sum()), a1_egfr=float((s[12:24] ** 2).sum()), a2_upcr=float((s[24:29] ** 2).sum()), a2_egfr=float((s[29:] ** 2).sum()))
    return dict(chi2=float(2 * b.cost), costs_by_start=costs, params={k: (v if isinstance(v, list) else float(v)) for k, v in d.items()}, parts=parts,
                a2_egfr_slope={a: float(np.polyfit(yrs, mo[a][0][12:, 1], 1)[0]) for a in ARMS}, D36={a: mo[a][1] for a in ARMS})
OUTF = 'results/r79_protection_route.json'; out = json.load(open(OUTF)) if os.path.exists(OUTF) else {}
for tag, cf, kt in (('Mn', False, 0.0), ('M0', True, 0.0), ('M1', False, 'free'), ('M2', True, 'free')):
    if tag in out: continue
    b, costs = fit(cf, kt, range(8)); out[tag] = summarise(b, cf, kt, costs); out[tag]['x'] = [float(z) for z in b.x]; print(tag, json.dumps({k: v for k, v in out[tag].items() if k != 'x'}), flush=True); json.dump(out, open(OUTF, 'w'), indent=1)
if 'profile' not in out:
    prof = []; warm = None
    for lk in np.linspace(-7, -3, 9):
        b, _ = fit(True, float(10 ** lk), range(20, 23), warm); warm = b.x; d = unpack(b.x, True, float(10 ** lk)); prof.append(dict(log10_ktox=float(lk), chi2=float(2 * b.cost), cni=float(d['cni']), kdam=float(d['kdam']), leak=float(d['leak']))); print('profile', prof[-1], flush=True)
    out['profile'] = prof; json.dump(out, open(OUTF, 'w'), indent=1)
c = np.array([q['chi2'] for q in out['profile']]); dof = 68 - 12; red = max(c.min() / dof, 1.0); ok = (c - c.min()) / red <= 3.84
v = dict(H1=bool(out['Mn']['chi2'] - min(out['M0']['chi2'], out['M1']['chi2']) > 3.84), H2=bool(abs(out['M0']['chi2'] - out['M1']['chi2']) < 3.84), H3=bool(ok[0]), H4=bool(out['M1']['params']['ktox'] > 3e-5))
out['verdict'] = v; out['profile_interval_scaled'] = dict(lo=float(np.linspace(-7, -3, 9)[ok].min()), hi=float(np.linspace(-7, -3, 9)[ok].max()), reduced_chi2=float(c.min() / dof)); print(v, out['profile_interval_scaled']); json.dump(out, open(OUTF, 'w'), indent=1)
