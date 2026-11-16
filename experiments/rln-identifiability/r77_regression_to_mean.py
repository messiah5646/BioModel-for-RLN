"""Round 77 -- prereg: project doc claude/prereg-round77.md (written before this ran)."""
import json, signal, numpy as np
from scipy.optimize import least_squares
rng = np.random.default_rng(0); N = 400000; Z = rng.standard_normal((4, N))
def sel_stats(mu, sT, se, thr=1.5):
    T = mu + sT * Z[0]; S, B0, F = (np.exp(T + se * Z[i]) for i in (1, 2, 3)); k = S >= thr; base = (S[k] + B0[k]) / 2
    return base.mean(), base.std(), F[k].mean() - base.mean(), k.mean()
A = []
for cv in (0.2, 0.3, 0.4, 0.5, 0.6):
    se = np.sqrt(np.log(1 + cv ** 2))
    for sd in (2.0, 2.5, 3.0):
        f = lambda x: [sel_stats(x[0], abs(x[1]), se)[0] - 3.87, sel_stats(x[0], abs(x[1]), se)[1] - sd]
        r = least_squares(f, [1.0, 0.5]); m, s, rtm, frac = sel_stats(r.x[0], abs(r.x[1]), se)
        A.append(dict(cv_within=cv, baseline_sd=sd, fit_mean=m, fit_sd=s, rtm_upcr=rtm, share_of_week2=rtm / -0.70, frac_selected=frac, reliability=float(r.x[1] ** 2 / (r.x[1] ** 2 + se ** 2))))
        print(f"cv {cv} sd {sd}: selected mean {m:.2f} sd {s:.2f} | RTM {rtm:+.3f} mg/mg = {100*rtm/-0.70:.0f}% of week-2 fall | selected {100*frac:.0f}% | log-reliability {A[-1]['reliability']:.2f}")
E_ = []
for cvw in (0.05, 0.08):
    for sdb in (20.0, 30.0):
        swi = cvw * 78; f = lambda x: [(lambda T: ((T + swi * Z[1])[(T + swi * Z[1]) > 45]).mean() - 77.4)(x[0] + sdb * Z[0])]
        r = least_squares(f, [78.0]); T = r.x[0] + sdb * Z[0]; S = T + swi * Z[1]; k = S > 45; rt = (T + swi * Z[3])[k].mean() - ((S + T + swi * Z[2]) / 2)[k].mean()
        E_.append(dict(cv_within=cvw, between_sd=sdb, rtm_egfr=float(rt))); print(f"eGFR cv {cvw} between-sd {sdb}: RTM {rt:+.2f} mL/min")
# ---- Test B: free steps on the placebo arm
import csv
from model import make_params, simulate, observe
D = list(csv.DictReader([l for l in open('../../data/real/aurora1_per_visit.csv') if not l.lstrip().startswith('#')]))
W = [2, 4, 8, 12, 16, 20, 24, 30, 36, 42, 48, 52]; TG = np.array([0.0] + [w * 7.0 for w in W]); n = 178
O = np.array([[float(r['upcr_change']), float(r['upcr_sd']) / np.sqrt(n), float(r['egfr_change']), float(r['egfr_sd']) / np.sqrt(n)] for r in D if r['arm'] == 'placebo'])
class TO(Exception): pass
def _al(*a): raise TO()
signal.signal(signal.SIGALRM, _al)
def curve(s_inf, kdam):
    p = make_params(0.0); p['p0'] = 3.87; p['G0'] = 77.4 / (1 - p['c_inf']); p['s_inf'] = s_inf; p['k_dam'] = kdam; p['k_tox'] = 1e-7
    Y = simulate(p, TG, rtol=1e-8, atol=1e-10); Ob = np.array([observe(y, p, t) for t, y in zip(TG, Y)]); return Ob[1:, :2] - Ob[0, :2]
def resid(x, steps):
    signal.setitimer(signal.ITIMER_REAL, 5.0)
    try: C = curve(np.exp(min(x[0], 0)), np.exp(x[1]))
    except Exception: signal.setitimer(signal.ITIMER_REAL, 0); return np.full(24, 1e3)
    signal.setitimer(signal.ITIMER_REAL, 0); su, sg = (x[2], x[3]) if steps else (0.0, 0.0)
    return np.r_[(C[:, 0] + su - O[:, 0]) / O[:, 1], (C[:, 1] + sg - O[:, 2]) / O[:, 3]]
B = {}
for steps in (False, True):
    best = None
    for s in range(6):
        g = np.random.default_rng(s); x0 = [np.log(g.uniform(0.2, 0.95)), np.log(10 ** g.uniform(-5, -3))] + ([g.normal(-0.5, 0.5), g.normal(2, 2)] if steps else [])
        r = least_squares(resid, x0, args=(steps,), method='lm', max_nfev=600)
        if best is None or r.cost < best.cost: best = r
    rr = resid(best.x, steps); C = curve(np.exp(min(best.x[0], 0)), np.exp(best.x[1]))
    B['steps' if steps else 'none'] = dict(chi2=float(2 * best.cost), chi2_upcr=float((rr[:12] ** 2).sum()), chi2_egfr=float((rr[12:] ** 2).sum()), s_inf=float(np.exp(min(best.x[0], 0))), k_dam=float(np.exp(best.x[1])),
                                         step_upcr=float(best.x[2]) if steps else 0.0, step_egfr=float(best.x[3]) if steps else 0.0, model_upcr=[round(float(v), 2) for v in C[:, 0]], model_egfr=[round(float(v), 1) for v in C[:, 1]])
    print('steps' if steps else 'none', {k: (round(v, 3) if isinstance(v, float) else v) for k, v in B['steps' if steps else 'none'].items()})
bound05 = min(a['rtm_upcr'] for a in A if a['cv_within'] == 0.5)
v = dict(H1=all(abs(a['rtm_upcr']) < 0.35 for a in A if a['cv_within'] <= 0.5), H2=all(e['rtm_egfr'] < 0.5 for e in E_), H3=B['steps']['chi2'] > 0.5 * B['none']['chi2'], H4=abs(B['steps']['step_upcr']) > abs(bound05), bound_cv05=bound05)
v = {k: (bool(x) if isinstance(x, (bool, np.bool_)) else float(x)) for k, x in v.items()}; print(v); json.dump(dict(testA_upcr=A, testA_egfr=E_, testB=B, verdict=v), open('results/r77_regression_to_mean.json', 'w'), indent=1)
