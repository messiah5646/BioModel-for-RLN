"""Round 78b -- RECHECK of round 78, requested by the user. Exploratory (not pre-registered); purpose is to try to break the result.
1 trivial baselines recomputed from the CSV, not copied from R50      2 saved parameters re-simulated: do the reported numbers reproduce
3 optimiser stability: holdout refit from new random starts           4 attack: drop week 2 from the fit -- does the step survive
5 second split never used before: fit weeks 2-36, predict 42-52       6 parameter-count penalty (AIC on chi2)"""
import json, numpy as np
ns = {}; exec(open('r78_shared_step.py').read().split("ALL = np.arange")[0], ns)
O1, ARMS, model, resid, unpack, least_squares = ns['O1'], ns['ARMS'], ns['model'], ns['resid'], ns['unpack'], ns['least_squares']
W = ns['W']; saved = json.load(open('results/r78_shared_step.json')); out = {}
def trivial(last, prev, outidx):
    e_p, e_l = [], []
    for a in ARMS:
        for col in (0, 2):
            for j in outidx:
                e_p.append(O1[a][last, col] - O1[a][j, col])
                sl = (O1[a][last, col] - O1[a][prev, col]) / (W[last] - W[prev]); e_l.append(O1[a][last, col] + sl * (W[j] - W[last]) - O1[a][j, col])
    return float(np.sqrt(np.mean(np.square(e_p)))), float(np.sqrt(np.mean(np.square(e_l))))
def rmse(mo, outidx):
    e = []
    for a in ARMS: e += list(mo[a][outidx, 0] - O1[a][outidx, 0]) + list(mo[a][outidx, 1] - O1[a][outidx, 2])
    return float(np.sqrt(np.mean(np.square(e))))
def fit(steps, idx, seeds):
    best = None; costs = []
    for s in seeds:
        g = np.random.default_rng(s); x0 = [np.log(g.uniform(0.2, 0.95)), np.log(10 ** g.uniform(-5, -3)), np.log(10 ** g.uniform(-8, -5)), g.normal(-1, 1), g.normal(0, 1), g.normal(0, 1)] + ([g.normal(-0.7, 0.4), g.normal(3, 1.5)] if steps else [])
        try: r = least_squares(resid, x0, args=(steps, idx), method='lm', max_nfev=500)
        except Exception: continue
        costs.append(float(2 * r.cost))
        if best is None or r.cost < best.cost: best = r
    return best, costs
OUT = np.array([7, 8, 9, 10]); FIT = np.arange(7)
# 1
out['1_trivial_recomputed'] = dict(zip(('persistence', 'linear'), trivial(6, 4, OUT))); print(1, out['1_trivial_recomputed'], flush=True)
# 2
def xfrom(d, steps):
    lg = lambda z: float(np.log(max(z, 1e-300) / max(1 - z, 1e-300)))
    return [np.log(max(d['s_inf'], 1e-300)), np.log(d['kdam']), np.log(max(d['ktox'], 1e-300)), lg(d['cni']), lg(d['leak']), lg(d['h'] / 0.15)] + ([d['su'], d['sgf']] if steps else [])
rep = {}
for tag, steps in (('none', False), ('steps', True)):
    f = saved['fits'][tag]; chi = float((resid(np.array(xfrom(f['params'], steps)), steps, np.arange(12)) ** 2).sum()); hr = rmse(model(np.array(xfrom(f['holdout_fit_params'], steps)), steps), OUT)
    rep[tag] = dict(chi2_reported=f['chi2'], chi2_recomputed=chi, holdout_reported=f['holdout_rmse'], holdout_recomputed=hr)
out['2_reproduce'] = rep; print(2, rep, flush=True)
# 3
b, costs = fit(True, FIT, range(100, 108)); out['3_new_starts'] = dict(holdout_rmse=rmse(model(b.x, True), OUT), fit_chi2_by_start=sorted(costs), step=[float(b.x[6]), float(b.x[7])]); print(3, out['3_new_starts'], flush=True)
# 4
b, _ = fit(True, np.arange(1, 7), range(200, 206)); out['4_drop_week2'] = dict(holdout_rmse=rmse(model(b.x, True), OUT), step=[float(b.x[6]), float(b.x[7])], week2_pred={a: [float(model(b.x, True)[a][0, 0]), float(model(b.x, True)[a][0, 1])] for a in ARMS}, week2_obs={a: [float(O1[a][0, 0]), float(O1[a][0, 2])] for a in ARMS}); print(4, out['4_drop_week2'], flush=True)
# 5
FIT2 = np.arange(9); OUT2 = np.array([9, 10, 11]); tp, tl = trivial(8, 7, OUT2); s5 = dict(persistence=tp, linear=tl)
for tag, steps in (('none', False), ('steps', True)):
    b, _ = fit(steps, FIT2, range(300, 306)); s5[tag] = rmse(model(b.x, steps), OUT2)
out['5_second_split_fit_2to36_predict_42to52'] = s5; print(5, s5, flush=True)
# 6
out['6_aic'] = dict(none=saved['fits']['none']['chi2'] + 2 * 6, steps=saved['fits']['steps']['chi2'] + 2 * 8); print(6, out['6_aic'])
json.dump(out, open('results/r78b_recheck.json', 'w'), indent=1)
