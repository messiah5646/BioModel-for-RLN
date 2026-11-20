"""Round 82 -- prereg: project doc claude/prereg-round82.md. External accuracy test on AURA-LV (real, unseen by the harness).
LESSON 7 -- the only accuracy that counts is on data the model never touched. AURORA was used to build the harness,
so accuracy on AURORA flatters it. AURA-LV is a different trial of the same drug: a fair, if small, exam."""
import csv, json, os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import vtrial as V
R = lambda f: list(csv.DictReader([l for l in open(os.path.join(HERE, '..', '..', 'data', 'real', f)) if not l.lstrip().startswith('#')]))
LV = {(r['arm'], int(r['week'])): r for r in R('aura_lv_per_visit.csv')}
A1 = {(r['arm'], int(r['week'])): float(r['upcr_change']) for r in R('aurora1_per_visit.csv')}
ARM = {'placebo': ('placebo', 'placebo'), 'treated': ('voclosporin_low_23.7mg_bid', 'treated')}
BASE_LV = {'placebo': float(LV[('placebo', 0)]['upcr_absolute']), 'treated': float(LV[('voclosporin_low_23.7mg_bid', 0)]['upcr_absolute'])}
BASE_A1 = {'placebo': 3.87, 'treated': 4.14}
obs = {(a, w): (float(LV[(ARM[a][0], w)]['upcr_change']), float(LV[(ARM[a][0], w)]['upcr_sd']) / np.sqrt(float(LV[(ARM[a][0], w)]['upcr_n']))) for a in ARM for w in (24, 48)}
ens = json.load(open(os.path.join(HERE, 'runs', 'ensemble_v2_boot8.json')))['kept']; spec = json.load(open(os.path.join(HERE, 'drugs', 'voclosporin_like.json')))
def harness(rescale):
    pred = {k: [] for k in obs}
    for m in ens:
        for a in ARM:
            p = V.base_params(m['s_inf'], m['k_dam'], m['k_tox'])
            if rescale: p['p0'] = BASE_LV[a]
            elif a == 'treated': p['p0'] = BASE_A1['treated']
            q = V.apply_drug(p, spec) if a == 'treated' else p
            d, _ = V.trajectory(q, [24, 48], m['step'])
            # the day-0 step was fitted on AURORA's UPCR scale; rescale it with baseline in the primary variant
            if rescale: d[:, 0] += m['step'][0] * (BASE_LV[a] / 3.87 - 1)
            elif a == 'treated': d[:, 0] += m['step'][0] * (BASE_A1['treated'] / 3.87 - 1)
            for i, w in enumerate((24, 48)): pred[(a, w)].append(float(d[i, 0]))
    return {k: (float(np.median(v)), float(min(v)), float(max(v))) for k, v in pred.items()}
P = {'harness_primary': harness(True), 'harness_secondary': harness(False)}
P['C1_copy_AURORA'] = {(a, w): (A1[(ARM[a][1], w)],) * 3 for a in ARM for w in (24, 48)}
P['C2_copy_AURORA_scaled'] = {(a, w): (A1[(ARM[a][1], w)] * BASE_LV[a] / BASE_A1[a],) * 3 for a in ARM for w in (24, 48)}
out = {'observed': {f'{a}_wk{w}': dict(change=v[0], se=v[1]) for (a, w), v in obs.items()}}
for name, pr in P.items():
    pts = {f'{a}_wk{w}': dict(pred=pr[(a, w)][0], range=pr[(a, w)][1:], err=pr[(a, w)][0] - obs[(a, w)][0], z=(pr[(a, w)][0] - obs[(a, w)][0]) / obs[(a, w)][1]) for (a, w) in obs}
    diffs = {}
    for w in (24, 48):
        dp = pr[('treated', w)][0] - pr[('placebo', w)][0]; do = obs[('treated', w)][0] - obs[('placebo', w)][0]; se = np.hypot(obs[('treated', w)][1], obs[('placebo', w)][1])
        diffs[f'wk{w}'] = dict(pred=dp, obs=do, err=dp - do, z=(dp - do) / se, sign_ok=bool(np.sign(dp) == np.sign(do)))
    out[name] = dict(points=pts, rmse_arms=float(np.sqrt(np.mean([p['err'] ** 2 for p in pts.values()]))), diffs=diffs, rmse_diffs=float(np.sqrt(np.mean([d['err'] ** 2 for d in diffs.values()]))))
out['C0_zero_difference'] = dict(rmse_diffs=float(np.sqrt(np.mean([(out['harness_primary']['diffs'][k]['obs']) ** 2 for k in ('wk24', 'wk48')]))))
h, c1, c2 = out['harness_primary'], out['C1_copy_AURORA'], out['C2_copy_AURORA_scaled']
out['verdict'] = dict(H1=h['rmse_arms'] < c1['rmse_arms'], H2=h['rmse_arms'] > c2['rmse_arms'],
                      H3=all(d['sign_ok'] and abs(d['z']) <= 2 for d in h['diffs'].values()), H4=h['rmse_diffs'] < out['C0_zero_difference']['rmse_diffs'])
json.dump(out, open(os.path.join(HERE, 'runs', 'external_aura_lv.json'), 'w'), indent=1)
print('observed', {k: round(v['change'], 2) for k, v in out['observed'].items()}, 'SE', {k: round(v['se'], 2) for k, v in out['observed'].items()})
for n in ('harness_primary', 'harness_secondary', 'C1_copy_AURORA', 'C2_copy_AURORA_scaled'):
    o = out[n]; print(f"{n:24s} arms RMSE {o['rmse_arms']:.3f} | diffs RMSE {o['rmse_diffs']:.3f} |", {k: f"{v['pred']:.2f} (z {v['z']:+.1f})" for k, v in o['points'].items()}, {k: f"{v['pred']:.2f} vs {v['obs']:.2f} z {v['z']:+.1f}" for k, v in o['diffs'].items()})
print('C0 diffs RMSE', round(out['C0_zero_difference']['rmse_diffs'], 3)); print(out['verdict'])
