"""Round 82 attack (pre-registration section 'Interpretation': a win over C2 must be attacked before reporting).
Attack A: the primary variant ALSO rescaled the fitted day-0 UPCR step by baseline ratio. The pre-registration said
'p0 set to each arm's baseline (the only change)'. So the step rescale was an unregistered extra; rerun without it.
Attack B: how big is the margin compared with measurement noise? Monte Carlo over the observed SEs."""
import json, os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
src = open(os.path.join(HERE, 'external_aura_lv.py')).read().split("P = {'harness_primary'")[0]; ns = {'__file__': os.path.join(HERE, 'external_aura_lv.py')}; exec(src, ns)
V, ens, spec, obs, ARM, BASE_LV = ns['V'], ns['ens'], ns['spec'], ns['obs'], ns['ARM'], ns['BASE_LV']
pred = {k: [] for k in obs}
for m in ens:
    for a in ARM:
        p = V.base_params(m['s_inf'], m['k_dam'], m['k_tox']); p['p0'] = BASE_LV[a]
        d, _ = V.trajectory(V.apply_drug(p, spec) if a == 'treated' else p, [24, 48], m['step'])
        for i, w in enumerate((24, 48)): pred[(a, w)].append(float(d[i, 0]))
med = {k: float(np.median(v)) for k, v in pred.items()}
o = json.load(open(os.path.join(HERE, 'runs', 'external_aura_lv.json')))
err = lambda pr: np.array([pr[k] - obs[k][0] for k in obs])
prim = {k: o['harness_primary']['points'][f'{k[0]}_wk{k[1]}']['pred'] for k in obs}
c2 = {k: o['C2_copy_AURORA_scaled']['points'][f'{k[0]}_wk{k[1]}']['pred'] for k in obs}
rm = lambda e: float(np.sqrt(np.mean(e ** 2)))
A = dict(pred={f'{a}_wk{w}': v for (a, w), v in med.items()}, rmse_arms=rm(err(med)),
         diffs={f'wk{w}': med[('treated', w)] - med[('placebo', w)] for w in (24, 48)})
# Attack B: resample the observed means within their SEs; how often does the primary harness beat C2?
rng = np.random.default_rng(0); se = np.array([obs[k][1] for k in obs]); y = np.array([obs[k][0] for k in obs])
P1, P2, PA = (np.array([d[k] for k in obs]) for d in (prim, c2, med)); wins = winsA = 0
for _ in range(20000):
    yy = y + rng.normal(0, se); wins += rm(P1 - yy) < rm(P2 - yy); winsA += rm(PA - yy) < rm(P2 - yy)
B = dict(share_primary_beats_C2=wins / 20000, share_unscaled_step_beats_C2=winsA / 20000)
out = dict(attack_A_no_step_rescale=A, attack_B_noise=B, C2_rmse=o['C2_copy_AURORA_scaled']['rmse_arms'], primary_rmse=o['harness_primary']['rmse_arms'])
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'runs', 'external_aura_lv_attack.json'), 'w'), indent=1)
