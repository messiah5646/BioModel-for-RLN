"""Round 83 -- prereg: project doc claude/prereg-round83.md.
LESSON 8 -- 'simulate a real drug from its properties, then compare with the trial' is only a fair test if the
properties come from somewhere OTHER than that trial's outcomes. Here the property is the dose ratio, and the
drug's base strengths were fitted on a different trial (AURORA)."""
import csv, json, os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import vtrial as V
R = lambda f: list(csv.DictReader([l for l in open(os.path.join(HERE, '..', '..', 'data', 'real', f)) if not l.lstrip().startswith('#')]))
LV = {(r['arm'], int(r['week'])): r for r in R('aura_lv_per_visit.csv')}
ens = json.load(open(os.path.join(HERE, 'runs', 'ensemble_v2_boot8.json')))['kept']
base = json.load(open(os.path.join(HERE, 'drugs', 'voclosporin_like.json')))
ARMS = {'placebo': ('placebo', None), 'low_23.7mg': ('voclosporin_low_23.7mg_bid', 1.0), 'high_39.5mg': ('voclosporin_high_39.5mg_bid', 39.5 / 23.7)}
def spec(scale): return dict(name=f'voclosporin x{scale:.3f}', effects=[dict(mechanism=e['mechanism'], strength=min(0.95, e['strength'] * scale)) for e in base['effects']])
out = {}
for a, (key, scale) in ARMS.items():
    b0 = float(LV[(key, 0)]['upcr_absolute']); preds = {24: [], 48: []}
    for m in ens:
        p = V.base_params(m['s_inf'], m['k_dam'], m['k_tox']); p['p0'] = b0
        q = V.apply_drug(p, spec(scale)) if scale else p
        d, _ = V.trajectory(q, [24, 48], m['step'])
        for i, w in enumerate((24, 48)): preds[w].append(float(d[i, 0]))
    rec = {}
    for w in (24, 48):
        r = LV[(key, w)]; o = float(r['upcr_change']); se = float(r['upcr_sd']) / np.sqrt(float(r['upcr_n'])); pm = float(np.median(preds[w]))
        rec[f'wk{w}'] = dict(pred=pm, pred_range=[min(preds[w]), max(preds[w])], obs=o, se=se, err=pm - o, z=(pm - o) / se,
                             pred_frac=pm / b0, obs_frac=o / b0)
    out[a] = dict(baseline=b0, **rec)
lo, hi = out['low_23.7mg'], out['high_39.5mg']
v = dict(H1=all(hi[f'wk{w}']['pred'] < lo[f'wk{w}']['pred'] for w in (24, 48)),
         H2=all(hi[f'wk{w}']['obs'] > lo[f'wk{w}']['obs'] for w in (24, 48)),
         H3=all(hi[f'wk{w}']['obs_frac'] > lo[f'wk{w}']['obs_frac'] for w in (24, 48)),
         H4=all(abs(hi[f'wk{w}']['z']) > 2 for w in (24, 48)))
out['verdict'] = v; out['rmse_all_6'] = float(np.sqrt(np.mean([out[a][f'wk{w}']['err'] ** 2 for a in ARMS for w in (24, 48)])))
json.dump(out, open(os.path.join(HERE, 'runs', 'existing_drugs_check.json'), 'w'), indent=1)
for a in ARMS:
    print(f"{a:12s} baseline {out[a]['baseline']:.2f} | " + ' | '.join(f"wk{w}: sim {out[a][f'wk{w}']['pred']:+.2f} real {out[a][f'wk{w}']['obs']:+.2f} (z {out[a][f'wk{w}']['z']:+.1f}; frac sim {out[a][f'wk{w}']['pred_frac']:.2f} real {out[a][f'wk{w}']['obs_frac']:.2f})" for w in (24, 48)))
print('RMSE over 6 arm-visits', round(out['rmse_all_6'], 3)); print(v)
