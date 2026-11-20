"""Round 84 -- prereg: project doc claude/prereg-round84.md.
LESSON 9 -- when a simple empirical model keeps beating a mechanistic one, build the empirical model properly
and make it the default. Model-based meta-analysis (MBMA) pools trials on a common scale (here: fraction of baseline
proteinuria remaining) and is tested by leaving one whole trial out at a time."""
import csv, json, os, sys, numpy as np
from scipy.optimize import least_squares
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import vtrial as V
R = lambda f: list(csv.DictReader([l for l in open(os.path.join(HERE, '..', '..', 'data', 'real', f)) if not l.lstrip().startswith('#')]))
def pts_aurora(arm, base, n):
    return [(int(r['week']), 1 + float(r['upcr_change']) / base, float(r['upcr_sd']) / np.sqrt(n) / base) for r in R('aurora1_per_visit.csv') if r['arm'] == arm]
def pts_lv(arm):
    rows = {int(r['week']): r for r in R('aura_lv_per_visit.csv') if r['arm'] == arm}; b = float(rows[0]['upcr_absolute'])
    return [(w, 1 + float(rows[w]['upcr_change']) / b, float(rows[w]['upcr_sd']) / np.sqrt(float(rows[w]['upcr_n'])) / b) for w in (24, 48)]
PL = {'AURORA1': pts_aurora('placebo', 3.87, 178), 'AURA-LV': pts_lv('placebo'),
      'CFZ533': [(round(int(r['day']) / 7), float(r['upcr_ratio']), float(r['upcr_sd']) / np.sqrt(float(r['upcr_n']))) for r in R('cfz533_per_visit.csv') if r['arm'] == 'placebo']}
TR = {'AURORA1': pts_aurora('treated', 4.14, 179), 'AURA-LV': pts_lv('voclosporin_low_23.7mg_bid')}
curve = lambda x, t: x[0] + (1 - x[0]) * np.exp(-np.asarray(t, float) / np.exp(x[1]))
def fit_curve(trials):
    P = [p for k in trials for p in PL[k]]; t, y, s = (np.array(c) for c in zip(*P))
    return least_squares(lambda x: (curve(x, t) - y) / s, [0.5, np.log(8.0)]).x
def copy_nearest(trials, w):
    return float(np.mean([min(PL[k], key=lambda p: abs(p[0] - w))[1] for k in trials]))
# mechanistic harness: fraction remaining for placebo and voclosporin-like arms, ensemble median
ens = json.load(open(os.path.join(HERE, 'runs', 'ensemble_v2_boot8.json')))['kept']; spec = json.load(open(os.path.join(HERE, 'drugs', 'voclosporin_like.json')))
weeks_all = sorted({p[0] for v in PL.values() for p in v}); mech = {'placebo': {}, 'treated': {}}
for arm, b in (('placebo', 3.87), ('treated', 4.14)):
    M = []
    for m in ens:
        p = V.base_params(m['s_inf'], m['k_dam'], m['k_tox']); p['p0'] = b
        d, _ = V.trajectory(V.apply_drug(p, spec) if arm == 'treated' else p, weeks_all, m['step']); M.append(1 + d[:, 0] / b)
    mech[arm] = dict(zip(weeks_all, np.median(np.array(M), 0)))
out = {'placebo': {}, 'drug': {}}
for held in PL:
    others = [k for k in PL if k != held]; x = fit_curve(others); rec = []
    for w, y, s in PL[held]:
        rec.append(dict(week=w, obs=y, se=s, E1=float(curve(x, w)), M=float(mech['placebo'][w]), C=copy_nearest(others, w)))
    rm = {k: float(np.sqrt(np.mean([(r[k] - r['obs']) ** 2 for r in rec]))) for k in ('E1', 'M', 'C')}
    out['placebo'][held] = dict(curve=dict(asymptote=float(x[0]), tau_weeks=float(np.exp(x[1]))), points=rec, rmse=rm,
                                max_absz={k: float(max(abs(r[k] - r['obs']) / r['se'] for r in rec)) for k in ('E1', 'M', 'C')})
for held, other in (('AURA-LV', 'AURORA1'), ('AURORA1', 'AURA-LV')):
    xo = fit_curve([k for k in PL if k != held]); rec = []
    for w in (24, 48):
        y, s = next((p[1], p[2]) for p in TR[held] if p[0] == w)
        tp = min(TR[other], key=lambda p: abs(p[0] - w)); pp = min(PL[other], key=lambda p: abs(p[0] - w)); ratio = tp[1] / pp[1]
        rec.append(dict(week=w, obs=y, se=s, E1_drug=float(curve(xo, w) * ratio), M=float(mech['treated'][w]), drug_ratio_from_other_trial=ratio))
    out['drug'][held] = dict(points=rec, rmse={k: float(np.sqrt(np.mean([(r[k] - r['obs']) ** 2 for r in rec]))) for k in ('E1_drug', 'M')})
P = out['placebo']; D = out['drug']
out['verdict'] = dict(H1=sum(P[k]['rmse']['E1'] < P[k]['rmse']['M'] for k in P) >= 2, H2=all(D[k]['rmse']['E1_drug'] < D[k]['rmse']['M'] for k in D),
                      H3=sum(P[k]['rmse']['E1'] < P[k]['rmse']['C'] for k in P) <= 1, H4=all(P['CFZ533']['max_absz'][k] > 2 for k in ('E1', 'M', 'C')))
json.dump(out, open(os.path.join(HERE, 'runs', 'mbma_loto.json'), 'w'), indent=1)
for k, v in P.items(): print('PLACEBO held-out', k, {m: round(e, 3) for m, e in v['rmse'].items()}, 'max|z|', {m: round(e, 1) for m, e in v['max_absz'].items()}, 'curve', {a: round(b, 2) for a, b in v['curve'].items()})
for k, v in D.items(): print('DRUG held-out', k, {m: round(e, 3) for m, e in v['rmse'].items()}, [(r['week'], round(r['obs'], 3), round(r['E1_drug'], 3), round(r['M'], 3)) for r in v['points']])
print(out['verdict'])
