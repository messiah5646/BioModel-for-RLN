"""Backtest gate for the harness. LESSON 6 -- before trusting a simulator on a NEW drug, run it on a drug whose
trial result you already know, and compare. Here: voclosporin vs AURORA 1 treated arm (weeks 2-52) and the
AURORA 2 eGFR slope (months 12-36). Honest label: the voclosporin strengths were fitted to AURORA in round 78,
so year 1 is IN-SAMPLE. The AURORA 2 slope was not used to fit them (round 78 already showed it is missed)."""
import csv, json, os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import vtrial as V
rows = lambda f: list(csv.DictReader([l for l in open(os.path.join(HERE, '..', '..', 'data', 'real', f)) if not l.lstrip().startswith('#')]))
A1 = rows('aurora1_per_visit.csv'); A2 = rows('aurora2_per_visit.csv')
T1 = np.array([[float(r['upcr_change']), float(r['upcr_sd']) / np.sqrt(179), float(r['egfr_change']), float(r['egfr_sd']) / np.sqrt(179)] for r in A1 if r['arm'] == 'treated'])
P1 = V.load_aurora()
def slope(arm):
    m = np.array([12, 18, 24, 30, 36]) / 12.0; e = np.array([[float(r['egfr_change']), float(r['egfr_se'])] for r in A2 if r['arm'] == arm])
    w = 1 / e[:, 1] ** 2; X = np.c_[np.ones(5), m]; C = np.linalg.inv(X.T @ (w[:, None] * X)); return float((C @ X.T @ (w * e[:, 0]))[1]), float(np.sqrt(C[1, 1]))
ens = json.load(open(os.path.join(HERE, 'runs', 'ensemble_v2_boot8.json')))['kept']; spec = json.load(open(os.path.join(HERE, 'drugs', 'voclosporin_like.json')))
weeks = V.FIT_WEEKS + V.A2_WEEKS   # months 12..36 (bug fixed: an earlier version used months 36..108)
runs = V.run_trial(spec, ens, [round(w, 2) for w in weeks])
T = np.array([r['drug'] for r in runs]); C = np.array([r['control'] for r in runs])
chi_t = [float((((t[:12, 0] - T1[:, 0]) / T1[:, 1]) ** 2).sum() + (((t[:12, 1] - T1[:, 2]) / T1[:, 3]) ** 2).sum()) for t in T]
yrs = np.array([12, 18, 24, 30, 36]) / 12.0
ms = {'placebo': [float(np.polyfit(yrs, c[12:, 1], 1)[0]) for c in C], 'treated': [float(np.polyfit(yrs, t[12:, 1], 1)[0]) for t in T]}
out = dict(year1_treated_chi2_median=float(np.median(chi_t)), year1_treated_chi2_range=[min(chi_t), max(chi_t)], n_points=24,
           wk52_diff_model=dict(dUPCR=float(np.median(T[:, 11, 0] - C[:, 11, 0])), dEGFR=float(np.median(T[:, 11, 1] - C[:, 11, 1]))),
           wk52_diff_observed=dict(dUPCR=float(T1[11, 0] - P1[11, 0]), dEGFR=float(T1[11, 2] - P1[11, 2])),
           a2_slope_observed={a: slope(a) for a in ('placebo', 'treated')},
           a2_slope_model={a: [float(np.median(v)), float(min(v)), float(max(v))] for a, v in ms.items()})
out['a2_slope_z'] = {a: (out['a2_slope_model'][a][0] - out['a2_slope_observed'][a][0]) / out['a2_slope_observed'][a][1] for a in ('placebo', 'treated')}
print(json.dumps(out, indent=1)); json.dump(out, open(os.path.join(HERE, 'runs', 'backtest_voclosporin.json'), 'w'), indent=1)
