"""P3 / Q58 analysis (push 3) per sle/PREREG-P3.md. Same estimators as run86.py.
Usage: python3 -I run_p3.py sri4_table.csv bicla_table.csv"""
import sys, json
import numpy as np, pandas as pd

def load(path):
    R = pd.read_csv(path)
    a, b, c, d = R.r_drug + .5, R.n_drug - R.r_drug + .5, R.r_ctrl + .5, R.n_ctrl - R.r_ctrl + .5
    R['lor'] = np.log(a * d / (b * c)); R['se'] = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d); return R
def iv(x, se): w = 1 / se ** 2; return float((w * x).sum() / w.sum())
def dl(x, se):
    w = 1 / se ** 2; m = (w * x).sum() / w.sum(); Q = (w * (x - m) ** 2).sum(); k = len(x)
    return max(0.0, (Q - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))
def primary(R): return R[R.week <= 52].sort_values('week').groupby(['trial', 'arm']).tail(1).reset_index(drop=True)
def score(df, m):
    e = df.obs - df[m]; w = 1 / df.se ** 2
    return round(float((w * e ** 2).sum() / w.sum()), 4), round(float((e ** 2 - df.se ** 2).mean()), 4)

def task_d(P, unit):
    rows = []
    for g in P[unit].unique():
        T, H = P[P[unit] != g], P[P[unit] == g]
        mu = iv(T.lor, T.se); t2 = dl(T.lor.values, T.se.values)
        for _, r in H.iterrows():
            anc = T[T.molecule == r.molecule]; lvl = 'molecule'
            if not len(anc): anc = T[T.target == r.target]; lvl = 'target'
            if len(anc):
                am = iv(anc.lor, anc.se); ase = 1 / np.sqrt((1 / anc.se ** 2).sum()); w = t2 / (t2 + ase ** 2)
                d3 = mu + w * (am - mu)
            else: d3, lvl = mu, 'none'
            rows.append(dict(trial=r.trial, arm=r.arm, molecule=r.molecule, anchor=lvl, obs=r.lor, se=r.se, D0=0.0, D1=mu, D3=d3, tau2=t2))
    return pd.DataFrame(rows)

def verdict(X):
    s0, s1, s3 = score(X, 'D0')[0], score(X, 'D1')[0], score(X, 'D3')[0]
    if s3 <= 0.9 * s0 and s3 <= 0.9 * s1: v = 'PASS (hierarchy generalises)'
    elif s3 >= s1: v = 'FAIL (D3 does not beat pooled mean)'
    else: v = 'PARTIAL'
    return dict(n_arms=len(X), n_trials=X.trial.nunique(), D0=score(X, 'D0'), D1=score(X, 'D1'), D3=score(X, 'D3'),
                D3_vs_D1=round(s3 / s1 - 1, 3), D3_vs_D0=round(s3 / s0 - 1, 3), verdict=v)

def task_c(R, P):
    C = []
    for (tr, arm), g in R.groupby(['trial', 'arm']):
        e = g[(g.week >= 12) & (g.week <= 26)]; l = g[(g.week >= 48) & (g.week <= 52)]
        if len(e) and len(l):
            e = e.sort_values('week').iloc[-1]; l = l.sort_values('week').iloc[-1]
            C.append(dict(trial=tr, arm=arm, cls=g.cls.iloc[0], te=e.week, tl=l.week, de=e.lor, see=e.se, dl=l.lor, sel=l.se))
    C = pd.DataFrame(C)
    if not len(C): return None, C
    pool = iv(P.lor, P.se); t2 = dl(P.lor.values, P.se.values)
    C['C0'] = C.de; C['C1'] = 0.0; C['C2'] = pool + (t2 / (t2 + C.see ** 2)) * (C.de - pool)
    def sc(df):
        out = {}
        for m in ['C0', 'C1', 'C2']:
            e = df.dl - df[m]; w = 1 / df.sel ** 2; out[m] = round(float((w * e ** 2).sum() / w.sum()), 4)
        out['n'] = len(df); out['gain_C2_vs_C1'] = round(1 - out['C2'] / out['C1'], 3) if out['C1'] > 0 else None
        out['gain_best_vs_C1'] = round(1 - min(out['C0'], out['C2']) / out['C1'], 3) if out['C1'] > 0 else None
        return out
    res = dict(overall=sc(C), tau2=round(t2, 4))
    res['overall']['verdict'] = 'early read-out fast-forwards' if res['overall']['gain_best_vs_C1'] >= 0.10 else 'does not fast-forward'
    for k in ['fast', 'slow']:
        res[k] = sc(C[C.cls == k]) if (C.cls == k).any() else None
    f, s = res['fast'], res['slow']
    if f and s and f['gain_C2_vs_C1'] is not None and s['gain_C2_vs_C1'] is not None:
        res['split_verdict'] = 'Supported' if (f['gain_C2_vs_C1'] > s['gain_C2_vs_C1'] and f['gain_C2_vs_C1'] > 0) else 'not supported'
    else: res['split_verdict'] = 'not testable (a class has no arms)'
    return res, C

out = {}
R = load(sys.argv[1]); P = primary(R)
Xt = task_d(P, 'trial'); Xm = task_d(P, 'molecule')
out['SRI4_LOTO_primary_anchored'] = verdict(Xt[Xt.anchor != 'none'])
out['SRI4_LOTO_all_arms'] = verdict(Xt)
out['SRI4_LOMO_anchored'] = verdict(Xm[Xm.anchor != 'none'])
out['SRI4_LOMO_all_arms'] = verdict(Xm)
out['anchor_levels_LOTO'] = Xt.anchor.value_counts().to_dict()
# one arm per trial (highest n_drug arm as the 'highest dose' proxy is ambiguous; use highest numeric dose in arm title)
def dose(a):
    import re
    m = re.findall(r'(\d+(?:\.\d+)?) ?(?:mg|mcg|µg)', a); return max(map(float, m)) if m else 0
P1 = P.assign(dz=P.arm.map(dose)).sort_values(['trial', 'dz']).groupby('trial').tail(1)
X1 = task_d(P1, 'trial'); out['SRI4_LOTO_one_arm_per_trial_anchored'] = verdict(X1[X1.anchor != 'none'])
out['Q58'], C = task_c(R, P)
# BICLA robustness
Rb = load(sys.argv[2]); Pb = primary(Rb); Xb = task_d(Pb, 'trial')
out['BICLA_LOTO_anchored'] = verdict(Xb[Xb.anchor != 'none']) if (Xb.anchor != 'none').any() else 'no anchored arms'
out['BICLA_LOTO_all_arms'] = verdict(Xb)
json.dump(out, open('p3_result.json', 'w'), indent=1)
Xt.to_csv('p3_taskD_loto.csv', index=False); Xm.to_csv('p3_taskD_lomo.csv', index=False); C.to_csv('p3_taskC.csv', index=False)
print(json.dumps(out, indent=1))
print(Xt[Xt.anchor != 'none'].round(3).to_string())
print(C.round(3).to_string())
