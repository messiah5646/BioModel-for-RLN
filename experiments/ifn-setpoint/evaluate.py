"""Round 96 evaluation exactly as PREREG.md. Reads scores.csv, writes result.json."""
import json, numpy as np, pandas as pd
from scipy.stats import spearmanr
d = pd.read_csv('scores.csv')
late_tp = {'CD3xCD19': 'M3', 'CART': 'post'}; base_tp = {'CD3xCD19': 'Base', 'CART': 'pre'}
recs = []
for (tr, dn), g in d.groupby(['trial', 'donor']):
    b = g[g.tp == base_tp[tr]].score.iloc[0]
    for tp in g.tp:
        if tp == base_tp[tr]: continue
        recs.append(dict(trial=tr, donor=dn, tp=tp, base=b, post=g[g.tp == tp].score.iloc[0],
                         late=(tp == late_tp[tr])))
R = pd.DataFrame(recs); R['ratio'] = R.post / R.base
out = {}
def lodo(target):
    T = R[target].copy(); rows = []
    for i, r in T.iterrows():
        tr_same = T[(T.trial == r.trial) & (T.donor != r.donor)]
        tr_all = T[T.donor != r.donor]
        p = dict(B0=r.base, B1=tr_same.post.mean(), M1=r.base * tr_same.ratio.median(),
                 M2=r.base * tr_all.ratio.median())
        rows.append({k: abs(v - r.post) / r.base for k, v in p.items()} | dict(trial=r.trial, donor=r.donor))
    return pd.DataFrame(rows)
L = lodo(R.late)
E = lodo((R.trial == 'CD3xCD19') & (R.tp == 'D28'))
mL = L[['B0', 'B1', 'M1', 'M2']].mean(); mE = E[['B0', 'B1', 'M1']].mean()
rho = {tr: spearmanr(R[R.late & (R.trial == tr)].base, R[R.late & (R.trial == tr)].post)[0] for tr in late_tp}
out['late_mean_rel_error'] = mL.round(4).to_dict(); out['early_D28_mean_rel_error'] = mE.round(4).to_dict()
out['late_median_rel_error'] = L[['B0', 'B1', 'M1', 'M2']].median().round(4).to_dict()
out['spearman_base_late'] = {k: round(v, 3) for k, v in rho.items()}
out['median_f'] = {'CD3xCD19_D28': round(R[(R.trial == 'CD3xCD19') & (R.tp == 'D28')].ratio.median(), 3),
                   'CD3xCD19_M3': round(R[(R.trial == 'CD3xCD19') & R.late].ratio.median(), 3),
                   'CART_post': round(R[(R.trial == 'CART') & R.late].ratio.median(), 3)}
out['H1'] = bool(mL.M1 <= 0.9 * mL.B1 and all(v >= 0.5 for v in rho.values()))
out['H1_error_part'] = bool(mL.M1 <= 0.9 * mL.B1)
out['H2'] = bool(mL.M1 <= 0.9 * mL.B0)
out['H3'] = bool(mL.M2 <= 1.10 * mL.M1)
out['per_donor_late'] = L.round(3).to_dict('records')
json.dump(out, open('result.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'per_donor_late'}, indent=1))
print(L.round(3).to_string())
