"""EXPLORATORY (not pre-registered): sensitivity of round-96 verdicts to the error scale and to the P2 outlier."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
d = pd.read_csv('scores.csv'); bt = {'CD3xCD19':'Base','CART':'pre'}; lt = {'CD3xCD19':'M3','CART':'post'}
R = []
for (tr, dn), g in d.groupby(['trial','donor']):
    R.append(dict(trial=tr, donor=dn, base=g[g.tp==bt[tr]].score.iloc[0], post=g[g.tp==lt[tr]].score.iloc[0]))
R = pd.DataFrame(R); R['ratio'] = R.post/R.base
def run(T, scale):
    e = {k: [] for k in ['B0','B1','M1','M2']}
    for i, r in T.iterrows():
        s = T[(T.trial==r.trial)&(T.donor!=r.donor)]; a = T[T.donor!=r.donor]
        p = dict(B0=r.base, B1=s.post.mean(), M1=r.base*s.ratio.median(), M2=r.base*a.ratio.median())
        sc = r.base if scale=='base' else T[T.trial==r.trial].base.std()
        for k in e: e[k].append(abs(p[k]-r.post)/sc)
    return {k: round(np.mean(v),3) for k,v in e.items()}
print('abs error / trial SD of baseline, all 18:', run(R,'sd'))
print('rel error, drop P2 (baseline 0.06):', run(R[R.donor!='P2'],'base'))
print('abs/SD, drop P2:', run(R[R.donor!='P2'],'sd'))
print('Spearman CD3xCD19 drop P2:', round(spearmanr(R[(R.trial=="CD3xCD19")&(R.donor!="P2")].base, R[(R.trial=="CD3xCD19")&(R.donor!="P2")].post)[0],3))
print('pooled Spearman base vs post (trial-ranked):', round(spearmanr(R.groupby("trial").base.rank(pct=True), R.groupby("trial").post.rank(pct=True))[0],3))
