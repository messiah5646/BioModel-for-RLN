"""Audit r120, exploratory (not pre-registered): sensitivity of round-119 verdict to error metric and to the low-baseline patients."""
import sys, pandas as pd, numpy as np
S = pd.read_csv(sys.argv[1]).pivot(index='patient', columns='tp', values='mono')
def ev(W, rel=True):
    r = []
    for p in W.index:
        b, o = W.loc[p, 'pre'], W.loc[p, 'post']; oth = W.drop(p)
        pr = dict(B0=b, B1=oth.post.mean(), M2=b*0.73)
        r.append({k: abs(v-o)/(b if rel else 1) for k, v in pr.items()})
    return pd.DataFrame(r, index=W.index).mean().round(3).to_dict()
print('prereg relative, all 5 :', ev(S))
print('absolute error, all 5  :', ev(S, False))
print('relative, drop P1      :', ev(S.drop('P1')))
print('relative, drop P1,P2   :', ev(S.drop(['P1','P2'])))
print('log-ratio error |log(pred/obs)|:')
r=[]
for p in S.index:
    b,o=S.loc[p,'pre'],S.loc[p,'post']; oth=S.drop(p)
    r.append(dict(B0=abs(np.log(b/o)),B1=abs(np.log(oth.post.mean()/o)),M2=abs(np.log(0.73*b/o))))
print(' ', pd.DataFrame(r).mean().round(3).to_dict())
print('CV of post scores', round(S.post.std()/S.post.mean(),2), 'CV of pre', round(S.pre.std()/S.pre.mean(),2))
