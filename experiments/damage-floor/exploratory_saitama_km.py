"""EXPLORATORY, not scored (Q16 PREREG): composite risk groups (CI>3 is one factor), not CI groups.
Digitised from Fig 2, Arthritis Res Ther 2026, doi 10.1186/s13075-026-03876-w (PMC13523271). Reading error ~+/-0.02."""
import numpy as np
from scipy.optimize import least_squares
t=np.array([3,6,12,24,36.])
C={'low':np.array([.66,.84,.90,.93,.93]),'moderate':np.array([.44,.56,.66,.69,.69]),'high':np.array([.18,.31,.37,.41,.41])}
def fit(groups,mode,tmask):
    # params: mode 'pool' [c,k]; 'floor' [k,c_g...]; 'rate' [c,k_g...]
    G=len(groups)
    def pred(p,g,tt):
        if mode=='pool': c,k=p
        elif mode=='floor': k=p[0]; c=p[1+g]
        else: c=p[0]; k=p[1+g]
        return c*(1-np.exp(-k*tt))
    def res(p): return np.concatenate([pred(p,g,t[tmask])-C[n][tmask] for g,n in enumerate(groups)])
    p0={'pool':[.6,.3],'floor':[.3]+[.6]*G,'rate':[.6]+[.3]*G}[mode]
    lo={'pool':[0,.01],'floor':[.01]+[0]*G,'rate':[0]+[.01]*G}[mode]
    hi={'pool':[1,5],'floor':[5]+[1]*G,'rate':[1]+[5]*G}[mode]
    r=least_squares(res,p0,bounds=(lo,hi))
    n=tmask.sum()*G; k=len(p0); sse=max((r.fun**2).sum(),1e-12)
    aicc=n*np.log(sse/n)+2*k+(2*k*(k+1)/(n-k-1) if n-k-1>0 else np.inf)
    return r.x,sse,aicc,lambda g,tt: pred(r.x,g,tt)
for groups in (['low','high'],['low','moderate','high']):
    print('== groups',groups)
    m=np.ones(5,bool)
    for mode in ('pool','floor','rate'):
        x,sse,aicc,_=fit(groups,mode,m); print(f'  {mode:5s} params={np.round(x,3)} SSE={sse:.4f} AICc={aicc:.1f}')
    # held-out: fit t<=6, predict t=36
    m6=t<=6
    out={}
    for mode in ('pool','floor','rate'):
        x,sse,aicc,pr=fit(groups,mode,m6)
        out[mode]=np.mean([abs(pr(g,36.)-C[n][-1]) for g,n in enumerate(groups)])
    out['carry']=np.mean([abs(C[n][1]-C[n][-1]) for n in groups])
    print('  held-out MAE at 36 mo (fit on 3,6 mo):',{k:round(v,3) for k,v in out.items()})
