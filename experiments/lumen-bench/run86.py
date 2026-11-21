import pandas as pd,numpy as np
R=pd.read_csv('resp_table.csv')
a,b,c,d=R.r_drug+.5,R.n_drug-R.r_drug+.5,R.r_ctrl+.5,R.n_ctrl-R.r_ctrl+.5
R['lor']=np.log(a*d/(b*c));R['se']=np.sqrt(1/a+1/b+1/c+1/d)
P=R[R.week<=104].sort_values('week').groupby(['trial','arm']).tail(1).reset_index(drop=True)
print(P[['trial','arm','molecule','cls','week','lor','se']].round(2).to_string())
def iv(x,se): w=1/se**2; return float((w*x).sum()/w.sum())
def dl(x,se):
    w=1/se**2; m=(w*x).sum()/w.sum(); Q=(w*(x-m)**2).sum(); k=len(x)
    return max(0.0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum()))
rows=[]
for tr in P.trial.unique():
    T,H=P[P.trial!=tr],P[P.trial==tr]
    mu=iv(T.lor,T.se); t2=dl(T.lor.values,T.se.values)
    for _,r in H.iterrows():
        cc=T[T.cls==r.cls]; d2=iv(cc.lor,cc.se) if len(cc) else mu
        anc=T[T.molecule==r.molecule]; lvl='molecule'
        if not len(anc): anc=T[T.target==r.target]; lvl='target'
        if len(anc):
            am=iv(anc.lor,anc.se); ase=1/np.sqrt((1/anc.se**2).sum()); w=t2/(t2+ase**2)
            d3=mu+w*(am-mu)
        else: d3=mu; lvl='none'; w=np.nan
        rows.append(dict(trial=tr,arm=r.arm,anchor=lvl,obs=r.lor,se=r.se,D0=0,D1=mu,D2=d2,D3=d3,w=w,tau2=t2))
X=pd.DataFrame(rows)
def score(df,m):
    e=df.obs-df[m]; w=1/df.se**2
    return float((w*e**2).sum()/w.sum()), float((e**2-df.se**2).mean())
for name,sub in [('PRIMARY: arms with anchor',X[X.anchor!='none']),('secondary: all arms',X)]:
    print('\n'+name,'n=',len(sub))
    for m in ['D0','D1','D2','D3']:
        s,nc=score(sub,m); print(f'  {m}: IV-MSE={s:.3f}  noise-corrected MSE={nc:.3f}')
print('\n',X.round(2).to_string())
# Task C'
C=[]
for (tr,arm),g in R.groupby(['trial','arm']):
    e=g[g.week<=26];l=g[g.week>=48]
    if len(e) and len(l):
        e=e.sort_values('week').iloc[-1]
        for _,x in l.iterrows(): C.append(dict(trial=tr,arm=arm,te=e.week,tl=x.week,de=e.lor,see=e.se,dl=x.lor,sel=x.se))
C=pd.DataFrame(C)
pool=iv(P.lor,P.se); t2=dl(P.lor.values,P.se.values)
C['C0']=C.de;C['C1']=0.0;C['C2']=pool+(t2/(t2+C.see**2))*(C.de-pool)
print("\nTASK C' rows",len(C),' tau2=',round(t2,3))
for m in ['C0','C1','C2']:
    e=C.dl-C[m];w=1/C.sel**2;print(f'  {m}: IV-MSE={(w*e**2).sum()/w.sum():.3f}  noise-corr={(e**2-C.sel**2).mean():.3f}')
print(C.round(2).to_string())
X.to_csv('r86_taskD.csv',index=False);C.to_csv('r86_taskC.csv',index=False)
