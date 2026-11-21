import pandas as pd,numpy as np
from scipy.optimize import least_squares
D=pd.read_csv('fr_table.csv');D['y']=np.log(D.fr)
TAU={'fast':1.0,'slow':6.9}
g=lambda t,c:1-np.exp(-np.asarray(t,float)/TAU[c])
def wrmse(rows):  # rows: list of (trial, err)
    df=pd.DataFrame(rows,columns=['trial','e'])
    per=df.groupby('trial').e.apply(lambda e:np.mean(e**2))
    return float(np.sqrt(per.mean())),per.apply(np.sqrt)
# ---- Task A
C=D[D.role=='ctrl'];trials=C.trial.unique();A={'A0':[],'A1':[],'A2':[]}
for tr in trials:
    tr_,te=C[C.trial!=tr],C[C.trial==tr]
    A['A0']+= [(tr,v-tr_.y.mean()) for v in te.y]
    k=-(tr_.y*tr_.week).sum()/(tr_.week**2).sum()
    A['A1']+= [(tr,v+k*t) for v,t in zip(te.y,te.week)]
    f=lambda p,t:np.log(p[0]+(1-p[0])*np.exp(-t/p[1]))
    p=least_squares(lambda p:f(p,tr_.week.values)-tr_.y.values,[0.4,10],bounds=([0.01,0.5],[0.99,500])).x
    A['A2']+= [(tr,v-f(p,t)) for v,t in zip(te.y,te.week)]
print('TASK A  control arm, LOTO RMSE (log FR)')
resA={m:wrmse(r) for m,r in A.items()}
for m,(s,per) in resA.items(): print(f'  {m}: {s:.3f}')
print(pd.DataFrame({m:per for m,(s,per) in resA.items()}).round(3).to_string())
# ---- build effects
E=[]
for (tr,arm),g_ in D[D.role=='drug'].groupby(['trial','arm']):
    for _,r in g_.iterrows():
        c=D[(D.trial==tr)&(D.role=='ctrl')&(D.week==r.week)]
        if len(c): E.append(dict(trial=tr,arm=arm,cls=r.cls,week=r.week,d=r.y-c.y.mean()))
E=pd.DataFrame(E);E.to_csv('effects.csv',index=False)
print('\nEFFECTS d=log(FR_drug/FR_ctrl)\n',E.round(3).to_string())
# ---- Task B
B={'B0':[],'B1':[],'M1':[],'M2':[]}
for tr in E.trial.unique():
    tr_,te=E[E.trial!=tr],E[E.trial==tr]
    for _,r in te.iterrows():
        B['B0'].append((tr,r.d));B['B1'].append((tr,r.d-tr_.d.mean()))
        cc=tr_[tr_.cls==r.cls]
        B['M1'].append((tr,r.d-(cc.d.mean() if len(cc) else tr_.d.mean())))
        if len(cc):
            gg=g(cc.week,r.cls);Ec=(gg*cc.d).sum()/(gg**2).sum()
        else: Ec=0
        B['M2'].append((tr,r.d-Ec*g(r.week,r.cls)))
print('\nTASK B  add-on effect, LOTO RMSE')
resB={m:wrmse(r) for m,r in B.items()}
for m,(s,per) in resB.items(): print(f'  {m}: {s:.3f}')
print(pd.DataFrame({m:per for m,(s,per) in resB.items()}).round(3).to_string())
# ---- Task C
Crow=[]
for (tr,arm),a in E.groupby(['trial','arm']):
    early=a[a.week<=26];late=a[a.week>=36]
    if not len(early) or not len(late): continue
    e=early.sort_values('week').iloc[-1]
    for _,l in late.iterrows():
        c=e.cls
        Crow.append(dict(trial=tr,arm=arm,te=e.week,tl=l.week,de=e.d,dl=l.d,
          C0=l.d-e.d,C1=l.d,C2=l.d-e.d*g(l.week,c)/g(e.week,c)))
Cd=pd.DataFrame(Crow)
print('\nTASK C  fast-forward early->late, rows:',len(Cd))
for m in ['C0','C1','C2']:
    s,per=wrmse(list(zip(Cd.trial,Cd[m])));print(f'  {m}: {s:.3f}')
Cd.round(3).to_csv('taskC.csv',index=False)
print(Cd.round(3)[['trial','arm','te','tl','de','dl']].to_string())
