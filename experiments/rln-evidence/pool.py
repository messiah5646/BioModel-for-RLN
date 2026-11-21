import pandas as pd,numpy as np
T=pd.read_csv('rln_table.csv')
def pool(k,n):
    k=np.asarray(k,float);n=np.asarray(n,float)
    y=np.log((k+.5)/(n-k+.5));v=1/(k+.5)+1/(n-k+.5);w=1/v
    m=(w*y).sum()/w.sum();Q=(w*(y-m)**2).sum();df=len(y)-1
    t2=max(0,(Q-df)/(w.sum()-(w**2).sum()/w.sum())) if df>0 else 0
    ws=1/(v+t2);mu=(ws*y).sum()/ws.sum();se=1/np.sqrt(ws.sum())
    f=lambda x:1/(1+np.exp(-x))
    I2=max(0,(Q-df)/Q) if Q>0 and df>0 else 0
    return f(mu),f(mu-1.96*se),f(mu+1.96*se),len(y),int(n.sum()),I2
rows=[]
for g,d in T.groupby('group'):
    for ep in ['cr','orr']:
        dd=d.dropna(subset=[ep])
        if len(dd): p,lo,hi,k,N,I2=pool(dd[ep],dd.n);rows.append((g,ep.upper(),k,N,round(p,2),round(lo,2),round(hi,2),round(I2,2),round(float(dd.month.mean()),1)))
R=pd.DataFrame(rows,columns=['group','endpoint','studies','patients','pooled','lo95','hi95','I2','mean_month'])
print(R.to_string(index=False))
R.to_csv('pooled.csv',index=False)
def comp(a,b,ep):
    x=R[(R.group==a)&(R.endpoint==ep)].iloc[0];y=R[(R.group==b)&(R.endpoint==ep)].iloc[0]
    sep=x.lo95>y.hi95 or y.lo95>x.hi95
    return f'{a} {x.pooled} [{x.lo95},{x.hi95}] vs {b} {y.pooled} [{y.lo95},{y.hi95}] -> {"SUPPORTED difference" if sep else "not distinguishable"}'
print('\nQ2 CR :',comp('anti_CD20','conventional_switch','CR'))
print('Q2 ORR:',comp('anti_CD20','conventional_switch','ORR'))
print('Q1: CAR_T not poolable under the rule; plasma_cell has 1 study ->',comp('plasma_cell','conventional_switch','CR'))
print('extra (not pre-registered): BAFF_axis vs conventional CR:',comp('BAFF_axis','conventional_switch','CR'))
