"""Score the rebuilt resp_table.csv against result86.txt (PREREG-Q17 P1-P3)."""
import pandas as pd,numpy as np,re,subprocess,sys
R=pd.read_csv('resp_table.csv')
a,b,c,d=R.r_drug+.5,R.n_drug-R.r_drug+.5,R.r_ctrl+.5,R.n_ctrl-R.r_ctrl+.5
R['lor']=np.log(a*d/(b*c));R['se']=np.sqrt(1/a+1/b+1/c+1/d)
txt=open('result86.txt').read().splitlines()
prim=[l.split() for l in txt[1:30]]
O=pd.DataFrame([[x[1],x[2],float(x[5]),float(x[6]),float(x[7])] for x in prim],columns=['trial','arm','week','lor','se'])
P=R[R.week<=104].sort_values('week').groupby(['trial','arm']).tail(1)
M=O.merge(P[['trial','arm','week','lor','se']],on=['trial','arm'],how='left',suffixes=('_86','_new'))
M['ok']=(M.week_86==M.week_new)&((M.lor_86-M.lor_new).abs()<=0.01)&((M.se_86-M.se_new).abs()<=0.01)
print(M.round(3).to_string());P1=bool(M.ok.all()) and len(P)==29;print('P1',P1,'(rebuilt arms',len(P),')')
i=[k for k,l in enumerate(txt) if l.startswith("TASK C'")][0]
C86=pd.DataFrame([l.split()[1:] for l in txt[i+5:i+18]],columns=['trial','arm','te','tl','de','see','dl','sel','C0','C1','C2']).iloc[:,:8]
for k in ['te','tl','de','see','dl','sel']: C86[k]=C86[k].astype(float)
rows=[]
for (tr,arm),g in R.groupby(['trial','arm']):
    e=g[g.week<=26];l=g[g.week>=48]
    if len(e) and len(l):
        e=e.sort_values('week').iloc[-1]
        for _,x in l.iterrows(): rows.append(dict(trial=tr,arm=arm,te=e.week,tl=x.week,de=e.lor,see=e.se,dl=x.lor,sel=x.se))
C=pd.DataFrame(rows);MC=C86.merge(C,on=['trial','arm','tl'],how='outer',suffixes=('_86','_new'))
MC['ok']=True
for k in ['te','de','see','dl','sel']: MC['ok']&=(MC[k+'_86']-MC[k+'_new']).abs()<=0.01
print(MC.round(2).to_string());P2=bool(MC.ok.all()) and len(C)==13;print('P2',P2,'(rebuilt rows',len(C),')')
out=subprocess.run([sys.executable,'-I','run86.py'],capture_output=True,text=True).stdout
open('resp_table/result86_rerun.txt','w').write(out)
num=lambda t:[float(x) for x in re.findall(r'(?:IV-MSE=)(-?[\d.]+)',t)]
a86,anew=num('\n'.join(txt)),num(out)
print('IV-MSE 86 ',a86);print('IV-MSE new',anew)
P3=len(a86)==len(anew)==11 and max(abs(x-y) for x,y in zip(a86,anew))<=0.002;print('P3',P3)
print('PASS' if P1 and P2 and P3 else 'FAIL')
