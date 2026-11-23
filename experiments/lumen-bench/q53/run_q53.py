"""Q53 / round 112 analysis per PREREG-Q53.md. Run from experiments/lumen-bench: python3 q53/run_q53.py"""
import sys,os
sys.path.insert(0,os.getcwd())
import numpy as np,pandas as pd
import lumen_forecast as LF

def iv(x,se): w=1/se**2; return float((w*x).sum()/w.sum())
def dl(x,se):
    w=1/se**2; m=(w*x).sum()/w.sum(); Q=(w*(x-m)**2).sum(); k=len(x)
    return max(0.0,(Q-(k-1))/(w.sum()-(w**2).sum()/w.sum()))

def rows(R):
    P=R[R.week<=104].sort_values('week').groupby(['trial','arm']).tail(1)
    pool=iv(P.lor,P.se); t2=dl(P.lor.values,P.se.values)
    C=[]
    for (tr,arm),g in R.groupby(['trial','arm']):
        e=g[g.week<=26]; l=g[(g.week>=48)&(g.week<=104)]
        if not (len(e) and len(l)): continue
        e=e.sort_values('week').iloc[-1]
        Rh=R[R.trial!=tr]
        f=LF.forecast(e.molecule,e.target,early=(e.week,e.r_drug,e.n_drug,e.r_ctrl,e.n_ctrl),R=Rh,holdout=tr)
        for _,x in l.iterrows():
            C.append(dict(trial=tr,arm=arm,n_drug=e.n_drug,te=e.week,tl=x.week,de=e.lor,see=e.se,dl=x.lor,sel=x.se,
                          C0=e.lor,C1=0.0,C2=pool+(t2/(t2+e.se**2))*(e.lor-pool),F=f['post_logOR']))
    return pd.DataFrame(C),pool,t2

def score(C,m,w=None):
    e=C.dl-C[m]; w=1/C.sel**2 if w is None else w
    return float((w*e**2).sum()/w.sum()), float((e**2-C.sel**2).mean())

def block(C,label,w=None):
    out={m:score(C,m,w) for m in ['C0','C1','C2','F']}
    print(f'\n{label}: rows={len(C)} arms={C.groupby(["trial","arm"]).ngroups} trials={C.trial.nunique()}')
    for m,(s,nc) in out.items(): print(f'  {m}: IV-MSE={s:.3f}  noise-corr={nc:.3f}')
    r=out['C0'][0]/out['C1'][0]; print(f'  C0/C1={r:.2f} (H1 pass if <=0.90): {"PASS" if r<=0.9 else "FAIL"}')
    return out

res={}
for name in ['resp_table.csv','resp_table_plus_pub.csv']:
    R=LF.load(name); C,pool,t2=rows(R)
    print('\n=====',name,f'pool={pool:.3f} tau2={t2:.3f} drift_var(all)={LF.drift_var(R):.3f}')
    res[name]=block(C,'ALL ROWS (primary)')
    if name=='resp_table_plus_pub.csv':
        Cx=C
        a=C.sort_values('tl').groupby(['trial','arm']).tail(1); ra=block(a,'(a) one row per arm, latest late visit')
        b=C[C.n_drug>=10]; rb=block(b,'(b) n_drug>=10')
        w=1/C.sel**2; w=w/w.groupby(C.trial).transform('sum'); rc=block(C,'(c) trial-balanced weights',w)
        o=res[name]; h1=o['C0'][0]<=0.9*o['C1'][0]
        rob=all(x['C0'][0]<=0.9*x['C1'][0] for x in [ra,rb,rc])
        h2=o['F'][0]<=o['C0'][0] and o['F'][0]<=0.9*o['C1'][0]
        print(f'\nH1 primary: {"PASS" if h1 else "FAIL"}; H1 robust (a,b,c): {"PASS" if rob else "FAIL"}; H2 forecaster: {"PASS" if h2 else "FAIL"}')
        # H3 drift by horizon
        def drift(D): return float(np.mean((D.dl-D.de)**2-D.see**2-D.sel**2)),len(D)
        s=drift(C[C.tl<=52]); l=drift(C[C.tl>=76])
        print(f'H3 drift var: 48-52 wk {s[0]:.3f} (n rows {s[1]}); 76-104 wk {l[0]:.3f} (n rows {l[1]}); ratio {l[0]/s[0] if s[0]>0 else float("nan"):.2f}')
        print(C.round(2).to_string())
        # leave-one-row diagnostics: which row drives C0 vs C1
        print('\nper-row IV contribution share (C0 err^2 * w):')
        w=1/C.sel**2; c0=(w*(C.dl-C.C0)**2); c1=(w*(C.dl-C.C1)**2)
        print(pd.DataFrame(dict(trial=C.trial,arm=C.arm,tl=C.tl,C0=(c0/w.sum()).round(3),C1=(c1/w.sum()).round(3))).to_string())
        for tr in C.trial.unique():
            D=C[C.trial!=tr]; print(f'  drop {tr:10s}: C0/C1={score(D,"C0")[0]/score(D,"C1")[0]:.2f}  F/C0={score(D,"F")[0]/score(D,"C0")[0]:.2f}')

# EXPLORATORY (not pre-registered): drop the molecule that carries the signal (voclosporin: AURA-LV + AURORA)
print('\nEXPLORATORY, not pre-registered:')
D=Cx[~Cx.trial.isin(['AURA-LV','AURORA'])]
block(D,'drop voclosporin trials (IV weights)')
w=1/D.sel**2; w=w/w.groupby(D.trial).transform('sum'); block(D,'drop voclosporin trials (trial-balanced)',w)
