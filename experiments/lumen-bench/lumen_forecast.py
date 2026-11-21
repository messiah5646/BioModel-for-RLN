"""LUMEN forecaster: predict a lupus nephritis add-on drug's complete-renal-response effect at ~1 year.

Evidence-anchored, as validated in rounds 85-87:
  prior  = same molecule > same target > all trials (inverse-variance, shrunk with between-trial tau^2)
  update = the drug's own early read-out (week <= 26), allowing early->late drift (estimated from trial data)
Mechanism class is NOT used for magnitude (leave-one-molecule-out: it does worse than the pooled average).
Research and simulation only. Not clinical guidance.
"""
import sys,json,argparse
import numpy as np,pandas as pd
DELAYED={'CD20','BAFF','CD40','CD80/86'}  # heuristic: NOBILITY/BLISS-LN separate late; not a validated rule
def load(path='resp_table.csv'):
    R=pd.read_csv(path)
    a,b,c,d=R.r_drug+.5,R.n_drug-R.r_drug+.5,R.r_ctrl+.5,R.n_ctrl-R.r_ctrl+.5
    R['lor']=np.log(a*d/(b*c));R['se']=np.sqrt(1/a+1/b+1/c+1/d);return R
def iv(x,se): w=1/se**2;return float((w*x).sum()/w.sum()),float(1/np.sqrt(w.sum()))
def dl(x,se):
    w=1/se**2;m=(w*x).sum()/w.sum();Q=(w*(x-m)**2).sum()
    return max(0.0,(Q-(len(x)-1))/(w.sum()-(w**2).sum()/w.sum()))
def drift_var(R):
    rows=[]
    for _,g in R.groupby(['trial','arm']):
        e,l=g[g.week<=26],g[g.week>=48]
        if len(e) and len(l):
            e=e.sort_values('week').iloc[-1];l=l.sort_values('week').iloc[0]
            rows.append((l.lor-e.lor,e.se**2+l.se**2))
    x=np.array(rows);return max(0.0,float(np.mean(x[:,0]**2-x[:,1])))
def forecast(molecule,target,early=None,R=None,holdout=None):
    R=load() if R is None else R
    P=R[R.week<=104].sort_values('week').groupby(['trial','arm']).tail(1)
    if holdout: P=P[P.trial!=holdout]   # for validation: remove the trial being predicted
    mu,_=iv(P.lor,P.se);t2=dl(P.lor.values,P.se.values)
    anc=P[P.molecule==molecule];lvl='molecule'
    if not len(anc): anc=P[P.target==target];lvl='target'
    if len(anc):
        am,ase=iv(anc.lor,anc.se);w=t2/(t2+ase**2)
        m0=mu+w*(am-mu);v0=t2*(1-w)+t2   # uncertainty about the drug's mean + a new trial's deviation
    else:
        lvl='none (pooled; mechanism class not used)';m0,v0=mu,t2+_**2
    out=dict(anchor=lvl,prior_logOR=round(m0,3),prior_sd=round(float(np.sqrt(v0)),3))
    if early:
        wk,rd,nd,rc,nc=early
        a,b,c,d=rd+.5,nd-rd+.5,rc+.5,nc-rc+.5
        le=np.log(a*d/(b*c));ve=1/a+1/b+1/c+1/d+drift_var(R)
        m=(m0/v0+le/ve)/(1/v0+1/ve);v=1/(1/v0+1/ve)
        out.update(early_logOR=round(float(le),3),post_logOR=round(float(m),3),post_sd=round(float(np.sqrt(v)),3))
        if target in DELAYED: out['warning']='target class separated late in NOBILITY/BLISS-LN; a week-24 read-out may understate the 1-year effect'
    else: m,v=m0,v0
    lo,hi=m-1.96*np.sqrt(v),m+1.96*np.sqrt(v)
    out.update(OR=round(float(np.exp(m)),2),OR_95=[round(float(np.exp(lo)),2),round(float(np.exp(hi)),2)],P_benefit=round(float(1-0.5*(1+__import__('math').erf((0-m)/np.sqrt(2*v)))),2))
    return out
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('molecule');ap.add_argument('target')
    ap.add_argument('--holdout');ap.add_argument('--early',nargs=5,type=float,metavar=('WEEK','R_DRUG','N_DRUG','R_CTRL','N_CTRL'))
    x=ap.parse_args();print(json.dumps(forecast(x.molecule,x.target,x.early,holdout=x.holdout),indent=1))
