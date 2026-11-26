"""LUMEN forecaster: predict a lupus nephritis add-on drug's complete-renal-response effect at ~1 year.

Evidence-anchored, as validated in rounds 85-87:
  prior  = same molecule > same target > all trials (inverse-variance, shrunk with between-trial tau^2)
  update = the drug's own early read-out (week <= 26), allowing early->late drift (estimated from trial data)
Mechanism class is NOT used for magnitude (leave-one-molecule-out: it does worse than the pooled average).
Round 116 (2026-11-23): across 40 non-renal SLE trials, molecule/target anchoring was 31% WORSE than the pooled mean
(Supported), and in LN it never cleared its pre-registered margin. The command-line default is therefore the pooled prior
(--anchor restores the old behaviour).
Round 127 (2026-11-23, RT-1 time-forward test, pre-registered c8e07af): the shrunk early update FAILED (lost to raw
carry-forward by >=10%); the CLI headline is now the pooled prior even when --early is given (--use-early-update restores it). The function keeps anchor=True as its default so that frozen forecasts and audits
(check_frozen.py, q53/run_q53.py) reproduce exactly.
Round 129 (2026-11-23, RT-2, PREREG d5b5f84): a control-arm meta-regression (posting year + steroid-stringent definition)
did NOT beat the pooled control rate (60 trials: INCONCLUSIVE 1.06x; LN 20 trials FAIL 1.17x; time-forward FAIL 1.38x).
LUMEN therefore outputs relative effects (log OR) only; control rate = pooled, not modelled. Pooled control-rate 95%
intervals covered 0.88 (LN 0.85) < 0.90, so any future absolute-rate output must widen them.
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
def forecast(molecule,target,early=None,R=None,holdout=None,anchor=True):
    R=load() if R is None else R
    P=R[R.week<=104].sort_values('week').groupby(['trial','arm']).tail(1)
    if holdout: P=P[P.trial!=holdout]   # for validation: remove the trial being predicted
    mu,_=iv(P.lor,P.se);t2=dl(P.lor.values,P.se.values)
    anc=P[P.molecule==molecule] if anchor else P.iloc[0:0];lvl='molecule'
    if not len(anc) and anchor: anc=P[P.target==target];lvl='target'
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
    ap.add_argument('--holdout');ap.add_argument('--anchor',action='store_true',help='use molecule/target anchoring (not default since round 116)');ap.add_argument('--early',nargs=5,type=float,metavar=('WEEK','R_DRUG','N_DRUG','R_CTRL','N_CTRL'))
    ap.add_argument('--use-early-update',action='store_true',help='make the shrunk early update the headline (off by default since round 127)')
    x=ap.parse_args()
    out=forecast(x.molecule,x.target,None,holdout=x.holdout,anchor=x.anchor)
    if x.early:
        upd=forecast(x.molecule,x.target,x.early,holdout=x.holdout,anchor=x.anchor)
        if x.use_early_update: out=upd
        else: out['early_update_EMERGING']={k:upd[k] for k in ('early_logOR','post_logOR','post_sd','OR','OR_95') if k in upd}
        out['early_note']=('Round 127 (RT-1, time-forward, PREREG c8e07af): the shrunk early update FAILED, losing to raw carry-forward '
                           'of the early read-out by >=10% (13 arms / 7 trials, voclosporin excluded; Emerging). Headline = pooled prior, BUT on the same arms '
                           'the pooled headline was itself worse than the update (IV-MSE 0.311 vs 0.259) and than carry-forward (0.146); '
                           'in LN without voclosporin (4 trials) null was best (0.279 vs pooled 0.626). NO early rule is validated '
                           '(reviewer R-20261124-1, round 132). '
                           'Round 141 (RT-4, regulatory documents, PREREG-RT4, UNDERPOWERED: 10 trials / 17 arms, only 1 RT-4-new trial scored): '
                           'pooled again lost, IV-MSE pooled 0.322 vs carry-forward 0.110 vs null 0.253; LN only (5 trials) pooled 0.497 vs null 0.299; '
                           'DEFAULT-FAIL pattern in all 9 leave-one-molecule-out variants. Rule 8: early-read-out rule settled; no further re-tests '
                           '(round 144, reviewer R-20261126-3).')
        # Rule 2b display (R-20261126-3, round 144): pooled loses to a trivial baseline by >=10% on every early subset,
        # so with --early no single headline is given; M0/M1/C0 are printed side by side. forecast() is unchanged.
        if not x.use_early_update:
            wk,rd,nd,rc,nc=x.early
            a,b,c,d=rd+.5,nd-rd+.5,rc+.5,nc-rc+.5
            le=float(np.log(a*d/(b*c)));se=float(np.sqrt(1/a+1/b+1/c+1/d))
            out['headline']='NONE: no validated rule; pooled lost to carry-forward on the all-trial early subsets (RT-1, RT-4) and to null in LN'
            out['early_side_by_side']={
                'M0_null':{'OR':1.0},
                'M1_pooled_prior':{'OR':out['OR'],'OR_95':out['OR_95']},
                'C0_carry_forward':{'OR':round(float(np.exp(le)),2),'OR_95_sampling_only':[round(float(np.exp(le-1.96*se)),2),round(float(np.exp(le+1.96*se)),2)],
                                    'note':'raw early read-out (week %g) carried to the late visit; early->late drift not added'%wk}}
    out['control_rate']='pooled, not modelled; absolute arm rates not output (round 129, RT-2 PREREG d5b5f84)'
    print(json.dumps(out,indent=1))
