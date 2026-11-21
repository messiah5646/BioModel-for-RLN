"""Prior-only LUMEN forecast from a (lor, se) arm table: identical maths to lumen_forecast.forecast() with early=None."""
import sys,json,math,numpy as np,pandas as pd
sys.path.insert(0,'../..');from lumen_forecast import iv,dl
def prior(molecule,target,R,holdout=None,maxweek=104):
    P=R[R.week<=maxweek].sort_values('week').groupby(['trial','arm']).tail(1)
    if holdout: P=P[P.trial!=holdout]
    mu,mse=iv(P.lor,P.se);t2=dl(P.lor.values,P.se.values)
    anc=P[P.molecule==molecule];lvl='molecule'
    if not len(anc): anc=P[P.target==target];lvl='target'
    if len(anc):
        am,ase=iv(anc.lor,anc.se);w=t2/(t2+ase**2);m=mu+w*(am-mu);v=t2*(1-w)+t2
    else: lvl='none';m,v=mu,t2+mse**2
    sd=math.sqrt(v)
    return dict(anchor=lvl,logOR=round(m,3),sd=round(sd,3),OR=round(math.exp(m),2),
      OR_95=[round(math.exp(m-1.96*sd),2),round(math.exp(m+1.96*sd),2)],P_benefit=round(0.5*(1+math.erf(m/(sd*math.sqrt(2)))),2),
      pooled_logOR=round(mu,3))
