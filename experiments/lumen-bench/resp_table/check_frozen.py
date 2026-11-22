"""P4 + downstream: reproduce arm_table_reconstructed.csv, r90_frozen.json and w1_frozen_2026-11-21.json from the rebuilt tables."""
import sys,json,numpy as np,pandas as pd
sys.path.insert(0,'.');sys.path.insert(0,'prospective/w1')
import lumen_forecast as L
from forecast_prior import prior
R=L.load('resp_table.csv');RP=L.load('resp_table_plus_pub.csv')
A=pd.read_csv('prospective/w1/arm_table_reconstructed.csv')
RPl=RP.sort_values('week').groupby(['trial','arm','week']).tail(1)
M=A.merge(RPl[['trial','arm','week','lor','se']],on=['trial','arm','week'],suffixes=('_94','_new'))
d=max((M.lor_94-M.lor_new).abs().max(),(M.se_94-M.se_new).abs().max())
print('P4 arm table: matched',len(M),'/',len(A),'rows; max |diff|',round(d,4))
F=json.load(open('prospective/r90_frozen.json'))
ins={'NOBILITY_obinutuzumab_wk104':('obinutuzumab','CD20',(52,22,63,14,62),'NOBILITY'),
     'TULIP-LN1_anifro_int_wk104':('anifrolumab','IFNAR',(52,23,51,15,49),'TULIP-LN1'),
     'TULIP-LN1_anifro_basic_wk104':('anifrolumab','IFNAR',(52,7,45,15,49),'TULIP-LN1')}
for k,(m,t,e,h) in ins.items():
    f=L.forecast(m,t,e,R=R,holdout=h);print('r90',k,'frozen',F[k]['post_logOR'],F[k]['post_sd'],'| rebuilt',f['post_logOR'],f['post_sd'])
W=json.load(open('prospective/w1/w1_frozen_2026-11-21.json'));mx=0
for n,v in W.items():
    f=prior(v['molecule'],v['target'],RP);mx=max(mx,abs(f['logOR']-v['logOR']),abs(f['sd']-v['sd']))
print('W1 13 frozen forecasts: max |diff| logOR/sd from full-precision rebuilt table',round(mx,4))
