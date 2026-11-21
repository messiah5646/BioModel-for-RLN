import json,pandas as pd;from forecast_prior import prior
R=pd.read_csv('arm_table_reconstructed.csv')
C=[('NCT05138133','IRIS anifrolumab ph3','anifrolumab','IFNAR',None),
('NCT06406205','QLG1074 (voclosporin) ph3 China','voclosporin','calcineurin',None),
('NCT05039619','obinutuzumab paediatric ph2','obinutuzumab','CD20',None),
('NCT07044115','MIL62 (anti-CD20) ph1/2','MIL62','CD20',None),
('NCT05126277','SIRIUS-LN ianalumab ph3','ianalumab','BAFF-R','BAFF'),
('NCT05201469','VIBRANT dazodalibep (VIB4920) ph2','dazodalibep','CD40L','CD40'),
('NCT04643470','zanubrutinib ph2','zanubrutinib','BTK',None),
('NCT05810948','efgartigimod ph2','efgartigimod','FcRn',None),
('NCT05268289','iptacopan ph2 (terminated)','iptacopan','factor B',None),
('NCT07363460','HSK39297 ph2','HSK39297','unknown',None),
('NCT07229742','SHR-2173 ph2','SHR-2173','unknown',None),
('NCT07038382','HLX79 + HLX01 ph2','HLX79','unknown',None),
('NCT07201129','cenerimod ph3','cenerimod','S1P1',None)]
out={}
for nct,name,mol,tgt,axis in C:
    f=prior(mol,tgt,R);f={k:(float(v) if hasattr(v,'item') else v) for k,v in f.items()};f.update(trial=name,molecule=mol,target=tgt)
    if axis: f['secondary_axis_anchor']={k:(float(v) if hasattr(v,'item') else v) for k,v in prior(mol,axis,R).items()};f['secondary_axis_anchor']['axis']=axis
    out[nct]=f;print(nct,name,f['anchor'],f['logOR'],f['sd'],f['OR'],f['OR_95'],f.get('secondary_axis_anchor',{}).get('OR'))
json.dump(out,open('w1_frozen_2026-11-21.json','w'),indent=1)
