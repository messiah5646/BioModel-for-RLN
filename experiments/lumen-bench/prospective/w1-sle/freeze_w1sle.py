"""RT-20261126-1 run 2: freeze prior-only forecasts for the 8 PREREG-W1-SLE trials with the SHIPPED lumen_forecast.py.
Run from experiments/lumen-bench:  python3 -I prospective/w1-sle/freeze_w1sle.py
No outcome data are read. Research and simulation only."""
import sys,os,json,hashlib,datetime,subprocess
sys.path.insert(0,os.getcwd())
import lumen_forecast as L
T={'sri4':'sle/sri4_table.csv','bicla':'sle/bicla_table.csv'}
C=[('NCT05039840','frexalimab APATURA ph2','frexalimab','CD40L','sri4',24),
   ('NCT06559163','obexelimab SunStone ph2','obexelimab','CD19/FcgRIIb','bicla',24),
   ('NCT05672576','cenerimod OPUS-1 ph3','cenerimod','S1P1','sri4',52),
   ('NCT05648500','cenerimod OPUS-2 ph3','cenerimod','S1P1','sri4',52),
   ('NCT05843643','upadacitinib SELECT-SLE ph3','upadacitinib','JAK','bicla',52),
   ('NCT05639114','ianalumab SIRIUS-SLE 1 ph3','ianalumab','BAFF-R','sri4',60),
   ('NCT05624749','ianalumab SIRIUS-SLE 2 ph3','ianalumab','BAFF-R','sri4',60),
   ('NCT06015230','GR1603 ph1b/2','GR1603','unknown','bicla',24)]
md5=lambda p:hashlib.md5(open(p,'rb').read()).hexdigest()
R={k:L.load(v) for k,v in T.items()}
out={'_meta':{'prereg':'experiments/lumen-bench/prospective/w1-sle/PREREG-W1-SLE.md','prereg_commit':'5515355',
  'frozen_utc':datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ'),
  'git_head':subprocess.run(['git','rev-parse','--short','HEAD'],capture_output=True,text=True).stdout.strip(),
  'lumen_forecast_md5':md5('lumen_forecast.py'),'tables_md5':{k:md5(v) for k,v in T.items()},
  'quantity':'log OR drug vs placebo, registered primary endpoint/visit, 0.5 cc; interval = m +/- 1.96 sd',
  'models':'M1 = forecast(anchor=False) [headline under test]; M0 = 0; M3 = forecast(anchor=True) [secondary]'}}
for nct,name,mol,tgt,tab,wk in C:
    m1=L.forecast(mol,tgt,R=R[tab],anchor=False);m3=L.forecast(mol,tgt,R=R[tab],anchor=True)
    for f in (m1,m3):
        mu,sd=f['prior_logOR'],f['prior_sd'];f['logOR_95']=[round(mu-1.96*sd,3),round(mu+1.96*sd,3)]
    out[nct]=dict(trial=name,molecule=mol,target=tgt,endpoint=('SRI-4' if tab=='sri4' else 'BICLA'),primary_week=wk,
                  table=T[tab],M0_null={'logOR':0.0},M1_pooled=m1,M3_anchored=m3,
                  C0_slot=('SRI-4 wk36 -> wk60 (descriptive, rule 8)' if mol=='ianalumab' else None))
    print(nct,name,tab,'M1',m1['prior_logOR'],m1['prior_sd'],m1['OR'],m1['OR_95'],'| M3',m3['anchor'],m3['OR'],m3['OR_95'])
json.dump(out,open('prospective/w1-sle/frozen_2026-11-26.json','w'),indent=1)
