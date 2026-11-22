"""Q17 / round 107: rebuild resp_table.csv (CRR responders per arm per visit) from the committed registry snapshot.
Usage: python3 -I build_resp_table.py <lupus_with_results.json.gz> <out.csv>
Selection rules are generic (regexes below); trial-specific entries are listed in ARMS (group-title -> arm label,
unavoidable labelling) and ALIAS (trials whose CRR measure is not titled 'complete'). See PREREG-Q17.md."""
import gzip,json,re,sys,math,csv
NCT={'ACCESS':'NCT00774852','CFZ533':'NCT03610516','PALIZADE':'NCT05781750','BELONG':'NCT00626197','AURA-LV':'NCT02141672',
 'ALXN2050':'NCT05097989','SANCTUARY':'NCT04564339','ORCHID-LN':'NCT04376827','TULIP-LN1':'NCT02547922','SELUNE':'NCT04181762',
 'BMS-LN':'NCT00430677','NOBILITY':'NCT02550652','AURORA':'NCT03021499','LUNAR':'NCT00282347','ATLAS':'NCT01499355',
 'BI655064':'NCT02770170','REGENCY':'NCT04221477','CALIBRATE':'NCT02260934','ALLURE':'NCT01714817','BLISS-LN':'NCT01639339'}
TARGET={'abatacept':'CD80/86','iscalimab':'CD40','zetomipzomib':'immunoproteasome','ocrelizumab':'CD20','voclosporin':'calcineurin',
 'ALXN2050':'factor D','ravulizumab':'C5','guselkumab':'IL-23','anifrolumab':'IFNAR','secukinumab':'IL-17A','obinutuzumab':'CD20',
 'rituximab':'CD20','BIIB023':'TWEAK','BI655064':'CD40','belimumab':'BAFF'}
FAST={'voclosporin','ravulizumab','ALXN2050','BIIB023'}           # PREREG-r86 labels
# trial -> list of (regex on group title, arm label, molecule); 'ctrl' marks the control arm
ARMS={'ACCESS':[(r'^Abatacept$','abatacept','abatacept'),(r'^Placebo$','ctrl',None)],
 'CFZ533':[(r'^CFZ533','iscalimab','iscalimab'),(r'^Placebo','ctrl',None)],
 'PALIZADE':[(r'30 mg','zeto_30','zetomipzomib'),(r'60 mg','zeto_60','zetomipzomib'),(r'^Placebo','ctrl',None)],
 'BELONG':[(r'^OCR 400','ocre_400','ocrelizumab'),(r'^OCR 1000','ocre_1000','ocrelizumab'),(r'^Placebo','ctrl',None)],
 'AURA-LV':[(r'Low Dose','vocl_low','voclosporin'),(r'High Dose','vocl_high','voclosporin'),(r'^Placebo$','ctrl',None)],
 'ALXN2050':[(r'^LN Cohort: ALXN2050 120','fD_120','ALXN2050'),(r'^LN Cohort: ALXN2050 180','fD_180','ALXN2050'),(r'^LN Cohort: Placebo$','ctrl',None)],
 'SANCTUARY':[(r'^LN Cohort: Ravulizumab','ravulizumab','ravulizumab'),(r'^LN Cohort: Placebo','ctrl',None)],
 'ORCHID-LN':[(r'^Guselkumab','guselkumab','guselkumab'),(r'^Placebo','ctrl',None)],
 'TULIP-LN1':[(r'Basic','anifro_basic','anifrolumab'),(r'Intensified','anifro_int','anifrolumab'),(r'^Placebo$','ctrl',None)],
 'SELUNE':[(r'^Secukinumab','secukinumab','secukinumab'),(r'^Placebo','ctrl',None)],
 'BMS-LN':[(r'^Abatacept 30/10','abata_30','abatacept'),(r'^Abatacept 10/10','abata_10','abatacept'),(r'^Placebo$','ctrl',None)],
 'NOBILITY':[(r'^OBINUTUZUMAB','obinutuzumab','obinutuzumab'),(r'^PLACEBO','ctrl',None)],
 'AURORA':[(r'^Voclosporin','voclosporin','voclosporin'),(r'^Placebo','ctrl',None)],
 'LUNAR':[(r'^Rituximab','rituximab','rituximab'),(r'^Placebo','ctrl',None)],
 'ATLAS':[(r'3 mg','biib_3','BIIB023'),(r'20 mg','biib_20','BIIB023'),(r'^Placebo','ctrl',None)],
 'BI655064':[(r'^120','bi_120','BI655064'),(r'^180','bi_180','BI655064'),(r'^240','bi_240','BI655064'),(r'^Placebo','ctrl',None)],
 'REGENCY':[(r'^Obinutuzumab','obinutuzumab','obinutuzumab'),(r'^Placebo','ctrl',None)],
 'CALIBRATE':[(r'\(RCB\)','belimumab','belimumab'),(r'\(RC\)','ctrl',None)],
 'ALLURE':[(r'^Abatacept IV','abatacept','abatacept'),(r'^Placebo IV','ctrl',None)],
 'BLISS-LN':[(r'^Belimumab','belimumab','belimumab'),(r'^Placebo','ctrl',None)]}
INC=re.compile(r'complete (renal )?(response|remission)|\bCRR\b|\bCR\b',re.I)
ALIAS={'AURORA':re.compile(r'adjudicated renal response at week|^number of participants with renal response at week',re.I)}  # AURORA's protocol 'renal response' = CRR
EXC=re.compile(r'time to|duration|number of months|incidence|sustain|maintain|through|within those|steroid|prednison|corticosteroid|'
               r'nephrotic|modified|mcrr|creatinine criteria|prr or|or prr|or crr|\(ri\) or|access definition|remaining|low dose',re.I)
CLS_OK=re.compile(r'^(CRR|CR\b|achieved response|responder|number of renal responders|week|day)',re.I)
CLS_BAD=re.compile(r'PRR|^PR\b|^NR|non-?responder|did not|PR -|NR -',re.I)
def weeks(txt):
    m=re.search(r'(week|day|month)s?\s*(\d+)',txt,re.I)
    if not m: return None
    u,v=m.group(1).lower(),int(m.group(2))
    return v if u=='week' else round(v/7) if u=='day' else round(v*52/12)
def half_up(x): return int(math.floor(x+0.5))
D={s['protocolSection']['identificationModule']['nctId']:s for s in json.load(gzip.open(sys.argv[1]))}
cand=[]   # (trial, group title, week, r, n, measure idx, type, title)
for t,n in NCT.items():
    for i,m in enumerate(D[n]['resultsSection']['outcomeMeasuresModule']['outcomeMeasures']):
        ttl=m['title']
        # correction 1 (after first comparison, reported): 'partial' is NOT a title exclusion; ranked CR/PR/NR outcomes
        # are kept and the PR/NR classes are dropped by CLS_BAD. Titles naming only partial response never match INC.
        if not (INC.search(ttl) or (t in ALIAS and ALIAS[t].search(ttl))) or EXC.search(ttl): continue
        cls_titles=[c.get('title','') for c in m.get('classes',[])]
        visit_in_class=any(weeks(c) for c in cls_titles)
        if re.search(r'\bduring\b',ttl,re.I) and not re.search(r'\bat (week|day|month)',ttl,re.I) and not visit_in_class: continue  # cumulative 'ever'
        unit=(m.get('unitOfMeasure') or '').lower();pct='percent' in unit
        if m.get('paramType') not in ('COUNT_OF_PARTICIPANTS','NUMBER'): continue
        g={x['id']:x['title'] for x in m.get('groups',[])}
        den={c['groupId']:int(c['value']) for d in m.get('denoms',[]) for c in d['counts']}
        for c in m.get('classes',[]):
            ct=c.get('title','')
            if ct and (CLS_BAD.search(ct) or not CLS_OK.search(ct)): continue
            wk=weeks(ct) or weeks(ttl) or weeks(m['timeFrame'])
            for cat in c['categories']:
                if cat.get('title') and CLS_BAD.search(cat['title']): continue
                for me in cat['measurements']:
                    v=me.get('value')
                    if v in (None,'NA'): continue
                    N=den[me['groupId']];r=half_up(float(v)*N/100) if pct else int(float(v))
                    cand.append((t,g[me['groupId']],wk,r,N,i,m['type'],ttl))
# map groups to arms; per (trial, arm, week) keep largest denominator, then PRIMARY, then first measure
best={}
for t,gt,wk,r,N,i,ty,ttl in cand:
    lab=next(((a,mol) for rx,a,mol in ARMS[t] if re.search(rx,gt)),None)
    if lab is None: continue
    k=(t,lab[0],wk);score=(N,ty=='PRIMARY',-i)
    if k not in best or score>best[k][0]: best[k]=(score,r,N,lab[1],i,ttl)
# drop visits where any arm's denominator is < 50% of that arm's largest denominator (incomplete-follow-up subsets)
maxN={};[maxN.__setitem__((t,a),max(maxN.get((t,a),0),v[2])) for (t,a,w),v in best.items()]
rows=[];dropped=[]
for (t,a,w),v in sorted(best.items()):
    if a=='ctrl': continue
    c=best.get((t,'ctrl',w))
    if c is None: continue
    if v[2]<0.5*maxN[(t,a)] or c[2]<0.5*maxN[(t,'ctrl')]: dropped.append((t,a,w,v[2],c[2]));continue
    mol=v[3];rows.append(dict(trial=t,nct=NCT[t],arm=a,molecule=mol,target=TARGET[mol],cls='fast' if mol in FAST else 'slow',week=w,
        r_drug=v[1],n_drug=v[2],r_ctrl=c[1],n_ctrl=c[2],measure_drug=v[4],measure_ctrl=c[4],source='ctgov'))
with open(sys.argv[2],'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(len(rows),'rows,',len({(r['trial'],r['arm']) for r in rows}),'arms,',len({r['trial'] for r in rows}),'trials; dropped (<50% denominators):',dropped)
