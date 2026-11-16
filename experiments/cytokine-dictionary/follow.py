import sys, numpy as np, pandas as pd
sys.argv=[sys.argv[0]]+sys.argv[1:]
d=pd.read_csv(sys.argv[1],index_col=0); d=d[(d.adj_p_value<0.05)&(d.log_fc.abs()>0.5)]
LY=["CD4_Naive_T_cell","CD8_Naive_T_cell","CD4_Memory_T_cell","CD8_Memory_T_cell","NK_CD56low","Naive_B_cell"]; MY=["CD14_Mono","CD16_Mono"]
INERT=["EPO","TPO","GDNF","PSPN","SCF","IL-5","HGF","VEGF","EGF"]
n=d.groupby(['cytokine','celltype']).size()
ne=lambda q,c:int(n.get((q,c),0))
inert=[ne(q,c) for q in INERT for c in LY+MY]; print('inert genes: median',np.median(inert),'95th',np.percentile(inert,95),'max',max(inert))
print('inert by cell',{c:[ne(q,c) for q in INERT] for c in LY+MY})
print('IL-18 present?', 'IL-18' in set(d.cytokine))
vec={k:g.set_index('gene').log_fc for k,g in d.groupby(['cytokine','celltype'])}
def cos(a,b):
    if a is None or b is None: return np.nan
    u=a.index.union(b.index); x,y=a.reindex(u).fillna(0).values,b.reindex(u).fillna(0).values; return float(x@y/np.sqrt((x@x)*(y@y)))
for q in ['IL-2','IL-15','IL-7','IL-21','IL-12']:
    ifn=d[(d.cytokine==q)&(d.gene=='IFNG')][['celltype','log_fc']].sort_values('log_fc',ascending=False).head(4)
    print(q,'IFNG induced in:',[(r.celltype,round(r.log_fc,1)) for r in ifn.itertuples()])
    for c in MY:
        sims={s:cos(vec[(q,c)],vec.get((s,c))) for s in d.cytokine.unique() if s!=q and (s,c) in vec and len(vec[(s,c)])>=20}
        top=sorted(sims.items(),key=lambda x:-x[1])[:6]
        print('   ',c,'cos to IFN-gamma %.2f'%sims.get('IFN-gamma',np.nan),'| top:',[(s,round(v,2)) for s,v in top])
# receptor chains measured in the 'negative' monocytes?
e=pd.read_csv(sys.argv[1],index_col=0).groupby(['gene','celltype']).logCPM.median()
for g in ['IL2RB','IL2RG','IL2RA','IL15RA','IL12RB1','IL12RB2','IL7R','IL21R','IFNGR1','IFNGR2','CSF2RA','CSF2RB','IL1R1']:
    print(g,{c:round(float(e[(g,c)]),1) for c in ['CD14_Mono','CD16_Mono','NK_CD56low','CD4_Naive_T_cell'] if (g,c) in e.index})
# GM-CSF -> lymphocytes: IL1 induction and similarity
x=d[(d.cytokine=='GM-CSF')&d.gene.isin(['IL1A','IL1B','IL6','TNF'])][['celltype','gene','log_fc']]; print(x.sort_values('log_fc',ascending=False).head(8).to_string())
