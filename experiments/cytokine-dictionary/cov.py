import pandas as pd, numpy as np, sys
d=pd.read_csv(sys.argv[1],index_col=0); ci=pd.read_excel(sys.argv[2],sheet_name='all_cytokines').set_index('name')
# is logCPM a property of (gene, celltype)?
g=d.groupby(['gene','celltype']).logCPM.agg(['std','count','median'])
print('logCPM sd within (gene,celltype), where n>1: median',g.loc[g['count']>1,'std'].median(),'95pct',g.loc[g['count']>1,'std'].quantile(.95))
expr=g['median'].unstack('celltype')    # genes x celltypes, NaN = never DE there
cyts=sorted(d.cytokine.unique()); print('cytokines not in info:',[c for c in cyts if c not in ci.index])
print(ci.loc[[c for c in cyts if c in ci.index],['gene','receptor gene']].to_string())
rows=[]
for q in cyts:
    if q not in ci.index or pd.isna(ci.loc[q,'receptor gene']): continue
    chains=[x.strip() for x in str(ci.loc[q,'receptor gene']).replace(';',',').split(',') if x.strip()]
    for c in expr.columns:
        v=[expr.loc[ch,c] if ch in expr.index else np.nan for ch in chains]
        known=[x for x in v if not np.isnan(x)]
        cls='neg' if any(x<3 for x in known) else ('pos' if len(known)==len(chains) else 'unknown')
        rows.append((q,c,cls,len(known),len(chains)))
r=pd.DataFrame(rows,columns=['cyt','ct','cls','known','chains'])
print(r.cls.value_counts().to_string()); print(pd.crosstab(r.ct,r.cls).to_string())
print('logCPM distribution of all rows:',d.logCPM.describe().to_string())
