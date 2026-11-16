import sys, json, numpy as np, pandas as pd
out = json.load(open("theories.json"))
t5 = pd.DataFrame(out["T5"]); t5["below_p05"] = t5.cosine < t5.null_p05
print(t5.groupby(["antagonist", "agonist"]).agg(cells=("cell", "size"), negative=("cosine", lambda x: int((x < 0).sum())), below_p05=("below_p05", "sum"),
      cos_min=("cosine", "min"), cos_median=("cosine", "median"), cos_max=("cosine", "max")).to_string())
raw = pd.read_csv(sys.argv[1], index_col=0)
print("num_DE_pbs_wells:", raw.num_DE_pbs_wells.describe().round(2).to_dict()); print(raw.num_DE_pbs_wells.value_counts().sort_index().head(14).to_dict())
d = raw[(raw.adj_p_value < 0.05) & (raw.log_fc.abs() > 0.5)]
# sign test on shared genes, antagonist vs agonist
for ant, ago in (("IL-1Ra", "IL-1-beta"), ("IL-1Ra", "IL-1-alpha"), ("IL-36Ra", "IL-36-alpha")):
    tot = opp = 0
    for c in ["CD14_Mono", "CD16_Mono", "cDC", "CD4_Naive_T_cell", "CD8_Naive_T_cell", "Treg", "NK_CD56low", "Naive_B_cell"]:
        a = d[(d.cytokine == ant) & (d.celltype == c)].set_index("gene").log_fc; g = d[(d.cytokine == ago) & (d.celltype == c)].set_index("gene").log_fc
        sh = a.index.intersection(g.index); tot += len(sh); opp += int((np.sign(a[sh]) != np.sign(g[sh])).sum())
    print(ant, "vs", ago, "shared significant genes", tot, "opposite sign", opp, "share", round(opp / max(tot, 1), 2))
# what do genes shared in opposite sign look like in monocytes?
a = d[(d.cytokine == "IL-1Ra") & (d.celltype == "CD14_Mono")].set_index("gene").log_fc; g = d[(d.cytokine == "IL-1-beta") & (d.celltype == "CD14_Mono")].set_index("gene").log_fc
sh = a.index.intersection(g.index); o = sh[np.sign(a[sh]) != np.sign(g[sh])]
print("CD14 mono, IL-1Ra down & IL-1b up:", sorted([x for x in o if a[x] < 0], key=lambda x: a[x])[:25])
print("CD14 mono, IL-1Ra up & IL-1b down:", sorted([x for x in o if a[x] > 0], key=lambda x: -a[x])[:15])
