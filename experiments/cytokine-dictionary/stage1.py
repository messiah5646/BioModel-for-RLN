"""Stage 1 on real data: do bystander effects exist in the human cytokine screen?

Input: the significant-gene table shipped in github.com/theislab/HumanCytokineDict
(edgeR, cytokine vs PBS per cell type, adj p < 0.1 and |log2FC| >= 0.5, 12 donors pooled).

The receptor lists below were written from textbook receptor biology BEFORE any
pair-level response in the table was read. They are a prior, not a measurement.
"""
import sys, json
import numpy as np, pandas as pd

LY = ["CD4_Naive_T_cell", "CD8_Naive_T_cell", "CD4_Memory_T_cell", "CD8_Memory_T_cell", "NK_CD56low", "Naive_B_cell"]
MY = ["CD14_Mono", "CD16_Mono"]
# cytokine -> (cell types expected to LACK the receptor, cell types expected to respond directly)
PRIOR = {
    "IL-2":   (MY, ["NK_CD56hi", "NK_CD56low", "Treg"]),
    "IL-15":  (MY, ["NK_CD56hi", "NK_CD56low", "CD8_Memory_T_cell"]),
    "IL-12":  (MY, ["NK_CD56hi", "NK_CD56low"]),
    "IL-18":  (MY + ["Naive_B_cell"], ["NK_CD56hi", "NK_CD56low", "MAIT"]),
    "IL-7":   (MY + ["NK_CD56low", "Naive_B_cell"], ["CD4_Naive_T_cell", "CD8_Naive_T_cell"]),
    "IL-21":  (MY, ["Naive_B_cell", "NK_CD56low"]),
    "M-CSF":  (LY, ["CD14_Mono"]),
    "IL-34":  (LY, ["CD14_Mono"]),
    "GM-CSF": (LY, ["CD14_Mono", "cDC"]),
    "G-CSF":  (LY, ["CD14_Mono"]),
    "C5a":    (LY, ["CD14_Mono"]),
    "IL-3":   (["CD4_Naive_T_cell", "CD8_Naive_T_cell", "NK_CD56low"], ["pDC", "CD14_Mono"]),
    "TSLP":   (["CD4_Naive_T_cell", "CD8_Naive_T_cell", "NK_CD56low", "Naive_B_cell"], ["cDC", "CD14_Mono"]),
    "CD40L":  (["CD4_Naive_T_cell", "CD8_Naive_T_cell", "CD4_Memory_T_cell", "CD8_Memory_T_cell", "NK_CD56low"], ["Naive_B_cell", "CD14_Mono"]),
    "BAFF":   (["CD14_Mono", "NK_CD56low", "CD4_Naive_T_cell", "CD8_Naive_T_cell"], ["Naive_B_cell"]),
    "APRIL":  (["CD14_Mono", "NK_CD56low", "CD4_Naive_T_cell", "CD8_Naive_T_cell"], ["Naive_B_cell"]),
}
# cytokines whose receptors are not expected on blood mononuclear cells: the null rate
INERT = ["EPO", "TPO", "GDNF", "PSPN", "SCF", "IL-5", "HGF", "VEGF", "EGF"]
PADJ, LFC, MIN_GENES, COS = 0.05, 0.5, 20, 0.3      # fixed before reading responses

d = pd.read_csv(sys.argv[1], index_col=0)
ci = pd.read_excel(sys.argv[2], sheet_name="all_cytokines").set_index("name")
d = d[(d.adj_p_value < PADJ) & (d.log_fc.abs() > LFC)]
cells = LY + MY + ["NK_CD56hi", "Treg", "MAIT", "cDC", "pDC"]
vec = {k: g.set_index("gene").log_fc for k, g in d.groupby(["cytokine", "celltype"])}
nde = lambda q, c: len(vec.get((q, c), ()))
energy = lambda q, c: float((vec[(q, c)] ** 2).sum()) if (q, c) in vec else 0.0
def cos(a, b):
    if a is None or b is None or len(a) == 0 or len(b) == 0: return np.nan
    u = a.index.union(b.index); x, y = a.reindex(u).fillna(0).values, b.reindex(u).fillna(0).values
    return float(x @ y / np.sqrt((x @ x) * (y @ y)))

out = {}
# measured contradiction check: where the table has logCPM for a receptor chain in a "negative" cell type
expr = d.groupby(["gene", "celltype"]).logCPM.median()
# ---------------- A: existence
by = [(q, c) for q, (neg, _) in PRIOR.items() for c in neg]
di = [(q, c) for q, (_, pos) in PRIOR.items() for c in pos]
nul = [(q, c) for q in INERT for c in LY + MY]
frac = lambda pairs: float(np.mean([nde(q, c) >= MIN_GENES for q, c in pairs]))
out["A"] = dict(n_bystander_pairs=len(by), share_bystander_responding=frac(by),
                n_direct_pairs=len(di), share_direct_responding=frac(di),
                n_inert_pairs=len(nul), share_inert_responding=frac(nul),
                median_genes=dict(bystander=float(np.median([nde(*p) for p in by])), direct=float(np.median([nde(*p) for p in di])),
                                  inert=float(np.median([nde(*p) for p in nul]))))
# ---------------- B and C: magnitude, and variation by cell type
rows = []
for q, (neg, pos) in PRIOR.items():
    nb, nd = [nde(q, c) for c in neg], [nde(q, c) for c in pos]
    eb, ed = [energy(q, c) for c in neg], [energy(q, c) for c in pos]
    rows.append(dict(cytokine=q, genes_bystander_median=float(np.median(nb)), genes_bystander_max=int(max(nb)),
                     genes_direct_median=float(np.median(nd)), size_ratio=float(np.sqrt(np.median(eb) / np.median(ed))) if np.median(ed) > 0 else np.nan,
                     bystander_by_cell={c: nde(q, c) for c in neg}, direct_by_cell={c: nde(q, c) for c in pos}))
out["B"] = rows
# ---------------- D: are bystander responses explained by a measured, induced cytokine?
lig = {s: [g.strip() for g in str(ci.loc[s, "gene"]).split(",")] for s in d.cytokine.unique() if s in ci.index}
drows = []
for q, c in by:
    if nde(q, c) < MIN_GENES: continue
    induced = {}
    for s, genes in lig.items():
        if s == q: continue
        hit = d[(d.cytokine == q) & d.gene.isin(genes) & (d.log_fc > LFC)]
        if len(hit): induced[s] = dict(sender=hit.sort_values("log_fc").iloc[-1].celltype, log2fc=float(hit.log_fc.max()))
    sims = {s: cos(vec[(q, c)], vec.get((s, c))) for s in lig if s != q and nde(s, c) >= MIN_GENES}
    cand = {s: v for s, v in sims.items() if s in induced}; ctrl = {s: v for s, v in sims.items() if s not in induced}
    best = max(cand, key=cand.get) if cand else None
    top = max(sims, key=sims.get) if sims else None
    # ambient-RNA control: similarity to the direct responders' own response to q
    amb = max([cos(vec[(q, c)], vec.get((q, r))) for r in PRIOR[q][1]] + [np.nan], key=lambda x: -1 if np.isnan(x) else x)
    drows.append(dict(cytokine=q, cell=c, genes=nde(q, c), n_induced=len(induced), best_mediator=best,
                      best_cos=cand.get(best) if best else None, mediator_sender=induced[best]["sender"] if best else None,
                      mediator_log2fc=induced[best]["log2fc"] if best else None,
                      top_overall=top, top_overall_cos=sims.get(top) if top else None, top_is_induced=bool(top in induced) if top else None,
                      ctrl_share_above=float(np.mean([v >= COS for v in ctrl.values()])) if ctrl else None,
                      cand_share_above=float(np.mean([v >= COS for v in cand.values()])) if cand else None,
                      cos_to_sender_response=None if np.isnan(amb) else amb))
D = pd.DataFrame(drows)
out["D"] = dict(n_significant_bystander=len(D),
                share_with_matching_induced_mediator=float((D.best_cos >= COS).mean()) if len(D) else None,
                share_top_match_is_induced=float(D.top_is_induced.mean()) if len(D) else None,
                mean_share_candidates_matching=float(D.cand_share_above.mean()) if len(D) else None,
                mean_share_controls_matching=float(D.ctrl_share_above.mean()) if len(D) else None,
                share_mediator_beats_sender_similarity=float((D.best_cos > D.cos_to_sender_response.fillna(-1)).mean()) if len(D) else None)
json.dump(out, open(sys.argv[3], "w"), indent=1)
D.to_csv(sys.argv[3].replace(".json", "_D.csv"), index=False)
print(json.dumps(out["A"], indent=1))
print(pd.DataFrame(rows)[["cytokine", "genes_bystander_median", "genes_bystander_max", "genes_direct_median", "size_ratio"]].to_string())
print(json.dumps(out["D"], indent=1))
pd.set_option("display.width", 250)
print(D[["cytokine", "cell", "genes", "n_induced", "best_mediator", "best_cos", "mediator_sender", "top_overall", "top_overall_cos", "top_is_induced", "cos_to_sender_response"]].round(2).to_string())
