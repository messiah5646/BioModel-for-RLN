"""
Is the literature-triage task (the only place a Jev-class decision model could
plug in) hard enough to need a learned model?

Round 41 barred Geneformer for losing to PCA. Same standard here: before admitting
a decision model, run the trivial comparator. For "does this paper's model include
compartment X as a state variable?", the trivial comparator is a keyword match.

IMPORTANT LIMITATION, stated up front: these are text snippets as recorded in our
own round docs (R37/R41/R42), not verbatim abstracts fetched now. n=9. This measures
the STRUCTURE of the task -- whether keyword presence determines the answer -- not
real-world accuracy. It cannot support an accuracy claim about any model.

The labels are our own hand determinations from rounds 37/40/42.
"""
import re, json

# (id, snippet as recorded in our round docs, hand label: does the MODEL carry
#  long-lived / CD19-negative plasma cells as a STATE VARIABLE?)
CASES = [
 ("Budu-Grajdeanu2010",
  "Four ODEs: immune complexes, pro-inflammatory mediators, damaged tissue D, "
  "anti-inflammatory mediators. D carries a decay term lumping repair, resolution "
  "and regeneration. Calibrated per-patient on four Ohio SLE Study patients.",
  False),
 ("Park2025",
  "Two compartments, blood and bone marrow. CAR-T split CD4/CD8 across four stages; "
  "reactive vs non-reactive CD19+ B cells. B-cell repopulation is zero-order "
  "naive-only. No kidney compartment; proteinuria driven by NETosis and IFN-alpha.",
  False),
 ("Volkova2024",
  "Seventeen nonlinear ODEs, 52 parameters, fitted to individual data from 254 "
  "anifrolumab patients. Readouts plasma IFN-alpha and IFNGS only. No organ compartment.",
  False),
 ("SimulationsPlusSLE",
  "SLE QSP model covering skin and joints, explicitly excluding lupus nephritis.",
  False),
 ("Schips2023",
  "Germinal centre agent-based model. Every somatic hypermutation can generate "
  "self-reactive clones; natural BCR redemption is insufficient. Plasma cell output "
  "is tracked as cells leave the germinal centre.",
  False),
 ("NETosis2024",
  "Modelling the impact of NETosis during the initial stage of SLE. ODE model of "
  "neutrophil extracellular traps and immune activation.",
  False),
 # THE TRAP. Round 42 retrieved this, it talks about long-lived plasma cells at
 # length, and it is NOT a mathematical model at all -- it is a review of MOUSE
 # models. A keyword matcher scores it as a hit and is wrong.
 ("PMC5410990_TRAP",
  "Modelling clinical systemic lupus erythematosus: similarities, differences and "
  "success stories. Cyclophosphamide does not delete long-lived plasma cells "
  "efficiently, and these cells persist in survival niches in the bone marrow and "
  "inflamed tissues. Plasmablasts, B cells, T follicular helper cells and myeloid "
  "cells are discussed.",
  False),
 # Our own model: genuinely has the compartment.
 ("OurModel_R40",
  "Six states: autoreactive B cells and short-lived plasma cells, both CD19+ and "
  "depletable; long-lived plasma cells, CD19-negative and not depletable; "
  "autoantibody titre; activity; monotone renal damage. dPl/dt = f_l*supp*B - d_Pl*Pl.",
  True),
 ("Nat Commun 2025 (experimental, not a model)",
  "Distinct autoreactive CD19-negative plasma cell subsets accumulate in "
  "lupus-prone mice. Flow cytometry and single-cell sequencing of splenic and "
  "bone marrow compartments.",
  False),
]

KEYWORDS = [r"long-lived plasma cell", r"CD19[- ]negative", r"CD19-", r"LLPC"]

def keyword_predict(text):
    t = text.lower()
    return any(re.search(k.lower(), t) for k in KEYWORDS)

print("TRIVIAL COMPARATOR: keyword match for the LLPC-as-state-variable question")
print(f"{'case':>34} | {'keyword':>8} | {'truth':>6} | result")
print("-"*72)
tp=fp=tn=fn=0
for cid, txt, truth in CASES:
    pred = keyword_predict(txt)
    if pred and truth: tp+=1; r="correct"
    elif pred and not truth: fp+=1; r="** FALSE POSITIVE **"
    elif not pred and truth: fn+=1; r="** FALSE NEGATIVE **"
    else: tn+=1; r="correct"
    print(f"{cid:>34} | {str(pred):>8} | {str(truth):>6} | {r}")
n=len(CASES)
print(f"\naccuracy {100*(tp+tn)/n:.1f}%  ({tp} TP, {tn} TN, {fp} FP, {fn} FN, n={n})")
print(f"precision {tp/max(tp+fp,1):.2f}  recall {tp/max(tp+fn,1):.2f}")
print()
print("WHY THE BASELINE FAILS, and it is not a tuning problem:")
print("  The question is not 'does the text mention X' but 'does the MODEL carry X")
print("  as a state variable'. Those come apart in two directions:")
print("   - a review discussing LLPC biology at length carries no model at all")
print("   - an experimental paper measuring CD19- plasma cells is not a model either")
print("  No keyword list separates these, because the discriminating information is")
print("  the grammatical role of the term, not its presence. Adding keywords cannot")
print("  fix a false positive caused by the term being present and correct.")
json.dump(dict(n=n,tp=tp,tn=tn,fp=fp,fn=fn,accuracy=(tp+tn)/n),
          open('results/keyword_baseline.json','w'), indent=1)
