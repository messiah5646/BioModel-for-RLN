# RESULT round 103 (Q38, K1): KPMP -> E-MTAB-13596 kidney immune mapping (descriptive)

PREREG: PREREG-K1-round103.md (ca04f95). Scripts: k1_prep.py, k1_map.py, k1_describe.py. Outputs: k1_val.json, k1_query_qc.csv,
k1_query_obs.csv.gz, k1_describe.json. GEO search: K1-GEO-CATALOGUE-round103.md.

## Held-out gate (KPMP, 15 of 77 donors held out, 7,393 immune cells)
| Method | Accuracy | Macro-F1 (15 immune labels) |
|---|---|---|
| Majority class (T) | 0.24 | 0.03 |
| kNN-15 on PCA-50 | 0.85 | 0.79 |
| scANVI (scArches query) | 0.85 | 0.83 |
Gate PASS (scANVI >= kNN - 0.05 and > majority): scANVI labels primary. The gain over kNN is in rare/proliferating labels
(cycMNP 0.90 vs 0.44, cycT 0.95 vs 0.79); common labels tie. ncMON is the weakest label for both (0.53).

## Query (7 CD45+ samples, 34,993 cells; 1,955 flagged by QC or scrublet; protected rule returned 14 flagged label groups (incl. plasma cells), 33,972 kept)
- scANVI is over-confident across labs: median max-probability 0.97-0.99 in every sample, yet 9-16 % of cells are "unresolved"
  (> half of 15 reference neighbours disagree), and scANVI agrees with kNN-PCA on only 75-87 % of cells. The same pattern in the
  nephrectomy control (10 % unresolved) shows this is a cross-lab effect, not disease.
- Reference gap caught by the marker check: cells labelled "ncMON" are CSF3R/FCGR3B/CXCR2/S100A8-high and FCGR3A-low, i.e.
  neutrophils. KPMP has no neutrophil label, so they were forced into the nearest myeloid class.
- AGN1 is 85 % proximal-tubule epithelium (sort failure; median mito 24 %); the non-immune catch-all labels prevented these from
  being called immune cells.
- LN (n = 1, class II): monocytes are 51 % of immune cells (AGN 1-23 %, control 9 %) and T cells 6 % (AGN 22-43 %, control 25 %);
  1,469 of 1,541 scANVI "MON" calls are also MON by kNN-PCA. LN myeloid cells are FCN1-high and CD163/C1QA/APOE/SPP1-low
  (infiltrating-monocyte-like, not resident macrophage). IFN6 score: myeloid +0.44 in LN vs -0.11 to -0.18 in all others; T/NK
  +0.22 vs -0.10 to -0.21.
- Round-91 interstitial non-responder programme in myeloid cells: LN 0.64, control 0.61, AGN 0.19-0.59. Not LN-specific here;
  it tracks myeloid composition more than disease. Cannot test the round-91 claim (one patient, no response label).

## Systematic search
No open 2024-2026 human LN kidney single-cell/nucleus series with > 1 LN patient. Best open LN kidney expression remains
bulk/ROI: GSE200306, GSE32591/2, and GSE263909 (GeoMx glomeruli, 33 LN vs 16 control ROIs, counts open; new candidate for Q19/Q34).
