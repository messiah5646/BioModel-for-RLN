# Round 103 (Q38): open 2024-2026 human LN kidney single-cell / spatial series

Search 2026-11-22. (1) data/real/catalogues/geo_lupus_series.csv (597 lupus series from the Q0 harvest), filtered for
kidney/renal/nephritis/glomerul/biopsy AND single-cell/nucleus/spatial words, 2024-2026: 12 hits. (2) GEO E-utilities (db=gds):
("lupus nephritis" OR lupus AND kidney/renal) AND (single cell/nucleus/snRNA/scRNA/spatial) AND gse AND Homo sapiens AND
PDAT 2024-2026: 17 hits. Sample characteristics read from GEO SOFT brief for each candidate.

| Accession | Date | What | LN kidney samples | Labels | Open matrix | Use |
|---|---|---|---|---|---|---|
| E-MTAB-13596 (ArrayExpress) | 2024-06 | CD45+ kidney immune scRNA, 10x 3' v3 | 1 (class II) + 5 ANCA-GN + 1 nephrectomy control | class, disease | yes (MTX) | K1 query (this round) |
| GSE263909 | 2024-06 | GeoMx glomerular ROIs (spatial, not single-cell) | 33 LN ROIs vs 16 control | disease only in SOFT | yes (Q3-normalised counts) | Q19/Q34 candidate: glomerular IFN vs control; check paper supplement for class/activity |
| GSE333529 | 2026-09 | Spatial, "skin-kidney neutrophil axis", healthy and lupus kidneys | <= 4 sections total | tissue only in SOFT | RAW tar | descriptive only |
| GSE294965 / GSE303481 | 2025 | Crescent formation, spatial + snRNA | 0 LN (ANCA-GN 11, control 4 in GSE303481) | ANCA-GN | h5ad | not LN; possible crescent comparator |
| GSE264393 / GSE264334 | 2024-04 | snRNA + Xenium method paper | 0 LN (control, IgAN) | - | RAW | not LN |
| GSE322527 | 2026-07 | urine vs blood T cells at LN flare | 0 kidney (urine) | flare | yes | Q42 |
| GSE254176, GSE306378, GSE298578/85773/85774, GSE287690 | 2025-26 | blood only | 0 kidney | - | - | not kidney |
| GSE328870, GSE334922-4 | 2026 | kidney transplant rejection | 0 LN | - | - | not LN |

Conclusion (Supported): in fully open repositories there is still NO 2024-2026 human LN kidney single-cell/nucleus dataset with
more than one LN patient. The only open LN kidney single-cell matrix found is the single class II patient in E-MTAB-13596.
AMP Phase I/II kidney scRNA remain controlled (ImmPort/Synapse). The largest open LN kidney expression resources remain bulk or
spatial-ROI (GSE200306 serial biopsies, GSE32591/2, GSE263909 GeoMx).
