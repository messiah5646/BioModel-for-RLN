# Q57: round 91 regenerated from open GEO files (2026-11-23, round 113)

Round 91 committed no script (audit round 111). `regenerate_r91.py` rebuilds it from GSE200306 (series matrix or the
two normalized supplementary files; both give identical numbers). Choices round 91 did not record are marked
RECONSTRUCTED in the script. Output: regen_matrix.txt, regen_suppl.txt.

| Claim (round 91) | Original | Regenerated | Reproduces? |
|---|---|---|---|
| H1 CD68 compartment test | not testable (no CD68) | not testable | yes |
| Glomerular macrophage score, CR > NR at Bx1 | p 0.199 (CR 16, NR 5) | p 0.196 (CR 20, NR 7) | yes (still not significant) |
| TI macrophage score, CR vs NR at Bx1 | no difference (p 0.38/0.64) | no difference (p 0.57/0.45) | yes |
| H2 TI: kidney PCs add to clinical | FAIL (-0.010) | FAIL (-0.036) | yes |
| H2 glom: kidney PCs add to clinical | "PASS by the letter" (+0.062), called no usable signal | FAIL (+0.000) | **no** -> drop the "pass" |
| H3 NR > CR after treatment, TI (descriptive) | CD99, MAP4K4, CSF1R, PTPN6, TGFB1, HLA-DRA, IFI16, LILRB1, LILRB2, C1S | CSF1R #1, MAP4K4 #2, PTPN6 #5, C1S #6, CD99 #8, LILRB1 #10, IFI16 #11, TGFB1 #13, LILRB2 #19, HLA-DRA #41 of 502 | mostly (6/10 in top 10; all in top 41) |
| H3 CR > NR | IL1RN, ENTPD1, RAG2, CFP, CASP10 | IL1RN, CFP, ENTPD1, CCL18, KCNJ2 | mostly (3/5) |

Sample counts differ (regenerated: more CR/NR patients matched, all 94 patient-biopsies labelled). The original medians
(~2.2) are on a different scale from the GEO log2 values (~4.6), so round 91 used an unrecorded normalisation.

**Decision (Q57):** cite round 91 in the preprint only via this script's numbers. Keep: "after treatment, NR tubulointerstitium
keeps a myeloid programme (CSF1R, LILRB1/2, C1S, PTPN6) with TGFB1 elevated", Emerging (NR n 9 vs CR n 24, one cohort,
descriptive, no independent replication: Q12 blocked). Drop: the glomerular H2 "pass".
