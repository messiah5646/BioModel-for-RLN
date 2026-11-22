# Q12 catalogue: independent kidney data to validate the round-91 H3 non-responder programme (2026-11-22)

Goal: an LN kidney expression dataset, independent of GSE200306 (Ohio State), with per-patient treatment-response
labels, ideally post-treatment biopsies. Searched GEO (E-utilities, "lupus nephritis" human series, 124 hits, all
kidney series with n>=15 read) and CELLxGENE (all collections mentioning lupus). No expression values were opened.

| Dataset | What it has | Response labels? | Usable? |
|---|---|---|---|
| GSE113342 (Ohio State, NanoString, 14 pts) | diagnosis vs repeat-flare biopsy, glom + TI | no (flare, not response) | no; same centre, likely overlapping patients |
| GSE127797 (Toronto, 54 biopsies; Almaani 2019) | glom + TI, ISN class only | no ("contact Dr Reich") | no (would need contacting someone: forbidden) |
| GSE239576 (Nagasaki, 21 pts, RNA-seq) | ISN class, "active LN" | no | no |
| GSE32591 / GSE104948 (ERCB) | glom + TI microarray | no | no |
| CELLxGENE | only blood lupus collections (Perez 2022, CD3xCD19) | - | no kidney |
| AMP Phase 1 kidney scRNA (Arazi 2019; Der 2019 tubular IFN vs 6-month response) | response labels exist | yes | ImmPort SDY997: requires a registered account |
| AMP Phase 2 kidney (Synapse syn64064827) | week-52 response | yes | controlled access (Synapse login + DUA) |

Verdict: BLOCKED. Every dataset with kidney expression AND response labels needs an account (ImmPort or Synapse).
A frozen test is ready to run the moment data arrive: score = mean z of the 10 round-91 NR-up genes (CD99, MAP4K4,
CSF1R, PTPN6, TGFB1, HLA-DRA, IFI16, LILRB1, LILRB2, C1S) minus mean z of the 5 CR-up genes (IL1RN, ENTPD1, RAG2,
CFP, CASP10), TI / tubular-interstitial cells, post-treatment where available else baseline. Pass: NR > CR,
one-sided Mann-Whitney p < 0.05, AND AUC beats baseline UPCR alone by >= 0.05. Must also beat 1,000 random
15-gene sets of the same panel (score AUC above the 95th percentile).
