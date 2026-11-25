# Round 137b open-data harvest (live session, 2026-11-25; Thomas: "do more, get all data")

Every file in `harvest_list.tsv` was fetched with `scripts/provenance/fetch.sh` into `data/real/harvest-r137/`
(gitignored, about 17 GB). Fetch times, sizes, md5 and HTTP Last-Modified are in `MANIFEST-harvest.tsv`; failures
in `harvest_failures.tsv`. Re-run: `bash experiments/unblock-2026-11-25/harvest.sh` (skips files already in the manifest).

## What was fetched and why
| Group | Files | Label fields (from round 137 metadata checks) | Queue status |
|---|---|---|---|
| Blood with per-donor nephritis labels | GSE135779 (matrix, genes, RAW 1.3 GB) + Nehar-Belaid 2020 MOESM3 (ST1b clinical); GSE226147 platelet counts; GSE211700 FPKM | SLEDAI, renal items, KIDNEY, Neph_all, class (GSE135779); per-item renal SLEDAI, progression status (GSE226147); LN vs non-LN (GSE211700) | Q2/Q14 question settled (rule 8, r99): stays closed |
| Serial blood in LN under treatment | GSE72747 / GSE72798 (matrix + RAW) | 10 LN, before / 3 / 6 months, SLEDAI, GFR | no consumer (IFN layer parked, r123) |
| Urine at biopsy with response | PLoS One 2022 S1 (pone.0263778.s002), PLoS One 2025 S1 (pone.0323334.s001) | MCP-1, EGF, UPCR, AI/CI, 6-month CR/PR (2022); cross-sectional (2025) | Q62 parked (no consumer) |
| Kidney spatial / single-nucleus | Danaher 2024 CosMx (analyses zip + cleaneddata.RData); VCAM1+ PT 2025 snRNA/snATAC/Visium (Zenodo 15868389) | class, glomerular lesion scores (Danaher); WHO class, SLEDAI (VCAM1) | descriptive; no response/chronicity labels |
| CAR-T / bispecific blood | STAR-T (Zenodo 15715577: proteomics, CyTOF, Visium HD, scRNA); CD3xCD19 bispecific CELLxGENE collection 6c147b03 (5 h5ad) | timepoints, dose | IFN-consumer rule applies (r123) |

## Deliberately NOT fetched here
- **RT-4 sources** (FDA/EMA reviews, PMC full texts in CATALOGUE.md Table A). PREREG-RT4 (67f2d26) requires run 3 to fetch
  them after the 1-h gap; fetching them in this harvest would break that.
- Anything behind a login, click-through, captcha/bot check or paywall (GSE200306 paper supplements, Danaher Table S7,
  AMP/ImmPort/Synapse, dbGaP). The policy forbids it and round 137 lists them for Thomas.

## Consequence for future tests
All files above were downloaded before any pre-registration. Any later test on them can be at best "12b-timed
(data committed/fetched before the PREREG; earlier looks not excludable)". Round 137 inspected metadata and label levels
only; no outcome relationship was analysed.
