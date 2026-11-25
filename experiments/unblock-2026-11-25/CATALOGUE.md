# Round 137: re-search for open alternatives to every blocked item (live session, 2026-11-25)

Thomas asked: "What's blocked, try to fix that, and try alternatives." Open-data policy unchanged (no accounts, no
clickwraps, no requests, no contact). Only metadata, column names, label levels and table locations were inspected.
**No outcome numbers were extracted or analysed**, so every route below can still be pre-registered (rule 12b).
Downloaded files stay out of the repo on purpose: a later run must re-fetch them with scripts/provenance/fetch.sh
after its PREREG.

## Verdict by blocked item

| Item | Was blocked on | Best open alternative found 2026-11-25 | Unblocked? | Decision at stake? |
|---|---|---|---|---|
| RT-3 / Q61 early read-out rule (UNDERPOWERED r136, 7 trials) | registries post only the late endpoint | **FDA multidiscipline reviews + EMA EPARs + open PMC full texts** give by-arm early (wk 12-26) responder data for trials already in the 60-trial tables but missing an early visit there | **Partly** (see table A) | **Yes**: the CLI `--early` rule (pooled vs carry-forward vs null) |
| Prospective scoring (13 frozen W1 trials) | no results posted | re-checked CT.gov API 2026-11-25: 0/13 have results; LNP023 (NCT05268289) now TERMINATED | No | (waits) |
| Q12 kidney response labels | AMP/ImmPort accounts | none qualifying; Danaher 2024 CosMx (Zenodo 13964258) has 8 pts, 7 responders / 1 refractory | No | - |
| Q19 / Q34 per-biopsy chronicity + mRNA | GSE200306 CI only as medians | GSE200306 paper (PMC9613357) supplement-2.xlsx may hold per-biopsy AI/CI; Danaher Table S7 (n ~10 biopsies) | **Needs a human browser** (bot challenge / 403) | weak: kidney IFN has no LUMEN consumer (r123) |
| Q9 urine flare | no open serial urine + flare data | none serial. Closest: PLoS One 2022 (PMC8912200) S1 file, 60 LN with outcome, urine MCP-1/EGF/UPCR at biopsy + 6-mo CR/PR + AI/CI per patient | No (not serial) | none today (no patient-level LUMEN layer) |
| Q16 / Q51b damage floor | no response split by chronicity | same PLoS One 2022 file: per-patient CI + 6-mo response + follow-up UPCR (one visit; cannot separate floor from rate) | Partly | none today (damage floor Hypothetical, no consumer; r126 gate) |
| Q2 Perez nephritis labels (dbGaP) | controlled | **GSE135779** (Nehar-Belaid 2020 scRNA, ~40 SLE) + Nat Immunol MOESM3 sheet ST1b: per-donor SLEDAI, renal items, `KIDNEY`, `Neph_all`, class; GSE226147 platelets (51 SLE, per-item renal SLEDAI) | **Yes, data are open** | **No**: blood IFN vs nephritis is settled (rule 8, r99); possible overlap with GSE65391 unverified |

## Table A: early read-out sources (locations only; no numbers recorded)

Trials already in the RT-1 60-trial set but with NO early visit there, for which an open source gives by-arm early data:

| Trial (NCT) | In 60 set, early there? | Open source | Location | Tabulated or figure |
|---|---|---|---|---|
| ACCESS LN (NCT00774852) | yes, target wk24 only | PMC4528976 | Table 2 (wk24), Table 4 (wk52 split by wk24 status) | tabulated (wk52 needs care) |
| BLISS-LN (NCT01639339) | yes, wk104 only | EMA Benlysta II/0080 variation report | Fig 19 p57/p59 (PERR/CRR by visit, %±SE, N per arm); Table 26 p55 | early = figure |
| BLISS-SC (NCT01484496) | yes, wk52 only | EMA Benlysta X/0046/G; PMC5434872 Fig 2A | Fig 7 p37 | figure |
| BLISS-NEA (NCT01345253) | yes, wk52 only | PMC5867402 | Fig 2A | figure |
| BLISS-76 (NCT00410384) | yes, wk52 only | PMC5007058 | Fig 2A | figure |
| PLUTO (NCT01649765) | yes, wk52 only | PMC7509523 | Fig 2A | figure |
| TULIP-2 (NCT02446899) | in SLE pipeline map | EMA Saphnelo EPAR | Fig 8 p87 / Fig 9 p88 (BICLA by visit, % with n) | figure |
| MUSE (NCT01438489) | in build_sle_table EXCLUDED set (definition rule) | FDA Saphnelo review 761123 | Tables 20-24 pp86-88 (SRI-4 wk24 with OCS taper; wk52) | tabulated; wk24 endpoint includes OCS taper -> definition check |
| LOTUS (NCT03517722) | EXCLUDED set (definition rule) | PMC9606504 | Fig 2A/2B (% to 1 dp with n) | figure, near-exact |
| telitacicept ph2b (NCT02885610) | RT-3 N06 (late only) | PMC10958275 | Fig 2B | figure |

Already early in RT-1 (sources would only cross-check): AURORA-1 / AURA-LV (FDA 213716 Tables 30/36/37/46/48;
voclosporin, excluded from the RT-1 primary), NOBILITY (PMC8762029 Fig 2), TULIP-1, CALIBRATE (PMC7839443 Table 3).
No early data found: REGENCY (EMA Gazyvaro report has wk50+ only), PAISLEY (wk32/48), LILAC, ianalumab.

### What this route can and cannot do
- **Can:** raise the early-read-out subset from 7 trials (RT-1 primary) to roughly 12-16, if figure-printed values are
  admitted under the RT-3 rule ("percentages with arm N converted to counts; values read off a plot NOT used unless
  the figure prints the number").
- **Cannot:** be "verifiably pre-registered". Every one of these trials' late outcomes is already in the analysed
  tables (r86, r112, r116, RT-1), so the best label is **12b-timed (late outcomes already seen; early values new)**.
  Rule 8 allows it only because the early visits are a new data type for these trials.
- Figure-only trials need a frozen rule for reading printed labels; trials whose figure shows bars without printed
  numbers are excluded, not digitised.

## Search records
Searches were run by three sub-agents (kidney; urine + blood-renal; regulatory) plus a CT.gov API re-check.
Kidney: ~25 GEO series, Zenodo/figshare, CELLxGENE (408 collections; lupus = blood only), ArrayExpress, PRIDE/CNGB.
Urine/blood: GEO "lupus AND urine" (200), 378 GEO lupus blood series characteristics parsed, PMC supplements.
Regulatory: Drugs@FDA reviews (213716, 761123, 125370), 6 EMA EPAR/variation reports, 30 PMC full texts.
Access failures: FDA .cfm index pages timed out/503; Benlysta-LN and Gazyva-LN FDA reviews not located;
science.org 403 and PMC bot challenge for two supplements (GSE200306 paper, Danaher Table S7); EMA 429 until spaced.
