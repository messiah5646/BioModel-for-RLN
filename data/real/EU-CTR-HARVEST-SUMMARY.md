# EU Clinical Trials Register — lupus nephritis per-visit endpoint harvest

Harvested 2026-11-02. Source register: https://www.clinicaltrialsregister.eu
All results pages are at `/ctr-search/trial/<EUDRACT>/results`.

## Enumeration

Search `?query=lupus+nephritis` returned **64 trials across 4 pages** (20/20/20/4).
Counting the "View results" markers across all four pages, **~40 trials carry a results
section**, not 16 as a prior session recorded. The 16 figure matches page 1 alone.
One trial (2005-004575-37, abatacept IM101075) is marked **"Removed from public view"**.

Note: not all 64 hits are lupus-nephritis trials. Several are SLE trials where nephritis is
a subpopulation or an exclusion criterion, and three page-4 hits (efalizumab psoriasis
27577, BeCarma cardiac MRI, EXTRAStudy ACE-inhibitor histopathology) are not LN trials at all.

## Trials checked (22 results pages opened)

| EudraCT | Trial / drug | Results? | Per-visit? | Visits | Endpoints with >1 post-baseline visit | Verdict for longitudinal fitting |
|---|---|---|---|---|---|---|
| 2017-003230-93 | CFZ533 / iscalimab, CCFZ533X2202 | yes | **yes** | **6** | UPCR ratio-to-baseline (Days 197,225,253,281,309,337); urine hyaline casts (~15 visits) | **HIGHEST VALUE new addition.** Only trial in the register with a renal endpoint at 6 post-baseline visits with mean+SD+per-visit n. Harvested. |
| 2015-002022-39 | Obinutuzumab, WA29748 (NOBILITY) | yes | **yes** | **15** | CD19+ B cells (baseline, wk2/4/12/24/52/104 + Bcfu months 6-48); PtGA VAS (6 visits) | **High value for a B-cell/PD submodel**, worthless for renal. Unit/title contradiction and duplicated rows. Harvested. |
| 2017-003101-17 | BI 655064 maintenance, 1293-0013 | yes | **yes** | 4 | SLEDAI total change (wk 12,26,42,52) | Moderate. Continuer/responder subset, no baseline anchor, n rises wk42→wk52. Harvested with warnings. |
| 2012-003364-51 | Voclosporin, AUR-VCS-2012-01 (AURA-LV) | yes | partial | 2 | UPCR change + SELENA-SLEDAI change (wk 24,48), with absolute values too | Moderate. Only 2 timepoints, but has both absolute and change values and a 3-arm dose contrast. Harvested. |
| 2012-000714-11 | Abatacept, IM101-291 | yes | partial | 2 (1 usable) | UPCR adjusted mean (D365,D729); UPCR median %; eGFR adjusted mean (D365,D729) | Low. Day 729 is all 99999 placeholders, so effectively one timepoint. Mixed SE/CI/IQR. Harvested as a trap record. |
| 2019-003211-57 | Secukinumab, CAIN457Q12301 (SELUNE) | yes | yes (PRO only) | 4 | FACIT-Fatigue, SF-36v2 PCS, LupusQoL (wk 12,24,36,52) | Low for renal: UPCR at wk52 only. Harvested for the PRO series and the completion-collapse trap. |
| 2016-004045-81 | Voclosporin AURORA 1 | yes | yes | — | — | Already harvested in a prior session. Not re-opened. |
| 2016-004046-28 | Voclosporin AURORA 2 | yes | yes | — | — | Already harvested in a prior session. Not re-opened. |
| 2011-004570-28 | Belimumab, BEL114054 (BLISS-LN) | yes | **no** | 1 | none | **Low.** Response rates only (PERR/CRR/ORR at wk 52 and 104). No continuous endpoint at any visit. Surprising given the trial's size (224+224). |
| 2015-001442-29 | Anifrolumab, D3461C00007 (TULIP-LN) | yes | **no** | 1 | none | Low. Primary is a Week-52 geometric mean ratio of UPCR only; the model includes visit terms but no per-visit values are printed. |
| 2015-001750-15 | BI 655064 induction, 1293.10 | yes | no | 2 (responder %) | CRR/PRR/MRR at wk 26 and 52 — proportions, not change scores | Low. No continuous endpoint reported. |
| 2018-004142-42 | BMS-986165, IM011-073 | yes | no | 1 | none | **Worthless.** Terminated for low enrolment; only 1 participant in the treated arm; no statistics performed. |
| 2018-003155-38 | Guselkumab, CNTO1959LUN2001 | yes | yes (labs only) | 2 | ~20 haematology/chemistry labs at Baseline, wk24, wk52 | Low. n = 0–4 per arm at wk52; 9999/99999 placeholders. No UPCR/eGFR/C3/C4/dsDNA numeric series. |
| 2020-005568-79 | Nipocalimab, 80202135LUN2001 | nominal | no | — | none | **No usable data.** Page carries only a summary-report attachment; no arms, no endpoints. Attachment label reads "EudraCT Cancelled Trial". |
| 2004-004917-41 | MMF vs IVC, WX17801 (ALMS induction) | yes | no | 1 | none | Low. Three Week-24 responder endpoints; comparisons use "mean difference (final values)", not change from baseline. No ALMS maintenance phase reported here. |
| 2006-005357-29 | Ocrelizumab, WA20500 (BELONG) | yes | nominal | many | SLEDAI-2K reduction, FACIT, SF-36, mBPI-SF, AUC of cGFR listed across wk 1–48 | **Worthless in practice.** Terminated early; most secondary endpoints show "not applicable" with 0 participants analysed. |
| 2010-020968-38 | CNTO 136 / sirukumab, CNTO136LUN2001 | yes | no | 1 | none | **Worthless.** Placebo arm n=4. All change endpoints at Week 24 only. |
| 2011-002159-32 | BIIB023, 211LE201 | yes | no | 1 | none | Low. All endpoints are Week-52 responder percentages. No continuous series. |
| 2016-003201-34 | Rituximab, U2970g | nominal | no | — | none | **No usable data.** Summary-report attachment only; no arms, no endpoints on the page. |
| 2010-018329-20 | Laquinimod, LN-LAQ-201 | nominal | no | — | none | **No usable data.** Summary-report attachment only; no arms, endpoints or sponsor field. |
| 2007-004518-15 | Atacicept + MMF, 28113 | nominal | no | — | none | **No usable data.** Summary-report attachment only. |
| 2021-005772-19 | Secukinumab OLE, CAIN457Q12301E1 | yes | no | 2 | CRR at Week 104E1 and Week 132 (proportions) | **Worthless and dangerous.** n=9/9 then 6/5; 0 of 31 completed; eGFR criterion anchored to the PARENT trial's "core Baseline". |
| 2022-001377-31 | Daxdilimab, HZNP-DAX-203 | yes | no | 1 | none | **Worthless.** Terminated; 0 of 19 completed; eGFR change at wk52 has "no data collected". |
| 2020-004047-86 | ALPN-101 / acazicolcept, AIS-A03 | yes | no | 1 | none | Low. SLE not LN. SLEDAI-2K change at Day 169 only. |
| 2019-000638-20 | Elsubrutinib/upadacitinib, M19-130 | yes | no | 1 | none | Low. SLE not LN. All endpoints at Week 24. |
| 2017-001489-53 | Ustekinumab, CNTO1275SLE3001 | yes | no | 1 | none | **Worthless.** SLE not LN; exactly one endpoint (SRI-4 at wk52) on the whole page. |
| 2019-000328-16 | AMG 570 / rozibafusp alfa, 20170588 | yes | yes (PRO only) | 5 | PROMIS-Fatigue, SF-36v2 PCS, LupusQoL, PtGA (wk 12,24,36,44,52) | Low for renal. SLE not LN; no renal or serology continuous endpoint. Not harvested to CSV. |
| 2018-000305-23 | BOS161721, BOS161721-02 | yes | yes (biomarker) | 4–9 | anti-dsDNA change (Days 30–270); pSTAT3 (D30,44,60,90); IL-21 signature (D15,90,180,270, all 99.999 placeholders) | Low–moderate. SLE not LN; Phase 1b n=5–9 per cohort. Worth a second look only for a serology submodel. |

## Trials with results NOT opened (triage reasons)

Register marks results available but not fetched this session, in rough priority order
for any follow-up: 2013-000594-69 (BIIB023 long-term extension 211LE202 — likely a
continuer subset of 211LE201), 2011-005856-32 (RITULUP, two rituximab regimens),
2008-002221-36 (atacicept + CYC 28779), 2006-002107-13 (CERL080A2420 EC-MPS, 6-month
open-label), 2005-003070-19 (ERL080A2410 EC-MPS), 2012-004893-25 (RITUXILUP CRO2035),
2015-005543-14 (belimumab after B-cell depletion CTU/2013/096), 2014-000488-42 (SYNBioSe),
2017-001619-35 (P1200_41 rituximab pathophysiology), 2007-001377-28 / 2021-001802-30 /
2015-003884-12 (not LN trials).

## Register-wide observations

1. **The register's dominant reporting style is the single-timepoint responder percentage.**
   Of 22 pages opened, only 6 carry any continuous endpoint at more than two post-baseline
   visits, and only 1 (CFZ533) carries a *renal* one.
2. **Four pages carry results in name only** (nipocalimab, rituximab U2970g, laquinimod,
   atacicept 28113): a PDF summary-report attachment with no arms, endpoints or tables in
   the HTML. The register's "results available" flag does not imply extractable data.
3. **Dispersion convention is not consistent, even inside one trial.** Abatacept IM101-291
   prints SE for UPCR, a 95% CI for eGFR and an IQR for median UPCR on the same page.
   CFZ533 prints a geometric mean + 95% CI for its primary UPCR and an arithmetic mean + SD
   for the same measure as a secondary.
4. **99999 / 9999999 / 99999999 / 9999 / 99.999 are all missing-data placeholders** in
   different trials. Only obinutuzumab prints a footnote explaining its code. Any automated
   numeric scrape of this register will silently ingest these as values.
5. **Baseline anchors are frequently absent.** CFZ533, abatacept, BI 655064 maintenance and
   SELUNE all print change or ratio scores with no absolute baseline for that endpoint,
   making cross-trial rescaling impossible.
