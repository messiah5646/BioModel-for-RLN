# Q9 catalogue: public longitudinal urine data for an LN flare early-warning test (2026-11-22)

Question (round 89 follow-up): does urine CD163 / CD206 (MRC1) rise BEFORE a lupus nephritis flare, where blood IFN did not?
Needed: patient-level, serial urine measurements with flare dates (or at least pre-flare and flare samples), openly downloadable.

## Searched
- NCBI GEO (gds, "lupus AND (urine OR urinary)", 60 series): no urine protein or urine-cell series with flare labels; hits are
  mouse kidney, cell lines, one urinary exosome miRNA set in pregnancy (GSE298542, n=11, no flares).
- Zenodo, Dryad (keyword searches): nothing relevant.
- PubMed (44 hits for urinary CD163/CD206 AND lupus) plus PMC data-availability statements of the longitudinal ones.

## What exists
| Source | Design | Flare labels | Patient-level data public? |
|---|---|---|---|
| Mejia-Vilet 2020 JASN (doi 10.1681/ASN.2019121285) | Mexican + OSU cohorts; uCD163 6 mo pre-flare -> flare -> 12 mo | Yes | No (figures only) |
| Renaudineau 2025 Rheumatology (doi 10.1093/rheumatology/keaf446) | 129 LN followed; uCD163 discriminated flares (>520 ng/mmol) | Yes | No |
| Gamez/Toulouse 2024 J Transl Autoimmun (doi 10.1016/j.jtauto.2024.100244) | diagnosis + follow-up | partly | "on request" |
| Fava 2024 JCI Insight, AMP Phase II urine Kiloplex (doi 10.1172/jci.insight.172569) | ~1200 urine proteins at biopsy, 3, 6, 12 mo, with response | No (response, not flare) | Synapse syn53124971: needs free Synapse login + ARK Portal clickwrap (checked via API: "unmet access requirements") |
| Anifrolumab TULIP-LN urine biomarkers 2026 A&R (doi 10.1002/art.70089) | trial, serial urine | No | AstraZeneca request only |
| Smith 2023 Rheumatology childhood LN (doi 10.1093/rheumatology/keac465) | serial | partly | not posted |

The JCI Insight "Supporting data values" file (downloaded from Europe PMC) holds only per-protein summary statistics
(AUC, p, logFC per contrast), not per-patient values.

## Verdict
No openly downloadable, patient-level, longitudinal urine dataset with flare labels exists that this project can use without
an account or a data request. The flare early-warning test cannot be run. Status: blocked (account needed).

The closest usable resource is AMP Phase II urine proteomics on Synapse (free account + click-through). It has serial urine
but labels treatment response, not flare, so it answers a neighbouring question (does urine CD163 at month 3 add to UPCR for
12-month response?). That test is pre-registered now in PREREG-AMP-urine.md so it is frozen before anyone sees the data.

Literature state (aggregate only): Supported that urine sCD163 tracks histological activity and is higher at flare than in
remission; Emerging that it rises before a flare (one cohort, Mejia-Vilet 2020, 6 months pre-flare -> flare).
