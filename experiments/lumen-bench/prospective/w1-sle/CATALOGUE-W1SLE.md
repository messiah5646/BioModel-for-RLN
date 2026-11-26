# RT-20261126-1 run 1: catalogue of unreported non-renal SLE RCTs (round 146, 2026-11-26)

**Registry metadata only. No outcome numbers were read or recorded.** This is not a PREREG; the PREREG (PREREG-W1-SLE) is run 2, >= 1 h after this push (rule 12b), followed by the freeze as a separate later commit.

Decision this can change: whether the shipped pooled-prior headline in `lumen_forecast.py` stays (pooled beats null by >= 10% on >= 5 prospectively scored SLE trials) or is replaced by null-centred intervals and dropped from the preprint (pooled loses by >= 10%). For run 1 itself: if >= 5 trials are eligible -> run 2 writes the PREREG and freeze; if < 5 -> mark the route UNDERPOWERED-at-catalogue and park it.

## Source
- CT.gov API v2, `query.cond=systemic lupus erythematosus`, status COMPLETED | ACTIVE_NOT_RECRUITING, phase 2 or 3, allocation RANDOMIZED: 183 studies (`raw/ctgov_sle_p1.json`, fetched via `fetch.sh`, md5 688da1d77368e35eb11aa16f71dcaa94, `MANIFEST.tsv`).
- PubMed E-utilities search by NCT ID for each registry-eligible SLE-disease-activity trial (`raw/pubmed/`), titles only.
- Web search (2026-11-26) for sponsor press releases / toplines; only yes/no and the source were recorded.

## Registry rules (`catalogue_w1sle.py`, deterministic)
1 no results posted on CT.gov (87 dropped) · 2 not already tabulated or seen by LUMEN (LN 20 / SLE 40 tables, RT-3, RT-4: 18 dropped) · 2b not in the W1 LN freeze (4) · 3 non-renal (13) · 3b not cutaneous-only (3) · 4 randomised · 5 a placebo / SoC arm (21 dropped) · 6 enrolled >= 20 · 7 primary endpoint a responder proportion (17 dropped). **20 passed.**

## Manual screen (`topline_screen.csv`)
- Endpoint not SRI-x / BICLA disease activity (vaccine response, LLDAS, continuous MEX-SLEDAI): 4 excluded.
- **Topline public before this catalogue: 6 excluded** (N-acetylcysteine 2012 paper; TULIP-SC; GLPG3667; envudeucitinib LUMUS; deucravacitinib POETYK SLE-1 and SLE-2, BMS release 2026-11-25).
- **Eligible: 9**, plus 1 pending (sirolimus NCT04582136: a news headline on mTOR inhibition in lupus was not opened; run 2 resolves it, exclude if unresolved).

| NCT | trial | primary | primary completion (registry) |
|---|---|---|---|
| NCT02633163 | MSC (MiSLE) | SRI wk24 | 2025-08 |
| NCT05039840 | frexalimab (APATURA) | SRI-4 wk24 | 2026-09 |
| NCT06559163 | obexelimab (SunStone) | BICLA wk24 | 2026-09 |
| NCT05672576 | cenerimod OPUS-1 | SRI-4 wk52 | 2026-12 |
| NCT05648500 | cenerimod OPUS-2 | SRI-4 wk52 | 2027-03 |
| NCT05843643 | upadacitinib SELECT-SLE | BICLA wk52 | 2027-03 |
| NCT05639114 | ianalumab SIRIUS-SLE 1 | SRI-4 wk60 | 2027-03 |
| NCT05624749 | ianalumab SIRIUS-SLE 2 | SRI-4 wk60 | 2027-04 |
| NCT06015230 | GR1603 ph1b/2 | BICLA wk24 | 2028-06 |
| (pending) NCT04582136 | sirolimus | SRI-4 wk24 | 2025-04 |

Carry-forward slots (rule 8 re-opening path): only the two ianalumab trials register an earlier responder visit (SRI-4 wk36 -> wk60).

## Notes for run 2
- Re-screen every eligible trial for a topline that appeared between this push and the PREREG push; exclude any that did.
- Same-molecule pairs (cenerimod x2, ianalumab x2) should be named in the PREREG with a one-molecule-per-pair sensitivity (cf. round 139).
- Accrual is slow: on registry dates, >= 5 could report by mid-2027 at the earliest; MiSLE (completed 2025-08) may never report.
- Scope gap: RECRUITING trials were not catalogued (route spec); most large SLE phase 3s still recruiting will report later anyway.
