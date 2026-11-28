# W3 run 1 (round 169, 2026-11-28): catalogue of unreported SLE/LN RCTs beyond the 21 W1/REC1 trials

**Registry metadata only. No outcome number was read or recorded** (studies with posted results are dropped before any other field is used; 20 of the 22 registry-passers are RECRUITING / NOT_YET_RECRUITING / ENROLLING_BY_INVITATION, so no outcome can exist; the other 2 are ACTIVE_NOT_RECRUITING with primary completion 2027-01 and COMPLETED-none).

Decision this can change (queue W3): whether any rule-8 settlement can ever reopen. If >= 10 eligible unreported trials beyond the 21 exist now -> PREREG-W3 + first batch freeze (freeze >= 1 h after the PREREG push, rule 12b); if < 10 -> record accrual and LUMEN states settled lines as 'closed for the foreseeable future'.

## Source
CT.gov API v2, `query.cond=lupus`, phase 2 or 3, allocation RANDOMIZED, two pulls via `scripts/provenance/fetch.sh` (`MANIFEST.tsv`):
open statuses (RECRUITING | NOT_YET_RECRUITING | ENROLLING_BY_INVITATION, 59 studies) and closed (ACTIVE_NOT_RECRUITING | COMPLETED, 200 studies). UNKNOWN / WITHDRAWN / TERMINATED not pulled.

## Registry rules (`catalogue_w3.py`, deterministic; = W1-SLE rules extended to renal trials and open statuses)
1 results posted (91 dropped) · 2 already tabulated/seen by LUMEN (18) · 2b in the 21-trial W1/REC1 freeze (19) · 2c W1-SLE manual exclusion 2026-11-26 (12) · 3 not SLE/LN (29) · 3b cutaneous-only (3) · 4 randomised interventional · 5 placebo/SoC arm (39 dropped) · 6 enrolled >= 20 · 7 responder primary (26 dropped). **22 passed.**

## Manual screen (`manual_screen.csv`)
5 excluded (multi-disease basket; SLE-ITP endpoint; undefined 'therapeutic benefit'; CLASI-70 cutaneous primary; belimumab-or-telitacicept arm). **Eligible: 17 (14 SLE: 11 SRI-4, 3 BICLA; 3 LN: 2 CRR, 1 KDIGO renal remission scored as CRR).**
Topline screen: none needed for the 20 open-status trials (no data can exist); the ACTIVE_NOT_RECRUITING dapagliflozin LN trial (NCT07678203) has primary completion 2027-01. No candidate reached primary completion before the analyst model's training cutoff (June 2026).

## Count vs the W3 gate
17 >= 10 -> **PREREG branch.**

## Accrual (descriptive, registry dates; Hypothetical as a projection)
- Eligible trials by first-posted year (21 frozen + 17 new; older years under-count because reported trials drop out): 2021 5, 2022 5, 2023 5, 2024 4, 2025 9, 2026 13 (to 28 Nov). Recent registration rate ~10/yr.
- Registry primary completion (upper bound on what could be scorable, before reporting lag and non-reporting): cumulative 9 by end 2026, 23 by end 2027, 36 by end 2028, 38 by end 2029 (21 + 17), plus ~10/yr new.
- Against POW1's pre-stated re-open n (round 167): O2 (n >= 20) earliest ~2028-29; CAL1 gap 0.20, O1, O3 (n >= 50) ~2030+; CAL1 gap 0.15 (n >= 100) mid-2030s. W1's 21 have 0/21 reported so far although 7 had registry primary completion by end 2026, so non-reporting/lag will push these later.
