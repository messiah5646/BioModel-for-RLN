# CURE1 source catalogue (round 163, 2026-11-28, scheduled run)

**Decision this can change:** whether CURE2 (cure-mixture vs long-remission PREREG) is feasible.
If >= 20 unique patients with open per-patient follow-up exist -> write CURE2 PREREG; if < 20 -> CURE2 blocked,
LUMEN prints only "remission duration observed, cure not assessable", track jumps to CURE4.

**Outcome blinding.** No per-patient relapse flag, relapse time, remission status or follow-up value was read or
extracted. What was looked at: PubMed titles, publication types, regex snippets of abstracts giving only "<n> patients"
and "follow-up ... <duration>" (search_meta.json), and for PMC open-access papers the table captions, header cells and
first-column row labels (no body values). Titles carry aggregate outcome words ("durable remission", "relapse of SLE"),
so CURE3 on these sources is "12b-timed", not "verifiably pre-registered" for aggregate direction (same convention as A135).

## Searches (reproducible)
- PubMed esearch, 2026-11-28: `("chimeric antigen receptor"[tiab] OR "CAR-T"[tiab] OR "CAR T"[tiab] OR "CAR T-cell"[tiab] OR blinatumomab[tiab] OR teclistamab[tiab] OR mosunetuzumab[tiab] OR "T cell engager"[tiab] OR "T-cell engager"[tiab] OR bispecific[tiab]) AND (lupus[tiab] OR SLE[tiab])` -> 393 hits (pubmed_hits_393.txt); 142 after dropping reviews/editorials/preclinical titles; 63 human-treatment candidates hand-picked by title (candidate_pmids.txt).
- ClinicalTrials.gov API v2, cond=lupus, intr=CAR-T / chimeric antigen receptor / blinatumomab / teclistamab / mosunetuzumab / T cell engager / CAR-NK -> 173 registered studies, **0 with posted results** (ctgov_reset_lupus.tsv).

## A. Sources with per-patient follow-up confirmed open (header/row-label level)
| id | therapy | cohort | n | max follow-up (abstract) | per-patient last-visit follow-up printed | per-patient remission status at last visit | per-patient time-to-relapse printed | drug-free status printed |
|---|---|---|---|---|---|---|---|---|
| PMID 40906188 / PMC12488764 (Z Rheumatol 2025, Erlangen long-term) | autologous CD19 CAR-T | Erlangen SLE (master record; subsumes PMID 34347960, 36109639, 38381673 SLE arm) | 11 | ~2.5 y (abstract) | yes ("Follow-up (Jahre)" row) | yes ("SLEDAI (letzte Visite)", "DORIS-Remission (letzte Visite)") | no column | not as a row (DORIS allows low-dose prednisolone/HCQ); text to check in CURE3 |
| PMID 41612420 / PMC12857108 (J Hematol Oncol 2026, extended f/u) | BCMA-CD19 compound CAR-T | Chinese cCAR phase 1 (master record; subsumes PMID 38777376, 41787704, 39690781, 41669912, likely 40549291) | 12 | not in abstract snippet (38777376 abstract: up to 46 months) | yes ("Last follow UP (months)") | yes ("Remission status") | no column | not as a column; text to check |

**Unique patients with open per-patient follow-up: 23** (11 + 12; different countries and products, no overlap expected).
Both give *current-status* data (status at a per-patient last visit), not event times.

## B. Open (PMC) sources, per-patient follow-up NOT confirmed at header level
| id | therapy | n | note |
|---|---|---|---|
| PMID 41501497 / PMC13004673 (Nat Med 2026, CASTLE basket) | CD19 CAR-T (Erlangen) | 24 (10 SLE) | tables aggregate; "Source Data" file present; overlap with 40906188 likely -> counted 0 |
| PMID 40386685 / PMC12083988 (relma-cel phase 1) | CD19 CAR-T | 8 (12 enrolled per snippet) | per-patient baseline table only; outcomes likely in figures/supplement |
| PMID 41541778 / PMC12799546 | CD19/CD22 bicistronic CAR-T, paediatric LN | 3 | baseline per-patient table only |
| PMID 40341741 / PMC12297247 | allogeneic CD19 CAR-T | ? | letter + supplement |
| PMID 42309556 / PMC13289224 | subcutaneous T-cell engager, phase Ib | 15 | aggregate tables |
| case reports PMID 41266084, 41920311, 42621160 (2), 40918126 | CD19 CAR-T | 1-2 each | narrative only |
| PMID 42876831 | CD19 CAR-T | 1 | **selected because of relapse** -> exclude from estimation (outcome-selected) |

## C. Not open full text (publisher only; supplements may be free but not checked)
NEJM 34347960, 38381673, 39047248, 39231352, 40961420; Lancet 38642568, 39615991, 41240964 (CAR-NK); Nat Med 36109639, 40866583, 40993243, 41644747, 42722888; ARD 38777376, 40450406, 40681435, 42321107, 42414128, 42680668; Med 40446794, 41856114, 42486100; A&R 41787704, 42098533; others in candidate_pmids.txt. Several overlap groups A.

## D. Historical comparators (conventional therapy), structure only
| id | design | per-patient? |
|---|---|---|
| PMID 35725295 WIN-Lupus (ARD 2022) | RCT, maintenance IS withdrawal vs continuation in LN | aggregate (KM) |
| PMID 39832909 (LSM 2025) | IS withdrawal after complete renal response | aggregate |
| PMID 39037915, 43015040, 38673484 | IS discontinuation cohorts in LN | aggregate |
| PMID 32870418, 25086082, 16078321 | SLE long-term / prolonged remission without treatment | aggregate |
| PMID 36117428, 11263789 | autologous HSCT for refractory LN / childhood SLE (earlier "reset" therapy, long follow-up) | to check |

## Verdict for the gate
23 unique patients >= 20 -> **CURE2 PREREG proceeds.** Constraint CURE2 must handle: the data are current-status
(in remission or not at one last visit per patient, follow-up mostly <= ~4 years), no relapse times; the PREREG must
use a current-status likelihood (P(relapse-free at t_i)) and say how re-remission after relapse is coded. Time-forward
split available: Erlangen 2025 (n 11) earlier, cCAR 2026 (n 12) later = held-out n 12 (>= 10 for Supported).
Expected (Hypothetical, before any data): with n 23 and <= 4 years, cure vs decreasing-hazard Weibull is unlikely to be
separable, so INCONCLUSIVE -> CURE4 is the probable decision-maker.
