# D-BICLA check (round 156, 2026-11-27)

Decision this can change: whether bicla_table.csv is rebuilt. If the vobarilizumab-ph2 wk48 rows are a parser/registry-copy
artefact of our build -> rebuild via build_sle_table.py with a documented rule; if they match the registry as posted -> keep the
table (md5 99eb279e) unchanged.

Source: ClinicalTrials.gov API v2, NCT02437890, fields resultsSection.outcomeMeasuresModule, fetched 2026-11-27 (this run).
Outcome "Number and Percentage of Subjects With mBICLA Response at Week 24 and Week 48" (denominators 62/64/62/62/62):

| arm | wk24 registry | wk48 registry | bicla_table wk24 / wk48 |
|---|---|---|---|
| Placebo | 28 | 28 | 28 / 28 |
| 75 mg q4w | 28 | 32 | 28 / 32 |
| 150 mg q4w | 22 | 22 | 22 / 22 |
| 150 mg q2w | 24 | 19 | 24 / 19 |
| 225 mg q2w | 22 | 22 | 22 / 22 |

Result: all 10 rows match the posted registry exactly. The identical wk24/wk48 counts (3 of 5 arms) are in the sponsor's
posting, not introduced by build_sle_table.py. No peer-reviewed full report exists (PubMed "vobarilizumab systemic lupus
erythematosus": 0 hits), so whether the sponsor copied wk24 into wk48 cannot be checked on open data.
Side note: the separate primary-endpoint record (wk24 only, MCP-Mod analysis) gives 29/28/24/24/23, 0-2 responders different
from the wk24 class above (different analysis set/imputation); the table consistently uses the two-visit measure.

Decision: keep bicla_table.csv unchanged (md5 99eb279e). No rebuild; W1-SLE frozen BICLA forecasts unaffected.
Label: Established that the table matches the registry (n 1 trial, 10 arm-visit rows); Hypothetical that the registry wk48 rows
are a sponsor copy (no open source can test it).
