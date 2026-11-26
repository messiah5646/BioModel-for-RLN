# Disclosure: run-2 timing (round 147)

The route asked run 2 to start >= 1 h after run 1. The run-1 catalogue commit 4ed6f82 is timed 16:27:28 UTC; the
PREREG-W1-SLE commit 5515355 is timed 17:27:07 UTC: **59 min 39 s, 21 s short of 1 h.**

Impact: none on blinding. The 1-h rule (12b) protects the gap between a PREREG and outcome data; here no outcome data
exist for any of the 8 trials (CT.gov hasResults = false for all, re-checked 2026-11-26 ~17:25 UTC; no topline found),
and scoring happens only in later W1 runs when results post, months after this freeze. The catalogue was not edited by
run 2 (md5 in the PREREG matches 4ed6f82). Not corrected by re-committing, to keep the record as it happened.
Reviewer may treat the catalogue -> PREREG step as "12b-timed minus 21 s".
