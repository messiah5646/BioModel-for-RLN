# PREREG Q17 (round 107): scripted rebuild of resp_table.csv (2026-11-22)

## Why
resp_table.csv (CRR counts per arm per visit, 20 add-on LN RCTs) was hand-extracted in round 86, never committed
(*.csv git-ignored) and lost. Round 94 re-entered the 29 primary rows at 2 dp from result86.txt. Counts and the
early visits were gone, so Task C' and lumen_forecast.py's early-update (drift_var) could not be rerun.

## What is built
`build_resp_table.py` reads data/real/ctgov/lupus_with_results.json.gz (committed registry snapshot) and writes
resp_table.csv with columns trial, nct, arm, molecule, target, cls, week, r_drug, n_drug, r_ctrl, n_ctrl, source.
Selection is written as code before any comparison: 20 trials by NCT id (PREREG-r86 list); outcome measures whose
title names complete renal response (CRR / complete renal remission / complete response) as protocol-defined,
excluding partial/overall/modified/sensitivity variants unless a trial has no other CRR measure; control = the
placebo arm (or standard-of-care arm); counts = COUNT_OF_PARTICIPANTS, or PERCENT x arm N rounded when only %.
Visit week parsed from the time frame. Mapping rules may need trial-specific entries (arm names); these are
written in the script with a comment, not tuned against result86 numbers.

## Pass rule (fixed now, before the rebuilt table is compared)
- P1: all 29 primary arms of result86.txt reproduced: same (trial, arm, week); |lor diff| <= 0.01 and
  |se diff| <= 0.01 (2-dp rounding tolerance).
- P2: all 13 Task C' rows reproduced (te, tl, de, see, dl, sel within 0.01).
- P3: run86.py on the rebuilt table gives D0-D3 IV-MSE (primary and secondary) and C0-C2 IV-MSE within 0.002
  of result86.txt.
- P4: the 29 registry rows of prospective/w1/arm_table_reconstructed.csv match (follows from P1). The 3 week-104
  rows (NOBILITY, TULIP-LN1 x2; round 90) are from publications, not the registry; they are appended as
  source=publication rows with counts if recoverable from r90 inputs, else lor/se only, and are not part of P1-P3.
PASS = P1 and P2 and P3. Any mismatch is listed row by row; a mismatch is not fixed by editing numbers, only by
a documented selection-rule correction, and every such correction is reported.
