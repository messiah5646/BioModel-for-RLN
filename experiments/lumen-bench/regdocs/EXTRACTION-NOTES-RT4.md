# RT-4 run 3 extraction notes (round 141, 2026-11-25) - committed BEFORE scoring

Judgement calls (stated before any score was computed):
1. **ACCESS wk52** = sum of printed strata (Table 4): only wk24 CR/PR patients continued; wk24 non-responders are
   non-responders at wk52. Population check: Table 4 wk24 CR counts 22 (aba) / 21 (ctrl) equal the existing wk24 row. Included.
   Design note: abatacept wk24 complete responders stopped abatacept + AZA at wk28 (disclosed; the row is kept as printed).
2. **LOTUS**: methods define SRI-4 as SLEDAI-2K >= 4 reduction, no new BILAG A / <= 1 new B, PGA < 10% worsening, without
   treatment failure; no OCS condition -> standard (PGA written as '<10%' rather than '<0.3 points'; judged equivalent).
   Fig 2B (wk24) uses the full FAS (n 293/193) while Fig 2A (wk52) uses the mFAS (173/116); the wk24 mFAS counts are
   derived by subtracting the printed post-hoc non-mFAS counts (55/120, 26/77). They equal the CT.gov wk24 counts on 173/116
   (A140 population check: difference 0 in both arms). Included. LOTUS early value was in committed CT.gov data (A140): S9 drops it.
3. **MUSE**: the FDA review prints SRI-4 without OCS condition at BOTH visits only under the 'Treatment Policy' estimand
   (ignores excess-OCS/immunosuppressant intercurrent events, missing = non-responder). Wk52 also has 'SRI-4 without OCS
   tapering' but wk24 does not, so the same-definition pair is Treatment Policy at both visits. PREREG section 4 requires
   only 'plain SRI-4 (no OCS-taper condition)' at both visits -> included literally; flagged because Treatment Policy is
   not the standard treatment-failure SRI-4 used in other rows. S9 drops MUSE.
4. FDA accessdata blocks the default curl user agent; the review was fetched with fetch.sh unchanged, with a browser
   user-agent supplied through CURL_HOME/.curlrc.
5. Figure-label values: LOTUS only (S7 drops it).

Before scoring, by count alone: trials with an early value newly extracted in RT-4 = LOTUS, MUSE (+ ACCESS, whose RT-4
row is the late wk52 visit that turns its existing wk24 row into an early visit) = 3 < 5 required (PREREG section 8).
The primary is therefore UNDERPOWERED by rule whatever the scores; scoring is still run and reported as descriptive.
