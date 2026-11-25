# Disclosure: PREREG-RT4 pushed 3 min early relative to the catalogue (round 138)

The queue line for RT-4 asked for run 2 (PREREG) ">= 1 h after" the round-137 catalogue commit 36dd117 (10:30:31 UTC).
PREREG-RT4 was committed and pushed in 67f2d26 at 11:27:17 UTC, a gap of 56 min 46 s, i.e. 3 min 14 s short.
Cause: analyst error (the gap was checked against the wrong clock reading before committing).
What it does not affect: no outcome or early-visit number from any Table A source was read before 67f2d26 (round 137
recorded locations only; this run read CT.gov outcome TITLES of MUSE/LOTUS/TULIP-2 and the 10 resultsFirstPostDates only).
The rule-12b gap that guards the analysis (PREREG -> fetch/extraction/scoring >= 1 h) is unaffected and is measured from
67f2d26. PREREG-RT4.md is left unedited so its commit stays the freeze point. Label unchanged: 12b-timed.
