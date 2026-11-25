# Disclosure: RT-4 run 3 fetch time (round 141, 2026-11-25)

The queue note set run 3 "due >= 14:27Z" (1 h after AMENDMENT-A140-RT4, 03df5c2, pushed 13:27:20Z).
The first fetch.sh call of run 3 was at **14:24:38Z**: 57 min 18 s after that amendment (2.7 min short of the queue's
1-h target), **2 h 57 min after PREREG-RT4** (67f2d26, 11:27:17Z) and 1 h 59 min after AMENDMENT-R20261125-RT4
(414a21e, 12:25:39Z). Cause: the run misread the clock (assumed 14:27 had passed). Rule 12b's 1-h gap is defined from
the PREREG push and is met; both amendments were pushed before any fetch (fetch_utc > amendment commit times, see
sources/MANIFEST.tsv). The A140 amendment adds descriptive sensitivities only (S9, label text) and was not edited after
13:27:20Z. Label unchanged: "12b-timed (late outcomes already seen; early values new for 8 of 10 Table A trials)".
