# P3 data step (round 115): table built after PREREG commit 64f5bd5; no scoring yet

`build_sle_table.py` on the committed CT.gov snapshot -> sri4_table.csv (119 rows, 80 arms, 40 trials) and
bicla_table.csv (83 rows, 53 arms, 26 trials). Onset class: fast 10, slow 17, other 13 trials.
Molecules with >1 trial (anchors): belimumab 5, baricitinib 3, anifrolumab 2, tabalumab 2, rigerimod 2.

## Parser fixes (allowed by the PREREG; before -> after, no previously parsed row changed)
1. "PBO" abbreviation recognised as placebo (PHOENYCS-GO, iberdomide).
2. "6 months" / "85 days" word order parsed (ACCESS abatacept-SLE).
3. A measure reporting SRI(4), SRI(5), SRI(6) as classes keeps the SRI(4) class (lulizumab).
4. A "0mg" placebo label is not a dose (lulizumab).
SRI-4: 106 rows / 36 trials -> 119 / 40. BICLA: 71 / 23 -> 83 / 26. Changed previously parsed rows: 0 and 0.
(An intermediate edit by the analyst broke the placebo rule for one lenabasum arm; caught by the before/after check
and fixed before this commit.)

## Excluded by the pre-registered rules (reported, not changed)
MUSE (NCT01438489) and daxdilimab (NCT04925934): responses conditioned on steroid taper. EMBRACE (NCT01632241): modified
SLEDAI. LOTUS ustekinumab ph3 (NCT03517722): measure titled "SRI-4 composite response", caught by the `composite response`
exclusion (an unintended casualty: SRI-4 is itself a composite; kept out as pre-registered). Brentuximab (NCT02533570):
title says "SRI response" without "4". TULIP-2 (NCT02446899) reports BICLA only, so it is in the BICLA table only.
