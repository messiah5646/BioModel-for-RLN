# Pre-registration, round 88: what rescues refractory lupus nephritis? (2026-11-21)
Written after reading titles only (PubMed, 2005-2026); no abstracts or results read yet.

## Inclusion
Human cohorts or trials of refractory / resistant / persistent / relapsing-refractory LN, >= 8 patients with LN,
reporting complete renal response (CR/CRR) and/or overall response (CR + PR) at 6-24 months (take the visit closest
to 12 months). From the abstract; if the abstract lacks counts, the study is listed as unusable, not guessed.
Excluded: reviews, case reports < 8 patients, non-LN SLE-only cohorts, mouse studies, first-line trials.

## Groups (by mechanism, fixed now)
conventional_switch: CNI/multitarget, MMF/MPS, cyclophosphamide, leflunomide, iguratimod, cyclosporine
anti_CD20: rituximab, obinutuzumab, ofatumumab (alone)
BAFF_axis: belimumab, telitacicept (alone, or after rituximab)
RTX_plus_BAFF: rituximab combined with belimumab/telitacicept
plasma_cell: bortezomib, daratumumab
CAR_T: any CAR-T (CD19, BCMA, CD19/BCMA, CD19/22), allogeneic CAR or NK cell products
HSCT: autologous stem cell transplant
other: low-dose IL-2, irinotecan, etc.

## Analysis
Random-effects pooled proportion per group (logit, DerSimonian-Laird), CR and overall response separately.
Q1: Does deep depletion (CAR_T or plasma_cell) have higher pooled CR than anti_CD20 AND than conventional_switch?
Q2: Is anti_CD20 pooled CR higher than conventional_switch?
"Supported" only if the 95% CIs of the pooled proportions do not overlap; otherwise "not distinguishable".
Also report follow-up length per group, because CAR-T cohorts are young and short.

## Known biases (stated now)
Single-arm, selected, publication bias toward positive series, heterogeneous definitions of refractory and of response.
Response after ANY switch in RLN is high (Lupus 2021: 28/45 CR), so conventional_switch is the comparator,
not zero. The result can rank options; it cannot prove one beats another.
