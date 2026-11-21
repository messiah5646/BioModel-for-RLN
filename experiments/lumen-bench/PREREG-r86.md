# Pre-registration, round 86: drug-anchored bridging on renal responder data (2026-11-21)

## Why
Round 85: mechanism class (fast/slow) does not predict add-on effect; trial-to-trial heterogeneity is ~0.29 log units
above a ~0.25 noise floor. Round 82 and round 85 task C both point the other way: a drug's OWN earlier data predict
its later data. Hypothesis: the useful simulator is anchored on the drug's own human data, not on mechanism.

## Honesty note
Responder counts were seen during extraction (/tmp/resp_dump.txt). Model forms, the anchor hierarchy, class labels,
endpoint rule and scoring rule are fixed here before the table is scored. Nothing tuned after.

## Data
CRR (complete renal response, protocol definition) counts per arm per visit, 20 placebo/active-controlled add-on LN
trials from ClinicalTrials.gov results (resp_table.csv). Excluded by rule: no responder counts (atacicept, daxdilimab),
n<5 single-arm or open-label (ACTHar, NKT-01, BMS-986165, etanercept, MISSION), maintenance design (NCT03385564),
AURORA 2 continuation (selected subset).
Effect = log odds ratio (drug vs control), 0.5 continuity correction, SE from counts.
Primary endpoint per arm = latest CRR visit <= 104 weeks.

## Labels (mechanism only, fixed now)
fast = renal-direct: voclosporin, ravulizumab (C5), ALXN2050 (factor D), BIIB023 (TWEAK, acts on renal parenchyma)
slow = immune-upstream: all others.
Anchor hierarchy: same molecule > same molecular target (CD20: rituximab/ocrelizumab/obinutuzumab; CD40:
iscalimab/BI 655064) > none.

## Task D: predict each arm's primary log OR, leave-one-trial-out
- D0: 0
- D1: inverse-variance (IV) pooled mean of all other trials' arms
- D2: IV mean of same mechanism class in other trials
- D3: anchored: IV mean of anchor arms in other trials, shrunk to D1 with w = tau2 / (tau2 + se_anchor^2);
      tau2 by DerSimonian-Laird on the training arms; no anchor -> D1
Score: IV-weighted MSE of (observed - predicted), and the noise-corrected MSE (mean of err^2 - se^2), on
(b) the subset of arms that HAVE an anchor [primary] and (a) all arms [secondary].
Claim "drug-anchored bridging works" if D3 beats D0, D1 and D2 by >= 10% on the primary score.

## Task C': fast-forward within an arm on responder data
Arms with CRR at <= 26 wk and at >= 48 wk. Predict late log OR from the latest early one.
- C0: carry forward
- C1: 0
- C2: shrunk carry-forward to the pooled mean (w from early SE and DL tau2)
Claim "early read-out fast-forwards" if C0 or C2 beats C1 by >= 10% (IV-weighted MSE).

## Content criterion
If D3 passes: the simulator's drug-effect layer should be anchored on the drug's own human data (Bayesian bridging).
If C' passes: early read-outs can shorten trials. If neither passes, the round names what is missing.
