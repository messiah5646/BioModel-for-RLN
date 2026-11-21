# Rounds 85-87: what transfers between lupus nephritis trials (2026-11-21)

## Scored as pre-registered
| Round | Test | Result | Verdict |
|---|---|---|---|
| 85 A | Control arm of unseen trial (proteinuria scale) | 2-param curve 0.515 vs constant 0.558 (-8%) | fail |
| 85 B | Add-on effect from mechanism class (proteinuria) | 0.377 vs pooled 0.379 | fail |
| 85 C | Cell-derived delay improves early->late projection | 0.434 vs carry-forward 0.432 | fail |
| 86 D | Anchored on same molecule/target (responders, 18 anchored arms) | IV-MSE 0.094 vs null 0.309, pooled 0.167, class 0.104 | fail by 0.4 points (beats null -70%, pooled -44%) |
| 86 C' | Fast-forward: carry week 24-26 effect to week 48-96 (13 rows) | IV-MSE 0.224 vs null 0.577 | **pass (-61%)** |
| 87 | Frozen predictions on unseen trial (Liu 2015 tacrolimus multitarget vs IVC, n=362) | observed log OR 0.90 (SE 0.23); class 0.59 closest, pooled/anchor 0.33, null 0 | class prediction closest (n=1 trial) |

## Exploratory (not pre-registered), and decisive for interpretation
Leave-one-molecule-out (29 arms): pooled 0.275, null 0.330, mechanism class 0.349.
Mechanism class carries no transferable information about a NEW molecule; its apparent value in 86 D and 87 came
from voclosporin's sibling trial (and tacrolimus shares voclosporin's target). Same conclusion as rounds 50-84.

## What transfers
1. The same trial's early read-out (week 24-26 CRR effect) -> its week 48-96 effect. Pre-registered pass.
   Fails for delayed-acting drugs: NOBILITY obinutuzumab log OR -0.10 at week 24, +0.60 at week 52.
2. The same molecule or the same molecular target in another trial (anchor). Large gain over null and pooled.
3. Not mechanism class alone. Not mechanistic simulation of magnitude.

## Consequence for the architecture
The drug-effect layer must be a hierarchical evidence model (molecule -> target -> global) on standardised
responder endpoints, updated by early human read-outs. Mechanism and cell layers are demoted from predicting
magnitude to: (a) deciding which evidence is exchangeable (same target and same depth of engagement),
(b) setting the read-out horizon (delayed B-cell / BAFF effects need >= 52 weeks), (c) choosing the early
biomarker (monocyte IFN at month 3 after B-cell depletion, replicated in 2 trials).
