# Pre-registration, round 91: kidney compartments and treatment response (GSE200306), 2026-11-21
Written after reading the GEO series description and characteristic field NAMES only; no expression values and no
per-sample clinical values seen.

Data: protocol biopsies at flare (Bx1) and after treatment (Bx2), 58 proliferative LN patients; laser-microdissected
glomeruli and tubulointerstitium (TI), NanoString 574 immune transcripts (normalised files from GEO);
clinical response group CR / PR / NR.

H1 (the compartment contradiction from the D3 RLN map; Emerging literature: glomerular CD68 -> better response,
J Autoimmun 2025; interstitial CD68 -> worse response, Front Immunol 2024):
  H1a: at Bx1, glomerular CD68 is higher in CR than NR.
  H1b: at Bx1, TI CD68 is higher in NR than CR.
  Each passes if the direction is right AND one-sided Mann-Whitney p < 0.10 (small n, stated in advance).
H2 (does kidney molecular state add to clinical baseline?): outcome CR vs not-CR, Bx1 only.
  M0 = logistic(baseline proteinuria, creatinine); M1 = M0 + first 3 PCs of TI expression; M2 = M0 + first 3 PCs of
  glomerular expression. Leave-one-out AUC. PCs fit inside each fold. Pass if M1 or M2 beats M0 by >= 0.05 AUC.
H3 (refractoriness in tissue): genes most different between NR at Bx2 and CR at Bx2 in TI, reported descriptively
  (top 10 by effect size); no pass rule. This feeds the RLN map, labelled Emerging.
