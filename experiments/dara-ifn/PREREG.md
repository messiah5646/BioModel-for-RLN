# PREREG P2b (= Q41, re-aimed by SYNTHESIS H1): does blood IFN follow the autoantibody pool? GSE294789, daratumumab

Committed BEFORE any expression matrix of GSE294789 is downloaded (rule 12a, push 1 of 3).

**Decision this can change:** whether LUMEN's IFN state is driven by the autoantibody pool (H1) or by the B-cell
compartment. If H1 passes -> the IFN set point gets the antibody/plasma-cell pool as its input and P6 is run next.
If H1 fails -> keep IFN as a B-lineage-depletion multiplier (round 96) and drop the "one slow clock" framing for IFN.
If inconclusive -> no model change; H1 stays Hypothetical.

## What was seen before this PREREG (disclosure)
1. GEO metadata (sample titles, files). Design differs from the synthesis doc and the earlier project draft:
   - Post-treatment samples are **week 9 and week 36**, not week 36 only.
   - **No library pools one patient's baseline with that patient's later sample.** Pooled libraries mix different
     patients (P5w0+P7w0, P1w36+P2w36, P5w9+P7w9, P5w36+P10w0, P7w36+P9w9, P8+P10w36). The planned within-library
     paired contrast does not exist.
2. The paper (Ostendorf et al., Nat Commun 2026, doi 10.1038/s41467-026-69112-w; NCT04810754), main text:
   - It already reports that the type I IFN signature in T cells **fell early** after treatment. The aggregate
     week-9 direction is therefore NOT blind. I found no week-36 IFN statement in the text.
   - **Dexamethasone premedication** with every dose through week 8 (20 mg, then 10 mg). Week 9 is one week after
     the last steroid dose. Glucocorticoids lower IFN signatures by themselves.
   - Dendritic-cell counts stable; NK cells fell; circulating CD19+ B cells stable; ASCs depleted.
3. Source data sheets **Figure 2a1/2a2 only** (per-patient anti-dsDNA, IU/ml). The Figure 3 sheets (which may hold
   IFN data) were NOT opened and will not be opened before the analysis commit.
   Anti-dsDNA ratio week 36 / week 0 for the patients with scRNA at both visits:
   P1 0.29, P2 0.07, P3 0.48, P5 0.42, P7 0.42, P9 1.36, P10 0.60.

## Assumptions to check in the analysis commit (report if violated)
- A1. GEO title numbering "Patient N" equals the paper's "SLE #N". Check: none possible from expression; stated as
  an assumption. If the authors' numbering differs, the R1 correlation is invalid -> label inconclusive.
- A2. Pooled libraries are split by **TCR clonotype matching** (patient-private CDR3 nucleotide sequences) against each
  patient's single-patient libraries. Cells whose clonotype is found in exactly one candidate patient's reference are
  assigned; others dropped. A patient-sample with < 50 assigned T cells is dropped. Hashtag (Antibody Capture) data,
  if present, are reported as a cross-check only.

## Score (frozen, same as rounds 96/98)
IFN6 = per-cell mean of log1p(CP10K) over IFI27, IFI44L, IFIT1, ISG15, RSAD2, SIGLEC1; cells with total counts > 500.
**T cells** = barcodes with a productive TCR contig and no productive BCR contig. **Memory B** = productive BCR, no TCR.
Per patient-sample score = mean over its T cells. Multiplier m = IFN6(week t) / IFN6(week 0).

## Primary analysis: week 36 vs week 0 (steroid-free by ~28 weeks)
Patients eligible: those with week 0 and week 36 T cells after A2 (expected up to 7: P1, P2, P3, P5, P7, P9, P10).

- **R1 (H1 test):** Spearman rho between log anti-dsDNA ratio (w36/w0) and log IFN6 multiplier m36, across patients.
  - **H1 supported** if rho >= 0.5 AND n >= 5.
  - **H1 fails** if rho <= 0.
  - Otherwise **inconclusive**. Exact permutation p reported, not used for the decision.
- **R2 (B-cell-compartment alternative):** daratumumab spares memory/circulating B cells, so a B-cell-driven IFN predicts
  median m36 near 1. **B-cell explanation fails** if median m36 <= 0.80 AND median memory-B cell fraction ratio
  (w36/w0, among sorted cells) >= 0.8. If the B-cell fraction falls below 0.8, R2 is inconclusive.
- **R3 (old rule, kept as secondary):** fraction of patients with m36 < 1. Reported; >= 60% is consistent with H1, but
  not decisive on its own.

## Secondary: week 9 (steroid-confounded; reported, never decisive)
m9 per patient, and Spearman with anti-dsDNA ratio w8/w0. If median m9 < median m36 by > 0.15, the early fall is
read as partly steroid- or drug-acute, not antibody-driven.

## QC sensitivity
Repeat R1 and R2 with IFN6 residualised on log total counts, log genes and mitochondrial fraction (per-cell linear
model fitted on all T cells). If the R1 verdict flips, report both and label the result inconclusive.

## Labels
Any outcome is at most **Emerging**: n <= 7, single-arm, refractory SLE (9/10 with renal history, not active LN
selection), dexamethasone co-treatment, partly unblinded (see disclosure).

## Not allowed after unblinding
Changing genes, cell definitions, visits, thresholds, or dropping patients except by A2's < 50-cell rule.

## Companion (descriptive, no rule): GSE263931 CAR-T, n = 7 — deferred to a separate PREREG.
