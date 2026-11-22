# PREREG Q11 (round 97): does a blood IFN / plasmablast signature predict MMF response beyond UPCR and SLEDAI?

Frozen 2026-11-22 BEFORE any response label or expression value was looked at.
What was seen before freezing: GEO series text, characteristic field NAMES, diagnosis / drug_group counts,
number of samples per patient and visit days (design only). Not seen: sri4_response values, any expression.

## Data
GSE224705 (PMID 38260685), whole blood, Affymetrix GPL13158, 448 samples, 183 people.
MMF: 130 samples / 44 patients. AZA: 54 samples / 20 patients.
IMPORTANT LIMIT, known in advance: the response label is **SRI-4** (global SLE response), not renal response,
and diagnosis is "SLE" for all patients (no LN flag). UPCR (mg/g) and SLEDAI are present per sample.

## Units and definitions
- Patient-level analysis. Baseline = each patient's earliest sample in that drug group.
- Outcome: if sri4_response is constant within a patient, use it; otherwise the label at the patient's LAST visit.
  Patients with no usable label are dropped (count reported).
- Features at baseline:
  - log10(UPCR + 1); SLEDAI (numeric).
  - IFN6 = mean of per-gene z-scores (z across all SLE baseline samples, both drugs) of IFI27, IFI44L, IFIT1, ISG15,
    RSAD2, SIGLEC1 (same 6 genes as rounds 89-96). Probe -> gene: GPL13158 annotation, highest-mean probe per gene.
  - PB4 = same construction for JCHAIN (IGJ), MZB1, TNFRSF17, TXNDC5.
- Model: L2 logistic regression (C = 1, features standardised inside each fold), leave-one-patient-out
  out-of-fold probabilities, AUC with responder = positive class.

## Hypotheses and pass rules
- **H1 (primary, MMF).** AUC(UPCR + SLEDAI + IFN6 + PB4) - AUC(UPCR + SLEDAI) >= 0.05 AND the paired
  patient-bootstrap (2,000 resamples of out-of-fold predictions) 95% lower bound of the difference > 0. Else FAIL.
- **H2 (replication, AZA).** Same rule in AZA patients.
- **H3 (trivial baselines, reported for every model).** Chance 0.50; SLEDAI alone; UPCR alone; IFN6 + PB4 alone.
- Sensitivity (reported, not a pass rule): MMF patients with baseline UPCR >= 500 mg/g (renal-active proxy).

## Expected
Non-responders are a minority (series text: 27 of 130 MMF samples), so ~10 MMF non-responders: CIs will be wide.
A FAIL here can mean "no added value" or "underpowered"; the bootstrap CI is reported so the two can be told apart.
SRI-4 requires a >= 4-point SLEDAI fall, so baseline SLEDAI is expected to carry signal by construction.

Research only, not clinical guidance.
