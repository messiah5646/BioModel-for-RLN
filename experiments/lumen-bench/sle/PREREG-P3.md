# PREREG P3 (+ Q58), round 115: does the LN forecaster's anchor rule generalise to non-renal SLE trials?

Committed before the SLE responder table is built (rule 12a push 1). Data step = `build_sle_table.py` (push 2).

**Decision this can change:** whether LUMEN's drug-effect layer uses the molecule -> target -> pooled anchor hierarchy
as a general rule (if it passes here too) or only as an LN-specific heuristic resting on n = 20 trials (if it fails).
Q58: whether "carry the early read-out forward" is kept as a default, and whether it depends on onset speed.

## Disclosure
Only outcome-measure TITLES, arm titles and designs were listed (55 candidate trials). No responder counts were printed.
However, the analyst (Claude) knows published headline results of several of these trials (e.g. belimumab BLISS,
anifrolumab TULIP, baricitinib SLE-BRAVE, tabalumab ILLUMINATE) from training. All rules below are mechanical and fixed
now, so that knowledge cannot steer choices, but the test is not blind in the strict sense.

## Data
`data/real/ctgov/lupus_with_results.json.gz` (committed snapshot). Include: allocation RANDOMIZED, a placebo-comparator
group, SLE (non-renal) population, results posted. Exclude by rule:
- trials whose active arms combine two investigational drugs or test two different drugs (NCT03656562 VAY736/CFZ533,
  NCT03978520 and NCT04451772 elsubrutinib/upadacitinib, NCT05306873 MMF/voclosporin);
- long-term extensions without a concurrent placebo arm (NCT03843125 baricitinib LTE);
- lupus nephritis trials (already in resp_table.csv).

**Endpoint (primary): SRI-4.** Measure title matches `SRI ?\(?4\)?|SRI-4|SLE Responder Index|Responder Index` and NOT
`[5-8]\)|SRI-?[5-8]|steroid|OCS|corticosteroid|prednis|sustain|time to|duration|interferon|IFN|test.?high|test.?low|
modified|mSRI|subgroup|stratum|composite response|with no|without|and OCG|dose ≤|dose <=|BICLA`.
**Endpoint (secondary analysis): BICLA**, title matches `BICLA|BILAG.{0,20}(based|Based) (Composite|Combined)` with the
same exclusions (except BICLA).
Counts: COUNT_OF_PARTICIPANTS taken as is; percentages converted with the arm denominator, half-up rounding
(as build_resp_table.py). Visit parsed from class title, then measure title, then time frame.
**Primary visit per arm = latest visit <= 52 weeks** with both arm and control present. Early visit (Q58) = latest
visit in 12-26 weeks. Denominator rule as build_resp_table.py (drop visits with < 50% of the arm's largest n).
Arm mapping: control = group title contains "placebo" and no dose pattern (`\d+(\.\d+)? ?(mg|mcg|µg|mg/kg)`) and no
molecule token; every other group = active arm of that trial's molecule. Groups matching `open-label|extension|\bLTE\b|long-term|safety follow|Part B|Part C|\(MEE\)|Cohort A|Phase 1b|pSS|
Sj[oö]gren|Rheumatoid|\(RA\)|Main Study\) /` are dropped (other diseases in basket trials, CLE-only cohorts,
extensions). Control groups titled "followed by" (placebo-to-drug crossover) are used only at visits <= 24 weeks.

## Labels (fixed now, by NCT)
Molecule and target per trial: see `LABELS` in build_sle_table.py, written in the same commit as this file.
Onset class for Q58 (mechanism, fixed now):
- **fast** = type I IFN axis (anifrolumab, sifalimumab, IFN-kinoid, litifilimab/BDCA2, daxdilimab/ILT7) and JAK/TYK2
  (baricitinib, deucravacitinib, brepocitinib/PF-06700841, upadacitinib).
- **slow** = B-lineage and T-B help: BAFF/APRIL (belimumab, tabalumab, atacicept, rozibafusp), CD20 (obinutuzumab),
  CD40L (dapirolizumab), BTK (evobrutinib, fenebrutinib, branebrutinib), CD80/86-ICOS (abatacept, acazicolcept),
  CD19 (obexelimab family, if present).
- **other** = all remaining (IL-12/23, IL-6, IL-2, cereblon, S1P, cannabinoid, RNase, peptide, CD30, IL-21, TLR7/8,
  ILT4 etc.). Not used in the Q58 split.

## Task D (as round 86), primary effect = log OR of SRI-4 at the primary visit
- D0 = 0; D1 = IV pooled mean of training arms; D3 = anchor (same molecule, else same target) IV mean shrunk to D1 with
  w = tau2 / (tau2 + se_anchor^2), tau2 DerSimonian-Laird on training arms; no anchor -> D1.
- **Primary: leave-one-trial-out (LOTO)**, scored on arms that have an anchor in other trials.
- Secondary: leave-one-molecule-out (LOMO): all trials of the molecule removed; anchors can only be same-target.
- Score: IV-weighted MSE of (observed - predicted); noise-corrected MSE reported.
- **Pass ("hierarchy generalises")**: LOTO D3 beats both D0 and D1 by >= 10% on anchored arms.
  **Fail**: D3 does not beat D1. Between: "partial". Report n trials / arms / anchored arms.
- Prior-art note: this is model-based meta-analysis (see claude/PRIOR-ART-trial-forecasting-2026-11-23.md); the claim
  tested is only whether the specific anchor hierarchy helps, not novelty.

## Task C' / Q58: early read-out -> week-52 log OR
Arms with an early (12-26 wk) and a primary (48-52 wk) SRI-4 visit. C0 carry-forward, C1 null (0), C2 shrunk
carry-forward to the pooled mean (w from early SE and tau2). IV-MSE.
- Overall: "early read-out fast-forwards" if C0 or C2 beats C1 by >= 10%.
- Split: report C0/C2 vs C1 separately for fast and slow arms. Prediction from the LN result (round 112, Hypothetical):
  carry-forward helps more for fast than slow mechanisms. **Supported** if the C2-vs-C1 gain is larger in fast than in
  slow and positive in fast; **not supported** otherwise. Small n per class is expected; label at most Emerging.

## Robustness (reported, not decisive)
BICLA instead of SRI-4 (trials that report BICLA); one arm per trial (highest dose) to remove within-trial correlation.

## Not allowed after the table is built
Changing regexes, labels, visit rules, exclusions or thresholds. Bugs in parsing may be fixed only if the fix is
reported with before/after counts (as round 107's correction 1).
