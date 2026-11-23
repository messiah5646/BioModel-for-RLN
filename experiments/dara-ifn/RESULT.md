# P2b RESULT: blood T-cell IFN does NOT follow the autoantibody pool after daratumumab (H1 not supported; weak evidence)

PREREG: experiments/dara-ifn/PREREG.md (commit 1e38b5c). Download: commit 9e7423a. Analysis: analyse.py -> out/.

## Pre-registered verdict
| Rule | Result | Verdict |
|---|---|---|
| R1 Spearman(log anti-dsDNA w36/w0, log IFN6 m36), n = 5 (P1, P3, P5, P7, P9) | rho = -0.30 (QC-adjusted -0.80); perm p 0.74 | **H1 fails** (rho <= 0); no QC flip |
| R2 B-cell explanation | median m36 = 1.10; B-fraction ratio 2.76 | **not rejected** (no IFN fall while memory B persist) |
| R3 fraction with m36 < 1 | 0.40 | below the 60% H1 bar |
| Week 9 (steroid-confounded, secondary) | median m9 = 0.97 (n = 6) | no clear early fall in this score |

## Deviations found during the analysis (all disclosed; none chosen after seeing IFN results)
1. **GEO labels of two pooled libraries contradict the TCR data.** TCR clonotype overlap identifies samples
   (same patient, other visit: 6-19% shared TRB CDR3; different patients: ~0%; out/label_audit_tcr_overlap_pct.csv).
   - GSM8919760, titled "P5 w36 + P10 w0", shares 0% with P10 w9 but 3-5% with the P5/P7 pools.
   - GSM8919763, titled "P7 w36 + P9 w9", shares 0% with P9 (w0, w9, w36) but 3.5% with P10 w9.
   So the pre-registered demultiplexing (which trusted the titles) assigns P5/P7 w36 cells unreliably, and P10 has
   no usable w36 sample. The pre-registered n = 5 therefore includes two unreliable patients.
2. **Not deposited / not served:** GSM8919775 (P8 + P10 w36) has no TCR file, so P10 drops by rule A2.
   GSM8919766 (P9 w9) h5 returns HTTP 403 from GEO (out/result.json `skipped_libraries`).
3. **P2 w36** got only 20 clonotype-matched T cells (P2 clonotypes recur weakly: 6% w0-w9), so it drops by A2.
   Hashtags in that library map 1:1 to clonotype-assigned patients (concordance 0.993), giving 193 P2 w36 T cells.

## Not pre-registered: reliable subset (P1, P3, P9 from clonotype-confirmed libraries; P2 via hashtag), n = 4
| Patient | anti-dsDNA w36/w0 | IFN6 m36 | m36 QC-adjusted | T cells w36 |
|---|---|---|---|---|
| 1 | 0.29 | 1.10 | 1.02 | 967 |
| 2 | 0.07 | 1.78 | 1.96 | 193 |
| 3 | 0.48 | 1.11 | 0.98 | 9,618 |
| 9 | 1.36 | 0.95 | 0.75 | 2,116 |
Spearman rho = -0.8 (QC-adjusted -1.0). The patient whose antibody fell most (P2) had the largest IFN rise; the one whose
antibody rose (P9) had the lowest IFN. This is the opposite of H1.

## Interpretation
- With the antibody pool reduced 50-90% in most patients and dexamethasone stopped ~28 weeks earlier, blood memory
  T-cell IFN6 at week 36 was not lower than baseline and did not track the antibody change. This fits round 89/96
  ("blood IFN is a patient set point") better than H1.
- It does not contradict the paper's early (week 9) fall: that window is confounded by dexamethasone, and the
  paper's claim is a cluster DEG analysis, not this per-patient 6-gene score.
- **Label: Emerging (against H1).** n = 4-5, single-arm refractory SLE, unverified assumption A1 (GEO patient numbers
  = paper SLE #N), sorted memory T cells only (no monocytes/pDCs, the main IFN producers and responders).

## Decision
- Do **not** build LUMEN's IFN state on the antibody pool (plan stop rule). Keep the round-96 multiplier model.
- H1 (SYNTHESIS) is downgraded from Hypothetical-with-a-test to **not supported for blood IFN**. The kidney-timing part
  of H1 (B-lineage drugs separate late) is untouched and is tested by the frozen VIBRANT forecast (H3).
- P6 (plasma-cell two-pool test) no longer depends on H1 and can run on its own question (refractoriness).
- Report the GEO label inconsistency to the data owners when Thomas does outreach (not contacted by the loop, rule 5).

## Round 114 self-audit (2026-11-23): two checks, one correction
1. **Assumption A1 (GEO "Patient N" = paper "SLE #N") is consistent with independent data** (n 5; audit r115 wording). CD38 mRNA in TCR+ CD8-negative cells at week 0
   ranks patients like the paper's flow cytometry (% CD38+ CD4 T, Fig 3b source data): Spearman 0.90, exact p 0.042,
   n = 5 (P3 lowest and P9 highest in both). `check_A1_cd38.py`, out/check_A1_cd38.txt. (CITE-seq CD38 protein could not
   be used: P9's week-0 library has a different 344-antibody panel.)
2. **Gene set is not driving the result.** The paper's own T-cell IFN genes (MX1, IFI44L, XAF1, ISG15, IFIT1) give the
   same week-36 answer on clonotype-confirmed samples: median m36 1.13, rho vs anti-dsDNA -0.80 (n 4).
   `check_geneset.py`, out/check_geneset.txt.
3. **Correction.** The table above says week 9 showed "no clear early fall" (median 0.97). That median included the
   mislabelled P5/P7 samples. On clonotype-confirmed samples (P1-P4) median m9 is 0.82-0.85: P3 and P4 fall ~35%,
   P1 and P2 do not. So there is a partial early fall (consistent with the paper) that is gone by week 36. That pattern
   fits a transient dexamethasone/acute-drug effect, not the antibody pool. The H1 verdict is unchanged.
