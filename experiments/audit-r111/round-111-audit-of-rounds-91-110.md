# Round 111: audit of rounds 91-110. The code reproduces exactly; 15 claim labels are downgraded

**Decision this can change:** which claims from rounds 91-110 go into the preprint (P4).
- If a claim passes (n ≥ 10, pre-registered, the number traces to committed code): it enters the preprint at its label.
- If it fails: it is downgraded to Emerging, or kept out of the preprint until its n is filled in.

Pre-registration: experiments/audit-r111/PREREG-audit-r111.md (commit ef03a17), committed before any re-run.

## Verdict

**The audit FAILS by its pre-registered rule.** Every number reproduces, but labels were wrong:

- 15 Supported/Established claims break rule 10 or rule 12 and are downgraded below.
- 19 more literature claims give no n, so they are held out of the preprint until their n is filled in.

The code side is clean.

## A. Reproduction (scratch copy; committed files untouched)

| Check | Result |
|---|---|
| A1 `build_resp_table.py` from `lupus_with_results.json.gz` | 43 rows / 29 arms / 20 trials. Byte-identical to the committed `resp_table.csv`. `resp_table_plus_pub.csv` (46 rows) is identical too. PASS |
| A2 `compare_q17.py` vs round-86 `result86.txt` | P1/P2/P3 True. All 11 IV-MSE values are identical (e.g. D3 0.094, C0 0.224). P4 arm table 32/32 rows, max diff 0.005. PASS |
| A3 `check_frozen.py` (round 90 + W1) | Round-90 frozen forecasts are reproduced exactly: NOBILITY 0.507/0.300, TULIP-LN int 0.392/0.248, basic 0.131/0.254. The 13 W1 forecasts reproduce with max diff 0.004. PASS |
| Spot checks | r92, r98, r100, r102, r104 and r105 numbers match their committed result files, checked at the extra decimals. The pooled set-point multiplier 0.73 sits in `setpoint.py`/`RESULT.md`; only the per-arm values (0.688/0.745) are in `result.json`. |

Reproduction problems found:

1. **Round 91 has no committed analysis script.** `kidney-compartments/` holds only PREREG, RESULT and `result.txt`, so its numbers cannot be regenerated.
2. **Round 102 protected-cell count.** The doc says 787 cells; `c4_result.json` says `protected_kept` 630. This is cosmetic and changes no verdict.

## B. Claim ledger, counts

| Item | Count |
|---|---|
| Supported/Established claims in rounds 91-110 | 67 |
| Flagged claims (lists overlap) | 34 |
| Statistical comparisons run | ≈156, over 13 scored rounds |
| Rounds with 0 tests (literature or bookkeeping) | 7 |

With ≈156 comparisons, about 8 nominal p < 0.05 hits are expected by chance. No round-91-110 claim rests on a single p-value without a pre-registered rule, but this is why the exploratory signals below stay Emerging. They are:

- ribosomal fraction (r102)
- AZA IFN (r97)
- composition 0.74 in C5 (r104)

**Pre-registration timing.** Six PREREG commits land only 1-3 minutes before their result commits (r91, r96, r97, r99, r105, r107). Commit order cannot prove the analysis had not already been run. The labels stay, but the preprint must say that pre-registration is evidenced by commit order only. From now on, push the PREREG and then a separate data-download commit before the analysis.

## C. Downgrades (rule 10/12) — authoritative from this round

| Round | Claim | Was | Now | Why |
|---|---|---|---|---|
| 93 S1 | Baseline clinical + kidney PCs at chance (GSE200306) | Supported | Emerging | Round 91 called this Emerging; no analysis script is committed |
| 94 S2 | Rebuilt arm table reproduces round-90 priors | Established | Emerging | n = 2 trials, not pre-registered. Superseded by r107/r111 (n = 20 trials) |
| 98 S4 | Embedding gains track processing; avoid them in per-library designs | Supported | Emerging | Secondary analysis with no pass rule; n = 8 donors |
| 103 S2 | scANVI over-confident across labs | Supported | Emerging | n = 7 query samples, one lab |
| 103 S3 | KPMP "ncMON" are neutrophils | Established | Emerging | n = 7 samples, one lab |
| 104 S1 | scVI does not transfer (0.90 → 0.53) | Supported | Emerging | n = 6 patients. The general "scVI ≠ better patient classifier" claim stays Supported on r100 (n = 10) and r102 (n = 96) |
| 104 S2 | IFN replicates its non-result at flare | Supported | Emerging | n = 6, not blind (r89 had already scored this set) |
| 105 S1 | Month-3 UPCR leaves no room for exposure | Supported | Emerging | Aggregate bound with unknown ρ; the doc itself caps it at Emerging |
| 106 S1 | Severe disease = lower ceiling, not slower rate | Supported | Emerging | Exploratory, not scored, one cohort, composite grouping |
| 107 S2 | TULIP-LN year-2 population selected | Supported | Emerging | n = 1 trial |
| 97 S2 / 109 S8 / 110 S6 | "Blood IFN is a slow patient set point" | Supported | **Emerging** | Rests on r89 (n = 6) and r96 (n = 18 donors; set-point persistence H1 **failed**, Spearman 0.49). P5 (ILLUMINATE, 1,760 pts) is the test that can settle it |
| 109 S4 | CD19 CAR-T resets the loop | Supported | Emerging | DORIS, n = 5, uncontrolled |
| 110 S2 | Long-lived plasma cells = refractory reservoir | Supported | Emerging | Daratumumab n = 2 cases; CD19+BCMA CAR-T 12/15 uncontrolled. P6 is the test |

**Held out of the preprint until n is filled in (label unchanged).** These are literature claims that state no n:

- r93 S3, S4, S7, S8
- r95 S4b, S5, S6
- r96 S1
- r101 S2
- r109 S2, S3, S5, S6, S7
- r110 S1, S3, S4

**Kept but qualified.** r108 S2, "kidney IFN falls as chronicity rises", stays Supported (n = 156 biopsies, AMP). It is one controlled-access study and has no open replication (Q19 blocked).

## D. Outside the window but load-bearing for P4

The round-90 frozen week-104 "PASS" (MSE 0.055 vs carry-forward 0.242) reproduces exactly. However, it rests on **n = 2 trials (3 arms)**, so by rule 12 it is **Emerging**, not Supported.

The forecaster's Supported base is the leave-one-trial-out benchmark (round 86, n = 20 trials, 29 arms), and that also reproduces exactly. This is why P3 (widening to non-renal SLE trials) matters before the preprint.

## E. Claims cleared for the preprint (pre-registered, code committed, reproduced, n ≥ 10)

| Claim | Label | n | Round |
|---|---|---|---|
| Trial-response table and leave-one-trial-out forecaster benchmark | Supported | 20 trials / 29 arms | r86, r107, r111 |
| Geneformer adds nothing over plain PCA for pre/post B-cell depletion | Supported | 11 donors | r92 |
| IFN6 survives QC adjustment | Supported | 11 donors | r98 |
| scVI ≤ PCA, pre/post depletion | Supported | 10 donors | r100 |
| scVI ≈ pseudobulk, SLE vs control | Supported | 96 donors | r102 |
| Blood IFN6 adds nothing to UPCR (+SLEDAI) for MMF response | Supported | 44 patients | r97 S1 |
| Blood IFN6 does not flag nephritis beyond extra-renal activity | Supported | 218 patients, 2 cohorts | r99 |

## Follow-ups (decision-gated)

1. **Rule change** (process, no new item): commit the PREREG, then the data download, then the analysis, in three separate pushes. This makes pre-registration verifiable.
2. **Commit a regeneration script for round 91**, or drop its numbers from the preprint. Decision: whether the "interstitial myeloid/TGF-β" programme is cited at all. Added as Q57.

Sources:

- Committed code and outputs on branch rounds-70-72: `experiments/lumen-bench/resp_table/*`, `experiments/ifn-setpoint/*`, `experiments/scvi-protocol/*_result.json`
- Round docs claude/round-91 … round-110 in the project
