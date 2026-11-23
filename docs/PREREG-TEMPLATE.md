# PREREG template (rules 1, 2, 2b, 12, 12a, 12b)

Copy into `experiments/<dir>/PREREG-<item>.md`, fill in, commit and push, then END the run (rule 12b).

1. **Item / round / date (UTC).**
2. **Decision this can change (rule 7):** <X>. If result A -> do <a>; if result B -> do <b>.
3. **What has already been seen** (data, tables, earlier rounds touching the same outcome).
4. **Data** (exact files / accessions; downloads via `scripts/provenance/fetch.sh`; frozen snapshot hash if no download).
5. **Models compared**, including ALL trivial baselines. Name the **strongest available trivial baseline**
   (carry-forward whenever the patient's or trial's own earlier value exists; else pooled/population mean).
6. **Metric and held-out scheme** (LOTO / time-forward / patient-held-out). Primary metric only one.
7. **Minimum n** below which the result is "underpowered", not a pass.
8. **Pass rule.** Model must beat the strongest trivial baseline by >= <margin>.
9. **FAIL rule (rule 2b, mandatory).** Losing to ANY trivial baseline by >= 10% on the primary metric = **FAIL**,
   not inconclusive. Anything between pass and FAIL = inconclusive.
10. **Sensitivities / secondaries** (labelled; none can turn a FAIL into a pass).
11. **Label rule** (rule 12): Supported/Established only with n stated; n < 10 -> Emerging at most.
12. Result JSON written through `scripts/provenance/stamp.py`.
