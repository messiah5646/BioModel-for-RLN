# Amendment A1 to PREREG-VC1 (round 142), committed before any real D28/M3 value was scored

PREREG-VC1 commit 0ea016d at 04:26:42 UTC. This amendment is committed at about 04:50 UTC, still before the 1-h mark (05:27).
Between the two commits, only code was tested, and only on **synthetic targets**. Every D28/M3 unit was rebuilt from the
same patient's *Base* cells plus a known injected shift (`run_vc1.py --synthetic neg|pos --quick`, 3 folds, 1 seed,
no inner tuning). No real post-treatment cell entered any score. Synthetic results are committed with this amendment.

## Why the change
The negative-control synthetic run (the shift is identical for every patient within a cell type, so no model should beat
B1 mean shift) showed a design defect in J1/J2 as written. "Base + decode(predicted latent delta)" cannot represent a
shift that lies outside the 50-component PCA span, while B1 represents it directly. The PREREG-as-written J2 scored 3.1x
worse than B1 there (0.0196 vs 0.0063). That is a property of the architecture, not of biology, and it would make the
real test a straw man.

## Change (models only; data, metric, folds, baselines, pass/FAIL rules and minimum n unchanged)
- **J1 and J2 become residual over the mean shift:** prediction = Base + B1 mean shift (training patients, same cell type)
  + decode(predicted patient-specific latent deviation). The latent models are trained on
  (z_target - z_source) minus the mean latent delta of the same (cell type, transition) in the training pairs.
  This is the VCC-2026 lesson (round 101): the model's job is to personalise the average response, not to relearn it.
- **J1orig and J2orig** (the PREREG-as-written forms) are kept and reported as secondary S10. They cannot ship.
- **Positive control (S11, code check, reported):** synthetic targets where the response size along baseline program v
  depends on the patient's own baseline score on a second program w. No trivial rule can copy this interaction.

## Code-check results (synthetic only, 3 folds, quick mode)
| Control | B1 | B3 (best trivial) | J1 | J2 | J2orig | noise ceiling |
|---|---|---|---|---|---|---|
| negative (no personal signal) | 0.00626 | 0.00626 | 0.00885 | 0.00672 | 0.01964 | 0.00518 |
| positive (baseline-dependent response) | 0.01356 | 0.01337 | 0.01268 | **0.01127** | 0.02526 | 0.00518 |

Reading: amended J2 does not invent personalised signal when none exists (7% worse than B1, the price of fitting noise).
It recovers a real personalised signal when one exists (16% better than the best trivial rule). So a real-data FAIL or
INCONCLUSIVE will mean there was no learnable personal signal at n 11, not that the model could not have found it.

## S6 (STAR-T transfer) is skipped, reason logged before the real run
The STAR-T `.h5mu` files (BD Rhapsody) carry only `Cell_Type_Experimental` (the sorted B/T fraction), not the fine
ontology labels the CD3xCD19 model is conditioned on. Mapping would need a new marker-based annotation. Doing that
without looking at post-treatment cells is a separate pre-registered task, not part of VC-1.
