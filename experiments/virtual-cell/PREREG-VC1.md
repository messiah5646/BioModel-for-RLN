# PREREG VC-1 (round 142): a learned latent "virtual-cell" predictor of month-3 blood cell state after CD3xCD19 therapy

1. **Item / round / date.** VC-1, round 142, live session 2026-11-26. Thomas: "can't we train anything with the existing
   data, at least a proof of concept" -> chose the virtual-cell predictor, pre-registered; "use all data augmentation and
   create the best course of action". Under rule 12b this file is committed and pushed first; training and scoring run
   **>= 1 h after this commit**. Before then, the only allowed operations on the data are building pseudobulk tables and
   testing code on permuted / synthetic targets. No model is scored against real D28 or M3 values before the 1-h mark.

2. **Decision this can change (rule 7).** Whether LUMEN gets a learned cell-state ("virtual cell") layer that predicts
   a patient's post-therapy blood state from their own baseline.
   - **PASS** -> build LUMEN virtual-cell module v0 around the winning latent predictor (J), labelled Emerging (n 11),
     with STAR-T transfer as its first external check.
   - **FAIL / INCONCLUSIVE** -> LUMEN's cell-state layer is the best trivial rule (copy the patient's baseline, add the
     average treatment shift, or shrink to the population level, whichever wins), and no deep model is shipped.
     The project records that learned latent dynamics add nothing at n 11 (same lesson as VCC 2026, round 101).

3. **What has already been seen (disclosure).**
   - The same 11 patients (A-319 CD3xCD19 bispecific phase 1, CELLxGENE collection 6c147b03, Nat Med 2026) were analysed in
     round 96 (blood IFN6 set point: Spearman 0.49, shared multiplier 0.73) and round 100 (pre vs M3 classification;
     scVI 0.91 vs PCA 0.94 AUC, composition 0.81). So the existence and rough size of a pre/post difference is known.
     Gene-level month-3 prediction from baseline was never run.
   - The h5ad files were fetched on 2026-11-25 (round 137b harvest, `experiments/unblock-2026-11-25/MANIFEST-harvest.tsv`),
     before this PREREG. On 2026-11-26 only obs metadata was read (donor x timepoint x cell-type counts).
   - Label ceiling: **"12b-timed (data on disk before the PREREG; earlier looks not excludable)"**.
   - Rule 8: embeddings as patient-level *classifiers* are settled. This is a different task (predicting a future state,
     scored against carry-forward), so it is not a re-test.

4. **Data.**
   - Primary: CD3xCD19 CELLxGENE h5ads (T cells, B cells, Monocytes, Myeloid, LDG; md5 in the harvest manifest), 11 SLE
     patients (P1-P9, P11, P12) x Base / D28 / M3, `raw.X` counts, `cell_type` labels as deposited, `dose`.
   - Unit = (patient, compartment file, cell_type, timepoint) pseudobulk: sum of raw counts -> CP10k -> log1p. A unit is
     **eligible** if it has >= 30 cells. A (patient, cell type) pair is **scored** if both its Base and M3 units are eligible.
   - Genes: the 2,000 most variable genes across **Base-timepoint** pseudobulks only (no outcome data), fixed for all folds,
     excluding mitochondrial, ribosomal (RPL/RPS), hemoglobin and sex-chromosome genes (XIST, RPS4Y1, DDX3Y, ...).
   - Composition: fraction of each cell type among cells of its compartment file, per patient x timepoint.
   - External (secondary only): STAR-T (Zenodo 15715577, 5 patients, BL ... M6, sorted B and T `.h5mu`).

5. **Models.** All fitted on training patients only (leave-one-patient-out, 11 folds).
   - **B0 carry-forward** (strongest trivial baseline, rule 2b): predicted M3 = the patient's own Base pseudobulk.
   - **B1 mean shift:** Base + mean (M3 - Base) of the same cell type across training patients (VCC "average perturbation").
   - **B2 population mean:** mean M3 of that cell type across training patients.
   - **B3 shrink to common level:** a x (Base + B1 shift) + (1 - a) x B2, one a per cell type in [0, 1] (grid 0.05),
     fitted by inner LOPO on the training patients (a = 1 is B1, a = 0 is B2; round-119 lesson).
   - **J1 linear latent predictor (ablation):** PCA encoder (50 components, fitted on training pseudobulks of all
     timepoints) -> ridge regression of z(M3) - z(Base) on [z(Base), cell-type one-hot, dose]; decode and add to Base.
   - **J2 LUMEN-JEPA (the challenger):** the same PCA encoder; an MLP predictor (2 hidden layers of 128, ReLU, dropout 0.2,
     weight decay 1e-3) that takes [z(source), cell-type embedding, dose, delta-t in days] and predicts z(target) - z(source).
     The loss is MSE in latent space, plus 0.1 x MSE after decoding to genes. Ensemble of 5 seeds. Prediction:
     Base + decode(predicted delta-z).
   - **Augmentation (J1 and J2, training patients only; never applied to the held-out patient's scored values):**
     (a) all three transitions as training pairs: Base->D28, D28->M3 and Base->M3, with delta-t as an input;
     (b) cell bootstrap: 20 pseudobulk replicates per eligible unit, resampling its cells with replacement, plus 10
     replicates from 50% subsamples, so the model sees realistic sampling noise;
     (c) Gaussian noise in the latent space (sd 0.1 x per-component sd) on source embeddings;
     (d) pair mixing: random source and target replicates of the same (patient, cell type).
     Hyperparameters above are fixed. There is no tuning on held-out patients; any inner tuning uses inner LOPO on the 10
     training patients, and only for J2's epoch count (grid 50/100/200).

6. **Metric and held-out scheme.** Leave-one-patient-out. **Primary metric:** for the held-out patient, the MSE over the
   2,000 genes between predicted and observed M3 pseudobulk, averaged with equal weight over the patient's scored cell
   types; then averaged over the 11 patients. Observed M3 uses all of the held-out patient's cells.
   - **Noise ceiling (reported, not a model):** split each held-out M3 unit's cells in half 20 times; the mean MSE between
     half-pseudobulks divided by 4 (variance of a half mean is 2x, of a half-half difference 4x, that of the full mean), is the floor no model can beat.

7. **Minimum n.** >= 9 patients with >= 3 scored cell types each, else UNDERPOWERED.

8. **Pass rule.** Let A(m) be the primary metric and best = min(A(B0), A(B1), A(B2), A(B3)).
   **PASS** if A(J2) <= 0.90 x best AND J2 beats the best baseline's per-patient MSE in >= 9 of 11 patients.

9. **FAIL rule (rule 2b).** A(J2) >= 1.10 x A of ANY trivial baseline (B0-B3) = **FAIL**. Anything between = INCONCLUSIVE.
   If J1 passes and J2 does not, the shipped candidate is J1 (same rules applied to J1); this is the only fallback.

10. **Secondaries (labelled; none can turn a FAIL into a PASS).**
    - S1 delta correlation: Pearson r between predicted and observed (M3 - Base) over genes, per scored unit (B0 has no
      delta; reported for B1, B3, J1, J2), and the VCC-style score (0 = B1 mean shift, 1 = noise ceiling).
    - S2 composition: MAE of predicted M3 cell-type fractions (B0 / B1 / B2 / J2 composition head).
    - S3 D28 as target (Base -> D28), same models.
    - S4 per-compartment results (T, B, monocyte, myeloid, LDG) and per-cell-type wins.
    - S5 without augmentation (J2-noaug) to show what augmentation buys.
    - S6 STAR-T external transfer (descriptive): models trained on all 11 CD3xCD19 patients predict STAR-T M3 (or nearest
      visit) from BL for T and B cell types whose labels can be mapped by a frozen marker rule. If the labels cannot be
      mapped without looking at post-treatment cells, S6 is skipped and the reason is logged.
    - S7 patient bootstrap (2,000 resamples, seed 142) of A(J2) / A(best).

11. **Labels (rule 12).** n = 11 patients, so **Emerging at most**, whatever the verdict, carrying the 12b-timed label.

12. **Output.** `experiments/virtual-cell/run_vc1.py` -> `vc1_result.json` via `scripts/provenance/stamp.py`, plus
    `vc1_per_patient.csv`. Pseudobulk builder: `build_pseudobulk.py` (committed before scoring).

Research and simulation only. Not clinical or drug-development guidance.

## Input md5 (frozen at PREREG)
- 5c869ad5331eb7a39834c71abca19b65  data/real/harvest-r137/CD3xCD19_Bcells.h5ad
- 0cf6282bce02d897c115abdec98d8910  data/real/harvest-r137/CD3xCD19_LDG.h5ad
- ee2e06cd8ddafb9ba515567894c483d6  data/real/harvest-r137/CD3xCD19_Monocytes.h5ad
- 5fe73c9630dcf744d53d8fe7a881449f  data/real/harvest-r137/CD3xCD19_Myeloid.h5ad
- dd61442c0fe93059517ea3b2d5bfd38b  data/real/harvest-r137/CD3xCD19_Tcells.h5ad
- ffae0e3ab1a47a56ce51bb1cc08d1163  data/real/harvest-r137/STAR-T_scRNA-seq.zip
- ca14be70b213b98bbc2b830e2b092329  scripts/provenance/stamp.py
- 4165bab8f72dee66fd8929e3edff191d  scripts/provenance/fetch.sh
