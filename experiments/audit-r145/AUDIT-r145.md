# Audit round 145 (rule 10): rounds 141-144

Scheduled run 2026-11-26 ~15:25 UTC. Clean worktree at 690a950 (branch rounds-70-72).

Decision this can change: which r141-r144 labels enter the preprint (P4) and whether rounds-70-72 is ready to merge.
If a load-bearing number does not reproduce or a label overstates n / pre-registration -> downgrade and correct before
the preprint carries it; if all reproduce -> rows enter preprint 3.2/3.5 at their current labels.

## Reproduction
| Check | Result |
|---|---|
| r141: `regdocs/run_rt4.py` in clean worktree -> rt4_result.json | identical except run_utc / git_head; primary 10 trials / 17 arms: M0 0.2531, M1 0.3219, M2 0.2277, C0 0.1096; DEFAULT-FAIL pattern; minimum n not met -> UNDERPOWERED; RT-4-new trials in primary: LOTUS only |
| r142/r143: 5 CD3xCD19 h5ads re-fetched via fetch.sh (MANIFEST-a145-refetch.tsv, 15:25-15:27Z) | all 5 md5 = PREREG-VC1 = r143 refetch |
| r143: `build_pseudobulk.py` | 1,100 units, 762 eligible, 437,558 cells, 237 Base units for HVG: same as r143 |
| r143: `vc1_descriptive.py` -> vc1_descriptive.json | 974 leaves, identical except run_utc / git_head (B2base 0.0356 vs B2 0.0321 at M3; D28 S3; dispersion pairs; without-P9 blocks) |
| r142: vc1_result.json means vs committed vc1_per_patient.csv | all 10 model means match (B0 0.0535, B1 0.0584, B2 0.0321, B3 0.0332, J2 0.0423, J2orig 0.0392, ceiling 0.0051); J2 beats B2 in 0/11; ratio 1.318 |
| r142: J2 re-trained, fold 1 (P1) only (`vc1_fold1_check.py`, same rng order as run_vc1.main; 444 s) | B2, J1, J2, J2noaug, J2orig per-patient MSE identical to committed row (rel. diff <= 3e-15), epochs 50 = 50. Full 11-fold re-train (70 min) not repeated |
| r144: CLI `lumen_forecast.py obinutuzumab CD20 --early 24 40 100 25 100` | headline "NONE: no validated rule ..."; M0/M1/C0 side by side (pooled OR 1.38 [0.75, 2.51]; carry-forward 1.98); RT-4 line present in early_note |
| r144: forecast() unchanged | a24cea8 diff touches only the `__main__` block; `resp_table/check_frozen.py` reproduces P4 arm table (32/32), r90 frozen (3/3) and W1 13 frozen forecasts (max diff 0.004) |
| PREREG-VC1, AMENDMENT-VC1-A1, PREREG-RT4 edited after commit? | no (one commit each: 0ea016d, 674d5fa, 67f2d26) |

## Tests counted, rounds 141-144
r141: 1 primary decision comparison (RT-4; 4 models) -> UNDERPOWERED; ~20 descriptive sensitivity sets (S1, S3, S5, S5b, S6, S7, S8, S8b x9, S9).
r142: 1 primary (VC-1 J2 vs best trivial baseline) -> FAIL; J1 secondary FAIL; S4 by compartment, S7 bootstrap, 2 synthetic controls.
r143: 4 LOPO tables + 3 dispersion pairs, all descriptive, post hoc. r144: 0 tests.
**2 primary decision comparisons, both negative (UNDERPOWERED, FAIL); ~30 descriptive.**

## Labels
| Claim | n | Pre-reg | Label | Audit |
|---|---|---|---|---|
| r141: RT-4 UNDERPOWERED; descriptive CF 0.110 < M2 < M0 < M1 | 10 trials / 17 arms | PREREG-RT4 + 2 amendments | Emerging, 12b-timed (late outcomes seen) | correct; reproduces |
| r141: "3 of 10 documents print usable early values" | 10 documents | - | plain count (r144 removed "Established") | correct |
| r142: VC-1 learned predictor FAILS vs population average | 11 patients | PREREG-VC1 + A1 | Emerging, 12b-timed | correct; **wording added:** A1 changed the challenger model 41 min before the real run, so under rule 12b as amended in r144 the clock would have restarted. A1 made J2 worse (0.0423 vs as-written 0.0392) and both forms FAIL, so no label or verdict change |
| r143: small visit-specific shift / "convergence" | 11 patients | post hoc | Hypothetical | correct (B2base/B2 1.11 reproduces; reviewer's 1.10 bar) |
| r144: CLI display + rule 12b amendment sentence | - | - | process | correct |
No claim in rounds 141-144 is labelled Supported or Established; nothing rests on n < 10 above Emerging. **No downgrades.**

## Notes
- The round-142 project file name still carries the withdrawn "converges" headline (kept for links; corrected header inside). Renaming would need a delete, which the loop may not do. Thomas may rename it.
- where-this-stands (project) was only through r135; brought through r145 this round. Preprint v0.1.9 adds the RT-3 row (3.2) and VC-1 row (3.5).

## Verdict
PASS on code (r141 and r143 reproduce exactly; r142 means and a re-trained fold reproduce exactly; r144 CLI and frozen
forecasts check), PASS on labels with one wording addition (VC-1 A1 timing). Branch rounds-70-72 is ready to merge (Thomas reviews).
