# AMENDMENT to PREREG-RT4 (reviewer flag R-20261125-1, round 139, 2026-11-25)

Written and pushed BEFORE any RT-4 run-3 fetch: `experiments/lumen-bench/regdocs/sources/` does not exist and no
`fetch.sh` call for RT-4 has been made at this commit. PREREG-RT4.md (67f2d26) is left unedited.

## Was the reviewer right? Yes (checked against the frozen tables, rule 13)
- Table A trials by molecule (`resp_table.csv`, `sle/sri4_table.csv`, `sle/bicla_table.csv`, CT.gov for MUSE/LOTUS,
  RT-3 for telitacicept): **belimumab 5** (BLISS-LN, BLISS-SC, BLISS-NEA [table name belimumab-NEA], BLISS-76, PLUTO),
  **anifrolumab 2** (TULIP-2, MUSE; TULIP-1 is also anifrolumab and already early in RT-1), abatacept 1 (ACCESS),
  ustekinumab 1 (LOTUS), telitacicept 1.
- RT-1 primary (7 trials: ALLEGORY, ALXN2050, BI655064, ORCHID-LN, SANCTUARY, fenebrutinib-ph2, rozibafusp-ph2) has 7
  distinct molecules. So the merged primary could carry 4-5 belimumab trials, and the ">= 5 RT-4-new trials" minimum
  (section 8) could be met with belimumab + anifrolumab alone. That is the round-112 voclosporin failure mode
  (signal resting on one molecule), and PREREG-RT4 has no molecule-level sensitivity (S5 resamples trials).

## Added (descriptive; models, metric, minimum n, decision rules and the shipped branch are UNCHANGED)
- **S8 molecule collapse.** Primary subset re-scored with all belimumab trials collapsed to one pseudo-trial and all
  anifrolumab trials to another: per-model IV-MSE where each collapsed molecule contributes the IV-weighted mean of its
  arms' squared errors with weight equal to the sum of its arms' weights divided by its number of trials (i.e. that
  molecule counts as one trial). Also **S8b leave-one-molecule-out**: primary IV-MSE for M0/M1/M2/C0 with each molecule
  removed in turn; report which section-9 rule each variant would select.
- **S5b molecule-level bootstrap.** 2,000 resamples of molecules (seed 139), descriptive 95% intervals of A(C0)/A(M1)
  and A(M2)/A(M1).
- **Label rule.** The round doc's first line, `rt4_result.json` and any CLI `--early` note state
  "k of n primary trials belimumab, j anifrolumab". If S8 or any S8b variant selects a different section-9 rule than
  the primary, they also state "rests on <molecule>", and the label stays Emerging with n per molecule.
- `rt4_result.json` carries keys `S8_molecule_collapse`, `S8b_lomo`, `S5b_molecule_bootstrap`.

## Disclosure added at the reviewer's request (R-20261125-2, applies to the RT-4 label)
Round-137 sub-agents opened the FDA reviews / EMA EPARs / PMC full texts to find table and figure locations (their
transcripts were not committed); earlier sight of early values and model training knowledge of published BLISS/TULIP
results cannot be excluded. Label stays "12b-timed (late outcomes already seen; early values new)", Emerging.

Research and simulation only. Not clinical or drug-development guidance.
