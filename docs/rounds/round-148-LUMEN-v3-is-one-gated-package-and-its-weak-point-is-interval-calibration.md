# Round 148 — LUMEN v3 is one gated package, and its weak point is interval calibration

*2026-11-27, live session with Thomas. The scheduled run had found nothing due (W1 is due 2026-11-28), and Thomas then
asked for three things: build the model architecture, research how to optimise it, and write the system design.*

**Decision this can change:** which predictive distribution the drug-effect layer ships. If O1a passes PREREG-O1 →
ship it as v3.1. If it fails or is inconclusive → keep M1 and close the interval question for open data.

## What was built

The package `lumen/` (commit 91876dc, branch `rounds-70-72`) can be installed with `pip install -e .`.

- **L0 provenance.** An md5 table for the forecaster and its four input tables. `lumen verify` reports all of them
  unchanged.
- **L1 measurement contracts.** `ArmCount` and `Observable` enforce the round-58 rules: scale, dispersion, per-visit
  n, population and missing codes. Run over the three committed tables, they found 0 violations in LN and SRI-4.
  In BICLA they flagged 2 vobarilizumab arms whose week-24 and week-48 counts are identical (control 28/62 at both
  visits). Shipped outputs use one row per arm, so this cannot change them; it is queued for a source check
  (D-BICLA).
- **L2 baselines and gate.** Null, pooled and carry-forward, scored by IV-MSE, log score and coverage. `gate.verdict`
  implements rule 2b as code. On RT-1's numbers it returns FAIL, as round 127 did.
- **L3 drug-effect layer.** It imports `experiments/lumen-bench/lumen_forecast.py` from its md5-frozen file and never
  copies it. Its output equals the legacy CLI's exactly (tested). With `--early` it still prints no headline.
- **L4 registry.** Eight modules are listed, each with its status, its label with n, and its re-open condition.
  **One is shipped:** the pooled prior. It cannot mark a module shipped without an n.
- **Tests.** 17 pass in a clean checkout. They reproduce the r90 and W1-SLE frozen forecasts exactly, and the legacy
  CLI's output (except for the early-read-out note, whose text is abridged in v3).

## What the optimisation research found

Details are in `docs/RESEARCH-optimising-LUMEN-2026-11-27.md`, drawn from 9 PubMed papers.

- **Learned models.** Simple linear baselines still match or beat deep perturbation models
  (Ahlmann-Eltze 2025, [doi:10.1038/s41592-025-02772-6](https://doi.org/10.1038/s41592-025-02772-6)). This agrees with
  VC-1 and rounds 92-104. **Established in the field; here Emerging, n 11.** No learned module is retried until open
  data reach ≥ 30 longitudinal patients.
- **Where gains can come from.** The shipped interval combines a fixed-effect mean, a DerSimonian-Laird tau² and a z
  quantile. All three are known to under-cover with few trials
  (Seide 2019, [doi:10.1186/s12874-018-0618-3](https://doi.org/10.1186/s12874-018-0618-3);
  Partlett & Riley 2017, [doi:10.1002/sim.7140](https://doi.org/10.1002/sim.7140)).
- **Synthetic check, not real data** (300 replicates per cell). With a true tau of 0.3-0.6, M1 covered 0.84-0.94 and
  REML + random-effects weights + a t quantile covered 0.94-0.99, against 0.95 nominal. **Hypothetical for LUMEN**:
  at its pool sizes (≥ 20 trials) the expected gain is only 1-3 coverage points.
- **Next candidates.** A robust MAP control-rate prior (O2) could let LUMEN print absolute rates honestly. Phase-2
  shrinkage (O3) needs a rule-8 check first.

## Pre-registered

**PREREG-O1** (4e8975a, pushed 04:15Z). O1a is compared with the shipped M1 and with a null of the same width, on
RT-1's time-forward units. The primary metric is the log score per trial. FAIL fires at 0.05 nats, or if IV-MSE is
10% worse than a baseline. The test is **12b-timed (committed data)**: every table has been analysed before. The
analysis runs ≥ 1 h later in a scheduled run. That run also freezes O1a/O1b forecasts for the 8 W1-SLE trials, as a
prospective record.

## Documents

- `docs/SYSTEM-DESIGN-v3.md` (also `claude/SYSTEM-DESIGN-v3.md` in the project) supersedes v2. Its new constraints
  C14-C18 are: no mechanism input, no early headline, relative effects only, a data-size entry condition for learned
  modules, and frozen maths imported rather than copied. It also specifies the module lifecycle and how the
  world-model slot would be tested.
- Comparisons this round: 0 on real outcomes; 1 synthetic simulation (4 interval methods × 9 cells). No labels changed.

## Was the request answered honestly?

Yes, within the project's rules. The architecture is the system that actually works plus slots that say why each
one is empty. It is not a new JEPA model trained on n = 11, which round 142 already showed fails. The system design
names the condition under which a world-model module could enter.

**Sources:** PubMed — the DOIs above plus Turner 2015 [doi:10.1002/sim.6381](https://doi.org/10.1002/sim.6381),
Röver 2021 [doi:10.1002/jrsm.1475](https://doi.org/10.1002/jrsm.1475), Lilienthal 2024
[doi:10.1002/jrsm.1685](https://doi.org/10.1002/jrsm.1685), Schmidli 2014
[doi:10.1111/biom.12242](https://doi.org/10.1111/biom.12242), Goteti 2024
[doi:10.1002/psp4.13083](https://doi.org/10.1002/psp4.13083). Repo commits 91876dc, 4e8975a, 04be263.

*Research and simulation only. Not clinical or drug-development guidance.*
