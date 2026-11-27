# LUMEN — System Design v3.0

**Status: canonical architecture spec, 2026-11-27 (round 148). Supersedes v2.0 (2026-11-02, `claude/SYSTEM-DESIGN.md`).**
Code: package `lumen/` on branch `rounds-70-72` (commit 91876dc). Evidence labels follow `claude/where-this-stands.md`
(round 145). This document adds no new empirical claims.

## 0. The one-paragraph version

LUMEN is a **gated forecasting system**. Modules can be added, but none produces headline output until it has beaten
trivial baselines on a pre-registered held-out test. Today exactly one module ships: the **drug-effect pooled prior**.
It gives a relative effect (log odds ratio) for a lupus or LN trial, with an interval. Every other layer (early
read-outs, control rates, patient state, cell state, flare, mechanism) exists as a registry entry with its verdict,
its n and the one thing that would re-open it. The JEPA / world-model aim is kept as a **slot with an entry condition**,
not as an architecture that is assumed to work.

### 0.1 What changed from v2, and why

| v2 (2026-11-02) | v3 | Reason |
|---|---|---|
| Purpose: identifiability instrument that can also simulate | Purpose: **trial forecaster with a module gate**; the mechanistic core is a supporting module | Rounds 85-145: only evidence pooling survived held-out tests |
| Layers L0-L5 around an ODE core | Layers L0-L5 around a **registry + gate** (§2) | Six mechanistic harnesses lost to 0-1 parameter rules |
| "Baselines are a first-class component" (principle) | Baselines + rule-2b gate are **code** (`lumen.baselines`, `lumen.gate`) | Principle was applied inconsistently until round 126 |
| Measurement contracts specified | Contracts **implemented** (`lumen.contracts`); they flag 2 identical-count rows in the BICLA table | Round 58 near-misses |
| Learned surrogates "last, only where they beat baselines" | `cell_state.vc1_learned` = **failed** entry; retry needs >= 30 longitudinal patients | VC-1, round 142 |
| No prospective layer | **L5 prospective record**: frozen forecasts with md5-pinned code | Rounds 94, 147 |

## 1. Design constraints (findings the architecture must not contradict)

C1-C13 of v2 stand. They cover the structural zero, baselines, the action entering the transition, prediction vs
parameter claims, the CD19⁻ reservoir as a variant, Tregs, the placebo arm, growable state, the missing
molecular-to-clinical bridge, contested reservoir identity, provenance, scale contracts, and structure buying
identifiability. New in v3:

| # | Finding | Consequence |
|---|---|---|
| C14 | The pooled mean beats molecule/target anchoring and mechanism class (r86, r87, r116; Supported, n 23 trials) | The drug-effect layer has **no mechanism input**; mechanism only decides exchangeability |
| C15 | No early read-out rule is validated; the question is closed (rule 8) | With an early read-out, the system prints **no headline**, just three models side by side |
| C16 | Control rates are not predictable from year or definition (RT-2; n 60 trials) | Output is relative effects only, until a gated control-rate module passes |
| C17 | Learned cell models lose to population averages at n ≈ 11 (VC-1) and in the field (Ahlmann-Eltze 2025) | Learned modules need a **data-size entry condition** before a test is even pre-registered |
| C18 | Frozen forecasts must reproduce byte-for-byte for audits | Shipped maths is **imported from an md5-frozen file, never copied**; any change is a new module |

## 2. Layers

```
  L5  PROSPECTIVE RECORD    frozen forecasts (r90, W1 LN x13, W1-SLE x8), md5 of code + tables, scored when results post
       ^
  L4  REGISTRY              every module: status, evidence label + n, consumes/produces, re-open condition
       ^                    only status 'shipped' may produce a headline
  L3  LAYERS                drug_effect (shipped) | control_rate | patient_state | cell_state | flare | mechanistic
       ^
  L2  BASELINES + GATE      M0 null, M1 pooled, C0 carry-forward; IV-MSE, log score, coverage; rule-2b verdict
       ^
  L1  MEASUREMENT           ArmCount (binary endpoints), Observable (scale, dispersion, per-visit n, population)
       ^
  L0  PROVENANCE            FROZEN md5 table, fetch.sh manifests, run stamps (run_utc, git_head, PREREG commit)
```

| Layer | Module | File |
|---|---|---|
| L0 | `provenance.FROZEN`, `verify`, `stamp` | `lumen/provenance.py` (+ `scripts/provenance/fetch.sh`) |
| L1 | `ArmCount`, `Observable`, `combinable`, `check_visits` | `lumen/contracts.py` |
| L2 | `null`, `carry_forward`, `iv_mse`, `log_score`, `coverage`; `gate.verdict` | `lumen/baselines.py`, `lumen/gate.py` |
| L3 | `drug_effect.forecast / prior` (wraps `experiments/lumen-bench/lumen_forecast.py`, md5 2c11b9a5) | `lumen/layers/drug_effect.py` |
| L3 (unshipped) | heterogeneity and robust-prior candidates | `lumen/candidates.py` |
| L4 | `REGISTRY`, `shipped()` | `lumen/registry.py` |
| L5 | frozen JSONs + freeze scripts | `experiments/lumen-bench/prospective/` |
| API | `lumen forecast | status | verify` | `lumen/cli.py` |

## 3. The module lifecycle (the core of v3)

```
 idea ──► candidate code (lumen/candidates.py or experiments/)   synthetic tests only; no real outcomes
      ──► decision line (rule 7)  ── none? ──► [parked: no decision at stake]
      ──► PREREG pushed  ──(>= 1 h, rule 12b)──►  fetch / analysis run  ──► gate.verdict
             PASS ──► registry status 'shipped' (label + n), new file frozen with its own md5; old outputs kept
             FAIL ──► registry status 'failed' + re-open condition; added to rule 8 if asked >= 2 times
             UNDERPOWERED / INCONCLUSIVE ──► stays a candidate; no re-test without new data
```

Invariants the code enforces:
1. `registry.Module` refuses status `shipped` without an n in its evidence string (rule 12).
2. `drug_effect.core()` refuses to load if `lumen_forecast.py` has drifted from its md5 (C18).
3. `gate.verdict` returns FAIL when a candidate loses to **any** baseline by >= the margin (rule 2b), so a failure
   cannot be reported as inconclusive.
4. Tests (`python -m pytest`) reproduce the r90 and W1-SLE frozen forecasts and the legacy CLI exactly.

## 4. Drug-effect layer (shipped)

Input: arm counts `(r_drug, n_drug, r_ctrl, n_ctrl)` at the latest visit ≤ 104 weeks, from one of three tables
(LN complete renal response, 20 trials; SLE SRI-4; SLE BICLA). Log OR uses a 0.5 continuity correction.

Output, without an early read-out: pooled log OR, SD = sqrt(DL tau² + Var(mean)), OR with a 95% interval, P(benefit).
With an early read-out: `headline = NONE`, and M0 / M1 / C0 are printed side by side. Control rate: "pooled, not
modelled".

Honest limit: a never-tested target gets the historical base rate, about OR 1.4 (0.75-2.5).

**Known weakness → O1.** The centre is a fixed-effect mean and tau² is DerSimonian-Laird. Both are known to
under-cover with few trials. In synthetic data, coverage was 0.84-0.94 against 0.95 nominal when tau ≥ 0.3. PREREG-O1
(4e8975a) tests REML + random-effects weights + a t quantile against M1 and M0 on the RT-1 time-forward units.

## 5. Gated slots (not shipped)

| Slot | Status | Entry / re-open condition |
|---|---|---|
| `drug_effect.early_update` | option-emerging, settled | W1: ≥ 5 frozen trials reporting both visits |
| `control_rate` | not implemented | a pre-registered model beats the pooled rate on coverage + log score (candidate O2) |
| `patient_state.ifn_setpoint` | standalone, no consumer | a module that needs IFN(t) before day 90 |
| `cell_state.*` (JEPA / virtual cell) | failed (VC-1) | open longitudinal single-cell data, ≥ 30 patients under one therapy |
| `flare` | clinical-only | a new data type (urine, kidney) on open data |
| `mechanistic_core` | supporting | beats persistence out-of-sample on levels |

### 5.1 The world-model slot, specified so it can be tested

If the data condition is met, a cell-state module must implement
`predict(baseline_state, therapy, horizon) -> state` together with a predictive spread. It is scored per patient
against **three** baselines: the patient's own baseline (carry-forward), the population average at the target visit,
and the others' average plus a mean therapy shift. The primary metric is held-out per-patient error. The module ships
only if it beats all three by ≥ 10%. VC-1's harness (`experiments/virtual-cell/run_vc1.py`) is the reference
implementation of this test.

## 6. Data flow for one forecast

```
lumen forecast MOL TARGET --table T
  -> provenance.verify(lumen_forecast.py, table T)        md5 or refuse
  -> lumen_forecast.load(T)                               log OR, SE per arm-visit
  -> lumen_forecast.forecast(anchor=False)                pooled prior (M1)
  -> output policy (C15, C16)                             headline / no headline, control rate refused
  -> JSON with a 'layer' block naming status and evidence
```

## 7. Build order from here

1. **W1** prospective scoring (weekly, from 2026-11-28): the only test that can validate M1 out of time.
2. **O1** analysis run (≥ 1 h after 4e8975a). PASS → `drug_effect_v31.py` frozen with its own md5; the registry
   gains a second shipped entry and M1 moves to `superseded`.
3. **O2** robust MAP control rate, only after O1, because it reuses the heterogeneity estimate.
4. Wire `lumen.contracts.check_visits` into the table builders as a warning: two BICLA rows (vobarilizumab
   wk 24 = wk 48) need a source check. This does not change shipped outputs, which use one row per arm.
5. Learned / world-model retry only when §5.1's data condition is met.

## 8. Known holes

- No outside reader (clinician, trial statistician, pharmacometrician). This is still the highest-value gap.
- Thin prospective record: 0 of the 21 W1 frozen forecasts (13 LN, 8 SLE) scored so far; r90 was scored on n 2 trials.
- Tregs, germinal centre and IFN are absent from the mechanistic core. They matter for biology, not for current
  forecasts.
- The method is model-based meta-analysis, not a new method (Goteti 2024). What is ours: the open LN table, the frozen
  prospective record, the negative transfer results, and the gate.

*Research and simulation only. Not clinical or drug-development guidance.*
