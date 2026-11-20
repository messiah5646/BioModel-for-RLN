# SLE / lupus nephritis in-silico trial harness

A setup for testing drug **mechanisms** on a virtual lupus-nephritis trial. You describe a drug as what it does
and how strongly; the harness simulates a control arm and a drug arm, compares them over up to three years, and
tells you whether a real trial of a given size could even see the difference.

**Research and simulation only. Not clinical or drug-development guidance.** It cannot tell you whether a drug
works in patients. It tells you what a mechanism of a chosen strength would do *if this model were right*, and how
much of that answer the real data leave undetermined.

## Quick start

```bash
cd experiments/sle-virtual-trial
python3 vtrial.py drugs/sglt2i_like.json                 # one drug
python3 vtrial.py drugs/*.json --weeks 156 --n 180       # several; 3-year follow-up; 180 patients per arm
python3 backtest.py                                       # the validation gate (run after vtrial.py once)
```

The first run builds the control-arm ensemble (about 2–3 minutes) and caches it in `runs/`. Each drug writes
`runs/<drug>/result.json` and `runs/<drug>/trajectories.png`.

## Write your own drug

```json
{"name": "my drug", "effects": [
  {"mechanism": "immunosuppression", "strength": 0.2},
  {"mechanism": "bcell_depletion", "strength": 0.25, "start_day": 0, "end_day": 28}
]}
```

| mechanism | biology it stands for | what it changes in the model |
|---|---|---|
| `immunosuppression` | sustained dampening of autoimmune drive (MMF-, CNI-, JAK-like) | lowers the inflammatory set-point `s_inf` by `strength` (0–1) |
| `antiproteinuric` | lowers protein leak without touching immunity (SGLT2i-, RAAS-like) | `leak_red = strength` |
| `haemodynamic` | reversible fall in filtration from vascular tone (a CNI side effect) | `h_drug = strength` |
| `bcell_depletion` | kills B cells and short-lived plasma cells for a window (anti-CD19/CD20, CAR-T-like) | kill rate per day over `start_day`–`end_day` |
| `activity_clearance` | speeds resolution of kidney inflammation directly | multiplies `k_clear` by `1 + strength` |

**Not representable** (the harness refuses them by name rather than guess): long-lived plasma-cell depletion
(anti-BCMA), interferon blockade (anifrolumab), regulatory T-cell expansion (low-dose IL-2), co-stimulation
blockade (abatacept, anti-CD40). The model has no state for these.

## How it works, in six lessons

**1. The disease model.** Six quantities change over time: B cells, short-lived and long-lived plasma cells,
autoantibody, kidney inflammation (activity) and irreversible kidney damage. What a trial measures, proteinuria
(UPCR) and filtration (eGFR), is computed from activity and damage. This is `model.py` from the project's
identifiability work.

**2. The control arm is fitted to real patients.** The placebo arm of AURORA 1 and 2 (voclosporin trial; both
arms got mycophenolate plus steroids) is the calibration target. It includes the *day-0 step* found in rounds 77–78:
both arms change abruptly within two weeks of starting background therapy, which the slow immune chain in the
model cannot produce on its own.

**3. An ensemble, not one answer.** Real data leave some parameters undetermined. The biggest is `k_tox`, how much
proteinuria itself damages the kidney: rounds 75, 76 and 79 found a ceiling and no floor. So the harness fits 40
control arms (5 values of `k_tox` across its allowed range, each with 8 noise-jittered copies of the data) and runs
every drug on all 40. The bands in the plots are that spread. If a conclusion changes across the `k_tox` range,
the summary says `verdict_depends_on_unidentified_k_tox: true`, which means: the data cannot settle this one.

**4. Strength is your input.** The model cannot tell you how strong a new drug is. That number has to come from
somewhere else (phase 1 biomarkers, a related drug, or a sweep over values).

**5. Power: could a real trial see it?** The predicted difference at week 52 is compared with the real
between-patient spread in AURORA (standard deviation of change), giving the chance that a two-arm trial of `--n`
patients per arm detects it at the usual 5% level. At zero effect this prints 0.025 by construction of the
formula, not 0.05.

**6. Backtest before you trust it.** `backtest.py` runs the voclosporin-like drug and compares it with what
AURORA actually found. Read its label: the voclosporin strengths were fitted to AURORA year 1 in round 78, so the
year-1 comparison is in-sample. The treated arm's 3-year eGFR slope was not used anywhere and is the real test.

## Current backtest result (2026-11-20)

| check | model | real | verdict |
|---|---|---|---|
| placebo add-on (negative control) | difference 0.00 everywhere | must be 0 | pass |
| week-52 UPCR difference, voclosporin vs placebo | −0.94 | −0.77 | in-sample, close |
| week-52 eGFR difference | −2.69 | −3.0 | in-sample, close |
| **treated-arm eGFR slope, months 12–36** | **−2.37 mL/min/yr** | **−1.00 ± 0.87** | **held out: z = −1.6, inside 2 SE but too steep** |
| placebo eGFR slope, months 12–36 | −3.29 | −3.17 ± 0.93 | in-sample (used in calibration) |

The held-out miss is in the same direction as round 79: the model under-produces the long-term protection seen
in the treated arm. Treat any 3-year eGFR projection for a drug as probably too pessimistic about protection.

## Results for the example drugs

Differences are drug minus control. Medians over the ensemble, range in brackets. Power assumes 180 per arm.

| drug | wk-52 UPCR | power | wk-52 eGFR | power | wk-156 eGFR | depends on `k_tox`? |
|---|---|---|---|---|---|---|
| voclosporin-like (backtest) | −0.94 [−1.18, −0.81] | 0.83 | −2.69 [−3.04, −1.86] | 0.40 | −0.94 [−1.65, +0.61] | **yes** |
| SGLT2i-like (illustrative strengths) | −0.93 [−1.13, −0.82] | 0.83 | −2.16 [−2.29, −1.34] | 0.28 | −1.85 [−2.19, +0.26] | **yes** |
| CD19 depletion, 4 weeks (illustrative) | −0.10 [−0.12, −0.09] | 0.05 | +0.88 [+0.66, +1.03] | 0.08 | +0.66 [+0.36, +0.86] | no |
| +30% immunosuppression (illustrative) | −0.89 [−1.21, −0.72] | 0.79 | +3.98 [+3.50, +4.90] | 0.71 | +7.81 [+6.86, +8.81] | no |

What to learn from this table:

- **Anti-proteinuric drugs are exactly where the model is blind.** Whether lowering proteinuria protects the kidney
  long-term depends on `k_tox`, which no available trial determines. The harness says so instead of picking an answer.
- **A one-off CD19 depletion barely moves the kidney in this model.** That is a property of the model's structure
  (it protects long-lived plasma cells, and autoantibody from them keeps driving activity), and that structure is
  contested (project `DATA-CATALOGUE.md` §5). Read it as "this model's assumption", not as evidence about CAR-T.
- **Sustained immunosuppression is the mechanism the model rewards most**, because activity drives damage in every
  ensemble member. That is again the model's structure speaking.

## What this setup is not

- Not a patient simulator: arm means only, no responder rates, no individual variation.
- Not general SLE: kidney only, calibrated on one trial's background therapy.
- Not validated predictively: the underlying model does not beat simple forecasting rules on real data
  (project rounds 50–53, 78, 79).
- Not a drug ranking. The ranking in the table above follows directly from assumptions you can read in `model.py`.

## Files

- `vtrial.py`: the harness (ensemble, drug mapping, simulation, power, report)
- `backtest.py`: validation against AURORA
- `drugs/*.json`: example drug specs
- `runs/`: cached ensemble and outputs

## Accuracy on a real trial the harness never saw (round 82, 2026-11-20)

AURA-LV (voclosporin 23.7 mg vs placebo, earlier trial, different patients), UPCR change at weeks 24 and 48.
Scored as pre-registered (project doc `prereg-round82.md`); `external_aura_lv.py`, `external_aura_lv_attack.py`.

| predictor | error on 4 arm values (RMSE, mg/mg) | treatment difference, wk 24 / 48 (real: −1.55 / −1.61) |
|---|---|---|
| **harness, as pre-registered** (baseline set to AURA-LV's) | **0.64** | −1.32 / −1.39 |
| copy AURORA's numbers | 0.93 | −1.15 / −0.91 |
| **copy AURORA's numbers, scaled by baseline** | **0.44** | −1.59 / −1.32 |
| harness + extra step rescale (NOT pre-registered) | 0.37 | −1.54 / −1.60 |

The harness loses to a zero-parameter rule that just rescales the earlier trial's results. A variant that wins
exists, but it uses a change I added in code without writing it down first, so it does not count. Typical error
per value is about 0.4–0.6 mg/mg, against measurement uncertainty of about 0.4 mg/mg per value.
