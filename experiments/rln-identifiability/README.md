# RLN identifiability experiments (rounds 40–44)

Mechanistic model of refractory lupus nephritis used **not to predict patients** but to
establish what measurement can and cannot identify. Every number in `docs/rounds/` is
reproduced by the scripts here.

## Files

| file | what it does |
|---|---|
| `model.py` | shared six-state model, observation operator, sensitivity helper. No import side effects. |
| `r40_feedback.py` | the structural zero, and whether real ODE feedback restores identifiability of the damage split |
| `r44_squeeze.py` | virtual population: placebo falsifier (r43) and the variability–power squeeze (r44) |
| `r44_analyze.py` | multi-seed analysis, incl. the threshold-artifact discriminator |
| `results/*.json` | raw outputs, checkpointed per cell |

Run from this directory: `python3 r40_feedback.py`. (Not `python3 -I` — isolated mode
drops the script directory from `sys.path`.)

## The model

Six states. `Pl` (long-lived plasma cells) is **CD19-negative**, so CD19-directed therapy
cannot reach it; `D` (renal damage) is **monotone**. Both are load-bearing — see below.

Depletion is **off by default** in `make_params`: the core carries no therapy, and
interventions belong in modules. `r40_feedback.py` switches it on (`dep_kill`) because
that analysis *is* a depletion study.

## Results these scripts reproduce

- **The structural zero.** At `phi=0` the apparent molecular sensitivity to `k_dam` scales
  as exactly `1/h` across five step sizes — the signature of a true zero, not a small
  effect. Reproduces round 40 independently.
- **Feedback does not rescue it at realistic strength.** Fisher marginal `sigma` on
  log `k_dam` improves only 1.01–1.23× at `phi=0` and ~2× at `phi=20`. Meanwhile
  **severity is worth 437×**. The information lives in the clinical channel.
- **The LLPC floor.** ~89% B-cell and 97% short-lived-plasma-cell depletion buys a **7.9%**
  autoantibody fall, because ~99.7% of production comes from the protected CD19⁻ pool.
  A B-cell count is a misleading surrogate.
- **The variability–power squeeze.** Cohen's d falls 0.488 → 0.234 as between-patient CV
  rises 0.40 → 1.00 (p = 5.4e-6). The *continuous* endpoint erodes too, so this is not a
  binary-threshold artifact.

## Two instruments, because one was not enough

The cosine-VIF is **not monotone in the observable set** — it reported gains below 1.0,
i.e. adding a channel appearing to destroy information, which is impossible. Fisher
marginal variance is monotone by construction and governs every verdict. The VIF code is
retained only to reproduce the structural zero.

## Numbers that must not be quoted about real trials

`treated()` is a **stylised** voclosporin-like arm, not fitted to AURORA. Its effect
(+9.8 pp at CV 0.70) is weaker than voclosporin's observed (+18 pp), so the absolute
required-N figures are a property of this drug model. The 4.5× escalation **with CV** is
the real result; the absolute N is not.

## Known limitations

No Tregs, no germinal centre, no T-cell help, no interferon. Synthetic throughout, one
operating point. The multiplicative-reshaping explanation for the squeeze was offered
after seeing the result and is **not separately tested** (additive vs multiplicative
jitter is the obvious falsifier, not run).

*Research and simulation only. Not clinical or drug-development guidance.*
