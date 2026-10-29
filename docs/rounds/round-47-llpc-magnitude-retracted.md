# Round 47 — the LLPC headline's magnitude is an assumption, not a result (2026-10-29)

Code: `experiments/rln-identifiability/r47_llpc_share.py`. Computed.

**Partial retraction of this project's most quotable number.** The structural claim
survives; the figure attached to it does not.

---

## The claim under attack

Rounds 40/44 and `docs/system-design.md` §4.2 report: *"~89% B-cell and 97%
short-lived-plasma-cell depletion buys a 7.9% autoantibody fall, because ~99.7% of
production comes from the protected CD19⁻ pool."*

It is the most clinically striking thing the project has produced, and the most likely to
be quoted. So it is the one most worth attacking.

**`model.py` sets `g_s = 0.4·d_Ab`, `g_l = 0.6·d_Ab`** — 60% of autoantibody production
from the long-lived pool at steady state. That split is **tier-C, assumed, never measured.**

## The real-world constraint I had not used

`lessons.md` records the signature dissociation seen in CD19 CAR-T for lupus: **anti-dsDNA
seroconverts while vaccine IgG holds.** Vaccine IgG is the classic long-lived-plasma-cell
product. So if autoantibody disappears while LLPC-derived protective antibody persists,
the **autoreactive** LLPC share must be *low* — not 60%. Our parameterisation is plausibly
backwards for the autoreactive compartment specifically.

## Result

Autoantibody fall at day 240 under CD19 depletion, sweeping the assumed share:

| `g_l` share | Ab fall @ d240 | production from protected pool | B-cell fall |
|---|---|---|---|
| 0.00 | **75.7%** | 0.0% | 89.6% |
| 0.05 | 55.2% | 93.0% | 89.6% |
| 0.10 | 42.5% | 96.6% | 89.6% |
| 0.20 | 27.6% | 98.4% | 89.6% |
| 0.40 | 13.7% | 99.4% | 89.6% |
| **0.60 (ours)** | **7.1%** | **99.7%** | 89.6% |
| 0.95 | **1.2%** | 100.0% | 89.6% |

**The autoantibody fall spans 1.2% to 75.7% as an assumed parameter moves.** Our 7.9%
figure is one point on that curve, and the curve covers nearly the whole range.

## What survives and what does not

**RETRACTED: the magnitude.** "Deep B-cell depletion buys only a ~8% autoantibody fall" is
not supported. It is an assumption wearing a result's clothes, and the real-world
seroconversion constraint argues the true value is at the *other* end — a large fall.

**SURVIVES, and more robustly than expected: the structural claim.** "After depletion the
protected pool dominates what remains of production" holds for *any* share above ~5%:
93.0% at `g_l` = 0.05, 96.6% at 0.10, ≥99% above 0.40. It is near-saturated because
depletion crushes the CD19⁺ sources by ~97% regardless, so even a small protected
contribution dominates the remainder.

**And the corrected version is a better story.** At a low autoreactive LLPC share, *both*
hold at once: autoantibody falls substantially **and** the residual is almost entirely
from the pool CD19 therapy cannot reach. That is deep short-term response followed by
relapse from an untouched reservoir — which matches the observed CAR-T lupus pattern far
better than our original "barely moves" claim, and it preserves the clinical point that a
**B-cell count is a misleading surrogate** while fixing the arithmetic behind it.

## Consequence

- `docs/system-design.md` §4.2 and `where-this-stands` must not quote 7.9%.
- `g_l`/`g_s` is now the model's highest-value unmeasured parameter: it swings the headline
  across its entire range. It should be **fitted to paired anti-dsDNA and vaccine-IgG
  titres**, which `lessons.md` already identified as the cheapest decisive measurement.
- The identifiability results of rounds 40/44 are **unaffected** — they concern
  `k_dam`/`k_tox` in the kidney block, which the humoral split does not touch (that is
  round 40's structural zero, in our favour for once).

## What this does not establish

- Single-parameter sweep at one operating point; no joint variation with `d_Pl`, which also
  sets the pool's persistence.
- The seroconversion constraint is taken from our own notes on the CAR-T literature and was
  **not re-verified at source this round**. It is the load-bearing external fact here, so it
  should be.
- No fit was performed. This shows the headline is sensitive to an assumption; it does not
  establish the correct value.

*Research and simulation only. Not clinical or drug-development guidance.*
