# Round 51 — my scale-mismatch argument was wrong; the real reason is structural (2026-10-29)

Computed. Code: `experiments/rln-identifiability/r51_scale_mismatch.py`.

**Question:** can any model in the NVIDIA drug-discovery catalogue improve this lupus
model? I had already answered "no, scale mismatch — the stack works at molecule scale and
our unknowns are population rate constants." **That argument is false, and the test
falsified it.** The answer is still no, for a stronger reason.

---

## 1. The argument I gave, and its falsification

Claim: the failing quantity (round 50's week 24→48 eGFR trend, understated ~4×) depends on
kidney-block rate constants, not on anything a molecular-scale model could ground, so
perfect molecular knowledge buys nothing.

Relative sensitivity of that quantity, |d log(trend)/d log(param)|, by group:

| group | max sensitivity |
|---|---|
| humoral / molecular | **5.65** (`d_Ab`), then `d_Pl` 4.93, `k_a` 4.47, `g_l` 2.99 |
| kidney block | **6.72** (`c_inf`), then `k_prog` 4.61, `k_tox` 2.25 |
| observation / pharmacology | 4.96 (`s_inf`) |

**Ratio kidney/humoral = 1.2×, against the 10× the argument needed.** Humoral parameters
are essentially as influential as kidney ones. The scale-mismatch argument fails.

## 2. Why it failed — round 40's zero is DIRECTIONAL, and I conflated the directions

Round 40 established that kidney *parameters* do not move molecular *states* under
feed-forward coupling, so an upstream observable carries exactly zero information about a
downstream degeneracy. **That is one direction, and it is not symmetric.**

In V2, activity is driven by autoantibody (`dA/dt = k_a·Ab − k_clear·A`) and the damage
driver depends on activity. So humoral **parameters** propagate strongly downstream even
though humoral **observables** are worthless for resolving the downstream damage split.

> **Upstream observables cannot identify downstream degeneracies. Upstream parameters
> strongly influence downstream predictions. These are different statements, and treating
> the first as licensing the second is an error.**

Worth recording as a general trap: this project's strongest result is easy to over-apply.

## 3. The real reason the answer is still no: the defect is structural

If the failing quantity is sensitive to `d_Ab`, could grounding `d_Ab` close the gap?
Sweeping it across 64×:

| `d_Ab` × | wk 24→48 eGFR | reaches observed −3.30? |
|---|---|---|
| 0.25 | −1.324 | no |
| 0.50 | **+0.940** | no |
| 1.00 | −0.295 | no |
| 2.00 | −1.030 | no |
| 4.00 | −1.311 | no |
| 8.00 | −1.435 | no |
| 16.00 | −1.495 | no |

**The trend saturates near −1.5 and never reaches −3.30.** So although the quantity is
sensitive to the most influential humoral parameter, **no value of it closes a 4× gap.**

> **The round-50 defect is structural, not a parameter-grounding problem.** Structure is
> fixed by changing equations, not by improving parameter estimates — and every model in
> the catalogue grounds parameters, not structure. That is the correct reason the answer is
> no, and it is stronger than the one I gave because it does not depend on a scale
> argument that turns out to be false.

## 4. An anomaly, flagged not tidied

The `d_Ab` response is **non-monotone**: 0.25× gives −1.324, 0.50× gives **+0.940**
(filtration *rising*), 1.00× gives −0.295. A monotone parameter should not reverse the
sign of a downstream trend twice. **Unexplained.** Project rule is to follow the anomaly
rather than tidy it; it is not followed here and is logged as an open item. It may be a
competition between the inflammatory eGFR factor (`c_inf·A`, reversible) and nephron loss
(`N`, monotone), which move filtration in opposite directions on different timescales —
that is a hypothesis, untested.

## 5. The one honest exception, and where the stack could contribute

**`d_Ab` is the single humoral parameter where a structure/affinity model has a real
route.** Autoantibody clearance is set by FcRn-mediated IgG recycling, which is a
protein–protein affinity problem of exactly the kind OpenFold3/Boltz-2 address. It is also
the most influential humoral parameter on the failing quantity. So the route is real — but
round 41's independent benchmarks put Boltz-2 affinity at r ≈ 0.42 with compression, and
§3 shows it cannot close this particular gap anyway.

**Where the stack could genuinely contribute to the project**, as opposed to this model:
`lessons.md` records that specificity is receptor-level, that an anti-dsDNA B cell and a
tetanus B cell are identical to any CAR, and that "no model invents a targeting handle."
CAAR-T — autoantigen as receptor — attacks exactly that, and RFdiffusion, ProteinMPNN and
RFD4-Proteina are binder-design tools for precisely this problem. **That extends the
project sideways into intervention design; it does not improve the dynamics model.**

## 6. What this does not establish

- One operating point, one variant (round 49's V2), one failing quantity. Sensitivities are
  local and linearised.
- §3 sweeps `d_Ab` alone. A joint move of several humoral parameters might reach −3.30; not
  tested, and if it did, the conclusion in §3 would need revisiting.
- "Structural, not parameter" is established for `d_Ab` and inferred for the rest from the
  sensitivity ranking. Not proven for every parameter.
- §5's FcRn route is a mechanism argument, not a test. Nothing was run against any
  structure or affinity model.

## 7. Next

1. **Follow the §4 anomaly** — a sign reversal in a monotone parameter is the kind of thing
   this project has twice turned into its biggest gains.
2. Round 50 §6 still stands as the top item: locate the crossover horizon.
3. If intervention design is ever in scope, §5 names the one place the catalogue is the
   right tool — and it is a different project from this one.

*Research and simulation only. Not clinical or drug-development guidance.*
