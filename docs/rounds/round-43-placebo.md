# Round 43 — the placebo falsifier: survived, but not on the terms I pre-registered (2026-10-28)

Pre-registration: `out/PREREG_round43.md`, written before any test number existed.
Code: `src/t6_placebo.py`, `src/t6b_pop.py`, `src/t6c_rest.py`. Seed 20261028.

**H1 passed, H2 passed, H3 FALSIFIED.** The architecture *can* produce a realistic
control-arm response rate and a realistic treatment separation — but only at a
between-patient variability larger than I declared plausible in advance. The headline is
therefore conditional, and the pre-registration is what makes that visible.

**The round also produced an unpredicted structural result (§4) that is more useful than
the falsifier itself.**

---

## 1. What was actually tested

`lessons.md` §8 names the falsifier: real placebo arms in lupus trials improve 33–44%, and
"a model with a stable disease attractor cannot reproduce that."

Sharpening it changed the test. A **deterministic** model gives a response *rate* of 0% or
100%; 33–44% is a statement about a **population**. Also, this model does not have a stable
attractor under control conditions — `supp(t)` encodes the background MMF/steroids a lupus
control arm actually receives. So the question became:

> Does a plausible between-patient parameter distribution produce a control-arm complete
> renal response rate in the target band, **while** keeping the treated arm higher?

Endpoint: AURORA-1-like CRR at week 52 — **UPCR ≤ 0.5 AND (eGFR ≥ 60 or no >20% fall)**.
Jittered parameters chosen **by sensitivity, not alphabetically** (lessons §7, the R24
error): `p0`, `d_Ab`, `k_clear`, `k_a`, `g_l`, `d_Pl`. Note `g_l` — long-lived plasma cell
antibody output — ranks 5th of 19, independent support for §4.2's compartment.

## 2. Result

n = 1200 per arm per cell, log-normal jitter, binomial 95% CI.

| CV | control CRR | treated CRR | diff (pp) | H1 | H2 |
|---|---|---|---|---|---|
| 0.10 | 0.0% [0.0, 0.0] | 0.2% [0.0, 0.4] | +0.2 | n | n |
| 0.20 | 0.8% [0.3, 1.2] | 6.8% [5.4, 8.3] | +6.1 | n | n |
| 0.30 | 4.4% [3.3, 5.6] | 17.8% [15.7, 20.0] | +13.4 | n | **Y** |
| 0.40 | 10.4% [8.6, 12.1] | 24.7% [22.3, 27.2] | +14.4 | n | **Y** |
| 0.50 | 16.8% [14.7, 18.9] | 30.0% [27.4, 32.6] | +13.3 | n | **Y** |
| **0.70** | **25.3% [22.8, 27.8]** | **36.9% [34.1, 39.6]** | **+11.6** | **Y** | **Y** |
| 1.00 | 33.9% [31.1, 36.7] | 42.9% [40.0, 45.8] | +8.9 | **Y** | n |

**At CV = 0.70 the model produces control 25.3% and treated 36.9%.** For reference,
AURORA-1 reported roughly 23% and 41%. That is close, and it was not tuned to those
numbers — the jitter set and thresholds were fixed in the pre-registration.

## 3. H3 is falsified, and I am not moving the threshold

H3 pre-declared "plausible" as **CV ≤ 50%**. At CV = 0.50 control CRR is 16.8%
[14.7, 18.9] — the CI excludes the 20% band edge. **H1 is reachable only at CV ≈ 0.70.**

One could argue 70% between-patient CV on immune rate constants is defensible; inter-patient
variability in immunology is genuinely large. **That argument is not available to me here.**
I declared 0.50 before seeing the table, and `lessons.md` §2 is explicit that choosing a
policy after seeing the comparison is how this project has gone wrong before. So:

> **The model survives the placebo falsifier, conditional on accepting ~70% between-patient
> variability — more than I pre-declared plausible.** Anyone rejecting that spread should
> read this round as the falsifier still standing.

## 4. The unpredicted result: a variability–power squeeze

Not predicted, not asked for, and more useful than the falsifier:

**As variability rises, the control response rate rises toward realism, but the treatment
effect shrinks** — +14.4 pp at CV 0.40, +13.3 at 0.50, +11.6 at 0.70, **+8.9 at 1.00.**

The two requirements pull against each other. The window where a realistic placebo rate and
a ≥10 pp separation coexist is narrow — roughly **CV 0.6–0.8**, and it closes entirely by
CV 1.0, where H2 fails.

**Why it matters.** The between-patient heterogeneity needed to reproduce a realistic lupus
placebo arm is itself large enough to erode the detectable treatment effect. This is a
statement about **trial power**, derived from a mechanism rather than assumed in a power
calculation: a population heterogeneous enough to behave like a real LN cohort is one in
which a true drug effect looks smaller. It also offers a mechanistic account of why lupus
trials fail so often, without invoking a bad drug.

**[HYPOTHESIS]** — one model structure, one operating point, one seed per cell. Worth a
dedicated round, and it feeds contribution A directly.

## 5. The trivial comparator — and for once the mechanism wins

Prereg step 3 required it. A 2-parameter exponential fitted to the control-arm mean UPCR
(r = 0.00462/day, c = 0.554, RMSE 0.107), with log-normal variability on **both** its own
parameters:

| CV | 0.10 | 0.30 | 0.50 | 0.70 | 1.00 |
|---|---|---|---|---|---|
| trivial-model CRR | 0.0% | 0.2% | 2.7% | 5.9% | **9.1%** |

**It never reaches the band — 9.1% at CV = 1.00, against the mechanistic model's 33.9%.**

This project's trivial rules usually tie or win (R15, R22, R24). Here they do not, and the
reason is structural: `UPCR = p0·A·N` is multiplicative in several jittered parameters, so
log-normal spread produces a far heavier left tail than a 2-parameter curve can. **This is
the first result in the project where the mechanistic structure earns its keep against the
trivial comparator on a population endpoint.** Stated narrowly: on *this* endpoint, not in
general.

## 6. Bugs and faults found

- **v1 killed at its own 300 s timeout (exit 143).** Two causes, both named in
  `lessons.md`: it inherited `rtol = 1e-11` from the finite-difference code (R36 lost a job
  to exactly this), and it printed only at the end, so the kill lost everything.
- **Fixed and verified, not assumed:** at `rtol = 1e-7` UPCR₅₂ differs from the tight solve
  by **2e-9** relative. Checkpointing now writes per CV, which is the only reason CV 0.1–0.4
  survived the kill.
- Importing `t1_feedback` re-executes its module-level sweep (~6 s per invocation). Noted,
  **not** silently patched, per the surgical-changes rule.

## 7. What this does not establish

- **A matching rate is not evidence the mechanism is right.** R24 showed near-identical
  trajectories from different hypotheses. This tests *reachability of a rate*, nothing more.
- Synthetic throughout. Illustrative parameters, one operating point, one model structure.
- **One seed per cell.** Binomial CIs are reported, but seed-to-seed variation in the jitter
  draw is unmeasured. Multi-seed repetition is outstanding.
- The treated arm is a stylised voclosporin-like effect (45% extra suppression + 25%
  haemodynamic UPCR drop). It was not fitted to AURORA, so the closeness in §2 is
  encouraging rather than evidential.
- CRR here uses two of AURORA-1's criteria; rescue medication and steroid-dose criteria are
  not modelled.

## 8. Next

1. **Multi-seed the CV 0.5–0.8 window.** Cheap, and it is the only thing standing between
   §4's squeeze and a real result.
2. **Follow the anomaly: the variability–power squeeze** (§4). Characterise where the
   window closes. It feeds contribution A and is the most interesting thing here.
3. Re-examine whether 70% CV is defensible **against external data**, not against my own
   intuition — virtual-population literature reports inter-patient spreads that could
   settle it.
4. Unchanged from R42: clear gate A (PMID 42006208) and G0.

*Research and simulation only. Not clinical or drug-development guidance.*
