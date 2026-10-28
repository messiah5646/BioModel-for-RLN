# Pre-registration — Round 44: is the variability–power squeeze real?
Written 2026-10-28 BEFORE any test number. Seeds 0-7 (declared now).

## The claim under attack
R43 §4 (unpredicted): as between-patient CV rises, control CRR rises toward realism while
the treated−control gap shrinks (+14.4 → +13.3 → +11.6 → +8.9 pp over CV 0.4→1.0).
R43 labelled it [HYPOTHESIS] on one seed per cell. This round tries to kill it.

## The alternative explanation I must rule out (the obvious one)
**Binary-threshold artifact.** CRR is a threshold on a continuous quantity. As CV → ∞,
the fraction of parameter space satisfying UPCR ≤ 0.5 is increasingly determined by the
spread rather than by treatment, so ANY two arms must converge. If so, the "squeeze" is
arithmetic, not biology, and carries no trial-design meaning.

**Discriminating test:** compute a CONTINUOUS effect size (standardised mean difference in
log UPCR₅₂) on the same populations.
- If the squeeze appears in the binary endpoint ONLY → threshold artifact. **Retract R43 §4.**
- If it appears in the continuous effect size TOO → the erosion is in the underlying
  distributions, not the threshold. R43 §4 survives.

## Hypotheses
- **H1.** Multi-seed (8 seeds, n=600/arm) confirms the binary squeeze: diff_pp at CV 0.8
  is below diff_pp at CV 0.4, with non-overlapping seed-SD intervals.
- **H2 (the real test).** Continuous standardised effect size (Cohen's d on log UPCR₅₂)
  ALSO declines over the same CV range. *Falsified if d is flat or rising* → §4 retracted
  as a threshold artifact.
- **H3.** Required sample size per arm (80% power, α=0.05, two-proportion) rises
  monotonically with CV across 0.4→1.0.

## Method, declared in advance
- CVs 0.40, 0.50, 0.60, 0.70, 0.80, 1.00; seeds 0–7; n=600/arm/cell.
- **Paired within seed**: control and treated arms drawn from the SAME seed, so the
  comparison is paired (lessons §2).
- Report **seed SD separately** from binomial CI.
- Checkpoint to disk per cell (R43's kill lost an unsaved run).
- rtol 1e-7, already verified against the tight solve at 2e-9 relative.

## Inclusion rule, declared before results
Exclude a simulated patient only on integrator failure or eGFR clamp. No other exclusion.

## What it cannot establish
Still synthetic, one model structure, one operating point. A surviving squeeze is a
property of THIS model, not a measured fact about lupus trials.
