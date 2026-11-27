# Round 149: phase has never been tested on the relative effect, so O3 is pre-registered

*2026-11-27 · scheduled run · queue item O3 · branch rounds-70-72 · PREREG 05328df (pushed 04:26Z)*

**Decision this can change:** whether LUMEN's drug-effect layer takes trial phase (2 vs 3) as a covariate on the log
odds ratio. If the phase model beats the shipped pooled mean by >= 10% time-forward -> phase enters drug_effect as a
v3.1 candidate (CLI `--phase`); if it fails or is inconclusive -> no phase term, and the question joins rule 8's
settled list.

## What this run did
1. **Took O3.** No PRIORITY item was todo; O1 is waiting for its 1-hour clock (PREREG 4e8975a, 04:15Z); O2 is gated
   behind O1; W1 is due 2026-11-28 (7 days after its last done note).
2. **Rule-8 check (is this already settled?).** Searched the repo and the project. RT-2 (round 129) tested posting
   year and definition stringency on control-arm rates only. Rounds 86, 116 and RT-1 tested molecule/target anchors
   and mechanism class on the effect. **Phase on the relative effect is untested**, so the item goes ahead.
3. **Pre-registered the test** (`experiments/lumen-v3/PREREG-O3.md`) and stopped, as rule 12b requires. Analysis runs
   in a later run, at least 1 hour after the push.

## The test, in plain words
Phase-2 trials often look better than the phase-3 trials that follow them (the classic case is anifrolumab). If that
holds in LUMEN's tables, a forecast for a new phase-2 trial should be pulled toward what phase 2 typically shows,
and a phase-3 forecast toward phase 3. P1 fits a straight line, log OR = a + b × (phase 3), with the same weights the
shipped pooled mean uses, so that without the phase term it is exactly the shipped model. Each trial is forecast
from trials posted earlier only (as in RT-1). Primary score: inverse-variance-weighted squared error, LN + SLE
combined. PASS needs 10% better than the pooled mean; losing to null or pooled by >= 10% is a FAIL (rule 2b).
At least 20 scored trials, with >= 6 of each phase, are needed or the result is UNDERPOWERED.

## What was seen before the PREREG
Only the phase counts (no outcomes): LN 20 trials (12 phase 2, 7 phase 3, 1 phase 2/3); SLE 40 trials (26 phase 2,
13 phase 3, 1 phase 1/2). The outcome tables are committed and were analysed in earlier rounds in other ways, so the
test is labelled **12b-timed (committed data; earlier looks not excludable)**. The analyst knows some public
results that hint at phase-2 optimism, which is why the question was asked.

## Labels
No new empirical claims this round. Phase-2 optimism in LUMEN's tables: **Hypothetical** until O3 runs.

## Next
O1 analysis (due >= 05:15Z), then O3 analysis (due >= 05:26Z), W1 on 2026-11-28.

Sources: PREREG-O3 (commit 05328df); RT-2 round 129; RT-1 round 127; claude/RESEARCH-optimising-LUMEN-2026-11-27.md §5.

*Research and simulation only. Not clinical or drug-development guidance.*
