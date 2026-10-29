# Literature triage

The recurring screening question behind gates A, B and G0: **does this paper's model
carry compartment X as a state variable?** Not "does the paper mention X."

`keyword_baseline.py` is the trivial comparator, run before any decision model is
considered (round 41 barred Geneformer for losing to PCA; same standard here).

## Result

Precision **0.33** on 9 hand-labelled cases from rounds 37/40/42 — accuracy 77.8%,
recall 1.00, two false positives. Both false positives really happened to us:

- **PMC5410990** discusses long-lived plasma cells at length and is a review of *mouse*
  models, carrying no mathematical model at all. Round 42 had to read it to find out.
- A 2025 *Nat Commun* paper measures autoreactive CD19⁻ plasma cells experimentally and
  models nothing.

**The failure is structural, not a tuning problem.** The discriminating information is the
grammatical role of the term, not its presence, so no keyword list separates a review that
discusses a cell type from a model that carries it. Both false positives arise from the
term being present *and correct*.

## Status

`n = 9 is not a benchmark.` This establishes that the task is non-lexical. It establishes
no accuracy figure for any model, and the labels are not independent — they are ours.

The open item is a labelled set of ~50–100 abstracts, seeded from the paper
determinations already in rounds 37/40/41/42. That set is worth building whether or not a
decision model is ever used: it turns contribution B's screening from unmeasured judgement
into a measured task.

See `docs/rounds/round-45-jev.md` for why a System One decision model cannot enter the
model architecture, and what would have to be true before one is admitted here.

*Research and simulation only. Not clinical or drug-development guidance.*
