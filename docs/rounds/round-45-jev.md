# Round 45 — Jev and the System One decision models: useful, but not where it looks (2026-10-29)

Mixed round: §4 is computed, everything else is retrieval and therefore hypothesis
(rule 1). Code: `src/t8_triage.py`. Raw: `out/t8_triage.json`.

**Verdict: a Jev-class model cannot enter the architecture at all — that is a category
error — but it addresses a real and currently blocking bottleneck in the research loop.
Jev itself is not admissible on present evidence, for a reason that is methodological
rather than about quality. The open-weights alternatives are the admissible route, and
neither can be tested from this container.**

---

## 1. What Jev is

TypeSafe AI's first "System One" model, released 15 Sept 2026; Diogo Almeida (co-inventor
of RLHF/InstructGPT) as CEO. It is **not an LLM**. You send unstructured state plus typed
questions and get typed decisions back:

- **Choice** — one of up to 255 options
- **Score** — an ordered rubric, 2–10 levels
- **Bool/Noul** — a yes/no with a probability in [0,1]

Each answer carries a confidence derived from the shape of the probability distribution.
Questions are evaluated in parallel and in isolation, non-autoregressively. Text-only,
64k context. Vendor list price $0.042/M input tokens with output unmetered, 70–500 ms
end-to-end, claimed 193.6× faster and 444.6× cheaper than frontier LLMs.

Training method is called RLCD (Reinforcement Learning for Calibrated Decisions).

## 2. It cannot enter the architecture — and this is the main point

The instinct is to put it at **L5 (experiment selection)** or in the **L4 provenance
gate**. Both are wrong, and not marginally:

- **L5 ranks experiments by expected information per unit cost, computed from the Fisher
  information matrix.** That is an exact calculation, and `argmax` over it is optimal
  given the model. Substituting a learned decision would replace an exact computation
  with an unverified guess. Strictly worse.
- **The L4 provenance gate is `if marginal_sigma < threshold: release`.** A comparison
  operator. Its entire value is that releases are deterministic and auditable — putting a
  probabilistic model there would destroy the property the gate exists to provide.
- **It cannot satisfy the §5.1 module contract.** A module needs `drift()`, internal
  state, and kinetics. Jev is text-only and has no dynamics. §5.2's NC3 law means a
  module without internal state is dose-equivalent anyway.

**Nothing in L1–L5 is bottlenecked on judgment.** Those layers are bottlenecked on
computation and on data access. A decision model helps with neither.

## 3. Where the real slot is

`lessons.md` §6: **"Audit is the bottleneck, not generation."** That is where this fits.

Right now **three gates are open and all three are literature-triage work**:
- **Gate A** — the *Kidney Int Rep* trial-design review, blocked (reCAPTCHA + publisher
  errors), gating the lead contribution.
- **Gate B** — the identifiability atlas needs a structured screen over many published
  lupus/LN models.
- **G0** — uncleared since round 40, gating every novelty claim about LLPC omission.

The recurring question in all of them is a **Noul**: *does this paper's model carry
compartment X as a state variable?* That is precisely the typed-decision shape.

## 4. Is the task hard enough to need a model? COMPUTED: yes, and the baseline failure is structural

Round 41 barred Geneformer for losing to PCA. Same standard applies, so the trivial
comparator comes first: a keyword match for `long-lived plasma cell | CD19-negative | LLPC`.

Hand-labelled set of 9 cases taken from our own round-37/40/42 determinations:

| | |
|---|---|
| accuracy | 77.8% (1 TP, 6 TN, **2 FP**, 0 FN) |
| **precision** | **0.33** |
| recall | 1.00 |

**Both false positives are cases that actually happened to us.**

- **PMC5410990.** Round 42 retrieved it, it discusses long-lived plasma cells at length
  ("cyclophosphamide does not delete long-lived plasma cells efficiently... persist in
  survival niches"), and it is **not a mathematical model at all** — it is a review of
  *mouse* models. I had to read it to find that out, and it could not clear G0.
- **The 2025 *Nat Commun* paper** on autoreactive CD19⁻ plasma cell subsets, flagged in
  round 42 as experimental rather than modelling.

**Why this is not a tuning problem.** The question is not "does the text mention X" but
"does the *model* carry X as a state variable." Those come apart in both directions: a
review can discuss the biology thoroughly while carrying no model, and an experimental
paper can measure the exact cell type while modelling nothing. The discriminating
information is the **grammatical role** of the term, not its presence. No keyword list
separates them, because the false positives arise from the term being present *and
correct*. Adding keywords cannot fix that.

So: at precision 0.33 we read ~3 papers per real hit. On a task where audit is the
bottleneck, a calibrated screen at even ~0.8 precision is a 2–3× reduction in reading
load. **That is the actual, quantified value proposition — and it is in the research
tooling, not in the model.**

## 5. Why Jev itself is not admissible on present evidence

Not a quality judgement. The disqualifier is the evaluation method:

> The vendor's accuracy evidence comes from internal workflow evals that **use the average
> of GPT-6 Astra and Fable 5.1 as reference probabilities, not ground-truth labels.**

That measures **agreement with other LLMs, not correctness.** It is `lessons.md` §5 in a
new domain — *"a model evaluated against the simulator that generated its training data
has measured nothing."* Our triage question has a real ground truth: whether a paper's
model carries a state variable is a fact, settled by reading the paper. LLM-agreement is
not evidence about it, and the vendor notes its own reference answers favour OpenAI and
Anthropic models.

Compounding:
- **No public benchmarks, deliberately.** The vendor FAQ asks the question and does not
  answer it. Round 41's complaint about missing simple baselines, again.
- **Calibration is claimed with no metrics and no data.** Calibration is the *one*
  property that would make a triage screen usable — you must know when to distrust it and
  escalate to a human read. Unverified.
- **"0% hallucination" is not empirical** — it is schema matching guaranteed by design.
  Type-safety is not correctness; a well-typed answer can be wrong.
- **Closed waitlist early access**, with a per-version "jaggedness" page. A closed,
  versioned API that can change under us conflicts with rule 1: a number must stay
  re-derivable from code in the repo.

## 6. The open-weights alternatives — the admissible route, unverified

Several exist, and they fix the reproducibility objection in principle:
`ikermoel/open-alternative-jev` (typed calibrated decisions from any open-weights LLM in
one forward pass, HF + vLLM, advertising "honest benchmarks"), `wfzyx/von` (sub-15 ms,
non-autoregressive, local), `allebee/jevk5` (Apache-2.0 weights and code). One aggregator
claims 54 open alternatives.

**Two honest caveats.** First, every source I found for these is SEO-grade content
(`madewithjev.com`, `respan.ai`, `flowtivity.ai`, `metatext.io`, DataCamp) and **I verified
none at source.** Round 37 recorded this exact failure: leads from a third-party GitHub
list went unverified, and the one that was checked "turned out weaker than its label."

Second, **neither can be tested from this container.** Tested per host, per
`lessons.md` §9:

| host | result |
|---|---|
| `pypi.org` | 200 |
| `raw.githubusercontent.com` | 301 |
| `github.com` | 400 |
| **`huggingface.co`** | **unreachable** |
| **`api.typesafe.ai`** | **unreachable** |

No weights, no API. Plus two cores and no GPU. So an empirical verdict is **blocked**,
not pending.

## 7. The gate before anything is admitted

**Build our own labelled triage set, then benchmark against the keyword baseline.** We
already have a seed: every paper determination in rounds 37/40/41/42 is a hand label, and
§4's trap cases are the discriminating examples. Target ~50–100 labelled abstracts.

That set is worth building **regardless of whether any decision model is ever used**: it
converts contribution B's screening from unmeasured judgement into a measured task, and
it is the only thing that could justify admitting a model under this project's own rules.

**n = 9 is not a benchmark.** §4 establishes that the task is structurally non-lexical. It
does **not** establish any accuracy figure for anything.

## 8. What this does not establish

- §1, §5, §6 are retrieval. No vendor claim was reproduced; most sources are secondary and
  several are content farms.
- §4's snippets are **paraphrases recorded in our own round docs, not verbatim abstracts
  fetched now.** It measures task structure, not real-world accuracy. The labels are my
  own, so the set is not independent.
- No decision model of any kind was run. There is no comparison between Jev, an open
  alternative, and an LLM — only between the keyword baseline and the ground truth.
- The claim that L5/L4 cannot benefit is an argument from their structure, not an
  experiment. It would be falsified by a decision problem in the loop that is genuinely
  judgement-bound rather than computation-bound; I did not find one, which is weak evidence.

## 9. Next

1. **Build the labelled triage set** (~50–100 abstracts from rounds 37/40/41/42 plus new
   retrieval). Decision-useful on its own, and it is the admission gate.
2. Verify two or three open-alternative repos **at source** before treating the ecosystem
   as real.
3. If a model is ever admitted, it screens and **ranks for human reading**; it never
   clears a gate by itself. G0 and gate B are novelty claims, and round 34's rule gives
   those one verdict after primary sources are read **by us**.
4. Unchanged and still ahead of this: additive-vs-multiplicative jitter, whether the
   R43/R44 CV window is empty, gate A.

*Research and simulation only. Not clinical or drug-development guidance.*
