# Workarounds for every standing blocker

Live document. Each row carries a **tested** status, not a plan. Routes marked FAILS were
attempted and failed; routes marked WORKS were attempted and succeeded this session.

---

## W1 — Holdout erosion (self-inflicted, now the most urgent)

**Blocker.** Rounds 48 and 49 have each scored against the same two held-out numbers.
A third test designed to fix round 49's miss would be fitting, not prediction.

**Workaround.** Re-harvest the **per-visit** series rather than summary endpoints. Round 35
established the EU Clinical Trials Register results sections carry full tabulated per-arm
per-visit endpoint tables with n and SD — 13 visits in year 1 per arm. Then fit on the
early visits and hold out the later ones, turning n = 2 into roughly 30 held-out points
from one source.

**Status.** Route **confirmed reachable** this session: the register search returned 64
lupus nephritis trials, **16 with results sections**. Per-visit table extraction not yet
re-run. This is the single highest-value unblocked action.

**Fallback.** Whole-trial holdout — fit one trial, predict another — using the other 15
trials with results.

## W2 — Gate A: the trial-design review (gates the lead contribution)

**Blocker.** PMC and PubMed serve reCAPTCHA; publisher fulltext returns client errors.

**Routes attempted and FAILED:** PMC fulltext (reCAPTCHA), PubMed abstract (reCAPTCHA),
publisher fulltext (client_error), trade summary (client_error), DOI resolver
(client_error), **publisher PDF (robots_disallowed)**.

**Routes still untried:** Europe PMC REST, Semantic Scholar API, CORE, scholar.archive.org,
DOAJ (the journal is open access, so a compliant full-text route should exist), the
journal's own issue-level table of contents, or a browser session with the desktop app.

**Status.** **Six routes failed.** Reclassified from "blocked pending retry" to **needs a
different class of access** — a browser session, or an institutional route.

**Reframe that removes the blocker's urgency:** the review answers *which endpoints trials
use and how much they vary*. Contribution A answers *which parameters are estimable from a
given schedule*. Those are different questions, so the review is a **prior-art check, not
a dependency**. Contribution A can be drafted with the gate explicitly open.

## W3 — G0: do published lupus models omit long-lived plasma cells?

**Blocker.** Uncleared since round 40. Round 42's only locatable review turned out to be
about mouse models and listed no state variables.

**Workaround.** Stop trying to clear it from reviews. **Read the primary models'
equations.**

**Status — WORKS.** One of four targets cleared at source this session: the Bull Math Biol
NETosis model's full four-equation system was retrieved from the publisher. Its states are
apoptotic material, autoantigen, activated macrophages, neutrophils — **no long-lived
plasma cells, no CD19-negative plasma cells, no regulatory T cells**, and the authors
explicitly state they omit B-cell antigen presentation and autoantibody production.

**Remaining three targets**, same route: the Biomath chronic-phase SLE model, the
integrative type-1-IFN SLE model, and the ODE-informed ML stratification paper. Two are
open access.

## W4 — Patient-level longitudinal data

**Blocker.** Vivli, NIDDK-CR, ImmPort controlled tier and UK Biobank all key on an
institutional signing official. Not technical.

**Workarounds, in cost order.** ImmPort registration is free with no stated affiliation
requirement and same-day (round 37). The population-scale LN atlas at Synapse
**syn64064827** (156 LN subjects, 538k kidney cells) is days, not months. UK Biobank is
closed to new applications until late 2026.

**Status.** Untested this session. **But the identifiability programme needs none of it** —
that is the whole point of `positioning.md`. Treat this as an upgrade, not a dependency.

## W5 — The BLISS-LN clinical study report

**Blocker.** Genuinely public at GSK's register, but the file URL exceeds WebFetch's
length limit (round 35).

**Workaround.** Fetch the register's trial page and follow whatever shorter link it
exposes, rather than constructing the file URL. Failing that, the
clinicalstudydatarequest route.

**Status.** Untested. Valuable because this agent has no haemodynamic proteinuria effect,
which removes the route-3 confound that §3 of round 49 is now stuck on.

## W6 — Decision-model weights (round 45)

**Blocker.** `huggingface.co` and `api.typesafe.ai` both unreachable; two cores, no GPU.

**Workaround.** None worth taking. `pypi.org` is reachable, but round 45's verdict already
turns on a labelled evaluation set we do not have, and building that set is the real gate.
**Deprioritised, not blocked.**

## W7 — Shell egress to data archives

**Blocker.** All archives refused at the gateway; only pypi and raw.githubusercontent.

**Standing workaround.** Test **per host and per tool** — WebFetch and the shell take
different routes, and this session confirmed it again: the EU register is unreachable by
curl and reachable by WebFetch. Never conclude a host is blocked from one tool's failure.

---

## The pattern worth keeping

Four of these seven dissolve rather than resolve: W2 and W4 stop being blockers once the
work is scoped so they are not dependencies, W6 is deprioritised, and W7 is a testing
discipline rather than an obstacle. **The two that genuinely bind are W1 and W5** — both
are data-harvesting jobs on a source already confirmed reachable, and W1 is the one that
decides whether round 49's next test counts as science.

*Research and simulation only. Not clinical or drug-development guidance.*
