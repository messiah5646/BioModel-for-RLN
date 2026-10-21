# Research state

> Living document. Last updated: 2026-10-22.

This file tracks where the BioModel-for-RLN thinking actually stands — distinct from the aspirational README. It is the in-repo mirror of the Notion project hub.

## Current stage

- Repo scaffolding only. No model or data code committed yet.
- Core premise fixed: a **JEPA-style latent predictive world model** of biological state, with a reasoning layer for proposing/ranking interventions.
- RLN is the chosen biological anchor domain.

## Open questions (to be filled in)

- State representation: what level (molecular / pathway / cell / tissue) and what modality the encoder consumes.
- Perturbation encoding: how interventions are represented in the model's action/condition space.
- Data: which real datasets ground the model, and the access plan for each.
- Evaluation: what counts as a successful *in silico* prediction, and the benchmark it is measured against.
- Novelty: what is genuinely new vs. existing virtual-cell / biology-foundation-model work.

## Notes

Detailed findings, prior-art checks, and decision history are maintained in the Notion hub and will be ported into `docs/` as they stabilise.
