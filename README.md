# BioModel-for-RLN

A research codebase for a **biology world model** — a JEPA-style predictive model of biological state — used to simulate human and natural biology, fast-forward candidate experiments, and pre-screen trials *in silico* before wet-lab testing. **RLN** is the current biological anchor domain.

> Status: **early scaffolding.** Structure and docs are being laid down; model and data code to follow.

## Idea

Learn a latent world model of biological state (cell / pathway / system level) that predicts how that state evolves under perturbations, in representation space rather than by reconstructing raw measurements (the JEPA premise). A reasoning layer on top proposes and ranks interventions, so that expensive trials can be prioritised — and many dead ends discarded — computationally first.

## Repository layout

```
src/biomodel/     # core package: model, encoders, predictors, objectives
data/             # datasets and loaders (large/raw data gitignored)
configs/          # experiment + model configs
experiments/      # experiment entrypoints and run scripts
notebooks/        # exploratory analysis
scripts/          # utility + data-prep scripts
docs/             # design notes, architecture, research-state write-ups
tests/            # unit / integration tests
```

## Documentation & logs

Design docs, architecture notes, and a dated build log are mirrored in Notion (project hub: *BioModel-for-RLN*). The `docs/` folder here holds the versioned, in-repo copies.

## Status

See [`docs/research-state.md`](docs/research-state.md) for where the thinking currently stands.
