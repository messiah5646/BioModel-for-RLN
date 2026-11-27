"""LUMEN v3 — one package for the lupus nephritis trial forecaster and its gated module slots.

Layers (see docs/SYSTEM-DESIGN-v3.md):
  L0  provenance      md5-frozen inputs, run stamps                      lumen.provenance
  L1  measurement     arm counts and observables with scale contracts    lumen.contracts
  L2  baselines+gate  null / pooled / carry-forward, proper scores, 2b   lumen.baselines, lumen.gate
  L3  layers          drug effect (shipped), others gated                lumen.layers
  L4  registry        what ships, at which evidence label and n          lumen.registry
The drug-effect maths is the frozen experiments/lumen-bench/lumen_forecast.py (md5 2c11b9a5), imported,
never copied or edited, so every frozen forecast and audit reproduces.
Research and simulation only. Not clinical or drug-development guidance.
"""
__version__ = "3.0.0"
