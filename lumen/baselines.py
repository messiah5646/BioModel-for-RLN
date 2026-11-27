"""L2 baselines and scores. Trivial baselines are a first-class component (constraint C2), not a benchmark afterthought.

Scores take arrays of held-out observations (obs log OR, its sampling SE) and a model's predictive mean and SD.
  iv_mse    : the project's primary metric since round 112 (weights 1/se_obs^2) — same as run_rt1.ivmse
  log_score : mean negative log predictive density of obs under N(mean, sd^2 + se_obs^2); proper, rewards honest width
  coverage  : share of obs inside mean +/- 1.96 sqrt(sd^2 + se_obs^2) — same as run_rt1.coverage
"""
from __future__ import annotations
import math
import numpy as np


def null(n: int):
    """M0: no effect, log OR 0. Predictive SD is not defined for a point null; give 0 (scores use se_obs only)."""
    return np.zeros(n), np.zeros(n)


def carry_forward(early_lor, early_se):
    """C0: the trial's own early read-out carried to the late visit, sampling SE only (no drift added)."""
    return np.asarray(early_lor, float), np.asarray(early_se, float)


def iv_mse(obs, se_obs, mean) -> float:
    obs, se_obs, mean = (np.asarray(x, float) for x in (obs, se_obs, mean))
    w = 1 / se_obs ** 2
    return float((w * (obs - mean) ** 2).sum() / w.sum())


def log_score(obs, se_obs, mean, sd) -> float:
    obs, se_obs, mean, sd = (np.asarray(x, float) for x in (obs, se_obs, mean, sd))
    v = sd ** 2 + se_obs ** 2
    return float(np.mean(0.5 * np.log(2 * math.pi * v) + (obs - mean) ** 2 / (2 * v)))


def coverage(obs, se_obs, mean, sd, z: float = 1.96) -> float:
    obs, se_obs, mean, sd = (np.asarray(x, float) for x in (obs, se_obs, mean, sd))
    return float((np.abs(obs - mean) <= z * np.sqrt(sd ** 2 + se_obs ** 2)).mean())
