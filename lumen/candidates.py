"""Optimisation candidates for the drug-effect layer. NOT SHIPPED. Nothing here is called by lumen.layers.

Each candidate is a drop-in replacement for the pooled prior (M1). It may enter lumen.layers only after a
pre-registered held-out test (PREREG-O1, experiments/lumen-v3/) returns PASS through lumen.gate.

Why these (docs/RESEARCH-optimising-LUMEN-2026-11-27.md):
  * M1 estimates tau^2 by DerSimonian-Laird, which is biased low with few or unbalanced trials, and frequentist
    prediction intervals under-cover in that regime (Partlett & Riley 2017, doi:10.1002/sim.7140;
    Seide, Roever & Friede 2019, doi:10.1186/s12874-018-0618-3).
  * M1's centre is the fixed-effect (1/se^2) mean, so a few large trials dominate; random-effects weights
    1/(se^2 + tau^2) are the standard choice when the target is a NEW trial.
  * Bayesian normal-normal models with an empirical or weakly informative tau prior give better coverage with few
    studies (Turner et al. 2015, doi:10.1002/sim.6381; Roever et al. 2021, doi:10.1002/jrsm.1475;
    Lilienthal et al. 2024, doi:10.1002/jrsm.1685).
  * A robust mixture (Schmidli et al. 2014, doi:10.1111/biom.12242) adds a vague component so a truly new
    mechanism is not over-shrunk to the historical mean.
"""
from __future__ import annotations
import math
import numpy as np

# ---- tau^2 estimators -------------------------------------------------------------------------------------------


def tau2_dl(y, se):
    """DerSimonian-Laird; identical to lumen_forecast.dl (the shipped estimator)."""
    y, v = np.asarray(y, float), np.asarray(se, float) ** 2
    w = 1 / v; m = (w * y).sum() / w.sum(); Q = (w * (y - m) ** 2).sum()
    return max(0.0, (Q - (len(y) - 1)) / (w.sum() - (w ** 2).sum() / w.sum()))


def tau2_pm(y, se, tol=1e-10, it=200):
    """Paule-Mandel: solve sum w_i (y_i - mu_w)^2 = k - 1 with w_i = 1/(v_i + tau^2)."""
    y, v = np.asarray(y, float), np.asarray(se, float) ** 2
    k = len(y); lo, hi = 0.0, max(1.0, 10 * float(np.var(y)))

    def f(t2):
        w = 1 / (v + t2); m = (w * y).sum() / w.sum()
        return (w * (y - m) ** 2).sum() - (k - 1)
    if f(0.0) <= 0:
        return 0.0
    while f(hi) > 0:
        hi *= 2
    for _ in range(it):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
        if hi - lo < tol:
            break
    return (lo + hi) / 2


def tau2_reml(y, se, it=500, tol=1e-10):
    """Restricted maximum likelihood by fixed-point iteration (Viechtbauer 2005 form)."""
    y, v = np.asarray(y, float), np.asarray(se, float) ** 2
    t2 = tau2_dl(y, se)
    for _ in range(it):
        w = 1 / (v + t2); m = (w * y).sum() / w.sum()
        new = max(0.0, ((w ** 2) * ((y - m) ** 2 - v)).sum() / (w ** 2).sum() + 1 / w.sum())
        if abs(new - t2) < tol:
            return new
        t2 = new
    return t2


# ---- predictive distributions for a NEW trial's log OR -----------------------------------------------------------


def predictive_frequentist(y, se, tau2_fn=tau2_reml, re_weights=True, t_quantile=True):
    """Return (mean, sd, df) of the predictive distribution for a new trial.

    sd^2 = tau^2 + Var(mu_hat). With t_quantile the 95% interval uses t_{k-2} (Higgins, Thompson & Spiegelhalter 2009)
    instead of 1.96; callers turn that into an equivalent normal sd via sd * t / 1.96 for scoring.
    """
    y, se = np.asarray(y, float), np.asarray(se, float)
    t2 = tau2_fn(y, se)
    w = 1 / (se ** 2 + t2) if re_weights else 1 / se ** 2
    mu = float((w * y).sum() / w.sum()); var_mu = float(1 / w.sum())
    sd = math.sqrt(t2 + var_mu)
    df = max(1, len(y) - 2)
    return mu, sd, (df if t_quantile else None)


def predictive_bayes(y, se, log_tau2_prior=(-2.56, 1.74), half_normal_tau=None, grid=400):
    """Normal-normal hierarchical model, flat prior on mu, tau prior on a grid; returns predictive (mean, sd).

    log_tau2_prior: log-normal prior on tau^2 (mean, sd of log tau^2). The default (-2.56, 1.74) is the value commonly
    quoted from Turner et al. 2015 for a general setting; CHECK the setting-specific row of their table (pharmacological
    vs placebo/control) before any PREREG fixes it. half_normal_tau: if given, use a half-normal(scale) prior on tau
    instead (Roever et al. 2021 suggest 0.5 for log OR as weakly informative).
    """
    y, v = np.asarray(y, float), np.asarray(se, float) ** 2
    taus = np.linspace(1e-4, 2.5, grid)
    if half_normal_tau is not None:
        lp = -0.5 * (taus / half_normal_tau) ** 2
    else:
        m, s = log_tau2_prior
        lt2 = np.log(taus ** 2)
        lp = -0.5 * ((lt2 - m) / s) ** 2 + np.log(2 / taus)   # density of tau from log-normal tau^2 (Jacobian 2/tau)
    ll = np.empty_like(taus); mus = np.empty_like(taus); vmus = np.empty_like(taus)
    for i, t in enumerate(taus):                               # marginal likelihood with mu integrated out (flat prior)
        w = 1 / (v + t ** 2); mu = (w * y).sum() / w.sum()
        mus[i], vmus[i] = mu, 1 / w.sum()
        ll[i] = 0.5 * np.log(vmus[i]) - 0.5 * np.sum(np.log(v + t ** 2)) - 0.5 * np.sum(w * (y - mu) ** 2)
    lpost = lp + ll; p = np.exp(lpost - lpost.max()); p /= p.sum()
    mean = float((p * mus).sum())
    second = float((p * (vmus + taus ** 2 + mus ** 2)).sum())   # E[Var | tau] + E[mean^2 | tau]
    return mean, math.sqrt(max(second - mean ** 2, 0.0))


def robust_mixture(mean, sd, w_vague=0.2, vague_sd=1.0):
    """Schmidli-style robustification: (1 - w) N(mean, sd^2) + w N(0, vague_sd^2). Returns the mixture's mean and sd,
    plus the components for exact log-score evaluation."""
    m = (1 - w_vague) * mean
    second = (1 - w_vague) * (sd ** 2 + mean ** 2) + w_vague * vague_sd ** 2
    return dict(mean=m, sd=math.sqrt(second - m ** 2), components=[(1 - w_vague, mean, sd), (w_vague, 0.0, vague_sd)])
