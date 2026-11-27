"""Candidates are checked on SYNTHETIC data only. Their performance on real trials is PREREG-O1's job."""
import numpy as np
from lumen import candidates as c
from lumen.layers.drug_effect import core


def _sim(k=400, tau=0.3, seed=1):
    rng = np.random.default_rng(seed)
    se = rng.uniform(0.15, 0.6, k)
    return 0.3 + rng.normal(0, tau, k) + rng.normal(0, se), se


def test_dl_identical_to_shipped():
    y, se = _sim(25)
    assert abs(c.tau2_dl(y, se) - core().dl(y, se)) < 1e-12


def test_estimators_recover_tau():
    y, se = _sim()
    for f in (c.tau2_dl, c.tau2_pm, c.tau2_reml):
        assert abs(np.sqrt(f(y, se)) - 0.3) < 0.06, f.__name__


def test_bayes_predictive_close_to_truth_with_many_trials():
    y, se = _sim()
    m, sd = c.predictive_bayes(y, se)
    assert abs(m - 0.3) < 0.05 and abs(sd - 0.3) < 0.07


def test_bayes_not_degenerate_when_dl_gives_zero():
    y, se = np.array([0.30, 0.31, 0.29]), np.array([0.3, 0.3, 0.3])   # DL tau^2 = 0 exactly
    assert c.tau2_dl(y, se) == 0
    _, sd_b = c.predictive_bayes(y, se)
    assert sd_b > np.sqrt(1 / (1 / se ** 2).sum())   # tau prior keeps between-trial spread in the predictive


def test_robust_mixture_widens():
    r = c.robust_mixture(0.3, 0.2)
    assert r["sd"] > 0.2 and abs(r["mean"] - 0.24) < 1e-12
