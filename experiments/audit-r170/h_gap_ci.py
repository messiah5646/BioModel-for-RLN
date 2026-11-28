"""Audit r170 check: bootstrap CI for the time-forward gap of CAL1 baseline H (r168 labels its over-prediction
Supported without a CI). Same units (trial means), same bootstrap as run_cal1.summarise. Descriptive audit check of a
label; not a re-test (rule 8): no decision rule is applied and nothing printed by LUMEN changes from it."""
import sys, pathlib, numpy as np, pandas as pd
R = pathlib.Path(__file__).resolve().parents[2]
U = pd.read_csv(R / 'experiments/lumen-v3/cal1_units.csv')
T = U.groupby('trial')[['H', 'E1', 'P1']].mean()
rng = np.random.default_rng(20261127); n = len(T)
H, E, P = (T[c].to_numpy() for c in ('H', 'E1', 'P1'))
g = [float((H[i] - E[i]).mean()) for i in (rng.integers(0, n, n) for _ in range(4000))]
print(f"n_trials {n}; H mean {H.mean():.4f}; obs {E.mean():.4f}; gap {H.mean()-E.mean():+.4f}; "
      f"95% CI [{np.percentile(g,2.5):+.4f}, {np.percentile(g,97.5):+.4f}]; P(gap>0.10) {np.mean(np.array(g)>=0.10):.2f}")
