"""Audit r125: regenerate round-121 exploratory bootstrap CI for c (no committed script produced q25_exploratory_ci.txt).
Usage: python3 -I r121_ci_recheck.py <dir with q25_MMF_patients.csv, q25_AZA_patients.csv>. Percentile bootstrap over patients."""
import sys, numpy as np, pandas as pd
D = sys.argv[1]
for drug in ['MMF', 'AZA']:
    P = pd.read_csv(f'{D}/q25_{drug}_patients.csv', index_col=0)
    for seed in [0, 1, 2]:
        rng = np.random.default_rng(seed); cs = []
        for _ in range(5000):
            d = P.iloc[rng.integers(0, len(P), len(P))]
            if d.R.nunique() < 2: continue
            A = np.column_stack([np.ones(len(d)), d['first'], d.R])
            cs.append(np.linalg.lstsq(A, d['last'].values, rcond=None)[0][2])
        lo, hi = np.percentile(cs, [2.5, 97.5])
        print(drug, 'seed', seed, 'c 95% CI', round(lo, 3), round(hi, 3), '| SD first', round(P['first'].std(), 3),
              '| corr first-last', round(P[['first', 'last']].corr().iloc[0, 1], 3), '| lo/SD', round(lo / P['first'].std(), 3))
