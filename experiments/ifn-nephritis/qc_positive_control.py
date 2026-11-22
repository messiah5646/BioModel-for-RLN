# QC only (not a pre-registered test): IFN6 must separate SLE from healthy if the scoring is right.
import sys, runpy, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
D = sys.argv[1]; sys.argv = ['x', D, '/tmp']
g = {}
src = open('run.py').read().split('num = lambda')[0]; exec(src, g)
for gse, col, sle, hc in [('GSE65391', 'disease state', 'SLE', 'Healthy'), ('GSE72326', 'group', 'SLE', 'Healthy control of SLE')]:
    X, M = g['read_matrix'](gse); E = g['ifn_illumina'](X); M = M[M[col].isin([sle, hc])]
    s = g['zscore_ifn'](E, M.index); y = (M[col] == sle).astype(int)
    print(gse, 'SLE vs healthy IFN6 AUC', round(roc_auc_score(y, s), 3), 'n', len(M))
X, M = g['read_matrix']('GSE99967'); IFN = g['IFN']; ENT = g['ENTREZ']
E = X.loc[[ENT[k] + '_at' for k in IFN]].astype(float).T; E.columns = IFN
s = ((E - E.mean()) / E.std()).mean(axis=1).loc[M.index]; y = (M['disease state'] == 'SLE').astype(int)
print('GSE99967 SLE vs control IFN6 AUC', round(roc_auc_score(y, s), 3), 'n', len(M))
