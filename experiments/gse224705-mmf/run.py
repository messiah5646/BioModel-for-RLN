"""Round 97 (Q11). Usage: python3 -I run.py <data_dir> <out_dir>. Follows PREREG.md exactly."""
import sys, gzip, json, numpy as np, pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score
D, O = sys.argv[1], sys.argv[2]
IFN = ['IFI27','IFI44L','IFIT1','ISG15','RSAD2','SIGLEC1']; PB = ['JCHAIN','IGJ','MZB1','TNFRSF17','TXNDC5']
ann = pd.read_csv(f'{D}/gpl.txt.gz', sep='\t', skiprows=27, usecols=['ID','Gene symbol'], dtype=str).dropna()
want = {p: s.split('///')[0] for p, s in zip(ann.ID, ann['Gene symbol']) if s.split('///')[0] in IFN + PB}  # TXNDC5 probe is TXNDC5///BLOC1S5
rows = {}
with gzip.open(f'{D}/sm.txt.gz', 'rt') as f:
    for l in f:
        if l.startswith('"ID_REF"'): cols = [c.strip('"') for c in l.rstrip('\n').split('\t')[1:]]
        elif l[0] == '"':
            p, *v = l.rstrip('\n').split('\t'); p = p.strip('"')
            if p in want: rows[p] = np.array(v, float)
X = pd.DataFrame(rows, index=cols).T; X['gene'] = X.index.map(want)
X = X.assign(m=X.drop(columns='gene').mean(axis=1)).sort_values('m').groupby('gene').tail(1).set_index('gene').drop(columns='m').T
if 'IGJ' in X and 'JCHAIN' not in X: X = X.rename(columns={'IGJ':'JCHAIN'})
pbg = [g for g in ['JCHAIN','MZB1','TNFRSF17','TXNDC5'] if g in X]
print('genes found', list(X.columns), 'PB used', pbg)
c = pd.read_csv(f'{D}/chars.csv'); c['days'] = pd.to_numeric(c.days_from_first_visit)
c['upcr'] = pd.to_numeric(c['protein_creatinine_ratio_urine_(mg/g)'], errors='coerce')
c['sledai'] = pd.to_numeric(c.sledai, errors='coerce')
out, pts = {}, []
for drug in ['MMF','AZA']:
    s = c[c.drug_group == drug].sort_values('days')
    const = s.groupby('patient_id').sri4_response.nunique().le(1).mean()
    lab = s.groupby('patient_id').sri4_response.last(); base = s.groupby('patient_id').head(1).set_index('patient_id')
    b = base.assign(y=lab); pts.append(b.assign(drug=drug))
    out[drug] = {'n_pat': len(b), 'frac_constant_label': float(const), 'labels': b.y.value_counts().to_dict()}
B = pd.concat(pts)
B['y'] = (B.y == 'YES').astype(int)
Xb = X.loc[B.gsm.values]; Xb.index = B.index
ref = X.loc[B.gsm.values]  # all SLE baseline samples, both drugs
z = (Xb - ref.mean()) / ref.std()
B['IFN6'] = z[IFN].mean(axis=1).values; B['PB4'] = z[pbg].mean(axis=1).values
B['lupcr'] = np.log10(B.upcr + 1)
B[['drug','y','lupcr','sledai','IFN6','PB4','upcr','days']].to_csv(f'{O}/patients.csv')
def oof(d, f):
    p = np.zeros(len(d)); Xv = d[f].values; y = d.y.values
    for i in range(len(d)):
        m = np.ones(len(d), bool); m[i] = False
        p[i] = make_pipeline(StandardScaler(), LogisticRegression(C=1.0)).fit(Xv[m], y[m]).predict_proba(Xv[i:i+1])[0, 1]
    return p
M = {'SLEDAI': ['sledai'], 'UPCR': ['lupcr'], 'IFN6+PB4': ['IFN6','PB4'], 'IFN6': ['IFN6'],
     'base: UPCR+SLEDAI': ['lupcr','sledai'], 'full: +IFN6+PB4': ['lupcr','sledai','IFN6','PB4']}
rng = np.random.default_rng(97)
def analyse(d, tag):
    d = d.dropna(subset=['lupcr','sledai','IFN6','PB4']); r = {'n': len(d), 'n_resp': int(d.y.sum()), 'n_nonresp': int((1-d.y).sum())}
    if min(r['n_resp'], r['n_nonresp']) < 3: print(tag, 'too few in a class', r); return r
    P = {k: oof(d, f) for k, f in M.items()}; r['auc'] = {k: round(roc_auc_score(d.y, v), 3) for k, v in P.items()}
    r['auc']['chance'] = 0.5
    y = d.y.values; a, b = P['full: +IFN6+PB4'], P['base: UPCR+SLEDAI']; diffs = []
    for _ in range(2000):
        i = rng.integers(0, len(d), len(d))
        if y[i].min() == y[i].max(): continue
        diffs.append(roc_auc_score(y[i], a[i]) - roc_auc_score(y[i], b[i]))
    dlt = r['auc']['full: +IFN6+PB4'] - r['auc']['base: UPCR+SLEDAI']
    r['delta'] = round(dlt, 3); r['delta_ci95'] = [round(float(np.percentile(diffs, q)), 3) for q in (2.5, 97.5)]
    r['pass'] = bool(dlt >= 0.05 and r['delta_ci95'][0] > 0)
    # raw (unfitted) direction checks
    r['raw_auc_neg_sledai'] = round(roc_auc_score(y, d.sledai), 3); r['raw_auc_IFN6'] = round(roc_auc_score(y, d.IFN6), 3)
    r['raw_auc_PB4'] = round(roc_auc_score(y, d.PB4), 3); r['raw_auc_lupcr'] = round(roc_auc_score(y, d.lupcr), 3)
    print(tag, json.dumps(r)); return r
res = {'label_structure': out, 'PB_genes_used': pbg + (['TXNDC5'] if 'TXNDC5' in X else []),
       'H1_MMF': analyse(B[B.drug == 'MMF'], 'MMF'), 'H2_AZA': analyse(B[B.drug == 'AZA'], 'AZA'),
       'sens_MMF_UPCR_ge500': analyse(B[(B.drug == 'MMF') & (B.upcr >= 0.5)], 'MMF UPCR>=0.5 g/g')}
json.dump(res, open(f'{O}/result.json', 'w'), indent=1)

# EXPLORATORY (not pre-registered): LOO-AUC is biased toward 0 for uninformative features; repeat with 5-fold x 50
from sklearn.model_selection import RepeatedStratifiedKFold, cross_val_predict
ex = {}
for drug in ['MMF','AZA']:
    d = B[B.drug == drug].dropna(subset=['lupcr','sledai','IFN6','PB4']); ex[drug] = {}
    for k, f in M.items():
        aucs = []
        for rep in range(50):
            from sklearn.model_selection import StratifiedKFold
            p = cross_val_predict(make_pipeline(StandardScaler(), LogisticRegression(C=1.0)), d[f].values, d.y.values,
                                  cv=StratifiedKFold(5, shuffle=True, random_state=rep), method='predict_proba')[:, 1]
            aucs.append(roc_auc_score(d.y, p))
        ex[drug][k] = [round(float(np.mean(aucs)), 3), round(float(np.percentile(aucs, 5)), 3), round(float(np.percentile(aucs, 95)), 3)]
    print('EXPL 5fold', drug, ex[drug])
res['exploratory_5fold_x50_mean_p5_p95'] = ex
res['deviations'] = ['UPCR field is in g/g scale despite name (median 0.11); sensitivity threshold applied as >=0.5 g/g',
                     'TXNDC5 probe annotated TXNDC5///BLOC1S5, used as TXNDC5']
json.dump(res, open(f'{O}/result.json', 'w'), indent=1)
