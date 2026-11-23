"""Q22 (round 119) analysis exactly as PREREG-Q22.md. Usage: python3 -I analyse.py DATADIR OUTDIR"""
import sys, glob, os, json, numpy as np, pandas as pd, scipy.io
from scipy.stats import spearmanr
G = ['IFI27', 'IFI44L', 'IFIT1', 'ISG15', 'RSAD2', 'SIGLEC1']
F_LATE = 0.73
src, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True)
rows, cells = [], []
for m in sorted(glob.glob(os.path.join(src, '*_matrix.mtx.gz'))):
    base = m[:-len('_matrix.mtx.gz')]; name = os.path.basename(base).split('_', 1)[1]  # PRE1 / POST1
    tp = 'pre' if name.startswith('PRE') else 'post'; pat = 'P' + name[-1]
    M = scipy.io.mmread(m).tocsc()
    feat = pd.read_csv(base + '_features.tsv.gz', sep='\t', header=None)
    sym = feat[1].astype(str).values
    if M.shape[0] != len(sym): M = M.T.tocsc()
    miss = [g for g in G + ['LYZ', 'CD3E'] if g not in set(sym)]
    if miss: sys.exit(f'STOP: missing genes {miss} in {name}')
    tot = np.asarray(M.sum(0)).ravel(); keep = tot > 500
    M = M[:, keep]; tot = tot[keep]
    def lg(g): return np.log1p(M[np.where(sym == g)[0][0], :].toarray().ravel() / tot * 1e4)
    ngen = np.asarray((M > 0).sum(0)).ravel()
    mt = np.where(pd.Series(sym).str.startswith('MT-').values)[0]
    mito = np.asarray(M[mt, :].sum(0)).ravel() / tot
    lyz = lg('LYZ'); mono = lyz > 2; tcell = (lg('CD3E') > 1) & ~mono
    sc = np.mean([lg(g) for g in G], axis=0)
    rows.append(dict(patient=pat, tp=tp, n_cells=int(keep.sum()), n_mono=int(mono.sum()), n_t=int(tcell.sum()),
                     mono=sc[mono].mean(), t=sc[tcell].mean()))
    cells.append(pd.DataFrame(dict(patient=pat, tp=tp, s=sc[mono], lc=np.log(tot[mono]), lg=np.log(ngen[mono]),
                                   mito=mito[mono])))
S = pd.DataFrame(rows); S.to_csv(os.path.join(out, 'scores.csv'), index=False); print(S.round(4).to_string())
C = pd.concat(cells)
X = np.c_[np.ones(len(C)), C.lc, C.lg, C.mito]; beta = np.linalg.lstsq(X, C.s.values, rcond=None)[0]
C['r'] = C.s - X @ beta + C.s.mean()  # residual re-centred on the pooled mean (keeps scores positive for ratios)
R_qc = C.groupby(['patient', 'tp']).r.mean().unstack()

def evaluate(W):  # W: index patient, columns pre/post
    rec = []
    for p in W.index:
        b, o = W.loc[p, 'pre'], W.loc[p, 'post']; others = W.drop(p)
        pred = dict(B0=b, B1=others.post.mean(), M2=b * F_LATE)
        rec.append({k: abs(v - o) / b for k, v in pred.items()} | dict(patient=p, ratio=o / b))
    E = pd.DataFrame(rec); m = E[['B0', 'B1', 'M2']].mean()
    med = E.ratio.median()
    if m.M2 <= 0.9 * m.B0 and m.M2 <= 0.9 * m.B1: v = 'PASS'
    elif med >= 0.90 and m.B0 <= m.M2: v = 'FAIL'
    else: v = 'INCONCLUSIVE'
    rng = np.random.default_rng(0)
    boot = [np.median(rng.choice(E.ratio.values, len(E))) for _ in range(5000)]
    return dict(mean_rel_error=m.round(4).to_dict(), median_ratio=round(med, 3),
                ratio_boot90=[round(float(np.percentile(boot, 5)), 3), round(float(np.percentile(boot, 95)), 3)],
                spearman_base_post=round(float(spearmanr(W.pre, W.post)[0]), 3), verdict=v,
                per_patient=E.round(3).to_dict('records'))
res = dict(primary_monocyte=evaluate(S.pivot(index='patient', columns='tp', values='mono')),
           secondary_tcell=evaluate(S.pivot(index='patient', columns='tp', values='t')),
           qc_residualised_monocyte=evaluate(R_qc))
res['qc_flip'] = res['qc_residualised_monocyte']['verdict'] != res['primary_monocyte']['verdict']
res['final_verdict'] = 'INCONCLUSIVE' if res['qc_flip'] else res['primary_monocyte']['verdict']
json.dump(res, open(os.path.join(out, 'result.json'), 'w'), indent=1)
for k in ['primary_monocyte', 'secondary_tcell', 'qc_residualised_monocyte']:
    print(k, {a: b for a, b in res[k].items() if a != 'per_patient'})
    print(pd.DataFrame(res[k]['per_patient']).to_string())
print('FINAL', res['final_verdict'], 'qc_flip', res['qc_flip'])
