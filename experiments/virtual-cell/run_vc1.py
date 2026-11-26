"""VC-1 (round 142): LUMEN-JEPA virtual-cell proof of concept, scored under PREREG-VC1.md.

  python3 run_vc1.py --synthetic   code test: every D28/M3 value is replaced by a transform of the patient's own Base
                                   cells (no real post-treatment value is read). Allowed before the 1-h mark.
  python3 run_vc1.py               real run (>= 1 h after the PREREG commit, rule 12b)
"""
import os, sys, json, argparse, time
import numpy as np, pandas as pd, scipy.sparse as sp
import torch, torch.nn as nn
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge

ROOT = os.path.dirname(os.path.abspath(__file__)); W = os.path.join(ROOT, 'work')
TDAYS = {'Base': 0.0, 'D28': 28.0, 'M3': 90.0}
NPC = 50; NBOOT = 20; NHALF = 10; NOISE = 0.1; PAIRS = 30; SEEDS = 5; EPOCH_GRID = (50, 100, 200)
torch.set_num_threads(2)


# ---------------------------------------------------------------- data
def load(synthetic):
    U = pd.read_csv(os.path.join(W, 'units.csv'))
    C = pd.read_csv(os.path.join(W, 'cells.csv'))
    X = sp.load_npz(os.path.join(W, 'cells_hvg.npz')).tocsr()
    U['ct'] = U.compartment + ':' + U.cell_type
    U['dose_num'] = U.dose.str.extract(r'([\d.]+)')[0].astype(float)
    cells_of = C.groupby('unit').indices
    lib = C.lib.to_numpy()
    rng = np.random.default_rng(142)
    if synthetic:  # fake post-treatment units from the same patient's Base cells; real D28/M3 cells never used
        base_cells = {(r.patient, r.ct): cells_of[r.unit] for r in U[U.timepoint == 'Base'].itertuples()}
        fake = {}
        for r in U[U.timepoint != 'Base'].itertuples():
            bc = base_cells.get((r.patient, r.ct))
            fake[r.unit] = None if bc is None else rng.choice(bc, size=len(bc), replace=True)
        for u, idx in fake.items():
            if idx is None: U.loc[U.unit == u, 'eligible'] = False
            else: cells_of[u] = idx
        U.loc[U.unit.isin(fake), 'n_cells'] = [len(cells_of[u]) if u in cells_of else 0 for u in U.unit[U.unit.isin(fake)]]
    return U, X, lib, cells_of


def pb(X, lib, idx, synth_shift=None):
    v = np.log1p(np.asarray(X[idx].sum(0)).ravel() / lib[idx].sum() * 1e4)
    return v if synth_shift is None else v + synth_shift


class Store:
    """Pseudobulks (full and augmented) per unit. In synthetic mode a fixed per-cell-type shift is added to D28/M3."""
    def __init__(s, U, X, lib, cells_of, synthetic):
        s.U, s.X, s.lib, s.cells_of = U.set_index('unit', drop=False), X, lib, cells_of
        s.full, s.rng = {}, np.random.default_rng(142)
        s.shift, s.scale = {}, {}
        if synthetic:
            r = np.random.default_rng(7)
            for ct in U.ct.unique(): s.shift[ct] = r.normal(0, 0.15, X.shape[1])
            if synthetic == 'pos':  # positive control: an interaction a trivial rule cannot copy. The size of the
                # response along baseline program v depends on the patient's own baseline score on a second program w
                # (v, w = top-2 principal axes of Base pseudobulks, so the signal is low-rank and biologically shaped).
                B = U[(U.timepoint == 'Base') & U.eligible]
                PB = np.vstack([pb(X, lib, cells_of[u]) for u in B.unit])
                Pc = PCA(n_components=2, random_state=0).fit(PB); v, w = Pc.components_
                sc = (PB - Pc.mean_) @ w
                for ct, g in B.assign(sc=sc).groupby('ct'):
                    zz = (g.sc - g.sc.mean()) / (g.sc.std() + 1e-9) if len(g) > 1 else g.sc * 0
                    for pt, z in zip(g.patient, zz): s.scale[(pt, ct)] = float(z)
                s.v = v * 4.0
        s.synthetic = synthetic

    def _sh(s, u):
        r = s.U.loc[u]
        if not s.synthetic or r.timepoint == 'Base': return None
        extra = s.v * s.scale.get((r.patient, r.ct), 0.0) if s.synthetic == 'pos' else 0.0
        return (s.shift[r.ct] + extra) * (0.5 if r.timepoint == 'D28' else 1.0)

    def get(s, u):
        if u not in s.full: s.full[u] = pb(s.X, s.lib, s.cells_of[u], s._sh(u))
        return s.full[u]

    def aug(s, u):
        if not hasattr(s, '_aug'): s._aug = {}
        if u in s._aug: return s._aug[u]
        idx = s.cells_of[u]; out = []
        for _ in range(NBOOT): out.append(pb(s.X, s.lib, s.rng.choice(idx, len(idx), replace=True), s._sh(u)))
        for _ in range(NHALF): out.append(pb(s.X, s.lib, s.rng.choice(idx, max(len(idx) // 2, 1), replace=False), s._sh(u)))
        s._aug[u] = np.vstack(out); return s._aug[u]

    def halves(s, u, n=20):
        idx = s.cells_of[u]; m = []
        for _ in range(n):
            p = s.rng.permutation(idx); h = len(p) // 2
            m.append(np.mean((pb(s.X, s.lib, p[:h], s._sh(u)) - pb(s.X, s.lib, p[h:], s._sh(u))) ** 2))
        return float(np.mean(m)) / 4.0


# ---------------------------------------------------------------- baselines
def unit_id(E, p, ct, t):
    r = E[(E.patient == p) & (E.ct == ct) & (E.timepoint == t)]
    return None if r.empty else int(r.unit.iloc[0])


def baseline_stats(S, E, pats, ct, tgt='M3'):
    deltas, tgts = [], []
    for p in pats:
        b, m = unit_id(E, p, ct, 'Base'), unit_id(E, p, ct, tgt)
        if m is not None: tgts.append(S.get(m))
        if b is not None and m is not None: deltas.append(S.get(m) - S.get(b))
    return (np.mean(deltas, 0) if deltas else None), (np.mean(tgts, 0) if tgts else None)


def fit_a(S, E, pats, ct, tgt='M3'):
    grid = np.round(np.arange(0, 1.0001, 0.05), 2); err = np.zeros(len(grid)); n = 0
    for q in pats:
        b, m = unit_id(E, q, ct, 'Base'), unit_id(E, q, ct, tgt)
        if b is None or m is None: continue
        sh, mu = baseline_stats(S, E, [x for x in pats if x != q], ct, tgt)
        if sh is None or mu is None: continue
        x, y = S.get(b), S.get(m); n += 1
        for i, a in enumerate(grid): err[i] += np.mean((a * (x + sh) + (1 - a) * mu - y) ** 2)
    return 1.0 if n == 0 else float(grid[np.argmin(err)])


# ---------------------------------------------------------------- latent models
class MLP(nn.Module):
    def __init__(s, d_in, d_out):
        super().__init__()
        s.f = nn.Sequential(nn.Linear(d_in, 128), nn.ReLU(), nn.Dropout(0.2), nn.Linear(128, 128), nn.ReLU(),
                            nn.Dropout(0.2), nn.Linear(128, d_out))
    def forward(s, x): return s.f(x)


def make_pairs(S, E, pats, cts, pca, augment, rng, pairs=PAIRS):
    src, dst, feat = [], [], []
    trans = [('Base', 'D28'), ('D28', 'M3'), ('Base', 'M3')] if augment else [('Base', 'M3')]
    ctix = {c: i for i, c in enumerate(cts)}
    for p in pats:
        for ct in cts:
            for a, b in trans:
                ua, ub = unit_id(E, p, ct, a), unit_id(E, p, ct, b)
                if ua is None or ub is None: continue
                if augment:
                    A, B = S.aug(ua), S.aug(ub)
                    ia, ib = rng.integers(0, len(A), pairs), rng.integers(0, len(B), pairs)
                    za, zb = pca.transform(A[ia]), pca.transform(B[ib])
                else:
                    za, zb = pca.transform(S.get(ua)[None]), pca.transform(S.get(ub)[None])
                for k in range(len(za)):
                    oh = np.zeros(len(cts)); oh[ctix[ct]] = 1
                    meta = [E.loc[E.unit == ua, 'dose_num'].iloc[0], TDAYS[a] / 90, (TDAYS[b] - TDAYS[a]) / 90]
                    src.append(za[k]); dst.append(zb[k]); feat.append(np.concatenate([oh, meta]))
    return np.array(src), np.array(dst), np.array(feat)


def train_mlp(zs, zt, f, sd, epochs, seed, augment, comps, checkpoints=None):
    torch.manual_seed(seed); np.random.seed(seed)
    Xin = np.hstack([zs / sd, f]).astype(np.float32); Y = ((zt - zs) / sd).astype(np.float32)
    net = MLP(Xin.shape[1], Y.shape[1]); opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=1e-3)
    Xt, Yt = torch.tensor(Xin), torch.tensor(Y); Ct = torch.tensor(comps.astype(np.float32)); sdt = torch.tensor(sd.astype(np.float32))
    nz = zs.shape[1]; snaps = {}
    for ep in range(1, epochs + 1):
        net.train(); perm = torch.randperm(len(Xt))
        for i in range(0, len(Xt), 256):
            j = perm[i:i + 256]; x = Xt[j].clone()
            if augment: x[:, :nz] += NOISE * torch.randn_like(x[:, :nz])
            pred = net(x); lat = ((pred - Yt[j]) ** 2).mean()
            gene = ((((pred - Yt[j]) * sdt) @ Ct) ** 2).mean()
            loss = lat + 0.1 * gene; opt.zero_grad(); loss.backward(); opt.step()
        if checkpoints and ep in checkpoints: snaps[ep] = {k: v.clone() for k, v in net.state_dict().items()}
    return net, snaps


def predict_mlp(net, z, f, sd):
    net.eval()
    with torch.no_grad():
        return z + net(torch.tensor(np.hstack([z / sd, f]).astype(np.float32))).numpy() * sd


# ---------------------------------------------------------------- one fold
def latent_fold(S, E, train, cts, rng, augment=True, tune=True, seeds=SEEDS, inner=True):
    units = [u for u in E.unit if E.loc[E.unit == u, 'patient'].iloc[0] in train]
    P = np.vstack([S.get(u) for u in units])
    pca = PCA(n_components=min(NPC, len(units) - 1), random_state=0).fit(P)
    zs, zt, f = make_pairs(S, E, train, cts, pca, augment, rng)
    sd = zs.std(0) + 1e-6
    mu = mean_delta(zs, zt, f, len(cts)); zt_res = zt - mu(f)
    ridge = Ridge(alpha=1.0).fit(np.hstack([zs / sd, f]), (zt_res - zs) / sd)
    ridge_orig = Ridge(alpha=1.0).fit(np.hstack([zs / sd, f]), (zt - zs) / sd)
    ep = 100
    if tune and inner:
        errs = {e: [] for e in EPOCH_GRID}
        for q in train:
            tr = [x for x in train if x != q]
            a, b, c = make_pairs(S, E, tr, cts, pca, augment, rng, pairs=10)
            if len(a) == 0: continue
            mi = mean_delta(a, b, c, len(cts))
            net, snaps = train_mlp(a, b - mi(c), c, sd, max(EPOCH_GRID), 0, augment, pca.components_, set(EPOCH_GRID))
            shi = {ct: baseline_stats(S, E, tr, ct)[0] for ct in cts}
            for ct in cts:
                ub, um = unit_id(E, q, ct, 'Base'), unit_id(E, q, ct, 'M3')
                if ub is None or um is None: continue
                z = pca.transform(S.get(ub)[None]); fe = feat_row(E, ub, ct, cts)
                for e in EPOCH_GRID:
                    net.load_state_dict(snaps[e])
                    dz = predict_mlp(net, z, fe, sd) - z
                    sh0 = shi[ct] if shi[ct] is not None else 0.0
                    errs[e].append(np.mean((S.get(ub) + sh0 + (dz @ pca.components_)[0] - S.get(um)) ** 2))
        ep = min(EPOCH_GRID, key=lambda e: np.mean(errs[e]) if errs[e] else np.inf)
    nets = [train_mlp(zs, zt_res, f, sd, ep, s, augment, pca.components_)[0] for s in range(seeds)]
    nets_orig = [train_mlp(zs, zt, f, sd, ep, s, augment, pca.components_)[0] for s in range(seeds)]
    return dict(pca=pca, sd=sd, ridge=ridge, ridge_orig=ridge_orig, nets=nets, nets_orig=nets_orig, epochs=ep, n_pairs=len(zs))


def mean_delta(zs, zt, f, nct):
    """Mean latent delta per (cell type, transition) over the training pairs; returns a function of the feature rows."""
    key = lambda F: [(int(np.argmax(r[:nct])), round(float(r[nct + 1]), 4), round(float(r[nct + 2]), 4)) for r in F]
    k = key(f); d = zt - zs; acc = {}
    for i, kk in enumerate(k): acc.setdefault(kk, []).append(d[i])
    m = {kk: np.mean(v, 0) for kk, v in acc.items()}
    return lambda F: np.vstack([m.get(kk, np.zeros(zs.shape[1])) for kk in key(F)])


def feat_row(E, ub, ct, cts):
    oh = np.zeros(len(cts)); oh[cts.index(ct)] = 1
    return np.concatenate([oh, [E.loc[E.unit == ub, 'dose_num'].iloc[0], 0.0, 1.0]])[None]


def latent_predict(L, S, E, ub, ct, cts, shift):
    """Amendment A1: J1/J2 = Base + B1 mean shift + decoded patient-specific latent deviation.
    J1orig/J2orig = the PREREG-as-written form (Base + decoded latent delta, no mean-shift residual)."""
    x = S.get(ub); z = L['pca'].transform(x[None]); f = feat_row(E, ub, ct, cts); sd = L['sd']; C = L['pca'].components_
    sh0 = shift if shift is not None else 0.0
    j1 = x + sh0 + ((L['ridge'].predict(np.hstack([z / sd, f])) * sd) @ C)[0]
    j2 = x + sh0 + ((np.mean([predict_mlp(n, z, f, sd) for n in L['nets']], 0) - z) @ C)[0]
    j1o = x + ((L['ridge_orig'].predict(np.hstack([z / sd, f])) * sd) @ C)[0]
    j2o = x + ((np.mean([predict_mlp(n, z, f, sd) for n in L['nets_orig']], 0) - z) @ C)[0]
    return j1, j2, j1o, j2o


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--synthetic', choices=['neg', 'pos'], default=None)
    ap.add_argument('--quick', action='store_true', help='code test only: 1 seed, no inner tuning, 3 folds')
    ap.add_argument('--target', default='M3'); args = ap.parse_args()
    t0 = time.time()
    U, X, lib, cells_of = load(args.synthetic)
    S = Store(U, X, lib, cells_of, args.synthetic)
    E = U[U.eligible].copy()
    pats = sorted(E.patient.unique(), key=lambda p: int(p[1:]))
    cts = sorted(E.ct.unique())
    rng = np.random.default_rng(142)
    rows, units_out = [], []
    folds = pats[:3] if args.quick else pats
    for p in folds:
        train = [q for q in pats if q != p]
        scored = [ct for ct in cts if unit_id(E, p, ct, 'Base') is not None and unit_id(E, p, ct, args.target) is not None]
        L = latent_fold(S, E, train, cts, rng, augment=True, tune=not args.quick, seeds=1 if args.quick else SEEDS)
        Ln = latent_fold(S, E, train, cts, rng, augment=False, tune=False, seeds=1 if args.quick else SEEDS)
        res = {m: [] for m in ['B0', 'B1', 'B2', 'B3', 'J1', 'J2', 'J2noaug', 'J1orig', 'J2orig', 'ceiling']}
        for ct in scored:
            ub, um = unit_id(E, p, ct, 'Base'), unit_id(E, p, ct, args.target)
            x, y = S.get(ub), S.get(um)
            sh, mu = baseline_stats(S, E, train, ct, args.target)
            b1 = x if sh is None else x + sh; b2 = x if mu is None else mu
            a = fit_a(S, E, train, ct, args.target) if sh is not None and mu is not None else 1.0
            b3 = a * b1 + (1 - a) * b2
            j1, j2, j1o, j2o = latent_predict(L, S, E, ub, ct, cts, sh)
            _, j2n, _, _ = latent_predict(Ln, S, E, ub, ct, cts, sh)
            preds = dict(B0=x, B1=b1, B2=b2, B3=b3, J1=j1, J2=j2, J2noaug=j2n, J1orig=j1o, J2orig=j2o)
            mse = {m: float(np.mean((v - y) ** 2)) for m, v in preds.items()}
            mse['ceiling'] = S.halves(um)
            d = y - x
            dc = {m: (float(np.corrcoef(v - x, d)[0, 1]) if np.std(v - x) > 0 else np.nan) for m, v in preds.items() if m != 'B0'}
            units_out.append(dict(patient=p, ct=ct, n_base=int(E.loc[E.unit == ub, 'n_cells'].iloc[0]),
                                  n_tgt=int(E.loc[E.unit == um, 'n_cells'].iloc[0]), a=a, no_shift=sh is None,
                                  **{f'mse_{m}': v for m, v in mse.items()}, **{f'dcorr_{m}': v for m, v in dc.items()}))
            for m in res: res[m].append(mse[m])
        rows.append(dict(patient=p, n_scored=len(scored), epochs=L['epochs'], n_pairs=L['n_pairs'],
                         **{m: float(np.mean(v)) if v else np.nan for m, v in res.items()}))
        print(p, len(scored), {m: round(float(np.mean(v)), 5) for m, v in res.items() if v}, f'{time.time()-t0:.0f}s', flush=True)
    R = pd.DataFrame(rows); D = pd.DataFrame(units_out)
    tag = (f'synthetic{args.synthetic}_' if args.synthetic else '') + ('quick_' if args.quick else '') + (args.target if args.target != 'M3' else '')
    R.to_csv(os.path.join(ROOT, f'vc1_per_patient{"_" + tag if tag else ""}.csv'), index=False)
    D.to_csv(os.path.join(ROOT, f'vc1_per_unit{"_" + tag if tag else ""}.csv'), index=False)
    out = summarise(R, D)
    out.update(mode=f'synthetic-{args.synthetic}' if args.synthetic else 'real', quick=args.quick, target=args.target,
               n_patients=int(len(R)), cell_types=len(cts), runtime_s=round(time.time() - t0))
    fn = os.path.join(ROOT, f'vc1_result{"_" + tag if tag else ""}.json')
    sys.path.insert(0, os.path.join(ROOT, '..', '..', 'scripts', 'provenance')); from stamp import stamp
    json.dump(stamp(out, prereg='experiments/virtual-cell/PREREG-VC1.md'), open(fn, 'w'), indent=1)
    print(json.dumps(out, indent=1))


def summarise(R, D):
    A = {m: float(R[m].mean()) for m in ['B0', 'B1', 'B2', 'B3', 'J1', 'J2', 'J2noaug', 'J1orig', 'J2orig', 'ceiling']}
    base = ['B0', 'B1', 'B2', 'B3']; best = min(base, key=lambda m: A[m])
    out = dict(A=A, best_baseline=best)
    elig = R[R.n_scored >= 3]
    out['n_eligible_patients'] = int(len(elig))
    for J in ['J2', 'J1']:
        wins = int((R[J] < R[best]).sum())
        ratio = A[J] / A[best]
        fail = any(A[J] >= 1.10 * A[m] for m in base)
        if len(elig) < 9: v = 'UNDERPOWERED'
        elif fail: v = 'FAIL'
        elif ratio <= 0.90 and wins >= 9: v = 'PASS'
        else: v = 'INCONCLUSIVE'
        out[J] = dict(ratio_to_best=ratio, wins_vs_best=wins, verdict=v)
    out['verdict'] = out['J2']['verdict'] if out['J2']['verdict'] == 'PASS' else (
        'PASS (J1 fallback)' if out['J1']['verdict'] == 'PASS' else out['J2']['verdict'])
    out['dcorr_median'] = {m: float(D[f'dcorr_{m}'].median()) for m in ['B1', 'B3', 'J1', 'J2', 'J2noaug', 'J2orig']}
    num = A['B1'] - np.array([A[m] for m in ['J2', 'J1', 'B3']]); den = A['B1'] - A['ceiling']
    out['vcc_score'] = dict(zip(['J2', 'J1', 'B3'], (num / den).tolist())) if den > 0 else None
    rng = np.random.default_rng(142); r = []
    for _ in range(2000):
        s = R.sample(len(R), replace=True, random_state=int(rng.integers(1e9)))
        r.append(s['J2'].mean() / s[best].mean())
    out['S7_boot_ratio_J2_best_95'] = [float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5))]
    out['S4_by_compartment'] = {c: {m: float(g[f'mse_{m}'].mean()) for m in ['B0', 'B1', 'B3', 'J2']}
                                for c, g in D.assign(comp=D.ct.str.split(':').str[0]).groupby('comp')}
    return out


if __name__ == '__main__':
    main()
