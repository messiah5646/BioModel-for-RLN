"""VC-1 descriptive follow-up (round 143, reviewer R-20261126-1). POST HOC, DESCRIPTIVE: cannot change the VC-1 FAIL.

Same pseudobulks, patients, eligible units and LOPO scoring as run_vc1.py (B0-B3 reuse its functions). Adds:
  B2base  training patients' mean BASE pseudobulk used as the prediction of the held-out patient's target.
          If B2base ~ B2, the 'population M3 average' wins because a pooled profile beats any single other library,
          not because month 3 is a shared post-therapy state.
  S3      the same baselines with D28 as the target (pre-registered in PREREG-VC1, not run in round 142).
  disp    noise-corrected between-patient dispersion per cell type at Base, D28, M3, on the SAME patients at each
          compared pair of timepoints: mean over genes of between-patient variance minus mean per-unit sampling
          variance (half-split estimate, as run_vc1.Store.halves).
Every block is run on all 11 patients and without P9 (P9's 'M3' library is named S14_M13 in obs.sample).

    python3 vc1_descriptive.py        (needs work/ from build_pseudobulk.py)
"""
import os, sys, json
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, ROOT)
from run_vc1 import load, Store, unit_id, baseline_stats, fit_a


def mean_of(S, E, pats, ct, t):
    v = [S.get(u) for u in (unit_id(E, p, ct, t) for p in pats) if u is not None]
    return np.mean(v, 0) if v else None


def lopo(S, E, pats, cts, tgt):
    rows, units = [], []
    for p in pats:
        train = [q for q in pats if q != p]
        res = {m: [] for m in ['B0', 'B1', 'B2', 'B2base', 'B3']}
        for ct in cts:
            ub, um = unit_id(E, p, ct, 'Base'), unit_id(E, p, ct, tgt)
            if ub is None or um is None: continue
            x, y = S.get(ub), S.get(um)
            sh, mu = baseline_stats(S, E, train, ct, tgt)
            mb = mean_of(S, E, train, ct, 'Base')
            b1 = x if sh is None else x + sh; b2 = x if mu is None else mu
            a = fit_a(S, E, train, ct, tgt) if sh is not None and mu is not None else 1.0
            pr = dict(B0=x, B1=b1, B2=b2, B2base=x if mb is None else mb, B3=a * b1 + (1 - a) * b2)
            for m, v in pr.items(): res[m].append(float(np.mean((v - y) ** 2)))
            units.append(dict(patient=p, ct=ct, **{m: res[m][-1] for m in res}))
        if res['B0']: rows.append(dict(patient=p, n=len(res['B0']), **{m: float(np.mean(v)) for m, v in res.items()}))
    R = pd.DataFrame(rows)
    A = {m: float(R[m].mean()) for m in ['B0', 'B1', 'B2', 'B2base', 'B3']}
    rng = np.random.default_rng(143); bt = []
    for _ in range(2000):
        s = R.sample(len(R), replace=True, random_state=int(rng.integers(1e9))); bt.append(s.B2base.mean() / s.B2.mean())
    return dict(n_patients=int(len(R)), A=A, ratio_B2base_over_B2=A['B2base'] / A['B2'],
                boot95_B2base_over_B2=[float(np.percentile(bt, 2.5)), float(np.percentile(bt, 97.5))],
                wins_B2_over_B2base=int((R.B2 < R.B2base).sum()), wins_B2_over_B0=int((R.B2 < R.B0).sum()),
                per_patient=R.round(5).to_dict('records'))


def dispersion(S, E, pats, cts, t1, t2):
    out = []
    for ct in cts:
        both = [p for p in pats if unit_id(E, p, ct, t1) is not None and unit_id(E, p, ct, t2) is not None]
        if len(both) < 3: continue
        d = {}
        for t in (t1, t2):
            us = [unit_id(E, p, ct, t) for p in both]
            M = np.vstack([S.get(u) for u in us])
            d[t] = float(M.var(0, ddof=1).mean() - np.mean([S.halves(u) for u in us]))
        out.append(dict(ct=ct, n=len(both), **d))
    D = pd.DataFrame(out)
    return dict(n_cell_types=int(len(D)), median={t1: float(D[t1].median()), t2: float(D[t2].median())},
                n_lower_at_later=int((D[t2] < D[t1]).sum()), per_ct=D.round(5).to_dict('records'))


def main():
    U, X, lib, cells_of = load(None)
    S = Store(U, X, lib, cells_of, None)
    E = U[U.eligible].copy()
    allp = sorted(E.patient.unique(), key=lambda p: int(p[1:]))
    cts = sorted(E.ct.unique())
    p9_sample = None
    out = dict(note='POST HOC, DESCRIPTIVE (reviewer R-20261126-1); cannot change the VC-1 FAIL', blocks={})
    for tag, pats in [('all', allp), ('noP9', [p for p in allp if p != 'P9'])]:
        b = {}
        for tgt in ['M3', 'D28']: b[f'lopo_{tgt}'] = lopo(S, E, pats, cts, tgt)
        S.rng = np.random.default_rng(143)
        b['disp_Base_M3'] = dispersion(S, E, pats, cts, 'Base', 'M3')
        b['disp_Base_D28'] = dispersion(S, E, pats, cts, 'Base', 'D28')
        b['disp_D28_M3'] = dispersion(S, E, pats, cts, 'D28', 'M3')
        out['blocks'][tag] = b
        print(tag, {k: (v['A'] if 'A' in v else v['median']) for k, v in b.items()}, flush=True)
    sys.path.insert(0, os.path.join(ROOT, '..', '..', 'scripts', 'provenance')); from stamp import stamp
    json.dump(stamp(out, prereg='experiments/virtual-cell/PREREG-VC1.md'),
              open(os.path.join(ROOT, 'vc1_descriptive.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
