"""Audit r145: re-train VC-1 J2 for the FIRST LOPO fold only (P1), exactly as run_vc1.main() does it
(same rng seed 142 consumed in the same order, so fold 1 is reproducible on its own), and compare with
the committed vc1_per_patient.csv row. Full 11-fold re-run (70 min) not repeated. Audit only; changes nothing."""
import os, sys, json, time
import numpy as np, pandas as pd
VC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'virtual-cell'); sys.path.insert(0, VC)
import run_vc1 as R
t0 = time.time()
U, X, lib, cells_of = R.load(None); S = R.Store(U, X, lib, cells_of, None)
E = U[U.eligible].copy(); pats = sorted(E.patient.unique(), key=lambda p: int(p[1:])); cts = sorted(E.ct.unique())
rng = np.random.default_rng(142); p = pats[0]; train = [q for q in pats if q != p]; tgt = 'M3'
scored = [ct for ct in cts if R.unit_id(E, p, ct, 'Base') is not None and R.unit_id(E, p, ct, tgt) is not None]
L = R.latent_fold(S, E, train, cts, rng, augment=True, tune=True, seeds=R.SEEDS)
Ln = R.latent_fold(S, E, train, cts, rng, augment=False, tune=False, seeds=R.SEEDS)
res = {m: [] for m in ['B2', 'J1', 'J2', 'J2noaug', 'J2orig']}
for ct in scored:
    ub, um = R.unit_id(E, p, ct, 'Base'), R.unit_id(E, p, ct, tgt); x, y = S.get(ub), S.get(um)
    sh, mu = R.baseline_stats(S, E, train, ct, tgt); b2 = x if mu is None else mu
    j1, j2, j1o, j2o = R.latent_predict(L, S, E, ub, ct, cts, sh); _, j2n, _, _ = R.latent_predict(Ln, S, E, ub, ct, cts, sh)
    for m, v in dict(B2=b2, J1=j1, J2=j2, J2noaug=j2n, J2orig=j2o).items(): res[m].append(float(np.mean((v - y) ** 2)))
new = {m: float(np.mean(v)) for m, v in res.items()}
old = pd.read_csv(os.path.join(VC, 'vc1_per_patient.csv')).set_index('patient').loc[p]
out = dict(patient=p, epochs_new=int(L['epochs']), epochs_committed=int(old.epochs),
           new=new, committed={m: float(old[m]) for m in new}, rel_diff={m: new[m] / float(old[m]) - 1 for m in new},
           runtime_s=round(time.time() - t0))
sys.path.insert(0, os.path.join(VC, '..', '..', 'scripts', 'provenance')); from stamp import stamp
json.dump(stamp(out, prereg='experiments/virtual-cell/PREREG-VC1.md'), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'vc1_fold1_check.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
