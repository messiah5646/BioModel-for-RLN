"""Round 73 -- prereg: project doc claude/prereg-round73.md (not blind for H1/H2; see that doc)."""
import json, numpy as np
from model import make_params, simulate, obs_vector, VISIT_GRID
from r40_feedback import PARS, NOISE, CHANNELS
SEV = [1, 2, 4, 6, 10, 14, 20, 30]; OBS = CHANNELS['S0']; FULL = PARS + ['G0']; K = FULL.index('k_dam')
def setup(sev):
    p = make_params(0.0); p['k_dam'] *= sev; p['k_tox'] *= sev; p['dep_kill'] = 0.25; return p
def fim(p):
    cols = []
    for nm in FULL:
        h = abs(p[nm]) * 1e-2; pp = dict(p); pp[nm] += h; pm = dict(p); pm[nm] -= h
        cols.append((obs_vector(pp, VISIT_GRID, OBS) - obs_vector(pm, VISIT_GRID, OBS)) / (2 * h) * p[nm])
    J = np.array(cols).T; sig = np.repeat(NOISE[OBS], len(VISIT_GRID)); return J.T @ np.diag(1 / sig ** 2) @ J
def marg(F, names):
    idx = [FULL.index(n) for n in names]; S = F[np.ix_(idx, idx)]
    rk = int(np.linalg.matrix_rank(S, tol=np.linalg.eigvalsh(S).max() * 1e-12))
    return float(np.sqrt(np.linalg.pinv(S, rcond=1e-12)[names.index('k_dam'), names.index('k_dam')])), rk == len(names)
def expo(D, S):
    b, a = np.polyfit(np.log(D), np.log(S), 1); pr = a + b * np.log(D); y = np.log(S)
    return float(b), float(1 - ((y - pr) ** 2).sum() / ((y - y.mean()) ** 2).sum())
rows = []
for sev in SEV:
    p = setup(sev); F = fim(p); D = float(simulate(p, VISIT_GRID)[-1, 5])
    s_alone = float(1 / np.sqrt(F[K, K])); s_known, ok1 = marg(F, PARS); s_free, ok2 = marg(F, FULL)
    loo = {}
    for nm in FULL:
        if nm == 'k_dam': continue
        s, _ = marg(F, [n for n in FULL if n != nm]); loo[nm] = s_free / s          # factor by which sigma drops if nm were known
    top = max(loo, key=loo.get)
    rows.append(dict(sev=sev, D=D, s_alone=s_alone, s_G0known=s_known, s_G0free=s_free, ratio=s_free / s_known, full_rank=bool(ok1 and ok2), loo_top=top, loo={k: round(v, 3) for k, v in loo.items()}))
    print(f"sev {sev:2d} D {D:.4f} alone {s_alone:9.3f} G0known {s_known:9.3f} G0free {s_free:9.3f} ratio {s_free/s_known:5.2f} top {top} ({loo[top]:.2f}x) k_tox {loo['k_tox']:.2f}x G0 {loo['G0']:.2f}x rank_ok {ok1 and ok2}")
D = np.array([r['D'] for r in rows]); bk, r2k = expo(D, [r['s_G0known'] for r in rows]); bf, r2f = expo(D, [r['s_G0free'] for r in rows])
v = dict(exp_G0known=bk, r2_G0known=r2k, exp_G0free=bf, r2_G0free=r2f,
         H1=all(r['ratio'] >= 1.5 for r in rows), H2=abs(bf - bk) < 0.3, H3=r2f >= 0.99,
         H4=all(r['loo_top'] == 'G0' for r in rows if r['sev'] <= 6) and all(r['loo_top'] == 'k_tox' for r in rows if r['sev'] >= 14))
print(v); json.dump(dict(rows=rows, verdict=v), open('results/r73_g0_assumption.json', 'w'), indent=1)
