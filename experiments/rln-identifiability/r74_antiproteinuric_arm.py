"""Round 74 -- prereg: project doc claude/prereg-round74.md (written before this ran; blind)."""
import json, numpy as np
from model import make_params, simulate, obs_vector, VISIT_GRID
from r40_feedback import PARS, NOISE, CHANNELS
SEV = [1, 2, 4, 6, 10, 14, 20, 30]; OBS = CHANNELS['S0']; NM = PARS + ['leak_red']; K = NM.index('k_dam'); n = len(PARS)
def setup(sev, L):
    p = make_params(0.0); p['k_dam'] *= sev; p['k_tox'] *= sev; p['dep_kill'] = 0.25; p['leak_red'] = L; p['drug_on'] = 0.0; p['drug_off'] = 1e9; return p
def fim(p):
    cols = []
    for nm in NM:
        if nm == 'leak_red':
            h = 0.005; pp = dict(p); pp[nm] += h; pm = dict(p); pm[nm] -= h
            if p[nm] == 0: cols.append(np.zeros(len(VISIT_GRID) * len(OBS))); continue
            cols.append((obs_vector(pp, VISIT_GRID, OBS) - obs_vector(pm, VISIT_GRID, OBS)) / (2 * h))
        else:
            h = abs(p[nm]) * 1e-2; pp = dict(p); pp[nm] += h; pm = dict(p); pm[nm] -= h
            cols.append((obs_vector(pp, VISIT_GRID, OBS) - obs_vector(pm, VISIT_GRID, OBS)) / (2 * h) * p[nm])
    J = np.array(cols).T; sig = np.repeat(NOISE[OBS], len(VISIT_GRID)); return J.T @ np.diag(1 / sig ** 2) @ J
def sig(F, m):
    S = F[:m, :m]; w = np.linalg.eigvalsh(S); ok = int(np.linalg.matrix_rank(S, tol=w.max() * 1e-12)) == m
    return float(np.sqrt(np.linalg.pinv(S, rcond=1e-12)[K, K])), ok
rows = []
for sev in SEV:
    Fc = fim(setup(sev, 0.0)); FL = {L: fim(setup(sev, L)) for L in (0.1, 0.3, 0.5)}
    cc, o1 = sig(2 * Fc, n); cl, o2 = sig(Fc + FL[0.3], n); clu, o3 = sig(Fc + FL[0.3], n + 1); ll, o4 = sig(2 * FL[0.3], n)
    dose = {L: cc / sig(Fc + FL[L], n)[0] for L in FL}
    D = float(simulate(setup(sev, 0.0), VISIT_GRID)[-1, 5]); DL = float(simulate(setup(sev, 0.3), VISIT_GRID)[-1, 5])
    r = dict(sev=sev, D_ctrl=D, D_leak=DL, s_CC=cc, s_CL=cl, s_CLu=clu, s_LL=ll, gain_CL=cc / cl, gain_CLu=cc / clu, gain_LL=cc / ll, dose={str(k): v for k, v in dose.items()}, rank_ok=bool(o1 and o2 and o3 and o4))
    rows.append(r); print(f"sev {sev:2d} D {D:.3f}/{DL:.3f} sCC {cc:8.3f} sCL {cl:7.3f} sCLu {clu:8.3f} sLL {ll:8.3f} | gain CL {cc/cl:6.2f} CLu {cc/clu:6.2f} LL {cc/ll:5.2f} | dose {dose[0.1]:.2f} {dose[0.3]:.2f} {dose[0.5]:.2f} rank {r['rank_ok']}")
v = dict(H1=all(r['gain_CL'] >= 3 for r in rows), H1_dead=any(r['gain_CL'] < 1.5 for r in rows), H2=all(r['gain_CLu'] >= 2 for r in rows),
         H3=all(r['gain_LL'] < 1.5 for r in rows), H4=all(r['dose']['0.1'] < r['dose']['0.3'] < r['dose']['0.5'] for r in rows))
print(v); json.dump(dict(rows=rows, verdict=v), open('results/r74_antiproteinuric_arm.json', 'w'), indent=1)
