"""Round 72 -- prereg in project doc claude/prereg-round72.md (written before this ran)."""
import json, numpy as np
from model import make_params, simulate, observe, VISIT_GRID, BASELINE, Y0
from r40_feedback import PARS, NOISE, CHANNELS
OBS = CHANNELS['S0']
def setup(sev):
    p = make_params(0.0); p['k_dam'] *= sev; p['k_tox'] *= sev; p['dep_kill'] = 0.25; return p
def ov(p, D0):
    y0 = Y0.copy(); y0[5] = D0; Y = simulate(p, VISIT_GRID, y0=y0)
    O = np.array([observe(y, p, t) for t, y in zip(VISIT_GRID, Y)]) / BASELINE
    return O[:, OBS].T.ravel()
def fim(p, D0, names):
    cols = []
    for nm in names:
        if nm == 'D0':
            h = 0.01
            d = (ov(p, D0 + h) - ov(p, D0)) / h if D0 == 0 else (ov(p, D0 + h) - ov(p, D0 - h)) / (2 * h)
        else:
            h = abs(p[nm]) * 1e-2; pp = dict(p); pp[nm] += h; pm = dict(p); pm[nm] -= h
            d = (ov(pp, D0) - ov(pm, D0)) / (2 * h) * p[nm]
        cols.append(d)
    J = np.array(cols).T
    # CV noise: sd proportional to the predicted value (r40 used baseline-standardised constant sd; kept identical)
    sig = np.repeat(NOISE[OBS], len(VISIT_GRID)); return J.T @ np.diag(1 / sig ** 2) @ J
def marg(F, names, target):
    w, _ = np.linalg.eigh(F); cond = float(w.max() / max(w.min(), 1e-300))
    C = np.linalg.pinv(F, rcond=1e-12); i = names.index(target)
    return float(np.sqrt(C[i, i])), cond, int(np.linalg.matrix_rank(F, tol=w.max() * 1e-12))
rows = []
for sev in (1, 6, 20):
    for D0 in (0.0, 0.2, 0.4):
        p = setup(sev); full = PARS + ['G0', 'D0']; F = fim(p, D0, full)
        def sub(names, chron=None):
            idx = [full.index(n) for n in names]; S = F[np.ix_(idx, idx)].copy()
            if chron: S[names.index('D0'), names.index('D0')] += 1 / chron ** 2
            return S
        r = dict(sev=sev, D0=D0, D_final=float(simulate(p, VISIT_GRID, y0=np.r_[Y0[:5], D0])[-1, 5]))
        for tag, names, ch in (('r', PARS, None), ('r2', PARS + ['G0'], None), ('a', PARS + ['D0'], None), ('b', full, None), ('c', full, 0.10)):
            S = sub(names, ch); s, cond, rk = marg(S, names, 'k_dam'); r['kdam_' + tag] = s; r['rank_' + tag] = f"{rk}/{len(names)}"; r['cond_' + tag] = cond
            if 'D0' in names: r['sD0_' + tag] = marg(S, names, 'D0')[0]
        r['infl_b'] = r['kdam_b'] / r['kdam_r2']; r['infl_c'] = r['kdam_c'] / r['kdam_r2']
        r['recovered'] = 1 - np.log(r['infl_c']) / np.log(r['infl_b']) if r['infl_b'] > 1.0001 else None
        # sd of chronicity observation needed to recover half
        need = None
        for sd in np.geomspace(0.5, 0.001, 60):
            S = sub(full, sd); ic = marg(S, full, 'k_dam')[0] / r['kdam_r2']
            if r['infl_b'] > 1.0001 and 1 - np.log(ic) / np.log(r['infl_b']) >= 0.5: need = float(sd); break
        r['sd_for_half'] = need; rows.append(r)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if not k.startswith('cond')}, flush=True)
H0 = all(r['sD0_a'] < 0.05 for r in rows); H1 = all(r['infl_b'] >= 2 for r in rows); H1weak = all(r['infl_b'] >= 1.2 for r in rows)
H2 = all((r['recovered'] or 0) >= 0.5 for r in rows); H3 = all(r['sD0_b'] >= 0.10 for r in rows if r['sev'] == 1)
v = dict(H0=H0, H1=H1, H1_at_least_weak=H1weak, H2=H2, H3=H3); print(v)
json.dump(dict(rows=rows, verdict=v), open('results/r72_baseline_damage.json', 'w'), indent=1)
