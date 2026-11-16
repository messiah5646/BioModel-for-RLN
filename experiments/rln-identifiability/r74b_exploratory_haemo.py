"""Round 74b -- EXPLORATORY, not pre-registered. Unfavourable case: the leak arm also carries a reversible
haemodynamic eGFR offset h_drug = 0.10, and both leak_red and h_drug are estimated. Also G0 estimated."""
import json, numpy as np
from model import make_params, obs_vector, VISIT_GRID
from r40_feedback import PARS, NOISE, CHANNELS
OBS = CHANNELS['S0']; NM = PARS + ['G0', 'leak_red', 'h_drug']; K = NM.index('k_dam')
def setup(sev, L, H):
    p = make_params(0.0); p['k_dam'] *= sev; p['k_tox'] *= sev; p['dep_kill'] = 0.25; p['leak_red'] = L; p['h_drug'] = H; p['drug_on'] = 0.0; p['drug_off'] = 1e9; return p
def fim(p):
    cols = []
    for nm in NM:
        ab = nm in ('leak_red', 'h_drug')
        if ab and p[nm] == 0: cols.append(np.zeros(len(VISIT_GRID) * 2)); continue
        h = 0.005 if ab else abs(p[nm]) * 1e-2; pp = dict(p); pp[nm] += h; pm = dict(p); pm[nm] -= h
        cols.append((obs_vector(pp, VISIT_GRID, OBS) - obs_vector(pm, VISIT_GRID, OBS)) / (2 * h) * (1 if ab else p[nm]))
    J = np.array(cols).T; s = np.repeat(NOISE[OBS], len(VISIT_GRID)); return J.T @ np.diag(1 / s ** 2) @ J
def sg(F, names):
    i = [NM.index(x) for x in names]; S = F[np.ix_(i, i)]; return float(np.sqrt(np.linalg.pinv(S, rcond=1e-12)[names.index('k_dam'), names.index('k_dam')]))
out = []
for sev in (1, 6, 20):
    Fc = fim(setup(sev, 0, 0)); Fl = fim(setup(sev, 0.3, 0.10))
    cc = sg(2 * Fc, PARS + ['G0']); a = sg(Fc + Fl, PARS + ['G0', 'leak_red', 'h_drug']); b = sg(Fc + Fl, PARS + ['G0'])
    out.append(dict(sev=sev, s_CC_G0free=cc, gain_all_unknown=cc / a, gain_drug_known=cc / b)); print(out[-1])
json.dump(out, open('results/r74b_exploratory_haemo.json', 'w'), indent=1)
