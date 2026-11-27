"""O3 (round 152): implements experiments/lumen-v3/PREREG-O3.md (commit 05328df) exactly.

Question: should LUMEN's drug_effect forecast shift by trial phase (phase 2 vs phase 3)?
Run from the repo root:  python3 -I experiments/lumen-v3/run_o3.py
Imports run_rt1.py unchanged (tables, DATE, CFG, late_arms, score_stratum; shipped iv/dl). Research and simulation only.
"""
import sys, os, json, math
import numpy as np, pandas as pd
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'experiments', 'lumen-bench', 'timeforward'))
from lumen import provenance, baselines as B, gate
ok = provenance.verify_all(); assert all(ok.values()), ok
import run_rt1 as R1   # unchanged
OUT = os.path.dirname(os.path.abspath(__file__))
PREREG = 'experiments/lumen-v3/PREREG-O3.md'

PHASES = {}
for s in R1.ct:
    p = s['protocolSection']
    PHASES[p['identificationModule']['nctId']] = p.get('designModule', {}).get('phases') or []

def ph3(nct, recode_23=False):
    """PREREG s5: ph3 = 1 if 'PHASE3' in phases (PHASE2/PHASE3 -> 1). Sensitivity (c): PHASE2/PHASE3 -> 0."""
    ps = PHASES.get(nct, [])
    if not ps: return None
    if recode_23 and 'PHASE2' in ps and 'PHASE3' in ps: return 0
    return 1 if 'PHASE3' in ps else 0

def wls(y, se, x, x0, tau2, re=False):
    """WLS of y on [1, x] (or [1] if x is None). Weights 1/se^2 (P1) or 1/(se^2+tau2) (P1-RE).
    Returns predictive mean, var = tau2 + x0'(X'WX)^-1 x0, coefficients and their covariance."""
    w = 1 / (se ** 2 + (tau2 if re else 0.0))
    X = np.ones((len(y), 1)) if x is None else np.column_stack([np.ones(len(y)), x])
    XtW = X.T * w; cov = np.linalg.inv(XtW @ X); beta = cov @ (XtW @ y)
    x0 = np.array([1.0]) if x is None else np.array([1.0, x0])
    return float(x0 @ beta), float(tau2 + x0 @ cov @ x0), beta, cov

def build(mode, recode_23=False):
    rows = []
    for st, R in [('LN', R1.LN), ('SLE', R1.SLE)]:
        sc = pd.DataFrame(R1.score_stratum(R, st, mode)); sc = sc[sc.subset == 'all']
        trials = R[['trial', 'nct']].drop_duplicates(); trials = trials.assign(date=trials.nct.map(R1.DATE))
        nct_of = dict(zip(trials.trial, trials.nct))
        for t, g in sc.groupby('trial'):
            tr_date = trials[trials.trial == t].date.iloc[0]
            tt = trials[trials.date < tr_date].trial if mode == 'tf' else trials[trials.trial != t].trial
            P = R1.late_arms(R[R.trial.isin(tt)], st).copy()
            P['ph3'] = [ph3(n, recode_23) for n in P.nct]
            x0 = ph3(nct_of[t], recode_23)
            if x0 is None or P.ph3.isna().any(): raise SystemExit('missing phase (none expected; PREREG s5)')
            y, se = P.lor.values.astype(float), P.se.values.astype(float)
            tau2 = R1.dl(y, se)
            # identity check (PREREG s6): intercept-only WLS == shipped M1
            m_chk, v_chk, _, _ = wls(y, se, None, None, tau2)
            assert abs(m_chk - g.M1.iloc[0]) < 1e-10 and abs(v_chk - g.v_M1.iloc[0]) < 1e-10, (t, m_chk, g.M1.iloc[0])
            tph = P.drop_duplicates('trial').ph3
            fallback = bool((tph == 1).sum() < 3 or (tph == 0).sum() < 3)
            if fallback:
                mP, vP, mR, vR, b = m_chk, v_chk, m_chk, v_chk, np.nan
            else:
                mP, vP, beta, _ = wls(y, se, P.ph3.values.astype(float), x0, tau2)
                mR, vR, _, _ = wls(y, se, P.ph3.values.astype(float), x0, tau2, re=True)
                b = float(beta[1])
            for _, r in g.iterrows():
                rows.append(dict(mode=mode, stratum=st, trial=t, arm=r.arm, ph3=x0, obs=r.obs, se_obs=r.se_obs,
                                 M0=0.0, M1=r.M1, v_M1=r.v_M1, P1=mP, v_P1=vP, P1RE=mR, v_P1RE=vR,
                                 b_train=b, fallback=fallback, n_train_trials=int(P.trial.nunique())))
    return pd.DataFrame(rows)

def summ(D, label):
    out = dict(label=label, n_arms=int(len(D)), n_trials=int(D.trial.nunique()),
               n_trials_no_fallback=int(D[~D.fallback].trial.nunique()),
               n_ph3_trials_no_fallback=int(D[(~D.fallback) & (D.ph3 == 1)].trial.nunique()),
               n_ph2_trials_no_fallback=int(D[(~D.fallback) & (D.ph3 == 0)].trial.nunique()))
    if not len(D): return out
    for m, v in [('M0', 'v_M1'), ('M1', 'v_M1'), ('P1', 'v_P1'), ('P1RE', 'v_P1RE')]:
        sd = np.sqrt(D[v].values)
        out[m] = dict(ivmse=round(B.iv_mse(D.obs, D.se_obs, D[m]), 4),
                      logscore=round(B.log_score(D.obs, D.se_obs, D[m], sd), 4),
                      cov95=round(B.coverage(D.obs, D.se_obs, D[m], sd), 3))
    g = gate.verdict(out['P1']['ivmse'], {'M0': out['M0']['ivmse'], 'M1': out['M1']['ivmse']}, reference='M1',
                     n_units=out['n_arms'], n_trials=out['n_trials'], margin=0.10)
    out['gate'] = dict(verdict=g.verdict, reasons=g.reasons)
    out['P1_over_M1'] = round(out['P1']['ivmse'] / out['M1']['ivmse'], 3) if out['M1']['ivmse'] else None
    return out

def verdict(s):
    """PREREG s8-s10: minimum n, then the gate."""
    if s['n_trials_no_fallback'] < 20 or s['n_ph3_trials_no_fallback'] < 6 or s['n_ph2_trials_no_fallback'] < 6:
        return 'UNDERPOWERED'
    return s['gate']['verdict']

def full_fit():
    """Descriptive (f): phase coefficient b (log OR, ph3 minus ph2) on all late arms per stratum."""
    out = {}
    for st, R in [('LN', R1.LN), ('SLE', R1.SLE)]:
        P = R1.late_arms(R, st).copy(); P['ph3'] = [ph3(n) for n in P.nct]
        y, se, x = P.lor.values.astype(float), P.se.values.astype(float), P.ph3.values.astype(float)
        tau2 = R1.dl(y, se); res = dict(n_arms=int(len(P)), n_trials=int(P.trial.nunique()),
                                        n_ph3_trials=int(P.drop_duplicates('trial').ph3.sum()), tau2_DL=round(float(tau2), 4))
        for nm, re in [('FE_weights', False), ('RE_weights', True)]:
            _, _, beta, cov = wls(y, se, x, 0, tau2, re=re); seb = math.sqrt(cov[1, 1])
            res[nm] = dict(intercept_ph2=round(float(beta[0]), 3), b_ph3_minus_ph2=round(float(beta[1]), 3),
                           ci95=[round(float(beta[1] - 1.96 * seb), 3), round(float(beta[1] + 1.96 * seb), 3)])
        out[st] = res
    return out

def freeze_w1():
    """Prospective record (s11): P1 forecasts for the 8 frozen W1-SLE trials, same training rows as O1's freeze."""
    sys.path.insert(0, os.path.join(ROOT, 'experiments', 'lumen-bench'))
    import lumen_forecast as L
    F = json.load(open(os.path.join(ROOT, 'experiments/lumen-bench/prospective/w1-sle/frozen_2026-11-26.json')))
    cat = pd.read_csv(os.path.join(ROOT, 'experiments/lumen-bench/prospective/w1-sle/candidates_w1sle.csv'))
    cph = dict(zip(cat.nct, cat.phase))
    out = {}
    for nct, f in F.items():
        if nct == '_meta': continue
        R = L.load(os.path.join(ROOT, 'experiments/lumen-bench', f['table']))
        P = R[R.week <= 104].sort_values('week').groupby(['trial', 'arm']).tail(1).copy()
        P['ph3'] = [ph3(n) for n in P.nct]
        phs = str(cph.get(nct, '')); x0 = 1 if 'PHASE3' in phs else 0
        y, se = P.lor.values.astype(float), P.se.values.astype(float); tau2 = L.dl(y, se)
        tph = P.drop_duplicates('trial').ph3; fb = bool((tph == 1).sum() < 3 or (tph == 0).sum() < 3)
        if fb: m, v, _, _ = wls(y, se, None, None, tau2)
        else: m, v, _, _ = wls(y, se, P.ph3.values.astype(float), x0, tau2)
        s = math.sqrt(v)
        out[nct] = dict(trial=f['trial'], table=f['table'], catalogue_phase=phs, ph3=x0, fallback=fb,
                        M1_frozen=[f['M1_pooled']['prior_logOR'], f['M1_pooled']['prior_sd']],
                        P1=dict(mean=round(m, 3), sd=round(s, 3), logOR_95=[round(m - 1.96 * s, 3), round(m + 1.96 * s, 3)]))
    return out

def main():
    tf, lo, tf23 = build('tf'), build('loto'), build('tf', recode_23=True)
    res = dict(prereg=PREREG, prereg_commit='05328df', md5_verify=ok)
    res['primary'] = summ(tf, 'time-forward, all subset, LN+SLE, IV-MSE')
    res['verdict'] = verdict(res['primary'])
    res['S_a_loto'] = summ(lo, 'LOTO'); res['S_a_loto']['verdict_by_rule_descriptive'] = verdict(res['S_a_loto'])
    res['S_b_LN'] = summ(tf[tf.stratum == 'LN'], 'time-forward LN only')
    res['S_b_SLE'] = summ(tf[tf.stratum == 'SLE'], 'time-forward SLE only')
    res['S_c_ph23_as_ph2'] = summ(tf23, 'time-forward, PHASE2/PHASE3 recoded as phase 2')
    res['S_f_full_fit'] = full_fit()
    nf = tf[~tf.fallback]
    res['S_g_by_heldout_phase'] = dict(ph2=summ(nf[nf.ph3 == 0], 'held-out phase-2 arms, no fallback'),
                                       ph3=summ(nf[nf.ph3 == 1], 'held-out phase-3 arms, no fallback'))
    res['fallback_trials_tf'] = sorted(tf[tf.fallback].trial.unique().tolist())
    res['b_train_range_tf'] = [round(float(nf.b_train.min()), 3), round(float(nf.b_train.max()), 3)] if len(nf) else None
    pd.concat([tf, lo, tf23.assign(mode='tf_ph23_as_ph2')]).to_csv(os.path.join(OUT, 'o3_predictions.csv'), index=False)
    res = provenance.stamp(res, prereg=PREREG)
    json.dump(res, open(os.path.join(OUT, 'o3_result.json'), 'w'), indent=1, default=str)
    fz = provenance.stamp(dict(note='Descriptive prospective record per PREREG-O3 s11; no outcomes read; PREREG-W1-SLE unchanged',
                                forecasts=freeze_w1()), prereg=PREREG)
    json.dump(fz, open(os.path.join(ROOT, 'experiments/lumen-bench/prospective/w1-sle/frozen_O3_2026-11-27.json'), 'w'), indent=1)
    print(json.dumps(res, indent=1, default=str))

if __name__ == '__main__':
    main()
