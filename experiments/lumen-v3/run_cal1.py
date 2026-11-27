"""CAL1 (round 161, RT-20261127-1 run 2): implements experiments/lumen-v3/PREREG-CAL1.md (commit 6c2dd0a) exactly.

Is the printed P_benefit calibrated? Display decision only; no model change (rule 8).
Run from the repo root:  python3 -I experiments/lumen-v3/run_cal1.py --selftest   (simulated units only)
                         python3 -I experiments/lumen-v3/run_cal1.py              (real run)
Imports run_rt1.py unchanged (LN, SLE, DATE, CFG, late_arms, score_stratum). Research and simulation only.
"""
import sys, os, json, math, io
import numpy as np, pandas as pd
from scipy.stats import norm
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'experiments', 'lumen-bench', 'timeforward'))
from lumen import provenance
OUT = os.path.dirname(os.path.abspath(__file__))
PREREG = 'experiments/lumen-v3/PREREG-CAL1.md'
SEED = 20261127; NBOOT = 2000; MARGIN = 0.10; MIN_TRIALS = 20


# ---------------------------------------------------------------- scoring core (shared by selftest and real run)
def add_preds(U):
    """U: one row per held-out arm with M1, v_M1, obs, se_obs, H, H2. Adds events and predictions (PREREG §5, §10)."""
    U = U.copy()
    sd = np.sqrt(U.v_M1 + U.se_obs ** 2)
    U['E1'] = (U.obs > 0).astype(float)
    U['P1'] = norm.cdf(U.M1 / sd)
    U['Dlat'] = norm.cdf(U.M1 / np.sqrt(U.v_M1))           # the printed P_benefit
    U['C05'] = 0.5
    U['E2'] = (U.obs / U.se_obs > 1.96).astype(float)
    U['P2'] = norm.cdf((U.M1 - 1.96 * U.se_obs) / sd)
    U['E2b'] = ((U.obs / U.se_obs).abs() > 1.96).astype(float)
    U['P2b'] = norm.cdf((U.M1 - 1.96 * U.se_obs) / sd) + norm.cdf((-U.M1 - 1.96 * U.se_obs) / sd)
    return U

def trial_table(U, preds, event):
    """Within-trial means of each prediction, of the event and of each squared error (unit = trial, PREREG §6)."""
    d = U[['trial', event] + preds].copy()
    for p in preds: d['b_' + p] = (d[p] - d[event]) ** 2
    return d.groupby('trial').mean()

def summarise(U, model, base, event, rng_seed=SEED, units='trial'):
    preds = [model] + base
    if units == 'arm':
        T = U[[event] + preds].copy()
        for p in preds: T['b_' + p] = (T[p] - T[event]) ** 2
    else:
        T = trial_table(U, preds, event)
    n = len(T)
    gap = float((T[model] - T[event]).mean())
    brier = {p: float(T['b_' + p].mean()) for p in preds}
    rng = np.random.default_rng(rng_seed)
    gb, db = [], []
    A = T.to_dict('list'); arr = {k: np.asarray(v) for k, v in A.items()}
    for _ in range(NBOOT):
        i = rng.integers(0, n, n)
        gb.append(float((arr[model][i] - arr[event][i]).mean()))
        db.append(float(arr['b_' + model][i].mean() - arr['b_' + base[0]][i].mean()))
    ci = lambda x: [round(float(np.percentile(x, 2.5)), 4), round(float(np.percentile(x, 97.5)), 4)]
    out = dict(n_units=int(n), units=units, mean_pred=round(float(T[model].mean()), 4),
               observed_fraction=round(float(T[event].mean()), 4), gap=round(gap, 4), gap_ci95=ci(gb),
               brier={p: round(v, 4) for p, v in brier.items()},
               brier_diff_vs_strongest=round(brier[model] - brier[base[0]], 4), brier_diff_ci95=ci(db),
               mean_pred_baselines={p: round(float(T[p].mean()), 4) for p in base})
    out['verdict'] = verdict(out, model, base)
    return out

def verdict(s, model, base):
    """PREREG §7-9. base[0] = strongest trivial baseline (H or H2), base[1] = C05."""
    if s['n_units'] < MIN_TRIALS: return 'UNDERPOWERED'
    b = s['brier']; g = s['gap']; lo, hi = s['gap_ci95']
    fail = abs(g) >= MARGIN or b[model] >= 1.10 * b[base[0]] or b[model] >= 1.10 * b[base[1]]
    if fail: return 'FAIL'
    if abs(g) < MARGIN and lo > -MARGIN and hi < MARGIN and b[model] <= b[base[0]] and b[model] <= b[base[1]]:
        return 'PASS'
    return 'INCONCLUSIVE'

def reliability(U, col, event='E1'):
    """D2: trials split into tertiles of mean-arm prediction."""
    T = U.groupby('trial')[[col, event]].mean().sort_values(col)
    T['tert'] = pd.qcut(T[col].rank(method='first'), 3, labels=['low', 'mid', 'high'])
    return {str(k): dict(n_trials=int(len(g)), mean_pred=round(float(g[col].mean()), 3),
                         observed=round(float(g[event].mean()), 3)) for k, g in T.groupby('tert', observed=True)}

def coverage95(U):
    half = 1.96 * np.sqrt(U.v_M1 + U.se_obs ** 2)
    return dict(n_arms=int(len(U)), coverage=round(float(((U.obs - U.M1).abs() <= half).mean()), 3))


# ---------------------------------------------------------------- self-test on simulated units
def selftest():
    rng = np.random.default_rng(1)
    rows = []
    for t in range(200):
        mu, tau = 0.3, 0.3; vm = tau ** 2 + 0.01
        theta = rng.normal(mu, tau)
        for a in range(2):
            se = rng.uniform(0.2, 0.5)
            rows.append(dict(trial=f't{t}', M1=mu, v_M1=vm, obs=rng.normal(theta, se), se_obs=se))
    U = pd.DataFrame(rows)
    H = norm.cdf(0.3 / np.sqrt(0.1 + 0.35 ** 2)); U['H'] = H; U['H2'] = 0.2
    U = add_preds(U)
    s = summarise(U, 'P1', ['H', 'C05'], 'E1')
    print('selftest calibrated sim:', json.dumps(s))
    assert abs(s['gap']) < 0.05, s                 # calibrated simulator -> gap ~ 0
    assert s['brier']['P1'] <= s['brier']['C05'] + 1e-9
    V = U.copy(); V['P1'] = 0.95                   # deliberately over-confident -> FAIL
    s2 = summarise(V, 'P1', ['H', 'C05'], 'E1'); assert s2['verdict'] == 'FAIL', s2
    s3 = summarise(U.iloc[:30], 'P1', ['H', 'C05'], 'E1'); assert s3['verdict'] == 'UNDERPOWERED', s3
    print('selftest OK: calibrated gap', s['gap'], '| over-confident ->', s2['verdict'], '| n 15 ->', s3['verdict'])


# ---------------------------------------------------------------- real run
def baselines_for(R, st, mode, R1):
    """H and H2 per scored trial, from the same training trials as run_rt1.score_stratum (PREREG §5)."""
    trials = R[['trial', 'nct']].drop_duplicates()
    trials = trials.assign(date=trials.nct.map(R1.DATE))
    out = {}
    for _, tr in trials.iterrows():
        tt = trials[trials.date < tr.date].trial if mode == 'tf' else trials[trials.trial != tr.trial].trial
        train = R[R.trial.isin(tt)]
        if train.trial.nunique() < 5: continue
        P = R1.late_arms(train, st)
        out[tr.trial] = dict(H=float((P.lor > 0).mean()), H2=float((P.lor / P.se > 1.96).mean()),
                             n_train_arms=int(len(P)))
    return out

def units(R, st, mode, R1):
    D = pd.DataFrame(R1.score_stratum(R, st, mode))
    D = D[D.subset == 'all'].copy()
    B = baselines_for(R, st, mode, R1)
    assert set(D.trial) == set(B), 'baseline/score trial sets differ'
    for k in ['H', 'H2', 'n_train_arms']: D[k] = D.trial.map(lambda t: B[t][k])
    return D

def main():
    ok = provenance.verify_all(); assert all(ok.values()), ok
    import run_rt1 as R1   # unchanged
    tf = pd.concat([units(R1.LN, 'LN', 'tf', R1), units(R1.SLE, 'SLE', 'tf', R1)], ignore_index=True)
    # assertion against the committed RT-1 predictions (PREREG §4)
    ref = pd.read_csv(os.path.join(ROOT, 'experiments/lumen-bench/timeforward/rt1_predictions.csv'))
    ref = ref[(ref['mode'] == 'tf') & (ref.subset == 'all')]
    m = tf.merge(ref, on=['stratum', 'trial', 'arm'], suffixes=('', '_ref'))
    assert len(m) == len(tf) == len(ref), (len(m), len(tf), len(ref))
    for c in ['M1', 'v_M1', 'obs', 'se_obs']:
        assert np.allclose(m[c], m[c + '_ref'], atol=1e-9, rtol=0), c
    tf = add_preds(tf)
    res = dict(prereg=PREREG, label='12b-timed (committed data; earlier looks not excludable)',
               n_trials=int(tf.trial.nunique()), n_arms=int(len(tf)),
               n_trials_by_stratum={s: int(g.trial.nunique()) for s, g in tf.groupby('stratum')},
               assertion='M1, v_M1, obs, se_obs equal rt1_predictions.csv (tf, all) to 1e-9: passed')
    mol = tf.groupby('trial').molecule.first().str.lower()
    res['molecule_note'] = dict(belimumab=int(mol.str.contains('belimumab').sum()),
                                anifrolumab=int(mol.str.contains('anifrolumab').sum()),
                                n_distinct_molecules=int(mol.nunique()))
    res['A_primary'] = summarise(tf, 'P1', ['H', 'C05'], 'E1')
    res['B_secondary_sig'] = summarise(tf, 'P2', ['H2', 'C05'], 'E2')
    res['B_prime_either_direction'] = dict(n_trials=int(tf.trial.nunique()),
                                           mean_pred=round(float(tf.groupby('trial').P2b.mean().mean()), 4),
                                           observed=round(float(tf.groupby('trial').E2b.mean().mean()), 4))
    res['D1_printed_P_benefit'] = summarise(tf, 'Dlat', ['H', 'C05'], 'E1')
    res['D2_reliability'] = dict(P1=reliability(tf, 'P1'), H=reliability(tf, 'H'), Dlat=reliability(tf, 'Dlat'))
    res['D3_by_stratum'] = {s: summarise(g, 'P1', ['H', 'C05'], 'E1') for s, g in tf.groupby('stratum')}
    lo = pd.concat([units(R1.LN, 'LN', 'loto', R1), units(R1.SLE, 'SLE', 'loto', R1)], ignore_index=True)
    lo = add_preds(lo)
    res['D4_loto'] = summarise(lo, 'P1', ['H', 'C05'], 'E1')
    res['D5_arm_level'] = summarise(tf, 'P1', ['H', 'C05'], 'E1', units='arm')
    res['D6_coverage95_M1'] = coverage95(tf)
    BI = R1.load(os.path.join(ROOT, 'experiments/lumen-bench/sle/bicla_table.csv')) if hasattr(R1, 'load') else None
    if BI is not None and BI.nct.map(R1.DATE).notna().all():
        bt = add_preds(units(BI, 'SLE', 'tf', R1))
        res['D7_bicla'] = summarise(bt, 'P1', ['H', 'C05'], 'E1')
        res['D7_bicla']['printed_P_benefit'] = summarise(bt, 'Dlat', ['H', 'C05'], 'E1')
    else:
        res['D7_bicla'] = 'not scored (BICLA trials lack resultsFirstPostDate in the committed snapshot)'
    res['decision_a'] = {'PASS': 'KEEP', 'FAIL': 'DROP', 'INCONCLUSIVE': 'KEEP with "calibration not established"',
                         'UNDERPOWERED': 'KEEP, "calibration not scored"'}[res['A_primary']['verdict']]
    res['decision_b'] = 'add P_sig_benefit (Emerging)' if res['B_secondary_sig']['verdict'] == 'PASS' else 'nothing added'
    res = provenance.stamp(res, prereg=PREREG)
    cols = ['stratum', 'trial', 'nct', 'date', 'arm', 'molecule', 'n_train_trials', 'n_train_arms', 'M1', 'v_M1',
            'obs', 'se_obs', 'E1', 'P1', 'Dlat', 'H', 'E2', 'P2', 'H2']
    buf = io.StringIO(); tf[cols].to_csv(buf, index=False)
    provenance.write_once_text(os.path.join(OUT, 'cal1_units.csv'), buf.getvalue())
    provenance.write_once(os.path.join(OUT, 'cal1_result.json'), res, default=str)
    print(json.dumps(res, indent=1, default=str))

if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else main()
