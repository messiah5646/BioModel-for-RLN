"""RT-20261123-1 (round 127): time-forward backtest of the shipped LN/SLE trial forecaster.

Implements experiments/lumen-bench/timeforward/PREREG-RT1.md (commit c8e07af) exactly.
Run from the repo root:  python3 -I experiments/lumen-bench/timeforward/run_rt1.py
Research and simulation only. Not clinical guidance.
"""
import sys, os, gzip, json
import numpy as np, pandas as pd
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'experiments', 'lumen-bench'))
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'provenance'))
from lumen_forecast import load, iv, dl, drift_var   # shipped functions, unchanged
from stamp import stamp

OUT = os.path.join(ROOT, 'experiments', 'lumen-bench', 'timeforward')

# --- data ---
ct = json.load(gzip.open(os.path.join(ROOT, 'data/real/ctgov/lupus_with_results.json.gz')))
DATE = {}
for s in ct:
    p = s['protocolSection']
    DATE[p['identificationModule']['nctId']] = p['statusModule'].get('resultsFirstPostDateStruct', {}).get('date')
LN = load(os.path.join(ROOT, 'experiments/lumen-bench/resp_table.csv')); LN['stratum'] = 'LN'
SLE = load(os.path.join(ROOT, 'experiments/lumen-bench/sle/sri4_table.csv')); SLE['stratum'] = 'SLE'
CFG = {  # (late max week, early min, early max, early-subset target min, max)
    'LN': dict(late_max=104, e_min=0, e_max=26, t_min=48, t_max=104),
    'SLE': dict(late_max=52, e_min=12, e_max=26, t_min=48, t_max=52)}

def last(g, lo, hi):
    g = g[(g.week >= lo) & (g.week <= hi)]
    return None if not len(g) else g.sort_values('week').iloc[-1]

def late_arms(R, st):
    c = CFG[st]
    return R[R.week <= c['late_max']].sort_values('week').groupby(['trial', 'arm']).tail(1)

def prior(train, st, molecule=None, target=None, anchor=False):
    P = late_arms(train, st)
    mu, semu = iv(P.lor, P.se); t2 = dl(P.lor.values, P.se.values)
    if not anchor:
        return mu, t2 + semu ** 2, 'pooled'
    anc = P[P.molecule == molecule]; lvl = 'molecule'
    if not len(anc): anc = P[P.target == target]; lvl = 'target'
    if len(anc):
        am, ase = iv(anc.lor, anc.se); w = t2 / (t2 + ase ** 2)
        return mu + w * (am - mu), t2 * (1 - w) + t2, lvl
    return mu, t2 + semu ** 2, 'none'

def n_drift_pairs(R):
    n = 0
    for _, g in R.groupby(['trial', 'arm']):
        if len(g[g.week <= 26]) and len(g[g.week >= 48]): n += 1
    return n

def score_stratum(R, st, mode):
    """mode 'tf' = time-forward (strictly earlier trials), 'loto' = leave-one-trial-out."""
    c = CFG[st]; rows = []
    trials = R[['trial', 'nct']].drop_duplicates()
    trials = trials.assign(date=trials.nct.map(DATE))
    assert trials.date.notna().all(), 'missing resultsFirstPostDate'
    for _, tr in trials.iterrows():
        if mode == 'tf':
            tt = trials[trials.date < tr.date].trial
        else:
            tt = trials[trials.trial != tr.trial].trial
        train = R[R.trial.isin(tt)]
        ntrain = train.trial.nunique()
        if ntrain < 5: continue   # burn-in
        mu, v0, _ = prior(train, st)
        dv = drift_var(train) if n_drift_pairs(train) >= 2 else None
        for arm, g in R[R.trial == tr.trial].groupby('arm'):
            mol, tgt = g.molecule.iloc[0], g.target.iloc[0]
            ma, va, lvl = prior(train, st, mol, tgt, anchor=True)
            base = dict(stratum=st, trial=tr.trial, nct=tr.nct, date=tr.date, arm=arm, molecule=mol,
                        n_train_trials=ntrain, mu=mu, v0=v0, anchor_lvl=lvl)
            L = last(g, 0, c['late_max'])
            if L is not None:   # all-arm (S3/S4)
                rows.append(dict(base, subset='all', target_week=int(L.week), obs=L.lor, se_obs=L.se,
                                 M0=0.0, M1=mu, M3=ma, v_M1=v0, v_M3=va))
            E = last(g, c['e_min'], c['e_max']); T = last(g, c['t_min'], c['t_max'])
            if E is not None and T is not None and E.week < T.week:
                r = dict(base, subset='early', target_week=int(T.week), early_week=int(E.week), obs=T.lor,
                         se_obs=T.se, M0=0.0, M1=mu, C0=E.lor, M3=ma, v_M1=v0, v_M3=va,
                         drift_var=dv, m2_scored=dv is not None)
                if dv is not None:
                    ve = E.se ** 2 + dv
                    r['M2'] = (mu / v0 + E.lor / ve) / (1 / v0 + 1 / ve); r['v_M2'] = 1 / (1 / v0 + 1 / ve)
                rows.append(r)
    return rows

def ivmse(d, m):
    d = d[d[m].notna()]; w = 1 / d.se_obs ** 2
    return float((w * (d.obs - d[m]) ** 2).sum() / w.sum())

def summarise_early(d, label):
    d = d[d.m2_scored == True]
    out = dict(label=label, n_arms=int(len(d)), n_trials=int(d.trial.nunique()),
               trials=sorted(d.trial.unique().tolist()))
    if not len(d): out['verdict'] = 'UNDERPOWERED (0 arms)'; return out
    for m in ['M0', 'M1', 'M2', 'C0', 'M3']: out['ivmse_' + m] = round(ivmse(d, m), 4)
    m2 = out['ivmse_M2']; best = min(out['ivmse_M0'], out['ivmse_M1'], out['ivmse_C0'])
    if out['n_arms'] < 6 or out['n_trials'] < 4: v = 'UNDERPOWERED'
    elif m2 >= 1.10 * best: v = 'FAIL'
    elif m2 <= 0.90 * out['ivmse_M1'] and m2 < 1.10 * out['ivmse_C0'] and m2 < 1.10 * out['ivmse_M0']: v = 'PASS'
    else: v = 'INCONCLUSIVE'
    out['verdict_by_rule'] = v
    out['M2_over_M1'] = round(m2 / out['ivmse_M1'], 3) if out['ivmse_M1'] else None
    return out

def coverage(d, m, vcol):
    d = d[d[m].notna()]
    half = 1.96 * np.sqrt(d[vcol] + d.se_obs ** 2)
    return dict(n=int(len(d)), coverage=round(float(((d.obs - d[m]).abs() <= half).mean()), 3))

def main():
    res = {}
    for mode in ['tf', 'loto']:
        rows = score_stratum(LN, 'LN', mode) + score_stratum(SLE, 'SLE', mode)
        D = pd.DataFrame(rows); D['mode'] = mode; res[mode] = D
    tf, lo = res['tf'], res['loto']
    E = tf[tf.subset == 'early']; A = tf[tf.subset == 'all']
    novoc = E[E.molecule.str.lower() != 'voclosporin']
    out = dict(prereg='experiments/lumen-bench/timeforward/PREREG-RT1.md')
    out['primary'] = summarise_early(novoc, 'time-forward early subset, LN+SLE, voclosporin excluded')
    out['early_subset_unscored_for_M2'] = E[E.m2_scored != True][['stratum', 'trial', 'arm']].to_dict('records')
    out['S1_with_voclosporin'] = summarise_early(E, 'time-forward early subset, LN+SLE, voclosporin included')
    out['S2_LN_novoc'] = summarise_early(novoc[novoc.stratum == 'LN'], 'LN only, voclosporin excluded')
    out['S2_LN_with_voc'] = summarise_early(E[E.stratum == 'LN'], 'LN only, voclosporin included')
    out['S2_SLE'] = summarise_early(E[E.stratum == 'SLE'], 'SLE only')
    s3 = {}
    for st, d in [('LN', A[A.stratum == 'LN']), ('SLE', A[A.stratum == 'SLE']), ('combined', A)]:
        s3[st] = dict(n_arms=int(len(d)), n_trials=int(d.trial.nunique()),
                      ivmse_M0=round(ivmse(d, 'M0'), 4), ivmse_M1=round(ivmse(d, 'M1'), 4),
                      ivmse_M3=round(ivmse(d, 'M3'), 4))
    out['S3_all_arm'] = s3
    out['S4_coverage'] = dict(M1_all_combined=coverage(A, 'M1', 'v_M1'),
                              M1_early_novoc=coverage(novoc[novoc.m2_scored == True], 'M1', 'v_M1'),
                              M2_early_novoc=coverage(novoc[novoc.m2_scored == True], 'M2', 'v_M2'))
    out['S4_W1_rule'] = ('pooled-only coverage >= 0.90 -> W1 criterion (a) cannot tell models apart; report (b) as primary'
                         if out['S4_coverage']['M1_all_combined']['coverage'] >= 0.90 else 'coverage < 0.90; W1 (a) stays')
    LE = lo[lo.subset == 'early']
    out['S5_loto_primary'] = summarise_early(LE[LE.molecule.str.lower() != 'voclosporin'], 'LOTO, same primary subset')
    out['burn_in_scored_trials'] = dict(LN=int(A[A.stratum == 'LN'].trial.nunique()), SLE=int(A[A.stratum == 'SLE'].trial.nunique()))
    pd.concat([tf, lo]).to_csv(os.path.join(OUT, 'rt1_predictions.csv'), index=False)
    out = stamp(out, prereg='experiments/lumen-bench/timeforward/PREREG-RT1.md')
    json.dump(out, open(os.path.join(OUT, 'rt1_result.json'), 'w'), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))

if __name__ == '__main__':
    main()
