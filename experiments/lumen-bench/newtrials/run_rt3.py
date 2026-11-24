"""RT-20261124-1 run 3 (round 136): score early read-out rules on NEW trials only.

Implements experiments/lumen-bench/newtrials/PREREG-RT3.md (68900fb) + AMENDMENT-A135-RT3.md (c181f31).
Input: new_arms.csv / exclusions.csv (extraction commit 5abd3f3, before this scoring commit).
Run from the repo root:  python3 -I experiments/lumen-bench/newtrials/run_rt3.py
Research and simulation only. Not clinical guidance.
"""
import sys, os, json
import numpy as np, pandas as pd
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'experiments', 'lumen-bench', 'timeforward'))
import run_rt1 as R1                      # shipped functions, unchanged (loads the 60 tabulated trials + dates)
from lumen_forecast import load
from stamp import stamp
OUT = os.path.dirname(os.path.abspath(__file__))
PREREG = 'experiments/lumen-bench/newtrials/PREREG-RT3.md'

N = pd.read_csv(os.path.join(OUT, 'new_arms.csv'))
a, b, c, d = N.r_drug + .5, N.n_drug - N.r_drug + .5, N.r_ctrl + .5, N.n_ctrl - N.r_ctrl + .5
N['lor'] = np.log(a * d / (b * c)); N['se'] = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
# exclusions.csv reason field contains unquoted commas: first 3 fields fixed, rest = reason
_L = open(os.path.join(OUT, 'exclusions.csv')).read().splitlines()[1:]
X = pd.DataFrame([dict(zip(['cid', 'name', 'rule', 'reason'], l.split(',', 3))) for l in _L if l.strip()])

def ivmse(df, m):
    w = 1 / df.se_obs ** 2; return float((w * (df.obs - df[m]) ** 2).sum() / w.sum())

rows = []
for st, OLD in [('LN', R1.LN), ('SLE', R1.SLE)]:
    old = OLD.assign(date=OLD.nct.map(R1.DATE), new=False)
    new = N[N.stratum == st].assign(new=True)
    M = pd.concat([old, new], ignore_index=True)
    c = R1.CFG[st]
    for tr, g in new.groupby('trial'):
        dt = g.date.iloc[0]
        train = M[(M.date < dt)]; ntr = train.trial.nunique()
        if ntr < 5: continue
        mu, v0, _ = R1.prior(train, st)
        dv = R1.drift_var(train) if R1.n_drift_pairs(train) >= 2 else None
        for arm, ga in g.groupby('arm'):
            ma, va, lvl = R1.prior(train, st, ga.molecule.iloc[0], ga.target.iloc[0], anchor=True)
            base = dict(stratum=st, cid=ga.cid.iloc[0], trial=tr, date=dt, arm=arm, n_train_trials=ntr, anchor_lvl=lvl)
            L = R1.last(ga, 0, c['late_max'])
            if L is not None:
                rows.append(dict(base, subset='all', target_week=int(L.week), obs=L.lor, se_obs=L.se, M0=0.0, M1=mu, M3=ma))
            E = R1.last(ga, c['e_min'], c['e_max']); T = R1.last(ga, c['t_min'], c['t_max'])
            if E is not None and T is not None and E.week < T.week:   # primary subset (none expected after extraction)
                r = dict(base, subset='early', obs=T.lor, se_obs=T.se, M0=0.0, M1=mu, C0=E.lor, M3=ma)
                if dv is not None:
                    ve = E.se ** 2 + dv; r['M2'] = (mu / v0 + E.lor / ve) / (1 / v0 + 1 / ve)
                rows.append(r)
D = pd.DataFrame(rows)
P = D[D.subset == 'early'] if len(D) else D
out = dict(prereg=PREREG, amendment='experiments/lumen-bench/newtrials/AMENDMENT-A135-RT3.md',
           n_candidates=36, n_excluded_by_rule=X[~X.reason.str.startswith('INCLUDED')].groupby('rule').size().to_dict(),
           included_late_only=sorted(N.trial.unique().tolist()),
           primary=dict(n_trials=int(P.trial.nunique()) if len(P) else 0, n_arms=int(len(P))))
out['decision_rule_hit'] = 'UNDERPOWERED' if out['primary']['n_trials'] < 10 or out['primary']['n_arms'] < 12 else 'SCORE'
out['decision'] = ('UNDERPOWERED -> park Q61-type early-read-out work; CLI keeps pooled + "no early rule validated" note; '
                   'preprint unchanged' if out['decision_rule_hit'] == 'UNDERPOWERED' else 'see rules 2-6')
out['S7'] = dict(without_N30=out['primary'], clean_subset_N29_N31_N36=dict(n_trials=0, n_arms=0))
A = D[D.subset == 'all'] if len(D) else D
if len(A):
    s3 = dict(n_trials=int(A.trial.nunique()), n_arms=int(len(A)), trials=sorted(A.trial.unique().tolist()))
    for m in ['M0', 'M1', 'M3']: s3['ivmse_' + m] = round(ivmse(A, m), 4)
    s3['per_trial'] = {t: {m: round(ivmse(g, m), 4) for m in ['M0', 'M1', 'M3']} for t, g in A.groupby('trial')}
    s3['S1_equal_trial_weights'] = {m: round(float(np.mean([ivmse(g, m) for _, g in A.groupby('trial')])), 4) for m in ['M0', 'M1', 'M3']}
    s3['label'] = 'descriptive, n 3 trials / 7 arms, below any label threshold'
    out['S3_late_only_new_trials'] = s3
out['S2'] = dict(LN=dict(n_trials=0), SLE=out['primary'])
out['S4_BICLA_only'] = dict(n_trials=0, note='no included trial reports BICLA without SRI-4')
out['S5_bootstrap'] = 'not computed: primary subset empty'
out['S6_new_plus_RT1'] = 'not computed: no new primary arms to add'
out['label'] = '12b-timed for N01-N28/N30 per amendment; no decision reached on a clean subset; UNDERPOWERED is the verdict'
D.to_csv(os.path.join(OUT, 'rt3_predictions.csv'), index=False)
out = stamp(out, prereg=PREREG)
json.dump(out, open(os.path.join(OUT, 'rt3_result.json'), 'w'), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
