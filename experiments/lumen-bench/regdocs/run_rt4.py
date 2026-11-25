"""RT-4 run 3 (round 141): score early read-out rules on the merged early subset.

Implements experiments/lumen-bench/regdocs/PREREG-RT4.md (67f2d26) + AMENDMENT-R20261125-RT4.md (414a21e)
+ AMENDMENT-A140-RT4.md (03df5c2). Input: rt4_rows.csv / rt4_exclusions.csv (extraction commit 1d2f8b0, before scoring).
Run from the repo root:  python3 -I experiments/lumen-bench/regdocs/run_rt4.py
Research and simulation only. Not clinical guidance.
"""
import sys, os, json
import numpy as np, pandas as pd
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'experiments', 'lumen-bench', 'timeforward'))
import run_rt1 as R1                      # shipped functions, unchanged
from stamp import stamp
OUT = os.path.dirname(os.path.abspath(__file__))
PREREG = 'experiments/lumen-bench/regdocs/PREREG-RT4.md'
LABEL = ('12b-timed (late outcomes already seen; early values new for 8 of 10 Table A trials; LOTUS/MUSE early values '
         'were in committed CT.gov data, unread). Round-137 sub-agents opened these documents to find locations '
         '(transcripts not committed); earlier sight of early values and model training knowledge of published '
         'results cannot be excluded.')

def lor(df):
    a, b, c, d = df.r_drug + .5, df.n_drug - df.r_drug + .5, df.r_ctrl + .5, df.n_ctrl - df.r_ctrl + .5
    df = df.copy(); df['lor'] = np.log(a * d / (b * c)); df['se'] = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d); return df

NEW = lor(pd.read_csv(os.path.join(OUT, 'rt4_rows.csv')))
TEL = lor(pd.read_csv(os.path.join(ROOT, 'experiments/lumen-bench/newtrials/new_arms.csv')))
TEL['cls'] = 'other'
R1.DATE['NCT02885610'] = '2023-12-21'          # telitacicept (PREREG section 4 dates)
RT4_TRIALS = sorted(NEW.trial.unique())
cols = ['trial', 'nct', 'arm', 'molecule', 'target', 'cls', 'week', 'r_drug', 'n_drug', 'r_ctrl', 'n_ctrl', 'lor', 'se']
LN = pd.concat([R1.LN[cols], NEW[NEW.stratum == 'LN'][cols]], ignore_index=True)
SLE = pd.concat([R1.SLE[cols], NEW[NEW.stratum == 'SLE'][cols], TEL[cols]], ignore_index=True)

def ivmse(d, m):
    w = 1 / d.se_obs ** 2; return float((w * (d.obs - d[m]) ** 2).sum() / w.sum())

def scores(d):
    return {m: round(ivmse(d, m), 4) for m in ['M0', 'M1', 'M2', 'C0']}

def rule(A, ALN=None, nln=0):
    """PREREG section 9 rules 2-6 (rule 1 = minimum n handled by caller)."""
    def lnok(m): return True if nln < 3 else ALN[m] < 1.10 * ALN['M0']
    cf = A['C0'] <= .9 * A['M1'] and A['C0'] < 1.1 * min(A['M0'], A['M2']) and lnok('C0')
    up = A['M2'] <= .9 * A['M1'] and A['M2'] < 1.1 * min(A['M0'], A['C0']) and lnok('M2')
    if cf and up: return 'CF-SWITCH' if A['C0'] <= A['M2'] else 'UPDATE-SWITCH'
    if cf: return 'CF-SWITCH'
    if up: return 'UPDATE-SWITCH'
    if A['M1'] >= 1.1 * min(A['M0'], A['C0']): return 'DEFAULT-FAIL'
    if A['M1'] <= min(A['M0'], A['M2'], A['C0']): return 'POOLED-WINS'
    return 'INCONCLUSIVE'

def summ(d, label):
    o = dict(label=label, n_trials=int(d.trial.nunique()), n_arms=int(len(d)), trials=sorted(d.trial.unique().tolist()))
    if not len(d): return o
    o['ivmse'] = scores(d)
    ln = d[d.stratum == 'LN']
    o['rule_selected'] = rule(o['ivmse'], scores(ln) if len(ln) else None, ln.trial.nunique())
    o['n_LN_trials'] = int(ln.trial.nunique())
    return o

def collapse(d, mols=('belimumab', 'anifrolumab')):
    """S8: each listed molecule counts as one trial (weight = its arms' weight sum / its n trials)."""
    out = {}
    for m in ['M0', 'M1', 'M2', 'C0']:
        num = den = 0.0
        for key, g in d.groupby(d.apply(lambda r: r.molecule if r.molecule in mols else r.trial, axis=1)):
            w = 1 / g.se_obs ** 2; k = g.trial.nunique()
            se2 = float((w * (g.obs - g[m]) ** 2).sum() / w.sum()); W = float(w.sum()) / k
            num += W * se2; den += W
        out[m] = round(num / den, 4)
    return out

def boot(d, by, seed, B=2000):
    rng = np.random.default_rng(seed); units = d[by].unique(); r1, r2 = [], []
    for _ in range(B):
        s = rng.choice(units, len(units), replace=True)
        bd = pd.concat([d[d[by] == u] for u in s]); A = scores(bd)
        r1.append(A['C0'] / A['M1']); r2.append(A['M2'] / A['M1'])
    q = lambda x: [round(float(np.percentile(x, 2.5)), 3), round(float(np.percentile(x, 97.5)), 3)]
    return {'C0_over_M1_95': q(r1), 'M2_over_M1_95': q(r2), 'units': int(len(units))}

def main():
    D = pd.DataFrame(R1.score_stratum(LN, 'LN', 'tf') + R1.score_stratum(SLE, 'SLE', 'tf'))
    E = D[(D.subset == 'early') & (D.m2_scored == True)].copy()
    novoc = E[E.molecule.str.lower() != 'voclosporin']
    P = novoc
    new_in_primary = sorted(set(P.trial) & set(RT4_TRIALS))
    out = dict(prereg=PREREG, amendments=['414a21e AMENDMENT-R20261125-RT4', '03df5c2 AMENDMENT-A140-RT4'],
               extraction_commit='1d2f8b0', label=LABEL)
    out['primary'] = summ(P, 'merged early subset, LN+SLE, M2-scored, voclosporin excluded')
    n_ok = P.trial.nunique() >= 12 and len(P) >= 18 and len(new_in_primary) >= 5
    out['primary']['rt4_new_trials_in_primary'] = new_in_primary
    out['primary']['minimum_n_met'] = bool(n_ok)
    out['decision'] = 'UNDERPOWERED' if not n_ok else out['primary']['rule_selected']
    out['molecule_count_label'] = '%d of %d primary trials belimumab, %d anifrolumab' % (
        P[P.molecule == 'belimumab'].trial.nunique(), P.trial.nunique(), P[P.molecule == 'anifrolumab'].trial.nunique())
    out['early_subset_not_M2_scored'] = D[(D.subset == 'early') & (D.m2_scored != True)][['stratum', 'trial', 'arm']].to_dict('records')
    out['S1_rt4_new_only'] = summ(P[P.trial.isin(RT4_TRIALS)], 'RT-4-extracted trials only (descriptive)')
    out['S2_LN'] = summ(P[P.stratum == 'LN'], 'LN only, voclosporin excluded')
    out['S2_SLE'] = summ(P[P.stratum == 'SLE'], 'SLE only')
    out['S2b_with_voclosporin'] = summ(E, 'voclosporin included')
    A = D[D.subset == 'all']
    out['S3_all_arm'] = {st: dict(n_trials=int(d.trial.nunique()), n_arms=int(len(d)),
                                  **{m: round(ivmse(d, m), 4) for m in ['M0', 'M1', 'M3']})
                         for st, d in [('LN', A[A.stratum == 'LN']), ('SLE', A[A.stratum == 'SLE']), ('combined', A)]}
    out['S4_BICLA'] = 'not run: TULIP-2 early BICLA not printed (figure not labelled); no RT-4 BICLA rows'
    out['S5_trial_bootstrap'] = boot(P, 'trial', 138)
    out['S5b_molecule_bootstrap'] = boot(P, 'molecule', 139)
    out['S6_equal_trial_weights'] = {m: round(float(np.mean([ivmse(g, m) for _, g in P.groupby('trial')])), 4)
                                     for m in ['M0', 'M1', 'M2', 'C0']}
    out['S7_no_figure_values'] = summ(P[P.trial != 'LOTUS'], 'without figure-label trials (LOTUS)')
    PL = P[P.stratum == 'LN']   # LN clause: no belimumab/anifrolumab arm is in the LN primary, so LN scores are unchanged
    out['S8_molecule_collapse'] = dict(ivmse=collapse(P), rule_selected=rule(collapse(P), scores(PL), PL.trial.nunique()))
    out['S8b_lomo'] = {m: summ(P[P.molecule != m], 'without ' + m) for m in sorted(P.molecule.unique())}
    out['S8b_lomo'] = {k: dict(n_trials=v['n_trials'], ivmse=v.get('ivmse'), rule_selected=v.get('rule_selected'))
                       for k, v in out['S8b_lomo'].items()}
    out['S9_clean_early'] = dict(primary=summ(P[~P.trial.isin(['LOTUS', 'MUSE'])], 'primary without LOTUS, MUSE'),
                                 S1=summ(P[P.trial.isin(set(RT4_TRIALS) - {'LOTUS', 'MUSE'})], 'S1 without LOTUS, MUSE'))
    out['RT1_reference'] = 'RT-1 primary (7 trials): C0 0.146 < M2 0.259 < M1 0.311 < M0 0.367'
    D.to_csv(os.path.join(OUT, 'rt4_predictions.csv'), index=False)
    out = stamp(out, prereg=PREREG)
    json.dump(out, open(os.path.join(OUT, 'rt4_result.json'), 'w'), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))

if __name__ == '__main__':
    main()
