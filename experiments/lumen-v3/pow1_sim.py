#!/usr/bin/env python3
"""POW1 (round 167): W1 detectability audit. SIMULATION ONLY - reads no outcome data (none exists: W1 0/21 scorable).

Design (fixed before the first run, committed before results):
  Inputs: the write-once REC1 freeze (21 W1 + W1-SLE trials: frozen P1, P_a, P_b, P_c, B, D, se_ant, arms) and the
  frozen O1/O2/O3 SLE records (8/5/8 trials). Nothing else.
  Simulated outcome: per active arm e ~ Bernoulli(p_true) (arms independent -> power is an UPPER bound, because real
  arms share a control). Larger prospective sets (n 30/50/100, 'W3') resample the 21 W1 trials with replacement.
  Sims: 2000 per cell; paired trial bootstrap 1000 resamples (PREREGs use 2000; 1000 changes CI ends by < 1%).
  Questions, each answered for n = 5, 10, 21 (and W3 sizes):
   Q1 CAL1-type gross calibration: power that the bootstrap 95% CI of g = mean(P1 - e) excludes 0 when the true
      rate is P1 - delta; minimum detectable delta (MDD) at 80% power. Also power of CAL1's own FAIL rule
      (|g| >= 0.10 and CI excludes 0 on that side), and the 0.80-vs-0.69 case.
   Q2 REC1 s9: P(SWITCH-eligible) for each competitor when it is the truth; P(any false SWITCH) and P(CLOSE)
      when B is the truth; P(CLOSE) when P1 is the truth.
   Q3 O1/O3/O2 and the CLI early rule (log-OR / logit scale, normal sampling with REC1 se_ant):
      O1: true predictive sd = k x M1 sd; power of PREREG-O1 s8 (a correctly specified model beats M1 by >= 0.05
          nats/trial with paired CI < 0); MDD k. O3: true phase-2 shift beta on log OR; power of P1 IV-MSE <= 0.90 x
          M1 (O3 s9); MDD beta. Generic paired IV-MSE ratio (CLI early rule, rule 2b 10% margin): error-sd ratio r
          giving 80% power. O2: control-rate logit bias b detectable by NLPD margin 0.05; MDD b.
      Each is also reported against its own PREREG minimum n (O1/O3/CAL1 >= 20, O2 >= 20 + >= 10 LN, CLI >= 5
      both-visit trials) and the number of frozen prospective records that exist.
Output: pow1_result.json via scripts/provenance stamp (lumen.provenance.stamp).
"""
import json, os, sys, math
import numpy as np
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
PROS = os.path.join(ROOT, 'experiments', 'lumen-bench', 'prospective')
NSIM, NBOOT = 2000, 1000
SIZES = [5, 10, 21, 30, 50, 100]
rng = np.random.default_rng(20261128)


def load():
    d = json.load(open(os.path.join(PROS, 'rec1', 'frozen_REC1_2026-11-28.json')))
    T = []
    for k, v in d.items():
        if k == '_meta':
            continue
        nc = max(1, len(v['control_arms']))
        m = max(1, v['n_arms'] - nc)
        T.append(dict(nct=k, table=v['table'], m=m, se=v['se_ant'], P1=v['P1'], a=v['P_a'], b=v['P_b'],
                      c=v['P_c'], B=v['B'], D=v['D_legacy_latent'], M1=v['M1'], v=v['v']))
    return T


def draw_idx(nT, n):
    if n <= nT:
        return np.array([rng.choice(nT, n, replace=False) for _ in range(NSIM)])
    return rng.integers(0, nT, (NSIM, n))


def sim_ebar(p, m):
    """p, m: (NSIM, n) per-trial true prob and active-arm count -> within-trial mean of Bernoulli arms."""
    out = np.zeros_like(p)
    for mm in np.unique(m):
        sel = m == mm
        out[sel] = rng.binomial(mm, np.clip(p[sel], 0, 1)) / mm
    return out


def boot_ci(x):
    """x: (NSIM, n) per-trial values -> percentile 95% CI of the mean, per sim."""
    nsim, n = x.shape
    idx = rng.integers(0, n, (NBOOT, n))
    lo = np.empty(nsim); hi = np.empty(nsim)
    for s in range(nsim):
        bm = x[s][idx].mean(1)
        lo[s], hi[s] = np.percentile(bm, [2.5, 97.5])
    return lo, hi


def mdd(powers, grid):
    for g, p in zip(grid, powers):
        if p >= 0.80:
            return g
    return None


def q1(T):
    P1 = np.array([t['P1'] for t in T]); D = np.array([t['D'] for t in T]); m = np.array([t['m'] for t in T])
    res = {}
    deltas = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
    for n in SIZES:
        idx = draw_idx(len(T), n)
        row = {}
        for dlt in deltas:
            ptrue = P1[idx] - dlt
            e = sim_ebar(ptrue, m[idx])
            g = P1[idx] - e
            lo, hi = boot_ci(g)
            gm = g.mean(1)
            row[f'{dlt:.2f}'] = dict(power_CI_excludes_0=float(np.mean(lo > 0)),
                                     power_CAL1_FAIL=float(np.mean((gm >= 0.10) & (lo > 0))))
        res[n] = dict(by_delta=row,
                      MDD_delta_80pct_CI=mdd([row[f'{x:.2f}']['power_CI_excludes_0'] for x in deltas], deltas),
                      MDD_delta_80pct_CAL1_FAIL=mdd([row[f'{x:.2f}']['power_CAL1_FAIL'] for x in deltas], deltas))
    # the research doc's rough case: forecast 0.80 vs true 0.69, one arm per trial
    case = {}
    for n in [5, 10, 21, 50]:
        e = rng.binomial(1, 0.69, (NSIM, n)).astype(float)
        g = 0.80 - e
        lo, hi = boot_ci(g)
        case[n] = float(np.mean(lo > 0))
    res['case_0.80_vs_0.69_one_arm'] = case
    return res


def rec1_rule(P, e, B, D):
    """P: dict name->(NSIM,n) frozen probs; e: (NSIM,n) within-trial means (approximation: Brier on trial means of
    arm Brier = use per-trial (P-e)^2 with e the arm-mean; exact for 1 arm). Returns dict of eligibility arrays."""
    brier = {k: ((v - e) ** 2) for k, v in P.items()}
    bB = brier['B']
    gap = {k: np.abs((v - e).mean(1)) for k, v in P.items()}
    out = {}
    for X in ['a', 'b', 'c']:
        dX = brier[X] - bB
        lo, hi = boot_ci(dX)
        mX, mB = brier[X].mean(1), bB.mean(1)
        elig = (mX <= 0.90 * mB) & (hi < 0) & (gap[X] < gap['B']) & (mX < brier['D'].mean(1)) & \
               (mX <= brier['C05'].mean(1))
        fail = (mX >= 1.10 * mB) | (gap[X] >= gap['B'] + 0.10)
        out[X] = (elig, fail)
    return out


def q2(T):
    keys = ['P1', 'a', 'b', 'c', 'B', 'D']
    A = {k: np.array([t[k] for t in T]) for k in keys}
    m = np.array([t['m'] for t in T])
    res = {}
    for n in [5, 10, 21, 50, 100]:
        idx = draw_idx(len(T), n)
        P = {k: A[k][idx] for k in keys}
        P['C05'] = np.full_like(P['B'], 0.5)
        row = {}
        for truth in ['a', 'b', 'c', 'B', 'P1']:
            e = sim_ebar(P[truth], m[idx])
            r = rec1_rule(P, e, P['B'], P['D'])
            anyel = r['a'][0] | r['b'][0] | r['c'][0]
            allfail = r['a'][1] & r['b'][1] & r['c'][1]
            row[f'truth_{truth}'] = dict(
                P_switch_eligible={X: float(r[X][0].mean()) for X in 'abc'},
                P_any_switch=float(anyel.mean()), P_close_all_fail=float(allfail.mean()),
                P_keep_open=float((~anyel & ~allfail).mean()))
        res[n] = row
    # how far apart are the frozen forecasts? expected Brier gap if truth = X
    sep = {}
    for X in 'abc':
        p = A[X]
        sep[X] = dict(mean_abs_diff_vs_B=float(np.mean(np.abs(p - A['B']))),
                      expected_brier_ratio_X_over_B_if_X_true=float(np.mean(p * (1 - p)) /
                                                                  np.mean(p * (1 - p) + (p - A['B']) ** 2)))
    res['separation'] = sep
    return res


def nlpd(y, mu, var):
    return 0.5 * np.log(2 * np.pi * var) + 0.5 * (y - mu) ** 2 / var


def q3(T):
    se = np.array([t['se'] for t in T]); M1 = np.array([t['M1'] for t in T]); v = np.array([t['v'] for t in T])
    res = {}
    # O1: true tau inflation k
    ks = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
    o1 = {}
    for n in [8, 21, 50, 100]:
        idx = draw_idx(len(T), n)
        s, mu, vv = se[idx], M1[idx], v[idx]
        row = {}
        for k in ks:
            theta = mu + np.sqrt(k ** 2 * vv) * rng.standard_normal(mu.shape)
            y = theta + s * rng.standard_normal(mu.shape)
            d = nlpd(y, mu, k ** 2 * vv + s ** 2) - nlpd(y, mu, vv + s ** 2)  # correct model minus M1
            lo, hi = boot_ci(d)
            row[f'{k}'] = float(np.mean((d.mean(1) <= -0.05) & (hi < 0)))
        o1[n] = dict(power_by_sd_inflation=row, MDD_k=mdd([row[f'{k}'] for k in ks], ks))
    res['O1_interval_model'] = o1
    # O3: phase-2 optimism beta (true log OR of phase-2 trials = M1 - beta); P1 knows beta; IV-MSE ratio <= 0.90
    betas = [0.0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0]
    o3 = {}
    for n in [8, 21, 50, 100]:
        idx = draw_idx(len(T), n)
        s, mu, vv = se[idx], M1[idx], v[idx]
        ph2 = rng.random(mu.shape) < 0.5
        row = {}
        for bta in betas:
            theta = mu - bta * ph2 + np.sqrt(vv) * rng.standard_normal(mu.shape)
            y = theta + s * rng.standard_normal(mu.shape)
            w = 1 / (s ** 2)
            mse1 = (w * (y - mu) ** 2).sum(1) / w.sum(1)
            mseP = (w * (y - (mu - bta * ph2)) ** 2).sum(1) / w.sum(1)
            row[f'{bta}'] = float(np.mean(mseP <= 0.90 * mse1))
        o3[n] = dict(power_by_beta=row, MDD_beta=mdd([row[f'{b}'] for b in betas], betas))
    res['O3_phase_term'] = o3
    # generic paired IV-MSE: forecaster A error sd = s0, B error sd = r*s0 (same true effects), B loses by >= 10%
    rs = [1.1, 1.25, 1.5, 2.0, 3.0]
    gen = {}
    for n in [5, 10, 21, 50]:
        idx = draw_idx(len(T), n)
        s, vv = se[idx], v[idx]
        row = {}
        for r in rs:
            tot0 = vv + s ** 2
            eA = np.sqrt(tot0) * rng.standard_normal(s.shape)
            eB = eA + np.sqrt((r ** 2 - 1) * tot0) * rng.standard_normal(s.shape)
            w = 1 / s ** 2
            mA = (w * eA ** 2).sum(1) / w.sum(1); mB = (w * eB ** 2).sum(1) / w.sum(1)
            row[f'{r}'] = float(np.mean(mB >= 1.10 * mA))
        gen[n] = dict(P_worse_by_10pct_detected=row, MDD_sd_ratio=mdd([row[f'{r}'] for r in rs], rs),
                      false_alarm_r1=None)
        eA = np.sqrt(vv + s ** 2) * rng.standard_normal(s.shape); eB = np.sqrt(vv + s ** 2) * rng.standard_normal(s.shape)
        w = 1 / s ** 2
        gen[n]['false_alarm_r1'] = float(np.mean((w * eB ** 2).sum(1) / w.sum(1) >= 1.10 * (w * eA ** 2).sum(1) / w.sum(1)))
    res['generic_IV_MSE_10pct_rule'] = gen
    # O2: control logit bias b, control arm n = n_a from se (approx: se_ctrl ~ se/sqrt(2)); predictive sd 0.6 logit
    bs = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5]
    o2 = {}
    for n in [5, 20, 50]:
        idx = draw_idx(len(T), n)
        sc = se[idx] / math.sqrt(2)
        psd = 0.6  # order of the frozen O2 95% interval half-width (~1.2 logit) / 1.96
        row = {}
        for b in bs:
            truth = b + psd * rng.standard_normal(sc.shape)
            y = truth + sc * rng.standard_normal(sc.shape)
            d = nlpd(y, b, psd ** 2 + sc ** 2) - nlpd(y, 0.0, psd ** 2 + sc ** 2)
            lo, hi = boot_ci(d)
            row[f'{b}'] = float(np.mean((d.mean(1) <= -0.05) & (hi < 0)))
        o2[n] = dict(power_by_logit_bias=row, MDD_bias=mdd([row[f'{b}'] for b in bs], bs))
    res['O2_control_rate'] = o2
    return res


def frozen_counts():
    c = {}
    for nm, f in [('O1', 'frozen_O1_2026-11-27.json'), ('O2', 'frozen_O2_2026-11-27.json'),
                  ('O3', 'frozen_O3_2026-11-27.json')]:
        c[nm] = len(json.load(open(os.path.join(PROS, 'w1-sle', f)))['forecasts'])
    c['REC1_and_CAL1'] = 21
    c['REC1_LN'] = 13
    return c


def main():
    T = load()
    out = dict(item='POW1', round=167, note='simulation only; no outcome data read', nsim=NSIM, nboot=NBOOT,
               seed=20261128, n_trials_frozen=len(T), frozen_counts=frozen_counts(),
               prereg_min_n={'CAL1': 20, 'O1': 20, 'O2': '20 combined + 10 LN', 'O3': '20 incl >=6 per phase',
                             'REC1': '5 look / 10 switch', 'CLI_early': '5 trials with both visits'},
               Q1_calibration=q1(T), Q2_REC1=q2(T), Q3=q3(T))
    sys.path.insert(0, ROOT)
    try:
        from lumen.provenance import stamp
        out = stamp(out, prereg=None)
    except Exception as ex:  # stamp signature differs -> record git head ourselves
        import subprocess
        out['git_head'] = subprocess.run(['git', '-C', ROOT, 'rev-parse', '--short', 'HEAD'], capture_output=True,
                                         text=True).stdout.strip()
        out['stamp_note'] = f'lumen.provenance.stamp not used: {ex}'
    json.dump(out, open(os.path.join(HERE, 'pow1_result.json'), 'w'), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str)[:200])


if __name__ == '__main__':
    main()
