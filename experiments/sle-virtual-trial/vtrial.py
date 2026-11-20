"""SLE / lupus-nephritis in-silico trial harness.

LESSON 0 -- what this file is. A "virtual trial" here means: take a mechanistic model of the disease, build a
CONTROL arm that matches real trial data, add a DRUG described as changes to the model's parameters, simulate both
arms forward, and compare. It is a hypothesis machine, not an oracle: it tells you what a mechanism of a given
strength WOULD do IF the model is right, and how sure you can be given what real data leave undetermined.

Pipeline (each step is a function below):
  1. load_aurora()        real AURORA 1 placebo-arm means (weeks 2-52), the calibration target
  2. build_ensemble()     many plausible control patients-on-average, all consistent with that real data
  3. apply_drug()         translate a drug spec (JSON) into parameter changes
  4. run_trial()          simulate control vs drug for every ensemble member
  5. power()              could a real trial of n per arm even detect the predicted difference?
  6. report()             numbers + plot + the caveats that travel with them

Research and simulation only. Not clinical or drug-development guidance.
"""
import csv, json, os, signal, sys, argparse
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import norm
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'rln-identifiability'))
from model import make_params, simulate, observe                      # noqa: E402  the project's six-state ODE model

DATA = os.path.join(HERE, '..', '..', 'data', 'real', 'aurora1_per_visit.csv')
FIT_WEEKS = [2, 4, 8, 12, 16, 20, 24, 30, 36, 42, 48, 52]
BASE_UPCR, BASE_EGFR, N_PLACEBO = 3.87, 77.4, 178                     # AURORA 1 placebo baseline (round 35)

# LESSON 1 -- k_tox (how much proteinuria itself damages the kidney) is NOT identifiable from real trial data:
# rounds 75, 76, 79 found a ceiling (~3e-5 /day) and no floor. So we never pick one value; every result is
# computed across this whole range, and if the answer changes across it we say so.
KTOX_GRID = [0.0, 1e-6, 3e-6, 1e-5, 3e-5]

class _Timeout(Exception): pass
def _alarm(*_): raise _Timeout()
signal.signal(signal.SIGALRM, _alarm)


def load_aurora():
    rows = list(csv.DictReader([l for l in open(DATA) if not l.lstrip().startswith('#')]))
    o = np.array([[float(r['upcr_change']), float(r['upcr_sd']), float(r['egfr_change']), float(r['egfr_sd'])]
                  for r in rows if r['arm'] == 'placebo' and int(r['week']) in FIT_WEEKS])
    return o                                                          # columns: dUPCR, SD, dEGFR, SD


def base_params(s_inf, k_dam, k_tox):
    p = make_params(0.0)
    p['p0'] = BASE_UPCR; p['G0'] = BASE_EGFR / (1 - p['c_inf'])
    p['s_inf'], p['k_dam'], p['k_tox'], p['k_prog'] = s_inf, k_dam, k_tox, 0.0
    p['dep_kill'] = 0.0
    return p


def trajectory(p, weeks, step):
    """Simulated change from baseline in UPCR and eGFR at the given weeks.
    LESSON 2 -- the 'step': rounds 77-78 found that both arms of AURORA change abruptly within two weeks of
    starting background therapy (steroid pulse is the prime suspect). The model's slow immune chain cannot make
    that, so a fitted one-off offset is added to every post-baseline visit. It is shared by both arms."""
    tg = np.array([0.0] + [w * 7.0 for w in weeks])
    signal.setitimer(signal.ITIMER_REAL, 15.0)
    try:
        Y = simulate(p, tg, rtol=1e-8, atol=1e-10)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    O = np.array([observe(y, p, t) for t, y in zip(tg, Y)])
    d = O[1:, :2] - O[0, :2]
    return d + np.array(step)[None, :], float(Y[-1, 5])               # also return final damage D


A2_MONTHS = [12, 18, 24, 30, 36]
A2_WEEKS = [round(m * 30.44 / 7.0, 2) for m in A2_MONTHS]


def load_aurora2():
    f = os.path.join(HERE, '..', '..', 'data', 'real', 'aurora2_per_visit.csv')
    rows = list(csv.DictReader([l for l in open(f) if not l.lstrip().startswith('#')]))
    return np.array([[float(r['upcr_change']), float(r['upcr_se']), float(r['egfr_change']), float(r['egfr_se'])]
                     for r in rows if r['arm'] == 'placebo'])


def build_ensemble(n_boot=8, seed=0, verbose=True):
    """LESSON 3 -- an ensemble, not a single fit. For each allowed k_tox and each bootstrap copy of the data
    (real means jittered by their standard errors), refit the free parameters to the AURORA placebo arm.
    The spread of the ensemble is the uncertainty that real data leave in the CONTROL arm.

    LESSON 3b -- calibrate on everything the control arm has: AURORA 1 (weeks 2-52) AND AURORA 2 (months 12-36).
    AURORA 2 is a completer subset with eGFR capped at 90, so it gets free level offsets and contributes SHAPE only
    (as rounds 76 and 79 did). The treated arm is never used here, which keeps the voclosporin 3-year slope a
    held-out test. (History, kept honest: an earlier version fitted year 1 only and appeared to fail the 3-year
    backtest badly. That "failure" was a bug in backtest.py, which evaluated months 36-108 instead of 12-36. With the
    bug fixed, the year-1-only version passes on placebo, z = -0.6, and is borderline on treated, z = -2.07. The
    3-year calibration was kept because it uses more real data, not because the first version was broken.)"""
    o = load_aurora(); o2 = load_aurora2()
    se_u = o[:, 1] / np.sqrt(N_PLACEBO); se_g = o[:, 3] / np.sqrt(N_PLACEBO)
    weeks = FIT_WEEKS + A2_WEEKS; n1 = len(FIT_WEEKS)
    rng = np.random.default_rng(seed); members = []
    for kt in KTOX_GRID:
        for b in range(n_boot):
            j = (lambda m, s: m + rng.normal(0, s)) if b else (lambda m, s: m)
            yu, yg = j(o[:, 0], se_u), j(o[:, 2], se_g); y2u, y2g = j(o2[:, 0], o2[:, 1]), j(o2[:, 2], o2[:, 3])
            def res(x):
                try:
                    d, _ = trajectory(base_params(np.exp(min(x[0], 0)), np.exp(x[1]), kt), weeks, (x[2], x[3]))
                except Exception:
                    return np.full(34, 1e3)
                a, c = d[:n1], d[n1:]
                return np.r_[(a[:, 0] - yu) / se_u, (a[:, 1] - yg) / se_g,
                             (c[:, 0] + x[4] - y2u) / o2[:, 1], (c[:, 1] + x[5] - y2g) / o2[:, 3]]
            best = None
            for s in range(3):
                g = np.random.default_rng(100 * s + b)
                x0 = [np.log(g.uniform(0.3, 0.9)), np.log(10 ** g.uniform(-4.5, -3.3)), g.normal(-1, .3), g.normal(3.5, 1), 0.0, 0.0]
                r = least_squares(res, x0, method='lm', max_nfev=400)
                if best is None or r.cost < best.cost: best = r
            m = dict(k_tox=kt, s_inf=float(np.exp(min(best.x[0], 0))), k_dam=float(np.exp(best.x[1])),
                     step=[float(best.x[2]), float(best.x[3])], chi2=float(2 * best.cost), boot=b)
            members.append(m)
            if verbose: print(f"  member k_tox={kt:.0e} boot={b}: chi2={m['chi2']:.1f} s_inf={m['s_inf']:.2f} k_dam={m['k_dam']:.2e}", flush=True)
    # keep members that fit their own (jittered) data acceptably: 28 dof, 95% point 41.3, x2 slack because the
    # published errors ignore visit-to-visit correlation (round 78 section 5)
    kept = [m for m in members if m['chi2'] <= 2 * 41.3]
    return kept, members


MECHANISMS = {
    # name: (what it means biologically, how it maps onto the model)
    'immunosuppression': ('sustained dampening of autoimmune drive (MMF-, CNI-, JAK-like)', 's_inf *= (1 - strength)'),
    'antiproteinuric': ('lowers glomerular protein leak without touching immunity (SGLT2i-, RAAS-like)', 'leak_red = strength'),
    'haemodynamic': ('reversible fall in filtration from vascular tone (CNI-like side effect)', 'h_drug = strength'),
    'bcell_depletion': ('kills B cells and short-lived plasma cells for a window (anti-CD19/CD20, CAR-T-like)',
                        'dep_kill = strength per day between start_day and end_day'),
    'activity_clearance': ('speeds resolution of renal inflammation directly (complement/anti-inflammatory-like)',
                           'k_clear *= (1 + strength)'),
}
UNSUPPORTED = {
    'plasma_cell_depletion': 'long-lived plasma cells (anti-BCMA, daratumumab) have no kill term in the model',
    'interferon_blockade': 'the model has no interferon state (anifrolumab-like drugs cannot be represented)',
    'treg_expansion': 'the model has no regulatory T cells (low-dose IL-2-like drugs cannot be represented)',
    'costimulation_blockade': 'no T-cell help or germinal centre (abatacept-, anti-CD40-like drugs cannot be represented)',
}


def apply_drug(p, spec):
    """LESSON 4 -- a drug is described by WHAT it does (mechanism) and HOW MUCH (strength). The strength is an
    input you choose; the model cannot tell you a new drug's strength. That number has to come from elsewhere
    (phase 1 biomarkers, a related drug, or a sweep)."""
    q = dict(p)
    for eff in spec['effects']:
        m, s = eff['mechanism'], float(eff['strength'])
        if m in UNSUPPORTED: raise ValueError(f"mechanism '{m}' unsupported: {UNSUPPORTED[m]}")
        if m == 'immunosuppression': q['s_inf'] = q['s_inf'] * (1 - s)
        # drug_on = 1 day, NOT 0: the baseline visit is measured before the first dose. Starting at day 0 made the
        # baseline itself drug-lowered and cancelled the effect in change-from-baseline (bug found by the first run).
        elif m == 'antiproteinuric': q['leak_red'] = s; q['drug_on'], q['drug_off'] = 1.0, 1e9
        elif m == 'haemodynamic': q['h_drug'] = s; q['drug_on'], q['drug_off'] = 1.0, 1e9
        elif m == 'bcell_depletion': q['dep_kill'] = s; q['dep_on'] = float(eff.get('start_day', 0)); q['dep_off'] = float(eff.get('end_day', 28))
        elif m == 'activity_clearance': q['k_clear'] = q['k_clear'] * (1 + s)
        else: raise ValueError(f"unknown mechanism '{m}'. Known: {list(MECHANISMS)}")
    return q


def run_trial(spec, ensemble, weeks):
    out = []
    for m in ensemble:
        p = base_params(m['s_inf'], m['k_dam'], m['k_tox'])
        c, Dc = trajectory(p, weeks, m['step']); t, Dt = trajectory(apply_drug(p, spec), weeks, m['step'])
        out.append(dict(k_tox=m['k_tox'], control=c.tolist(), drug=t.tolist(), D_control=Dc, D_drug=Dt))
    return out


def power(delta, sd, n):
    """LESSON 5 -- a predicted difference is only useful if a real trial could see it. Two-arm comparison of
    means, two-sided alpha 0.05: power = Phi(|delta| / (sd*sqrt(2/n)) - 1.96)."""
    return float(norm.cdf(abs(delta) / (sd * np.sqrt(2.0 / n)) - 1.96))


def report(spec, runs, weeks, n_per_arm, outdir):
    os.makedirs(outdir, exist_ok=True)
    o = load_aurora(); sd52_u, sd52_g = o[-1, 1], o[-1, 3]
    W = np.array(weeks); diff = np.array([np.array(r['drug']) - np.array(r['control']) for r in runs])   # members x weeks x 2
    i52 = list(weeks).index(52) if 52 in weeks else -1
    d52u, d52g = diff[:, i52, 0], diff[:, i52, 1]; dend_g = diff[:, -1, 1]
    by_kt = {}
    for kt in KTOX_GRID:
        sel = [j for j, r in enumerate(runs) if r['k_tox'] == kt]
        if sel: by_kt[f'{kt:.0e}'] = dict(dUPCR_wk52=float(np.median(d52u[sel])), dEGFR_wk52=float(np.median(d52g[sel])),
                                          dEGFR_end=float(np.median(dend_g[sel])))
    end_vals = [v['dEGFR_end'] for v in by_kt.values()]
    k_sensitive = (max(end_vals) - min(end_vals)) > 0.5 * max(1.0, abs(np.median(dend_g)))
    summ = dict(drug=spec['name'], members=len(runs), n_per_arm=n_per_arm,
                week52=dict(dUPCR_median=float(np.median(d52u)), dUPCR_range=[float(d52u.min()), float(d52u.max())],
                            dEGFR_median=float(np.median(d52g)), dEGFR_range=[float(d52g.min()), float(d52g.max())],
                            power_UPCR=power(np.median(d52u), sd52_u, n_per_arm), power_EGFR=power(np.median(d52g), sd52_g, n_per_arm)),
                end=dict(week=int(W[-1]), dEGFR_median=float(np.median(dend_g)), dEGFR_range=[float(dend_g.min()), float(dend_g.max())]),
                by_k_tox=by_kt, verdict_depends_on_unidentified_k_tox=bool(k_sensitive),
                caveats=['Arm means only: no responder rates, no individual patients.',
                         'Control arm calibrated to AURORA 1 placebo (MMF + steroids background, lupus nephritis). Not general SLE.',
                         'Drug strength is an input, not a prediction.',
                         'The model does not beat trivial forecasting rules on real data (rounds 50-53, 78, 79).',
                         'No interferon, Tregs, T-cell help or long-lived plasma-cell kill in the model.',
                         'Research and simulation only. Not clinical or drug-development guidance.'])
    json.dump(dict(spec=spec, summary=summ, weeks=weeks, runs=runs), open(os.path.join(outdir, 'result.json'), 'w'), indent=1)
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for k, (lab, obs_col) in enumerate((('UPCR change (mg/mg)', 0), ('eGFR change (mL/min/1.73m²)', 2))):
            C = np.array([r['control'] for r in runs])[:, :, k]; T = np.array([r['drug'] for r in runs])[:, :, k]
            for M, col, nm in ((C, '#6b7280', 'control (background therapy)'), (T, '#2563eb', spec['name'])):
                ax[k].fill_between(W, M.min(0), M.max(0), color=col, alpha=.2); ax[k].plot(W, np.median(M, 0), color=col, label=nm)
            ax[k].scatter(FIT_WEEKS, o[:, obs_col], s=12, color='k', label='AURORA 1 placebo (real)', zorder=3)
            ax[k].set_xlabel('week'); ax[k].set_title(lab); ax[k].axhline(0, lw=.5, color='k')
        ax[0].legend(fontsize=7); fig.suptitle('In-silico trial. Bands = ensemble range across unidentified parameters. Research only.', fontsize=9)
        fig.tight_layout(); fig.savefig(os.path.join(outdir, 'trajectories.png'), dpi=130)
    except Exception as e:
        print('plot skipped:', e)
    return summ


def main():
    ap = argparse.ArgumentParser(description='SLE/LN in-silico trial harness (research only)')
    ap.add_argument('drug', nargs='+', help='drug spec JSON file(s)')
    ap.add_argument('--weeks', type=int, default=156, help='follow-up length in weeks')
    ap.add_argument('--n', type=int, default=180, help='patients per arm for the power estimate')
    ap.add_argument('--boot', type=int, default=8)
    ap.add_argument('--out', default=os.path.join(HERE, 'runs'))
    a = ap.parse_args()
    ens_file = os.path.join(a.out, f'ensemble_v2_boot{a.boot}.json')
    if os.path.exists(ens_file):
        kept = json.load(open(ens_file))['kept']
    else:
        print('Building control-arm ensemble from real AURORA 1 placebo data ...', flush=True)
        kept, allm = build_ensemble(a.boot); os.makedirs(a.out, exist_ok=True)
        json.dump(dict(kept=kept, all=allm), open(ens_file, 'w'), indent=1)
    print(f'{len(kept)} ensemble members kept', flush=True)
    weeks = sorted(set(FIT_WEEKS + list(range(60, a.weeks + 1, 8)) + [a.weeks]))
    for f in a.drug:
        spec = json.load(open(f)); runs = run_trial(spec, kept, weeks)
        s = report(spec, runs, weeks, a.n, os.path.join(a.out, os.path.splitext(os.path.basename(f))[0]))
        print(json.dumps({k: v for k, v in s.items() if k != 'caveats'}, indent=1), flush=True)


if __name__ == '__main__':
    main()
