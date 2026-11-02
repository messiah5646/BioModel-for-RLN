"""
R61: Stage 1's gate on the B-cell axis. Prereg: claude/prereg-round61.md.

Fit weeks 2-24 (8 points, both arms), predict week 52 (2 points, never used).
Four free parameters, declared before fitting: kill, t_off, d_B, s_scale.

Observable is percent change from baseline in CD19+ B cells. The CSV's baseline row is
in absolute cells/uL and is NEVER used as a data point (round 58 TRAP 1); B(0) is the
model's own initial condition.

Obinutuzumab is anti-CD20: it depletes B. It does not touch Pl (CD20-negative) and this
endpoint does not measure Ab, so only the B compartment is integrated here.

SOLVER NOTE: a first version used solve_ivp with max_step=5 and did not finish in 12
minutes -- the depletion window's discontinuity at t_off forced the step size down. The
B equation is LINEAR, so it is solved in closed form per segment here instead. This is
exact, not an approximation, and is verified against solve_ivp at the fitted parameters
in the residual-check block at the end.
"""
import csv, json, math
import numpy as np
from scipy.optimize import least_squares

VISITS = {'week_2': 14.0, 'week_4': 28.0, 'week_12': 84.0,
          'week_24': 168.0, 'week_52': 364.0}      # days
FIT = ['week_2', 'week_4', 'week_12', 'week_24']
HELD = 'week_52'
ARMS = ['obinutuzumab_1000mg', 'placebo']
TAU_S = 45.0          # background immunosuppression time constant, from model.py


def load():
    with open('../../data/real/obinutuzumab_nobility_per_visit.csv') as f:
        rows = list(csv.DictReader([l for l in f if not l.lstrip().startswith('#')]))
    D = {}
    for r in rows:
        if r['visit'] not in VISITS or r['unit_note'] != 'percent_change':
            continue
        if r['bcell_value'] == 'not reported' or r['bcell_sd'] == 'not reported':
            continue
        D.setdefault(r['arm'], {})[r['visit']] = dict(
            t=VISITS[r['visit']], v=float(r['bcell_value']),
            sd=float(r['bcell_sd']), n=int(r['bcell_n']))
    return D


D = load()


def _P(t, a, d_B, s_inf):
    """Particular solution of B' + aB = d_B*(s_inf + (1-s_inf)e^{-t/TAU_S})."""
    lam = 1.0 / TAU_S
    denom = a - lam
    if abs(denom) < 1e-9:
        denom = math.copysign(1e-9, denom if denom != 0 else 1.0)
    return d_B * s_inf / a + d_B * (1.0 - s_inf) * math.exp(-t * lam) / denom


def _seg(B0, t0, t1, a, d_B, s_inf):
    """Exact linear-ODE step from t0 to t1 at constant kill rate."""
    return math.exp(-a * (t1 - t0)) * (B0 - _P(t0, a, d_B, s_inf)) + _P(t1, a, d_B, s_inf)


def bcell_curve(arm, x, ts):
    """B(t) with B(0)=1. x = [log kill, t_off, log d_B, logit s_scale]."""
    kill, t_off = math.exp(x[0]), x[1]
    d_B = math.exp(x[2])
    s_inf = 1.0 / (1.0 + math.exp(-x[3]))
    treated = (arm == 'obinutuzumab_1000mg')
    t_off = max(0.0, t_off)

    out = []
    for t in ts:
        B, cur = 1.0, 0.0
        if treated and t_off > 0.0:
            b1 = min(t, t_off)
            B = _seg(B, 0.0, b1, d_B + kill, d_B, s_inf)
            cur = b1
        if t > cur:
            B = _seg(B, cur, t, d_B, d_B, s_inf)
        out.append(B)
    return np.array(out)


def residuals(x):
    r = []
    for arm in ARMS:
        ts = [D[arm][v]['t'] for v in FIT]
        B = bcell_curve(arm, x, ts)
        for v, b in zip(FIT, B):
            o = D[arm][v]
            sem = o['sd'] / math.sqrt(o['n'])
            r.append((100.0 * (b - 1.0) - o['v']) / sem)
    return r


best = None
for seed in range(60):
    rng = np.random.default_rng(seed)
    x0 = np.array([np.log(10 ** rng.uniform(-2, 0.5)), rng.uniform(20, 300),
                   np.log(10 ** rng.uniform(-3, -1)), rng.normal(0, 1.5)])
    try:
        res = least_squares(residuals, x0, method='lm', max_nfev=20000)
    except Exception:
        continue
    if best is None or res.cost < best.cost:
        best = res

kill, t_off = math.exp(best.x[0]), best.x[1]
d_B, s_inf = math.exp(best.x[2]), 1 / (1 + math.exp(-best.x[3]))
print("R61 -- Stage 1 gate on the B-cell axis (prereg: prereg-round61.md)")
print(f"fit cost {best.cost:.4f} on 8 points | kill={kill:.4f}/d  t_off={t_off:.1f}d  "
      f"d_B={d_B:.5f}/d  s_inf={s_inf:.3f}")

fit_chi2 = 2.0 * best.cost
print(f"\nSTOPPING RULE (prereg sec 6): chi2 on fit points = {fit_chi2:.3f}, dof = 4")
if fit_chi2 > 80.0:
    print("  -> fit did NOT reproduce the fit points. Round stops here.")
    raise SystemExit(0)
print("  -> fit reproduces the fit points; proceeding to held-out comparison.")

print(f"\n{'arm':>22} | {'visit':>8} | {'observed':>9} | {'model':>9} | {'resid/SEM':>9}")
print("-" * 70)
for arm in ARMS:
    ts = [D[arm][v]['t'] for v in FIT]
    B = bcell_curve(arm, best.x, ts)
    for v, b in zip(FIT, B):
        o = D[arm][v]
        sem = o['sd'] / math.sqrt(o['n'])
        print(f"{arm:>22} | {v:>8} | {o['v']:9.3f} | {100*(b-1):9.3f} | "
              f"{(100*(b-1)-o['v'])/sem:9.2f}")

errs = {}
for arm in ARMS:
    o = D[arm][HELD]
    model = 100.0 * (bcell_curve(arm, best.x, [o['t']])[0] - 1.0)
    pers = D[arm]['week_24']['v']
    v1, v2 = D[arm]['week_12']['v'], D[arm]['week_24']['v']
    t1, t2 = D[arm]['week_12']['t'], D[arm]['week_24']['t']
    lin = v2 + (v2 - v1) / (t2 - t1) * (o['t'] - t2)
    mean4 = float(np.mean([D[arm][v]['v'] for v in FIT]))
    errs[arm] = dict(obs=o['v'], sem=o['sd'] / math.sqrt(o['n']),
                     model_pred=model, model=abs(model - o['v']),
                     T0=abs(pers - o['v']), T1=abs(lin - o['v']), T2=abs(mean4 - o['v']))

print(f"\nHELD OUT week 52 (2 points, never used in fitting)")
print(f"{'arm':>22} | {'observed':>9} | {'model':>9} | {'|e| model':>10} | "
      f"{'|e| T0':>8} | {'|e| T1':>8} | {'|e| T2':>8} | {'1 SEM':>7}")
print("-" * 100)
for arm in ARMS:
    e = errs[arm]
    print(f"{arm:>22} | {e['obs']:9.3f} | {e['model_pred']:9.3f} | {e['model']:10.3f} | "
          f"{e['T0']:8.3f} | {e['T1']:8.3f} | {e['T2']:8.3f} | {e['sem']:7.3f}")

rmse = {k: math.sqrt(np.mean([errs[a][k] ** 2 for a in ARMS]))
        for k in ('model', 'T0', 'T1', 'T2')}
print(f"\npooled RMSE  model {rmse['model']:.3f} | persistence(T0) {rmse['T0']:.3f} | "
      f"linear(T1) {rmse['T1']:.3f} | armmean(T2) {rmse['T2']:.3f}")

per_arm_win = all(errs[a]['model'] < errs[a]['T0'] for a in ARMS)
h1 = (rmse['model'] < rmse['T0']) and per_arm_win
print(f"\nH1 (model beats persistence on pooled RMSE *and* in BOTH arms): "
      f"{'SUPPORTED' if h1 else 'FALSIFIED'}")
if (rmse['model'] < rmse['T0']) and not per_arm_win:
    print("  NOTE: pooled win WITHOUT a per-arm win -> reported as a WASH per prereg sec 5.")
t = errs['obinutuzumab_1000mg']
h2 = t['model'] < t['sem']
print(f"H2 (treated-arm error within 1 SEM = {t['sem']:.3f} pp): "
      f"{'SUPPORTED' if h2 else 'FALSIFIED'}  (model {t['model']:.3f} pp = "
      f"{t['model']/t['sem']:.2f} SEM; persistence was {t['T0']/t['sem']:.2f} SEM)")

# ---- solver cross-check: closed form vs solve_ivp at the fitted parameters ----
from scipy.integrate import solve_ivp
def _num(arm, t):
    treated = (arm == 'obinutuzumab_1000mg')
    def rhs(tt, y):
        k = kill if (treated and 0.0 <= tt <= t_off) else 0.0
        s = s_inf + (1 - s_inf) * math.exp(-tt / TAU_S)
        return [d_B * s - d_B * y[0] - k * y[0]]
    sol = solve_ivp(rhs, (0, t), [1.0], t_eval=[t], rtol=1e-10, atol=1e-13,
                    max_step=2.0)
    return sol.y[0][-1]
print("\nsolver cross-check (closed form vs solve_ivp) at week 52:")
for arm in ARMS:
    cf = bcell_curve(arm, best.x, [364.0])[0]
    nm = _num(arm, 364.0)
    print(f"  {arm:>22}: closed {cf:.9f} | numeric {nm:.9f} | diff {abs(cf-nm):.2e}")

json.dump(dict(x=[float(v) for v in best.x], fit_cost=float(best.cost),
               heldout=errs, rmse=rmse, H1=bool(h1), H2=bool(h2)),
          open('results/r61_bcell_gate.json', 'w'), indent=1)
print("\nsaved -> results/r61_bcell_gate.json")
