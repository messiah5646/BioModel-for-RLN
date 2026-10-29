"""
Round 46. Prereg: docs/rounds/prereg_round46.md (written before any test number).

T1 additive vs multiplicative jitter -- falsifies round 44's own mechanism explanation
T2 is the viable CV window empty?
T3 can enrichment on BASELINE-MEASURABLE quantities recover effect size?
"""
import json, os
import numpy as np
from model import make_params, simulate, observe, treated

JIT = ['p0', 'd_Ab', 'k_clear', 'k_a', 'g_l', 'd_Pl']
TG = np.array([0., 364.])
RES = 'results/r46'


def base_params():
    p = make_params(phi=1.0); p['k_dam'] *= 10; p['k_tox'] *= 10
    return p


def patient(p):
    """Returns baseline UPCR, baseline eGFR, wk52 CRR, wk52 log UPCR, valid."""
    try:
        Y = simulate(p, TG, rtol=1e-7, atol=1e-9)
    except RuntimeError:
        return None
    O = np.array([observe(y, p) for y in Y])
    u0, g0, u, g = O[0, 0], O[0, 1], O[-1, 0], O[-1, 1]
    if g <= 0.05 * p['G0'] + 1e-9:
        return None
    crr = bool(u <= 0.5 and (g >= 60.0 or g >= 0.8 * g0))
    return u0, g0, crr, float(np.log(max(u, 1e-12)))


def draw(cv, rng, base, mode='mult'):
    """One patient's parameter vector. mode='add' matches the variance of 'mult'."""
    p = dict(base)
    if mode == 'mult':
        sd = np.sqrt(np.log(1 + cv ** 2))
        for nm in JIT:
            p[nm] = base[nm] * np.exp(rng.normal(0, sd) - 0.5 * sd * sd)
    else:
        # additive with the SAME variance as the multiplicative case: a lognormal
        # with CV has sd = cv*mean, so draw normal(mean, cv*mean) and floor at >0
        for nm in JIT:
            v = base[nm] + rng.normal(0, cv * abs(base[nm]))
            p[nm] = max(v, 1e-9 * abs(base[nm]))
    return p


def cohen_d(lu_c, lu_t):
    sp = np.sqrt((np.var(lu_c, ddof=1) + np.var(lu_t, ddof=1)) / 2)
    return (np.mean(lu_c) - np.mean(lu_t)) / sp if sp > 0 else 0.0


def req_n(p1, p2, alpha=0.05, power=0.8):
    from scipy.stats import norm
    if abs(p1 - p2) < 1e-9:
        return np.inf
    z = norm.ppf(1 - alpha / 2) + norm.ppf(power)
    return z ** 2 * (p1 * (1 - p1) + p2 * (1 - p2)) / (p1 - p2) ** 2


def run_cell(cv, seed, base, n, mode='mult'):
    """Virtual twins: same patient under control and treatment."""
    rng = np.random.default_rng(2000 + seed)
    rows = []
    for _ in range(n):
        p = draw(cv, rng, base, mode)
        rc = patient(p)
        rt = patient(treated(p))
        if rc is None or rt is None:
            continue
        rows.append((rc[0], rc[1], rc[2], rc[3], rt[2], rt[3]))
    return np.array(rows, dtype=float)


def summarise(A):
    """A columns: u0, g0, crr_c, logu_c, crr_t, logu_t"""
    pc, pt = A[:, 2].mean(), A[:, 4].mean()
    return dict(n=len(A), control=float(pc), treated=float(pt),
                diff_pp=float((pt - pc) * 100),
                cohen_d=float(cohen_d(A[:, 3], A[:, 5])),
                req_n=float(req_n(pc, pt)))


def main():
    os.makedirs('results', exist_ok=True)
    base = base_params()

    # ---------- T1: additive vs multiplicative ----------
    path = f'{RES}_t1.json'
    t1 = json.load(open(path)) if os.path.exists(path) else []
    done = {(r['mode'], r['cv'], r['seed']) for r in t1}
    for mode in ['mult', 'add']:
        for cv in [0.40, 0.70, 1.00]:
            for seed in range(4):
                if (mode, cv, seed) in done:
                    continue
                A = run_cell(cv, seed, base, 500, mode)
                s = summarise(A); s.update(mode=mode, cv=cv, seed=seed)
                t1.append(s); json.dump(t1, open(path, 'w'), indent=1)
                print(f"T1 {mode} cv={cv:.2f} s={seed} d={s['cohen_d']:.3f} "
                      f"diff={s['diff_pp']:+.1f}pp n={s['n']}", flush=True)

    # ---------- T2: is the window empty? ----------
    path = f'{RES}_t2.json'
    t2 = json.load(open(path)) if os.path.exists(path) else []
    done = {(r['cv'], r['seed']) for r in t2}
    for cv in [0.50, 0.55, 0.60, 0.65, 0.70, 0.75]:
        for seed in range(8):
            if (cv, seed) in done:
                continue
            A = run_cell(cv, seed, base, 500, 'mult')
            s = summarise(A); s.update(cv=cv, seed=seed)
            t2.append(s); json.dump(t2, open(path, 'w'), indent=1)
            print(f"T2 cv={cv:.2f} s={seed} ctrl={s['control']*100:.1f}% "
                  f"diff={s['diff_pp']:+.1f}pp", flush=True)

    # ---------- T3: enrichment on baseline-measurable quantities ----------
    path = f'{RES}_t3.json'
    if not os.path.exists(path):
        A = run_cell(0.70, 0, base, 4000, 'mult')
        u0, g0 = A[:, 0], A[:, 1]
        uq = np.quantile(u0, [1/3, 2/3]); uqq = np.quantile(u0, [0.25, 0.75])
        um = np.median(u0)
        RULES = {
            'none (all patients)':            np.ones(len(A), bool),
            'baseline UPCR bottom tertile':   u0 <= uq[0],
            'baseline UPCR middle tertile':   (u0 > uq[0]) & (u0 <= uq[1]),
            'baseline UPCR top tertile':      u0 > uq[1],
            'baseline UPCR bottom quartile':  u0 <= uqq[0],
            'baseline UPCR top quartile':     u0 > uqq[1],
            'baseline eGFR >= 60':            g0 >= 60.0,
            'baseline eGFR < 60':             g0 < 60.0,
            'UPCR<=med & eGFR>=60':           (u0 <= um) & (g0 >= 60.0),
            'UPCR<=med & eGFR<60':            (u0 <= um) & (g0 < 60.0),
            'UPCR>med & eGFR>=60':            (u0 > um) & (g0 >= 60.0),
            'UPCR>med & eGFR<60':             (u0 > um) & (g0 < 60.0),
        }
        out = []
        for name, m in RULES.items():
            if m.sum() < 50:
                out.append(dict(rule=name, n=int(m.sum()), skipped=True)); continue
            s = summarise(A[m]); s.update(rule=name, frac_eligible=float(m.mean()))
            out.append(s)
            print(f"T3 {name:32s} n={s['n']:5d} ctrl={s['control']*100:5.1f}% "
                  f"diff={s['diff_pp']:+5.1f}pp d={s['cohen_d']:.3f} "
                  f"reqN={s['req_n']:7.0f}", flush=True)
        json.dump(out, open(path, 'w'), indent=1)
    print("ALLDONE")


if __name__ == '__main__':
    main()
