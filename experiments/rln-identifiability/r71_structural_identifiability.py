"""Round 71: local structural identifiability of the FULL six-state RLN model (experiments/rln-identifiability/model.py,
commit 2e76e7c), smooth version (no clamps, no drug window), by exact Taylor-series sensitivity rank mod a prime.
Outputs' Taylor coefficients at t=0 are differentiated wrt every unknown (dual numbers); generic rank of that Jacobian
= number of locally identifiable parameter combinations. A parameter is identifiable iff deleting its column drops rank.
Null vectors are reported in log-coordinates (v_i * theta_i), where scaling symmetries show up as small integers."""
import random, sys, json
from fractions import Fraction
Pm = (1 << 61) - 1
inv = lambda a: pow(a % Pm, Pm - 2, Pm)

class S:                                    # truncated series with gradient; c[k] = [value, d/dtheta_1, ...]
    K = 0; NP = 0
    def __init__(s, c=None): s.c = c if c is not None else [[0] * (S.NP + 1) for _ in range(S.K + 1)]
    @staticmethod
    def const(v, idx=None):
        o = S(); o.c[0][0] = v % Pm
        if idx is not None: o.c[0][idx + 1] = 1
        return o
    def __add__(a, b):
        if not isinstance(b, S): b = S.const(b)
        return S([[(x + y) % Pm for x, y in zip(ra, rb)] for ra, rb in zip(a.c, b.c)])
    __radd__ = __add__
    def __neg__(a): return S([[(-x) % Pm for x in r] for r in a.c])
    def __sub__(a, b): return a + (-(b if isinstance(b, S) else S.const(b)))
    def __rsub__(a, b): return S.const(b) - a
    def __mul__(a, b):
        if not isinstance(b, S): b = S.const(b)
        o = S(); n = S.NP + 1
        for k in range(S.K + 1):
            row = [0] * n
            for j in range(k + 1):
                ra, rb = a.c[j], b.c[k - j]; a0, b0 = ra[0], rb[0]
                if a0 == 0 and b0 == 0 and not any(ra) and not any(rb): continue
                row[0] += a0 * b0
                for i in range(1, n): row[i] += ra[i] * b0 + a0 * rb[i]
            o.c[k] = [x % Pm for x in row]
        return o
    __rmul__ = __mul__
    def integ(a, x0):                       # x0 + int_0^t a
        o = S(); o.c[0] = list(x0.c[0])
        for k in range(S.K): iv = inv(k + 1); o.c[k + 1] = [x * iv % Pm for x in a.c[k]]
        return o

def rank_mod(M):
    M = [r[:] for r in M]; rk = 0; piv = []
    for col in range(len(M[0])):
        pr = next((r for r in range(rk, len(M)) if M[r][col] % Pm), None)
        if pr is None: continue
        M[rk], M[pr] = M[pr], M[rk]; iv = inv(M[rk][col]); M[rk] = [x * iv % Pm for x in M[rk]]
        for r in range(len(M)):
            if r != rk and M[r][col] % Pm:
                f = M[r][col]; M[r] = [(x - f * y) % Pm for x, y in zip(M[r], M[rk])]
        piv.append(col); rk += 1
    return rk, piv, M

def ratrec(a, bound=40):                    # small rational reconstruction mod Pm
    for d in range(1, bound + 1):
        n = a * d % Pm
        if n <= bound: return Fraction(n, d)
        if Pm - n <= bound: return Fraction(-(Pm - n), d)
    return None

RATE = ['d_B','d_Ps','d_Pl','d_Ab','k_clear','s_B','f_s','f_l','g_s','g_l','k_a','k_dam','k_tox','p0','G0','c_inf','s_inf','tau_s','phi','k_prog']
ICS = ['B0','Ps0','Pl0','Ab0','A0','D0']

def jac(unknown, experiments, outputs, K, seed, tied=False):
    """experiments: list of dicts(dep=bool, depkill_unknown=bool). Returns Jacobian rows and names."""
    random.seed(seed); names = list(unknown)
    if any(e.get('dep') for e in experiments) and 'dep_kill' not in names and experiments[0].get('depkill_unknown', False): names.append('dep_kill')
    S.K = K; S.NP = len(names); val = {n: random.randint(2, 10**6) for n in RATE + ICS + ['dep_kill']}
    def P(n):
        return S.const(val[n], names.index(n)) if n in names else S.const(val[n])
    p = {n: P(n) for n in RATE + ICS + ['dep_kill']}
    if tied:                                # model.py normalisation: s_B=d_B, f_s=d_Ps, f_l=d_Pl, k_a=k_clear, g_s+g_l=d_Ab, ICs=1, D0=0
        p['s_B'], p['f_s'], p['f_l'], p['k_a'] = p['d_B'], p['d_Ps'], p['d_Pl'], p['k_clear']; p['g_l'] = p['d_Ab'] - p['g_s']
        for n in ICS: p[n] = S.const(0 if n == 'D0' else 1)
    rows = []; itau = None
    for e in experiments:
        k = p['dep_kill'] if e.get('dep') else S.const(0)
        x0 = [p[n] for n in ICS] + [S.const(1)]; x = list(x0)
        # 1/tau_s as a series constant with gradient: d(1/tau)/dtau = -1/tau^2
        it = S.const(inv(val['tau_s']));
        if 'tau_s' in names: it.c[0][names.index('tau_s') + 1] = (-inv(val['tau_s']) ** 2) % Pm
        for _ in range(K):
            B, Ps, Pl, Ab, A, D, E = x; N = 1 - D; s = p['s_inf'] + (1 - p['s_inf']) * E; fb = 1 + p['phi'] * D
            f = [p['s_B'] * s * fb - p['d_B'] * B - k * B, p['f_s'] * s * B - p['d_Ps'] * Ps - k * Ps, p['f_l'] * s * B - p['d_Pl'] * Pl,
                 p['g_s'] * Ps + p['g_l'] * Pl - p['d_Ab'] * Ab, p['k_a'] * Ab - p['k_clear'] * A,
                 (p['k_dam'] * A + p['k_tox'] * p['p0'] * A * N + p['k_prog']) * N, -(it * E)]
            x = [fi.integ(x0i) for fi, x0i in zip(f, x0)]
        B, Ps, Pl, Ab, A, D, E = x; N = 1 - D
        obs = dict(UPCR=p['p0'] * A * N, eGFR=p['G0'] * N * (1 - p['c_inf'] * A), Ab=Ab, B=B)
        for o in outputs:
            for kk in range(K + 1): rows.append(obs[o].c[kk][1:])
    return rows, names, val

def analyse(tag, unknown, experiments, outputs, K=None, tied=False):
    res = []
    for seed in (0, 1):
        ranks = {}
        for KK in ((K, K + 3) if K else (len(unknown) + 6, len(unknown) + 9)):
            rows, names, val = jac(unknown, experiments, outputs, KK, seed, tied); rk, piv, R = rank_mod(rows); ranks[KK] = rk
        ident = []
        for i, n in enumerate(names):
            sub = [r[:i] + r[i + 1:] for r in rows]; ident.append(rank_mod(sub)[0] == rk - 1)
        # null space in log coordinates
        free = [c for c in range(len(names)) if c not in piv]; nulls = []
        for fc in free:
            v = [0] * len(names); v[fc] = 1
            for r, pc in enumerate(piv): v[pc] = (-R[r][fc]) % Pm
            lv = [vi * val[n] % Pm for vi, n in zip(v, names)]; nz = [i for i, x in enumerate(lv) if x]
            sc = inv(lv[fc]); lv = [x * sc % Pm for x in lv]
            rr = {names[i]: ratrec(lv[i]) for i in nz}
            nulls.append({k_: (str(v_) if v_ is not None else 'nonrational') for k_, v_ in rr.items()})
        res.append(dict(seed=seed, n_unknown=len(names), rank_by_order=ranks, non_identifiable=[n for n, ok in zip(names, ident) if not ok], null_log=nulls))
    out = dict(tag=tag, outputs=outputs, experiments=experiments, tied=tied, runs=res); print(json.dumps(out), flush=True); return out

if __name__ == '__main__':
    E1 = [dict(dep=False)]; E12 = [dict(dep=False), dict(dep=True, depkill_unknown=False)]
    ALL4 = ['UPCR', 'eGFR', 'Ab', 'B']; CL = ['UPCR', 'eGFR']
    TIED = ['d_B','d_Ps','d_Pl','d_Ab','k_clear','g_s','k_dam','k_tox','p0','G0','c_inf','s_inf','tau_s','phi','k_prog']
    out = []
    out.append(analyse('A tied normalisation, 4 outputs, natural course', TIED, E1, ALL4, tied=True))
    out.append(analyse('B tied, clinical outputs only', TIED, E1, CL, tied=True))
    out.append(analyse('C tied, 4 outputs, + depletion experiment (kill known)', TIED, E12, ALL4, tied=True))
    out.append(analyse('D free rates, ICs known, 4 outputs', RATE, E1, ALL4))
    out.append(analyse('E free rates + unknown ICs, 4 outputs', RATE + ICS, E1, ALL4))
    out.append(analyse('F free rates, ICs known, 4 outputs, + depletion', RATE, E12, ALL4))
    json.dump(out, open('ident_full.json', 'w'), indent=1)
