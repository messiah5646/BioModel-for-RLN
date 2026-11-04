"""
R67: WHY is the exponent -2? Decomposing R66's scaling law into its sources.

R66 established sigma(log k_dam) = c * D_final^(-1.989), r2 = 0.9997, validated out of
sample to 3.3% median error. It left the exponent UNEXPLAINED and offered a guess:

  "[HYPOTHESIS, not tested] the extra power plausibly comes from the k_dam/k_tox
   degeneracy -- R38 measured their correlation at r = -0.982, and the marginal sigma
   of a near-collinear pair degrades faster than either parameter's own sensitivity.
   One power from sensitivity, one from conditioning."

This tests that guess by splitting sigma into two factors that multiply:

  sigma_marginal  = sigma_alone  x  inflation
  (all others estimated)  (all others KNOWN)   (price of not knowing them)

  sigma_alone = 1/sqrt(F_kk)            -- k_dam's own Fisher information
  inflation   = sigma_marginal / sigma_alone

Taking logs, the exponents ADD. So if the guess is right:
  exponent(sigma_alone) ~= -1   and   exponent(inflation) ~= -1.

PRE-REGISTRATION, before computing. Contamination: I know the total is -1.989 and I
proposed the 1+1 split, so the SPLIT is the hypothesis under test, not a prediction of
the total.

  H1. exponent(sigma_alone) is in [-1.25, -0.75].   Expected: ~-1, the trivial
      "bigger effect, better estimate" scaling.
  H2. exponent(inflation) is in [-1.25, -0.75], and the two exponents sum to the
      total within 0.05.  Expected: ~-1.
  H3. k_tox is the single largest contributor to the inflation. Operationalised as a
      leave-one-out: moving each nuisance parameter from "estimated" to "known" one at
      a time, dropping k_tox reduces the inflation MORE than dropping any other.

  Falsifiers, stated so the round can fail:
   - H1 outside its band -> the sensitivity half is NOT the trivial scaling, and R66's
     framing of "-1 would be trivial" was wrong.
   - H2 outside its band -> the 1+1 decomposition is wrong, whatever H1 does.
   - H3 fails -> the degeneracy is NOT k_tox-driven and the inflation comes from the
     wider nuisance set, which would make the mechanism a different one.

A pass on all three turns R66's empirical regularity into a structural account.
A failure on H2 leaves the law standing and the explanation dead -- the law is
out-of-sample validated and does not depend on this round.
"""
import json
import numpy as np
from model import make_params, simulate, obs_vector, VISIT_GRID
from r40_feedback import PARS, NOISE, CHANNELS

SEV = [1, 2, 4, 6, 10, 14, 20, 30]
PHI = 0.0
OBS = CHANNELS['S0']            # clinical only, as in the headline claim
K = PARS.index('k_dam')
T = PARS.index('k_tox')


def fim(p, channels, rel_h=1e-2):
    """Identical construction to r40_feedback.marginal_sigma, but returns F."""
    cols = []
    for nm in PARS:
        h = abs(p[nm]) * rel_h
        pp = dict(p); pp[nm] += h
        pm = dict(p); pm[nm] -= h
        d = (obs_vector(pp, VISIT_GRID, channels)
             - obs_vector(pm, VISIT_GRID, channels)) / (2 * h) * p[nm]
        cols.append(d)
    J = np.array(cols).T
    sig = np.repeat(NOISE[channels], len(VISIT_GRID))
    return J.T @ np.diag(1.0 / sig ** 2) @ J


def setup(sev, phi=PHI):
    p = make_params(phi)
    p['k_dam'] *= sev
    p['k_tox'] *= sev
    p['dep_kill'] = 0.25
    return p


def marg(F, idx):
    """Marginal sigma on parameter K when the parameters in idx are jointly estimated."""
    sub = F[np.ix_(idx, idx)]
    pos = idx.index(K)
    return float(np.sqrt(np.linalg.pinv(sub, rcond=1e-12)[pos, pos]))


def expo(D, S):
    b, a = np.polyfit(np.log(D), np.log(S), 1)
    pred = np.exp(a + b * np.log(D))
    r2 = 1 - np.sum((np.log(S) - np.log(pred)) ** 2) / \
        np.sum((np.log(S) - np.log(S).mean()) ** 2)
    return float(b), float(np.exp(a)), float(r2)


print("R67 -- decomposing R66's -2 exponent. Prereg in this file's docstring.\n")
rows = []
for sev in SEV:
    p = setup(sev)
    F = fim(p, OBS)
    D = float(simulate(p, VISIT_GRID)[-1, 5])
    s_all = marg(F, list(range(len(PARS))))          # everything estimated
    s_alone = float(1.0 / np.sqrt(F[K, K]))          # everything else KNOWN
    s_pair = marg(F, [K, T])                         # only k_tox unknown
    rows.append(dict(sev=sev, D=D, s_all=s_all, s_alone=s_alone, s_pair=s_pair,
                     infl=s_all / s_alone, infl_pair=s_pair / s_alone))

print(f"{'sev':>4} {'D_final':>8} {'sig_alone':>10} {'sig_pair':>10} {'sig_all':>10}"
      f" {'infl':>8} {'infl_pair':>10}")
print("-" * 68)
for r in rows:
    print(f"{r['sev']:4d} {r['D']:8.4f} {r['s_alone']:10.4f} {r['s_pair']:10.4f}"
          f" {r['s_all']:10.4f} {r['infl']:8.2f} {r['infl_pair']:10.2f}")

D = np.array([r['D'] for r in rows])
b_all, c_all, r2_all = expo(D, np.array([r['s_all'] for r in rows]))
b_alone, c_alone, r2_alone = expo(D, np.array([r['s_alone'] for r in rows]))
b_infl, c_infl, r2_infl = expo(D, np.array([r['infl'] for r in rows]))
b_pair, _, r2_pair = expo(D, np.array([r['infl_pair'] for r in rows]))

print(f"\nexponents in D_final (8 severities, phi=0):")
print(f"  sigma_all      (total)         b = {b_all:+.3f}   r2 = {r2_all:.5f}")
print(f"  sigma_alone    (sensitivity)   b = {b_alone:+.3f}   r2 = {r2_alone:.5f}")
print(f"  inflation      (conditioning)  b = {b_infl:+.3f}   r2 = {r2_infl:.5f}")
print(f"  inflation from k_tox alone     b = {b_pair:+.3f}   r2 = {r2_pair:.5f}")
print(f"\n  sum of parts  {b_alone:+.3f} {b_infl:+.3f} = {b_alone + b_infl:+.3f}"
      f"   vs total {b_all:+.3f}   (residual {abs(b_alone + b_infl - b_all):.4f})")

h1 = -1.25 <= b_alone <= -0.75
h2 = (-1.25 <= b_infl <= -0.75) and abs(b_alone + b_infl - b_all) < 0.05
print(f"\n  H1 (sensitivity exponent in [-1.25,-0.75]): "
      f"{'SUPPORTED' if h1 else 'FALSIFIED'}  (b = {b_alone:+.3f})")
print(f"  H2 (inflation exponent in [-1.25,-0.75] and parts sum to total): "
      f"{'SUPPORTED' if h2 else 'FALSIFIED'}  (b = {b_infl:+.3f})")

# ---------------- H3: which nuisance parameter drives the inflation ----------------
print("\nH3 -- leave-one-out: move each nuisance parameter from ESTIMATED to KNOWN")
p = setup(10)
F = fim(p, OBS)
full = marg(F, list(range(len(PARS))))
alone = float(1.0 / np.sqrt(F[K, K]))
loo = []
for j, nm in enumerate(PARS):
    if j == K:
        continue
    idx = [i for i in range(len(PARS)) if i != j]
    s = marg(F, idx)
    loo.append((nm, s, (full - s) / (full - alone)))
loo.sort(key=lambda x: -x[2])
print(f"  (at sev=10: sigma all-estimated {full:.4f}, all-others-known {alone:.4f},"
      f" inflation {full/alone:.1f}x)")
print(f"  {'dropped param':>12} {'sigma':>10} {'% of inflation removed':>24}")
for nm, s, frac in loo:
    print(f"  {nm:>12} {s:10.4f} {frac:23.1%}")
h3 = loo[0][0] == 'k_tox'
print(f"\n  H3 (k_tox is the largest single contributor): "
      f"{'SUPPORTED' if h3 else 'FALSIFIED'}  (largest = {loo[0][0]},"
      f" {loo[0][2]:.1%})")

json.dump(dict(rows=rows, b_all=b_all, b_alone=b_alone, b_infl=b_infl, b_pair=b_pair,
               r2_all=r2_all, r2_alone=r2_alone, r2_infl=r2_infl,
               loo=[dict(param=n, sigma=s, frac=f) for n, s, f in loo],
               H1=bool(h1), H2=bool(h2), H3=bool(h3)),
          open('results/r67_exponent_decomposition.json', 'w'), indent=1)
print("\nsaved -> results/r67_exponent_decomposition.json")
