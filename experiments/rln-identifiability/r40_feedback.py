"""
Round 40/41: does a REAL dynamic feedback path restore identifiability of the
(k_dam, k_tox) damage split?

Round 40 proved that under feed-forward coupling an upstream observable carries
EXACTLY zero information about a downstream degeneracy, and correctly refused to
count its own s_B-shift emulation as a feedback test. This writes the feedback into
the ODE (phi) and measures it two ways.

Two instruments, because one was not enough:
  * cosine-VIF = 1/(1-cos^2). NOT monotone in the observable set -- it reports
    "gains" below 1.0, i.e. adding a channel appearing to destroy information, which
    is impossible. Retained only to reproduce round 40's structural zero.
  * Fisher marginal variance. Monotone by construction. This one governs.
"""
import json
import numpy as np
from model import make_params, simulate, obs_vector, sens, VISIT_GRID

CHANNELS = {'S0': [0, 1], 'S1': [0, 1, 2], 'S2': [0, 1, 3], 'S3': [0, 1, 2, 3]}
PARS = ['k_dam', 'k_tox', 'd_Pl', 'f_l', 'dep_kill', 'c_inf', 'p0', 'k_clear', 'd_Ab']
NOISE = np.array([0.15, 0.08, 0.20, 0.25])      # CV: UPCR, eGFR, autoAb, B-cell


def vif_pair(p, channels, rel_h=1e-2):
    sa = sens(p, 'k_dam', VISIT_GRID, channels, rel_h)
    sb = sens(p, 'k_tox', VISIT_GRID, channels, rel_h)
    c = min(abs(float(sa @ sb) / (np.linalg.norm(sa) * np.linalg.norm(sb))), 1.0)
    return c, 1.0 / max(1.0 - c * c, 1e-300)


def marginal_sigma(p, channels, rel_h=1e-2):
    """Marginal sigma on LOG parameters from the full FIM with nuisance parameters."""
    cols = []
    for nm in PARS:
        h = abs(p[nm]) * rel_h
        pp = dict(p); pp[nm] += h
        pm = dict(p); pm[nm] -= h
        d = (obs_vector(pp, VISIT_GRID, channels)
             - obs_vector(pm, VISIT_GRID, channels)) / (2 * h) * p[nm]
        cols.append(d)
    J = np.array(cols).T
    sig = np.repeat(NOISE[channels], len(VISIT_GRID))   # obs_vector is channel-major
    F = J.T @ np.diag(1.0 / sig ** 2) @ J
    return np.sqrt(np.clip(np.diag(np.linalg.pinv(F, rcond=1e-12)), 0, None))


def main():
    print("STEP-SIZE CONVERGENCE at phi=0 (expect the structural zero)")
    for rel_h in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]:
        p = make_params(0.0); p['dep_kill'] = 0.25
        c, v = vif_pair(p, [0, 1, 2, 3], rel_h)
        s = sens(p, 'k_dam', VISIT_GRID, [2], rel_h)
        print(f"  rel_h={rel_h:7.0e} |cos|={c:.10f} VIF={v:11.1f} "
              f"max|dAb/dk_dam|={np.max(np.abs(s)):.3e}")
    print("  -> the apparent sensitivity scales as exactly 1/h: a true zero, not an effect.\n")

    print("FISHER MARGINAL SIGMA on log k_dam (monotone instrument; this one governs)")
    print(f"{'sev':>4} {'phi':>5} {'D_fin':>7} | {'S0':>9} {'S3':>9} | {'gain':>6}")
    out = []
    for sev in [1, 10, 20]:
        for phi in [0.0, 5.0, 20.0]:
            p = make_params(phi); p['k_dam'] *= sev; p['k_tox'] *= sev
            p['dep_kill'] = 0.25          # this analysis is a CD19-depletion study
            Df = simulate(p, VISIT_GRID)[-1, 5]
            s0 = marginal_sigma(p, [0, 1]); s3 = marginal_sigma(p, [0, 1, 2, 3])
            print(f"{sev:4d} {phi:5.1f} {Df:7.4f} | {s0[0]:9.3f} {s3[0]:9.3f} | "
                  f"{s0[0]/s3[0]:6.3f}")
            out.append(dict(sev=sev, phi=phi, D_final=float(Df),
                            sigma_S0=float(s0[0]), sigma_S3=float(s3[0]),
                            gain=float(s0[0] / s3[0])))
    json.dump(out, open('results/r40_feedback.json', 'w'), indent=1)
    worst = max(r["sigma_S0"] for r in out)
    best = min(r["sigma_S0"] for r in out)
    best_gain = max(r["gain"] for r in out)
    print(f"\n  Severity is worth {worst/best:.0f}x on sigma(k_dam); "
          f"the observable set at most {best_gain:.2f}x.")


if __name__ == '__main__':
    main()
