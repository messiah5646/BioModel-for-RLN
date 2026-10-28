"""
Shared RLN model core for rounds 40-44.

Six states: autoreactive B cells and short-lived plasma cells (both CD19+ and
depletable), long-lived plasma cells (CD19-, NOT depletable), autoantibody titre,
reversible inflammatory activity, and monotone irreversible renal damage.

Two structural commitments that the rounds depend on:
  * `Pl` is CD19-negative, so a CD19-directed therapy cannot reach it. This is what
    sets the autoantibody floor and the modelled mechanism of serological relapse.
  * `D` is monotone (dD/dt >= 0). The only published lupus model with a damage
    variable gives it a decay term, making it reversible.

`phi` is the damage -> antigen exposure -> B-cell activation feedback. phi=0 is the
feed-forward case; with phi=0 the kidney parameters move the molecular states by
exactly zero, which is round 40's structural-zero result.

Normalisation: every molecular compartment sits at 1.0 at the untreated steady state,
and the observation gains are set to AURORA-like baselines (UPCR 4.17, eGFR 78.2).

NOTE: this module has NO module-level side effects. An earlier version of the round-40
script printed and swept at import time, which cost ~6 s on every import and was a real
drag on the round-43 population runs. Scripts here guard everything behind __main__.
"""
import numpy as np
from scipy.integrate import solve_ivp

Y0 = np.array([1., 1., 1., 1., 1., 0.])
STATES = ('B', 'Ps', 'Pl', 'Ab', 'A', 'D')
BASELINE = np.array([4.17, 78.2, 1.0, 1.0])     # UPCR, eGFR, autoAb, B-cell


def make_params(phi=0.0):
    """Reference parameterisation. Steady state is exact; eGFR never hits its clamp."""
    d_B, d_Ps, d_Pl, d_Ab, k_clear = 0.030, 0.10, 0.0019, 0.033, 0.10
    return dict(
        d_B=d_B, d_Ps=d_Ps, d_Pl=d_Pl, d_Ab=d_Ab, k_clear=k_clear,
        s_B=d_B, f_s=d_Ps, f_l=d_Pl,            # -> B* = Ps* = Pl* = 1
        g_s=0.4 * d_Ab, g_l=0.6 * d_Ab,         # -> Ab* = 1
        k_a=k_clear,                            # -> A*  = 1
        kill=0.25,
        k_dam=6.0e-5, k_tox=1.2e-5,             # the degenerate pair
        p0=4.17, G0=92.0, c_inf=0.15,           # UPCR* = 4.17, eGFR* = 78.2
        s_inf=0.30, tau_s=45.0,
        phi=phi,
        # CD19-directed depletion window. OFF by default: the core carries no therapy,
        # and interventions are supplied as modules. The round-40/41 identifiability
        # analysis switches it on, because that analysis IS a depletion study -- the
        # long-lived-pool decay rate d_Pl is only identifiable from the depletion nadir.
        dep_kill=0.0, dep_on=180.0, dep_off=240.0,
    )


def supp(t, p):
    """Background immunosuppression (MMF/steroid-like): 1.0 -> s_inf."""
    return p['s_inf'] + (1.0 - p['s_inf']) * np.exp(-t / p['tau_s'])


def rhs(t, y, p):
    B, Ps, Pl, Ab, A, D = y
    Dc = min(max(D, 0.0), 0.95)
    N = 1.0 - Dc                       # surviving nephron fraction
    s = supp(t, p)
    fb = 1.0 + p['phi'] * Dc           # damage -> antigen -> B-cell activation
    # CD19-directed depletion reaches B and Ps only. Pl is CD19-negative.
    k = p['dep_kill'] if p['dep_on'] <= t <= p['dep_off'] else 0.0
    return [
        p['s_B'] * s * fb - p['d_B'] * B - k * B,
        p['f_s'] * s * B - p['d_Ps'] * Ps - k * Ps,
        p['f_l'] * s * B - p['d_Pl'] * Pl,
        p['g_s'] * Ps + p['g_l'] * Pl - p['d_Ab'] * Ab,
        p['k_a'] * Ab - p['k_clear'] * A,
        (p['k_dam'] * A + p['k_tox'] * p['p0'] * A * N) * N,   # monotone
    ]


def observe(y, p):
    """UPCR, eGFR, autoantibody, B-cell. UPCR is attenuated by nephron loss (the fold)."""
    B, Ps, Pl, Ab, A, D = y
    N = 1.0 - np.clip(D, 0.0, 0.95)
    return np.array([
        p['p0'] * A * N,
        p['G0'] * N * max(1.0 - p['c_inf'] * A, 0.05),
        Ab,
        B,
    ])


def simulate(p, tgrid, rtol=1e-11, atol=1e-13, y0=None):
    """Tight tolerances by default (needed for finite-difference sensitivity).

    For population runs use rtol=1e-7: verified to change UPCR_52 by 2e-9 relative
    against the tight solve. Round 43 v1 inherited 1e-11 here and died at its own
    300 s timeout, which is the same failure round 36 recorded.
    """
    sol = solve_ivp(rhs, (0.0, tgrid[-1]), Y0 if y0 is None else y0, args=(p,),
                    t_eval=tgrid, method='LSODA', rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y.T


def obs_matrix(p, tgrid, **kw):
    return np.array([observe(y, p) for y in simulate(p, tgrid, **kw)])


def obs_vector(p, tgrid, channels, **kw):
    """Baseline-standardised, channel-major flattening (used by the sensitivity code)."""
    O = obs_matrix(p, tgrid, **kw) / BASELINE
    return O[:, channels].T.ravel()


def sens(p, pname, tgrid, channels, rel_h=1e-2, **kw):
    """Central finite difference on a relative step.

    rel_h=1e-2 is converged here: |cos| between the k_dam and k_tox sensitivity
    vectors is stable to the 9th decimal across rel_h 1e-1 .. 1e-5. Smaller steps
    amplify solver roundoff -- the apparent molecular sensitivity at phi=0 scales as
    exactly 1/h, which is the signature of a true zero.
    """
    h = abs(p[pname]) * rel_h
    pp = dict(p); pp[pname] += h
    pm = dict(p); pm[pname] -= h
    return (obs_vector(pp, tgrid, channels, **kw)
            - obs_vector(pm, tgrid, channels, **kw)) / (2 * h)


def treated(p, cni_supp=0.45, haemo=0.25):
    """Stylised voclosporin-like arm: extra immunosuppression plus a reversible
    haemodynamic UPCR drop (route 3, which rounds 35/38 found the data demands).

    NOT fitted to AURORA. Its effect size (+9.8 pp at CV 0.70) is weaker than
    voclosporin's observed (+18 pp), so absolute required-N figures derived from it
    must not be quoted about real trials.
    """
    q = dict(p)
    q['s_inf'] = p['s_inf'] * (1.0 - cni_supp)
    q['p0'] = p['p0'] * (1.0 - haemo)
    return q


# AURORA-like visit grid: 13 visits in year 1, 6-monthly to month 36,
# plus fortnightly sampling around a day 180-240 depletion window.
VISIT_GRID = np.unique(np.concatenate([
    np.arange(0, 365, 28.0),
    np.arange(365, 1096, 182.0),
    np.arange(170, 261, 14.0),
]))
