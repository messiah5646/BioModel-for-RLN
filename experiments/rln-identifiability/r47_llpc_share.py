"""
R47: attack round 40/44's most striking claim -- "89% B-cell depletion buys a 7.9%
autoantibody fall, because 99.7% of production comes from the protected CD19- pool."

Is that a RESULT or a consequence of an ASSUMED parameter?

model.py sets g_s=0.4*d_Ab, g_l=0.6*d_Ab, i.e. 60% of autoantibody production comes
from the long-lived pool at steady state. That split is tier-C (assumed), never measured.

THE REAL-WORLD CONSTRAINT I had not used. lessons.md records the signature dissociation
seen in CD19 CAR-T for lupus: *anti-dsDNA seroconverts while vaccine IgG holds*. Vaccine
IgG is the classic long-lived-plasma-cell product, so if AUTOantibody disappears while
LLPC-derived protective antibody persists, the autoreactive LLPC share must be LOW -- not
60%. Our parameterisation may be backwards for the autoreactive compartment specifically.

This sweeps the share and asks which values are compatible with observed seroconversion.
"""
import numpy as np
from model import make_params, simulate, observe

TG = np.linspace(0, 400, 120)


def ab_fall(gl_share, kill=0.25):
    """Autoantibody fall at day 240 under CD19 depletion (day 180-240)."""
    p = make_params(phi=0.0)
    p['g_s'] = (1.0 - gl_share) * p['d_Ab']
    p['g_l'] = gl_share * p['d_Ab']
    p['dep_kill'] = kill
    Yd = simulate(p, TG, rtol=1e-9, atol=1e-11)
    p0 = dict(p); p0['dep_kill'] = 0.0
    Yn = simulate(p0, TG, rtol=1e-9, atol=1e-11)
    i = int(np.argmin(np.abs(TG - 240)))
    ab_n, ab_d = Yn[i, 3], Yd[i, 3]
    # share of production from the protected pool at day 240, under depletion
    share = p['g_l'] * Yd[i, 2] / (p['g_s'] * Yd[i, 1] + p['g_l'] * Yd[i, 2])
    return 100 * (1 - ab_d / ab_n), 100 * share, 100 * (1 - Yd[i, 0] / Yn[i, 0])


print("Sensitivity of the headline to an ASSUMED parameter")
print("(seroconversion in CD19 CAR-T lupus implies a LARGE autoantibody fall)\n")
print(f"{'g_l share':>10} | {'Ab fall @d240':>14} | {'prod. from Pl':>14} | {'B-cell fall':>12}")
print("-" * 62)
for s in [0.0, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 0.95]:
    f, sh, bf = ab_fall(s)
    flag = "  <-- our assumed value" if abs(s - 0.60) < 1e-9 else ""
    print(f"{s:10.2f} | {f:13.1f}% | {sh:13.1f}% | {bf:11.1f}%{flag}")

print("\nThe 7.9% figure is NOT robust: it is monotone in an assumed parameter and")
print("spans essentially the whole range as that parameter moves. The DIRECTION")
print("(the protected pool's share rises after depletion) is structural; the")
print("MAGNITUDE is an assumption wearing a result's clothes.")
