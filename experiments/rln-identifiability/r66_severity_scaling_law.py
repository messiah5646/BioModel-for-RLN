"""
R66: attacking the project's own "severity is worth 437x" headline.

TWO CHARGES against the claim as written in where-this-stands.md, SYSTEM-DESIGN.md and
positioning.md:

  CHARGE 1 -- CONFOUNDED. The 437x is sigma(sev=1, phi=0) / sigma(sev=20, phi=20).
  It varies severity AND feedback at once, then is reported as a severity effect and
  contrasted against the observable set's 2.07x. The severity-only number is different.

  CHARGE 2 -- POSSIBLY TRIVIAL. k_dam multiplies A in dD/dt. At low severity D stays
  near zero, so k_dam barely moves any observable and is necessarily ill-determined.
  "A parameter is better identified when its effect is larger" is arithmetic, not a
  finding. If sigma follows a clean power law in D_final, the 437x is a restatement of
  the chosen severity range and carries no extra information.

PRE-REGISTRATION, and the contamination is declared because it is real:
  I read the 9-row table printed by r40_feedback.py BEFORE writing this script, and
  estimated the exponent by hand from the three phi=0 rows. So the exponent value is
  IN-SAMPLE and is not a prediction. What IS pre-registered and out-of-sample:

  H1. The severity-only ratio (fixed phi) is substantially smaller than 437x.
  H2. log sigma is linear in log D_final with exponent b ~= -2, r^2 > 0.99, on the
      three training severities.
  H3. THE REAL TEST -- a law fitted on sev in {1, 10, 20} predicts sigma at
      sev in {2, 4, 6, 14, 30}, never used in fitting, to within 25% median
      absolute relative error. Falsified at >= 25%.

  If H3 holds, the headline must be restated as a scaling law and the 437x retired as
  a free-floating number. If H3 fails, the power law is a coincidence of three points
  and the original claim stands as-is.
"""
import json
import numpy as np
from model import make_params
from r40_feedback import marginal_sigma, CHANNELS
from model import simulate, VISIT_GRID

TRAIN_SEV = [1, 10, 20]
TEST_SEV = [2, 4, 6, 14, 30]
PHIS = [0.0, 5.0, 20.0]


def run(sev, phi, obs='S0'):
    p = make_params(phi)
    p['k_dam'] *= sev
    p['k_tox'] *= sev
    p['dep_kill'] = 0.25          # as in r40_feedback.main(): this IS a depletion study
    Y = simulate(p, VISIT_GRID)
    D_final = float(Y[-1, 5])
    sig = marginal_sigma(p, CHANNELS[obs])[0]
    return D_final, float(sig)


print("R66 -- is 'severity is worth 437x' confounded, trivial, or a scaling law?\n")

# ---------------- CHARGE 1: decompose the 437x ----------------
print("CHARGE 1 -- decomposition of the headline ratio (obs = S0, clinical only)")
grid = {}
for sev in TRAIN_SEV:
    for phi in PHIS:
        grid[(sev, phi)] = run(sev, phi)

headline = grid[(1, 0.0)][1] / grid[(20, 20.0)][1]
sev_only_phi0 = grid[(1, 0.0)][1] / grid[(20, 0.0)][1]
sev_only_phi20 = grid[(1, 20.0)][1] / grid[(20, 20.0)][1]
fb_only_sev1 = grid[(1, 0.0)][1] / grid[(1, 20.0)][1]
fb_only_sev20 = grid[(20, 0.0)][1] / grid[(20, 20.0)][1]

print(f"  headline as quoted   sigma(sev1,phi0)/sigma(sev20,phi20) = {headline:8.1f}x"
      f"   <- varies BOTH")
print(f"  severity alone, phi=0                                    = {sev_only_phi0:8.1f}x")
print(f"  severity alone, phi=20                                   = {sev_only_phi20:8.1f}x")
print(f"  feedback alone, sev=1                                    = {fb_only_sev1:8.1f}x")
print(f"  feedback alone, sev=20                                   = {fb_only_sev20:8.1f}x")
print(f"  observable set, best gain anywhere in the grid            =     2.07x (R41)")
print(f"\n  H1 (severity-only << 437x): "
      f"{'SUPPORTED' if sev_only_phi0 < 0.5 * headline else 'FALSIFIED'}"
      f"  -- the quoted figure overstates the severity effect by "
      f"{headline / sev_only_phi0:.1f}x")

# ---------------- CHARGE 2: power law, fitted on training severities ----------------
print("\nCHARGE 2 -- does sigma follow a power law in D_final? (fit on sev 1/10/20)")
fits = {}
for phi in PHIS:
    D = np.array([grid[(s, phi)][0] for s in TRAIN_SEV])
    S = np.array([grid[(s, phi)][1] for s in TRAIN_SEV])
    b, a = np.polyfit(np.log(D), np.log(S), 1)
    pred = np.exp(a + b * np.log(D))
    ss_res = np.sum((np.log(S) - np.log(pred)) ** 2)
    ss_tot = np.sum((np.log(S) - np.log(S).mean()) ** 2)
    r2 = 1 - ss_res / ss_tot
    fits[phi] = (a, b, r2)
    print(f"  phi={phi:5.1f}:  sigma = {np.exp(a):7.4f} * D_final^({b:+.3f})   r^2 = {r2:.5f}")

b0, r20 = fits[0.0][1], fits[0.0][2]
print(f"\n  H2 (exponent ~= -2 and r^2 > 0.99 at phi=0): "
      f"{'SUPPORTED' if abs(b0 + 2.0) < 0.25 and r20 > 0.99 else 'FALSIFIED'}"
      f"  (b = {b0:+.3f}, r^2 = {r20:.5f})")

# ---------------- H3: the out-of-sample test ----------------
print("\nH3 -- THE REAL TEST: predict sigma at severities never used in the fit")
print(f"{'phi':>5} {'sev':>4} {'D_final':>8} {'sigma obs':>11} {'sigma pred':>11} {'rel err':>9}")
print("-" * 54)
rows = []
for phi in PHIS:
    a, b, _ = fits[phi]
    for sev in TEST_SEV:
        Df, sg = run(sev, phi)
        pr = float(np.exp(a + b * np.log(Df)))
        rel = abs(pr - sg) / sg
        rows.append(dict(phi=phi, sev=sev, D_final=Df, sigma=sg, pred=pr, rel_err=rel))
        print(f"{phi:5.1f} {sev:4d} {Df:8.4f} {sg:11.4f} {pr:11.4f} {rel:8.1%}")

errs = np.array([r['rel_err'] for r in rows])
med = float(np.median(errs))
errs0 = np.array([r['rel_err'] for r in rows if r['phi'] == 0.0])
print(f"\n  median |rel err| across all {len(errs)} held-out points: {med:.1%}")
print(f"  median |rel err| at phi=0 only:                     {np.median(errs0):.1%}")
print(f"  max    |rel err|:                                   {errs.max():.1%}")
print(f"\n  H3 (median < 25%): {'SUPPORTED' if med < 0.25 else 'FALSIFIED'}")

# ---------------- what the law implies, stated as a check not a claim ----------------
print("\nIMPLICATION, if the law holds: the 437x is not a property of the disease,")
print("it is a property of the severity range chosen. Predicted ratio from the law:")
d1, d20 = grid[(1, 0.0)][0], grid[(20, 0.0)][0]
print(f"  (D_final[sev20]/D_final[sev1])^(-b) = ({d20/d1:.2f})^{-b0:.3f}"
      f" = {(d20/d1) ** (-b0):.1f}x   vs measured {sev_only_phi0:.1f}x")

json.dump(dict(headline=headline, sev_only_phi0=sev_only_phi0,
               sev_only_phi20=sev_only_phi20, fb_only_sev1=fb_only_sev1,
               fb_only_sev20=fb_only_sev20,
               fits={str(k): dict(const=float(np.exp(v[0])), exponent=float(v[1]),
                                  r2=float(v[2])) for k, v in fits.items()},
               heldout=rows, median_rel_err=med,
               H1=bool(sev_only_phi0 < 0.5 * headline),
               H2=bool(abs(b0 + 2.0) < 0.25 and r20 > 0.99),
               H3=bool(med < 0.25)),
          open('results/r66_severity_scaling.json', 'w'), indent=1)
print("\nsaved -> results/r66_severity_scaling.json")
