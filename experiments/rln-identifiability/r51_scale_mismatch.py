"""
R51: can molecular-scale models inform the quantity our model gets WRONG?

Round 50's defect is specific: the model understates the observed week 24->48 eGFR
decline ~4x and loses to linear extrapolation at short horizon. The question is whether
anything a molecular-scale model (structure, affinity, sequence) could ground is even
CONNECTED to that quantity.

Round 40's structural zero gives the machinery: if a parameter does not move a quantity,
no amount of precision about that parameter helps. Apply it in the other direction --
rank every parameter by its influence on the FAILING quantity, grouped by the scale a
molecular model could plausibly inform.
"""
import numpy as np
from model import make_params, simulate, observe

TG = np.array([0.0, 24 * 7.0, 48 * 7.0])

GROUPS = {
    'humoral / molecular (what a protein or sequence model could touch)':
        ['s_B', 'd_B', 'f_s', 'd_Ps', 'f_l', 'd_Pl', 'g_s', 'g_l', 'd_Ab', 'k_a', 'dep_kill'],
    'kidney block (population rate constants)':
        ['k_dam', 'k_tox', 'k_prog', 'c_inf'],
    'observation / pharmacology':
        ['p0', 'G0', 's_inf', 'tau_s', 'leak_red', 'h_drug'],
}


def base():
    p = make_params(phi=0.0)
    p['p0'], p['G0'] = 4.14, 78.3 / (1 - p['c_inf'])
    p['k_dam'], p['k_prog'], p['k_tox'] = 0.0, 1.03e-4, 0.2 * 1.03e-4   # round 49 V2
    p['s_inf'] = 0.536
    p['drug_on'], p['drug_off'] = 1.0, 1e9
    return p


def late_trend(p):
    """the quantity round 50 gets wrong: week 24 -> 48 change in eGFR."""
    Y = simulate(p, TG, rtol=1e-10, atol=1e-12)
    O = np.array([observe(y, p, t) for t, y in zip(TG, Y)])
    return O[2, 1] - O[1, 1]


p0 = base()
ref = late_trend(p0)
print(f"FAILING QUANTITY: week 24->48 eGFR change")
print(f"  model {ref:+.3f}   observed (treated) -3.30, (placebo) -1.70   -> understated ~4x\n")
print("Relative sensitivity  |d log(trend) / d log(param)|, grouped by scale")
print("(a parameter with ~0 cannot inform this quantity at ANY precision -- round 40)\n")
rows = []
for gname, params in GROUPS.items():
    print(f"  {gname}")
    gmax = 0.0
    for nm in params:
        if nm not in p0 or p0[nm] == 0:
            print(f"    {nm:12s}      (zero/absent in this variant)")
            continue
        h = abs(p0[nm]) * 1e-3
        pp = dict(p0); pp[nm] += h
        pm = dict(p0); pm[nm] -= h
        s = (late_trend(pp) - late_trend(pm)) / (2 * h) * p0[nm] / abs(ref)
        gmax = max(gmax, abs(s))
        print(f"    {nm:12s} {abs(s):10.3e}")
        rows.append((gname, nm, abs(s)))
    print(f"    -> group max {gmax:.3e}\n")

hum = max(s for g, n, s in rows if g.startswith('humoral'))
kid = max(s for g, n, s in rows if g.startswith('kidney'))
obs = max(s for g, n, s in rows if g.startswith('observation'))
print(f"GROUP MAXIMA   humoral {hum:.3e} | kidney {kid:.3e} | observation {obs:.3e}")
print(f"kidney / humoral ratio = {kid/max(hum,1e-300):.1f}x")
print()
print("VERDICT on grounding the failing quantity with a molecular-scale model:")
if kid > 10 * hum:
    print("  The failing quantity is dominated by kidney-block rate constants.")
    print("  Humoral parameters are weaker by >10x, so perfect molecular-scale knowledge")
    print("  cannot fix round 50's defect. This is round 40's structural result applied")
    print("  in reverse: precision about a parameter the quantity barely depends on buys")
    print("  nothing.")
else:
    print("  Humoral parameters are NOT negligible for this quantity -- molecular grounding")
    print("  could in principle help. The scale-mismatch argument fails.")


def why_and_reachable():
    """Why the structural-zero argument did not transfer, and is the gap reachable?"""
    print("\n" + "="*70)
    print("WHY THE SCALE-MISMATCH ARGUMENT FAILED -- round 40's zero is DIRECTIONAL")
    print("  Round 40: kidney PARAMETERS do not move molecular STATES (feed-forward).")
    print("  That is the upstream-observable direction. The reverse is not symmetric:")
    print("  in V2 activity is driven by autoantibody (dA/dt = k_a*Ab - k_clear*A), and")
    print("  the damage driver depends on activity. So humoral PARAMETERS strongly")
    print("  influence downstream quantities even though humoral OBSERVABLES carry no")
    print("  information about the downstream damage split. I conflated the directions.")
    print()
    print("IS THE ROUND-50 GAP REACHABLE BY THE TOP HUMORAL PARAMETER?")
    print("  d_Ab (autoantibody clearance) is the most influential humoral parameter.")
    print(f"  {'d_Ab x':>8} | {'wk24->48 eGFR':>14} | reaches observed -3.30?")
    print("  " + "-"*52)
    p = base(); ref_dab = p['d_Ab']
    hit = None
    for mult in [0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]:
        q = dict(p); q['d_Ab'] = ref_dab * mult
        t = late_trend(q)
        if t <= -3.30 and hit is None:
            hit = mult
        print(f"  {mult:8.2f} | {t:+14.3f} | {'YES' if t <= -3.30 else 'no'}")
    print()
    if hit:
        print(f"  Reachable at d_Ab x{hit:.2f} -- so the defect CAN be closed by this")
        print(f"  parameter alone, and grounding it is worth real effort.")
    else:
        print("  NOT reachable across a 64x range of d_Ab. So although the failing")
        print("  quantity is SENSITIVE to d_Ab, no physiologically plausible value closes")
        print("  the gap -- the defect is structural, not a parameter-grounding problem.")


if __name__ == '__main__' and '--why' in __import__('sys').argv:
    why_and_reachable()
