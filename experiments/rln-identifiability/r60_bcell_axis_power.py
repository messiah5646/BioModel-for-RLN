"""
R60: can the Stage 1 gate be tested on the B-cell axis with AGGREGATE trial data?

Stage 1's gate (ROADMAP.md) requires a model of the B/Ps/Pl/Ab axis to beat persistence
and linear extrapolation OUT OF SAMPLE on levels. Round 58 harvested the only multi-visit
B-cell series in the EU CTR: NOBILITY's CD19+ percent-change, 5 usable randomised visits.

This is NOT a model fit. It asks the prior question: is there enough signal in this data
for ANY model to be distinguishable from persistence? If the per-visit standard error
exceeds the total excursion of the series, the gate is untestable here regardless of what
model is proposed -- an identifiability/power statement, not a modelling one.

PRE-REGISTRATION, written before computing (contamination declared: the CSV was read in
round 58 and its values are visible in the header traps, so this is not blind):

  H1. The treated arm's total excursion across weeks 2-52 is SMALLER than the median
      per-visit standard error. Predicted: TRUE -> the arm is a saturated floor and
      persistence is near-exact by construction.
  H2. The placebo arm's signal-to-noise ratio (excursion / median SEM) is < 2.
      Predicted: TRUE.
  H3. Therefore no model can demonstrate an advantage over persistence on this data at
      conventional thresholds. Operationalised: the RMSE improvement a model would need
      to clear 1 SEM of the held-out point exceeds the persistence RMSE itself.

  Falsifier: if either arm's excursion exceeds 2x its median SEM, the data carries usable
  signal and H3 fails -- a real Stage 1 test becomes possible and should be run.

Usable visits: weeks 2, 4, 12, 24, 52 only. Weeks 104+ are a separate "Bcfu" follow-up
cohort with n=0-5, non-monotone n, and month-12 values byte-identical to week 104.
The baseline row is in absolute cells/uL while all others are percent change -- it is
EXCLUDED, never differenced against (round 58 TRAP 1).
"""
import csv
import math

VISITS = {'week_2': 2, 'week_4': 4, 'week_12': 12, 'week_24': 24, 'week_52': 52}
PATH = '../../data/real/obinutuzumab_nobility_per_visit.csv'


def load():
    with open(PATH) as f:
        rows = list(csv.DictReader([l for l in f if not l.lstrip().startswith('#')]))
    out = {}
    for r in rows:
        if r['visit'] not in VISITS:
            continue                       # excludes baseline (different unit) and Bcfu
        if r['unit_note'] != 'percent_change':
            continue
        if r['bcell_value'] == 'not reported' or r['bcell_sd'] == 'not reported':
            continue
        out.setdefault(r['arm'], []).append(
            (VISITS[r['visit']], float(r['bcell_value']),
             float(r['bcell_sd']), int(r['bcell_n'])))
    for a in out:
        out[a].sort()
    return out


D = load()
print("R60 -- is the Stage 1 gate testable on NOBILITY's B-cell series?")
print("Usable randomised visits only (weeks 2-52); baseline and Bcfu excluded.\n")

summary = {}
for arm, rows in D.items():
    sems = [sd / math.sqrt(n) for _, _, sd, n in rows]
    vals = [v for _, v, _, _ in rows]
    exc = max(vals) - min(vals)
    med_sem = sorted(sems)[len(sems) // 2]
    snr = exc / med_sem
    summary[arm] = dict(excursion=exc, med_sem=med_sem, snr=snr, vals=vals, sems=sems)
    print(f"{arm}")
    print(f"  {'week':>6} | {'value':>9} | {'SD':>9} | {'n':>3} | {'SEM':>8}")
    for (w, v, sd, n), s in zip(rows, sems):
        print(f"  {w:>6} | {v:9.3f} | {sd:9.3f} | {n:>3} | {s:8.3f}")
    print(f"  total excursion {exc:.2f} pp | median SEM {med_sem:.2f} pp"
          f" | SNR {snr:.2f}\n")

# ---- the pre-declared trivial comparator, scored on the one held-out point ----
print("Held-out test: fit on weeks 2-24, predict week 52 (1 point per arm).")
print(f"{'arm':>22} | {'persistence':>11} | {'linear':>9} | {'observed':>9} | {'1 SEM':>7}")
print("-" * 70)
for arm, rows in D.items():
    fitw = [r for r in rows if r[0] <= 24]
    held = [r for r in rows if r[0] == 52][0]
    pers = fitw[-1][1]
    (w1, v1, _, _), (w2, v2, _, _) = fitw[-2], fitw[-1]
    lin = v2 + (v2 - v1) / (w2 - w1) * (52 - w2)
    sem52 = held[2] / math.sqrt(held[3])
    summary[arm]['pers_err'] = abs(pers - held[1])
    summary[arm]['lin_err'] = abs(lin - held[1])
    summary[arm]['sem52'] = sem52
    print(f"{arm:>22} | {pers:11.3f} | {lin:9.3f} | {held[1]:9.3f} | {sem52:7.3f}")

print()
for arm in D:
    s = summary[arm]
    print(f"{arm}: |persistence error| {s['pers_err']:.3f} pp, "
          f"|linear error| {s['lin_err']:.3f} pp, 1 SEM at wk52 {s['sem52']:.3f} pp")
    print(f"   persistence error is {s['pers_err']/s['sem52']:.2f} SEM "
          f"-> {'WITHIN' if s['pers_err'] < s['sem52'] else 'OUTSIDE'} noise")

print("\n--- verdicts ---")
t = summary['obinutuzumab_1000mg']
p = summary['placebo']
print(f"H1 (treated excursion < median SEM): "
      f"{'SUPPORTED' if t['excursion'] < t['med_sem'] else 'FALSIFIED'} "
      f"({t['excursion']:.2f} vs {t['med_sem']:.2f})")
print(f"H2 (placebo SNR < 2): "
      f"{'SUPPORTED' if p['snr'] < 2 else 'FALSIFIED'} (SNR {p['snr']:.2f})")
gate = all(summary[a]['pers_err'] < summary[a]['sem52'] for a in D)
print(f"H3 (persistence already within 1 SEM in BOTH arms -> no room for a model to "
      f"show an advantage): {'SUPPORTED' if gate else 'FALSIFIED'}")
