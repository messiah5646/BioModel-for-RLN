"""LUMEN patient-state layer: blood IFN as a therapy-lowered patient set point (round 96).

IFN_i(t) = S_i * f(t)
  S_i  : the patient's own baseline score (patient state; round 89: ~54% of variance is patient identity)
  f(t) : therapy multiplier after deep B-cell depletion, f(0) = 1.
         Fitted values (round 96, median post/base, monocyte 6-gene score):
           CD3xCD19 day 28: 0.98 (no fall yet), CD3xCD19 month 3: 0.75, CD19 CAR-T post: 0.69.
         One shared late multiplier (F_LATE) serves both mechanisms (H3 passed).
  Shape: no change until LAG_DAYS, then exponential approach to F_LATE. LAG and HALF_LIFE are bounded, not
  identified: only D28 and M3 exist, and CAR-T post has no recorded day.
Not a flare variable (round 89). Research use only.
"""
import numpy as np
F_LATE = 0.73        # pooled median of the 18 late post/base ratios (0.732); 0.69 CAR-T, 0.75 CD3xCD19
LAG_DAYS = 28        # f(28) ~ 1 observed
HALF_LIFE_DAYS = 14  # placeholder within the range consistent with f(90) ~ 0.75 (any half-life <= ~20 d fits)

def f(t_days, f_late=F_LATE, lag=LAG_DAYS, half_life=HALF_LIFE_DAYS, depleted=True):
    t = np.asarray(t_days, float)
    if not depleted:
        return np.ones_like(t)
    k = np.log(2) / half_life
    return np.where(t <= lag, 1.0, f_late + (1 - f_late) * np.exp(-k * (t - lag)))

def ifn(t_days, baseline, **kw):
    return baseline * f(t_days, **kw)

if __name__ == '__main__':
    print({t: round(float(f(t)), 3) for t in [0, 28, 60, 90, 180]})
