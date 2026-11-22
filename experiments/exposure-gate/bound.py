"""Q15 Part B: binormal bound on incremental AUC of an exposure marker over month-3 proteinuria (PREREG.md)."""
import numpy as np, json
from scipy.stats import norm
d_of = lambda A: np.sqrt(2) * norm.ppf(A)
A_of = lambda d: norm.cdf(d / np.sqrt(2))
def dauc(A1, A2, rho):
    d1, d2 = d_of(A1), d_of(A2)
    d12 = np.sqrt((d1**2 + d2**2 - 2*rho*d1*d2) / (1 - rho**2))
    return A_of(d12) - A1
def pooled_from_iqr(m1, q1a, q1b, n1, m0, q0a, q0b, n0, log=False):
    f = np.log if log else (lambda x: x)
    s1 = (f(q1b) - f(q1a)) / 1.35; s0 = (f(q0b) - f(q0a)) / 1.35
    sp = np.sqrt(((n1-1)*s1**2 + (n0-1)*s0**2) / (n1+n0-2))
    return A_of(abs(f(m1) - f(m0)) / sp)
# month-3 24h urine protein, Ren Fail 2026 (PMC13757897), n=317 with 3-month data: CRR 0.54 (0.19-1.31) n236, non-CRR 2.75 (1.30-4.97) n81
U = {"raw_IQR (primary)": pooled_from_iqr(0.54, .19, 1.31, 236, 2.75, 1.30, 4.97, 81),
     "log_IQR": pooled_from_iqr(0.54, .19, 1.31, 236, 2.75, 1.30, 4.97, 81, log=True),
     "Youden_2x2 (1.3 g)": ((173/236) + (62/81)) / 2,
     "AMP month-3 proteinuria (round 93)": 0.80}
# exposure markers
sp = np.sqrt((116*21.7**2 + 66*16.2**2) / 182)
E = {"MPA-AUC ROC, LSM 2026 paediatric n65 (primary)": 0.70,
     "MPA-AUC meta-analysis 16 studies, mean/SD (secondary)": A_of((51.4-30.3)/sp),
     "LUNAR complete peripheral B-cell depletion 2x2 (wk78 window)": ((25/27) + (13/41)) / 2}
rhos = [-0.2, 0, 0.2, 0.4, 0.6]
out = {"UPCR_AUC": U, "exposure_AUC": E, "grid": {}}
print("month-3 UPCR AUCs:", {k: round(v, 3) for k, v in U.items()})
print("exposure AUCs:", {k: round(v, 3) for k, v in E.items()})
for uk, A1 in U.items():
    for ek, A2 in E.items():
        row = {r: round(float(dauc(A1, A2, r)), 3) for r in rhos}
        out["grid"][f"{uk} | {ek}"] = row
        print(f"{uk:36s} | {ek[:52]:52s} | " + "  ".join(f"rho{r:+.1f}:{v:+.3f}" for r, v in row.items()))
# verdict per PREREG, on primary UPCR
A1 = U["raw_IQR (primary)"]; verdict = {}
for ek, A2 in E.items():
    g = {r: dauc(A1, A2, r) for r in np.linspace(-0.2, 0.6, 81)}
    if all(v >= 0.05 for r, v in g.items() if 0 <= r <= 0.4 + 1e-9): v = "B-PASS"
    elif any(v >= 0.05 for r, v in g.items() if r <= 0.4 + 1e-9): v = "B-CONTINGENT"
    else: v = "B-FAIL"
    rmax = max([r for r, x in g.items() if x >= 0.05], default=None)
    verdict[ek] = (v, None if rmax is None else round(float(rmax), 2))
print("VERDICT (primary UPCR):", verdict)
# what standalone exposure AUC would be needed for dAUC>=0.05 at rho=0 and rho=0.3
for r in (0, 0.3):
    need = next(a for a in np.arange(0.5, 0.99, 0.001) if dauc(A1, a, r) >= 0.05)
    print(f"exposure AUC needed at rho={r}: {need:.3f}")
out["verdict"] = verdict
json.dump(out, open(__file__.replace("bound.py", "result.json"), "w"), indent=1, default=float)
