"""R44 analysis. Seed SD reported SEPARATELY from binomial CI (lessons section 2)."""
import numpy as np, json
from scipy import stats
R=json.load(open('results/r44_squeeze.json'))
cvs=sorted({r['cv'] for r in R})
print(f"{'CV':>5} | {'control %':>14} | {'treated %':>14} | {'diff pp':>14} | "
      f"{'Cohen d':>13} | {'req N/arm':>11}")
print("-"*92)
agg={}
for cv in cvs:
    g=[r for r in R if r['cv']==cv]
    if not g: continue
    c=np.array([x['control'] for x in g])*100
    t=np.array([x['treated'] for x in g])*100
    d=np.array([x['diff_pp'] for x in g])
    cd=np.array([x['cohen_d'] for x in g])
    rn=np.array([x['req_n'] for x in g])
    agg[cv]=dict(c=c,t=t,d=d,cd=cd,rn=rn,n=len(g))
    print(f"{cv:5.2f} | {c.mean():6.1f} +-{c.std(ddof=1):5.2f} | {t.mean():6.1f} +-{t.std(ddof=1):5.2f} | "
          f"{d.mean():+6.1f} +-{d.std(ddof=1):5.2f} | {cd.mean():6.3f} +-{cd.std(ddof=1):5.3f} | "
          f"{np.median(rn):11.0f}   (+- = seed SD, {len(g)} seeds)")

lo,hi=min(cvs),max(cvs)
print("\n--- H1: binary squeeze (diff_pp falls) ---")
if lo in agg and hi in agg:
    a,b=agg[lo]['d'],agg[hi]['d']
    tt=stats.ttest_ind(a,b,equal_var=False)
    print(f"  diff_pp CV{lo:.2f} = {a.mean():+.2f} +-{a.std(ddof=1):.2f} vs "
          f"CV{hi:.2f} = {b.mean():+.2f} +-{b.std(ddof=1):.2f}")
    print(f"  Welch t = {tt.statistic:+.2f}, p = {tt.pvalue:.2e}  -> "
          f"{'H1 SUPPORTED' if (b.mean()<a.mean() and tt.pvalue<0.05) else 'H1 NOT SUPPORTED'}")

print("\n--- H2 (THE DISCRIMINATOR): does the CONTINUOUS effect size also fall? ---")
if lo in agg and hi in agg:
    a,b=agg[lo]['cd'],agg[hi]['cd']
    tt=stats.ttest_ind(a,b,equal_var=False)
    sl=stats.linregress([r['cv'] for r in R],[r['cohen_d'] for r in R])
    print(f"  Cohen d CV{lo:.2f} = {a.mean():.3f} +-{a.std(ddof=1):.3f} vs "
          f"CV{hi:.2f} = {b.mean():.3f} +-{b.std(ddof=1):.3f}")
    print(f"  Welch t = {tt.statistic:+.2f}, p = {tt.pvalue:.2e}")
    print(f"  slope of d vs CV = {sl.slope:+.4f} (p={sl.pvalue:.2e}, r2={sl.rvalue**2:.3f})")
    if b.mean()<a.mean() and tt.pvalue<0.05:
        print("  -> H2 SUPPORTED: erosion is in the DISTRIBUTIONS, not the threshold.")
        print("     R43 section 4 SURVIVES.")
    else:
        print("  -> H2 FALSIFIED: continuous effect size does not fall.")
        print("     The squeeze is a BINARY-THRESHOLD ARTIFACT. R43 section 4 RETRACTED.")

print("\n--- H3: required sample size per arm rises with CV ---")
med=[(cv,np.median(agg[cv]['rn'])) for cv in cvs if cv in agg]
for cv,m in med: print(f"  CV {cv:.2f}: median required N/arm = {m:.0f}")
if len(med)>1:
    mono=all(med[i+1][1]>=med[i][1] for i in range(len(med)-1))
    print(f"  monotonic increase: {'YES -> H3 supported' if mono else 'NO -> H3 not supported'}")
    print(f"  fold change {med[0][0]:.2f}->{med[-1][0]:.2f}: {med[-1][1]/med[0][1]:.2f}x")
