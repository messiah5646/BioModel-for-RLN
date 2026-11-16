"""Round 70 pilot: local structural identifiability of the damage equation ALONE (BUILD-PLAN track A step 2).
dD/dt = (k_dam*A + k_tox*p0*A*N*leak + k_prog)*N, N = 1-D; output y = D.
Method: rank of Jacobian of (y', y'', ...) wrt parameters at random rational points (observability-rank test).
Case 1: A constant, known. Case 2: A(t) known, dA/dt = -lam*A (lam known). Case 3: A constant but UNKNOWN."""
import sympy as sp, random
kd,kt,p0,lk,kp,A,D,lam=sp.symbols('k_dam k_tox p0 leak k_prog A D lam',positive=True)
N=1-D; f=(kd*A+kt*p0*A*N*lk+kp)*N
def rank(dA,pars,n=7):
    ys=[f]
    for _ in range(n-1): ys.append(sp.diff(ys[-1],D)*f+sp.diff(ys[-1],A)*dA)
    J=sp.Matrix(ys).jacobian(pars); rs=[]
    for s in range(3):
        random.seed(s); sub={v:sp.Rational(random.randint(2,40),random.randint(41,90)) for v in (kd,kt,p0,lk,kp,A,D,lam)}
        rs.append(J.subs(sub).rank())
    return rs
P=[kd,kt,p0,lk,kp]
print("case1 A const known, 5 params: rank",rank(0,P))
print("case2 A decaying known, 5 params: rank",rank(-lam*A,P))
print("case3 A const unknown, 6 unknowns: rank",rank(0,P+[A]))
c=sp.symbols('c',positive=True); f=(kd*A+c*A*N+kp)*N
print("reparam c=k_tox*p0*leak; case1 [k_dam,c,k_prog]:",rank(0,[kd,c,kp]),"; case2:",rank(-lam*A,[kd,c,kp]))
