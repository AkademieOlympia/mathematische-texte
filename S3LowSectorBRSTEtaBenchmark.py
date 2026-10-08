# -*- coding: utf-8 -*-
"""
S3LowSectorBRSTEtaBenchmark.py

Analytic S^3 benchmark for the low G2 sectors V_14 and V_27 of the
Freed-Moore / BRST combination

    eta_R - 2 eta_FP - eta_NK.

This is a MATHEMATICAL BENCHMARK ONLY.
It is not a physical M-theory prediction.

External spectrum on a round S^3 of radius L:
    lambda_{n,±} = ±(n+3/2)/L
    multiplicity m_n = (n+1)(n+2).

For one internal 2x2 block
    [[s q, d/ell],
     [d/ell, 0]]
define
    rho = L/ell,
    y   = q L,
    C   = s y / 2,
    Delta = d rho.

The corresponding 11D dimensionless eigenvalues over an external
eigenvalue epsilon*a_n, a_n=n+3/2, are
    C ± sqrt((epsilon*a_n + C)^2 + Delta^2).

Away from zero crossings, zeta/eta continuation gives the local
polynomial
    eta_local(0) =
        C * (1 + 2 Delta^2 - (16/3) C^2)
      = (s y/2) * (1 + 2 d^2 rho^2 - (4/3)s^2 y^2).

Finite sign-crossing corrections are then added exactly:
    +2 m_n sign(C)
for each (n, epsilon) for which both eigenvalues have the same sign,
equivalently
    C^2 > (epsilon*a_n + C)^2 + Delta^2.

Low-sector input:
V14:
  R  : (s,d^2)=(-15,8), (9,8)
  FP : (5/3,12)
  NK : (-15,12)

V27:
  R  : (-15,28/3), (1,28/3), (9,28/3)
  FP : (-1,40/3)
  NK : (9,40/3)

Benchmark point:
    rho = 1/2
    x = q ell = sqrt(2/5)
    y = q L = x/2
This inherits the magnitude L=ell/2 only as a same-scale S^3 stress test.
"""

from __future__ import annotations
import sympy as sp

Q = sp.Rational
rho, y = sp.symbols("rho y", positive=True, real=True)

def eta_local(s, d2):
    s = sp.Rational(s)
    d2 = sp.Rational(d2)
    C = s*y/2
    Delta2 = d2*rho**2
    return sp.simplify(C * (1 + 2*Delta2 - sp.Rational(16,3)*C**2))

def crossing_correction(s, d2, rho_val, y_val, nmax=100):
    s = sp.Rational(s)
    d2 = sp.Rational(d2)
    C = sp.simplify(s*y_val/2)
    Delta2 = sp.simplify(d2*rho_val**2)
    corr = sp.Integer(0)
    events = []
    for n in range(nmax+1):
        a = sp.Rational(2*n+3, 2)
        mult = (n+1)*(n+2)
        for eps in (+1, -1):
            R2 = sp.simplify((eps*a + C)**2 + Delta2)
            diff = sp.N(C**2 - R2, 50)
            if diff > 0:
                step = 2*mult*(1 if C > 0 else -1)
                corr += step
                events.append({
                    "n": n,
                    "epsilon": eps,
                    "multiplicity": mult,
                    "correction": int(step),
                })
    return sp.Integer(corr), events

R14 = [(-15,8),(9,8)]
FP14 = (Q(5,3),12)
NK14 = (-15,12)

R27 = [(-15,Q(28,3)),(1,Q(28,3)),(9,Q(28,3))]
FP27 = (-1,Q(40,3))
NK27 = (9,Q(40,3))

B14_local = sp.factor(sum(eta_local(*b) for b in R14)
                      - 2*eta_local(*FP14) - eta_local(*NK14))
B27_local = sp.factor(sum(eta_local(*b) for b in R27)
                      - 2*eta_local(*FP27) - eta_local(*NK27))

x = sp.sqrt(Q(2,5))
rho0 = Q(1,2)
y0 = x/2

def sector_correction(R, FP, NK):
    corr = 0
    details = {}
    for i,b in enumerate(R):
        c,e = crossing_correction(*b, rho0, y0)
        corr += c
        details[f"R_{i}"] = {"corr": int(c), "events": e}
    c,e = crossing_correction(*FP, rho0, y0)
    corr -= 2*c
    details["FP"] = {"raw_corr": int(c), "weighted_corr": int(-2*c), "events": e}
    c,e = crossing_correction(*NK, rho0, y0)
    corr -= c
    details["NK"] = {"raw_corr": int(c), "weighted_corr": int(-c), "events": e}
    return sp.Integer(corr), details

corr14, details14 = sector_correction(R14, FP14, NK14)
corr27, details27 = sector_correction(R27, FP27, NK27)

B14_bench = sp.simplify(B14_local.subs({rho:rho0, y:y0}) + corr14)
B27_bench = sp.simplify(B27_local.subs({rho:rho0, y:y0}) + corr27)

eta_14 = sp.simplify(14*B14_bench)
eta_27 = sp.simplify(27*B27_bench)
eta_total = sp.simplify(eta_14 + eta_27)

print("B14_local(rho,y) =", B14_local)
print("B27_local(rho,y) =", B27_local)
print()
print("benchmark rho=1/2, x=q ell=sqrt(2/5), y=x/2")
print("14 crossing correction =", corr14)
print("27 crossing correction =", corr27)
print()
print("B14 benchmark =", B14_bench, "≈", sp.N(B14_bench, 16))
print("B27 benchmark =", B27_bench, "≈", sp.N(B27_bench, 16))
print("dim-weighted eta_14 =", eta_14, "≈", sp.N(eta_14, 16))
print("dim-weighted eta_27 =", eta_27, "≈", sp.N(eta_27, 16))
print("LOW-SECTOR BRST ETA =", eta_total, "≈", sp.N(eta_total, 16))
print("phase (pi/4)*eta ≈", sp.N(sp.pi*eta_total/4, 16))
print()
print("This is not the full G2 KK tower and not a physical M-theory eta value.")
