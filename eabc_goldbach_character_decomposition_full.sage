# SageMath 10.5
# EABC-Goldbach-Bias: exakte Charakterzerlegung bis NMAX
# Basis: 1, chi_3, chi_4, chi_12=chi_3*chi_4
#
# E=1 mod 12, A=5 mod 12, B=7 mod 12, C=11 mod 12.
# Für Primzahlen p>3:
# 1_E=(1+chi3)(1+chi4)/4
# 1_A=(1-chi3)(1+chi4)/4
# 1_B=(1+chi3)(1-chi4)/4
# 1_C=(1-chi3)(1-chi4)/4

from sage.all import *
import numpy as np
import csv

NMAX = 10^6
NMAX = int(NMAX)

def chi3(n):
    r = n % 3
    if r == 1: return 1
    if r == 2: return -1
    return 0

def chi4(n):
    r = n % 4
    if r == 1: return 1
    if r == 3: return -1
    return 0

# FFT length
L = 1
while L <= 2*NMAX:
    L <<= 1

basis_names = ["1","chi3","chi4","chi12"]
basis = {k: np.zeros(L, dtype=np.float64) for k in basis_names}

for p in prime_range(5, NMAX+1):
    p = int(p)
    c3 = chi3(p)
    c4 = chi4(p)
    basis["1"][p] = 1
    basis["chi3"][p] = c3
    basis["chi4"][p] = c4
    basis["chi12"][p] = c3*c4

F = {k: np.fft.rfft(v) for k,v in basis.items()}

# C_ab(N)=sum_{p+q=N} a(p)b(q), ordered
C = {}
for a in basis_names:
    for b in basis_names:
        C[(a,b)] = np.rint(
            np.fft.irfft(F[a]*F[b], n=L)[:NMAX+1]
        ).astype(np.int64)

sig = {
    "E":{"1":1,"chi3":1,"chi4":1,"chi12":1},
    "A":{"1":1,"chi3":-1,"chi4":1,"chi12":-1},
    "B":{"1":1,"chi3":1,"chi4":-1,"chi12":-1},
    "C":{"1":1,"chi3":-1,"chi4":-1,"chi12":1},
}

def Rord(X,Y):
    ans = np.zeros(NMAX+1, dtype=np.float64)
    for a in basis_names:
        for b in basis_names:
            ans += sig[X][a]*sig[Y][b]*C[(a,b)]
    return ans/16.0

comparisons = {
    0:(("E","C"),("A","B")),
    2:(("E","E"),("B","B")),
    6:(("E","A"),("B","C")),
    10:(("A","A"),("C","C")),
}

print("=== EABC Goldbach character audit ===")
print("NMAX =", NMAX)
print()

for r,(pair1,pair2) in comparisons.items():
    X,Y = pair1
    U,V = pair2
    start = r if r >= 4 else r+12
    Ns = np.arange(start, NMAX+1, 12, dtype=np.int64)
    Ns = Ns[(Ns>=4) & ((Ns % 2)==0)]

    R1 = Rord(X,Y)[Ns]
    R2 = Rord(U,V)[Ns]
    D = R1-R2
    T = R1+R2
    mask = T>0

    print("N mod 12 =", r, ":", X+Y, "-", U+V)
    print(" mean D =", float(np.mean(D)))
    print(" global relative D =", float(np.sum(D)/np.sum(T)))
    print(" mean relative D (T>0) =", float(np.mean(D[mask]/T[mask])))
    print()

# Exakte vereinfachte Differenzformeln (wegen C_ab=C_ba):
#
# D_0  = R_EC-R_AB
#      = 1/4 ( C_{1,chi12} - C_{chi3,chi4} )
#
# D_2  = R_EE-R_BB
#      = 1/4 ( C_{1,chi4} + C_{1,chi12}
#              + C_{chi3,chi4} + C_{chi3,chi12} )
#
# D_6  = R_EA-R_BC
#      = 1/4 ( C_{1,chi4} - C_{chi3,chi12} )
#
# D_10 = R_AA-R_CC
#      = 1/4 ( C_{1,chi4} - C_{1,chi12}
#              - C_{chi3,chi4} + C_{chi3,chi12} )
#
# Diese Identitäten sind algebraisch exakt für p,q>3.
