import numpy as np
from collections import Counter
from scipy.linalg import eig

# --------------------------------------------
# 1. Primzahlsieb
# --------------------------------------------

N = 10_000_000


def prime_sieve(n):

    sieve = np.ones(n + 10, dtype=bool)

    sieve[:2] = False

    for p in range(2, int(np.sqrt(n)) + 1):

        if sieve[p]:
            sieve[p*p:n+10:p] = False

    return sieve


print("Erzeuge Primzahlen...")

isprime = prime_sieve(N)

# --------------------------------------------
# 2. Primvierlinge finden
# --------------------------------------------

print("Suche Primvierlinge...")

quadruplets = []

for p in range(5, N - 8):

    if (
        isprime[p]
        and isprime[p+2]
        and isprime[p+6]
        and isprime[p+8]
    ):
        quadruplets.append(p)

print("Gefundene Vierlinge:", len(quadruplets))

# --------------------------------------------
# EABC-Klassen
# --------------------------------------------

def eabc_class(x):

    r = x % 12

    if r == 1:
        return "E"

    if r == 5:
        return "A"

    if r == 7:
        return "B"

    if r == 11:
        return "C"

    return "."

orientations = []

for p in quadruplets:

    pattern = "".join(
        eabc_class(p+h)
        for h in [0,2,6,8]
    )

    orientations.append(pattern)

print("\nOrientierungen:")
print(Counter(orientations))

# --------------------------------------------
# 4. Mittelpunkt-Test
# --------------------------------------------

bad_centers = []

for p in quadruplets:

    m = p + 4

    if m % 3 != 0:
        bad_centers.append(m)

print("\nFehlerhafte Mittelpunkte:", len(bad_centers))

# --------------------------------------------
# 5. mod-210 Analyse
# --------------------------------------------

res210 = Counter()

for p in quadruplets:
    res210[p % 210] += 1

print("\nmod 210:")

for r, c in sorted(res210.items()):
    print(r, c)

# --------------------------------------------
# 6. Quaternionische Phasen
# --------------------------------------------

phase = {
    "ABCE":  1j,
    "CEAB": -1j
}

qsum = 0j

for o in orientations:
    qsum += phase[o]

qnorm = abs(qsum) / len(orientations)

print("\nQuaternionische Phase")
print("Summe =", qsum)
print("Normiert   =", qnorm)

# --------------------------------------------
# 7. Chiraler Bias
# --------------------------------------------

abce = orientations.count("ABCE")
ceab = orientations.count("CEAB")

bias = (abce - ceab) / (abce + ceab)

print("\nChiraler Bias")
print("ABCE =", abce)
print("CEAB =", ceab)
print("Bias =", bias)

# --------------------------------------------
# 8. Hardy-Littlewood-Skalierung
# --------------------------------------------

expected = N / (np.log(N)**4)

C_emp = len(quadruplets) / expected

print("\nHardy-Littlewood")
print("x/log(x)^4 =", expected)
print("C_emp      =", C_emp)

# ============================================
# 9. Quaternionischer EABC-Operator
# ============================================

print("\nBaue EABC-Operator...")

# Quaternionische Phasen
i = 1j

# Zyklische chirale Kopplung
# E -> A -> B -> C -> E

D = np.array([

    [0,  1,  0, -i],
    [1,  0,  1,  0 ],
    [0,  1,  0,  1 ],
    [i,  0,  1,  0 ]

], dtype=complex)

print("\nEABC-Operator:")
print(D)

# --------------------------------------------
# 10. Spektrum
# --------------------------------------------

vals, vecs = eig(D)

print("\nEigenwerte:")

for v in vals:
    print(v)

# --------------------------------------------
# 11. Unitäre Rotationsmatrix
# --------------------------------------------

T = np.array([

    [0,1,0,0],
    [0,0,1,0],
    [0,0,0,1],
    [1,0,0,0]

], dtype=complex)

print("\nT^4 =")

print(np.linalg.matrix_power(T,4))

tvals, _ = eig(T)

print("\nEigenwerte von T:")

for v in tvals:
    print(v)

# ============================================
# 12. Liouville-Funktion
# ============================================

print("\nBerechne Liouville-Funktion...")

spf = np.arange(N+10)

for i2 in range(2, int(np.sqrt(N))+1):

    if spf[i2] == i2:

        for j in range(i2*i2, N+10, i2):

            if spf[j] == j:
                spf[j] = i2


def liouville(n):

    omega = 0

    while n > 1:

        p = spf[n]

        while n % p == 0:
            n //= p
            omega += 1

    return -1 if omega % 2 else 1

# --------------------------------------------
# 13. Liouville-Quaternion-Korrelation
# --------------------------------------------

corr = 0j

for p, o in zip(quadruplets, orientations):

    corr += liouville(p) * phase[o]

corr_norm = abs(corr) / len(quadruplets)

print("\nLiouville-Korrelation")
print("Korrelation =", corr)
print("Normiert    =", corr_norm)

# ============================================
# 14. Fazit
# ============================================

print("\n==============================")
print("FAZIT")
print("==============================")

if qnorm < 0.01:
    print("Quaternionische Phasen löschen sich nahezu aus.")

if abs(bias) < 0.01:
    print("Chirale Balance numerisch stabil.")

if corr_norm < 0.05:
    print("Liouville-Oszillation weitgehend entkoppelt.")
