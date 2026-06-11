import sympy as sp
from collections import Counter

N = 10**7

primes = list(sp.primerange(2, N + 10))
P = set(primes)

offsets = [0, 2, 6, 8]

quads = []
witness_values = []

for p in primes:
    if p + 8 > N:
        break

    S = sum((p + h) in P for h in offsets)
    W = 3.5 - S

    if p % 12 in (5, 11):
        witness_values.append(W)

    if W < 0:
        quads.append(p)

def orientation(p):
    if p % 12 == 5:
        return "ABCE"
    elif p % 12 == 11:
        return "CEAB"
    return "other"

print("N =", N)
print("Primvierlinge:", len(quads))
print("Orientierungen:", Counter(orientation(p) for p in quads))
print("Residuen mod 210:", Counter(p % 210 for p in quads))

neg = sum(w < 0 for w in witness_values)
print("Negative Witness-Treffer:", neg)

print("Erste Vierlinge:")
for p in quads[:20]:
    print(p, (p, p+2, p+6, p+8), orientation(p), "mod210=", p % 210)