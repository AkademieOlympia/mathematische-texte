import numpy as np

N = 2_000_000

# -------- Primzahlsieb --------
sieve = np.ones(N+10, dtype=bool)
sieve[:2] = False

for i in range(2, int(N**0.5)+1):
    if sieve[i]:
        sieve[i*i:N+1:i] = False

primes = np.nonzero(sieve)[0]


# -------- Familienfunktion --------
def family(p):
    r = p % 12
    if r == 1:
        return 0  # E
    if r == 5:
        return 1  # A
    if r == 7:
        return 2  # B
    if r == 11:
        return 3  # C
    return -1


# -------- Vierlinge suchen --------
quadruples = []

for p in primes:
    if p+8 > N:
        break
    
    if sieve[p+2] and sieve[p+6] and sieve[p+8]:
        
        quadruples.append([
            p,
            p+2,
            p+6,
            p+8,
            family(p),
            family(p+2),
            family(p+6),
            family(p+8)
        ])

quadruples = np.array(quadruples, dtype=np.int64)

print("Anzahl Vierlinge:", len(quadruples))


# -------- speichern --------
np.save("prime_quadruples_EABC_2M.npy", quadruples)

print("Datei gespeichert: prime_quadruples_EABC_2M.npy")