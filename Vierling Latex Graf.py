import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------
# 1. Primzahlen & Vierlinge
# -----------------------------------------
def sieve(n):
    prime_mask = np.ones(n+1, dtype=bool)
    prime_mask[0:2] = False
    for i in range(2, int(np.sqrt(n))+1):
        if prime_mask[i]:
            prime_mask[i*i:n+1:i] = False
    return np.where(prime_mask)[0]

def find_quadruplets(prime_numbers, limit_count=2000):
    prime_set = set(prime_numbers)
    prime_quadruplets = []
    for p in prime_numbers:
        if (p+2 in prime_set and
            p+6 in prime_set and
            p+8 in prime_set):
            prime_quadruplets.append((p, p+2, p+6, p+8))
            if len(prime_quadruplets) >= limit_count:
                break
    return prime_quadruplets

# -----------------------------------------
# Daten
# -----------------------------------------
primes = sieve(500000)
quadruplets = find_quadruplets(primes, 2000)

starts = np.array([q[0] for q in quadruplets])
spacings = np.diff(starts)
spacings = spacings / np.mean(spacings)

# -----------------------------------------
# GUE Vergleich
# -----------------------------------------
s = np.linspace(0,3,200)
wigner = (np.pi/2)*s*np.exp(-np.pi*s**2/4)

# -----------------------------------------
# FIGURE 1: Histogramm
# -----------------------------------------
plt.figure(figsize=(6,4))
plt.hist(spacings, bins=50, density=True, alpha=0.6, label="Vierlinge")
plt.plot(s, wigner, 'r-', label="Wigner (GUE)")
plt.legend()
plt.title("Spacing-Verteilung")
plt.xlabel("s")
plt.ylabel("Dichte")
plt.tight_layout()
plt.savefig("spacing_histogram.png", dpi=300)
plt.close()

# -----------------------------------------
# FIGURE 2: Level Repulsion
# -----------------------------------------
hist, bins = np.histogram(spacings, bins=80, range=(0,2), density=True)
centers = (bins[:-1] + bins[1:]) / 2

plt.figure(figsize=(6,4))
plt.plot(centers, hist, label="Vierlinge")
plt.plot(s, wigner, 'r--', label="Wigner")
plt.legend()
plt.title("Level Repulsion")
plt.xlabel("s")
plt.ylabel("P(s)")
plt.tight_layout()
plt.savefig("level_repulsion.png", dpi=300)
plt.close()

print("Fertig: PNG-Dateien im aktuellen Ordner erzeugt.")