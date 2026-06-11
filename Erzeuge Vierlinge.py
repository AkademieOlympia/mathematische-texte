import csv
import numpy as np

LIMIT = 200000000   # Grenze erhöhen für mehr Vierlinge

def sieve(n):
    sieve = np.ones(n+1, dtype=bool)
    sieve[:2] = False
    for i in range(2, int(n**0.5)+1):
        if sieve[i]:
            sieve[i*i:n+1:i] = False
    return sieve

print("Erzeuge Primzahlsieb...")
prime = sieve(LIMIT)

print("Suche Vierlinge...")

quartets = []

for p in range(5, LIMIT-8):
    if prime[p] and prime[p+2] and prime[p+6] and prime[p+8]:
        quartets.append((p,p+2,p+6,p+8))

print("Gefunden:", len(quartets))

print("Schreibe CSV...")

with open("prime_quadruplets.csv","w",newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["p","p2","p6","p8"])
    writer.writerows(quartets)

print("Fertig.")