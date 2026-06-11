import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

gamma = np.load("zeros6.npy")[:200]

# primes
def primes(n):
    sieve = np.ones(n+1, bool)
    sieve[:2] = False
    for i in range(2,int(n**0.5)+1):
        if sieve[i]:
            sieve[i*i::i]=False
    return np.nonzero(sieve)[0]

p = primes(20000)[3:]

x = 2*np.log(p)

E = np.cos(np.outer(gamma,x)).sum(axis=0)
E /= np.max(np.abs(E))

phi = x%(2*np.pi)

plt.scatter(phi,E,s=3)
plt.xlabel("phase = log(p²) mod 2π")
plt.ylabel("amplitude")
plt.savefig("Potenzen_Prim.png", dpi=150)
print("Plot gespeichert: Potenzen_Prim.png")