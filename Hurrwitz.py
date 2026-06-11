import numpy as np
from sympy import primerange

# --------------------
# Primzahlvierlinge
# --------------------

def prime_quadruples(limit=2000):
    primes = list(primerange(2, limit))
    P = set(primes)
    
    quads = []
    
    for p in primes:
        if all(q in P for q in [p, p+2, p+6, p+8]):
            quads.append([p,p+2,p+6,p+8])
    
    return quads


V = prime_quadruples(5000)

# --------------------
# EABC Klassifikation mod 12
# --------------------

def family(n):
    
    r = n % 12
    
    if r==1:
        return 0
        
    if r==5:
        return 1
        
    if r==7:
        return 2
        
    if r==11:
        return 3
        
    return -1


# --------------------
# Distanzfunktionen
# --------------------

def d_log(v,w):
    
    return np.linalg.norm(
        np.log(v)-np.log(w)
    )


def d_fam(v,w):
    
    return sum(
        family(a)!=family(b)
        for a,b in zip(v,w)
    )


def d_mod12(v,w):
    
    return sum(
        min(
            abs(a%12-b%12),
            12-abs(a%12-b%12)
        )
        for a,b in zip(v,w)
    )


# --------------------
# Kernel
# --------------------

def weight(v,w,
           sigma=1.5,
           lam=1,
           mu=0.2):
    
    return np.exp(
        -d_log(v,w)**2/sigma**2
    ) * np.exp(
        -lam*d_fam(v,w)
    ) * np.exp(
        -mu*d_mod12(v,w)
    )


# --------------------
# Matrix
# --------------------

N = len(V)

W = np.zeros((N,N))

for i in range(N):
    for j in range(N):
        
        if i!=j:
            W[i,j] = weight(V[i],V[j])

D = np.diag(W.sum(axis=1))

L = D - W

# --------------------
# Spektrum
# --------------------

eigvals = np.linalg.eigvalsh(L)

print(eigvals[:20])