import sympy as sp
import numpy as np

# ---------- Twins ----------
def twin_primes(N):
    twins = []
    for p in range(5, N):
        if sp.isprime(p) and sp.isprime(p+2):
            twins.append(p)
    return twins

# ---------- Klassifikation ----------
def twin_type(p):
    r = p % 12
    if r == 5:
        return "AB"   # (5,7)
    elif r == 11:
        return "CE"   # (11,1)
    else:
        return "other"

# ---------- Quaternion ----------
def Q_twin(p):
    return np.array([p, p+2, 0, 0], dtype=float)  # reduziert

# Quaternion algebra
def quat_mult(a,b):
    w1,x1,y1,z1 = a
    w2,x2,y2,z2 = b
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])

def quat_inv(q):
    conj = np.array([q[0],-q[1],-q[2],-q[3]])
    return conj / np.dot(q,q)

# ---------- BM-Rotation ----------
u = np.array([1,1,1,1],dtype=float)/2
u_inv = quat_inv(u)

def transform(Qv):
    return quat_mult(quat_mult(u, Qv), u_inv)

# ---------- Analyse ----------
N = 200000
twins = twin_primes(N)

count_AB = 0
count_CE = 0

after_AB = 0
after_CE = 0

patterns = []

for p in twins:
    t = twin_type(p)
    
    if t == "AB":
        count_AB += 1
    elif t == "CE":
        count_CE += 1
    
    # Transformation
    Qv = Q_twin(p)
    Qr = transform(Qv)
    
    # Projektion (nur erste beiden Komponenten)
    p_new = int(sp.nextprime(abs(Qr[0])))
    
    t2 = twin_type(p_new)
    
    if t2 == "AB":
        after_AB += 1
    elif t2 == "CE":
        after_CE += 1
    
    patterns.append(t2)

# ---------- Ergebnisse ----------
print("Original:")
print("AB:", count_AB, "CE:", count_CE)

print("\nNach BM-Transformation:")
print("AB:", after_AB, "CE:", after_CE)

# ---------- Clusteranalyse ----------
def cluster_measure(data, window):
    c = 0
    for i in range(len(data)-window):
        w = data[i:i+window]
        c += len(w) - len(set(w))
    return c/len(data)

for W in [50,100,200]:
    print("Cluster W=",W,":", cluster_measure(patterns,W))