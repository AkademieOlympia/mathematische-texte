import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# Partition mod 12
def classify(p):
    r = p % 12
    if r == 1: return "E"
    if r == 5: return "A"
    if r == 7: return "B"
    if r == 11: return "C"
    return None

def delta_theta(n):
    fac = sp.factorint(n)
    logE = 0
    logABC = 0
    
    for p, k in fac.items():
        if p <= 3:
            continue
        c = classify(p)
        if c == "E":
            logE += k * np.log(p)
        else:
            logABC += k * np.log(p)
    
    if logE == 0 and logABC == 0:
        return None
    
    delta = logABC - logE
    theta = np.arctan(np.exp(delta))
    return delta, theta

N = 100000
thetas = []
deltas = []

for n in range(1, N+1):
    if np.gcd(n, 6) != 1:
        continue
    res = delta_theta(n)
    if res:
        d, t = res
        deltas.append(d)
        thetas.append(t)

plt.hist(thetas, bins=100)
plt.title("Theta-Verteilung")
plt.show()