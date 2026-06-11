import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

# --- Daten laden ---
df = pd.read_csv("Vier bis 4700.csv")

# Primzahlen
p = np.array(sorted(set(df.iloc[:,0].values)))[:120]

# Familien bestimmen (mod 12)
def fam(n):
    r = n % 12
    if r == 1: return "E"
    if r == 5: return "A"
    if r == 7: return "B"
    if r == 11: return "C"

F = [fam(x) for x in p]

# --- Zyklische Ordnung ---
order = ["E","A","B","C"]
index = {c:i for i,c in enumerate(order)}

n = len(p)
T = np.zeros((n,n), dtype=complex)

for i in range(n):
    for j in range(max(0,i-1), min(n,i+2)):
        if i != j:
            k = (index[F[j]] - index[F[i]]) % 4
            phi = 0.3 * k + 0.05 * np.log(p[j] - p[i])
            theta = np.log(p[j]) - np.log(p[i]) + phi
            T[i,j] = np.exp(1j * theta)

D = (T + T.conj().T)/2
D = D / np.max(np.abs(D))

eig = np.linalg.eigvalsh(D)

# --- Spacings ---
s = np.diff(np.sort(eig))
s = s / np.mean(s)

# --- Kennzahlen ---
p_small = np.mean(s < 0.2)
moment2 = np.mean(s**2)

print("P(s<0.2):", p_small)
print("<s^2>:", moment2)