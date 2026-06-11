import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Kein GUI – Plot wird nur gespeichert
import matplotlib.pyplot as plt

# Nullstellen laden
gamma = np.load("zeros6.npy")[:200]

# Vierlinge
df = pd.read_csv("prime_quadruplets.csv")

p1 = df["p"].values
p2 = df["p2"].values
p3 = df["p6"].values
p4 = df["p8"].values

# Produktphase (log-Summe statt log(Produkt) um int64-Overflow zu vermeiden)
x = np.log(p1.astype(float)) + np.log(p2.astype(float)) + np.log(p3.astype(float)) + np.log(p4.astype(float))

# Interferenz
E = np.cos(np.outer(gamma, x)).sum(axis=0)

E /= np.max(np.abs(E))

phi = x % (2*np.pi)

plt.scatter(phi, E, s=3)
plt.xlabel("phase = log(p1 p2 p3 p4) mod 2π")
plt.ylabel("resonance amplitude")
plt.title("Resonance spectrum of prime quadruplets")
plt.savefig("Vierlinge_Eigenfrequenzen.png", dpi=150)
print("Plot gespeichert: Vierlinge_Eigenfrequenzen.png")