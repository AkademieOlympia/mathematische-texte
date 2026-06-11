import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

# --- Zeta ---
gamma = np.sort(gamma)
s_zeta = np.diff(gamma)
s_zeta /= np.mean(s_zeta)

# --- Bamberg ---
eig = np.sort(eigvals)
s_b = np.diff(eig)
s_b /= np.mean(s_b)

# --- Histogramm ---
plt.hist(s_zeta, bins=60, density=True, alpha=0.5, label="Zeta")
plt.hist(s_b, bins=60, density=True, alpha=0.5, label="Bamberg")
plt.legend()
plt.title("Zeta vs Bamberg")
plt.show()

# --- KS-Test ---
stat, pval = ks_2samp(s_zeta, s_b)

print("KS-Statistik:", stat)
print("p-Wert:", pval)