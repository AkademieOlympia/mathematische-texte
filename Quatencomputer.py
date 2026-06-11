import numpy as np
from scipy.linalg import expm, eigh
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Zielzahl: 221 = 13 * 17
N_target = 221

# modulo-60 Klassen
states = [1,7,11,13,17,19,23,29,
          31,37,41,43,47,49,53,59]

N = len(states)

# Hamiltonmatrix
H = np.zeros((N,N), dtype=complex)

for i,a in enumerate(states):
    for j,b in enumerate(states):
        # Reelle Sinus-Terme sind antisymmetrisch und würden beim Hermitisieren verschwinden.
        chirality = 0
        for k in range(1, 8):
            chirality += (1 / k) * np.sin(2 * np.pi * k * (a - b) / 60)

        modular = np.exp(
            2j * np.pi * pow(int(a), int(b), N_target) / N_target
        )

        H[i,j] = 1j * chirality + 0.5 * modular

# Hermitisieren
H = 0.5*(H + H.conj().T)

# Eigenwerte
evals, evecs = eigh(H)

print("Eigenwerte:")
print(evals)

# Zeitentwicklung
t = 1.0
U = expm(-1j * H * t)

# Anfangszustand
psi0 = np.zeros(N, dtype=complex)
psi0[0] = 1.0

# Entwicklung
psi_t = U @ psi0


# Wahrscheinlichkeiten
prob = np.abs(psi_t)**2

plt.plot(states, prob, 'o-')
plt.title("Chirale Resonanzentwicklung")
plt.xlabel("Modulo-60 Zustand")
plt.ylabel("Wahrscheinlichkeit")
plt.grid(True)
plt.savefig("quatencomputer_resonanzentwicklung.png", dpi=150, bbox_inches="tight")
plt.close()