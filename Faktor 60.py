import numpy as np
from scipy.linalg import eigh

# Zielzahl: 221 = 13 * 17
N_target = 221

# Zustandsraum
states = np.arange(2, 120)
M = len(states)

# Komplexer Hamiltonian
H = np.zeros((M, M), dtype=complex)

for i, a in enumerate(states):
    for j, b in enumerate(states):

        # Chirale Moden
        chirality = 0
        for k in range(1, 8):
            chirality += (1 / k) * np.sin(2 * np.pi * k * (a - b) / 60)

        # Echte modulare Dynamik
        modular = np.exp(
            2j * np.pi * pow(int(a), int(b), N_target) / N_target
        )

        H[i, j] = chirality + 0.5 * modular

# Hermitisieren
H = 0.5 * (H + H.conj().T)

# Eigenwerte
evals, evecs = eigh(H)

print("Spektrum:")
print(evals[-20:])