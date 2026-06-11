import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# Dimension
N = 64

# Zufallsseed
np.random.seed(42)

# Chirale Fourier-artige Kopplung
def chiral_matrix(N):

    H = np.zeros((N,N), dtype=complex)

    for a in range(N):
        for b in range(N):

            val = 0

            for k in range(1,8):

                val += (
                    (1/k)
                    * np.sin(
                        2*np.pi*k*(a-b)/N
                    )
                )

            H[a,b] = val

    return H

# Vier quaternionische Komponenten
HR = chiral_matrix(N)

HI = 0.3 * np.random.randn(N,N)
HJ = 0.3 * np.random.randn(N,N)
HK = 0.3 * np.random.randn(N,N)

# Hermitisieren
HR = 0.5*(HR + HR.conj().T)
HI = 0.5*(HI + HI.T)
HJ = 0.5*(HJ + HJ.T)
HK = 0.5*(HK + HK.T)

# Effektiver komplexer Hamiltonian
H = HR + 1j*HI + 0.5*HJ + 0.5*HK

# Eigenwerte/-vektoren
evals, evecs = eigh(H)

# Kohärenzmaße
coherences = []

for n in range(N):

    psi = evecs[:,n]

    visible = np.real(psi)

    C = (
        np.linalg.norm(visible)**2
        /
        np.linalg.norm(psi)**2
    )

    coherences.append(C)

coherences = np.array(coherences)

print("Mittlere Sichtbarkeit:")
print(np.mean(coherences))

print("Min/Max:")
print(np.min(coherences), np.max(coherences))

# Histogramm
plt.hist(coherences, bins=20)

plt.xlabel("Sichtbare Kohärenz")
plt.ylabel("Anzahl Eigenzustände")

plt.title("Quaternionische Sichtbarkeitsverteilung")

plt.show()