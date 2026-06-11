import matplotlib.pyplot as plt
import numpy as np

def oest_entropy(n):
    # Simulation der OEST-Härte: Hoch bei Rest 7 mod 8, niedrig bei 32-Glattheit
    base = 1.0
    if n % 32 == 0: base = 0.2  # Glättung
    if n % 8 == 7: base = 1.8   # A-Therme Instabilität
    return base * (1 + 0.1 * np.sin(n/10)) # Rauschen der K-Gruppe

x = np.arange(1, 1001)
y = [oest_entropy(n) for n in x]

plt.figure(figsize=(12, 6))
plt.stem(x, y, linefmt='b-', markerfmt='bo', basefmt='r-')
plt.title("OEST-Strukturanalyse: Zähmung der A-Therme durch 32-Glattheit")
plt.xlabel("Natürliche Zahl n")
plt.ylabel("Strukturelle Entropie H(n)")
plt.grid(True)
plt.savefig("OEST_Verifizierung_10000.png")
plt.show()