import numpy as np
import matplotlib.pyplot as plt

# Konstanten
c_val = 299792458       # Lichtgeschwindigkeit [m/s]
kB_val = 1.380649e-23   # Boltzmann-Konstante [J/K]
e_val = 1.602176634e-19 # Elementarladung [C]

# Von-Klitzing-Konstante R_K = h/e² → daraus h = e²·R_K
R_K = 25812.80745  # Von-Klitzing [Ω], festgelegt
h = (e_val**2) * R_K

def planck_spectrum(f, T):
    """Spektrale Energiedichte u(ν,T) der Weißstrahlung, nur aus R_K (und e, c, k_B)."""
    exponent = (h * f) / (kB_val * T)
    with np.errstate(over='ignore'):
        factor = (8 * np.pi * h * f**3) / (c_val**3)
        distribution = 1.0 / (np.exp(exponent) - 1.0)
    return factor * distribution

frequencies = np.linspace(1e12, 1000e12, 1000)  # 1 … 1000 THz
temperatures = [3000, 4000, 5000, 6000]  # K

plt.figure(figsize=(10, 6))
for T in temperatures:
    u = planck_spectrum(frequencies, T)
    plt.plot(frequencies / 1e12, u * 1e15, label=f'T = {T} K')

plt.title(r'Weißstrahlung: Plancksche Strahlungskurven aus der Von-Klitzing-Konstante $R_K$')
plt.xlabel('Frequenz [THz]')
plt.ylabel(r'Energiedichte $u(\nu, T)$ [$10^{-15}\,\mathrm{J\,s/m^3}$]')
plt.legend()
plt.grid(True, which="both", ls="-", alpha=0.3)
plt.tight_layout()
plt.savefig('planck_bamberg_model.png')
plt.show()