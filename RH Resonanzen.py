"""
RH-Resonanzen: Dichte-Modulation durch Riemann-Nullstellen
und Einfluss auf die Plancksche Strahlungskurve (Bamberg-Modell).
"""
import numpy as np
import matplotlib.pyplot as plt

# Riemann-Nullstellen (Imaginärteile): Sage oder Fallback
try:
    from sage.all import zeta_zeros  # pylint: disable=import-error
    zeros = [float(z.imag()) for z in zeta_zeros()[:30]]
except (ImportError, ModuleNotFoundError):
    # Fallback: erste 30 Nullstellen (ungefähr)
    zeros = [
        14.1347, 21.0220, 25.0109, 30.4249, 32.9351, 37.5862, 40.9187, 43.3271,
        48.0052, 49.7738, 52.9703, 56.4462, 59.3470, 60.8318, 65.1125, 67.0798,
        69.5464, 72.0672, 75.7047, 77.1448, 79.3374, 82.9104, 84.7355, 87.4253,
        88.8091, 92.4919, 94.6513, 95.8706, 98.8312, 101.3179
    ]

# Konstanten für Planck-Spektrum (aus R_K)
c_val = 299792458
kB_val = 1.380649e-23
e_val = 1.602176634e-19
R_K = 25812.80745
h = (e_val**2) * R_K

def planck_spectrum(f, temp):
    """Spektrale Energiedichte u(ν,T) [J·s/m³]."""
    f = np.asarray(f)
    exponent = (h * f) / (kB_val * temp)
    with np.errstate(over='ignore'):
        factor = (8 * np.pi * h * f**3) / (c_val**3)
        distribution = 1.0 / (np.exp(exponent) - 1.0)
    return factor * distribution

# --- 1. Dichte-Modulation durch Riemann-Nullstellen ---
n_pts = 1000
x = np.linspace(0, 50, n_pts)
density_mod = np.zeros(n_pts)
for gamma in zeros:
    density_mod += np.cos(gamma * np.log(x + 1.1))

# --- 2. Planck-Kurve und modulierte Feinstruktur ---
temp_k = 5000  # K
freqs = np.linspace(1e12, 1000e12, n_pts)
planck = planck_spectrum(freqs, temp_k)

# Skala log(freq) auf [0, 50] abbilden wie x
ln_f = np.log(freqs)
x_freq = (ln_f - ln_f.min()) / (ln_f.max() - ln_f.min()) * 50
density_at_freq = np.interp(x_freq, x, density_mod)
# Modulation auf ±2.5 % begrenzen
mod_norm = (density_at_freq - density_at_freq.min()) / (np.ptp(density_at_freq) + 1e-12)
planck_mod = planck * (1 + 0.05 * (mod_norm - 0.5))

# --- 3. Darstellung ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(x, density_mod, color='purple', lw=1.2)
ax1.set_title("Dichte-Modulation durch Riemann-Nullstellen")
ax1.set_xlabel("Skalierte Frequenz / Skala")
ax1.set_ylabel("Resonanz-Amplitude")
ax1.grid(True, alpha=0.3)

ax2.plot(freqs / 1e12, planck_mod * 1e15, label="Moduliert (RH-Feinstruktur)")
ax2.plot(freqs / 1e12, planck * 1e15, 'r--', alpha=0.6, lw=1, label="Klassisch Planck")
ax2.set_title("Feinstruktur der Strahlung (Bamberg-Modell)")
ax2.set_xlabel("Frequenz [THz]")
ax2.set_ylabel(r"$u(\nu,T)$ [$10^{-15}\,\mathrm{J\,s/m^3}$]")
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("RH_Resonanzen.png", dpi=150)
plt.show()
