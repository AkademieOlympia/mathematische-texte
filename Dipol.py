import numpy as np
import matplotlib.pyplot as plt

# 1. Gitterdefinition in der komplexen s-Ebene: s = x + iy
x = np.linspace(-6, 12, 500)
y = np.linspace(-8, 8, 500)
X, Y = np.meshgrid(x, y)
S = X + 1j * Y

# Singularitäts-Vermeidung direkt auf den Kernen
S[np.abs(S - 8.0) < 1e-2] = 8.0 + 1e-2 * 1j
S[np.abs(S + 2.0) < 1e-2] = -2.0 + 1e-2 * 1j

# 2. Lokales Zwei-Defekt-Modell omega_loc(s) = -1/(s-8) + 1/(s+2)
omega_loc = -1.0 / (S - 8.0) + 1.0 / (S + 2.0)

# Vektorfeld-Komponenten der Phasenströmung / logarithmischen Ableitung
U = np.real(omega_loc)
V = np.imag(omega_loc)

# Betrag für logarithmische Intensiätsdarstellung
mag = np.sqrt(U**2 + V**2)
mag_clipped = np.clip(mag, 1e-4, 50.0)

# 3. Plot-Generierung
plt.figure(figsize=(11, 8), dpi=300)

# Stromlinien-Darstellung der Strömung
plt.streamplot(
    x, y, U, V, 
    color=np.log(mag_clipped), 
    cmap='plasma', 
    density=1.8, 
    linewidth=1.0, 
    arrowsize=1.2
)

# Positionen der Ladungen (Pol und Nullstelle)
plt.scatter([8.0], [0.0], color='red', s=140, zorder=5, label='Pol (s = 8, q = -1)')
plt.scatter([-2.0], [0.0], color='cyan', s=140, zorder=5, label='Nullstelle (s = -2, q = +1)')

# Schutzzone R0 = 10 um den Hauptpol (s = 8)
protection_circle = plt.Circle((8.0, 0.0), 10.0, color='lime', fill=False, linestyle='--', linewidth=1.5, label='Schutzzone R0 = 10')
plt.gca().add_patch(protection_circle)

plt.axvline(0, color='gray', linestyle=':', alpha=0.6)
plt.axhline(0, color='gray', linestyle=':', alpha=0.6)

plt.title('Lokales Zwei-Defekt-Dipolfeld omega_loc(s) im Spektralraum', fontsize=13, fontweight='bold', pad=15)
plt.xlabel('Re(s)', fontsize=11)
plt.ylabel('Im(s)', fontsize=11)
plt.xlim(-6, 12)
plt.ylim(-8, 8)
plt.legend(loc='upper right', framealpha=0.9)
plt.grid(True, alpha=0.25)

plt.tight_layout()
plt.savefig('ckmrv_dipole_field_flow.png', dpi=300)
print("[SUCCESS] Dipolfeld-Phasenströmung erfolgreich gerendert und als ckmrv_dipole_field_flow.png gespeichert.")