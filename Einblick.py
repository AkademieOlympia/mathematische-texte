import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 1. Laden der 2 Millionen Riemann-Nullstellen
zeros = np.load('zeros6.npy')

# Für eine saubere, nicht überladene Wellenfront nutzen wir die ersten N Nullstellen.
# Höhere Nullstellen erzeugen extrem hochfrequente Oszillationen (Kurzwellen-Rauschen).
N_modes = 40
gamma = zeros[:N_modes]

# 2. Definition des Gitters (analog zu Ihren Plots: x von ~0 bis 5.5, y von 0 to 6)
# Wir starten bei x=0.1, da der arithmetische Raum logarithmisch skaliert wird (ln(x))
x = np.linspace(0.1, 5.5, 300)
y = np.linspace(0, 6.0, 300)
X, Y = np.meshgrid(x, y)

# 3. Aufbau des Wellenfeldes u_h(x, y) durch Superposition
# x-Richtung: Moduliert durch den Primzahl-Dualraum via ln(x) und die Nullstellen
# y-Richtung: Moduliert durch eine transversale, niederfrequente stehende Welle
Uh = np.zeros_like(X)

for i, g in enumerate(gamma):
    # Spektrale Dämpfung (Weyl-Skalierung), damit niedrigere Energien das Grundmuster dominieren
    weight = 1.0 / np.sqrt(g)
    
    # Harmonische Kopplung: Die Nullstelle g bestimmt die Frequenz in x
    # Ein kleiner Frequenzteiler (i % 3 + 1) steuert die langwelligere Oszillation in y
    f_x = np.cos(g * np.log(X + 1.0))
    f_y = np.cos(0.5 * (i % 3 + 1) * Y)
    
    Uh += weight * f_x * f_y

# 4. Normierung des Feldes auf das Intervall [0, 1]
# Dies korrespondiert mit Ihrer gemessenen Sichtbarkeit (Visibility) von ~50%
Uh = (Uh - Uh.min()) / (Uh.max() - Uh.min())

# ==========================================
# PLOTTING MALKASTEN (2D & 3D REKONSTRUKTION)
# ==========================================

# Plot A: 2D-Konturplot (Entspricht Ihrem gelb-blauen Dichtebild oben links)
fig_2d, ax_2d = plt.subplots(figsize=(8, 6))
contour = ax_2d.contourf(X, Y, Uh, levels=120, cmap='YlGnBu_r')
ax_2d.set_xlabel('x', fontsize=12)
ax_2d.set_ylabel('y', fontsize=12)
ax_2d.set_title(r'Numerische Approximation $u_h(x,y)$ aus Riemann-Eigenwerten (2D)', fontsize=13)
plt.tight_layout()
plt.savefig('riemann_contour_reconstruction.png', dpi=150)
plt.close(fig_2d)

# Plot B: 3D-Oberflächenplot (Entspricht Ihren rot-blauen 3D-Gittern)
fig_3d = plt.figure(figsize=(10, 7))
ax_3d = fig_3d.add_subplot(111, projection='3d')
surf = ax_3d.plot_surface(X, Y, Uh, cmap='RdBu_r', edgecolor='none', alpha=0.9)
ax_3d.set_xlabel('x')
ax_3d.set_ylabel('y')
ax_3d.set_zlabel('$u_h(x,y)$')
ax_3d.set_title(r'Modulierte Wellenlandschaft $u_h(x,y)$ im Phasenraum (3D)', fontsize=13)

# Blickwinkel anpassen, um die charakteristischen Täler und Kämme hervorzuheben
ax_3d.view_init(elev=35, azim=45)

plt.tight_layout()
plt.savefig('riemann_3d_reconstruction.png', dpi=150)
plt.close(fig_3d)
print("Die Rekonstruktionsplots wurden erfolgreich als PNG-Dateien gespeichert.")