import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def get_hurwitz_units():
    """Die 24 Hurwitz-Einheiten: 8 Lipschitz + 16 Halbe."""
    units = []
    # 8 Lipschitz-Einheiten (±1 auf einer Achse)
    for i in range(4):
        for s in (-1, 1):
            p = [0, 0, 0, 0]
            p[i] = s
            units.append(p)
    # 16 Hurwitz-Halbe: (±1/2, ±1/2, ±1/2, ±1/2)
    for a in (-1, 1):
        for b in (-1, 1):
            for c in (-1, 1):
                for d in (-1, 1):
                    units.append([a/2, b/2, c/2, d/2])
    return np.array(units)

def project_4d_to_3d(points_4d, w_camera=2.0):
    """ Projiziert 4D-Punkte in den 3D-Raum (Stereografische Projektion) """
    projected_points = []
    for p in points_4d:
        x, y, z, w = p
        # Formel für die Projektion: w ist die vierte Dimension
        factor = 1 / (w_camera - w)
        projected_points.append([x * factor, y * factor, z * factor])
    return np.array(projected_points)

# 1. Holen der 24 Hurwitz-Einheiten (aus dem vorherigen Schritt)
units_4d = get_hurwitz_units() 

# 2. Projektion nach 3D
points_3d = project_4d_to_3d(units_4d)

# 3. Plotten
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Farben: Lipschitz-Einheiten (Achsen) vs. Hurwitz-Halbe (Ecken)
colors = ['red' if np.all(np.mod(p, 1) == 0) else 'blue' for p in units_4d]

ax.scatter(points_3d[:, 0], points_3d[:, 1], points_3d[:, 2], 
           c=colors, s=100, edgecolors='black')

ax.set_title("3D-Projektion der 24 Hurwitz-Einheiten (Raumzeit-Keim)")
plt.savefig("Littelwood_Hurwitz.png", dpi=150)
print("Plot gespeichert: Littelwood_Hurwitz.png")