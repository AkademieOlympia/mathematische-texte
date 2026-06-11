# Visualisierung der Hurwitz-Quaternionen-Gitter (S_7 Basis)
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def visualisiere_hurwitz_raumzeit(p_limit=7, r=5):
    """Visualisiert Hurwitz-Gitterpunkte (Lipschitz + Halbzahlig) als 3D-Scatter."""
    # r: Radius der Betrachtung (r=5 → ~2300 Punkte)
    values = list(range(-r, r + 1))
    half_values = [x + 0.5 for x in range(-r, r)]
    points = []
    
    # Ganzzahlige Punkte (Lipschitz)
    for a in values:
        for b in values:
            for c in values:
                points.append((a, b, c))
                
    # Halbzahlig (Vervollständigung zum Hurwitz-Gitter)
    for a in half_values:
        for b in half_values:
            for c in half_values:
                points.append((a, b, c))

    # 3. 3D-Plot (Projektion der 4D-Struktur auf 3D)
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    xs, ys, zs = zip(*points)
    # Färbung basierend auf der Norm (arithmetische Last)
    norms = [math.sqrt(x**2 + y**2 + z**2) for x, y, z in points]
    
    scatter = ax.scatter(xs, ys, zs, c=norms, cmap='plasma', s=8, alpha=0.5)
    
    ax.set_title(f"Hurwitz-Gitter Projektion (Basis S_{p_limit})")
    ax.set_xlabel('Real (a)')
    ax.set_ylabel('i-Anteil (b)')
    ax.set_zlabel('j-Anteil (c)')
    plt.colorbar(scatter, label='Arithmetische Distanz (Norm)')
    
    plt.show()

# Visualisierung starten
visualisiere_hurwitz_raumzeit(7)