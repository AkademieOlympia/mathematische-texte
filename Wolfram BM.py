import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def quaternion_automaton_step(current_grid_size, resolution=4):
    """
    Simuliert die Evolution der Primzahl-Dichte in einem wachsenden 
    quaternionischen Gitter (vereinfachtes 2D-Projektionsmodell).
    resolution: Punkte pro Einheit (4 = 4x feiner als Gitter mit Schritt 1)
    """
    # Feineres Gitter: mehr Punkte pro Einheit
    step = 1.0 / resolution
    x = np.arange(-current_grid_size, current_grid_size + step/2, step)
    y = np.arange(-current_grid_size, current_grid_size + step/2, step)
    X, Y = np.meshgrid(x, y)
    
    # Norm im Hurwitz-Sinn (vereinfacht: x^2 + y^2 + z^2 + w^2)
    # Wir betrachten einen 2D-Schnitt für die Visualisierung
    Norms = X**2 + Y**2
    Norms_int = np.round(Norms).astype(int)  # Rundung für feines Gitter
    
    # Regel: Punkt ist aktiv (1), wenn Norm eine Primzahl ist
    def is_prime_simple(n):
        if n < 2: return False
        for i in range(2, int(np.sqrt(n)) + 1):
            if n % i == 0: return False
        return True
    
    v_is_prime = np.vectorize(is_prime_simple)
    Active_Points = v_is_prime(Norms_int)
    
    return Active_Points

def quaternion_automaton_step_3d(current_grid_size, resolution=2):
    """
    3D-Version: Volle Norm x²+y²+z² im Hurwitz-Sinn.
    resolution niedriger halten (2D wächst mit n², 3D mit n³).
    """
    step = 1.0 / resolution
    x = np.arange(-current_grid_size, current_grid_size + step/2, step)
    y = np.arange(-current_grid_size, current_grid_size + step/2, step)
    z = np.arange(-current_grid_size, current_grid_size + step/2, step)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    Norms = X**2 + Y**2 + Z**2
    Norms_int = np.round(Norms).astype(int)
    
    def is_prime_simple(n):
        if n < 2: return False
        for i in range(2, int(np.sqrt(n)) + 1):
            if n % i == 0: return False
        return True
    
    v_is_prime = np.vectorize(is_prime_simple)
    Active = v_is_prime(Norms_int)
    
    # Koordinaten der aktiven (Primzahl-)Punkte
    idx = np.where(Active)
    return X[idx], Y[idx], Z[idx]

# Wir simulieren das Wachstum über 3 "Zeitstufen"
plt.figure(figsize=(15, 5))
for i, size in enumerate([10, 20, 40]):
    grid = quaternion_automaton_step(size, resolution=4)
    plt.subplot(1, 3, i+1)
    plt.imshow(grid, cmap='magma', interpolation='nearest')
    plt.title(f"Evolution Stufe {i+1} (R={size})")
    plt.axis('off')

plt.tight_layout()
plt.savefig('quaternion_nks_evolution.png', dpi=150)
print("Plot gespeichert: quaternion_nks_evolution.png")

# 3D-Version: Primzahl-Punkte im vollen 3D-Raum (kleinere Radien wegen n³-Skalierung)
fig3d = plt.figure(figsize=(15, 5))
for i, size in enumerate([6, 10, 14]):
    xp, yp, zp = quaternion_automaton_step_3d(size, resolution=2)
    ax = fig3d.add_subplot(1, 3, i+1, projection='3d')
    ax.scatter(xp, yp, zp, c=zp, cmap='magma', s=1, alpha=0.6)
    ax.set_title(f"3D Evolution Stufe {i+1} (R={size})")
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
plt.tight_layout()
plt.savefig('quaternion_nks_evolution_3d.png', dpi=150)
print("Plot gespeichert: quaternion_nks_evolution_3d.png")

plt.show()