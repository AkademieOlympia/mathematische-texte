import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Gitter definieren (entsprechend den Achsen Ihrer Plots)
x = np.linspace(0, 5.5, 200)
y = np.linspace(0, 6.0, 200)
X, Y = np.meshgrid(x, y)

# Mathematische Struktur des Wellenfeldes u(x, y)
# Die langwellige Modulation in y-Richtung gekoppelt mit der Hochfrequenz in x-Richtung
f_lang = np.cos(0.5 * Y)
f_hoch = np.sin(4.0 * X) * np.cos(1.5 * X)

# Superposition und Kopplung (erzeugt die charakteristischen Augen/Rauten)
U = 0.5 * (1.0 + f_lang * f_hoch)

# Visualisierung: 3D-Oberflächenplot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Color-Map analog zu Ihren rot-blauen Plots (z.B. 'RdBu_r' oder 'twilight')
surf = ax.plot_surface(X, Y, U, cmap='RdBu_r', edgecolor='none', alpha=0.9)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('u_h(x,y)')
ax.set_title('Analytische Rekonstruktion des Interferenzmusters')
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()