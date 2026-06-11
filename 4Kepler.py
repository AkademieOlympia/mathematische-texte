import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ============================================
# Gedankenexperiment:
# 1) Startellipse (1. Primvierling)
# 2) Breakout am e-Punkt (Fokusposition x=f)
# 3) Gerade durch die großen Zahlen 11 und 13
# 4) Landung in der nächsten (größeren) Ellipse
# ============================================

# Erste Ellipse (kleiner Primvierling)
a1 = 13.0
b1 = 4.0
x01, y01 = 0.0, 0.0
f1 = np.sqrt(a1**2 - b1**2)
e1 = f1 / a1

# Zweite Ellipse (nächster Primvierling), deutlich darunter
# so gewählt, dass die horizontale Breakout-Gerade y=0 die Ellipse rechts schneidet
a2 = 24.0
b2 = 12.0
x02, y02 = 0.0, -10.0

print("Ellipse 1: a,b,f,e =", a1, b1, f1, e1)
print("Ellipse 2: a,b,center =", a2, b2, (x02, y02))

# Parameterisierung
theta = np.linspace(0, 2 * np.pi, 1200)
x1 = x01 + a1 * np.cos(theta)
y1 = y01 + b1 * np.sin(theta)
x2 = x02 + a2 * np.cos(theta)
y2 = y02 + b2 * np.sin(theta)

# Vierlingspunkte
quad1_x = np.array([-7, -5, 11, 13], dtype=float)
quad1_y = np.full_like(quad1_x, y01, dtype=float)

quad2_x = np.array([-13, -11, 17, 19], dtype=float)
quad2_y = np.full_like(quad2_x, y02, dtype=float)

# Breakout-Punkt am e-Punkt/Fokus der ersten Ellipse
breakout = np.array([f1, y01], dtype=float)

# Gerade durch die großen Zahlen 11 und 13 -> y = 0
# Landepunkt: rechter Schnitt von y=0 mit Ellipse 2
rhs = 1.0 - ((0.0 - y02) ** 2) / (b2**2)
if rhs < 0:
    raise ValueError("Konfiguration ungueltig: y=0 schneidet Ellipse 2 nicht.")
x_land = x02 + a2 * np.sqrt(rhs)
landing = np.array([x_land, 0.0], dtype=float)
print("Breakout (e-Punkt):", breakout)
print("Landung in Ellipse 2:", landing)

fig, ax = plt.subplots(figsize=(11, 8))
ax.set_xlim(-30, 30)
ax.set_ylim(-24, 10)
ax.set_aspect("equal")
ax.set_title("EABC Gedankenexperiment: Breakout-Linie von Ellipse 1 in Ellipse 2")

# Hilfslinien
ax.axvline(0, color="gray", linestyle="--", linewidth=1, alpha=0.6)
ax.axhline(y01, color="gray", linestyle=":", linewidth=1, alpha=0.6)
ax.axhline(y02, color="gray", linestyle=":", linewidth=1, alpha=0.6)

# Ellipsen
ax.plot(x1, y1, lw=2, label="Ellipse 1 (Startvierling)")
ax.plot(x2, y2, lw=2, label="Ellipse 2 (naechster Vierling)")

# Mittelpunkte und Fokus
ax.scatter([x01, x02], [y01, y02], marker="x", s=90, label="Mittelpunkte")
ax.scatter([breakout[0]], [breakout[1]], marker="D", s=90, label="Breakout e-Punkt")

# Vierlingspunkte
ax.scatter(quad1_x, quad1_y, s=90, label="1. Primvierling")
ax.scatter(quad2_x, quad2_y, s=90, label="2. Primvierling")

for x, y, lbl in zip(quad1_x, quad1_y, ["p-7", "p-5", "p+11", "p+13"]):
    ax.text(x, y + 0.4, lbl, ha="center", fontsize=10)
for x, y, lbl in zip(quad2_x, quad2_y, ["q-13", "q-11", "q+17", "q+19"]):
    ax.text(x, y + 0.4, lbl, ha="center", fontsize=10)

# Markiere die beiden "großen Zahlen" der ersten Reihe
ax.scatter([11, 13], [0, 0], s=120, marker="s", label="große Zahlen 11,13")

# Animationsobjekte
particle1, = ax.plot([], [], "o", ms=8)
trail1, = ax.plot([], [], lw=1.2)
trail1_x, trail1_y = [], []

ray_line, = ax.plot([], [], lw=2.2)
ray_head, = ax.plot([], [], "o", ms=7)

# Frames: erst Bewegung auf Ellipse 1, dann Breakout-Linie
n_ellipse = len(theta)
n_ray = 220
n_total = n_ellipse + n_ray


def update(frame):
    if frame < n_ellipse:
        t = theta[frame]
        px = x01 + a1 * np.cos(t)
        py = y01 + b1 * np.sin(t)
        particle1.set_data([px], [py])
        trail1_x.append(px)
        trail1_y.append(py)
        trail1.set_data(trail1_x, trail1_y)

        # Ray noch aus
        ray_line.set_data([], [])
        ray_head.set_data([], [])
    else:
        # Partikel am Breakout-Punkt "einfrieren"
        particle1.set_data([breakout[0]], [breakout[1]])
        trail1.set_data(trail1_x, trail1_y)

        # Gerade von breakout -> landing animiert ausfahren
        u = (frame - n_ellipse) / max(1, n_ray - 1)
        u = np.clip(u, 0.0, 1.0)
        current = breakout + u * (landing - breakout)
        ray_line.set_data([breakout[0], current[0]], [breakout[1], current[1]])
        ray_head.set_data([current[0]], [current[1]])

    return particle1, trail1, ray_line, ray_head


ani = FuncAnimation(fig, update, frames=n_total, interval=20, blit=True)
ax.legend(loc="upper left")

try:
    plt.show()
except KeyboardInterrupt:
    plt.close("all")
    print("\nAnimation per Ctrl+C beendet.")