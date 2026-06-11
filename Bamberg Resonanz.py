# SageMath 10.5: Prüfung der "Bamberger" Gitter-Resonanz
# Wir simulieren eine Projektion von euklidischen Koordinaten auf eine S3-Sphäre
try:
    from sage.all import (
        RootSystem,
        cos,
        line,
        pi,
        scatter_plot,
        sin,
        sqrt,
        vector,
        zeta,
    )
except ModuleNotFoundError as exc:
    raise SystemExit(
        "SageMath wird fuer dieses Skript benoetigt. "
        "Bitte starte es mit 'sage -python \"Bamberg Resonanz.py\"'."
    ) from exc
import numpy as np

def projektion_bamberg(n_max):
    # Summe der Kehrwerte der Quadrate bis n_max (lokale Dichte)
    lokale_dichte = sum(1/n^2 for n in range(1, n_max + 1))
    fehler = abs(zeta(2) - lokale_dichte)
    return lokale_dichte.n(), fehler.n()

# Beispiel für n=7 (deine Iterationsstufe)
dichte, diff = projektion_bamberg(7)
print(f"Dichte bei Iteration 7: {dichte}")
print(f"Abweichung zum idealen pi^2/6: {diff}")

# Zusatz fuer die Energiedoku: E8-Projektion mit 5-facher Symmetrie
E8 = RootSystem(['E', 8]).ambient_space()
roots = list(E8.roots())

def get_e8_projection():
    coords = [np.array(vector(r), dtype=float) for r in roots]
    phi = float((1 + sqrt(5)) / 2)
    v1 = np.array([1, phi, 0, 1, phi, 0, 1, phi], dtype=float)
    v2 = np.array([phi, 0, 1, phi, 0, 1, phi, 0], dtype=float)

    projected = []
    for c in coords:
        x = float(np.dot(c, v1))
        y = float(np.dot(c, v2))
        projected.append((x, y))
    return projected

points = get_e8_projection()
p = scatter_plot(points, marker='o', size=15, facecolor='blue', edgecolor='cyan')

for i in range(5):
    angle = float(i * 2 * pi / 5)
    p += line([(0, 0), (5 * cos(angle), 5 * sin(angle))], color='gold', thickness=1, alpha=0.5)

p.show(title="E8-Projektion (Bamberg-Resonanz / 5-ary Seed)", aspect_ratio=1, frame=False, axes=False)

print(f"Anzahl der projizierten Knoten: {len(points)}")
print(f"Zeta(2) Dichte-Faktor: {float((pi**2) / 6)}")