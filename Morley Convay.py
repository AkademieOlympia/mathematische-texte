# Visualisierung: Conway-Kreis und Morley-Dreieck (ohne Sage – NumPy + Matplotlib)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.collections import LineCollection

def draw_conway_and_morley(A_coords=(0, 0), B_coords=(6, 1), C_coords=(2, 5)):
    A = np.array(A_coords, dtype=float)
    B = np.array(B_coords, dtype=float)
    C = np.array(C_coords, dtype=float)

    # Seitenlängen: a gegenüber A, b gegenüber B, c gegenüber C
    a = np.linalg.norm(B - C)
    b = np.linalg.norm(A - C)
    c = np.linalg.norm(A - B)
    s = (a + b + c) / 2

    # Inkreismittelpunkt I und Radius r
    I = (a * A + b * B + c * C) / (a + b + c)
    area = np.sqrt(s * (s - a) * (s - b) * (s - c))
    r = area / s

    # Conway-Kreis: R_c = sqrt(r^2 + s^2)
    R_c = np.sqrt(r**2 + s**2)

    def angle(v1, v2):
        n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if n1 < 1e-12 or n2 < 1e-12:
            return 0.0
        cos_a = np.clip(np.dot(v1, v2) / (n1 * n2), -1, 1)
        return np.arccos(cos_a)

    def get_trisector_direction(origin, target1, target2, fraction):
        v1 = target1 - origin
        v2 = target2 - origin
        total = angle(v1, v2)
        cross = v1[0] * v2[1] - v1[1] * v2[0]
        alpha = total * fraction if cross >= 0 else -total * fraction
        rot = np.array([[np.cos(alpha), -np.sin(alpha)], [np.sin(alpha), np.cos(alpha)]])
        return rot @ v1

    def intersect(p1, v1, p2, v2):
        # p1 + t*v1 = p2 + u*v2  =>  t*v1 - u*v2 = p2 - p1
        M = np.column_stack([v1, -v2])
        try:
            t, _ = np.linalg.solve(M, p2 - p1)
            return p1 + t * v1
        except np.linalg.LinAlgError:
            return p1

    # Morley-Dreieck: Schnittpunkte der Winkeldreiteilenden
    vAB_1 = get_trisector_direction(A, B, C, 1/3)
    vAC_1 = get_trisector_direction(A, C, B, 1/3)
    vBA_1 = get_trisector_direction(B, A, C, 1/3)
    vBC_1 = get_trisector_direction(B, C, A, 1/3)
    vCA_1 = get_trisector_direction(C, A, B, 1/3)
    vCB_1 = get_trisector_direction(C, B, A, 1/3)

    M1 = intersect(A, vAB_1, B, vBA_1)
    M2 = intersect(B, vBC_1, C, vCB_1)
    M3 = intersect(C, vCA_1, A, vAC_1)

    # Plot
    fig, ax = plt.subplots(figsize=(8, 8))

    # Basis-Dreieck (kräftig, damit es auf dem großen Plot sichtbar bleibt)
    ax.plot([A[0], B[0], C[0], A[0]], [A[1], B[1], C[1], A[1]], 'k-', lw=2.5, label='Basis-Dreieck', zorder=3)

    # Conway-Kreis (zuerst zeichnen, damit er hinter dem Dreieck liegt)
    circ = Circle(I, R_c, fill=False, color='blue', lw=2.5, label='Conway-Kreis', zorder=1)
    ax.add_patch(circ)

    # Morley-Dreieck
    morley = Polygon([M1, M2, M3], fill=True, facecolor='red', alpha=0.3, edgecolor='red', lw=2, label='Morley-Dreieck', zorder=2)
    ax.add_patch(morley)

    ax.set_aspect('equal')
    ax.legend(loc='upper right')
    ax.set_title('Conway-Kreis und Morley-Dreieck (#Energiedoku)')
    ax.grid(True, alpha=0.3)

    # Achsen so, dass der gesamte Conway-Kreis sichtbar ist (R_c >> Dreieck)
    margin = max(R_c * 0.05, 0.5)
    ax.set_xlim(I[0] - R_c - margin, I[0] + R_c + margin)
    ax.set_ylim(I[1] - R_c - margin, I[1] + R_c + margin)

    return fig, ax

if __name__ == "__main__":
    fig, ax = draw_conway_and_morley()
    plt.savefig("Morley_Conway.png", dpi=150)
    print("Plot gespeichert: Morley_Conway.png")
