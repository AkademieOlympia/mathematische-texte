import os

import matplotlib.pyplot as plt
import numpy as np


# Zentrale Konfiguration fuer schnelle Vergleiche.
STEPS = 9
Q0 = np.array([5.0, 1.0, 2.0, 4.0], dtype=float)
CASE_VECTORS = {
    "zyklisch": [0.5, 0.5, 0.5, 0.5],
    "x-achse": [0.8, 0.6, 0.0, 0.0],
    "gemischt": [0.6, 0.2, 0.7, 0.3],
}
OUTPUT_FILENAME = "zeitkristall_vergleich.png"


def quat_mult(q, p):
    w1, x1, y1, z1 = q
    w2, x2, y2, z2 = p
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2
    ])


def quat_conj(q):
    return np.array([q[0], -q[1], -q[2], -q[3]])


def norm(q):
    return float(np.sum(q * q))


def unit(q):
    q = np.array(q, dtype=float)
    return q / np.linalg.norm(q)


def projectors(q):
    _, A, B, C = q
    return {
        "A-B": A - B,
        "B-C": B - C,
        "C-A": C - A,
    }


def evolve(Q, U, steps=6):
    xs = [Q.copy()]
    for _ in range(steps):
        Q = quat_mult(U, quat_mult(Q, quat_conj(U)))
        Q = np.round(Q, 12)
        xs.append(Q)
    return xs


def detect_period(orbit, tol=1e-8):
    q0 = orbit[0]
    for k in range(1, len(orbit)):
        if np.allclose(orbit[k], q0, atol=tol, rtol=0):
            return k
    return None


def analyze_case(Q0, U, steps=12):
    orbit = evolve(Q0, U, steps=steps)
    orbit_arr = np.array(orbit, dtype=float)
    proj = [projectors(q) for q in orbit]
    return {
        "orbit": orbit_arr,
        "norms": np.array([norm(q) for q in orbit_arr]),
        "period": detect_period(orbit_arr),
        "A-B": np.array([p["A-B"] for p in proj]),
        "B-C": np.array([p["B-C"] for p in proj]),
        "C-A": np.array([p["C-A"] for p in proj]),
    }


def plot_cases(results, output_path):
    n = len(results)
    fig, axes = plt.subplots(n, 2, figsize=(12, 4 * n), sharex=False)
    if n == 1:
        axes = np.array([axes])

    for row, (name, data) in enumerate(results.items()):
        orbit = data["orbit"]
        xs = np.arange(len(orbit))

        ax_components = axes[row, 0]
        ax_components.plot(xs, orbit[:, 0], marker="o", label="w")
        ax_components.plot(xs, orbit[:, 1], marker="o", label="A")
        ax_components.plot(xs, orbit[:, 2], marker="o", label="B")
        ax_components.plot(xs, orbit[:, 3], marker="o", label="C")
        ax_components.set_title(f"{name}: Komponenten")
        ax_components.set_xlabel("Schritt")
        ax_components.set_ylabel("Wert")
        ax_components.grid(True, alpha=0.3)
        ax_components.legend()

        ax_proj = axes[row, 1]
        ax_proj.plot(xs, data["A-B"], marker="o", label="A-B")
        ax_proj.plot(xs, data["B-C"], marker="o", label="B-C")
        ax_proj.plot(xs, data["C-A"], marker="o", label="C-A")
        ax_proj.plot(xs, data["norms"], linestyle="--", color="black", label="Norm")
        ax_proj.set_title(f"{name}: Projektoren und Norm")
        ax_proj.set_xlabel("Schritt")
        ax_proj.set_ylabel("Wert")
        ax_proj.grid(True, alpha=0.3)
        ax_proj.legend()

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


cases = {name: unit(vec) for name, vec in CASE_VECTORS.items()}

results = {}
for name, U in cases.items():
    data = analyze_case(Q0, U, steps=STEPS)
    results[name] = data
    print(f"\nFall: {name}")
    print(f"Periode: {data['period']}")
    for i, q in enumerate(data["orbit"]):
        print(f"{i} {q} Norm = {data['norms'][i]} Proj = {projectors(q)}")

output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), OUTPUT_FILENAME)
plot_cases(results, output_path)
print(f"\nPlot gespeichert unter: {output_path}")