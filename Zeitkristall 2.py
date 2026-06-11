import numpy as np
import math
from itertools import product


# ============================================================
# Quaternionen-Grundfunktionen
# ============================================================

def quat_mult(q, p):
    w1, x1, y1, z1 = q
    w2, x2, y2, z2 = p
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ], dtype=float)


def quat_conj(q):
    return np.array([q[0], -q[1], -q[2], -q[3]], dtype=float)


def quat_norm(q):
    return float(np.sum(q*q))


def quat_normalize(q, tol=1e-15):
    n = np.linalg.norm(q)
    if n < tol:
        raise ValueError("Quaternion hat fast Norm 0 und kann nicht normalisiert werden.")
    return q / n


# ============================================================
# Projektionen / Observablen
# ============================================================

def projectors(q):
    _, A, B, C = q
    return {
        "A-B": float(A - B),
        "B-C": float(B - C),
        "C-A": float(C - A),
    }


# ============================================================
# Dynamik
# ============================================================

def evolve(q0, u, steps=100, rounding_digits=12):
    xs = [q0.copy()]
    q = q0.copy()

    uc = quat_conj(u)

    for _ in range(steps):
        q = quat_mult(u, quat_mult(q, uc))
        q = np.round(q, rounding_digits)
        xs.append(q)

    return xs


# ============================================================
# Periodensuche / Rückkehrdiagnostik
# ============================================================

def detect_period(orbit, tol=1e-8, max_period=None):
    """
    Sucht die kleinste Periode relativ zum Startzustand orbit[0].
    """
    q0 = orbit[0]
    upper = len(orbit) - 1 if max_period is None else min(max_period, len(orbit)-1)

    for k in range(1, upper + 1):
        if np.allclose(orbit[k], q0, atol=tol, rtol=0):
            return k
    return None


def min_return_distance(orbit):
    """
    Minimaler Abstand eines späteren Zustands zum Startzustand.
    """
    q0 = orbit[0]
    dists = [np.linalg.norm(q - q0) for q in orbit[1:]]
    return float(min(dists)) if dists else 0.0


def norm_drift(orbit):
    """
    Maximale Abweichung der Norm vom Startwert.
    """
    n0 = quat_norm(orbit[0])
    return float(max(abs(quat_norm(q) - n0) for q in orbit))


# ============================================================
# Klassifikation
# ============================================================

def classify_orbit(orbit, tol_period=1e-8):
    period = detect_period(orbit, tol=tol_period)
    drift = norm_drift(orbit)
    min_dist = min_return_distance(orbit)

    if period == 1:
        orbit_type = "stationär"
    elif period is not None and period > 1:
        orbit_type = f"diskreter Zeitkristall (Periode {period})"
    else:
        orbit_type = "Schalenorbit / Präzession / gemischt"

    return {
        "orbit_type": orbit_type,
        "period": period,
        "max_norm_drift": drift,
        "min_return_distance": min_dist,
    }


# ============================================================
# Operatoren erzeugen
# ============================================================

def rotation_quaternion(axis, theta):
    """
    Quaternion für Rotation um Achse 'x', 'y', 'z' im imaginären Unterraum.
    theta in Radiant.
    """
    half = theta / 2.0
    c = math.cos(half)
    s = math.sin(half)

    if axis == "x":
        return quat_normalize(np.array([c, s, 0.0, 0.0], dtype=float))
    elif axis == "y":
        return quat_normalize(np.array([c, 0.0, s, 0.0], dtype=float))
    elif axis == "z":
        return quat_normalize(np.array([c, 0.0, 0.0, s], dtype=float))
    else:
        raise ValueError(f"Unbekannte Achse: {axis}")


def special_operators():
    """
    Einige bewusst gewählte Operatoren.
    """
    ops = []

    # der bekannte "zyklische" Fall
    ops.append(("zyklisch_symmetrisch", quat_normalize(np.array([0.5, 0.5, 0.5, 0.5], dtype=float))))

    # Achsenrotationen mit Winkeln, die oft endliche Ordnung liefern
    for axis in ["x", "y", "z"]:
        for frac_name, frac in [("pi_2", 0.5), ("pi_3", 1/3), ("pi_4", 0.25), ("pi_5", 0.2)]:
            theta = math.pi * frac
            name = f"{axis}_rot_{frac_name}"
            ops.append((name, rotation_quaternion(axis, theta)))

    # ein gemischter Operator
    ops.append(("gemischt_1", quat_normalize(np.array([1.0, 2.0, 3.0, 4.0], dtype=float))))
    ops.append(("gemischt_2", quat_normalize(np.array([2.0, -1.0, 1.0, 3.0], dtype=float))))

    return ops


# ============================================================
# Startzustände
# ============================================================

def initial_states():
    return [
        ("q0_ref", np.array([5.0, 1.0, 2.0, 4.0], dtype=float)),
        ("q0_stationaer_test", np.array([5.0, 1.0, 1.0, 1.0], dtype=float)),
        ("q0_partial_deg", np.array([5.0, 2.0, 2.0, 1.0], dtype=float)),
        ("q0_generic_1", np.array([5.0, 3.0, 7.0, 11.0], dtype=float)),
        ("q0_generic_2", np.array([5.0, -2.0, 4.0, 1.0], dtype=float)),
    ]


# ============================================================
# Bericht
# ============================================================

def print_short_orbit(orbit, max_lines=6):
    limit = min(len(orbit), max_lines)
    for i in range(limit):
        q = orbit[i]
        print(f"    n={i:2d}  Q={q}  N={quat_norm(q):.10f}  Proj={projectors(q)}")
    if len(orbit) > limit:
        print("    ...")


def run_study(steps=120, tol_period=1e-8, show_examples=True):
    ops = special_operators()
    states = initial_states()

    results = []

    for op_name, u in ops:
        for state_name, q0 in states:
            orbit = evolve(q0, u, steps=steps)
            info = classify_orbit(orbit, tol_period=tol_period)

            row = {
                "operator": op_name,
                "state": state_name,
                "period": info["period"],
                "orbit_type": info["orbit_type"],
                "max_norm_drift": info["max_norm_drift"],
                "min_return_distance": info["min_return_distance"],
            }
            results.append(row)

            print("=" * 90)
            print(f"Operator: {op_name} | Zustand: {state_name}")
            print(f"Typ: {info['orbit_type']}")
            print(f"Periode: {info['period']}")
            print(f"Max. Normdrift: {info['max_norm_drift']:.3e}")
            print(f"Min. Rückkehrabstand: {info['min_return_distance']:.6e}")

            if show_examples:
                print_short_orbit(orbit, max_lines=5)

    return results


# ============================================================
# Zusammenfassung
# ============================================================

def summarize_results(results):
    print("\n" + "#" * 90)
    print("ZUSAMMENFASSUNG")
    print("#" * 90)

    counts = {}
    for row in results:
        key = row["orbit_type"]
        counts[key] = counts.get(key, 0) + 1

    for k, v in sorted(counts.items(), key=lambda x: x[0]):
        print(f"{k:40s} : {v}")

    print("\nFälle mit endlicher Periode:")
    for row in results:
        if row["period"] is not None:
            print(
                f"  Operator={row['operator']:20s} "
                f"Zustand={row['state']:20s} "
                f"Periode={row['period']}"
            )


if __name__ == "__main__":
    results = run_study(steps=150, tol_period=1e-8, show_examples=False)
    summarize_results(results)