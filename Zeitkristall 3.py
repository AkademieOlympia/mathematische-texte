import numpy as np
import math
import csv
from pathlib import Path


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


def quat_norm_sq(q):
    return float(np.sum(q * q))


def quat_normalize(q, tol=1e-15):
    n = np.linalg.norm(q)
    if n < tol:
        raise ValueError("Quaternion mit fast Norm 0 kann nicht normalisiert werden.")
    return q / n


# ============================================================
# Projektionen
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

def evolve(q0, u, steps=300, rounding_digits=12):
    xs = [q0.copy()]
    q = q0.copy()
    uc = quat_conj(u)

    for _ in range(steps):
        q = quat_mult(u, quat_mult(q, uc))
        q = np.round(q, rounding_digits)
        xs.append(q)

    return xs


# ============================================================
# Diagnose
# ============================================================

def detect_period(orbit, tol=1e-8, max_period=None):
    q0 = orbit[0]
    upper = len(orbit) - 1 if max_period is None else min(max_period, len(orbit) - 1)

    for k in range(1, upper + 1):
        if np.allclose(orbit[k], q0, atol=tol, rtol=0):
            return k
    return None


def norm_drift(orbit):
    n0 = quat_norm_sq(orbit[0])
    return float(max(abs(quat_norm_sq(q) - n0) for q in orbit))


def min_return_distance(orbit):
    q0 = orbit[0]
    return float(min(np.linalg.norm(q - q0) for q in orbit[1:]))


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
# Operatoren
# ============================================================

def rotation_quaternion(axis, theta):
    half = theta / 2.0
    c = math.cos(half)
    s = math.sin(half)

    if axis == "x":
        return quat_normalize(np.array([c, s, 0.0, 0.0], dtype=float))
    if axis == "y":
        return quat_normalize(np.array([c, 0.0, s, 0.0], dtype=float))
    if axis == "z":
        return quat_normalize(np.array([c, 0.0, 0.0, s], dtype=float))

    raise ValueError(f"Unbekannte Achse: {axis}")


def special_operators():
    ops = []

    # Triadischer Sonderfall
    ops.append({
        "operator_name": "zyklisch_symmetrisch",
        "family": "triadisch",
        "axis": "-",
        "angle_label": "-",
        "angle_value": None,
        "expected_period": 3,
        "u": quat_normalize(np.array([0.5, 0.5, 0.5, 0.5], dtype=float))
    })

    # Achsenrotationen theta = pi/k
    for axis in ["x", "y", "z"]:
        for k in [2, 3, 4, 5, 6, 7, 8]:
            theta = math.pi / k
            ops.append({
                "operator_name": f"{axis}_rot_pi_{k}",
                "family": "achsial",
                "axis": axis,
                "angle_label": f"pi/{k}",
                "angle_value": theta,
                "expected_period": 2 * k,
                "u": rotation_quaternion(axis, theta)
            })

    # Gemischte Operatoren
    ops.append({
        "operator_name": "gemischt_1",
        "family": "gemischt",
        "axis": "-",
        "angle_label": "-",
        "angle_value": None,
        "expected_period": None,
        "u": quat_normalize(np.array([1.0, 2.0, 3.0, 4.0], dtype=float))
    })
    ops.append({
        "operator_name": "gemischt_2",
        "family": "gemischt",
        "axis": "-",
        "angle_label": "-",
        "angle_value": None,
        "expected_period": None,
        "u": quat_normalize(np.array([2.0, -1.0, 1.0, 3.0], dtype=float))
    })

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
# Studie
# ============================================================

def run_study(steps=500, tol_period=1e-8):
    ops = special_operators()
    states = initial_states()
    results = []

    for op in ops:
        u = op["u"]

        for state_name, q0 in states:
            orbit = evolve(q0, u, steps=steps)
            info = classify_orbit(orbit, tol_period=tol_period)

            row = {
                "operator_name": op["operator_name"],
                "family": op["family"],
                "axis": op["axis"],
                "angle_label": op["angle_label"],
                "angle_value": op["angle_value"],
                "expected_period": op["expected_period"],
                "state_name": state_name,
                "period_found": info["period"],
                "orbit_type": info["orbit_type"],
                "max_norm_drift": info["max_norm_drift"],
                "min_return_distance": info["min_return_distance"],
                "matches_expected_period": (
                    op["expected_period"] == info["period"]
                    if op["expected_period"] is not None
                    else None
                ),
            }
            results.append(row)

    return results


# ============================================================
# Ausgabe
# ============================================================

def print_summary(results):
    print("\n" + "#" * 100)
    print("ZUSAMMENFASSUNG")
    print("#" * 100)

    counts = {}
    for row in results:
        key = row["orbit_type"]
        counts[key] = counts.get(key, 0) + 1

    for k, v in sorted(counts.items()):
        print(f"{k:42s} : {v}")

    print("\n" + "#" * 100)
    print("ACHSIALE TESTS: Erwartete vs. gefundene Perioden")
    print("#" * 100)

    for row in results:
        if row["family"] == "achsial":
            print(
                f"{row['operator_name']:12s} | "
                f"{row['state_name']:18s} | "
                f"erwartet={str(row['expected_period']):>4s} | "
                f"gefunden={str(row['period_found']):>4s} | "
                f"ok={row['matches_expected_period']}"
            )

    print("\n" + "#" * 100)
    print("GEMISCHTE OPERATOREN")
    print("#" * 100)

    for row in results:
        if row["family"] == "gemischt":
            print(
                f"{row['operator_name']:10s} | "
                f"{row['state_name']:18s} | "
                f"Periode={row['period_found']} | "
                f"MinDist={row['min_return_distance']:.6e} | "
                f"Normdrift={row['max_norm_drift']:.3e}"
            )


def save_csv(results, path):
    fieldnames = [
        "operator_name",
        "family",
        "axis",
        "angle_label",
        "angle_value",
        "expected_period",
        "state_name",
        "period_found",
        "orbit_type",
        "max_norm_drift",
        "min_return_distance",
        "matches_expected_period",
    ]

    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


# ============================================================
# Hauptprogramm
# ============================================================

if __name__ == "__main__":
    out_dir = Path("/Users/thomashoffbauer/Desktop/Mathematische Texte")
    csv_path = out_dir / "zeitkristall_studie.csv"

    results = run_study(steps=600, tol_period=1e-8)
    print_summary(results)
    save_csv(results, csv_path)

    print(f"\nCSV gespeichert unter: {csv_path}")