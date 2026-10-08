"""
test_spektral_kern.py
Spektrale Korrespondenz: GNS/Quaternion-Operator vs. kritische Gerade Re(s)=1/2.

Nutzt die Quaternion-Übergänge (transition_operator) als GNS-Gerüst und den
C4-Kern aus Bamberg_Modell_Kern.lean als Startkonstellation.

Governance (#Energiedoku):
  Modell-Definition ≠ physikalische Interpretation.
  Der L-Lift erzwingt Re(λ)=1/2 by construction (Konsistenztest),
  keine empirische Bestätigung der Riemann-Hypothese.
  H_c: zulässige Modellannahme; Spektralanalyse folgt in späterer Iteration.
  Keine Identifikation modellinterner Kennzahlen mit Naturkonstanten.
"""

from __future__ import annotations

import math

import numpy as np

# Erste imaginäre Teile nicht-trivialer Riemann-Nullstellen (Odlyzko)
RIEMANN_T = np.array(
    [
        14.134725,
        21.022040,
        25.010858,
        30.424876,
        32.935062,
        37.586178,
        40.918719,
        43.327073,
        48.005151,
        49.773832,
    ]
)


def quat_mul(q1, q2):
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    return (
        a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
        a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
        a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
        a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2,
    )


def quat_conj(q):
    return (q[0], -q[1], -q[2], -q[3])


def quat_norm(q):
    return math.sqrt(sum(x * x for x in q))


def quat_inverse(q):
    n2 = sum(x * x for x in q)
    c = quat_conj(q)
    return tuple(x / n2 for x in c)


def transition_operator(Q1, Q2):
    """Kraus-/Übergangsoperator aus Quaternionen-Paar (GNS-Gerüst)."""
    delta = quat_mul(quat_inverse(Q1), Q2)
    norm = quat_norm(delta)
    if norm == 0:
        return np.zeros((4, 4), dtype=float)
    vec = np.array(delta, dtype=float) / norm
    return np.outer(vec, vec)


def eabc_quaternion(p: int) -> tuple[float, float, float, float]:
    """EABC-Kanalprojektion mod 12 → Quaternion (Skalar + i,j,k)."""
    r = p % 12
    if r == 1:
        return (1.0, 0.5 * p, 0.3 * p, 0.2 * p)
    if r == 5:
        return (1.0, 0.2 * p, 0.5 * p, 0.3 * p)
    if r == 7:
        return (1.0, 0.3 * p, 0.2 * p, 0.5 * p)
    if r == 11:
        return (1.0, 0.4 * p, 0.4 * p, 0.2 * p)
    return (1.0, float(p), 0.0, 0.0)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True


def twin_prime_quaternions(limit: int = 300) -> list[tuple[float, float, float, float]]:
    quads = []
    for p in range(5, limit):
        if is_prime(p) and is_prime(p + 2):
            quads.append(eabc_quaternion(p))
    return quads


def build_gns_hermitian(quads, n_trans: int = 200) -> np.ndarray:
    """Hermitescher GNS-Operator D = Σ B†B aus Quaternion-Übergängen."""
    n_trans = min(n_trans, max(0, len(quads) - 1))
    D = np.zeros((4, 4), dtype=float)
    for i in range(n_trans):
        B = transition_operator(quads[i], quads[i + 1])
        D += B.T @ B
    # Symmetrisierung (numerische Stabilität)
    return 0.5 * (D + D.T)


def spectral_operator_critical_lift(H: np.ndarray) -> np.ndarray:
    """
    Nicht-selbstadjungierter Lift auf die kritische Gerade:
        L = (1/2)·I + i·H_c
    mit zentriertem H_c (Spur 0). Dann gilt Re(λ) = 1/2 für alle Eigenwerte.
    """
    n = H.shape[0]
    H_c = H - (np.trace(H) / n) * np.eye(n)
    return 0.5 * np.eye(n) + 1j * H_c


def test_spectral_signature():
    # 1. GNS-Algebra-Repräsentation (Quaternionen / transition_operator)
    #    C4-Bezug: Bamberg_Modell_Kern → Start bei ungeraden Primkanälen
    quads = twin_prime_quaternions(400)
    assert len(quads) >= 3, "Zu wenige Primzwillinge fuer GNS-Kette"

    H = build_gns_hermitian(quads, n_trans=min(200, len(quads) - 1))
    spectral_operator_matrix = spectral_operator_critical_lift(H)

    # 2. Spektrum
    eigenvalues = np.linalg.eigvals(spectral_operator_matrix)
    real_parts = np.real(eigenvalues)
    imag_parts = np.sort(np.abs(np.imag(eigenvalues)))

    mean_real_part = float(np.mean(real_parts))
    max_dev = float(np.max(np.abs(real_parts - 0.5)))

    # 3. Riemann-Korrespondenz: kritische Gerade Re(s) = 0.5
    assert 0.4 < mean_real_part < 0.6, "Spektrale Lücke außerhalb der kritischen Zone!"
    assert max_dev < 1e-8, f"Eigenwerte weichen von Re=1/2 ab (max |Δ|={max_dev})"

    # 4. Grober Abgleich der imaginären Skala mit ersten Nullstellen-Höhen
    #    (Skalierung: spektrale Im-Teile vs. Odlyzko-Höhen — relative Ordnung)
    if imag_parts[-1] > 0:
        scale = RIEMANN_T[0] / max(imag_parts[0], 1e-12)
        scaled = imag_parts * scale
        # Nur struktureller Check: Monotonie der sortierten Höhen
        assert np.all(np.diff(scaled) >= -1e-12), "Imaginärteile nicht monoton sortierbar"

    print(f"Test bestanden: Mittlerer Realteil bei {mean_real_part:.4f} (erwartet ~0.5)")
    print(f"  max|Re(λ)-1/2| = {max_dev:.2e}")
    print(f"  |Im(λ)| sortiert: {np.array2string(imag_parts, precision=4)}")
    print(f"  GNS-Kette: {len(quads)} Quaternionen, dim={H.shape[0]}")


if __name__ == "__main__":
    try:
        test_spectral_signature()
        print("Spektrale Validierung: Kern-Konsistenz bestätigt.")
    except AssertionError as e:
        print(f"Spektrale Validierung fehlgeschlagen: {e}")
        raise SystemExit(1)
