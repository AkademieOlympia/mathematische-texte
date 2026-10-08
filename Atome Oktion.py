#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EABC Framework - Octonionic Periodic Table Mapping (Z = 1 to 93)
Abbildung der Elemente Z=1..93 auf Oktonionen-Vektoren in O ≅ R^8.

Ohne Sage: Oktonionen als 8D-Vektoren (e0..e7); Addition und Norm^2 reichen
für dieses Mapping (keine Multiplikation nötig).
"""

from __future__ import annotations

import numpy as np

# Standardbasis von O ≅ R^8: e[0] = 1, e[1]..e[7] imaginär
e = [np.eye(8, dtype=float)[i] for i in range(8)]


def get_quantum_numbers(Z: int) -> tuple[int, int, int | float, float]:
    """
    Liefert vereinfachte Valenz-Quantenzahlen (n, l, ml, ms) für Z = 1..93.
    """
    if Z <= 2:
        n, l = 1, 0
    elif Z <= 4:
        n, l = 2, 0
    elif Z <= 10:
        n, l = 2, 1
    elif Z <= 12:
        n, l = 3, 0
    elif Z <= 18:
        n, l = 3, 1
    elif Z <= 20:
        n, l = 4, 0
    elif Z <= 30:
        n, l = 4, 2  # d-block
    elif Z <= 36:
        n, l = 4, 1
    elif Z <= 38:
        n, l = 5, 0
    elif Z <= 48:
        n, l = 5, 2
    elif Z <= 54:
        n, l = 5, 1
    elif Z <= 56:
        n, l = 6, 0
    elif Z <= 71:
        n, l = 4, 3  # 14 Lanthanoide (f-block)
    elif Z <= 80:
        n, l = 5, 2
    elif Z <= 86:
        n, l = 6, 1
    elif Z <= 88:
        n, l = 7, 0
    elif Z <= 93:
        n, l = 5, 3  # Actinoide bis Neptunium
    else:
        n, l = 7, 3

    ms = 0.5 if (Z % 2 == 1) else -0.5
    ml = (Z % (2 * l + 1)) - l if l > 0 else 0
    return n, l, ml, ms


def octonion_vector_for_element(Z: int) -> np.ndarray:
    """Konstruiert den oktonionischen Vektor O(Z) ∈ R^8."""
    n, l, ml, ms = get_quantum_numbers(Z)

    vec = n * e[0]
    num_units = min(2 * l + 1, 7)
    for k in range(1, num_units + 1):
        coeff = (ml + ms) * (1.0 / k)
        vec = vec + coeff * e[k]
    return vec


def format_octonion(vec: np.ndarray) -> str:
    """Kompakte Darstellung a0 + a1 e1 + ... (Nullterme weglassen)."""
    parts: list[str] = []
    for i, a in enumerate(vec):
        if abs(a) < 1e-12:
            continue
        if i == 0:
            parts.append(f"{a:.4g}")
        else:
            sign = "+" if a >= 0 and parts else ""
            parts.append(f"{sign}{a:.4g}*e{i}" if sign else f"{a:.4g}*e{i}")
    return " ".join(parts) if parts else "0"


# Demonstration für Schlüssel-Elemente
elements = [(1, "H"), (6, "C"), (26, "Fe"), (57, "La"), (92, "U"), (93, "Np")]

print("=" * 65)
print("OKTONIONISCHE VEKTOREN FÜR DIE ELEMENTE (EABC-MODELL)")
print("=" * 65)
for Z, sym in elements:
    vec = octonion_vector_for_element(Z)
    norm_sq = float(np.dot(vec, vec))
    print(f"Z = {Z:2d} ({sym:2s}) | O(Z) = {format_octonion(vec)}")
    print(f"         | coords = {np.array2string(vec, precision=4, suppress_small=True)}")
    print(f"         | Norm^2 = {norm_sq:.4f}")
    print("-" * 65)
