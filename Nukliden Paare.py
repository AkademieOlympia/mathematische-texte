"""
EABC Nuklid-zu-Hurwitz-Mapping.

Hurwitz-Quaternionen H = D1 ∪ D2 verbieten gemischte Koordinaten:
  D1 — alle vier Komponenten ganzzahlig
  D2 — alle vier Komponenten in Z + 1/2
"""

from __future__ import annotations

from fractions import Fraction
from typing import Literal


Sector = Literal["D1 (Ganzzahlig)", "D2 (Halbganzzahlig)"]


def half_shift(x: int) -> Fraction:
    """Projiziert x ∈ ℤ strikt nach ℤ + 1/2."""
    if x % 2 != 0:
        return Fraction(x, 2)
    return Fraction(x + 1, 2)


def is_hurwitz(q: tuple[Fraction, Fraction, Fraction, Fraction]) -> bool:
    """True genau dann, wenn q ∈ D1 oder q ∈ D2 (keine Mischkoordinaten)."""
    dens = {c.denominator for c in q}
    if dens == {1}:
        return True
    if dens == {2} and all(c.numerator % 2 != 0 for c in q):
        return True
    return False


def map_nuclide_to_hurwitz(Z: int, N: int):
    """
    Mappt ein Nuklidpaar (Z, N) auf die Hurwitz-Klassen D1/D2.

    Even-Even → D1 via ganzzahlige Halbierung.
    Sonst → D2 via komponentenweise half_shift (alle in ℤ + 1/2).
    """
    A = Z + N
    I = N - Z

    if Z % 2 == 0 and N % 2 == 0:
        q = (
            Fraction(Z, 2),
            Fraction(N, 2),
            Fraction(A, 2),
            Fraction(I, 2),
        )
        sector: Sector = "D1 (Ganzzahlig)"
    else:
        q = (
            half_shift(Z),
            half_shift(N),
            half_shift(A),
            half_shift(I),
        )
        sector = "D2 (Halbganzzahlig)"

    if not is_hurwitz(q):
        raise ValueError(f"Abbildung verletzt Hurwitz-Zugehörigkeit: q={q}")

    norm = sum(c * c for c in q)
    return sector, q, norm


def _fmt_q(q: tuple[Fraction, Fraction, Fraction, Fraction]) -> str:
    return "(" + ", ".join(str(float(c) if c.denominator == 2 else int(c)) for c in q) + ")"


nuclides = [
    (8, 8, "C12/O16-Typ (Even-Even)"),
    (1, 0, "Proton p (Odd-A)"),
    (1, 1, "Deuterium d (Odd-Odd)"),
    (28, 50, "Ni-78 (Doppelt Magisch)"),
]

print("--- Korrigiertes EABC Nuklid-zu-Hurwitz-Mapping ---")
for Z, N, label in nuclides:
    sec, q_vec, nrm = map_nuclide_to_hurwitz(Z, N)
    print(
        f"Nuklid (Z={Z:2d}, N={N:2d}) [{label:26s}] -> {sec} | "
        f"q = {_fmt_q(q_vec)} | Norm = {float(nrm):.2f} | Hurwitz = {is_hurwitz(q_vec)}"
    )
