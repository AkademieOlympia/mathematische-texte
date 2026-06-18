#!/usr/bin/env python3
"""
Universelle euklidische Hebung — Minimal-Stub.

Kanonsiche Definitionen: collatz_eabc_euklidische_hebung.md

Demonstriert:
  - einen euklidischen Schritt E_Z(x,y) = (y, r) in Z
  - einen euklidischen Schritt in Z[i] (Gaußsche Ganzzahlen)
  - Norm-Defekt D(q) = q - Pi(q) für Hurwitz-Quaternionen

Kein Collatz-Beweis. Kein vollständiger gcd über alle Algebren.
"""

from __future__ import annotations

from dataclasses import dataclass


# ---------------------------------------------------------------------------
# Z — klassischer euklidischer Schritt
# ---------------------------------------------------------------------------


def euclidean_step_z(x: int, y: int) -> tuple[int, int]:
    """E_Z(x, y) = (y, r) mit x = q*y + r und |r| < |y|."""
    if y == 0:
        raise ValueError("y must be nonzero")
    q = x // y
    r = x - q * y
    if abs(r) >= abs(y) and r != 0:
        # Anpassung für negative Reste (klassische Variante)
        if abs(r - y) < abs(r):
            q += 1
            r -= y
        elif abs(r + y) < abs(r):
            q -= 1
            r += y
    return y, r


def defect_z(x: int) -> int:
    """D_R(x) = x - Pi(x); Pi = Identität auf Z."""
    return 0


def epsilon_z(x: int) -> int:
    """ε(x) = N(x) - N(Π(x)); auf Λ_Z = Z stets 0."""
    return 0


# ---------------------------------------------------------------------------
# Z[i] — Gaußsche Ganzzahlen (a + b*i)
# ---------------------------------------------------------------------------

Gauss = tuple[int, int]  # (re, im)


def g_add(a: Gauss, b: Gauss) -> Gauss:
    return (a[0] + b[0], a[1] + b[1])


def g_sub(a: Gauss, b: Gauss) -> Gauss:
    return (a[0] - b[0], a[1] - b[1])


def g_mul(a: Gauss, b: Gauss) -> Gauss:
    ar, ai = a
    br, bi = b
    return (ar * br - ai * bi, ar * bi + ai * br)


def g_conj(a: Gauss) -> Gauss:
    return (a[0], -a[1])


def g_norm_sq(a: Gauss) -> int:
    return a[0] * a[0] + a[1] * a[1]


def g_round_nearest(z: Gauss) -> Gauss:
    """Pi_C: nächste Gaußsche Ganzzahl (Koordinaten-Runden mit lokaler Korrektur)."""
    re, im = z
    candidates = []
    for dr in (-1, 0, 1):
        for di in (-1, 0, 1):
            q = (int(round(re)) + dr, int(round(im)) + di)
            candidates.append(q)
    return min(candidates, key=lambda q: g_norm_sq(g_sub(z, q)))


def g_inv(a: Gauss) -> Gauss:
    n = g_norm_sq(a)
    if n == 0:
        raise ValueError("division by zero")
    c = g_conj(a)
    return (c[0] // n if c[0] % n == 0 else c[0] / n, c[1] // n if c[1] % n == 0 else c[1] / n)


def g_inv_exact(a: Gauss) -> tuple[Gauss, int]:
    """y^{-1} als (conjugate, norm_sq) für exakte Restrechnung."""
    n = g_norm_sq(a)
    if n == 0:
        raise ValueError("division by zero")
    return g_conj(a), n


def epsilon_zi(x: Gauss) -> int:
    """
    ε(x) = N(x) - N(Π(x)) — lokaler Normdefekt.

    Für x ∈ Z[i] (Gitterpunkt): Π(x) = x, also ε(x) = 0.
    Für rationale Koordinaten vor Rundung: N(x) - N(Π(x)) ≥ 0.
    """
    pi_x = g_round_nearest(x)
    return max(0, g_norm_sq(x) - g_norm_sq(pi_x))


def euclidean_step_zi(x: Gauss, y: Gauss) -> tuple[Gauss, Gauss]:
    """E_C(x,y) = (y, r) mit N(r) < N(y) (Gauß-Euklid, ein Schritt)."""
    z_conj, n = g_inv_exact(y)
    # z = x * y^{-1}  →  (x * conj(y)) / N(y)  als rationale Koordinaten
    prod = g_mul(x, z_conj)
    z_rat = (prod[0] / n, prod[1] / n)
    q = g_round_nearest(z_rat)
    r = g_sub(x, g_mul(q, y))
    return y, r


# ---------------------------------------------------------------------------
# Hurwitz-Quaternionen — Defekt D(q) = q - Pi(q)
# ---------------------------------------------------------------------------

Hurwitz = tuple[int, int, int, int]  # Koeffizienten * 2 (halbe Gitterpunkte)


def _h_from_halves(h: Hurwitz) -> tuple[float, float, float, float]:
    return tuple(c / 2 for c in h)  # type: ignore[return-value]


def h_norm_sq(h: Hurwitz) -> int:
    return sum(c * c for c in h) // 4  # N mit halben Koeffizienten


def h_is_hurwitz(h: Hurwitz) -> bool:
    a, b, c, d = h
    all_int = all(x % 2 == 0 for x in h)
    all_half = all(x % 2 == 1 for x in h)
    return (all_int or all_half) and (a + b + c + d) % 2 == 0


def h_round_to_hurwitz(coords: tuple[float, float, float, float]) -> Hurwitz:
    """Pi_H: nächster Hurwitz-Punkt (Brute-Force über Kandidaten)."""
    best: Hurwitz | None = None
    best_d = 10**18
    for a2 in range(int(round(2 * coords[0])) - 2, int(round(2 * coords[0])) + 3):
        for b2 in range(int(round(2 * coords[1])) - 2, int(round(2 * coords[1])) + 3):
            for c2 in range(int(round(2 * coords[2])) - 2, int(round(2 * coords[2])) + 3):
                for d2 in range(int(round(2 * coords[3])) - 2, int(round(2 * coords[3])) + 3):
                    h = (a2, b2, c2, d2)
                    if not h_is_hurwitz(h):
                        continue
                    dist = sum((a2 / 2 - coords[i]) ** 2 for i, coords_i in enumerate(coords))
                    if dist < best_d:
                        best_d = dist
                        best = h
    if best is None:
        raise ValueError("no Hurwitz candidate")
    return best


@dataclass(frozen=True, slots=True)
class HurwitzDefect:
    q: Hurwitz
    pi_q: Hurwitz
    defect_halves: Hurwitz

    @property
    def norm_sq(self) -> int:
        d = self.defect_halves
        return sum(c * c for c in d) // 4


def hurwitz_defect(q: Hurwitz) -> HurwitzDefect:
    """D_H(q) = q - Pi(q) für q in H_H (Koeffizienten als Halbzahlen * 2)."""
    coords = _h_from_halves(q)
    pi_q = h_round_to_hurwitz(coords)
    defect = tuple(q[i] - pi_q[i] for i in range(4))
    return HurwitzDefect(q=q, pi_q=pi_q, defect_halves=defect)


def _demo() -> None:
    print("=== Z ===")
    x, y = 1071, 462
    y1, r1 = euclidean_step_z(x, y)
    print(f"E_Z({x}, {y}) = ({y1}, {r1}), |r| < |y|: {abs(r1) < abs(y)}")

    print("\n=== Z[i] ===")
    gx, gy = (7, 1), (3, 1)  # 7+i, 3+i
    y2, r2 = euclidean_step_zi(gx, gy)
    print(f"E_C({gx}, {gy}) = ({y2}, {r2}), N(r) < N(y): {g_norm_sq(r2) < g_norm_sq(gy)}")

    print("\n=== Hurwitz Defekt ===")
    # 2+i in Hurwitz-Kodierung: (4, 2, 0, 0) halbe Koeffizienten
    q = (5, 3, 3, 3)  # 5/2 + 3/2 i + 3/2 j + 3/2 k  (Summe ganzzahlig)
    d = hurwitz_defect(q)
    print(f"q = {q}, Pi(q) = {d.pi_q}, D(q) = {d.defect_halves}, N(D)^2 ~ {d.norm_sq}")


if __name__ == "__main__":
    _demo()
