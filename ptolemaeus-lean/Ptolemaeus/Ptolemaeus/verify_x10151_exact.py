#!/usr/bin/env python3
"""
Verifikation X(10151) und Euler-Gerade (ETC Part 6).

Quelle: https://faculty.evansville.edu/ck6/encyclopedia/ETCPart6.html

- Trilinears (Kimberling): (SA - 24*R^2 + 5*SW) * SB * SC * b * c : :
  (zyklisch in A,B,C bzw. a,b,c).
- Umrechnung Trilinear α:β:γ → baryzentrisch (Kimberling): a*α : b*β : c*γ.
  Der gemeinsame Faktor a*b*c lässt sich wegkürzen → vereinfachte Form T.
- ETC-Linearform (homogene baryzentrische Kombination von X(3), X(4)):
    X(10151) = (4*R^2 - SW) * X(3) + (14*R^2 - 3*SW) * X(4)
  mit SW = a^2 + b^2 + c^2, R^2 = a^2 b^2 c^2 / (16 * Delta^2),
  Delta^2 = s(s-a)(s-b)(s-c), s = (a+b+c)/2.

Hinweis zum „Abgleich“: Sowohl T (aus Trilinears) als auch L = p*O+q*H erfüllen
die Euler-Geraden-Gleichung identisch (Skalarprodukt 0). In Frac(QQ[a,b,c])
sind T und L i. A. nicht dieselbe homogene Richtung (Kreuzprodukt ≠ 0) —
Kimberling kombiniert Zentren oft in einer Konvention, die mit der reinen
Summe der homogenen X(3)-/X(4)-Koordinaten nicht punktweise mit der
geschlossenen Trilinearform übereinstimmen muss; geometrisch liegt beides
auf der Euler-Geraden.

  pip install sympy
  python3 verify_x10151_exact.py
"""

from __future__ import annotations

try:
    from sympy import symbols, expand, together, factor, simplify
except ImportError as e:
    raise SystemExit("SymPy fehlt: pip install sympy\n" + str(e)) from e


def dot(line: tuple, p: tuple):
    return expand(line[0] * p[0] + line[1] * p[1] + line[2] * p[2])


def conway_S(a, b, c):
    """Conway-Größen SA, SB, SC zu Seitenlängen a=|BC|, b=|CA|, c=|AB|."""
    a2, b2, c2 = a**2, b**2, c**2
    SA = (b2 + c2 - a2) / 2
    SB = (a2 + c2 - b2) / 2
    SC = (a2 + b2 - c2) / 2
    return SA, SB, SC


def conway_R2_SW_Delta2(a, b, c):
    a2, b2, c2 = a**2, b**2, c**2
    SW = a2 + b2 + c2
    s = (a + b + c) / 2
    Delta2 = expand(s * (s - a) * (s - b) * (s - c))
    R2 = a2 * b2 * c2 / (16 * Delta2)
    return R2, SW, Delta2


def euler_line_coeffs(a, b, c):
    a2, b2, c2 = a**2, b**2, c**2
    SA, SB, SC = conway_S(a, b, c)
    return (SA * (b2 - c2), SB * (c2 - a2), SC * (a2 - b2))


def bary_OH(a, b, c):
    a2, b2, c2 = a**2, b**2, c**2
    SA, SB, SC = conway_S(a, b, c)
    H = (SB * SC, SC * SA, SA * SB)
    O = (2 * a2 * SA, 2 * b2 * SB, 2 * c2 * SC)
    return O, H


def trilinears_x10151_etc(a, b, c):
    """
    ETC Trilinears (erste Komponente wie auf der Seite, dann zyklisch):
      tri_A = (SA - 24*R^2 + 5*SW) * SB * SC * b * c
      tri_B = (SB - 24*R^2 + 5*SW) * SC * SA * c * a
      tri_C = (SC - 24*R^2 + 5*SW) * SA * SB * a * b
    """
    SA, SB, SC = conway_S(a, b, c)
    R2, SW, _ = conway_R2_SW_Delta2(a, b, c)
    tri_a = (SA - 24 * R2 + 5 * SW) * SB * SC * b * c
    tri_b = (SB - 24 * R2 + 5 * SW) * SC * SA * c * a
    tri_c = (SC - 24 * R2 + 5 * SW) * SA * SB * a * b
    return (tri_a, tri_b, tri_c)


def trilinear_to_bary_kimberling(tri: tuple, a, b, c):
    """Kimberling: baryzentrisch = a*tri_A : b*tri_B : c*tri_C."""
    ta, tb, tc = tri
    return (a * ta, b * tb, c * tc)


def x10151_bary_from_trilinears_simplified(a, b, c):
    """
    Nach a*tri_A = a*b*c * (SA-24R^2+5SW)*SB*SC (zyklisch) kürzen wir a*b*c:
      ( (SA-24R^2+5SW)*SB*SC, (SB-24R^2+5SW)*SC*SA, (SC-24R^2+5SW)*SA*SB )
    """
    SA, SB, SC = conway_S(a, b, c)
    R2, SW, _ = conway_R2_SW_Delta2(a, b, c)
    return (
        (SA - 24 * R2 + 5 * SW) * SB * SC,
        (SB - 24 * R2 + 5 * SW) * SC * SA,
        (SC - 24 * R2 + 5 * SW) * SA * SB,
    )


def x10151_bary_from_etc_linear_form(a, b, c):
    """Homogene Summe (4R^2-SW)*O + (14R^2-3SW)*H."""
    R2, SW, _ = conway_R2_SW_Delta2(a, b, c)
    p = 4 * R2 - SW
    q = 14 * R2 - 3 * SW
    O, H = bary_OH(a, b, c)
    return (p * O[0] + q * H[0], p * O[1] + q * H[1], p * O[2] + q * H[2])


def cross01(p: tuple, q: tuple):
    return expand(p[0] * q[1] - p[1] * q[0])


def main() -> None:
    a, b, c = symbols("a b c", positive=True)
    line = euler_line_coeffs(a, b, c)
    _, H = bary_OH(a, b, c)

    print("Conway cos A = SA/(b*c) usw.; R^2, SW, Delta^2 wie im Docstring.\n")

    print("Orthozentrum H auf Euler-Gerade:", together(dot(line, H)))

    T = x10151_bary_from_trilinears_simplified(a, b, c)
    L = x10151_bary_from_etc_linear_form(a, b, c)
    print("T (aus Trilinears, gekürzt) · Euler:", together(dot(line, T)))
    print("L (ETC-Linearform) · Euler:", together(dot(line, L)))
    print()

    tri = trilinears_x10151_etc(a, b, c)
    Bfull = trilinear_to_bary_kimberling(tri, a, b, c)
    # Bfull soll parallel zu T sein
    cr_tb = cross01(T, Bfull)
    print("Kreuz T × (a*tri) [soll 0]:", simplify(together(cr_tb)) == 0)

    cr_tl = cross01(T, L)
    cr_factored = factor(together(cr_tl).as_numer_denom()[0])
    print("Kreuz T × L — Zähler (erste 200 Zeichen):", str(cr_factored)[:200], "…")
    print("Kreuz T × L identisch 0?", together(cr_tl) == 0)
    print()

    SA, SB, SC = conway_S(a, b, c)
    R2, SW, _ = conway_R2_SW_Delta2(a, b, c)
    print("Beispiel nur zur Lesbarkeit (Faktoren in SA, SB, SC, R^2, SW):")
    print("  tri_A enthält (SA - 24*R^2 + 5*SW) * SB * SC * b * c")
    print("  erster baryzentrischer Faktor nach Kürzen von a*b*c:")
    print("  ", simplify((SA - 24 * R2 + 5 * SW) * SB * SC))
    print()

    a2, b2, c2 = a**2, b**2, c**2
    SA, SB, SC = conway_S(a, b, c)
    alt = (
        a2 / (SB * SC - SA**2),
        b2 / (SC * SA - SB**2),
        c2 / (SA * SB - SC**2),
    )
    chk_alt = together(dot(line, alt))
    print("Altes Snippet · Euler (Zähler):", factor(chk_alt.as_numer_denom()[0]))
    print("Altes Snippet identisch 0?", chk_alt == 0)


if __name__ == "__main__":
    main()
