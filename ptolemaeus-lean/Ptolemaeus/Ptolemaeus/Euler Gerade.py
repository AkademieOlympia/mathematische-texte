#!/usr/bin/env python3
"""
Euler-Gerade: Umkreismittelpunkt O, Schwerpunkt G und Höhenschnittpunkt H
liegen auf einer Geraden.

Ursprünglich SageMath (`from sage.all import *`). Diese Version nutzt nur
die Standardbibliothek (exakte Brüche) — kein pip nötig.

Start: python3 "Euler Gerade.py"
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class Pt:
    x: Fraction
    y: Fraction

    @staticmethod
    def from_ints(x: int, y: int) -> Pt:
        return Pt(Fraction(x), Fraction(y))


def centroid(a: Pt, b: Pt, c: Pt) -> Pt:
    return Pt((a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3)


def circumcenter(a: Pt, b: Pt, c: Pt) -> Pt:
    """Umkreismittelpunkt über zwei Mittelsenkrechte (allgemein, keine Degeneration)."""
    # Mittelsenkrechte zu AB: |p-A|² = |p-B|²  ⇔  2(b-a)·(p - M_ab) = 0
    # also (b-a)·p = (b-a)·M_ab — Normale der Mittelsenkrechten = Seitenvektor (b-a)
    d_abx, d_aby = b.x - a.x, b.y - a.y
    m_abx = (a.x + b.x) / 2
    m_aby = (a.y + b.y) / 2
    n_abx, n_aby = d_abx, d_aby

    d_acx, d_acy = c.x - a.x, c.y - a.y
    m_acx = (a.x + c.x) / 2
    m_acy = (a.y + c.y) / 2
    n_acx, n_acy = d_acx, d_acy

    # Lineare Gleichung n·p = n·m für beide Linien -> 2x2 lösen
    # n_abx * ox + n_aby * oy = n_abx*m_abx + n_aby*m_aby
    r1 = n_abx * m_abx + n_aby * m_aby
    r2 = n_acx * m_acx + n_acy * m_acy
    det = n_abx * n_acy - n_aby * n_acx
    if det == 0:
        raise ValueError("Punkte kollinear — kein Umkreis.")
    ox = (r1 * n_acy - r2 * n_aby) / det
    oy = (n_abx * r2 - n_acx * r1) / det
    return Pt(ox, oy)


def orthocenter(a: Pt, b: Pt, c: Pt) -> Pt:
    """Höhenschnittpunkt: Schnitt zweier Höhen (Senkrecht zu BC durch A, senkrecht zu AC durch B)."""
    # Höhe durch A: (p-a)·(c-b) = 0  (Normale der Höhengeraden = Richtung BC)
    d_bcx, d_bcy = c.x - b.x, c.y - b.y
    n_ax, n_ay = d_bcx, d_bcy
    d_acx, d_acy = c.x - a.x, c.y - a.y
    n_bx, n_by = d_acx, d_acy

    # n·p = n·(Stützpunkt)
    r1 = n_ax * a.x + n_ay * a.y
    r2 = n_bx * b.x + n_by * b.y
    det = n_ax * n_by - n_ay * n_bx
    if det == 0:
        raise ValueError("Punkte kollinear.")
    hx = (r1 * n_by - r2 * n_ay) / det
    hy = (n_ax * r2 - n_bx * r1) / det
    return Pt(hx, hy)


def cross_z(o: Pt, p: Pt, q: Pt) -> Fraction:
    """Z-Komponente von (p-o) × (q-o); 0 genau dann, wenn O,P,Q kollinear."""
    return (p.x - o.x) * (q.y - o.y) - (p.y - o.y) * (q.x - o.x)


def check_euler_line(a: Pt, b: Pt, c: Pt) -> tuple[Pt, Pt, Pt, bool]:
    o = circumcenter(a, b, c)
    g = centroid(a, b, c)
    h = orthocenter(a, b, c)
    collinear = cross_z(o, g, h) == 0
    return o, g, h, collinear


def check_x10151_note() -> str:
    return (
        "Hinweis: Kimberling X(10151) und volle baryzentrische Algebra sind "
        "in SageMath am einfachsten; hier ist die Euler-Gerade (O–G–H) exakt "
        "mit Bruchrechnung verifiziert."
    )


def main() -> None:
    print("Euler-Gerade — Python (Standardbibliothek, exakte Brüche)\n")

    a = Pt.from_ints(0, 0)
    b = Pt.from_ints(5, 0)
    c = Pt.from_ints(2, 4)

    o, g, h, ok = check_euler_line(a, b, c)
    print("Dreieck A=(0,0), B=(5,0), C=(2,4)")
    print(f"  Umkreismittelpunkt O = ({o.x}, {o.y})")
    print(f"  Schwerpunkt        G = ({g.x}, {g.y})")
    print(f"  Höhenschnittpunkt  H = ({h.x}, {h.y})")
    print(f"  Kollinear (O, G, H)? {ok}")
    print(f"  Kreuzprodukt (soll 0 sein): {cross_z(o, g, h)}")
    print("\n" + check_x10151_note())


if __name__ == "__main__":
    main()
