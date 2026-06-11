# SageMath 10.9 — Quaternionischer Differentialoperator für das Bamberger Modell
#
# Terminal:
#   sage "Mag Vector.py"
#   sage -python "Mag Vector.py"
#
# Laufende Sage-Sitzung (Notebook / sage:):
#   load("/Pfad/zu/Mag Vector.py")     einmal
#   attach("/Pfad/zu/Mag Vector.py")   neu laden nach jeder Dateiänderung
#
# Reines „python“ ohne Sage: nicht unterstützt (siehe Import unten).

import sys

try:
    from sage.all import SR, QuaternionAlgebra, diff, var  # type: ignore[import-untyped]
except ImportError:
    print(
        "SageMath ist nicht geladen (Modul sage fehlt).\n"
        "Aufruf z.B.:\n"
        '  sage "Mag Vector.py"\n'
        "oder:\n"
        '  sage -python "Mag Vector.py"\n'
        "https://www.sagemath.org/",
        file=sys.stderr,
    )
    sys.exit(2)

# 1. Symbolische Variablen für Raum und Zeit
t, x, y, z, c = var("t x y z c")

# 2. Hamilton-Quaternionen über dem symbolischen Ring
Q = QuaternionAlgebra(SR, -1, -1)
i, j, k = Q.gens()


def _simplify_quaternion_element(F):
    """Vereinfacht F best effort (API je nach Sage-Version unterschiedlich)."""
    try:
        return F.canonical_form()
    except AttributeError:
        pass
    try:
        return F.reduce()
    except AttributeError:
        pass
    try:
        basis = F.parent().basis()
        return sum(comp.simplify_full() * b for comp, b in zip(F.coefficient_tuple(), basis))
    except Exception:
        return F


def apply_quaternion_operator(Phi, Ax, Ay, Az):
    """
    Wendet D = (1/c * d/dt, i*d/dx, j*d/dy, k*d/dz) auf Psi = Phi/c + A an.

    Phi: Skalarpotential (Ausdruck in t,x,y,z)
    Ax, Ay, Az: Komponenten des Vektorpotentials
    """
    Psi = (Phi / c) + Ax * i + Ay * j + Az * k

    D_t = (1 / c) * diff(Psi, t)
    D_x = i * diff(Psi, x)
    D_y = j * diff(Psi, y)
    D_z = k * diff(Psi, z)

    F = D_t + D_x + D_y + D_z
    res = _simplify_quaternion_element(F)

    print("### Ergebnisse der quaternionischen Feld-Ableitung ###")
    print(f"Gesamt-Feldtensor F: {res}")
    print("-" * 30)

    return res


if __name__ == "__main__":
    # Test-Szenario (Potenz mit **; Sage-Preparser akzeptiert alternativ ^ im Notebook)
    test_Phi = x**2 + y**2
    test_Ax = -y * (c / 2)
    test_Ay = x * (c / 2)
    test_Az = 0

    _ = apply_quaternion_operator(test_Phi, test_Ax, test_Ay, test_Az)
