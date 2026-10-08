#!/usr/bin/env python3
"""
Minimaler numerischer Test des zyklengekoppelten Floquet-Transferoperators
für q=2 und die ersten fünf arithmetischen Schalen.

Modellannahme:
    A e_m = sigma_3(m) e_m
    K_{F,r}(z)e_m =
        kappa * omega^r /
        ((sigma_3(m)+mu)^2 * (sigma_3(m)-z)) e_m

mit omega=-1, r=0,1.

Für einen Modus pro Schale gilt:
    D_r(z) = prod_m [1 - kappa*omega^r /
                     ((sigma_3(m)+mu)^2*(sigma_3(m)-z))]

Die Nullstellen sind daher exakt
    z_{m,r} = sigma_3(m)
              - kappa*omega^r/(sigma_3(m)+mu)^2.

Wichtig:
- Dies ist ein endliches, diagonales Toy-Modell.
- Es beweist keine L-Funktions- oder Modularitätsidentität.
- Ohne zusätzliche nichtselbstadjungierte/komplexe Kopplung liegen die
  Resonanznullstellen auf der reellen Achse.
"""

import math
import cmath

try:
    # Optional SageMath
    from sage.all import divisors
    def sigma3(n):
        return int(sum(d**3 for d in divisors(n)))
except Exception:
    def divisors_py(n):
        out = []
        for d in range(1, int(math.isqrt(n)) + 1):
            if n % d == 0:
                out.append(d)
                if d*d != n:
                    out.append(n // d)
        return out

    def sigma3(n):
        return sum(d**3 for d in divisors_py(n))


def shell_data(M=5):
    return [(m, sigma3(m)) for m in range(1, M + 1)]


def k_eigenvalue(a_m, z, r, mu=1.0, kappa=1.0):
    omega = -1.0
    return kappa * (omega ** r) / ((a_m + mu)**2 * (a_m - z))


def sector_determinant(z, r, M=5, mu=1.0, kappa=1.0):
    D = 1.0 + 0.0j
    for _, a_m in shell_data(M):
        D *= (1.0 - k_eigenvalue(a_m, z, r, mu=mu, kappa=kappa))
    return D


def full_determinant(z, M=5, mu=1.0, kappa=1.0):
    return (
        sector_determinant(z, 0, M=M, mu=mu, kappa=kappa)
        * sector_determinant(z, 1, M=M, mu=mu, kappa=kappa)
    )


def exact_sector_roots(r, M=5, mu=1.0, kappa=1.0):
    omega = -1.0
    roots = []
    for m, a_m in shell_data(M):
        z = a_m - kappa * (omega ** r) / (a_m + mu)**2
        roots.append((m, a_m, z))
    return roots


if __name__ == "__main__":
    M = 5
    mu = 1.0
    kappa = 1.0

    print("Erste fünf Schalen:")
    for m, a in shell_data(M):
        print(f"m={m:2d}  sigma_3(m)={a:6d}")

    for r in (0, 1):
        print(f"\nSektor r={r}:")
        for m, a, z in exact_sector_roots(r, M=M, mu=mu, kappa=kappa):
            shift = z - a
            print(
                f"m={m:2d}  Pol={a:12.8f}  "
                f"Nullstelle={z:16.12f}  Verschiebung={shift:+.6e}"
            )

    # Plausibilitätscheck: Determinante an den analytischen Nullstellen
    print("\nNumerischer Nullstellentest:")
    for r in (0, 1):
        for m, a, z in exact_sector_roots(r, M=M, mu=mu, kappa=kappa):
            val = sector_determinant(z, r, M=M, mu=mu, kappa=kappa)
            print(f"r={r}, m={m}: |D_r(z_root)| = {abs(val):.3e}")
