"""
CollatzXi — Größen Ξ, Ξ_* am ersten A_K>1-Austritt.

    U_K = 1 + sum_{j=1}^{K-1} (1 - 2^{-s_{j-1}}) A_j
    C_K = U_K - 2^{-s_{K-1}} A_K   (strikt C_K < U_K)
    Ξ   = n0 (A_K - 1) / K
    Ξ_* = n0 (A_K - 1) / U_K

Ξ_* ≥ 1 ist keine unabhängige Majorante: aus
    U_K = A_K (n_K + 2^{-s_{K-1}}) - n_0
folgt am Crossing A_K > 1 (n_0, n_K ungerade)
    Ξ_* ≥ 1  ⟺  n_K < n_0.
Die unabhängige Präfix-Majorante ist M_K in CollatzMajorant.py.

Epistemik: exakte Diagnose. Census ist kein Collatz-Beweis.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from CollatzBallotTree import actual_block


@dataclass
class ExitCert:
    n0: int
    K: int
    L: int
    R: int
    S: int
    Q: int
    r_exit: int
    s_exit: int
    A: Fraction
    U: Fraction
    C: Fraction
    nK: int
    Xi: Fraction
    Xi_star: Fraction

    @property
    def descended(self) -> bool:
        return self.nK < self.n0

    @property
    def certified(self) -> bool:
        return self.Xi_star >= 1


def first_exit_cert(n0: int, max_blocks: int = 400) -> ExitCert:
    """Erster Block mit A_K > 1. n0 ungerade > 1."""
    if n0 <= 1 or n0 % 2 == 0:
        raise ValueError("n0 must be odd and > 1")
    n = n0
    R = Q = B = 0
    U = Fraction(1)
    A_prev: Fraction | None = None
    s_prev: int | None = None
    for k in range(1, max_blocks + 1):
        r, s, nxt = actual_block(n)
        B = 3**r * B + (3**r - 2**r) * (1 << Q)
        R += r
        Q += r + s
        A = Fraction(1 << Q, 3**R)
        C = Fraction(B, 3**R)
        if k >= 2:
            assert A_prev is not None and s_prev is not None
            U += (1 - Fraction(1, 1 << s_prev)) * A_prev
        C_from_U = U - A / Fraction(1 << s)
        if C_from_U != C:
            raise RuntimeError(f"C vs U mismatch n0={n0} k={k}: {C} != {C_from_U}")
        expected = A * nxt
        if n0 + C != expected:
            raise RuntimeError(
                f"affine mismatch n0={n0} k={k}: {n0}+{C} != {A}*{nxt}"
            )
        if A > 1:
            Am1 = A - 1
            return ExitCert(
                n0=n0,
                K=k,
                L=k - 1,
                R=R,
                S=Q - R,
                Q=Q,
                r_exit=r,
                s_exit=s,
                A=A,
                U=U,
                C=C,
                nK=nxt,
                Xi=n0 * Am1 / k,
                Xi_star=n0 * Am1 / U,
            )
        A_prev = A
        s_prev = s
        n = nxt
        if n == 1:
            break
    raise RuntimeError(f"no A>1 exit for n0={n0} within {max_blocks} blocks")


def census_xi_star(X: int, keep_smallest: int = 25) -> dict:
    """Ungerade n in (1, X]. Exakte Brüche, Speicher O(keep_smallest)."""
    import heapq

    n_odd = 0
    n_xi_lt1 = 0
    n_star_lt1: list[ExitCert] = []
    n_not_descended: list[ExitCert] = []
    # Max-Heap der keep kleinsten Ξ* via negative keys
    small: list[tuple[float, int, ExitCert]] = []
    min_star = None
    arg_min = None
    for n0 in range(3, X + 1, 2):
        cert = first_exit_cert(n0)
        n_odd += 1
        if cert.Xi < 1:
            n_xi_lt1 += 1
        if cert.Xi_star < 1:
            n_star_lt1.append(cert)
        if not cert.descended:
            n_not_descended.append(cert)
        val = float(cert.Xi_star)
        if min_star is None or val < min_star:
            min_star = val
            arg_min = n0
        item = (-val, n0, cert)
        if len(small) < keep_smallest:
            heapq.heappush(small, item)
        elif val < -small[0][0]:
            heapq.heappushpop(small, item)
        if n0 % 1_000_000 == 1:
            print(f"    … n0={n0} minΞ*={min_star:.6f}", flush=True)
    smallest = sorted(((-neg, n0, c) for neg, n0, c in small), key=lambda t: t[0])
    return {
        "X": X,
        "n_odd": n_odd,
        "n_Xi_lt1": n_xi_lt1,
        "n_Xistar_lt1": n_star_lt1,
        "n_not_descended": n_not_descended,
        "min_Xistar": min_star,
        "argmin_Xistar": arg_min,
        "smallest": smallest,
    }


def fmt_cert(c: ExitCert) -> str:
    return (
        f"n0={c.n0} K={c.K} L={c.L} A-1={float(c.A-1):.6f} "
        f"U={float(c.U):.6f} Ξ={float(c.Xi):.6f} Ξ*={float(c.Xi_star):.6f} "
        f"nK={c.nK} desc={c.descended} cert={c.certified}"
    )


def main() -> None:
    print("=== Drei Ξ<1-Fälle, exakte Brüche ===", flush=True)
    reported = {
        27: (15, 0.504389, 2.765027, 2.736262, 23),
        31: (14, 0.974824, 2.485655, 5.490521, 23),
        63: (12, 0.420958, 3.026151, 1.669280, 61),
    }
    for n0, (K, xi, U, xs, nK) in reported.items():
        c = first_exit_cert(n0)
        print(" ", fmt_cert(c), flush=True)
        print(
            f"    A={c.A}  U={c.U}  C={c.C}  Ξ*={c.Xi_star}",
            flush=True,
        )
        assert c.K == K, (n0, c.K, K)
        assert c.nK == nK
        assert c.C < c.U
        assert c.descended and c.certified
        assert abs(float(c.Xi) - xi) < 1e-6
        assert abs(float(c.U) - U) < 1e-5
        assert abs(float(c.Xi_star) - xs) < 1e-5
        # elementares Zertifikat
        assert c.n0 * (c.A - 1) >= c.U
        assert c.n0 + c.C == c.A * c.nK
        assert c.A * c.nK < c.n0 * c.A
    print("  CSV/Identitäten OK", flush=True)

    for X in (10**5, 10**6, 10**7):
        print(f"\n=== Census Ξ_* , ungerade n ≤ {X} ===", flush=True)
        out = census_xi_star(X)
        print(
            f"  odd={out['n_odd']}  Ξ<1: {out['n_Xi_lt1']}  "
            f"Ξ*<1: {len(out['n_Xistar_lt1'])}  "
            f"nK≥n0: {len(out['n_not_descended'])}",
            flush=True,
        )
        print(
            f"  min Ξ* = {out['min_Xistar']:.6f} at n={out['argmin_Xistar']}",
            flush=True,
        )
        if out["n_Xistar_lt1"]:
            print("  COUNTEREXAMPLES Ξ*<1:", flush=True)
            for c in out["n_Xistar_lt1"][:10]:
                print("   ", fmt_cert(c), flush=True)
        print("  kleinste Ξ*:", flush=True)
        for val, n0, c in out["smallest"][:12]:
            print(
                f"    n={n0} Ξ*={val:.6f} Ξ={float(c.Xi):.6f} "
                f"K={c.K} A-1={float(c.A-1):.6f} U={float(c.U):.6f} nK={c.nK}",
                flush=True,
            )


if __name__ == "__main__":
    main()
