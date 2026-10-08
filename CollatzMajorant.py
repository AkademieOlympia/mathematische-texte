"""
Reine Präfix-Majorante M_K am ersten A_K>1-Austritt.

    M_K := K - sum_{j=0}^{K-2} 2^{-s_j}
         = 1 + sum_{j=1}^{K-1} (1 - 2^{-s_{j-1}})

Im sicheren Präfix A_j ≤ 1 gilt U_K ≤ M_K, und M_K hängt nicht von n_K
und nicht von den tatsächlichen A_j ab, nur von K und den s-Werten der
sicheren Blöcke.

Wegen 1 - 2^{-s} ≥ 1/2 für s ≥ 1 gilt außerdem
    M_K ≥ (K + 1)/2.
Das ist eine von der Blockstruktur unabhängige Unterschranke: wenn
n_0(A_K-1) < (K+1)/2, kann M_K den Startwert nicht zertifizieren.

Ξ_* ≥ 1 ist keine unabhängige Majorante: aus der Identität
    U_K = A_K (n_K + 2^{-s_{K-1}}) - n_0
folgt am Crossing A_K > 1
    Ξ_* ≥ 1  ⟺  n_0 ≥ n_K + 2^{-s_{K-1}}  ⟺  n_K < n_0
(da n_0, n_K ungerade).

Ganzzahliger Ballot-Rest, ohne Float:
    m_j = floor(log2(3^{R_j})) - Q_j = (3^{R_j}).bit_length() - 1 - Q_j.
Auf dem sicheren Präfix ist m_j ≥ 0, und m_j = 0 genau dann, wenn
    2^{Q_j} ≤ 3^{R_j} < 2^{Q_j+1},  also A_j ∈ (1/2, 1].
Es gilt d_j = {a R_j} + m_j und daher
    2^{-{a R_j}} = 2^{m_j} A_j.
M_K^{(2)} mit 2^{-{a R_j}} überschätzt deshalb jeden Term um 2^{m_j}.
Die „verfeinerte“ Summe mit den tatsächlichen m_j und {a R_j} ist identisch U_K.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import log2

from CollatzBallotTree import actual_block
from CollatzXi import first_exit_cert


def m_ballot(R: int, Q: int) -> int:
    """floor(R log2 3) - Q, exakt. m ≥ 0 ⇔ 2^Q ≤ 3^R."""
    return (3 ** R).bit_length() - 1 - Q


@dataclass
class MajorantCert:
    n0: int
    K: int
    L: int
    A: Fraction
    U: Fraction
    M: Fraction
    overflow: Fraction  # n0(A-1)
    nK: int
    s_safe: tuple[int, ...]
    r_safe: tuple[int, ...]
    r_exit: int
    s_exit: int
    R_exit: int
    Q_exit: int

    @property
    def M_floor(self) -> Fraction:
        """Universelle Unterschranke (K+1)/2."""
        return Fraction(self.K + 1, 2)

    @property
    def certified_M(self) -> bool:
        return self.overflow >= self.M

    @property
    def certified_floor(self) -> bool:
        return self.overflow >= self.M_floor

    @property
    def slack_bits(self) -> list[float]:
        """λ_j = R_j log2 3 - Q_j auf dem sicheren Präfix (A_j = 2^{-λ_j})."""
        out: list[float] = []
        R = Q = 0
        ln3 = log2(3)
        for r, s in zip(self.r_safe, self.s_safe):
            R += r
            Q += r + s
            out.append(R * ln3 - Q)
        return out


def first_majorant(n0: int, max_blocks: int = 400) -> MajorantCert:
    c = first_exit_cert(n0, max_blocks=max_blocks)
    n = n0
    r_safe: list[int] = []
    s_safe: list[int] = []
    for _ in range(c.L):
        r, s, nxt = actual_block(n)
        r_safe.append(r)
        s_safe.append(s)
        n = nxt
    r_exit, s_exit, nK = actual_block(n)
    M = Fraction(c.K)
    for sj in s_safe:
        M -= Fraction(1, 1 << sj)
    overflow = n0 * (c.A - 1)
    return MajorantCert(
        n0=n0,
        K=c.K,
        L=c.L,
        A=c.A,
        U=c.U,
        M=M,
        overflow=overflow,
        nK=nK,
        s_safe=tuple(s_safe),
        r_safe=tuple(r_safe),
        r_exit=r_exit,
        s_exit=s_exit,
        R_exit=c.R,
        Q_exit=c.Q,
    )


def _certifies_int(n0: int, Q: int, R: int, s_safe: list[int], K: int) -> bool:
    """n0(A-1) ≥ M, ganzzahlig."""
    if not s_safe:
        # M = 1, A = 2^Q / 3^R
        return n0 * ((1 << Q) - 3**R) >= 3**R
    smax = max(s_safe)
    den = 1 << smax
    Mnum = K * den - sum(den >> sj for sj in s_safe)
    return n0 * ((1 << Q) - 3**R) * den >= 3**R * Mnum


def census_majorant(X: int) -> dict:
    """Ungerade n in (1, X]. Fast path, exakte Integer-Arithmetik."""
    fails: list[int] = []
    fails_floor: list[int] = []
    n_odd = 0
    for n0 in range(3, X + 1, 2):
        n = n0
        R = Q = 0
        s_safe: list[int] = []
        for k in range(1, 401):
            r = (n + 1 & -(n + 1)).bit_length() - 1
            x = n
            for _ in range(r):
                x = (3 * x + 1) // 2
            s = (x & -x).bit_length() - 1
            nxt = x >> s
            R += r
            Q += r + s
            if (1 << Q) > 3**R:
                n_odd += 1
                if not _certifies_int(n0, Q, R, s_safe, k):
                    fails.append(n0)
                half_num = k + 1
                # overflow >= (K+1)/2  iff  2 n0 (2^Q - 3^R) >= (K+1) 3^R
                if 2 * n0 * ((1 << Q) - 3**R) < half_num * (3**R):
                    fails_floor.append(n0)
                break
            s_safe.append(s)
            n = nxt
        else:
            raise RuntimeError(f"no exit n0={n0}")
        if n0 % 2_000_000 == 1:
            print(f"    … n0={n0} M-fails={fails} floor-fails={fails_floor}", flush=True)
    return {
        "X": X,
        "n_odd": n_odd,
        "M_fails": fails,
        "floor_fails": fails_floor,
    }


def fmt(m: MajorantCert) -> str:
    return (
        f"n0={m.n0} K={m.K} n0(A-1)={float(m.overflow):.6f} "
        f"U={float(m.U):.6f} M={float(m.M):.6f} "
        f"(K+1)/2={float(m.M_floor):.3f} "
        f"M-cert={m.certified_M} floor-cert={m.certified_floor} "
        f"s={m.s_safe} exit=({m.r_exit},{m.s_exit}) nK={m.nK}"
    )


def main() -> None:
    print("=== Tautologie Ξ_* ⟺ Abstieg ===", flush=True)
    for n0 in (3, 7, 9, 27, 31, 63, 71, 95):
        c = first_exit_cert(n0)
        taut = (c.n0 >= c.nK + Fraction(1, 1 << c.s_exit)) == (c.nK < c.n0)
        equiv = (c.Xi_star >= 1) == (c.nK < c.n0)
        print(
            f"  n0={n0} Ξ*≥1={c.Xi_star >= 1} nK<n0={c.nK < c.n0} "
            f"equiv={equiv} taut-lemma={taut}",
            flush=True,
        )
        assert equiv and taut

    print("\n=== Majorante M_K, kritische Starts ===", flush=True)
    for n0 in (27, 31, 63, 71, 7, 9):
        m = first_majorant(n0)
        print(" ", fmt(m), flush=True)
        sl = m.slack_bits
        if sl:
            print(
                f"    min λ={min(sl):.4f} bit  max λ={max(sl):.4f} bit  "
                f"λ_last={sl[-1]:.4f}  s-counts={ {s: m.s_safe.count(s) for s in sorted(set(m.s_safe)) }}",
                flush=True,
            )
        # M^{(2)} mit tatsächlichem λ (Diagnose, noch keine unabhängige δ_j)
        M2 = Fraction(1)
        R = Q = 0
        s_prev = None
        n = n0
        for j, (r, s) in enumerate(zip(m.r_safe, m.s_safe), start=1):
            R += r
            Q += r + s
            A = Fraction(1 << Q, 3**R)
            if s_prev is None:
                pass
            else:
                pass
            if j >= 1:
                # U-Formel: erster Term j=1 hat Faktor (1-2^{-s_0}) A_1, s_0 = s of block 1? 
                # Code: U starts 1, then for k>=2 adds (1-2^{-s_prev})*A_prev
                # so first contribution is (1-2^{-s_1}) A_1 after second block...
                # U_K = 1 + sum_{j=1}^{K-1} (1-2^{-s_{j-1}}) A_j
                # For j=1, s_{0} is s of the block that produced A_1, i.e. first block.
                # Index in paper: s_{j-1} for term A_j. First term j=1 uses s_0 = s of block 1.
                # After first block, A_1 exists; U still 1 until we add at k=2.
                # Term for A_1 uses s of block 1. Yes s_safe[0] with A after block 1.
                pass
        # rebuild M2 = 1 + sum (1-2^{-s_{j-1}}) * 2^{-λ_j} = 1 + sum (1-2^{-s}) A_j = U
        print(f"    exact U vs M gap: M-U={float(m.M - m.U):.6f}", flush=True)
        print(
            f"    need M2 ≤ overflow={float(m.overflow):.6f}; "
            f"exact U={float(m.U):.6f} already {float(m.overflow/m.U):.3f}× darunter",
            flush=True,
        )

    for X in (10**5, 10**6, 10**7):
        print(f"\n=== Census M_K, ungerade n ≤ {X} ===", flush=True)
        out = census_majorant(X)
        print(
            f"  odd={out['n_odd']}  M-fails={out['M_fails']}  "
            f"floor-fails (n0(A-1)<(K+1)/2)={out['floor_fails']}",
            flush=True,
        )


if __name__ == "__main__":
    main()
