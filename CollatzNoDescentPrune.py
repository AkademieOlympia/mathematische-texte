"""
CollatzNoDescentPrune — affine Zylinder, Residuen rho, Branch-and-Bound.

Block (r, s): r ungerade Syracuse-Schritte, danach s extra Halbierungen
(Odd-to-Odd-Operator U aus Collatz reloaded):

    U(n) = (3^r (n+1) - 2^r) / 2^{r+s}

Affine Identität nach j Blöcken:

    2^{Q_j} n_j = 3^{R_j} n_0 + B_j,    D_j = 2^{Q_j} - 3^{R_j}.

No-Descent n_j >= n_0:
    D_j > 0  =>  n_0 <= B_j / D_j
    D_j < 0  =>  automatisch für alle n_0 >= 1 (B_j > 0).

Kanonischer Zylinder: exaktes (r, s) legt n_0 mod 2^{Q+1} fest.
Ist D>0 und rho * D > B, ist der Zylinder leer für Integer-No-Descent.

Epistemik: deterministische Diagnose. Kein Collatz-Beweis.
Der expandierende Zweig D<0 (z. B. konstantes (2,1)) wird nie gepruned.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass


def inv_odd(a: int, mod: int) -> int:
    return pow(a, -1, mod)


def apply_block(n: int, r: int, s: int) -> int | None:
    """U-Block, oder None wenn (r,s) nicht exakt passt."""
    if n % 2 == 0 or n <= 0:
        return None
    if ((n + 1) & - (n + 1)).bit_length() - 1 != r:  # v2(n+1)
        return None
    x = n
    for _ in range(r):
        if x % 2 == 0:
            return None
        x = (3 * x + 1) // 2
    v = (x & -x).bit_length() - 1 if x else 0
    if v != s:
        return None
    return x >> s


def v2(n: int) -> int:
    return (n & -n).bit_length() - 1 if n else 0


def actual_block(n: int) -> tuple[int, int, int]:
    r = v2(n + 1)
    x = n
    for _ in range(r):
        x = (3 * x + 1) // 2
    s = v2(x)
    return r, s, x >> s


@dataclass(frozen=True)
class PrefixState:
    step: int
    block: tuple[int, int]
    R: int
    Q: int
    D: int
    B: int
    rho: int
    modulus: int
    pruned: bool
    lift: int | None


def update_B(B: int, Q: int, r: int) -> int:
    return 3**r * B + (3**r - 2**r) * (1 << Q)


def rho_from_prefix(blocks: list[tuple[int, int]]) -> tuple[int, int]:
    """Kanonisches Residuum mod 2^{Q+1} durch inverse Blockkette von y=1."""
    Q = sum(r + s for r, s in blocks)
    M = 1 << (Q + 1)
    y = 1
    for r, s in reversed(blocks):
        y = ((1 << r) * inv_odd(pow(3, r, M), M) * (1 + (1 << s) * y) - 1) % M
    return y, M


def analyze_block_path(blocks: list[tuple[int, int]]) -> list[PrefixState]:
    R = Q = B = 0
    prefix: list[tuple[int, int]] = []
    results: list[PrefixState] = []
    prev_rho, prev_M = 0, 1

    for idx, (r, s) in enumerate(blocks):
        prefix.append((r, s))
        B = update_B(B, Q, r)
        R += r
        Q += r + s
        D = (1 << Q) - 3**R
        rho, M = rho_from_prefix(prefix)
        pruned = D > 0 and rho * D > B
        lift = None
        if idx > 0 and M % prev_M == 0:
            # rho_{K+1} = rho_K + ell * 2^{Q_K+1}, hier prev_M = 2^{Q_prev+1}
            delta = rho - prev_rho
            if delta % prev_M == 0:
                lift = delta // prev_M
        results.append(
            PrefixState(
                step=idx + 1,
                block=(r, s),
                R=R,
                Q=Q,
                D=D,
                B=B,
                rho=rho,
                modulus=M,
                pruned=pruned,
                lift=lift,
            )
        )
        prev_rho, prev_M = rho, M
        if pruned:
            break
    return results


def verify_affine_on_integer(n0: int, max_blocks: int = 12) -> list[dict]:
    """Prüft 2^Q n_j = 3^R n0 + B längs des echten U-Orbits von n0."""
    n = n0
    blocks = []
    rows = []
    for _ in range(max_blocks):
        r, s, nxt = actual_block(n)
        blocks.append((r, s))
        st = analyze_block_path(blocks)[-1]
        lhs = (1 << st.Q) * nxt
        rhs = 3**st.R * n0 + st.B
        rows.append(
            {
                "n": n,
                "block": (r, s),
                "n_next": nxt,
                "R": st.R,
                "Q": st.Q,
                "D": st.D,
                "B": st.B,
                "rho": st.rho,
                "modulus": st.modulus,
                "n0_mod": n0 % st.modulus,
                "affine_ok": lhs == rhs,
                "in_cylinder": n0 % st.modulus == st.rho,
                "no_descent": nxt >= n0,
                "pruned": st.pruned,
                "lift": st.lift,
            }
        )
        n = nxt
        if n == 1:
            break
    return rows


def bfs_no_descent(
    max_depth: int = 6,
    r_max: int = 4,
    s_max: int = 4,
) -> dict:
    """
    Branch-and-Bound auf Blockwörtern.
    Expandiert nur ungeprunte Präfixe. Zählt D<0 (nie prune) vs D>0.
    """
    alphabet = [(r, s) for r in range(1, r_max + 1) for s in range(1, s_max + 1)]
    q: deque[list[tuple[int, int]]] = deque([[]])
    visited = 0
    pruned = 0
    surviving_leaves = 0
    expanding = 0  # last D < 0
    contracting_alive = 0  # D>0 but not pruned
    examples_expanding: list[list[tuple[int, int]]] = []
    examples_pruned: list[list[tuple[int, int]]] = []
    examples_alive: list[list[tuple[int, int]]] = []

    while q:
        prefix = q.popleft()
        if len(prefix) >= max_depth:
            if prefix:
                st = analyze_block_path(prefix)[-1]
                surviving_leaves += 1
                if st.D < 0:
                    expanding += 1
                    if len(examples_expanding) < 8:
                        examples_expanding.append(prefix)
                else:
                    contracting_alive += 1
                    if len(examples_alive) < 8:
                        examples_alive.append(prefix)
            continue
        for blk in alphabet:
            word = prefix + [blk]
            visited += 1
            st = analyze_block_path(word)[-1]
            if st.pruned:
                pruned += 1
                if len(examples_pruned) < 8:
                    examples_pruned.append(word)
                continue
            q.append(word)

    return {
        "max_depth": max_depth,
        "alphabet": len(alphabet),
        "visited_prefixes": visited,
        "pruned": pruned,
        "surviving_leaves": surviving_leaves,
        "expanding_leaves_Dneg": expanding,
        "contracting_alive_leaves": contracting_alive,
        "prune_ratio": pruned / visited if visited else None,
        "examples_expanding": examples_expanding,
        "examples_pruned": examples_pruned,
        "examples_alive_contracting": examples_alive,
    }


def fmt_states(states: list[PrefixState]) -> str:
    lines = []
    for st in states:
        lines.append(
            f"  j={st.step} block={st.block} R={st.R} Q={st.Q} D={st.D} "
            f"B={st.B} rho={st.rho} mod={st.modulus} lift={st.lift} pruned={st.pruned}"
        )
    return "\n".join(lines)


def main() -> None:
    print("=== Affine-Check n0=27 ===")
    for row in verify_affine_on_integer(27, max_blocks=8):
        print(
            f"  {row['n']} -{row['block']}-> {row['n_next']} | "
            f"affine={row['affine_ok']} cyl={row['in_cylinder']} "
            f"n0%M={row['n0_mod']} rho={row['rho']} D={row['D']} "
            f"nd={row['no_descent']} pruned={row['pruned']} lift={row['lift']}"
        )

    print("\n=== Konstantes Wort (2,1)^8 ===")
    path21 = [(2, 1)] * 8
    states21 = analyze_block_path(path21)
    print(fmt_states(states21))
    print("  D<0 für alle j?", all(st.D < 0 for st in states21))
    print("  je gepruned?", any(st.pruned for st in states21))
    print("  rho vs 2^{3j+1}-5 (Fixpunkt -5):")
    for st in states21:
        fp = (1 << (3 * st.step + 1)) - 5
        print(f"    j={st.step}: rho={st.rho}  2^{{{3*st.step+1}}}-5={fp}  equal={st.rho==fp}")

    print("\n=== Integer 123 = 2^7-5, erzwungene (2,1)-Folge bis Bruch ===")
    n = 123
    n0 = 123
    for j in range(1, 8):
        r, s, nxt = actual_block(n)
        print(f"  {n} block=({r},{s}) -> {nxt}  no_desc={nxt>=n0}")
        if (r, s) != (2, 1):
            print("  bricht vom Wort (2,1) weg")
            break
        n = nxt

    print("\n=== BFS No-Descent-Pruning r,s<=3, Tiefe 5 ===")
    report = bfs_no_descent(max_depth=5, r_max=3, s_max=3)
    keep = {
        k: report[k]
        for k in (
            "max_depth",
            "alphabet",
            "visited_prefixes",
            "pruned",
            "surviving_leaves",
            "expanding_leaves_Dneg",
            "contracting_alive_leaves",
            "prune_ratio",
        )
    }
    print(" ", keep)
    print("  expandierende Beispiele:", report["examples_expanding"][:5])
    print("  geprunte Beispiele:", report["examples_pruned"][:5])
    print("  kontrahierend lebend:", report["examples_alive_contracting"][:5])


if __name__ == "__main__":
    main()
