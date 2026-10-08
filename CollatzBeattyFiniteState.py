"""
Finite-State Three-Shift — Diagnose der Beatty-Wand r=2, S_K = floor(α K),
α = log2(9/4), exakt über bit_length(3^{2K}).

Kein Collatz-Beweis. Ziel: maschinenprüfbare Fakten für
  (1) L_{K,q} und Sturmian-Faktoren,
  (2) welche Restklassen q ∈ {1,5,6} Beatty-Schritte überleben,
  (3) ob ℓ_K auf der kanonischen Bahn eine Funktion eines endlichen Faktors ist,
  (4) ob ein Null-Lift-Schwanz auftreten kann.

Der kanonische Rest ρ_K realisiert den Präfix der Länge K per Konstruktion.
Deshalb überlebt seine Restklasse jeden endlichen Shift — ein Automat, der
*alle* Zustände bei q∈{1,5,6} töten will, kann die Wand selbst nicht enthalten.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import islice

from CollatzBallotTree import actual_block, block_residue, extend_node, Node, v2


def S_wall(K: int) -> int:
    """floor(α K) = floor(K log2(9/4)), exakt."""
    if K <= 0:
        return 0
    R = 2 * K
    return (3**R).bit_length() - 1 - R


def s_letter(t: int) -> int:
    """s_t ∈ {1,2}, t ≥ 1."""
    return S_wall(t) - S_wall(t - 1)


def L_inc(K: int, q: int) -> int:
    """Q-Zuwachs über q Wandblöcke ab Index K (nach K Blöcken, nächste q)."""
    return 2 * q + S_wall(K + q) - S_wall(K)


def window(K: int, q: int) -> tuple[int, ...]:
    return tuple(s_letter(K + 1 + i) for i in range(q))


def factors_of_length(n: int, span: int = 4000) -> list[tuple[int, ...]]:
    seen: dict[tuple[int, ...], int] = {}
    for K in range(0, span):
        w = window(K, n)
        if w not in seen:
            seen[w] = K
    return sorted(seen, key=lambda w: seen[w])


def follow_blocks(n: int, letters: tuple[int, ...]) -> bool:
    cur = n
    for s in letters:
        r, ss, nxt = actual_block(cur)
        if r != 2 or ss != s:
            return False
        cur = nxt
    return True


def survivors_mod(letters: tuple[int, ...], bits: int) -> list[int]:
    """Ungerade n < 2^bits, die das Wortstück (2,s) exactly folgen."""
    mod = 1 << bits
    out = []
    for n in range(1, mod, 2):
        if follow_blocks(n, letters):
            out.append(n)
    return out


def composed_affine(letters: tuple[int, ...]) -> tuple[int, int]:
    """2^L n_q = 9^q n + 5 * sum_{m<q} 9^{q-1-m} 2^{L_m}."""
    L = 0
    Ssum = 0  # the sum with 5 later
    # rebuild: after m steps, 2^{L_m} n_m = 9^m n + 5 * Σ_{j<m} 9^{m-1-j} 2^{L_j}
    # S_m = Σ_{j<m} 9^{m-1-j} 2^{L_j}
    # recurrence: S_0=0; 2^{L_{m+1}} n' = 9 n_m + 5, L_{m+1}=L_m+2+s
    # S_{m+1} = 9 S_m + 2^{L_m}  with wait:
    # 2^{L+d} n' = 9 (9^m n + 5 S_m) / 2^L * 2^L  + 5 * 2^L? 
    # 2^{L+d} n' = 9 * (9^m n + 5 S_m) + 5 * 2^L = 9^{m+1} n + 5 (9 S_m + 2^L)
    q = len(letters)
    Ssum = 0
    for s in letters:
        Ssum = 9 * Ssum + (1 << L)
        L += 2 + s
    return L, 5 * Ssum


def canonical_lifts(Kmax: int) -> list[dict]:
    rows = []
    prev: Node | None = None
    path: list[tuple[int, int]] = []
    for K in range(1, Kmax + 1):
        s = s_letter(K)
        path.append((2, s))
        rho, M = block_residue(path)
        if prev is None:
            node = Node(
                path=((2, s),),
                R=2,
                S=s,
                Q=2 + s,
                rho=rho,
                modulus=M,
                ell=None,
                zero_lift_run=0,
            )
            ell = None
        else:
            node = extend_node(prev, 2, s)
            assert node is not None
            ell = node.ell
        rows.append(
            {
                "K": K,
                "s": s,
                "L1": 2 + s,
                "ell": ell,
                "rho": node.rho,
                "Q": node.Q,
                "zrun": node.zero_lift_run,
                "window6": window(max(K - 6, 0), 6) if K >= 6 else None,
            }
        )
        prev = node
    return rows


def ell_by_factor(rows: list[dict], n: int) -> dict[tuple[int, ...], set[int]]:
    """ℓ_{K+1} gruppiert nach s-Faktor der Länge n bei K (Start der Extension)."""
    groups: dict[tuple[int, ...], set[int]] = defaultdict(set)
    for i in range(len(rows) - 1):
        K = rows[i]["K"]  # after K blocks; next lift is ell at K+1
        w = window(K, n)
        ell = rows[i + 1]["ell"]
        if ell is not None:
            groups[w].add(ell)
    return dict(groups)


def q_shift_report(q: int, span: int = 400) -> dict:
    facs = factors_of_length(q, span=span)
    bits_needed = []
    surviving = []
    for w in facs:
        L, _ = composed_affine(w)
        bits = L + 1
        bits_needed.append(bits)
        surv = survivors_mod(w, bits)
        surviving.append((w, L, bits, len(surv), surv[:8]))
    return {
        "q": q,
        "n_factors": len(facs),
        "sturmian_bound": q + 1,
        "bits": bits_needed,
        "rows": surviving,
    }


def main() -> None:
    print("=== L_{K,q} exakt, K < 2000 ===", flush=True)
    for q in (1, 4, 5, 6, 51):
        vals = [L_inc(K, q) for K in range(0, 2000 - q)]
        print(f"  q={q:2d}  L in {min(vals)}..{max(vals)}", flush=True)

    print("\n=== Sturmian-Faktoren ===", flush=True)
    for n in (1, 2, 5, 6, 12, 52):
        facs = factors_of_length(n)
        print(f"  Länge {n:2d}: {len(facs)} Faktoren (Morse–Hedlund {n+1})", flush=True)

    print("\n=== Affine 2^L n' = 9^q n + C, Check q=1,6 ===", flush=True)
    for q in (1, 6):
        w = window(0, q)
        L, C = composed_affine(w)
        # n= ρ after 0 blocks: start odd. Use a number that follows the window if any
        bits = L + 1
        surv = survivors_mod(w, bits)
        print(f"  q={q} window={w} L={L} C={C} survivors_mod_2^{bits}={len(surv)}", flush=True)
        if surv:
            n0 = surv[0]
            n = n0
            for s in w:
                r, ss, n = actual_block(n)
                assert r == 2 and ss == s
            # n is n_q; check  (1<<L)*n == 9^q * n0 + C
            lhs = (1 << L) * n
            rhs = 9 ** q * n0 + C
            print(f"    identity on n0={n0}: {lhs == rhs}  lhs={lhs} rhs={rhs}", flush=True)

    print("\n=== q-Shift: Restklassen, die das Wortstück folgen ===", flush=True)
    for q in (1, 5, 6):
        rep = q_shift_report(q, span=800)
        print(
            f"  q={q}: {rep['n_factors']} Faktoren (Soll {rep['sturmian_bound']}), "
            f"bits={min(rep['bits'])}..{max(rep['bits'])}",
            flush=True,
        )
        for w, L, bits, nsurv, sample in rep["rows"]:
            print(
                f"    {w}  L={L}  M={bits}  survivors={nsurv}/{1 << (bits - 1)}  e.g. {sample[:4]}",
                flush=True,
            )

    print("\n=== Kanonische Bahn: ℓ und Faktor → ℓ ===", flush=True)
    rows = canonical_lifts(120)
    zeros = [r["K"] for r in rows if r["ell"] == 0]
    max_z = max(r["zrun"] for r in rows)
    print(
        f"  K≤120  #ℓ=0: {len(zeros)} at {zeros}  max zero-run={max_z}",
        flush=True,
    )
    for n in (1, 2, 6, 12, 24):
        g = ell_by_factor(rows, n)
        n_unique = sum(1 for v in g.values() if len(v) == 1)
        n_amb = sum(1 for v in g.values() if len(v) > 1)
        amb_ex = [(w, sorted(v)) for w, v in g.items() if len(v) > 1][:3]
        print(
            f"  Faktorlänge {n:2d}: {len(g)} Fenster, eindeutig ℓ: {n_unique}, "
            f"mehrdeutig: {n_amb}  z.B. {amb_ex}",
            flush=True,
        )

    print("\n=== Kanonisches ρ folgt dem Präfix (Sanity) ===", flush=True)
    path = [(2, s_letter(t)) for t in range(1, 21)]
    rho, _M = block_residue(path)
    n = rho
    ok = True
    for r_exp, s_exp in path:
        r, s, n = actual_block(n)
        if (r, s) != (r_exp, s_exp):
            ok = False
            break
    print(f"  ρ_20={rho} folgt 20 Wandblöcke: {ok}", flush=True)
    print(
        "  ⇒ die Restklasse von ρ_K überlebt jeden endlichen Three-Shift.",
        flush=True,
    )


if __name__ == "__main__":
    main()
