"""
CollatzBallotTree — prefix-constrained Baum im Regime A_K = 2^{Q}/3^{R} ≤ 1.

Exakte Prefix-Bedingung:  2^{Q_j} ≤ 3^{R_j}  für alle Präfixe j
(äquivalent zu S_j ≤ a R_j mit a = log2(3/2), ohne Float).

Messgrößen je Tiefe K:
  min ρ(w),  max Null-Lift-Serie,  min θ = ρ / 2^{Q+1}.

Unit-Test: (2,1)^K  ⇒  ρ_K = 2^{3K+1} − 5,  ℓ_K = 7,  θ_K → 1.

Epistemik: Diagnose des dauerhaft expansiven Regimes. Kein Collatz-Beweis.
Selbst Unbeschränktheit von ρ auf Ballot-Ästen schließt oszillierende
A_K ≶ 1 nicht aus.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


def inv_odd(a: int, mod: int) -> int:
    return pow(a, -1, mod)


def v2(n: int) -> int:
    return (n & -n).bit_length() - 1 if n else 0


def actual_block(n: int) -> tuple[int, int, int]:
    r = v2(n + 1)
    x = n
    for _ in range(r):
        x = (3 * x + 1) // 2
    s = v2(x)
    return r, s, x >> s


def prefix_safe(R: int, S: int) -> bool:
    """2^{R+S} ≤ 3^R, ganzzahlig."""
    return (1 << (R + S)) <= 3**R


def block_residue(blocks: list[tuple[int, int]]) -> tuple[int, int]:
    """Kanonischer Startrest ρ mod 2^{Q+1}."""
    Q = sum(r + s for r, s in blocks)
    M = 1 << (Q + 1)
    y = 1
    for r, s in reversed(blocks):
        y = ((1 << r) * inv_odd(pow(3, r, M), M) * (1 + (1 << s) * y) - 1) % M
    return y, M


@dataclass
class Node:
    path: tuple[tuple[int, int], ...]
    R: int
    S: int
    Q: int
    rho: int | None
    modulus: int
    ell: int | None
    zero_lift_run: int

    @property
    def theta(self) -> float | None:
        if self.rho is None:
            return None
        return self.rho / self.modulus


def extend_node(node: Node, r: int, s: int) -> Optional[Node]:
    R_new = node.R + r
    S_new = node.S + s
    Q_new = R_new + S_new
    if not prefix_safe(R_new, S_new):
        return None
    new_path = node.path + ((r, s),)
    rho_new, M_new = block_residue(list(new_path))
    if not node.path:
        ell = None
        zrun = 0
    else:
        rho_old = node.rho
        M_old = 1 << (node.Q + 1)
        if rho_new % M_old != rho_old:
            raise RuntimeError(
                f"cylinder mismatch: rho_new={rho_new} not ≡ {rho_old} mod {M_old}"
            )
        ell = (rho_new - rho_old) // M_old
        if not (0 <= ell < (1 << (r + s))):
            raise RuntimeError(f"lift out of range: ell={ell}, r+s={r+s}")
        zrun = node.zero_lift_run + 1 if ell == 0 else 0
    return Node(
        path=new_path,
        R=R_new,
        S=S_new,
        Q=Q_new,
        rho=rho_new,
        modulus=M_new,
        ell=ell,
        zero_lift_run=zrun,
    )


def unit_test_21(K: int = 8) -> None:
    path: list[tuple[int, int]] = []
    prev: Node | None = None
    for k in range(1, K + 1):
        path.append((2, 1))
        rho, M = block_residue(path)
        expect = (1 << (3 * k + 1)) - 5
        if rho != expect:
            raise AssertionError(f"(2,1)^{k}: rho={rho} ≠ {expect}")
        if M != 1 << (3 * k + 1):
            raise AssertionError(f"(2,1)^{k}: M={M}")
        if not prefix_safe(2 * k, k):
            raise AssertionError(f"(2,1)^{k} violates prefix_safe")
        if prev is not None:
            child = extend_node(prev, 2, 1)
            assert child is not None
            if child.ell != 7:
                raise AssertionError(f"(2,1)^{k}: ell={child.ell} ≠ 7")
            prev = child
        else:
            prev = Node(path=((2, 1),), R=2, S=1, Q=3, rho=rho, modulus=M, ell=None, zero_lift_run=0)
    print(f"UNIT (2,1)^{K}: rho_K=2^{{3K+1}}-5, ell=7  OK", flush=True)


def integer_follows(n0: int, path: tuple[tuple[int, int], ...]) -> bool:
    n = n0
    for r, s in path:
        ar, as_, nxt = actual_block(n)
        if (ar, as_) != (r, s):
            return False
        n = nxt
    return True


def ballot_depth_of_integer(n0: int, max_blocks: int = 40) -> dict:
    """Wie tief folgt der echte U-Orbit von n0 der Ballot-Bedingung?"""
    n = n0
    R = S = 0
    path = []
    for k in range(1, max_blocks + 1):
        r, s, nxt = actual_block(n)
        if not prefix_safe(R + r, S + s):
            return {
                "n0": n0,
                "ballot_blocks": k - 1,
                "first_break": (r, s),
                "path": path,
                "n_at_break": n,
                "n_next": nxt,
                "descended": nxt < n0,
            }
        path.append((r, s))
        R += r
        S += s
        n = nxt
        if n == 1:
            return {
                "n0": n0,
                "ballot_blocks": k,
                "first_break": None,
                "path": path,
                "n_at_break": n,
                "n_next": n,
                "descended": n < n0,
                "hit_one": True,
            }
    return {
        "n0": n0,
        "ballot_blocks": max_blocks,
        "first_break": None,
        "path": path,
        "capped": True,
        "n_at_break": n,
        "descended": n < n0,
    }


def census_prefix_tree(max_K: int, max_r: int, max_s: int) -> list[int]:
    """Zählt Ballot-Knoten je Tiefe ohne Residuen (nur R,S)."""
    from collections import Counter

    state: Counter[tuple[int, int]] = Counter({(0, 0): 1})
    counts = []
    for _K in range(1, max_K + 1):
        nxt: Counter[tuple[int, int]] = Counter()
        for (R, S), c in state.items():
            for r in range(1, max_r + 1):
                for s in range(1, max_s + 1):
                    if prefix_safe(R + r, S + s):
                        nxt[(R + r, S + s)] += c
        state = nxt
        counts.append(sum(state.values()))
    return counts


def search_prefix_tree(max_K: int = 6, max_r: int = 6, max_s: int = 4) -> tuple[list[dict], list[Node]]:
    root = Node(path=(), R=0, S=0, Q=0, rho=None, modulus=2, ell=None, zero_lift_run=0)
    level = [root]
    reports: list[dict] = []

    for K in range(1, max_K + 1):
        nxt: list[Node] = []
        for node in level:
            for r in range(1, max_r + 1):
                for s in range(1, max_s + 1):
                    child = extend_node(node, r, s)
                    if child is not None:
                        nxt.append(child)
        level = nxt
        if not level:
            reports.append({"K": K, "nodes": 0, "ended": True})
            print(f"Baum endet bei K = {K}")
            break

        min_rho_node = min(level, key=lambda x: x.rho)  # type: ignore[arg-type, return-value]
        max_zero_node = max(level, key=lambda x: x.zero_lift_run)
        min_theta_node = min(level, key=lambda x: x.theta)  # type: ignore[arg-type, return-value]
        n_zero = sum(1 for x in level if x.ell == 0)
        # Extremal: kleinstes ρ / 2^{Q+1} und kleinstes ρ relativ zur (2,1)-Benchmark
        bench = (1 << (3 * K + 1)) - 5
        slower = [x for x in level if x.rho < bench]  # type: ignore[operator]

        zero_examples = []
        for x in level:
            if x.ell == 0 and x.rho and x.rho > 1:
                follows = integer_follows(x.rho, x.path)
                zero_examples.append(
                    {
                        "rho": x.rho,
                        "run": x.zero_lift_run,
                        "path": x.path,
                        "follows": follows,
                    }
                )
                if len(zero_examples) >= 8:
                    break

        rep = {
            "K": K,
            "nodes": len(level),
            "min_rho": min_rho_node.rho,
            "min_rho_path": min_rho_node.path,
            "min_rho_Q": min_rho_node.Q,
            "min_theta": min_theta_node.theta,
            "min_theta_path": min_theta_node.path,
            "max_zero_run": max_zero_node.zero_lift_run,
            "max_zero_rho": max_zero_node.rho,
            "max_zero_path": max_zero_node.path,
            "n_zero_lifts": n_zero,
            "n_slower_than_21": len(slower),
            "zero_examples": zero_examples,
        }
        reports.append(rep)
        print(f"\nK = {K}", flush=True)
        print(f"  Knoten: {rep['nodes']}", flush=True)
        print(f"  min rho = {rep['min_rho']}  bei {rep['min_rho_path']}", flush=True)
        print(f"  min theta = {rep['min_theta']:.6f}  bei {rep['min_theta_path']}", flush=True)
        print(
            f"  längste Null-Lift-Serie = {rep['max_zero_run']}  "
            f"rho = {rep['max_zero_rho']}  bei {rep['max_zero_path']}",
            flush=True,
        )
        print(
            f"  ℓ=0 in dieser Tiefe: {n_zero};  ρ kleiner als (2,1)^K: {len(slower)}",
            flush=True,
        )
        if zero_examples:
            print(f"  ℓ=0, ρ>1 Beispiele: {zero_examples[:4]}", flush=True)

    return reports, level  # type: ignore[return-value]


def beam_search(
    max_K: int = 10,
    max_r: int = 8,
    max_s: int = 4,
    beam: int = 3000,
) -> list[dict]:
    """
    Behält je Tiefe die beam kleinsten ρ-Knoten plus alle ℓ=0.
    Untere Schranke an min ρ, falls der Extremalast im Beam bleibt.
    """
    root = Node(path=(), R=0, S=0, Q=0, rho=None, modulus=2, ell=None, zero_lift_run=0)
    level = [root]
    reports: list[dict] = []
    for K in range(1, max_K + 1):
        nxt: list[Node] = []
        for node in level:
            for r in range(1, max_r + 1):
                for s in range(1, max_s + 1):
                    child = extend_node(node, r, s)
                    if child is not None:
                        nxt.append(child)
        if not nxt:
            print(f"Beam-Baum endet bei K={K}", flush=True)
            break
        nxt.sort(key=lambda x: (x.rho, x.theta))
        zeros = [x for x in nxt if x.ell == 0]
        kept = []
        seen = set()
        for x in zeros + nxt:
            key = (x.path, x.rho)
            if key in seen:
                continue
            seen.add(key)
            kept.append(x)
            if len(kept) >= beam and x not in zeros:
                # nach dem Auffüllen mit Null-Lifts stoppen, sobald Beam voll
                if len(kept) >= beam + len(zeros):
                    break
        # hart: alle zeros plus bis zu `beam` kleinsten rho
        small = nxt[:beam]
        merged = {id(x): x for x in zeros + small}
        level = list(merged.values())
        min_rho_node = min(nxt, key=lambda x: x.rho)  # type: ignore[arg-type, return-value]
        min_theta_node = min(nxt, key=lambda x: x.theta)  # type: ignore[arg-type, return-value]
        max_zero_node = max(nxt, key=lambda x: x.zero_lift_run)
        bench = (1 << (3 * K + 1)) - 5
        n_zero = len(zeros)
        n_slow = sum(1 for x in nxt if x.rho is not None and x.rho < bench)
        rep = {
            "K": K,
            "nodes_full_level": len(nxt),
            "nodes_kept": len(level),
            "min_rho": min_rho_node.rho,
            "min_rho_path": min_rho_node.path,
            "min_theta": min_theta_node.theta,
            "min_theta_path": min_theta_node.path,
            "max_zero_run": max_zero_node.zero_lift_run,
            "max_zero_rho": max_zero_node.rho,
            "max_zero_path": max_zero_node.path,
            "n_zero_lifts": n_zero,
            "n_slower_than_21": n_slow,
        }
        reports.append(rep)
        print(
            f"K={K} generated={len(nxt)} kept={len(level)} "
            f"minρ={rep['min_rho']} minθ={rep['min_theta']:.6f} "
            f"max0={rep['max_zero_run']} ρ0={rep['max_zero_rho']} "
            f"ℓ=0:{n_zero} slower21:{n_slow}",
            flush=True,
        )
        print(f"  minρ-Pfad {rep['min_rho_path']}", flush=True)
        if max_zero_node.zero_lift_run:
            print(f"  max0-Pfad {rep['max_zero_path']}", flush=True)
        # Beam nur die gehaltenen Knoten weiterexpandieren
        level = small if not zeros else list({id(x): x for x in zeros + small}.values())
    return reports


def m_of_K(n_max: int = 200_000, max_blocks: int = 80) -> list[tuple[int, int, int]]:
    """m(K) = min { n0 odd > 1 : ballot_depth(n0) ≥ K }."""
    depths: dict[int, int] = {}
    for n0 in range(3, n_max, 2):
        depths[n0] = ballot_depth_of_integer(n0, max_blocks=max_blocks)["ballot_blocks"]
    maxd = max(depths.values()) if depths else 0
    rows = []
    for K in range(0, maxd + 1):
        cands = [n for n, d in depths.items() if d >= K]
        if not cands:
            break
        m = min(cands)
        rows.append((K, m, depths[m]))
    return rows


def main() -> None:
    unit_test_21(8)
    print("L(n)-Census: python3 CollatzL.py", flush=True)


if __name__ == "__main__":
    main()
