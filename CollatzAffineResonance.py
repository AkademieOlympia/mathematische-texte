"""
CollatzAffineResonance — affine Collatz-Klassen statt Paritäts-Gaps.

Naive Gap-Resonanz auf Paritätsmustern kollidiert stark (zu grob).
Die relevante Struktur ist die affine Transformation

    T^k(n) = (3^a n + c) / 2^k

mit Invariante

    R(n) = ( 3^a/2^k ,  c/2^k ,  |c|/(2^k (1+3^a)) )

Analogie zur Gap-Resonanz R(δ)=(ρ, Σ): hier trägt (A,B,γ) die
Information der affinen Klasse, nicht das reine Paritätswort.

Claim-boundary: Diagnostik / Atlas; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path


def collatz_step(x: int) -> int:
    if x % 2 == 0:
        return x // 2
    return 3 * x + 1


def affine_orbit(n: int, steps: int) -> tuple[int, int, int, int, int]:
    """
    Elementare Collatz-Schritte; liefert (a, k, c, T^k(n), x_end) mit

        T^k(n) = (3^a * n + c) / 2^k  ∈ ℤ

    Update-Regeln:
      gerade:  k ← k+1
      ungerade: a ← a+1,  c ← 3c + 2^k
    """
    if n <= 0:
        raise ValueError("n must be positive")
    a = 0
    k = 0
    c = 0
    x = n
    for _ in range(steps):
        if x % 2 == 0:
            k += 1
            x //= 2
        else:
            # x ↦ 3x+1 = (3^{a+1} n + (3c + 2^k)) / 2^k
            c = 3 * c + (1 << k)
            a += 1
            x = 3 * x + 1
    # Identität prüfen
    num = (3**a) * n + c
    den = 1 << k
    if num % den != 0:
        raise RuntimeError(f"affine identity failed for n={n}, steps={steps}")
    tval = num // den
    if tval != x:
        raise RuntimeError(f"orbit mismatch n={n}: affine={tval} orbit={x}")
    return a, k, c, tval, x


def resonance(n: int, steps: int = 50) -> tuple[Fraction, Fraction, Fraction, int, int, int]:
    """
    R = (A, B, γ) mit
      A = 3^a / 2^k
      B = c / 2^k
      γ = |c| / (2^k (1 + 3^a)) = |B| / (1 + 3^a)
    """
    a, k, c, _t, _x = affine_orbit(n, steps)
    A = Fraction(3**a, 1 << k)
    B = Fraction(c, 1 << k)
    gamma = Fraction(abs(c), (1 << k) * (1 + 3**a))
    return A, B, gamma, a, k, c


def resonance_key(n: int, steps: int = 50) -> tuple[Fraction, Fraction, Fraction]:
    A, B, gamma, *_ = resonance(n, steps)
    return A, B, gamma


def first_join_step(n: int, m: int, max_steps: int = 20_000) -> int | None:
    """Erste gemeinsame Vorwärtsorbit-Zahl; None falls nicht im Budget."""
    sn: set[int] = set()
    sm: set[int] = set()
    x, y = n, m
    for _ in range(max_steps):
        sn.add(x)
        sm.add(y)
        if x in sm or y in sn:
            return _
        x = collatz_step(x)
        y = collatz_step(y)
    return None


def test_atlas(N: int = 10_000, steps: int = 50) -> dict:
    atlas: dict[tuple[Fraction, Fraction, Fraction], list[int]] = defaultdict(list)
    meta: dict[tuple[Fraction, Fraction, Fraction], tuple[int, int, int]] = {}
    for n in range(1, N):
        A, B, gamma, a, k, c = resonance(n, steps)
        key = (A, B, gamma)
        atlas[key].append(n)
        meta[key] = (a, k, c)

    collisions = sum(1 for v in atlas.values() if len(v) > 1)
    multi = {k: v for k, v in atlas.items() if len(v) > 1}

    # Kollisionen = gleiche Restklasse mod 2^k?
    residue_pure = 0
    residue_mixed = 0
    # γ redundant? (A,B) bestimmt Klasse
    ab_only: dict[tuple[Fraction, Fraction], list[int]] = defaultdict(list)
    for key, vals in atlas.items():
        ab_only[(key[0], key[1])].extend(vals)

    join_examples: list[dict] = []
    for key, vals in multi.items():
        a, k, c = meta[key]
        mods = {n % (1 << k) for n in vals}
        if len(mods) == 1:
            residue_pure += 1
        else:
            residue_mixed += 1
        if len(join_examples) < 6 and len(vals) >= 2:
            j = first_join_step(vals[0], vals[1], max_steps=5000)
            join_examples.append(
                {
                    "n": vals[:8],
                    "a": a,
                    "k": k,
                    "c": c,
                    "mod": vals[0] % (1 << k),
                    "join_step": j,
                }
            )

    return {
        "N": N,
        "steps": steps,
        "classes": len(atlas),
        "collisions": collisions,
        "singleton_classes": len(atlas) - collisions,
        "collision_residue_pure": residue_pure,
        "collision_residue_mixed": residue_mixed,
        "classes_AB_only": len(ab_only),
        "gamma_redundant": len(ab_only) == len(atlas),
        "join_examples": join_examples,
        "sample_classes": [
            {
                "R": [str(x) for x in k],
                "n": v[:10],
                "size": len(v),
            }
            for k, v in list(atlas.items())[:12]
        ],
    }


def test_monotonicity(N: int = 5000, steps: int = 50) -> dict:
    """
    Monotonie-Sonden entlang n und entlang der Bahn.
    """
    # 1) γ(n) vs n: Anteil strikter Abstiege / Anstiege bei n → n+1
    gammas = []
    As = []
    for n in range(1, N):
        A, _B, gamma, *_ = resonance(n, steps)
        gammas.append(gamma)
        As.append(A)

    up = down = flat = 0
    for i in range(len(gammas) - 1):
        if gammas[i + 1] > gammas[i]:
            up += 1
        elif gammas[i + 1] < gammas[i]:
            down += 1
        else:
            flat += 1

    # 2) Entlang einer Bahn: A schrumpft oft (mehr /2 als *3), aber nicht immer
    orbit_shrink = 0
    orbit_checked = 0
    for n in range(2, min(N, 200)):
        A0, *_ = resonance(n, steps)
        x = n
        for _ in range(20):
            x = collatz_step(x)
            if x == 1:
                break
        A1, *_ = resonance(x, steps)
        orbit_checked += 1
        if A1 < A0:
            orbit_shrink += 1

    # 3) Fixschritte: größere a bei gleichem k ⇒ größeres A
    # Stichprobe gleicher k
    by_k: dict[int, list[tuple[int, Fraction]]] = defaultdict(list)
    for n in range(1, N):
        A, _B, _g, a, k, _c = resonance(n, steps)
        by_k[k].append((a, A))

    monotone_k_buckets = 0
    checked_k_buckets = 0
    for k, pairs in by_k.items():
        if len(pairs) < 4:
            continue
        checked_k_buckets += 1
        pairs_sorted = sorted(pairs, key=lambda t: t[0])
        ok = all(
            pairs_sorted[i][1] <= pairs_sorted[i + 1][1]
            for i in range(len(pairs_sorted) - 1)
        )
        if ok:
            monotone_k_buckets += 1

    return {
        "gamma_consecutive": {"up": up, "down": down, "flat": flat},
        "orbit_A_shrink_rate": orbit_shrink / max(orbit_checked, 1),
        "orbit_checked": orbit_checked,
        "A_monotone_in_a_given_k": {
            "buckets_checked": checked_k_buckets,
            "buckets_monotone": monotone_k_buckets,
        },
    }


def test_same_class_structure(N: int = 10_000, steps: int = 50) -> dict:
    """
    Frage 1: Erzeugt dieselbe affine Klasse verschiedene n?
    → ja genau die Kollisionsklassen; charakterisiere sie.
    """
    atlas: dict[tuple[Fraction, Fraction, Fraction], list[int]] = defaultdict(list)
    meta: dict[tuple[Fraction, Fraction, Fraction], tuple[int, int, int]] = {}
    for n in range(1, N):
        A, B, gamma, a, k, c = resonance(n, steps)
        key = (A, B, gamma)
        atlas[key].append(n)
        meta[key] = (a, k, c)

    # Innerhalb einer Klasse: gelten sie dieselbe (a,k,c)?
    # (sollte ja — R bestimmt A,B,γ und damit oft (a,k) bis auf Ambiguität)
    class_stats = []
    for key, vals in atlas.items():
        if len(vals) < 2:
            continue
        a, k, c = meta[key]
        # Prüfe: n ≡ n' mod 2^k / gcd(3^a, 2^k)=2^k  ⇒ gleiche Bahnspitze?
        # Affine: T^k(n)=(3^a n+c)/2^k gleich ⇔ 3^a (n-m) ≡ 0 mod 2^k
        # da 3^a ungerade: n≡m mod 2^k
        mods = {n % (1 << k) for n in vals}
        class_stats.append(
            {
                "a": a,
                "k": k,
                "c": c,
                "size": len(vals),
                "residue_classes_mod_2^k": len(mods),
                "sample": vals[:8],
            }
        )

    # Sortiere größte Kollisionen
    class_stats.sort(key=lambda d: -d["size"])
    return {
        "multi_classes": len(class_stats),
        "largest": class_stats[:10],
        "all_single_residue": all(d["residue_classes_mod_2^k"] == 1 for d in class_stats),
    }


def scan_horizons(N: int, step_list: list[int]) -> list[dict]:
    """Vergleich über Horizonte: Kollisionen erscheinen sobald N ≳ 2^k."""
    rows = []
    for steps in step_list:
        atlas = test_atlas(N, steps)
        struct = test_same_class_structure(N, steps)
        rows.append(
            {
                "steps": steps,
                "classes": atlas["classes"],
                "collisions": atlas["collisions"],
                "residue_pure": atlas["collision_residue_pure"],
                "residue_mixed": atlas["collision_residue_mixed"],
                "gamma_redundant": atlas["gamma_redundant"],
                "all_single_residue": struct["all_single_residue"],
                "largest": struct["largest"][:3],
            }
        )
    return rows


def main() -> int:
    N = 10_000
    # Kurze Horizonte: echte Kollisionen; steps=50: N ≪ 2^50 ⇒ trivial injektiv
    step_list = [5, 8, 10, 12, 16, 20, 50]

    print("=" * 72)
    print("AFFINE COLLATZ-RESONANZ  —  R = (A, B, γ)")
    print("=" * 72)
    print(f"N={N}")
    print("A=3^a/2^k,  B=c/2^k,  γ=|c|/(2^k(1+3^a))")
    print("Hinweis: bei festem k und N < 2^k oft trivial injektiv.")

    rows = scan_horizons(N, step_list)
    print("\n[1] Horizont-Scan")
    print(
        f"  {'steps':>5}  {'Klassen':>8}  {'Koll':>6}  "
        f"{'res-pure':>8}  {'res-mix':>7}  {'γ red.':>6}"
    )
    for r in rows:
        print(
            f"  {r['steps']:5d}  {r['classes']:8d}  {r['collisions']:6d}  "
            f"{r['residue_pure']:8d}  {r['residue_mixed']:7d}  "
            f"{'yes' if r['gamma_redundant'] else 'NO':>6}"
        )

    # Detail: Horizont mit vielen Kollisionen, aber nicht degeneriert
    detail_steps = next(
        (r["steps"] for r in rows if 0 < r["collisions"] < r["classes"]),
        10,
    )
    atlas = test_atlas(N, detail_steps)
    struct = test_same_class_structure(N, detail_steps)
    mono = test_monotonicity(min(N, 5000), detail_steps)

    print(f"\n[2] Detail bei steps={detail_steps}")
    print(f"  Klassen              = {atlas['classes']}")
    print(f"  Kollisionen          = {atlas['collisions']}")
    print(f"  Kollision res-pure   = {atlas['collision_residue_pure']}")
    print(f"  Kollision res-mixed  = {atlas['collision_residue_mixed']}")
    print(f"  γ redundant zu (A,B)? {atlas['gamma_redundant']}")

    print("\n[3] Größte Kollisionsklassen")
    for d in struct["largest"][:6]:
        print(
            f"  size={d['size']:4d}  (a,k,c)=({d['a']},{d['k']},{d['c']})  "
            f"residues={d['residue_classes_mod_2^k']}  sample={d['sample']}"
        )

    print("\n[4] Join-Beispiele (Vorwärtsorbit zweier Kollisionspartner)")
    for ex in atlas["join_examples"]:
        print(
            f"  n={ex['n']}  (a,k)=({ex['a']},{ex['k']})  "
            f"mod={ex['mod']}  join@{ex['join_step']}"
        )

    print("\n[5] Monotonie-Sonden")
    print(f"  γ(n+1) vs γ(n): {mono['gamma_consecutive']}")
    print(f"  A schrumpft entlang 20er-Bahn: {mono['orbit_A_shrink_rate']:.3f}")
    print(f"  A monoton in a | fix k: {mono['A_monotone_in_a_given_k']}")

    print("\n[6] Fazit")
    any_coll = any(r["collisions"] > 0 for r in rows)
    any_mix = any(r["residue_mixed"] > 0 for r in rows)
    if any_coll and not any_mix:
        print(
            "  Paritäts-Gaps waren zu grob; affine R trägt (a,k,c).\n"
            "  R(n)=R(m) ⇔ gleiches (a,k,c) ⇔ n≡m (mod 2^k) im Scan.\n"
            "  γ ist zu (A,B) redundant. steps=50/N=10^4: trivial injektiv."
        )
    elif any_mix:
        print(
            "  Kollisionen mit gemischten Restklassen — unerwartet;\n"
            "  affine Kodierung nicht rein 2-adisch."
        )
    else:
        print("  Keine Kollisionen in den getesteten Horizonten.")

    out = Path(__file__).resolve().parent / "collatz_affine_resonance.json"
    payload = {
        "claim": (
            "Affine Collatz resonance R=(3^a/2^k, c/2^k, |c|/(2^k(1+3^a))); "
            "parity gaps alone are too coarse; not a Collatz proof"
        ),
        "horizons": rows,
        "detail_steps": detail_steps,
        "atlas": atlas,
        "structure": struct,
        "monotonicity": {
            "gamma_consecutive": mono["gamma_consecutive"],
            "orbit_A_shrink_rate": mono["orbit_A_shrink_rate"],
            "A_monotone_in_a_given_k": mono["A_monotone_in_a_given_k"],
        },
        "definition": {
            "A": "3^a / 2^k",
            "B": "c / 2^k",
            "gamma": "|c| / (2^k (1+3^a))",
            "identity": "T^k(n) = (3^a n + c) / 2^k",
        },
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
