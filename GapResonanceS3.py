"""
Test S3 — Spektrale Trennung von Triplet-Gap-Klassen (k=3).

Vermutung:
  R(δ1)=R(δ2)  ⟺  δ2=δ1 oder δ2=S(δ1)
  mit S(0,a,b)=(0,b−a,b).

Zusätzlich Vermutung 2:
  δ_m=(0,1,m) hat streng monotones ρ_m ⇒ unendlich viele Resonanzklassen.

Claim-boundary: Gap-Geometrie; keine Primzahlspezifik.
"""

from __future__ import annotations

import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


def gap_matrix(delta: np.ndarray) -> np.ndarray:
    d = np.asarray(delta, dtype=float)
    k = len(d)
    M = np.zeros((k, k), dtype=float)
    for i in range(k):
        for j in range(k):
            if i != j:
                M[i, j] = 1.0 / abs(d[i] - d[j])
    return M


def raw_signature(delta) -> tuple[float, np.ndarray]:
    d = np.asarray(delta, dtype=float)
    k = 3
    M = gap_matrix(d)
    ones = np.ones(k)
    rho = float(np.linalg.norm(M @ ones) / math.sqrt(k))
    P = np.eye(k) - np.outer(ones, ones) / k
    PMP = 0.5 * (P @ M @ P + (P @ M @ P).T)
    eigs = np.sort(np.linalg.eigvalsh(PMP))[::-1]
    return rho, eigs


def sig_close(s1, s2, atol: float = 1e-9) -> bool:
    r1, e1 = s1
    r2, e2 = s2
    return abs(r1 - r2) < atol and bool(np.allclose(e1, e2, atol=atol, rtol=0.0))


def mirror_ab(a: int, b: int) -> tuple[int, int]:
    return b - a, b


def run_s3(bmax: int = 1000, atol: float = 1e-9) -> dict:
    print("=" * 72)
    print(f"TEST S3  —  k=3, 0<a<b≤{bmax},  atol={atol}")
    print("=" * 72)

    # (i) Spiegelinvarianz: R(δ) ≈ R(Sδ)
    n = 0
    mirror_fails = 0
    for a, b in itertools.combinations(range(1, bmax + 1), 2):
        n += 1
        sa, sb = mirror_ab(a, b)
        if not sig_close(raw_signature((0, a, b)), raw_signature((0, sa, sb)), atol):
            mirror_fails += 1
            if mirror_fails <= 3:
                print(f"  mirror fail: ({a},{b}) vs ({sa},{sb})")

    print(f"[i] Spiegelinvarianz: fails={mirror_fails}/{n}")

    # (ii) Trennung: verschiedene Spiegelorbits ⇒ nicht allclose
    # Kanonische Vertreter: lexmin(δ, Sδ)
    reps = []
    for a, b in itertools.combinations(range(1, bmax + 1), 2):
        sa, sb = mirror_ab(a, b)
        if (a, b) <= (sa, sb):
            reps.append((a, b, raw_signature((0, a, b))))

    # Bucket nach gerundetem ρ, dann allclose zwischen Orbits
    buckets: dict[float, list] = defaultdict(list)
    for a, b, s in reps:
        buckets[round(s[0], 6)].append((a, b, s))

    collisions = 0
    collision_examples = []
    for lst in buckets.values():
        for i in range(len(lst)):
            for j in range(i + 1, len(lst)):
                if sig_close(lst[i][2], lst[j][2], atol=atol):
                    collisions += 1
                    if len(collision_examples) < 5:
                        collision_examples.append((lst[i][:2], lst[j][:2]))

    print(f"[ii] Orbit-Kollisionen: {collisions}")
    if collision_examples:
        for ex in collision_examples:
            print(f"     {ex}")

    # (iii) Strukturzähler mit Rundung (Info)
    round_buckets: dict = defaultdict(list)
    for a, b in itertools.combinations(range(1, bmax + 1), 2):
        rho, eigs = raw_signature((0, a, b))
        key = (round(rho, 9), tuple(round(float(x), 9) for x in eigs))
        round_buckets[key].append((a, b))
    pairs = selfm = frag = 0
    for mems in round_buckets.values():
        S = sorted(set(mems))
        if len(S) == 1:
            a, b = S[0]
            if mirror_ab(a, b) == (a, b):
                selfm += 1
            else:
                frag += 1  # Rundungsfragmentierung
        elif len(S) == 2:
            (a1, b1), (a2, b2) = S
            if mirror_ab(a1, b1) == (a2, b2):
                pairs += 1
            else:
                frag += 1
        else:
            frag += 1

    n_orbits = len(reps)
    print(
        f"[iii] Orbits={n_orbits}  (erwartet: self={bmax//2} + pairs="
        f"{(n - bmax//2)//2})"
    )
    print(
        f"      round9-Buckets: pairs={pairs} self={selfm} "
        f"frag_or_odd={frag}  #R_round9={len(round_buckets)}"
    )

    passed = mirror_fails == 0 and collisions == 0
    print()
    print(f"VERDICT S3: {'PASS' if passed else 'FAIL'}")
    if passed:
        print(
            "  R klassifiziert Triplet-Gap-Geometrien bis auf Spiegelung "
            f"(geprüft für b≤{bmax}, allclose atol={atol})."
        )
    return {
        "bmax": bmax,
        "atol": atol,
        "n_deltas": n,
        "n_orbits": n_orbits,
        "mirror_fails": mirror_fails,
        "orbit_collisions": collisions,
        "pass": passed,
    }


def run_infinity_family(m_max: int = 500) -> dict:
    print()
    print("=" * 72)
    print(f"TEST ∞  —  Familie δ_m=(0,1,m), m=2..{m_max}")
    print("=" * 72)
    rhos = []
    for m in range(2, m_max + 1):
        rho, _ = raw_signature((0, 1, m))
        rhos.append(rho)
    mono_dec = all(rhos[i] > rhos[i + 1] for i in range(len(rhos) - 1))
    unique = len(set(round(r, 12) for r in rhos)) == len(rhos)
    print(f"  streng monoton fallend: {mono_dec}")
    print(f"  alle ρ verschieden (12 Digits): {unique}")
    print(f"  ρ(2)={rhos[0]:.12f}  ρ({m_max})={rhos[-1]:.12f}")
    passed = mono_dec and unique
    print(f"VERDICT ∞: {'PASS' if passed else 'FAIL'}")
    if passed:
        print(
            "  ⇒ unendlich viele Resonanzklassen "
            "(injektive Familie δ_m=(0,1,m))."
        )
    return {
        "m_max": m_max,
        "mono_decreasing": mono_dec,
        "all_unique": unique,
        "rho_2": rhos[0],
        "rho_mmax": rhos[-1],
        "pass": passed,
    }


def main() -> int:
    s3 = run_s3(1000)
    inf = run_infinity_family(500)
    out = Path(__file__).resolve().parent / "gap_resonance_s3.json"
    payload = {
        "claim": (
            "Resonance signature classifies triplet gap geometries up to "
            "reflection; infinitely many classes via (0,1,m); not about primes"
        ),
        "S3": s3,
        "infinity_family": inf,
        "conjecture_1": (
            "R(δ1)=R(δ2) iff δ2=δ1 or δ2=S(δ1) for δ=(0,a,b)"
        ),
        "conjecture_2": "|SpecGap_k|=∞ for each fixed k≥3",
    }
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0 if s3["pass"] and inf["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
