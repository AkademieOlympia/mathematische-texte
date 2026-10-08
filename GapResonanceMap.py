"""
Resonanzklassen — Gap-Geometrie als eigenständiges Objekt.

Signatur:
  R(δ) = (ρ(δ), Σ(δ))
  ρ = |M_δ 1| / √k
  Σ = spec(P M_δ P),  P = I − 11ᵀ/k

Resonanzäquivalenz: R(δ1)=R(δ2)
Natürlicher Raum: G_k = R^k / ⟨1⟩
Spiegelung/Translation → gleiche Signatur

Keine Primzahlspezifik.
"""

from __future__ import annotations

import itertools
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def gap_matrix(delta: np.ndarray, eps: float = 1e-30) -> np.ndarray:
    d = np.asarray(delta, dtype=float)
    k = len(d)
    M = np.zeros((k, k), dtype=float)
    for i in range(k):
        for j in range(k):
            if i != j:
                M[i, j] = 1.0 / max(abs(d[i] - d[j]), eps)
    return M


def center(delta: np.ndarray) -> np.ndarray:
    d = np.asarray(delta, dtype=float)
    return d - d.mean()


def signature(
    delta: np.ndarray,
    *,
    rho_digits: int = 12,
    eig_digits: int = 12,
    tol: float = 1e-12,
) -> tuple[float, tuple[float, ...]]:
    """
    R(δ) = (ρ, Σ) mit gerundeten Komponenten für exakte Dict-Keys.
    Σ = alle Eigenwerte von PMP (inkl. ~0 für ⟨1⟩), sortiert absteigend.
    """
    d = np.asarray(delta, dtype=float)
    k = len(d)
    M = gap_matrix(d)
    ones = np.ones(k)
    rho = float(np.linalg.norm(M @ ones) / math.sqrt(k))
    P = np.eye(k) - np.outer(ones, ones) / k
    PMP = 0.5 * (P @ M @ P + (P @ M @ P).T)
    eigs = np.sort(np.linalg.eigvalsh(PMP))[::-1]
    # numerisches Nullen glätten
    eigs = np.where(np.abs(eigs) < tol, 0.0, eigs)
    rho_r = round(rho, rho_digits)
    eig_r = tuple(round(float(x), eig_digits) for x in eigs)
    return rho_r, eig_r


def resonance_type(sig: tuple[float, tuple[float, ...]], tol: float = 1e-12) -> str:
    rho, _ = sig
    return "I_resonant" if abs(rho) > tol else "II_centered"


def G(*offsets: int) -> tuple[int, ...]:
    """G(a1,...,am) → δ=(0,a1,...,am) als Int-Tuple (kanonische Darstellung)."""
    return (0, *offsets)


def normalize_delta_repr(delta: tuple[int, ...] | list[int]) -> tuple[int, ...]:
    """Translation auf min=0, sortiert die absolute Repräsentation nicht um."""
    d = np.asarray(delta, dtype=int)
    return tuple(int(x) for x in (d - d.min()))


def mirror_normalize(delta: tuple[int, ...]) -> tuple[int, ...]:
    """
    Kanonischer Vertreter unter Translation + Spiegelung:
      δ ↦ δ - min(δ), und min( lex(δ), lex(max-δ) ).
    """
    d = np.asarray(delta, dtype=int)
    d = d - d.min()
    m = int(d.max())
    refl = m - d
    refl = refl - refl.min()
    a = tuple(sorted(int(x) for x in d))
    b = tuple(sorted(int(x) for x in refl))
    # speichere Positionsform mit 0 zuerst in Original-Reihenfolge der sortierten unique?
    # Für Katalog: behalte sortierte Positionsmenge als Menge
    return min(a, b)


# ---------------------------------------------------------------------------
# Kataloge
# ---------------------------------------------------------------------------

HAND_CATALOG: list[tuple[str, tuple[int, ...]]] = [
    ("Triplet", G(2, 6)),
    ("Vierling", G(2, 6, 8)),
    ("G(2,4,6,8)", G(2, 4, 6, 8)),
    ("G(2,6,8,12)", G(2, 6, 8, 12)),
    ("G(4,6,10,12)", G(4, 6, 10, 12)),
    ("G(2,4)", G(2, 4)),
    ("G(2,4,6)", G(2, 4, 6)),
    ("G(2,4,6,10)", G(2, 4, 6, 10)),
    ("G(2,6,10,12)", G(2, 6, 10, 12)),
    ("G(6,10,12,18)", G(6, 10, 12, 18)),
    ("G(2,8,12,14)", G(2, 8, 12, 14)),
    ("G(4,10,12,18)", G(4, 10, 12, 18)),
    ("G(2,6,10)", G(2, 6, 10)),
    ("G(0,2,6) raw", (0, 2, 6)),
]


def integer_gap_catalog(k: int, max_coord: int) -> list[tuple[int, ...]]:
    """
    Alle strikt aufsteigenden δ ∈ {0=δ0 < δ1 < … < δ_{k-1} ≤ max_coord}.
    (= admissible Positionsformen mit paarweise verschiedenen Koordinaten)
    """
    if k < 2:
        return []
    out = []
    for combo in itertools.combinations(range(1, max_coord + 1), k - 1):
        out.append((0, *combo))
    return out


@dataclass
class ClassBucket:
    signature: tuple[float, tuple[float, ...]]
    members: list[tuple[str, tuple[int, ...]]]
    type: str


def classify_catalog(
    items: list[tuple[str, tuple[int, ...]]],
) -> dict[tuple[float, tuple[float, ...]], ClassBucket]:
    buckets: dict[tuple[float, tuple[float, ...]], ClassBucket] = {}
    for name, delta in items:
        sig = signature(np.array(delta, dtype=float))
        if sig not in buckets:
            buckets[sig] = ClassBucket(
                signature=sig, members=[], type=resonance_type(sig)
            )
        buckets[sig].members.append((name, delta))
    return buckets


def count_signatures_for_k(k: int, max_coord: int) -> dict:
    """Zählt |{R(δ)}| für alle aufsteigenden Integer-δ in Dimension k."""
    cats = integer_gap_catalog(k, max_coord)
    sigs = set()
    type_I = 0
    type_II = 0
    for delta in cats:
        sig = signature(np.array(delta, dtype=float))
        sigs.add(sig)
        if resonance_type(sig) == "I_resonant":
            type_I += 1
        else:
            type_II += 1
    # nach Spiegelungs-/Translations-Normalform
    mirror_forms = {mirror_normalize(d) for d in cats}
    return {
        "k": k,
        "max_coord": max_coord,
        "n_deltas": len(cats),
        "n_signatures": len(sigs),
        "n_mirror_forms": len(mirror_forms),
        "n_type_I_points": type_I,
        "n_type_II_points": type_II,
        "compression": (
            len(cats) / len(sigs) if sigs else float("nan")
        ),
    }


def main() -> int:
    print("=" * 78)
    print("RESONANZKLASSEN  —  R(δ)=(ρ, Σ),  Äquivalenz R(δ1)=R(δ2)")
    print("=" * 78)
    print(
        "Raum: G_k = R^k/⟨1⟩\n"
        "Invarianzen: Translation δ↦δ+c1, Spiegelung (empirisch/geometrisch)\n"
        "claim-boundary: Spektralklassen der Gap-Geometrie; keine Primheit."
    )

    # --- Handkatalog ---
    print("\n[1] HANDKATALOG — Zusammenführung nach Signatur")
    print("-" * 78)
    buckets = classify_catalog(HAND_CATALOG)
    # sortiere nach ρ absteigend
    ordered = sorted(buckets.items(), key=lambda kv: -kv[0][0])

    for sig, bucket in ordered:
        rho, eigs = sig
        print(f"\nSignature:")
        print(f"  rho = {rho}")
        print(f"  spec = {eigs}")
        print(f"  type = {bucket.type}")
        print(f"Class ({len(bucket.members)} members):")
        for name, delta in bucket.members:
            print(f"  {delta}   # {name}")

    multi = [b for b in buckets.values() if len(b.members) > 1]
    print("\n[2] NICHTTRIVIALE ÄQUIVALENZKLASSEN")
    if not multi:
        print("  (keine Mehrfachklassen im Handkatalog)")
    else:
        for b in multi:
            deltas = [d for _, d in b.members]
            print(f"  rho={b.signature[0]}  →  {deltas}")

    # --- Finiteness-Scan ---
    print("\n[3] SIGNATURZÄHLUNG  |{R(δ)}| für Integer-δ (aufsteigend)")
    print("-" * 78)
    print(
        f"{'k':>3}  {'L':>4}  {'#δ':>8}  {'#R':>8}  "
        f"{'#mirror':>8}  {'#δ/#R':>8}  {'TypII':>6}"
    )
    finiteness_rows = []
    for k in (3, 4, 5):
        for L in (8, 12, 16, 20):
            if L < k:
                continue
            row = count_signatures_for_k(k, L)
            finiteness_rows.append(row)
            print(
                f"{row['k']:3d}  {row['max_coord']:4d}  {row['n_deltas']:8d}  "
                f"{row['n_signatures']:8d}  {row['n_mirror_forms']:8d}  "
                f"{row['compression']:8.3f}  {row['n_type_II_points']:6d}"
            )

    print(
        "\n  Lesart: wächst #R mit L≈unbeschränkt, oder saturiert?\n"
        "  (#mirror = Translations-/Spiegelungs-Normalformen; "
        "#R = Spektralsignaturen)"
    )

    # Detail: k=3, L=20 — Liste aller Signaturen
    print("\n[4] VOLLER SIGNATURATLAS  k=3, L≤12")
    print("-" * 78)
    cats3 = [(f"δ={d}", d) for d in integer_gap_catalog(3, 12)]
    b3 = classify_catalog(cats3)
    print(f"  #δ={len(cats3)}  #Signaturen={len(b3)}")
    for sig, bucket in sorted(b3.items(), key=lambda kv: -kv[0][0]):
        rho, eigs = sig
        # kompakte Mitgliederliste
        mems = [d for _, d in bucket.members]
        show = mems[:6]
        more = f" …(+{len(mems)-6})" if len(mems) > 6 else ""
        print(f"  ρ={rho:.8f}  |class|={len(mems):3d}  e.g. {show}{more}")
        print(f"           Σ={eigs}")

    # Typ-II Suche im Integer-Katalog
    print("\n[5] TYP-II IM INTEGER-KATALOG (ρ=0)")
    found_II = False
    for k, L in ((3, 20), (4, 16), (5, 14)):
        for d in integer_gap_catalog(k, L):
            sig = signature(np.array(d, dtype=float))
            if resonance_type(sig) == "II_centered":
                print(f"  gefunden: k={k} δ={d}  Σ={sig[1]}")
                found_II = True
                break
        if found_II:
            break
    if not found_II:
        print(
            "  Keine Typ-II-Klasse mit paarweise verschiedenen Integer-"
            "Koordinaten gefunden.\n"
            "  (ρ=0 erfordert M1=0 — bei positivem Gap-Graph typischerweise "
            "nur nach expliziter Zentrierung des Operators.)"
        )

    # JSON
    out = Path(__file__).resolve().parent / "gap_resonance_classes.json"
    payload = {
        "claim": (
            "Resonance signature R=(rho,spec(PMP)) classifies gap geometry; "
            "not primes"
        ),
        "hand_catalog_classes": [
            {
                "rho": sig[0],
                "spec": list(sig[1]),
                "type": bucket.type,
                "members": [
                    {"name": n, "delta": list(d)} for n, d in bucket.members
                ],
            }
            for sig, bucket in ordered
        ],
        "finiteness_scan": finiteness_rows,
        "open_question": (
            "Is |{R(δ): δ∈Z^k admissible}| finite for fixed k?"
        ),
    }
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    print("\n[6] VERDICT")
    print(
        "  Resonanzklasse = Faser von R. "
        "Translation/Spiegelung fallen zusammen, wo spektral sichtbar."
    )
    print(
        "  Offene Frage: Endlichkeit von |{R(δ)}| bei festem k — "
        "Scan zeigt Wachstum mit L (noch keine Saturierung)."
    )
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
