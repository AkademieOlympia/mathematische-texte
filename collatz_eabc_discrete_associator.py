#!/usr/bin/env python3
"""
Diskreter EABC-Assoziator auf V₄ = {E,A,B,C} — Negativkontrolle.

Kanonsiche Theorie: collatz_eabc_discrete_associator.md
                   collatz_eabc_holonomie.md  (Korrektur: Holonomie ≠ V₄-𝔞)
                   collatz_eabc_kommutator_assoziator.md

Φ(X,Y) aus eabc_from_lean / EABC.lean:
  Restklassen {1,5,7,11} mod 12, Φ(X,Y) = classOf(residue(X)·residue(Y)).
  V₄ ≅ Klein-Vierergruppe → Φ assoziativ → naive 𝔞 ≡ 0 (Theorem).

Echter EABC-Defekt: Γ((xy)z) - Γ(x(yz)) auf Trägern (Oktanion, Vierlinge ω) —
siehe collatz_eabc_holonomie_test.py.

𝔞_naiv(X,Y,Z) = sgn( Φ(Φ(X,Y),Z) vs Φ(X,Φ(Y,Z)) ) ∈ {-1,0,+1}  (stets 0 auf V₄).

Ausführung:
    python3 collatz_eabc_discrete_associator.py
    python3 collatz_eabc_discrete_associator.py --prime-limit 5000
"""

from __future__ import annotations

import argparse
import json
from math import isqrt
from pathlib import Path
from typing import Any, Literal

from eabc_from_lean import (
    Chirality,
    EClass,
    chirality_order,
    class_of,
    is_prime_quadruplet,
    q,
    residue,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_discrete_associator.json"

ALL_CLASSES: tuple[EClass, ...] = (EClass.E, EClass.A, EClass.B, EClass.C)
ABC_CLASSES: tuple[EClass, ...] = (EClass.A, EClass.B, EClass.C)
AssociatorSign = Literal[-1, 0, 1]
BracketPref = Literal["left", "right"]


def phi(x: EClass, y: EClass) -> EClass:
    """Φ: V₄×V₄→V₄ via mod-12-Multiplikation der Lean-Restklassen."""
    prod = (residue(x) * residue(y)) % 12
    out = class_of(prod)
    if out is None:
        raise ValueError(f"Φ({x.value},{y.value}) fällt aus V₄ (prod mod 12 = {prod})")
    return out


def phi_left(x: EClass, y: EClass, z: EClass) -> EClass:
    return phi(phi(x, y), z)


def phi_right(x: EClass, y: EClass, z: EClass) -> EClass:
    return phi(x, phi(y, z))


def associator_sign(x: EClass, y: EClass, z: EClass) -> AssociatorSign:
    """𝔞(X,Y,Z) ∈ {-1,0,+1}; bei V₄-Assoziativität stets 0."""
    left = phi_left(x, y, z)
    right = phi_right(x, y, z)
    if left is right:
        return 0
    # Fallback falls Klassenobjekte divergieren (sollte auf V₄ nicht vorkommen)
    order = {EClass.E: 0, EClass.A: 1, EClass.B: 2, EClass.C: 3}
    if order[left] > order[right]:
        return 1
    return -1


def multiplication_table() -> dict[str, dict[str, str]]:
    return {
        x.value: {y.value: phi(x, y).value for y in ALL_CLASSES}
        for x in ALL_CLASSES
    }


def associator_table(classes: tuple[EClass, ...]) -> dict[str, dict[str, dict[str, int]]]:
    return {
        x.value: {
            y.value: {z.value: associator_sign(x, y, z) for z in classes}
            for y in classes
        }
        for x in classes
    }


def prove_v4_klein_associativity() -> dict[str, Any]:
    """
    Expliziter V₄-Beweis: {1,5,7,11} mod 12 ist Klein-Vierergruppe, daher assoziativ.

    Schritte:
      1. Abschluss unter Multiplikation mod 12
      2. E=1 neutral
      3. Jedes Nicht-E-Element ist selbstinvers (a·a=E)
      4. Assoziativität: alle 4³ Tripel (oder Gruppentheorem)
    """
    residues = {c: residue(c) for c in ALL_CLASSES}
    closure_ok = True
    for x in ALL_CLASSES:
        for y in ALL_CLASSES:
            prod = (residues[x] * residues[y]) % 12
            if class_of(prod) is None:
                closure_ok = False

    identity_ok = all(phi(EClass.E, c) is c and phi(c, EClass.E) is c for c in ALL_CLASSES)
    inv_ok = all(phi(c, c) is EClass.E for c in ABC_CLASSES)

    # Klein-Isomorphie: A↦(1,0), B↦(0,1), C↦(1,1) in (Z/2)²
    z2_map = {EClass.E: (0, 0), EClass.A: (1, 0), EClass.B: (0, 1), EClass.C: (1, 1)}

    def z2_add(u: tuple[int, int], v: tuple[int, int]) -> tuple[int, int]:
        return ((u[0] + v[0]) % 2, (u[1] + v[1]) % 2)

    iso_ok = True
    for x in ALL_CLASSES:
        for y in ALL_CLASSES:
            xy = phi(x, y)
            if z2_add(z2_map[x], z2_map[y]) != z2_map[xy]:
                iso_ok = False

    assoc = check_associativity(ALL_CLASSES)
    return {
        "group": "Klein four-group V4 ≅ (Z/2)²",
        "closure_mod12": closure_ok,
        "identity_E": identity_ok,
        "self_inverse_nonE": inv_ok,
        "isomorphism_to_Z2_squared": iso_ok,
        "associative": assoc["associative"],
        "triples_tested": assoc["triples_tested"],
        "counterexample_count": assoc["counterexample_count"],
        "naive_associator_all_zero": assoc["associative"],
        "verdict": (
            "V₄ ist Klein-Gruppe; Φ assoziativ; naive 𝔞(X,Y,Z)≡0 für alle Tripel. "
            "Echter EABC-Defekt = projektionsbasierte Holonomie (collatz_eabc_holonomie.md)."
        ),
    }


def check_associativity(classes: tuple[EClass, ...] = ALL_CLASSES) -> dict[str, Any]:
    counterexamples: list[dict[str, str]] = []
    for x in classes:
        for y in classes:
            for z in classes:
                if phi_left(x, y, z) is not phi_right(x, y, z):
                    counterexamples.append(
                        {
                            "x": x.value,
                            "y": y.value,
                            "z": z.value,
                            "left": phi_left(x, y, z).value,
                            "right": phi_right(x, y, z).value,
                        }
                    )
    return {
        "domain": [c.value for c in classes],
        "triples_tested": len(classes) ** 3,
        "associative": len(counterexamples) == 0,
        "counterexample_count": len(counterexamples),
        "counterexamples": counterexamples[:8],
    }


def bracket_preference(chirality: Chirality) -> BracketPref:
    """Heuristik: ABCE ↔ linke Klammer ((·)·), CEAB ↔ rechte Klammer (·(··))."""
    return "left" if chirality is Chirality.ABCE else "right"


def quadruplet_chirality(p: int) -> Chirality:
    mod = p % 12
    if mod == 5:
        return Chirality.ABCE
    if mod == 11:
        return Chirality.CEAB
    raise ValueError(f"Kein Vierlingsstart: p={p} (mod 12 = {mod})")


def quadruplet_word(p: int) -> str:
    return "".join(class_of(n).value for n in q(p))


def chirality_quadruplet_report(limit: int) -> dict[str, Any]:
    """Lean-Test-Vierlinge: ABCE/CEAB und Klammerpräferenz (Heuristik)."""
    rows: list[dict[str, Any]] = []
    abce_count = 0
    ceab_count = 0
    for p in range(5, limit - 7, 2):
        if p % 12 not in (5, 11):
            continue
        if not is_prime_quadruplet(p):
            continue
        chi = quadruplet_chirality(p)
        pref = bracket_preference(chi)
        order = chirality_order(chi)
        a, b, c = order[0], order[1], order[2]
        left_val = phi_left(a, b, c).value
        right_val = phi_right(a, b, c).value
        if chi is Chirality.ABCE:
            abce_count += 1
        else:
            ceab_count += 1
        rows.append(
            {
                "p": p,
                "word": quadruplet_word(p),
                "chirality": chi.value,
                "bracket_preference": pref,
                "triple_abc": [a.value, b.value, c.value],
                "phi_left_abc": left_val,
                "phi_right_abc": right_val,
                "associator_sign": associator_sign(a, b, c),
                "preferred_matches_result": left_val == right_val,
            }
        )
    return {
        "limit": limit,
        "quadruplet_count": len(rows),
        "abce_count": abce_count,
        "ceab_count": ceab_count,
        "heuristic": (
            "ABCE → linke Klammer Φ(Φ(A,B),C); CEAB → rechte Klammer Φ(A,Φ(B,C)); "
            "auf V₄ liefern beide dasselbe Ergebnis (Assoziativität)."
        ),
        "samples": rows[:12],
        "all_associator_zero": all(r["associator_sign"] == 0 for r in rows),
    }


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def first_n_primes(n: int) -> list[int]:
    primes: list[int] = []
    candidate = 2
    while len(primes) < n:
        if _is_prime(candidate):
            primes.append(candidate)
        candidate += 1
    return primes


def prime_eabc_classes(n_primes: int) -> list[dict[str, Any]]:
    """Erste n_primes mit gültiger EABC-Klasse (p mod 12 ∈ {1,5,7,11})."""
    rows: list[dict[str, Any]] = []
    candidate = 2
    while len(rows) < n_primes:
        if _is_prime(candidate):
            cls = class_of(candidate)
            if cls is not None:
                rows.append({"index": len(rows) + 1, "p": candidate, "class": cls.value})
        candidate += 1
    return rows


def prime_associator_mean(n_primes: int) -> dict[str, Any]:
    """
    𝒜(N) = (1/M) Σ 𝔞_E(k) über aufeinanderfolgende Prim-EABC-Tripel.

    k indexiert das mittlere Prim p_k; Tripel (class(p_{k-1}), class(p_k), class(p_{k+1})).
    Nur Primzahlen in den vier EABC-Restklassen.
    """
    rows = prime_eabc_classes(n_primes)
    classes = [EClass(r["class"]) for r in rows]
    primes = [r["p"] for r in rows]
    signs: list[int] = []
    triples: list[dict[str, Any]] = []
    for k in range(1, len(primes) - 1):
        x, y, z = classes[k - 1], classes[k], classes[k + 1]
        s = associator_sign(x, y, z)
        signs.append(s)
        if len(triples) < 10:
            triples.append(
                {
                    "k": k + 1,
                    "primes": [primes[k - 1], primes[k], primes[k + 1]],
                    "classes": [x.value, y.value, z.value],
                    "associator_sign": s,
                }
            )
    mean = sum(signs) / len(signs) if signs else 0.0
    nonzero = sum(1 for s in signs if s != 0)
    return {
        "n_primes": n_primes,
        "triples_count": len(signs),
        "A_N": mean,
        "nonzero_count": nonzero,
        "stable_nonzero": nonzero > 0 and abs(mean) > 1e-9,
        "verdict": (
            "𝒜(N) ≡ 0 — Φ assoziativ auf V₄; kein diskreter Klammerdefekt in der Primfolge."
            if nonzero == 0
            else f"𝒜(N) = {mean:.6f} mit {nonzero} nichttrivialen Tripeln"
        ),
        "sample_triples": triples,
    }


def run(prime_limit: int = 10_000, n_primes: int = 200, output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    v4_proof = prove_v4_klein_associativity()
    assoc_full = check_associativity(ALL_CLASSES)
    assoc_abc = check_associativity(ABC_CLASSES)
    report: dict[str, Any] = {
        "meta": {
            "module": "collatz_eabc_discrete_associator.py",
            "theory": "collatz_eabc_discrete_associator.md",
            "holonomy_correction": "collatz_eabc_holonomie.md",
            "phi_definition": (
                "Φ(X,Y) = classOf(residue(X)·residue(Y) mod 12); "
                "V₄ ≅ (ℤ/12ℤ)×{1,5,7,11} unter Multiplikation; E neutral."
            ),
            "lean_sources": ["EABC.lean (residue, classOf)", "eabc_from_lean.py"],
        },
        "v4_klein_proof": v4_proof,
        "multiplication_table_V4": multiplication_table(),
        "associativity": {
            "full_V4": assoc_full,
            "ABC_subtriple": assoc_abc,
            "honest_verdict": (
                "Φ ist auf ganz V₄ assoziativ (Klein-Vierergruppe). "
                "𝔞 ≡ 0 für alle 4³ Tripel; auch auf {A,B,C}³."
                if assoc_full["associative"]
                else "Φ ist NICHT assoziativ auf V₄."
            ),
        },
        "associator_tables": {
            "full_V4_sign": associator_table(ALL_CLASSES),
            "ABC_subtriple_sign": associator_table(ABC_CLASSES),
        },
        "chirality_quadruplets": chirality_quadruplet_report(prime_limit),
        "prime_associator": prime_associator_mean(n_primes),
        "epistemic_labels": {
            "phi_V4": "Definition",
            "associativity_V4": "Theorem",
            "ABCE_CEAB_bracket_link": "Heuristik",
            "A_N_prime_sequence": "Experiment",
        },
    }
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Diskreter EABC-Assoziator auf V₄")
    parser.add_argument("--prime-limit", type=int, default=10_000, help="Obergrenze Vierlingssuche")
    parser.add_argument("--n-primes", type=int, default=200, help="Anzahl Primzahlen für 𝒜(N)")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(prime_limit=args.prime_limit, n_primes=args.n_primes, output=args.output)
    assoc = report["associativity"]["full_V4"]
    print("=== Diskreter EABC-Assoziator (V₄) ===")
    print(report["associativity"]["honest_verdict"])
    print(f"Assoziativ auf V₄: {assoc['associative']}  ({assoc['triples_tested']} Tripel)")
    quad = report["chirality_quadruplets"]
    print(f"Vierlinge bis {quad['limit']}: {quad['quadruplet_count']}  "
          f"(ABCE {quad['abce_count']}, CEAB {quad['ceab_count']})")
    pa = report["prime_associator"]
    print(f"𝒜({pa['n_primes']}) = {pa['A_N']:.6f}  —  {pa['verdict']}")
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
