"""
CollatzPathResonance — globale Pfadresonanz auf kritischen Affinblöcken.

Lokale Resonanz Λ(ω)=a log2(3)−m ist klassifiziert.
Globale Frage:

  ∃ zulässiger Pfad (ω_i) mit limsup (1/N) Σ Λ(ω_i) ≥ 0 ?

Für kritische Blöcke (a/m ≈ p_c) ist jedes Λ_i ≈ 0 bzw. leicht positiv;
ein Operatorzyklus nur aus solchen Blöcken hätte

  Λ_path = Σ Λ_i = k·Λ_0 ,  A_cycle = A_0^k .

R_C(ω) = (3^a, 2^m, B_ω)  — rationale Geometrie jenseits von Λ.

Claim-boundary: Pfad-/Zyklusdiagnostik; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
import random
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path


LOG2_3 = math.log2(3)
P_CRIT = math.log(2) / math.log(3)


@dataclass(frozen=True)
class AffinBlock:
    """Ein Paritätsblock ω mit exakter Affinform."""

    omega: str
    m: int
    a: int
    b: int  # T = (3^a n + b) / 2^m
    residue: int  # n ≡ residue (mod 2^m)

    @property
    def A(self) -> Fraction:
        return Fraction(3**self.a, 1 << self.m)

    @property
    def B(self) -> Fraction:
        return Fraction(self.b, 1 << self.m)

    @property
    def Lambda(self) -> float:
        return self.a * LOG2_3 - self.m

    @property
    def R_C(self) -> tuple[int, int, Fraction]:
        """Rationale Klasse (3^a, 2^m, B)."""
        return 3**self.a, 1 << self.m, self.B


def affine_from_omega(omega: tuple[int, ...]) -> tuple[int, int, Fraction, Fraction]:
    A = Fraction(1)
    B = Fraction(0)
    a = 0
    for e in omega:
        if e == 0:
            A /= 2
            B /= 2
        else:
            A *= Fraction(3, 2)
            B = (3 * B + 1) / 2
            a += 1
    m = len(omega)
    b = int(B * (1 << m))
    return a, b, A, B


def lift_residue(omega: tuple[int, ...]) -> int | None:
    candidates = {0}
    cur_mod = 1
    for j in range(len(omega)):
        new_mod = cur_mod * 2
        new_cand = set()
        for r0 in candidates:
            for half in (0, 1):
                r = r0 + half * cur_mod
                x = r
                ok = True
                for t in range(j + 1):
                    bit = 0 if x % 2 == 0 else 1
                    if bit != omega[t]:
                        ok = False
                        break
                    x = x // 2 if bit == 0 else (3 * x + 1) // 2
                if ok:
                    new_cand.add(r % new_mod)
        if not new_cand:
            return None
        candidates = new_cand
        cur_mod = new_mod
    return min(candidates)


def matches_word(n: int, omega: tuple[int, ...]) -> bool:
    x = n
    for e in omega:
        bit = 0 if x % 2 == 0 else 1
        if bit != e:
            return False
        x = x // 2 if e == 0 else (3 * x + 1) // 2
    return True


def make_block(omega: tuple[int, ...]) -> AffinBlock | None:
    a, b, _A, _B = affine_from_omega(omega)
    r = lift_residue(omega)
    if r is None:
        return None
    return AffinBlock(
        omega="".join(str(e) for e in omega),
        m=len(omega),
        a=a,
        b=b,
        residue=r,
    )


def compose_AB(blocks: list[AffinBlock]) -> tuple[Fraction, Fraction, float]:
    """(A,B,Λ_path) der Komposition blocks[0] dann blocks[1] … (links zuerst)."""
    A = Fraction(1)
    B = Fraction(0)
    Lam = 0.0
    for bl in blocks:
        # T_new = T_bl ∘ T_old  ⇒  A' = A_bl * A,  B' = A_bl * B + B_bl
        A = bl.A * A
        B = bl.A * B + bl.B
        Lam += bl.Lambda
    return A, B, Lam


def cycle_candidate(A: Fraction, B: Fraction) -> int | None:
    if A == 1:
        return None
    n = -B / (A - 1)
    if n.denominator == 1 and n > 0:
        return int(n)
    return None


def concat_bits(blocks: list[AffinBlock]) -> tuple[int, ...]:
    bits: list[int] = []
    for bl in blocks:
        bits.extend(int(c) for c in bl.omega)
    return tuple(bits)


def verify_cycle(n: int, blocks: list[AffinBlock]) -> bool:
    return matches_word(n, concat_bits(blocks)) and _apply_blocks(n, blocks) == n


def _apply_blocks(n: int, blocks: list[AffinBlock]) -> int:
    x = n
    for bl in blocks:
        for e in bl.omega:
            x = x // 2 if e == 0 else (3 * x + 1) // 2
    return x


def iter_words_with_a(m: int, a: int):
    for idxs in combinations(range(m), a):
        w = [0] * m
        for i in idxs:
            w[i] = 1
        yield tuple(w)


def build_critical_blocks(m: int, a: int) -> list[AffinBlock]:
    blocks = []
    for omega in iter_words_with_a(m, a):
        bl = make_block(omega)
        if bl is not None:
            blocks.append(bl)
    return blocks


def scan_k1_cycles(blocks: list[AffinBlock]) -> list[dict]:
    found = []
    A0 = blocks[0].A if blocks else Fraction(1)
    for bl in blocks:
        A, B, Lam = compose_AB([bl])
        n = cycle_candidate(A, B)
        if n is None:
            continue
        ok = verify_cycle(n, [bl])
        rec = {
            "k": 1,
            "n": n,
            "verified": ok,
            "A": str(A),
            "B": str(B),
            "Lambda_path": Lam,
            "mean_Lambda": Lam,
            "omegas": [bl.omega],
            "R_C": [list(bl.R_C[:2]) + [str(bl.R_C[2])]],
        }
        if ok or len(found) < 5:
            found.append(rec)
        if ok:
            # behalte alle verifizierten
            pass
    verified = [r for r in found if r["verified"]]
    return verified if verified else []


def scan_path_cycles(
    blocks: list[AffinBlock],
    k: int,
    samples: int,
    rng: random.Random,
) -> tuple[list[dict], dict]:
    """Zufällige k-Block-Pfade; Zykluskandidaten + Stats zu |A-1|, Λ_path."""
    verified = []
    candidates = []
    min_abs_A1 = None
    min_abs_rec = None
    max_mean_Lam = None
    max_mean_rec = None
    A_ge1 = 0
    checked = 0

    if not blocks:
        return [], {}

    for _ in range(samples):
        path = [rng.choice(blocks) for _ in range(k)]
        A, B, Lam = compose_AB(path)
        checked += 1
        mean = Lam / k
        absA1 = abs(float(A - 1))
        if A >= 1:
            A_ge1 += 1
        if min_abs_A1 is None or absA1 < min_abs_A1:
            min_abs_A1 = absA1
            min_abs_rec = {
                "k": k,
                "A": str(A),
                "abs_A_minus_1": absA1,
                "Lambda_path": Lam,
                "mean_Lambda": mean,
                "omegas": [b.omega for b in path],
            }
        if max_mean_Lam is None or mean > max_mean_Lam:
            max_mean_Lam = mean
            max_mean_rec = {
                "k": k,
                "A": str(A),
                "Lambda_path": Lam,
                "mean_Lambda": mean,
                "omegas": [b.omega for b in path],
            }
        n = cycle_candidate(A, B)
        if n is None:
            continue
        ok = verify_cycle(n, path)
        rec = {
            "k": k,
            "n": n,
            "verified": ok,
            "A": str(A),
            "B": str(B),
            "Lambda_path": Lam,
            "mean_Lambda": mean,
            "omegas": [b.omega for b in path],
        }
        candidates.append(rec)
        if ok:
            verified.append(rec)

    stats = {
        "k": k,
        "samples": checked,
        "A_ge1_count": A_ge1,
        "A_ge1_rate": A_ge1 / checked if checked else None,
        "rational_cycle_candidates": len(candidates),
        "verified_cycles": len(verified),
        "min_abs_A_minus_1": min_abs_rec,
        "max_mean_Lambda": max_mean_rec,
    }
    return verified, stats


def scan_repeated_block_cycles(blocks: list[AffinBlock], k_max: int = 6) -> list[dict]:
    """ω^k — derselbe Block k-mal (starke Resonanzkandidaten)."""
    found = []
    for bl in blocks:
        for k in range(1, k_max + 1):
            path = [bl] * k
            A, B, Lam = compose_AB(path)
            n = cycle_candidate(A, B)
            if n is None:
                continue
            if verify_cycle(n, path):
                found.append(
                    {
                        "k": k,
                        "n": n,
                        "verified": True,
                        "A": str(A),
                        "B": str(B),
                        "Lambda_path": Lam,
                        "mean_Lambda": Lam / k,
                        "omegas": [bl.omega],
                        "mode": "repeat",
                    }
                )
    return found


def b_diversity(blocks: list[AffinBlock]) -> dict:
    """Λ gleich, B verschieden — analog M_B vs spec(PMP)."""
    Bs = [bl.B for bl in blocks]
    uniq = len(set(Bs))
    return {
        "blocks": len(blocks),
        "unique_B": uniq,
        "unique_b_int": len({bl.b for bl in blocks}),
        "B_min": str(min(Bs)) if Bs else None,
        "B_max": str(max(Bs)) if Bs else None,
        "sample_R_C": [
            {"omega": bl.omega, "R_C": [bl.R_C[0], bl.R_C[1], str(bl.R_C[2])]}
            for bl in blocks[:5]
        ],
    }


def run_layer(m: int, a: int, rng: random.Random, samples_k2: int, samples_k3: int) -> dict:
    print(f"\n--- Layer m={m}, a={a}, dens={a/m:.6f}, Λ={a*LOG2_3-m:+.8f} ---")
    blocks = build_critical_blocks(m, a)
    print(f"  Blöcke realisierbar: {len(blocks)} / binom={math.comb(m,a)}")
    if not blocks:
        return {"m": m, "a": a, "blocks": 0}

    A0 = blocks[0].A
    Lam0 = blocks[0].Lambda
    print(f"  A_0={A0}  Λ_0={Lam0:+.8f}  typ={'III' if A0>1 else 'I' if A0<1 else 'II'}")

    div = b_diversity(blocks)
    print(f"  B-Diversität: unique_B={div['unique_B']}  (Λ fix, Geometrie in B)")

    # k=1 voll
    v1 = []
    cand1 = 0
    min_abs = None
    min_abs_bl = None
    for bl in blocks:
        A, B, Lam = compose_AB([bl])
        absA1 = abs(float(A - 1))
        if min_abs is None or absA1 < min_abs:
            min_abs = absA1
            min_abs_bl = bl
        n = cycle_candidate(A, B)
        if n is not None:
            cand1 += 1
            if verify_cycle(n, [bl]):
                v1.append(
                    {
                        "k": 1,
                        "n": n,
                        "A": str(A),
                        "B": str(B),
                        "Lambda_path": Lam,
                        "omega": bl.omega,
                    }
                )
    print(f"  k=1: rationale Kandidaten={cand1}  verifiziert={len(v1)}")
    if min_abs_bl is not None:
        print(f"  min|A-1|={min_abs:.6e}  (A konstant = {A0})")

    # Wiederholungen ω^k
    # zu teuer über alle blocks × k — sample + alle falls klein
    rep_found = []
    sample_rep = blocks if len(blocks) <= 5000 else rng.sample(blocks, 5000)
    for bl in sample_rep:
        for k in range(2, 5):
            path = [bl] * k
            A, B, Lam = compose_AB(path)
            n = cycle_candidate(A, B)
            if n is not None and verify_cycle(n, path):
                rep_found.append(
                    {
                        "k": k,
                        "n": n,
                        "A": str(A),
                        "Lambda_path": Lam,
                        "omega": bl.omega,
                        "mode": "repeat",
                    }
                )
    print(f"  ω^k (k=2..4, sample={len(sample_rep)}): verifiziert={len(rep_found)}")

    v2, st2 = scan_path_cycles(blocks, 2, samples_k2, rng)
    print(
        f"  k=2 sample={st2['samples']}: A≥1 rate={st2['A_ge1_rate']:.4f}  "
        f"rat.cand={st2['rational_cycle_candidates']}  verified={st2['verified_cycles']}"
    )
    v3, st3 = scan_path_cycles(blocks, 3, samples_k3, rng)
    print(
        f"  k=3 sample={st3['samples']}: A≥1 rate={st3['A_ge1_rate']:.4f}  "
        f"rat.cand={st3['rational_cycle_candidates']}  verified={st3['verified_cycles']}"
    )

    # gemischte Dichten: später in main

    return {
        "m": m,
        "a": a,
        "density": a / m,
        "Lambda_0": Lam0,
        "A_0": str(A0),
        "n_blocks": len(blocks),
        "B_diversity": div,
        "k1_candidates": cand1,
        "k1_verified": v1,
        "repeat_verified": rep_found,
        "k2": st2,
        "k2_verified": v2,
        "k3": st3,
        "k3_verified": v3,
        "min_abs_A_minus_1": min_abs,
        "_blocks": blocks,  # intern für Cache
    }


def mixed_sign_paths(
    blocks_pos: list[AffinBlock],
    blocks_neg: list[AffinBlock],
    samples: int,
    rng: random.Random,
    k: int = 4,
) -> dict:
    """
    Pfade mit gemischten Λ-Vorzeichen — kann mittlere Verstärkung ~0 halten,
    ohne dass jeder Block expansiv ist.
    """
    if not blocks_pos or not blocks_neg:
        return {"samples": 0}
    verified = []
    mean_pos = 0
    mean_ge0 = 0
    best_mean = None
    best = None
    for _ in range(samples):
        path = []
        for i in range(k):
            path.append(rng.choice(blocks_pos if rng.random() < 0.5 else blocks_neg))
        A, B, Lam = compose_AB(path)
        mean = Lam / k
        if mean > 0:
            mean_pos += 1
        if mean >= 0:
            mean_ge0 += 1
        if best_mean is None or abs(mean) < abs(best_mean):
            best_mean = mean
            best = {
                "mean_Lambda": mean,
                "Lambda_path": Lam,
                "A": str(A),
                "omegas": [b.omega for b in path],
                "Lambdas": [b.Lambda for b in path],
            }
        n = cycle_candidate(A, B)
        if n is not None and verify_cycle(n, path):
            verified.append(
                {
                    "n": n,
                    "A": str(A),
                    "mean_Lambda": mean,
                    "omegas": [b.omega for b in path],
                }
            )
    return {
        "k": k,
        "samples": samples,
        "mean_Lambda_pos_rate": mean_pos / samples,
        "mean_Lambda_ge0_rate": mean_ge0 / samples,
        "closest_mean_to_zero": best,
        "verified_cycles": verified,
    }


def main() -> int:
    rng = random.Random(19)
    print("=" * 72)
    print("COLLATZ-PFADRESONANZ  —  kritische Blöcke → Operatorzyklen")
    print("=" * 72)
    print(f"p_c={P_CRIT:.12f}")
    print("Frage: ∃ Pfad mit limsup (1/N) Σ Λ(ω_i) ≥ 0 als Operatorzyklus?")

    layers = [
        (19, 12),  # Λ≈+0.01955  gefährlichster Einzelblock
        (19, 11),  # Λ≈−1.565? wait 11*log2(3)-19
        (20, 13),
        (20, 12),
        (21, 13),
    ]
    # correct a=11 for m=19
    results = []
    block_cache: dict[tuple[int, int], list[AffinBlock]] = {}

    for m, a in layers:
        # sample sizes scale with cost
        sk2 = 80_000 if math.comb(m, a) > 100_000 else 150_000
        sk3 = 40_000 if math.comb(m, a) > 100_000 else 80_000
        # building all blocks for large binom is heavy — cap
        nbin = math.comb(m, a)
        if nbin > 120_000:
            print(f"\n--- Layer m={m}, a={a} (Stichprobe, binom={nbin}) ---")
            # sample blocks
            idxs_list = []
            seen = set()
            while len(idxs_list) < 40_000:
                idxs = tuple(sorted(rng.sample(range(m), a)))
                if idxs in seen:
                    continue
                seen.add(idxs)
                idxs_list.append(idxs)
            blocks = []
            for idxs in idxs_list:
                w = [0] * m
                for i in idxs:
                    w[i] = 1
                bl = make_block(tuple(w))
                if bl:
                    blocks.append(bl)
            block_cache[(m, a)] = blocks
            print(f"  sample blocks: {len(blocks)}")
            A0 = blocks[0].A
            Lam0 = blocks[0].Lambda
            print(f"  A_0={A0}  Λ_0={Lam0:+.8f}")
            v2, st2 = scan_path_cycles(blocks, 2, sk2, rng)
            v3, st3 = scan_path_cycles(blocks, 3, sk3, rng)
            # k1 on sample
            v1 = []
            cand1 = 0
            for bl in blocks:
                A, B, Lam = compose_AB([bl])
                n = cycle_candidate(A, B)
                if n is not None:
                    cand1 += 1
                    if verify_cycle(n, [bl]):
                        v1.append({"n": n, "A": str(A), "omega": bl.omega})
            print(f"  k=1: cand={cand1} verified={len(v1)}")
            print(
                f"  k=2: A≥1 rate={st2['A_ge1_rate']:.4f} verified={st2['verified_cycles']}"
            )
            print(
                f"  k=3: A≥1 rate={st3['A_ge1_rate']:.4f} verified={st3['verified_cycles']}"
            )
            results.append(
                {
                    "m": m,
                    "a": a,
                    "mode": "sample",
                    "n_blocks": len(blocks),
                    "Lambda_0": Lam0,
                    "A_0": str(A0),
                    "k1_verified": v1,
                    "k2": st2,
                    "k2_verified": v2,
                    "k3": st3,
                    "k3_verified": v3,
                    "B_diversity": b_diversity(blocks),
                }
            )
        else:
            layer = run_layer(m, a, rng, sk2, sk3)
            block_cache[(m, a)] = layer.pop("_blocks")
            results.append(layer)

    # Gemischte Pfade m=19: a=12 (pos) + a=11 (neg)
    print("\n--- Gemischte Pfade m=19 (a=12 ∪ a=11) ---")
    pos = block_cache.get((19, 12))
    neg = block_cache.get((19, 11))
    if pos is None:
        pos = build_critical_blocks(19, 12)
        block_cache[(19, 12)] = pos
    if neg is None:
        neg = build_critical_blocks(19, 11)
        block_cache[(19, 11)] = neg
    print(f"  |pos|={len(pos)} Λ≈{pos[0].Lambda:+.5f}  |neg|={len(neg)} Λ≈{neg[0].Lambda:+.5f}")
    mix = mixed_sign_paths(pos, neg, samples=100_000, rng=rng, k=4)
    print(
        f"  k=4 mix: meanΛ>0 rate={mix['mean_Lambda_pos_rate']:.4f}  "
        f"meanΛ≥0 rate={mix['mean_Lambda_ge0_rate']:.4f}  "
        f"verified_cycles={len(mix['verified_cycles'])}"
    )
    if mix.get("closest_mean_to_zero"):
        c = mix["closest_mean_to_zero"]
        print(
            f"  closest meanΛ→0: {c['mean_Lambda']:+.6f}  A={c['A']}  "
            f"Λs={[round(x,4) for x in c['Lambdas']]}"
        )

    # Fall-Diagnose
    any_verified = any(
        (r.get("k1_verified") or r.get("k2_verified") or r.get("k3_verified") or r.get("repeat_verified"))
        for r in results
    ) or bool(mix.get("verified_cycles"))

    print("\n[Struktureller Satz — Syracuse]")
    print(
        "  Für jedes Syracuse-Wort gilt B_ω ≥ 0 (Induktion), B=0 ⇒ ω=0^m ⇒ A<1.\n"
        "  Komposition: A_path = Π A_i > 0, B_path ≥ 0.\n"
        "  Ferner Λ_path = log2(A_path), also meanΛ ≥ 0 ⇔ A_path ≥ 1.\n"
        "  Fixpunkt n* = −B/(A−1): wenn A>1 und B≥0 ⇒ n* ≤ 0.\n"
        "  A=1 unmöglich (3^a ≠ 2^m).\n"
        "  ⇒ Es gibt keinen positiven ganzzahligen Operatorzyklus mit A_cycle ≥ 1."
    )
    print("\n[Fall-Diagnose]")
    print(f"  Empirie verified_cycle={any_verified} (erwartet False)")
    print(
        "  Fall A ist damit für geschlossene Operatorzyklen bewiesen, nicht nur gescannt.\n"
        "  Offen bleibt die Nicht-Zyklus-Frage: unendliche zulässige Pfade mit\n"
        "  limsup meanΛ ≥ 0 (ohne Rückkehr) — das ist die eigentliche Collatz-Form."
    )
    print(
        "  Gemischte endliche Pfade können meanΛ≥0 haben; sie liefern aber keinen\n"
        "  positiven Fixpunkt solange A_path≥1 und B_path≥0."
    )

    out = Path(__file__).resolve().parent / "collatz_path_resonance.json"
    payload = {
        "claim": (
            "Path resonance on critical Collatz affine blocks; "
            "search operator cycles with A>=1; not a proof"
        ),
        "p_crit": P_CRIT,
        "question": "exists admissible path with limsup mean Lambda >= 0 as operator cycle?",
        "layers": results,
        "mixed_m19": mix,
        "any_verified_cycle": any_verified,
        "theorem_syracuse": {
            "statement": (
                "No positive integer operator cycle with A_cycle>=1: "
                "B>=0, A>1 => n*=-B/(A-1)<=0; A=1 impossible"
            ),
            "open": "infinite non-cyclic admissible paths with limsup mean Lambda >= 0",
        },
        "interpretation": {
            "case_A": "proved for operator cycles (Syracuse B>=0)",
            "case_B": "ruled out for A>1, B>=0",
            "case_C": "finite paths with A~1 remain as near-resonances",
            "R_C": "(3^a, 2^m, B_omega)",
        },
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
