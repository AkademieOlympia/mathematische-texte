"""
CollatzResonanceSpectrum — Resonanzatlas im Operatorraum der Paritätswörter.

Statt Zahlen n≤N: alle (bzw. kritische) Wörter ω∈{0,1}^m.

Syracuse-Konvention (ein /2 pro Buchstabe):
  0: x ↦ x/2
  1: x ↦ (3x+1)/2

  T_ω(n) = (3^a n + b_ω) / 2^m
  A_ω = 3^a / 2^m
  Λ(ω) = log2(A) = a log2(3) − m

Typen:
  I   dissipativ  A<1  (Λ<0)
  II  kritisch    A=1  (Λ=0)  — für a,m∈ℕ unmöglich (3^a≠2^m)
  III expansiv    A>1  (Λ>0)

ρ_ω = |b| / |2^m − 3^a|

Claim-boundary: Spektrum der Affinklassen; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path


LOG2_3 = math.log2(3)
P_CRIT = math.log(2) / math.log(3)  # log2/log3 ≈ 0.630929


@dataclass(frozen=True)
class WordResonance:
    m: int
    a: int
    b: int
    omega: str
    A: Fraction
    Lambda: float
    rho: Fraction | None
    typ: str
    realizable: bool
    residue: int | None  # n ≡ residue (mod 2^m) falls realisierbar


def affine_offset_syracuse(omega: tuple[int, ...] | str) -> tuple[int, int, Fraction, Fraction]:
    """Liefert (a, b, A, B) mit T=(3^a n+b)/2^m, B=b/2^m."""
    if isinstance(omega, str):
        bits = tuple(int(c) for c in omega)
    else:
        bits = omega
    m = len(bits)
    A = Fraction(1)
    B = Fraction(0)
    a = 0
    for e in bits:
        if e == 0:
            A /= 2
            B /= 2
        else:
            A *= Fraction(3, 2)
            B = (3 * B + 1) / 2
            a += 1
    # B = b / 2^m
    b = int(B * (1 << m))
    assert A == Fraction(3**a, 1 << m)
    return a, b, A, B


def Lambda_of(a: int, m: int) -> float:
    return a * LOG2_3 - m


def type_of(A: Fraction) -> str:
    if A < 1:
        return "I"
    if A == 1:
        return "II"
    return "III"


def rho_of(b: int, a: int, m: int) -> Fraction | None:
    den = (1 << m) - 3**a
    if den == 0:
        return None
    return Fraction(abs(b), abs(den))


def _matches_word(n: int, omega: tuple[int, ...]) -> bool:
    x = n
    for e in omega:
        bit = 0 if x % 2 == 0 else 1
        if bit != e:
            return False
        if e == 0:
            x //= 2
        else:
            x = (3 * x + 1) // 2
    return True


def realizable_residue(omega: tuple[int, ...]) -> int | None:
    """
    Restklasse mod 2^m mit Syracuse-Wort ω, oder None.
    Inkrementeller Lift — O(#Kandidaten · m²), typisch eindeutig.
    """
    return _lift_residue(omega)


def _lift_residue(omega: tuple[int, ...]) -> int | None:
    """Lift der Restklasse Schritt für Schritt mod 2,4,...,2^m."""
    if not omega:
        return 0
    candidates = {0}  # mod 1
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
                    b = 0 if x % 2 == 0 else 1
                    if b != omega[t]:
                        ok = False
                        break
                    x = x // 2 if b == 0 else (3 * x + 1) // 2
                if ok:
                    new_cand.add(r % new_mod)
        if not new_cand:
            return None
        candidates = new_cand
        cur_mod = new_mod
    return min(candidates)


def analyze_word(omega: tuple[int, ...], check_real: bool = True) -> WordResonance:
    m = len(omega)
    a, b, A, _B = affine_offset_syracuse(omega)
    Lam = Lambda_of(a, m)
    typ = type_of(A)
    rho = rho_of(b, a, m)
    if check_real:
        res = realizable_residue(omega)
        real = res is not None
    else:
        res, real = None, True
    return WordResonance(
        m=m,
        a=a,
        b=b,
        omega="".join(str(e) for e in omega),
        A=A,
        Lambda=Lam,
        rho=rho,
        typ=typ,
        realizable=real,
        residue=res,
    )


def iter_all_words(m: int):
    for mask in range(1 << m):
        yield tuple((mask >> i) & 1 for i in range(m))


def iter_words_with_a(m: int, a: int):
    """Alle Wörter Länge m mit genau a Einsen."""
    for idxs in combinations(range(m), a):
        w = [0] * m
        for i in idxs:
            w[i] = 1
        yield tuple(w)


def spectrum_full(m: int, check_real: bool = True) -> dict:
    types = Counter()
    lam_by_a: dict[int, float] = {}
    closest = []  # kleinste |Λ|
    expansive_real = []
    critical_near = []
    realizable = 0
    total = 1 << m

    best_abs = None
    for omega in iter_all_words(m):
        wr = analyze_word(omega, check_real=check_real)
        types[wr.typ] += 1
        lam_by_a[wr.a] = wr.Lambda
        if wr.realizable:
            realizable += 1
        absL = abs(wr.Lambda)
        if best_abs is None or absL < best_abs - 1e-15:
            best_abs = absL
            closest = [wr]
        elif abs(absL - best_abs) < 1e-12:
            if len(closest) < 12:
                closest.append(wr)
        if wr.typ == "III" and wr.realizable and len(expansive_real) < 20:
            expansive_real.append(wr)
        if absL < 0.15 and wr.realizable and len(critical_near) < 30:
            critical_near.append(wr)

    # Histogramm über a (Λ hängt nur von a ab)
    hist = []
    for a in range(m + 1):
        count = math.comb(m, a)
        Lam = Lambda_of(a, m)
        hist.append(
            {
                "a": a,
                "count": count,
                "Lambda": Lam,
                "density": a / m if m else 0.0,
                "typ": type_of(Fraction(3**a, 1 << m)),
            }
        )

    return {
        "m": m,
        "total": total,
        "realizable": realizable,
        "types": dict(types),
        "histogram_by_a": hist,
        "closest_to_zero": [_wr_pub(w) for w in closest],
        "expansive_realizable_sample": [_wr_pub(w) for w in expansive_real],
        "near_critical_sample": [_wr_pub(w) for w in sorted(critical_near, key=lambda w: abs(w.Lambda))[:15]],
        "min_abs_Lambda": best_abs,
        "a_crit_real": P_CRIT * m,
    }


def spectrum_critical_slice(m: int, max_enum: int = 200_000, sample_cap: int = 4000) -> dict:
    """
    a ∈ Umgebung von p_c·m.
    Volle Enumeration nur wenn C(m,a) ≤ max_enum, sonst Stichprobe.
    """
    a_floor = int(math.floor(P_CRIT * m))
    a_ceil = int(math.ceil(P_CRIT * m))
    as_ = sorted(set([a_floor, a_ceil, max(0, a_floor - 1), min(m, a_ceil + 1)]))

    rows = []
    closest = None
    closest_abs = None

    for a in as_:
        Lam = Lambda_of(a, m)
        n_words = math.comb(m, a)
        words_iter = list(iter_words_with_a(m, a)) if n_words <= max_enum else None
        if words_iter is None:
            # Stichprobe: zufällige Indexmengen
            import random

            rng = random.Random(m * 1009 + a)
            sample = []
            seen = set()
            while len(sample) < min(sample_cap, n_words):
                idxs = tuple(sorted(rng.sample(range(m), a)))
                if idxs in seen:
                    continue
                seen.add(idxs)
                w = [0] * m
                for i in idxs:
                    w[i] = 1
                sample.append(tuple(w))
            words_iter = sample
            mode = "sample"
        else:
            mode = "full"

        real_a = 0
        local_best = None
        b_vals = []
        for omega in words_iter:
            wr = analyze_word(omega, check_real=True)
            if wr.realizable:
                real_a += 1
                b_vals.append(wr.b)
            absL = abs(wr.Lambda)
            if closest_abs is None or absL < closest_abs:
                closest_abs = absL
                closest = wr
            if local_best is None or abs(wr.Lambda) < abs(local_best.Lambda):
                local_best = wr
            elif (
                local_best is not None
                and abs(abs(wr.Lambda) - abs(local_best.Lambda)) < 1e-15
                and wr.realizable
                and not local_best.realizable
            ):
                local_best = wr

        rows.append(
            {
                "a": a,
                "Lambda": Lam,
                "density": a / m,
                "words_binom": n_words,
                "inspected": len(words_iter),
                "mode": mode,
                "realizable_inspected": real_a,
                "typ": type_of(Fraction(3**a, 1 << m)),
                "best": _wr_pub(local_best) if local_best else None,
                "b_mean_abs": (sum(abs(b) for b in b_vals) / len(b_vals)) if b_vals else None,
            }
        )

    return {
        "m": m,
        "a_targets": as_,
        "p_crit": P_CRIT,
        "by_a": rows,
        "closest": _wr_pub(closest) if closest else None,
        "min_abs_Lambda": closest_abs,
    }


def _wr_pub(wr: WordResonance | None) -> dict | None:
    if wr is None:
        return None
    return {
        "omega": wr.omega,
        "m": wr.m,
        "a": wr.a,
        "b": wr.b,
        "A": str(wr.A),
        "Lambda": wr.Lambda,
        "rho": str(wr.rho) if wr.rho is not None else None,
        "typ": wr.typ,
        "realizable": wr.realizable,
        "residue": wr.residue,
        "density": wr.a / wr.m if wr.m else None,
    }


def main() -> int:
    print("=" * 72)
    print("COLLATZ-RESONANZSPEKTRUM  —  Operatorraum Ω_m = {0,1}^m")
    print("=" * 72)
    print(f"p_c = log2/log3 = {P_CRIT:.12f}")
    print("Λ = a log2(3) − m    Typ I/II/III nach A ≶ 1")
    print("Hinweis: Typ II (A=1) ist für a,m>0 unmöglich (3^a ≠ 2^m).")

    full_ms = list(range(1, 17))  # 2^16 = 65536
    full_results = []
    print("\n[1] Vollständiger Atlas m=1..16")
    print(f"  {'m':>3}  {'total':>7}  {'real':>7}  {'I':>7}  {'III':>7}  {'min|Λ|':>10}  {'a*':>4}")
    for m in full_ms:
        # Realisierbarkeit für m≤16 ok; Scan 2^m
        check = m <= 16
        res = spectrum_full(m, check_real=check)
        full_results.append(res)
        a_star = min(res["histogram_by_a"], key=lambda h: abs(h["Lambda"]))["a"]
        print(
            f"  {m:3d}  {res['total']:7d}  {res['realizable']:7d}  "
            f"{res['types'].get('I',0):7d}  {res['types'].get('III',0):7d}  "
            f"{res['min_abs_Lambda']:10.6f}  {a_star:4d}"
        )

    # Detail m=16
    detail = full_results[-1]
    print("\n[2] Detail m=16 — Λ-Histogramm nach a")
    for h in detail["histogram_by_a"]:
        mark = " <-- nearest" if abs(h["Lambda"] - detail["min_abs_Lambda"]) < 1e-12 else ""
        if h["a"] in range(8, 13) or mark:
            print(
                f"  a={h['a']:2d}  dens={h['density']:.4f}  Λ={h['Lambda']:+.6f}  "
                f"typ={h['typ']}  #={h['count']}{mark}"
            )

    print("\n[3] Wörter nächst Λ=0 (m=16)")
    for w in detail["closest_to_zero"][:6]:
        print(
            f"  ω={w['omega']}  a={w['a']}  Λ={w['Lambda']:+.6f}  "
            f"A={w['A']}  ρ={w['rho']}  real={w['realizable']}  r={w['residue']}"
        )

    print("\n[4] Expandierende realisierbare Wörter (Sample m=16)")
    if detail["expansive_realizable_sample"]:
        for w in detail["expansive_realizable_sample"][:8]:
            print(
                f"  ω={w['omega']}  a={w['a']}  Λ={w['Lambda']:+.6f}  "
                f"A={w['A']}  real={w['realizable']}"
            )
    else:
        print("  keine (bei m=16: expansiv ⇒ a≥11, Dichte≥0.6875)")

    # Kritische Schnitte für größere m
    slice_ms = [18, 20, 22, 24, 26, 28, 30]
    slice_results = []
    print("\n[5] Kritischer Schnitt a≈p_c·m")
    print(f"  {'m':>3}  {'a-set':>16}  {'inspected':>10}  {'min|Λ|':>10}  {'typ*':>4}")
    for m in slice_ms:
        res = spectrum_critical_slice(m)
        slice_results.append(res)
        best = res["closest"]
        typ = best["typ"] if best else "?"
        inspected = sum(row["inspected"] for row in res["by_a"])
        print(
            f"  {m:3d}  {str(res['a_targets']):>16}  {inspected:10d}  "
            f"{res['min_abs_Lambda']:10.6f}  {typ:>4}"
        )

    print("\n[6] Nächste kritische Affinklassen (Slices)")
    for res in slice_results:
        w = res["closest"]
        if not w:
            continue
        print(
            f"  m={res['m']:2d}  a={w['a']}  dens={w['density']:.6f}  "
            f"Λ={w['Lambda']:+.8f}  A={w['A']}  ω={w['omega'][:40]}"
            f"{'…' if len(w['omega'])>40 else ''}  real={w['realizable']}"
        )

    # Globale Nähe zu Resonanz über m
    print("\n[7] Beste Annäherung an Λ=0 über m=1..30 (nur a-Optimum)")
    best_global = None
    for m in range(1, 31):
        a_opt = min(range(m + 1), key=lambda a: abs(Lambda_of(a, m)))
        Lam = Lambda_of(a_opt, m)
        rec = {"m": m, "a": a_opt, "Lambda": Lam, "density": a_opt / m, "A": str(Fraction(3**a_opt, 1 << m))}
        if best_global is None or abs(Lam) < abs(best_global["Lambda"]):
            best_global = rec
        if m in (1, 2, 5, 12, 17, 29, 30) or abs(Lam) < 0.02:
            print(
                f"  m={m:2d}  a={a_opt:2d}  dens={a_opt/m:.6f}  "
                f"Λ={Lam:+.8f}  A={rec['A']}"
            )
    print(f"  GLOBAL BEST: m={best_global['m']} a={best_global['a']} "
          f"Λ={best_global['Lambda']:+.8f}  A={best_global['A']}")

    print("\n[8] Fazit")
    print(
        "  Λ hängt nur von (a,m) ab — Wortstruktur steckt in b_ω und ρ.\n"
        "  Typ II (A=1) tritt nie auf. Kritische Wörter haben a/m≈p_c.\n"
        "  Collatz-Vermutung ⇔ kein unendlicher Pfad mit limsup Dichte ≥ p_c.\n"
        "  Spektrum = Klassifikation der Affinoperatoren, nicht Primselektion."
    )

    out = Path(__file__).resolve().parent / "collatz_resonance_spectrum.json"
    payload = {
        "claim": (
            "Collatz affine resonance spectrum over parity words; "
            "Lambda=a*log2(3)-m; types I/II/III; not a proof"
        ),
        "p_crit": P_CRIT,
        "full_atlas_m": full_ms,
        "full": [
            {
                "m": r["m"],
                "total": r["total"],
                "realizable": r["realizable"],
                "types": r["types"],
                "min_abs_Lambda": r["min_abs_Lambda"],
                "closest_to_zero": r["closest_to_zero"][:5],
                "histogram_by_a": r["histogram_by_a"],
            }
            for r in full_results
        ],
        "critical_slices": slice_results,
        "best_global_Lambda": best_global,
        "definition": {
            "A": "3^a / 2^m",
            "Lambda": "a*log2(3) - m",
            "rho": "|b| / |2^m - 3^a|",
            "convention": "syracuse",
        },
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
