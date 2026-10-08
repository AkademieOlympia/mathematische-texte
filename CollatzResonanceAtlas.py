"""
CollatzResonanceAtlas — bahnbezogene affine Resonanzklassen.

Entwicklung:
  Parität → Gap → Affine Klasse → Resonanz → Lyapunov-Kandidat

Kanonischer Horizont: erste echte Abstiegzeit
  τ(n) = min{ k > 0 : T^k(n) < n }

Affine Form nach τ Schritten:
  T^τ(n) = (3^a n + c) / 2^τ = A n + B

Bahninvariante:
  C(n) = (a, τ, ρ)  mit  ρ = |c| / (2^τ + 3^a)

Lyapunov-Kandidat:
  V(n) = −log A = τ log 2 − a log 3
  (A < 1 ⇔ a/τ < log2/log3)

Zyklustest: für jedes Paritätswort π der Länge k
  n = −B/(A−1)  (falls A≠1 und n∈ℕ)

Claim-boundary: strukturelle Diagnostik; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path


LOG2 = math.log(2)
LOG3 = math.log(3)
CRIT = LOG2 / LOG3  # ≈ 0.6309


def collatz_step(x: int) -> int:
    """Elementarer Collatz-Schritt (Diagnostik)."""
    return x // 2 if x % 2 == 0 else 3 * x + 1


def syracuse_step(x: int) -> tuple[int, int]:
    """
    Syracuse-Schritt (ein /2 pro Schritt):
      gerade: x ↦ x/2
      ungerade: x ↦ (3x+1)/2
    Liefert (x_neu, bit) mit bit∈{0,1}.
    """
    if x % 2 == 0:
        return x // 2, 0
    y = 3 * x + 1
    assert y % 2 == 0
    return y // 2, 1


def first_descent(n: int, max_steps: int = 100_000) -> tuple[int, int, int, int, list[int]]:
    """
    τ(n) = min{k>0 : T^k(n) < n} in der Syracuse-Map.

    Affine Form:
      T^τ(n) = (3^a n + c) / 2^τ = A n + B
    mit A = 3^a/2^τ, B = c/2^τ.
    """
    if n <= 1:
        return 0, 0, 0, n, []

    a = 0
    tau = 0
    c = 0
    x = n
    bits: list[int] = []
    for _ in range(max_steps):
        if x % 2 == 0:
            # /2: c unverändert relativ zu Nenner·2
            bits.append(0)
            tau += 1
            x //= 2
        else:
            # (3x+1)/2: a+=1, c ← (3c + 2^τ)/2  im Ganzzahlbild
            # vor dem /2: c' = 3c + 2^τ, dann Nenner 2^{τ+1}
            bits.append(1)
            c = 3 * c + (1 << tau)
            a += 1
            tau += 1
            x = (3 * x + 1) // 2
        if x < n:
            num = (3**a) * n + c
            den = 1 << tau
            assert num % den == 0 and num // den == x
            assert tau == len(bits)
            return a, tau, c, x, bits
    raise RuntimeError(f"no descent for n={n} within {max_steps}")


def rho_from_affine(a: int, k: int, c: int) -> Fraction:
    """ρ = |c| / (2^k + 3^a)."""
    return Fraction(abs(c), (1 << k) + 3**a)


def A_from_affine(a: int, k: int) -> Fraction:
    return Fraction(3**a, 1 << k)


def V_from_affine(a: int, k: int) -> float:
    """V = −log A = k log 2 − a log 3."""
    return k * LOG2 - a * LOG3


def C_of(n: int) -> tuple[int, int, Fraction, float, Fraction]:
    """C(n)=(a,τ,ρ) plus V und A."""
    a, tau, c, _x, _bits = first_descent(n)
    rho = rho_from_affine(a, tau, c)
    A = A_from_affine(a, tau)
    V = V_from_affine(a, tau)
    return a, tau, rho, V, A


# ---------------------------------------------------------------------------
# Affine class from parity pattern (Syracuse-style combined bit as in sketch)
# ---------------------------------------------------------------------------

def affine_from_pattern_syracuse(bits: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    """
    Bit=0: /2; Bit=1: (3x+1)/2  (wie User-Skizze).
    Jeder Buchstabe zählt als ein 'k'-Schritt der kombinierten Map.
    """
    A = Fraction(1)
    B = Fraction(0)
    for e in bits:
        if e == 0:
            A /= 2
            B /= 2
        else:
            A *= Fraction(3, 2)
            B = (3 * B + 1) / 2
    return A, B


def affine_from_pattern_elementary(bits: tuple[int, ...]) -> tuple[Fraction, Fraction, int, int, int]:
    """
    Elementare Bits: 0=/2, 1=3n+1 (ohne sofortiges /2).
    Liefert (A,B,a,k,c) mit A=3^a/2^k, B=c/2^k.
    """
    a = 0
    k = 0
    c = 0
    for e in bits:
        if e == 0:
            k += 1
        else:
            c = 3 * c + (1 << k)
            a += 1
    A = Fraction(3**a, 1 << k) if k >= 0 else Fraction(3**a)
    B = Fraction(c, 1 << k) if k > 0 else Fraction(c)
    return A, B, a, k, c


def possible_cycle_from_AB(A: Fraction, B: Fraction) -> int | None:
    if A == 1:
        return None
    n = -B / (A - 1)
    if n.denominator == 1 and n > 0:
        return int(n)
    return None


def verify_cycle(n: int, bits: tuple[int, ...], mode: str) -> bool:
    """Prüfe, ob das Bitmuster von n aus einen Zyklus der Länge len(bits) erzeugt."""
    x = n
    for e in bits:
        if mode == "syracuse":
            if e == 0:
                if x % 2 != 0:
                    return False
                x //= 2
            else:
                if x % 2 == 0:
                    return False
                x = 3 * x + 1
                if x % 2 != 0:
                    return False
                x //= 2
        else:  # elementary
            if e == 0:
                if x % 2 != 0:
                    return False
                x //= 2
            else:
                if x % 2 == 0:
                    return False
                x = 3 * x + 1
    return x == n


def scan_cycles(k_max: int, mode: str = "syracuse") -> dict:
    """Alle Paritätswörter bis Länge k_max; rationale Fixpunkte."""
    by_k: list[dict] = []
    all_found: list[dict] = []
    for k in range(1, k_max + 1):
        found = []
        for bits in product([0, 1], repeat=k):
            if mode == "syracuse":
                A, B = affine_from_pattern_syracuse(bits)
                a = sum(bits)
                kk = k
                # c aus B = c/2^k mit Syracuse-Zählung: jedes Bit ein /2
                c = int(B * (1 << k))
            else:
                A, B, a, kk, c = affine_from_pattern_elementary(bits)
                if kk != k:
                    # elementary: k = Anzahl der Nullen, nicht Wortlänge
                    pass
            n = possible_cycle_from_AB(A, B)
            if n is None:
                continue
            ok = verify_cycle(n, bits, mode)
            rec = {
                "n": n,
                "bits": "".join(str(b) for b in bits),
                "a": a,
                "k": k if mode == "syracuse" else kk,
                "A": str(A),
                "B": str(B),
                "verified": ok,
            }
            found.append(rec)
            if ok:
                all_found.append(rec)
        by_k.append(
            {
                "k": k,
                "candidates": len(found),
                "verified": sum(1 for r in found if r["verified"]),
                "samples": found[:5],
            }
        )
    # bekannter 4-2-1 Zyklus hervorheben
    trivial = [r for r in all_found if r["n"] in (1, 2, 4)]
    nontrivial = [r for r in all_found if r["n"] not in (1, 2, 4)]
    return {
        "mode": mode,
        "k_max": k_max,
        "by_k": by_k,
        "verified_total": len(all_found),
        "trivial_421": trivial[:20],
        "nontrivial_verified": nontrivial[:20],
        "nontrivial_count": len(nontrivial),
    }


# ---------------------------------------------------------------------------
# Atlas über n: Klassen (a,τ,c)
# ---------------------------------------------------------------------------

def build_atlas(N: int, max_steps: int = 100_000) -> dict:
    classes: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    expanding: list[dict] = []  # A ≥ 1
    V_neg: list[dict] = []  # V ≤ 0 ⇔ A ≥ 1
    C_atlas: dict[tuple[int, int, str], list[int]] = defaultdict(list)

    A_ge1 = 0
    descent_fail = 0

    for n in range(2, N):
        try:
            a, tau, c, x, bits = first_descent(n, max_steps=max_steps)
        except RuntimeError:
            descent_fail += 1
            continue
        classes[(a, tau, c)].append(n)
        rho = rho_from_affine(a, tau, c)
        A = A_from_affine(a, tau)
        V = V_from_affine(a, tau)
        C_atlas[(a, tau, str(rho))].append(n)
        if A >= 1:
            A_ge1 += 1
            if len(expanding) < 30:
                expanding.append(
                    {
                        "n": n,
                        "a": a,
                        "tau": tau,
                        "c": c,
                        "A": str(A),
                        "V": V,
                        "rho": str(rho),
                        "T": x,
                        "bits": "".join(str(b) for b in bits),
                    }
                )
            if V <= 0 and len(V_neg) < 30:
                V_neg.append({"n": n, "a": a, "tau": tau, "V": V, "A": str(A)})

    # Stabilität expandierender Klassen: bleibt A≥1 bei Iteration der Abstiegskarte?
    stable_expanding = []
    for rec in expanding:
        n0 = rec["n"]
        chain = [n0]
        x = n0
        unstable = False
        for _ in range(20):
            a, tau, c, x, _bits = first_descent(x)
            chain.append(x)
            if A_from_affine(a, tau) < 1:
                unstable = True
                break
            if x == chain[0]:
                break
        if not unstable and len(chain) > 1:
            stable_expanding.append({"start": n0, "chain": chain[:12], "A0": rec["A"]})

    # Statistik V
    Vs = []
    ratios = []
    for n in range(2, min(N, 5000)):
        a, tau, c, _x, _b = first_descent(n)
        Vs.append(V_from_affine(a, tau))
        ratios.append(a / tau if tau else 0.0)

    # größte Klassen
    largest = sorted(
        (
            {
                "a": a,
                "tau": tau,
                "c": c,
                "A": str(A_from_affine(a, tau)),
                "rho": str(rho_from_affine(a, tau, c)),
                "V": V_from_affine(a, tau),
                "size": len(ns),
                "sample": ns[:8],
            }
            for (a, tau, c), ns in classes.items()
        ),
        key=lambda d: -d["size"],
    )

    # C-Kollisionen vs (a,τ,c)-Klassen
    C_collisions = sum(1 for v in C_atlas.values() if len(v) > 1)

    return {
        "N": N,
        "num_affine_classes": len(classes),
        "num_C_classes": len(C_atlas),
        "C_collisions": C_collisions,
        "A_ge1_count": A_ge1,
        "descent_fail": descent_fail,
        "expanding_samples": expanding,
        "stable_expanding": stable_expanding[:10],
        "largest_classes": largest[:15],
        "V_stats": {
            "mean": sum(Vs) / len(Vs) if Vs else None,
            "min": min(Vs) if Vs else None,
            "max": max(Vs) if Vs else None,
            "frac_V_pos": sum(1 for v in Vs if v > 0) / len(Vs) if Vs else None,
            "frac_ratio_below_crit": sum(1 for r in ratios if r < CRIT) / len(ratios)
            if ratios
            else None,
        },
        "tau_hist": dict(Counter(tau for (_, tau, _) in classes.keys())),
    }


def lyapunov_probe(N: int = 2000) -> dict:
    """Prüfe V(T_desc(n)) < V(n) entlang der Abstiegskette (nicht elementar)."""
    descents = 0
    strict = 0
    weak = 0
    fails = 0
    examples_fail = []
    for n in range(2, N):
        a0, t0, _c0, x1, _b0 = first_descent(n)
        V0 = V_from_affine(a0, t0)
        if x1 <= 1:
            continue
        a1, t1, _c1, _x2, _b1 = first_descent(x1)
        V1 = V_from_affine(a1, t1)
        descents += 1
        if V1 < V0:
            strict += 1
        elif V1 <= V0:
            weak += 1
        else:
            fails += 1
            if len(examples_fail) < 12:
                examples_fail.append(
                    {"n": n, "x1": x1, "V0": V0, "V1": V1, "a0": a0, "t0": t0, "a1": a1, "t1": t1}
                )
    return {
        "checked": descents,
        "strict_decrease": strict,
        "weak_decrease": weak,
        "increase": fails,
        "fail_rate": fails / descents if descents else None,
        "fail_examples": examples_fail,
    }


def main() -> int:
    N = 10_000
    k_max = 16

    print("=" * 72)
    print("COLLATZ-RESONANZ-ATLAS  —  Bahninvariante C = (a, τ, ρ)")
    print("=" * 72)
    print(f"N={N},  ρ=|c|/(2^τ+3^a),  V=τ log2 − a log3,  crit=log2/log3={CRIT:.6f}")

    atlas = build_atlas(N)
    print("\n[1] Affine Klassen (a,τ,c)")
    print(f"  Klassen (a,τ,c)     = {atlas['num_affine_classes']}")
    print(f"  Klassen C=(a,τ,ρ)   = {atlas['num_C_classes']}")
    print(f"  C-Kollisionen       = {atlas['C_collisions']}")
    print(f"  A≥1 (expandierend)  = {atlas['A_ge1_count']}")
    print(f"  Descent-Fails       = {atlas['descent_fail']}")

    print("\n[2] V-Statistik (n<5000)")
    for k, v in atlas["V_stats"].items():
        if isinstance(v, float):
            print(f"  {k:24s} = {v:.6f}")
        else:
            print(f"  {k:24s} = {v}")

    print("\n[3] Größte Resonanzklassen")
    for d in atlas["largest_classes"][:8]:
        print(
            f"  size={d['size']:4d}  (a,τ,c)=({d['a']},{d['tau']},{d['c']})  "
            f"A={d['A']}  V={d['V']:.4f}  sample={d['sample']}"
        )

    print("\n[4] Expandierende Abstiege A≥1")
    if not atlas["expanding_samples"]:
        print("  keine")
    else:
        for rec in atlas["expanding_samples"][:10]:
            print(
                f"  n={rec['n']:5d}  a/τ={rec['a']}/{rec['tau']}  "
                f"A={rec['A']}  V={rec['V']:.4f}  bits={rec['bits']}"
            )
    print(f"  stabile Expanding-Ketten: {len(atlas['stable_expanding'])}")
    for rec in atlas["stable_expanding"][:5]:
        print(f"    start={rec['start']} chain={rec['chain']}")

    print("\n[5] Lyapunov-Sonde V entlang Abstiegskarte")
    lya = lyapunov_probe(min(N, 3000))
    print(f"  checked={lya['checked']}  strict↓={lya['strict_decrease']}  "
          f"↑={lya['increase']}  fail_rate={lya['fail_rate']}")
    if lya["fail_examples"]:
        print("  Gegenbeispiele V(Tτ(n))≮V(n):")
        for ex in lya["fail_examples"][:6]:
            print(
                f"    n={ex['n']}→{ex['x1']}  V:{ex['V0']:.4f}→{ex['V1']:.4f}  "
                f"(a,τ)=({ex['a0']},{ex['t0']})→({ex['a1']},{ex['t1']})"
            )

    print("\n[6] Zyklussuche (Syracuse-Bits, k≤{})".format(k_max))
    cyc_s = scan_cycles(k_max, mode="syracuse")
    for row in cyc_s["by_k"]:
        if row["k"] <= 12 or row["verified"]:
            print(
                f"  k={row['k']:2d}  candidates={row['candidates']:5d}  "
                f"verified={row['verified']:3d}"
            )
    print(f"  verified total={cyc_s['verified_total']}  "
          f"nontrivial={cyc_s['nontrivial_count']}")
    if cyc_s["nontrivial_verified"]:
        print("  nontrivial samples:")
        for r in cyc_s["nontrivial_verified"][:8]:
            print(f"    n={r['n']} bits={r['bits']} A={r['A']}")
    else:
        print("  keine verifizierten nichttrivialen Zyklen")
    if cyc_s["trivial_421"]:
        print(f"  trivial 4-2-1 Treffer: {len(cyc_s['trivial_421'])} (Auszug)")
        for r in cyc_s["trivial_421"][:4]:
            print(f"    n={r['n']} bits={r['bits']}")

    print("\n[7] Zyklussuche (elementare Bits, k≤{})".format(min(k_max, 14)))
    cyc_e = scan_cycles(min(k_max, 14), mode="elementary")
    print(
        f"  verified total={cyc_e['verified_total']}  "
        f"nontrivial={cyc_e['nontrivial_count']}"
    )
    for r in cyc_e["nontrivial_verified"][:6]:
        print(f"    n={r['n']} bits={r['bits']} A={r['A']}")

    print("\n[8] Strukturelle Beobachtung")
    print(
        "  Bei Syracuse-Affine ist stets c≥0, also B≥0.\n"
        "  Erster Abstieg T^τ(n)<n erzwingt daher A_τ<1\n"
        "  (sonst An+B≥n). Expandierende Erstabstiege existieren nicht.\n"
        "  V ist positiv auf Erstabstiegen, aber keine Lyapunov-Funktion\n"
        "  entlang der Abstiegskette (Übergang in gerade Klassen hebt V)."
    )

    print("\n[9] Fazit")
    print(
        "  Kanonischer Horizont τ macht (a,τ,c) bahnbezogen.\n"
        "  Atlas = Collatz-Gap-Klassen G={(a,τ,c)}.\n"
        "  Pro Paritätswort ≤1 rationaler Zyklus-Kandidat n=−B/(A−1);\n"
        "  bis k=16 nur der 4-2-1-Zyklus verifiziert."
    )

    out = Path(__file__).resolve().parent / "collatz_resonance_atlas.json"
    payload = {
        "claim": (
            "Collatz resonance atlas via first-descent affine class (a,τ,c); "
            "Lyapunov candidate V=-log A; cycle search n=-B/(A-1); not a proof"
        ),
        "crit": CRIT,
        "atlas": atlas,
        "lyapunov": lya,
        "cycles_syracuse": {
            "k_max": cyc_s["k_max"],
            "verified_total": cyc_s["verified_total"],
            "nontrivial_count": cyc_s["nontrivial_count"],
            "by_k": [
                {"k": r["k"], "candidates": r["candidates"], "verified": r["verified"]}
                for r in cyc_s["by_k"]
            ],
            "nontrivial_verified": cyc_s["nontrivial_verified"],
            "trivial_421_sample": cyc_s["trivial_421"][:10],
        },
        "cycles_elementary": {
            "k_max": cyc_e["k_max"],
            "verified_total": cyc_e["verified_total"],
            "nontrivial_count": cyc_e["nontrivial_count"],
            "nontrivial_verified": cyc_e["nontrivial_verified"],
        },
        "definition": {
            "tau": "min{k>0: T^k(n)<n}",
            "rho": "|c|/(2^tau + 3^a)",
            "C": "(a, tau, rho)",
            "V": "tau*log2 - a*log3 = -log A",
            "cycle": "n = -B/(A-1)",
        },
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
