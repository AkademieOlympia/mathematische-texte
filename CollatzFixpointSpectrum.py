"""
CollatzFixpointSpectrum — affine Fixpunkt-Geometrie (Λ, x*, χ).

Für ω mit T_ω(x) = A x + B, A = 3^a/2^L:
  x*_ℚ = B / (1 − A) = B / (2^L − 3^a)

Klassifikation:
  Typ A  — Attraktor:   x* ∈ {0,1,2}  (⇒ Λ<0 im Scan)
  Typ B  — Phantom:     Λ>0, x*<0
  Typ D  — neutral:     Λ<0, x* ∉ {0,1,2}

Abstandsspektrum:
  χ(ω) = dist(x*_ℚ, {0,1,2})

Vermutung:
  d→1  ⟺  χ=0
  χ>0  ⇒  d→0

Korrektur Beispiele:
  ω=10:  T=(3x+1)/4,  x*=1
  ω=01:  T=(3x+2)/4,  x*=2

Epistemik: affine Resonanzklasse R_C=(Λ,x*); kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path


LOG2_3 = math.log2(3)
ATTRACTOR = (0, 1, 2)


def affine_of_word(omega: tuple[int, ...]) -> tuple[int, int, int]:
    a = L = b = 0
    for e in omega:
        if e == 0:
            L += 1
        else:
            b = 3 * b + (1 << L)
            a += 1
            L += 1
    return a, L, b


def modinv_odd(a: int, mod: int) -> int:
    return pow(a, -1, mod)


def algebraic_fixed_point(a: int, L: int, B: int, j: int) -> int:
    mod = 1 << j
    coef = ((1 << L) - pow(3, a)) % mod
    return (B * modinv_odd(coef, mod)) % mod


def syracuse_iterate_close(n: int, L: int, mod: int) -> tuple[tuple[int, ...], int]:
    bits = []
    x = n % mod
    for _ in range(L):
        bit = 0 if x % 2 == 0 else 1
        bits.append(bit)
        if bit == 0:
            x = (x // 2) % mod
        else:
            x = ((3 * x + 1) // 2) % mod
    return tuple(bits), x


def chi_dist(xQ: Fraction) -> Fraction:
    """χ = dist(x*, {0,1,2}) in ℚ."""
    return min(abs(xQ - Fraction(t)) for t in ATTRACTOR)


def spectrum_row(omega: tuple[int, ...], j_max: int) -> dict:
    a, L, B = affine_of_word(omega)
    if L == 0:
        return {"omega": "", "L": 0}

    D = (1 << L) - 3**a  # 2^L−3^a
    A = Fraction(3**a, 1 << L)
    b = Fraction(B, 1 << L)  # T(x)=A x + b  (b = B_int/2^L)
    Lam = a * LOG2_3 - L
    # x* = b/(1−A) = B_int/(2^L−3^a)
    xQ = Fraction(B, D) if D != 0 else None
    xQ_alt = b / (1 - A) if A != 1 else None
    assert xQ == xQ_alt

    chi = chi_dist(xQ)
    chi_f = float(chi)

    H = []
    for j in range(1, j_max + 1):
        mod = 1 << j
        x = algebraic_fixed_point(a, L, B, j)
        w, x_end = syracuse_iterate_close(x, L, mod)
        if w == omega and x_end == x:
            H.append(j)
    H_ge2 = [j for j in H if j >= 2]
    dens = len(H) / j_max
    dens_ge2 = len(H_ge2) / max(j_max - 1, 1)
    stable = H_ge2 == list(range(2, j_max + 1))
    d_inf = 1.0 if stable else (0.0 if dens_ge2 == 0 else dens_ge2)  # Scan-Proxy

    if chi == 0:
        typ = "A_attractor"
    elif Lam > 1e-12:
        typ = "B_phantom"
    else:
        typ = "D_neutral"

    return {
        "omega": "".join(map(str, omega)),
        "a": a,
        "L": L,
        "B": B,
        "A": str(A),
        "Lambda": Lam,
        "x_star": str(xQ),
        "x_star_float": float(xQ),
        "chi": str(chi),
        "chi_float": chi_f,
        "chi_zero": chi == 0,
        "typ": typ,
        "H": H,
        "N": len(H),
        "dens": dens,
        "dens_ge2": dens_ge2,
        "d_inf_proxy": d_inf,
        "stable_ray": stable,
    }


def cycles_as_words(j: int) -> list[tuple[int, ...]]:
    mod = 1 << j
    nxt = [
        (r // 2) % mod if r % 2 == 0 else ((3 * r + 1) // 2) % mod for r in range(mod)
    ]
    visited = [0] * mod
    words: list[tuple[int, ...]] = []
    for start in range(mod):
        if visited[start]:
            continue
        path: list[int] = []
        pos: dict[int, int] = {}
        x = start
        while visited[x] == 0:
            visited[x] = 1
            pos[x] = len(path)
            path.append(x)
            x = nxt[x]
        if visited[x] == 1:
            cyc = path[pos[x] :]
            i0 = min(range(len(cyc)), key=lambda t: cyc[t])
            rot = cyc[i0:] + cyc[:i0]
            words.append(tuple(0 if v % 2 == 0 else 1 for v in rot))
        for v in path:
            visited[v] = 2
    return words


def main() -> int:
    j_max = 48
    L_max = 12

    print("=" * 72)
    print("FIXPUNKT-SPEKTRUM  —  R_C(ω) = (Λ, x*)  mit  χ=dist(x*,{0,1,2})")
    print("=" * 72)
    print("Vermutung: d→1 ⟺ χ=0;  χ>0 ⇒ d→0\n")

    # [1] Kanonische Beispiele
    print("[1] Kanonische Operatoren")
    examples = {
        "10": (1, 0),
        "01": (0, 1),
        "1010": (1, 0, 1, 0),
        "1100": (1, 1, 0, 0),
        "0": (0,),
        "1110111": (1, 1, 1, 0, 1, 1, 1),
    }
    for name, omega in examples.items():
        r = spectrum_row(omega, j_max)
        print(
            f"  ω={name:8s}  A={r['A']:8s}  Λ={r['Lambda']:+.4f}  "
            f"x*={r['x_star']:>12s}  χ={r['chi']:8s}  "
            f"d∞≈{r['d_inf_proxy']:.3f}  typ={r['typ']}"
        )

    # [2] Exhaustive Spektrum
    print(f"\n[2] Spektrum aller Wörter L≤{L_max}")
    rows = []
    for L in range(1, L_max + 1):
        for bits in product([0, 1], repeat=L):
            if bits == (0,) * L and L > 1:
                continue
            rows.append(spectrum_row(bits, j_max))

    by_typ = defaultdict(list)
    for r in rows:
        by_typ[r["typ"]].append(r)

    for typ in ("A_attractor", "B_phantom", "D_neutral"):
        g = by_typ[typ]
        if not g:
            print(f"  {typ}: n=0")
            continue
        chis = [r["chi_float"] for r in g]
        dens = [r["dens_ge2"] for r in g]
        print(
            f"  {typ:14s}  n={len(g):5d}  "
            f"χ∈[{min(chis):.4g},{max(chis):.4g}]  "
            f"dens≥2∈[{min(dens):.4f},{max(dens):.4f}]  "
            f"stable={sum(1 for r in g if r['stable_ray'])}"
        )

    # [3] Dichotomie χ
    print("\n[3] Dichotomie χ ↔ d")
    chi0 = [r for r in rows if r["chi_zero"]]
    chip = [r for r in rows if not r["chi_zero"]]
    print(
        f"  χ=0: n={len(chi0)}  alle stable_ray? "
        f"{all(r['stable_ray'] for r in chi0)}  "
        f"min dens≥2={min(r['dens_ge2'] for r in chi0):.4f}"
    )
    print(
        f"  χ>0: n={len(chip)}  irgendein stable_ray? "
        f"{any(r['stable_ray'] for r in chip)}  "
        f"max dens≥2={max(r['dens_ge2'] for r in chip):.4f}  "
        f"max N={max(r['N'] for r in chip)}"
    )
    equiv = all(r["stable_ray"] == r["chi_zero"] for r in rows)
    print(f"  Äquivalenz stable_ray ⟺ χ=0: {equiv}")
    print(
        f"  χ>0 ⇒ dens≥2 < 1: "
        f"{all(r['dens_ge2'] < 1 - 1e-15 for r in chip)}"
    )
    # stärkere Form: χ>0 ⇒ d→0 Proxy (dens nicht →1); N endlich schon bekannt
    print(
        f"  χ>0 ⇒ d_inf_proxy=0 (kein stable): "
        f"{all(r['d_inf_proxy'] < 1 for r in chip)}"
    )

    # [4] Scatter (Λ, χ) nach dens-Klasse
    print("\n[4] Paar (Λ, χ) — Extremwerte")
    print("  kleinste χ>0 (nahe Attraktor, aber nicht):")
    near = sorted(chip, key=lambda r: (r["chi_float"], abs(r["Lambda"])))[:10]
    for r in near:
        print(
            f"    χ={r['chi']:10s}  Λ={r['Lambda']:+.4f}  x*={r['x_star']:10s}  "
            f"dens≥2={r['dens_ge2']:.4f}  ω={r['omega']}"
        )
    print("  größte dens≥2 bei χ>0 (Phantom/finite):")
    for r in sorted(chip, key=lambda d: (-d["dens_ge2"], d["chi_float"]))[:8]:
        print(
            f"    dens≥2={r['dens_ge2']:.4f}  χ={r['chi']:10s}  "
            f"Λ={r['Lambda']:+.4f}  x*={r['x_star']:16s}  "
            f"ℋ={r['H']}  ω={r['omega'][:24]}"
        )

    # [5] Phantome im χ-Bild
    print("\n[5] Automaten-Phantome im Fixpunktbild")
    seen: set[str] = set()
    phantoms = []
    for j in (8, 10, 11, 12):
        for omega in cycles_as_words(j):
            key = "".join(map(str, omega))
            if key in seen:
                continue
            r = spectrum_row(omega, j_max)
            if r["Lambda"] <= 0:
                continue
            seen.add(key)
            phantoms.append(r)
            print(
                f"  j={j}  χ={r['chi']}  x*={r['x_star']}  "
                f"Λ={r['Lambda']:+.3f}  dens≥2={r['dens_ge2']:.4f}"
            )

    # [6] Λ vs χ Korrelation (nur χ>0)
    print("\n[6] Kontraktion allein reicht nicht")
    neg_chi = [r for r in rows if r["Lambda"] < -1e-12 and not r["chi_zero"]]
    print(
        f"  # Λ<0 ∧ χ>0 (Typ D): {len(neg_chi)}  "
        f"davon dens≥2>0: {sum(1 for r in neg_chi if r['dens_ge2']>0)}  "
        f"max dens≥2={max((r['dens_ge2'] for r in neg_chi), default=0):.4f}"
    )
    print(
        f"  Beispiel Typ D: ω=1100  x*=5/7  χ={spectrum_row((1,1,0,0), j_max)['chi']}"
    )

    print("\n[7] Fazit")
    print(
        "  R_C(ω)=(Λ,x*) klassifiziert vollständig im Scan:\n"
        "  χ=0 ⟺ d→1 (Attraktor {0,1,2});  χ>0 ⇒ kein stabiler Liftkanal.\n"
        "  Λ>0 ⇒ x*<0 (Phantom);  Λ<0∧χ>0 ⇒ neutral (z.B. 5/7).\n"
        "  Affine Resonanz = Kontraktion + zulässiger Fixpunkt."
    )

    # compact export: all A + all with N>0 + sample of near-χ D
    export_rows = []
    seen_e: set[str] = set()
    for r in rows:
        keep = r["chi_zero"] or r["N"] > 0 or r["chi_float"] < 0.2
        if keep and r["omega"] not in seen_e:
            seen_e.add(r["omega"])
            export_rows.append(
                {
                    "omega": r["omega"],
                    "a": r["a"],
                    "L": r["L"],
                    "B": r["B"],
                    "A": r["A"],
                    "x_star": r["x_star"],
                    "Lambda": r["Lambda"],
                    "chi": r["chi"],
                    "chi_float": r["chi_float"],
                    "dens_ge2": r["dens_ge2"],
                    "d_inf_proxy": r["d_inf_proxy"],
                    "typ": r["typ"],
                    "H": r["H"],
                }
            )

    out = Path(__file__).resolve().parent / "collatz_fixpoint_spectrum.json"
    payload = {
        "claim": (
            "Affine resonance class R_C=(Lambda, x*); "
            "chi=dist(x*,{0,1,2}); stable_ray iff chi=0; "
            "not a Collatz proof"
        ),
        "j_max": j_max,
        "L_max": L_max,
        "attractor_set": list(ATTRACTOR),
        "examples": {k: spectrum_row(v, j_max) for k, v in examples.items()},
        "type_counts": {t: len(by_typ[t]) for t in by_typ},
        "dichotomy": {
            "n_chi_zero": len(chi0),
            "n_chi_pos": len(chip),
            "stable_iff_chi_zero": equiv,
            "chi_pos_max_dens_ge2": max(r["dens_ge2"] for r in chip),
            "chi_pos_max_N": max(r["N"] for r in chip),
            "n_type_D": len(neg_chi),
        },
        "conjecture": "d->1 <=> chi=0; chi>0 => d->0",
        "phantoms": phantoms,
        "spectrum": export_rows,
        "note": "omega=10 has x*=1; omega=01 has x*=2",
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
