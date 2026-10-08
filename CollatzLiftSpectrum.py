"""
CollatzLiftSpectrum — Kartierung (Λ(ω), density(ℋ(ω))).

Struktur:
  ω → T_ω(n) = (3^a n + B)/2^L
  Cons(j) ⇔ word_L(x_ω(j))=ω ∧ S^L(x)≡x (mod 2^j)
  ℋ(ω) = { j : Cons(j) }          (Liftmenge, nicht monoton)
  N_ω(J) = |ℋ ∩ [1,J]|
  dens_ω(J) = N_ω(J)/J

Typen:
  A — stabile Kanäle: Λ<0 und dens→1 (ω=10 im Scan)
  B — Phantomkanäle:  Λ>0 und dens→0 (|ℋ|<∞, ggf. sporadisch)

Vermutung (stärker als Zyklusformulierung):
  Λ>0  ⇒  |ℋ(ω)| < ∞  (2-adische Maßnullmenge der konsistenten Ebenen)
  nicht zwingend Max-Formulierung; sporadische ℋ={12,25,37} erlaubt.

Epistemik: Automaten-/Lift-Spektrum; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from itertools import product
from pathlib import Path


LOG2_3 = math.log2(3)
P_CRIT = math.log(2) / math.log(3)


def modinv_odd(a: int, mod: int) -> int:
    return pow(a, -1, mod)


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


def algebraic_fixed_point(a: int, L: int, B: int, j: int) -> int:
    mod = 1 << j
    coef = ((1 << L) - pow(3, a)) % mod
    return (B * modinv_odd(coef, mod)) % mod


def lift_set(omega: tuple[int, ...], j_max: int) -> dict:
    """ℋ(ω)∩[1,J] und Dichte-Profil."""
    a, L, B = affine_of_word(omega)
    if L == 0:
        return {"omega": "", "a": 0, "L": 0, "Lambda": 0.0, "H": [], "N": 0, "dens": 0.0}

    Lam = a * LOG2_3 - L
    H: list[int] = []
    for j in range(1, j_max + 1):
        mod = 1 << j
        x = algebraic_fixed_point(a, L, B, j)
        w, x_end = syracuse_iterate_close(x, L, mod)
        if w == omega and x_end == x:
            H.append(j)

    N = len(H)
    dens = N / j_max
    # Asymptotische Dichte ab j≥2 (j=1 oft zu klein für echte Zyklen)
    H_ge2 = [j for j in H if j >= 2]
    span_ge2 = max(j_max - 1, 1)
    dens_ge2 = len(H_ge2) / span_ge2
    # cumulative density curve at checkpoints
    checkpoints = [j for j in (8, 16, 24, 32, 40, 48, j_max) if j <= j_max]
    dens_curve = {}
    for J in checkpoints:
        dens_curve[J] = sum(1 for j in H if j <= J) / J
        dens_curve[f"ge2@{J}"] = sum(1 for j in H if 2 <= j <= J) / max(J - 1, 1)

    gaps = [H[i + 1] - H[i] for i in range(len(H) - 1)] if len(H) >= 2 else []
    contiguous = N > 0 and H == list(range(H[0], H[-1] + 1))
    full_ray = H == list(range(1, j_max + 1))
    # stabiler Kanal: Cons auf allen Ebenen j=2..J (unendlicher Lift im Scan)
    stable_ray = H_ge2 == list(range(2, j_max + 1))

    return {
        "omega": "".join(map(str, omega)),
        "a": a,
        "L": L,
        "B": B,
        "mu": a / L,
        "Lambda": Lam,
        "Lambda_sign": "pos" if Lam > 1e-12 else ("neg" if Lam < -1e-12 else "zero"),
        "H": H,
        "N": N,
        "dens": dens,
        "dens_ge2": dens_ge2,
        "dens_curve": dens_curve,
        "max_H": max(H) if H else 0,
        "min_H": min(H) if H else 0,
        "gaps": gaps,
        "max_gap": max(gaps) if gaps else 0,
        "contiguous": contiguous,
        "full_ray": full_ray,
        "stable_ray": stable_ray,
        "sporadic": N >= 2 and not contiguous,
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


def classify(rec: dict, j_max: int) -> str:
    """Typ A / B / C (kontrahierend, endlicher Lift) / D (leer)."""
    if rec["N"] == 0:
        return "D_empty"
    if rec.get("stable_ray") or rec["full_ray"]:
        return "A_stable"
    if rec["Lambda_sign"] == "pos":
        return "B_phantom"
    return "C_contract_finite"


def main() -> int:
    j_max = 48
    L_max = 12

    print("=" * 72)
    print("LIFT-SPEKTRUM  —  (Λ(ω), dens(ℋ(ω)))")
    print("=" * 72)
    print(f"Cons ⇔ word=ω ∧ S^L(x)≡x;  dens = |ℋ∩[1,{j_max}]|/{j_max}")
    print(f"dens≥2 = |ℋ∩[2,{j_max}]|/({j_max}-1)  — j=1 oft trivial inkonsistent")
    print("Vermutung: Λ>0 ⇒ dens→0; stabile Kanäle ⇒ dens≥2→1\n")

    # [1] Kontrollkanäle
    print("[1] Kontrollkanäle (Typ A / Null)")
    controls = {
        "10": lift_set((1, 0), j_max),
        "0": lift_set((0,), j_max),
        "1100": lift_set((1, 1, 0, 0), j_max),
    }
    for name, rec in controls.items():
        print(
            f"  ω={name:6s}  Λ={rec['Lambda']:+.4f}  N={rec['N']:2d}  "
            f"dens={rec['dens']:.4f}  dens≥2={rec['dens_ge2']:.4f}  "
            f"stable_ray={rec['stable_ray']}  sporadic={rec['sporadic']}  "
            f"typ={classify(rec, j_max)}"
        )
        curve = " ".join(
            f"{J}:{rec['dens_curve'][J]:.2f}"
            for J in sorted(k for k in rec["dens_curve"] if isinstance(k, int))
        )
        print(f"    dens_curve  {curve}")

    # [2] Automaten-Phantome
    print("\n[2] Phantomkanäle aus G_j (Typ B)")
    phantoms: list[dict] = []
    seen: set[str] = set()
    for j in (8, 10, 11, 12):
        for omega in cycles_as_words(j):
            key = "".join(map(str, omega))
            if key in seen:
                continue
            rec = lift_set(omega, j_max)
            if rec["Lambda_sign"] != "pos":
                continue
            seen.add(key)
            rec["source_j"] = j
            phantoms.append(rec)
            Hshow = rec["H"] if len(rec["H"]) <= 10 else rec["H"][:8] + ["…"]
            print(
                f"  j={j}  ω={rec['omega'][:36]}{'…' if len(rec['omega'])>36 else ''}  "
                f"Λ={rec['Lambda']:+.3f}  N={rec['N']}  dens={rec['dens']:.4f}  "
                f"ℋ={Hshow}  sporadic={rec['sporadic']}"
            )

    # [3] Exhaustive short words — Spektrum
    print(f"\n[3] Spektrum aller Wörter L≤{L_max}, a≥1  (J={j_max})")
    rows: list[dict] = []
    for L in range(1, L_max + 1):
        for bits in product([0, 1], repeat=L):
            if bits.count(1) == 0:
                continue
            rows.append(lift_set(bits, j_max))

    by_sign: dict[str, list[dict]] = defaultdict(list)
    by_typ: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_sign[r["Lambda_sign"]].append(r)
        by_typ[classify(r, j_max)].append(r)

    def dens_stats(group: list[dict]) -> dict:
        if not group:
            return {"n": 0}
        ds = [r["dens"] for r in group]
        d2 = [r["dens_ge2"] for r in group]
        Ns = [r["N"] for r in group]
        return {
            "n": len(group),
            "dens_min": min(ds),
            "dens_max": max(ds),
            "dens_mean": sum(ds) / len(ds),
            "dens_ge2_max": max(d2),
            "dens_ge2_mean": sum(d2) / len(d2),
            "N_max": max(Ns),
            "n_stable_ray": sum(1 for r in group if r["stable_ray"]),
            "n_full_ray": sum(1 for r in group if r["full_ray"]),
            "n_sporadic": sum(1 for r in group if r["sporadic"]),
            "n_empty": sum(1 for r in group if r["N"] == 0),
        }

    for sign in ("pos", "neg", "zero"):
        st = dens_stats(by_sign[sign])
        print(
            f"  Λ {sign:3s}: n={st.get('n',0):5d}  dens∈[{st.get('dens_min',0):.4f},"
            f"{st.get('dens_max',0):.4f}]  mean={st.get('dens_mean',0):.4f}  "
            f"dens≥2_max={st.get('dens_ge2_max',0):.4f}  "
            f"N_max={st.get('N_max',0)}  stable_ray={st.get('n_stable_ray',0)}  "
            f"sporadic={st.get('n_sporadic',0)}  empty={st.get('n_empty',0)}"
        )

    print("\n  Typverteilung:")
    for typ in ("A_stable", "B_phantom", "C_contract_finite", "D_empty"):
        g = by_typ[typ]
        st = dens_stats(g)
        print(
            f"    {typ:20s}  n={st.get('n',0):5d}  "
            f"dens_mean={st.get('dens_mean',0):.4f}  N_max={st.get('N_max',0)}"
        )

    # [4] Spektrallücke: dens vs Λ
    print("\n[4] Spektrallücke dens|Λ>0 vs dens|Λ<0")
    pos = by_sign["pos"]
    neg = by_sign["neg"]
    dens_pos_max = max((r["dens"] for r in pos), default=0.0)
    dens_pos_ge2_max = max((r["dens_ge2"] for r in pos), default=0.0)
    dens_neg_stable = [r for r in neg if r["stable_ray"]]
    dens_neg_max_unstable = max(
        (r["dens_ge2"] for r in neg if not r["stable_ray"]), default=0.0
    )
    print(f"  max dens bei Λ>0:           {dens_pos_max:.4f}")
    print(f"  max dens≥2 bei Λ>0:         {dens_pos_ge2_max:.4f}  (soll →0)")
    print(f"  #Λ<0 mit stable_ray dens≥2=1: {len(dens_neg_stable)}")
    print(f"  max dens≥2 bei Λ<0, ¬stable:  {dens_neg_max_unstable:.4f}")
    print(
        f"  Trennung: alle Λ>0 haben dens≥2<1? "
        f"{all(r['dens_ge2'] < 1.0 - 1e-15 for r in pos)}"
    )
    print(
        f"  Λ>0 ⇒ |ℋ|<{j_max}: "
        f"{all(r['N'] < j_max for r in pos)}  "
        f"(max N={max((r['N'] for r in pos), default=0)})"
    )

    # [5] Top dens bei Λ>0 (worst phantoms) und full Λ<0
    print("\n[5] Extreme Paare (Λ, dens)")
    print("  höchste dens≥2 bei Λ>0:")
    for r in sorted(pos, key=lambda d: (-d["dens_ge2"], -d["N"], d["L"]))[:8]:
        print(
            f"    dens≥2={r['dens_ge2']:.4f}  N={r['N']:2d}  Λ={r['Lambda']:+.3f}  "
            f"L={r['L']:2d}  sporadic={r['sporadic']}  "
            f"ℋ={r['H'] if len(r['H'])<=8 else r['H'][:6]+['…']}"
            f"  ω={r['omega']}"
        )
    print("  stabile Liftkanäle Λ<0 (dens≥2=1, Cons auf j=2..J):")
    for r in sorted(dens_neg_stable, key=lambda d: (d["L"], d["omega"]))[:12]:
        print(
            f"    dens≥2=1  Λ={r['Lambda']:+.4f}  L={r['L']:2d} a={r['a']:2d}  "
            f"ω={r['omega']}"
        )
    if len(dens_neg_stable) > 12:
        print(f"    … +{len(dens_neg_stable)-12} weitere")

    # [6] Selbstähnlichkeit: dens_curve für Phantome vs Kontrolle
    print("\n[6] Lift-Selbstähnlichkeit N(J)/J")
    print("  Kontrolle ω=10:")
    c10 = controls["10"]
    for J in sorted(k for k in c10["dens_curve"] if isinstance(k, int)):
        d = c10["dens_curve"][J]
        d2 = c10["dens_curve"][f"ge2@{J}"]
        print(f"    J={J:2d}  dens={d:.4f}  dens≥2={d2:.4f}")
    if phantoms:
        print("  Phantome (mittlere dens≥2-Kurve):")
        keys = sorted(k for k in phantoms[0]["dens_curve"] if isinstance(k, int))
        for J in keys:
            vals = [p["dens_curve"][f"ge2@{J}"] for p in phantoms]
            print(
                f"    J={J:2d}  mean dens≥2={sum(vals)/len(vals):.4f}  "
                f"max={max(vals):.4f}"
            )

    # [7] Fazit
    print("\n[7] Fazit")
    print(
        "  Resonanzgröße: (Λ, dens(ℋ)); ℋ nicht monoton (sporadische Phantome).\n"
        "  Λ>0 ⇒ |ℋ|<∞, dens→0 — Maßnullmenge konsistenter Ebenen.\n"
        "  Stabile Kanäle: dens≥2→1 (ω=10 fehlt nur j=1 trivial).\n"
        "  Λ<0 allein ⇒ nicht dens→1 — die meisten Wörter haben ℋ=∅."
    )

    # Scatter samples for JSON (compact)
    scatter = [
        {
            "omega": r["omega"],
            "Lambda": r["Lambda"],
            "dens": r["dens"],
            "N": r["N"],
            "L": r["L"],
            "a": r["a"],
            "typ": classify(r, j_max),
            "sporadic": r["sporadic"],
            "H": r["H"] if r["N"] <= 12 else r["H"][:12] + ["…"],
        }
        for r in sorted(rows, key=lambda d: (-d["dens"], d["Lambda"]))
        if r["N"] > 0 or r["Lambda_sign"] == "pos"
    ]
    # keep manageable: all pos + all full + top dens
    scatter_out = []
    seen_sc: set[str] = set()
    for r in rows:
        keep = (
            r["Lambda_sign"] == "pos"
            or r["stable_ray"]
            or r["sporadic"]
            or (r["N"] >= 4 and r["dens_ge2"] >= 0.05)
        )
        if keep and r["omega"] not in seen_sc:
            seen_sc.add(r["omega"])
            scatter_out.append(
                {
                    "omega": r["omega"],
                    "Lambda": r["Lambda"],
                    "dens": r["dens"],
                    "dens_ge2": r["dens_ge2"],
                    "N": r["N"],
                    "L": r["L"],
                    "a": r["a"],
                    "typ": classify(r, j_max),
                    "sporadic": r["sporadic"],
                    "stable_ray": r["stable_ray"],
                    "H": r["H"],
                    "dens_curve": {
                        str(k): v for k, v in r["dens_curve"].items()
                    },
                }
            )

    out = Path(__file__).resolve().parent / "collatz_lift_spectrum.json"
    payload = {
        "claim": (
            "Lift spectrum: pair (Lambda, dens(H)); Lambda>0 => finite lift set "
            "(density->0); stable channels dens->1; H not necessarily monotone; "
            "not a Collatz proof"
        ),
        "j_max": j_max,
        "L_max": L_max,
        "p_crit": P_CRIT,
        "conjecture": {
            "Lambda_pos": "Lambda>0 => |H(omega)| < infinity (dens->0)",
            "stable": "attractors: Lambda<0 and dens_ge2->1 (Cons on j=2..J)",
            "note": "Lambda<0 alone does not imply dens->1; j=1 often trivial fail",
        },
        "controls": controls,
        "phantoms": phantoms,
        "sign_stats": {k: dens_stats(v) for k, v in by_sign.items()},
        "type_stats": {k: dens_stats(v) for k, v in by_typ.items()},
        "gap": {
            "max_dens_Lambda_pos": dens_pos_max,
            "max_dens_ge2_Lambda_pos": dens_pos_ge2_max,
            "n_Lambda_neg_stable_ray": len(dens_neg_stable),
            "max_N_Lambda_pos": max((r["N"] for r in pos), default=0),
            "all_Lambda_pos_N_lt_J": all(r["N"] < j_max for r in pos),
        },
        "scatter": scatter_out,
        "stable_ray_words": [r["omega"] for r in dens_neg_stable],
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
