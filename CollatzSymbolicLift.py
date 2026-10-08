"""
CollatzSymbolicLift — symbolische Lift-Obstruktion (nicht Valuation).

Korrektur: v_2(3^a−2^L)=0 gilt für jedes a≥1, L≥1 (Parität),
also keine spezielle Resonanz-Obstruktion.

Lift hat zwei Stufen:
  1) Algebraischer Fixpunkt der Affinform T_ω:
       (2^L − 3^a) n ≡ B_ω  (mod 2^j)
     Da 2^L−3^a ungerade: immer eindeutige Lösung x_ω(j).
  2) Automaten-Konsistenz (beide Bedingungen):
       word_L(x_ω(j)) = ω  und  S^L(x_ω(j)) ≡ x_ω(j)  (mod 2^j)
     Nur Wort-Match reicht nicht: algebraischer FP kann Paritäten
     treffen, ohne dass der Orbit im Automaten schließt.

h(ω) = max { j : Cons(j) }

Vermutung:
  Λ(ω) > 0  ⇒  h(ω) < ∞
  (kein stabiler Paritätslift / keine unendliche Zylinderkette)

Kontrollklasse {1,2}: ω=(1,0), Λ<0, h=∞ (im Scan).

Epistemik: symbolische Selbstkonsistenz; kein Collatz-Beweis.
"""

from __future__ import annotations

import json
import math
from itertools import product
from pathlib import Path


LOG2_3 = math.log2(3)
P_CRIT = math.log(2) / math.log(3)


def modinv_odd(a: int, mod: int) -> int:
    """Inverses von ungeradem a modulo 2^k (mod=2^k)."""
    # pow built-in works for odd a mod 2^k
    return pow(a, -1, mod)


def affine_of_word(omega: tuple[int, ...]) -> tuple[int, int, int]:
    """(a, L, B) mit T(n)=(3^a n+B)/2^L."""
    a = 0
    L = 0
    b = 0
    for e in omega:
        if e == 0:
            L += 1
        else:
            b = 3 * b + (1 << L)
            a += 1
            L += 1
    return a, L, b


def syracuse_word(n: int, L: int, mod: int | None = None) -> tuple[int, ...]:
    """Paritätswort der Länge L; n ggf. als Rest mod 2^m repräsentieren."""
    bits = []
    x = n if mod is None else n % mod
    for _ in range(L):
        if mod is not None:
            x %= mod
        bit = 0 if x % 2 == 0 else 1
        bits.append(bit)
        if bit == 0:
            x //= 2
        else:
            x = (3 * x + 1) // 2
        if mod is not None:
            x %= mod
    return tuple(bits)


def syracuse_iterate_close(n: int, L: int, mod: int) -> tuple[tuple[int, ...], int]:
    """Wort + Endzustand nach L Automaten-Schritten mod 2^j."""
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
    """
    x = B * (2^L − 3^a)^{-1}  mod 2^j.
    Äquivalent: (2^L − 3^a) x ≡ B (mod 2^j).
    """
    mod = 1 << j
    coef = ((1 << L) - pow(3, a)) % mod
    # coef ≡ -3^a ≡ odd (da 3^a odd)
    if coef % 2 == 0:
        # should not happen for a≥0,L≥1 with meaningful cycles; handle a=0
        coef = (coef + mod) % mod
    inv = modinv_odd(coef, mod)
    return (B * inv) % mod


def consistency_height(
    omega: tuple[int, ...],
    j_max: int = 48,
) -> dict:
    """
    Konsistenz Cons(j): word_L(x_ω(j)) = ω.

    Wichtig: Cons ist NICHT notwendig ein Anfangsintervall —
    Phantome können bei kleinem j scheitern und bei j=source gelten.

    h(ω)  = max { j ≤ j_max : Cons(j) }   (0 falls nie)
    h_ray = max H s.t. Cons(j) ∀ j∈[j_lo, H] für größtes Intervall am oberen Ende
    """
    a, L, B = affine_of_word(omega)
    if L == 0:
        return {"omega": "", "h": 0, "a": 0, "L": 0, "Lambda": 0.0}

    Lam = a * LOG2_3 - L
    mu = a / L if L else 0.0
    eta = mu - 0.5

    diff = (1 << L) - 3**a
    nu = 0
    t = abs(diff)
    while t % 2 == 0 and t > 0:
        t //= 2
        nu += 1

    ok_js: list[int] = []
    details_fail = None
    last_ok = None
    for j in range(1, j_max + 1):
        mod = 1 << j
        x = algebraic_fixed_point(a, L, B, j)
        w, x_end = syracuse_iterate_close(x, L, mod)
        # Echte Automaten-Konsistenz: Wort + Schließen S^L(x)≡x
        ok = (w == omega) and (x_end == x)
        if ok:
            ok_js.append(j)
            last_ok = {"j": j, "x": x}
        elif details_fail is None:
            details_fail = {
                "j": j,
                "x": x,
                "x_end": x_end,
                "closed": x_end == x,
                "word": "".join(map(str, w)),
                "expected": "".join(map(str, omega)),
                "word_match": w == omega,
            }

    h = max(ok_js) if ok_js else 0
    # längster Strahl nach oben: von h abwärts solange konsistent
    h_lo = h
    while h_lo - 1 in ok_js:
        h_lo -= 1
    contiguous_top = h >= 1 and set(range(h_lo, h + 1)).issubset(set(ok_js))

    return {
        "omega": "".join(map(str, omega)),
        "a": a,
        "L": L,
        "B": B,
        "mu": mu,
        "eta": eta,
        "Lambda": Lam,
        "nu_abs_diff": nu,
        "h": h,
        "h_lo": h_lo if ok_js else 0,
        "ok_js": ok_js,
        "n_ok": len(ok_js),
        "h_finite": h < j_max,
        "contiguous_top": contiguous_top,
        "last_ok": last_ok,
        "first_fail_in_scan": details_fail,
        "Lambda_pos": Lam > 1e-12,
        "eta_pos": eta > 1e-15,
    }


def cycles_as_words(j: int) -> list[tuple[int, ...]]:
    """Extrahiere Zykluswörter (eine Rotation pro Zyklus: Start = min node)."""
    mod = 1 << j
    nxt = []
    for r in range(mod):
        if r % 2 == 0:
            nxt.append((r // 2) % mod)
        else:
            nxt.append(((3 * r + 1) // 2) % mod)

    visited = [0] * mod
    words = []
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
            # word from smallest node
            i0 = min(range(len(cyc)), key=lambda t: cyc[t])
            rot = cyc[i0:] + cyc[:i0]
            bits = []
            for v in rot:
                bits.append(0 if v % 2 == 0 else 1)
            words.append(tuple(bits))
        for v in path:
            visited[v] = 2
    return words


def main() -> int:
    print("=" * 72)
    print("SYMBOLISCHER LIFT  —  h(ω) = max konsistente 2^j-Höhe")
    print("=" * 72)
    print("x_ω = B(2^L−3^a)^{-1} mod 2^j;  Cons ⇔ word=ω ∧ S^L(x)≡x")
    print("Korrektur: v_2(3^a−2^L)=0 immer (a≥1) — nicht die Obstruktion.\n")

    j_max = 48

    # [1] Sanity valuation
    print("[1] Valuation-Sanity (kein Resonanzsignal)")
    for omega in [(1, 0), (1, 1, 0, 0), (1,) * 5 + (0,) * 5]:
        a, L, B = affine_of_word(omega)
        diff = abs((1 << L) - 3**a)
        nu = 0
        t = diff
        while t % 2 == 0 and t:
            t //= 2
            nu += 1
        print(f"  ω={''.join(map(str,omega))}  a={a} L={L}  ν={nu}  (erwartet 0 falls a≥1)")

    # [2] Kontrollklasse 1-2
    print("\n[2] Kontrollklasse ω=10 (Zyklus 1↔2)")
    ctrl = consistency_height((1, 0), j_max=j_max)
    print(
        f"  Λ={ctrl['Lambda']:+.4f}  h={ctrl['h']}  "
        f"finite={ctrl['h_finite']}  last_ok={ctrl['last_ok']}"
    )

    # [3] Zero cycle word "0"
    print("\n[3] Kontrollklasse ω=0 (Fixpunkt 0)")
    z = consistency_height((0,), j_max=j_max)
    print(f"  Λ={z['Lambda']:+.4f}  h={z['h']}  finite={z['h_finite']}")

    # [4] Phantom cycle words from automata
    print("\n[4] Phantom-Zyklen aus G_j (symbolische Höhe h)")
    phantom_js = [8, 10, 11, 12]
    phantom_rows = []
    for j in phantom_js:
        for omega in cycles_as_words(j):
            rec = consistency_height(omega, j_max=j_max)
            if not rec["eta_pos"] and not rec["Lambda_pos"]:
                continue
            # also keep if it's the known phantom (Lambda>0)
            if not rec["Lambda_pos"]:
                continue
            rec["source_j"] = j
            phantom_rows.append(rec)
            print(
                f"  source j={j}  ω={rec['omega'][:40]}{'…' if len(rec['omega'])>40 else ''}  "
                f"L={rec['L']} a={rec['a']}  Λ={rec['Lambda']:+.3f}  "
                f"h={rec['h']}  ok_js={rec['ok_js'][:12]}{'…' if len(rec['ok_js'])>12 else ''}"
            )
            # Konsistenz genau auf source-Ebene?
            print(
                f"    Cons({j})={j in rec['ok_js']}  "
                f"Cons({j}+1)={ (j+1) in rec['ok_js']}  "
                f"x@h={rec['last_ok']}"
            )

    # Dedup by omega
    by_omega = {r["omega"]: r for r in phantom_rows}

    # [5] Random / exhaustive short words with Λ>0
    print("\n[5] Alle Wörter L≤12 mit Λ>0 — Höhenstatistik")
    lam_pos = []
    lam_neg = []
    for L in range(1, 13):
        for bits in product([0, 1], repeat=L):
            if bits.count(1) == 0:
                continue
            rec = consistency_height(bits, j_max=j_max)
            if rec["Lambda_pos"]:
                lam_pos.append(rec)
            else:
                lam_neg.append(rec)

    h_pos = [r["h"] for r in lam_pos]
    h_neg_sample = [r["h"] for r in lam_neg if r["h"] < j_max]  # finite ones
    h_neg_inf = sum(1 for r in lam_neg if r["h"] >= j_max)
    print(f"  #Wörter Λ>0 (L≤12, a≥1): {len(lam_pos)}")
    if h_pos:
        print(
            f"  h|Λ>0:  min={min(h_pos)}  max={max(h_pos)}  "
            f"mean={sum(h_pos)/len(h_pos):.2f}  alle endlich? {max(h_pos)<j_max}"
        )
    print(f"  #Wörter Λ≤0: {len(lam_neg)}  davon h≥{j_max}: {h_neg_inf}")

    print("\n  größte h bei Λ>0:")
    for r in sorted(lam_pos, key=lambda d: (-d["h"], d["L"]))[:10]:
        print(
            f"    h={r['h']:2d}  n_ok={r['n_ok']:2d}  L={r['L']:2d} a={r['a']:2d}  "
            f"Λ={r['Lambda']:+.3f}  ω={r['omega']}"
        )

    # [6] Implication check
    print("\n[6] Vermutung Λ>0 ⇒ h<∞ (kein unendlicher symbolischer Lift)")
    all_finite = all(r["h"] < j_max for r in lam_pos)
    phantoms_finite = all(r["h"] < j_max for r in by_omega.values())
    # stärkere Form: kein Cons(j) für j > L + C o.ä.; hier: max h
    max_h_pos = max(h_pos) if h_pos else 0
    print(f"  alle Λ>0-Wörter L≤12: h<{j_max} ? {all_finite}  (max h={max_h_pos})")
    print(f"  Automaten-Phantome:     h<{j_max} ? {phantoms_finite}")
    if by_omega:
        for w, r in by_omega.items():
            print(
                f"    Phantom h={r['h']} ok_js={r['ok_js']}  "
                f"(Konsistenz nur auf endlicher Träger-Menge)"
            )
    print(
        f"  Kontroll ω=10: h={ctrl['h']} ok_js=[1..{j_max}] "
        f"⇒ unendlicher Lift im Scan"
    )

    print("\n[7] Fazit")
    print(
        "  Valuation ν=0 ist universell — nicht die Sterilisierungsursache.\n"
        "  Algebraischer FP x_ω existiert stets in Z/2^j (ungerader Koeffizient).\n"
        "  Engpass: Automaten-Konsistenz word(x)=ω UND S^L(x)≡x.\n"
        "  Bei Λ>0 oft: Paritäten passen kurz, aber der Orbit schließt nicht\n"
        "  auf höherem Level — daher endliches h. Attraktor ω=10: h=∞ im Scan."
    )

    out = Path(__file__).resolve().parent / "collatz_symbolic_lift.json"
    payload = {
        "claim": (
            "Symbolic lift obstruction: algebraic FP always exists (odd coef), "
            "consistency word(x)=omega fails at finite h when Lambda>0; "
            "valuation is not the obstruction; not a Collatz proof"
        ),
        "p_crit": P_CRIT,
        "j_max": j_max,
        "valuation_note": "v_2(3^a-2^L)=0 for all a>=1,L>=1",
        "control_10": ctrl,
        "control_0": z,
        "automaton_phantoms": list(by_omega.values()),
        "words_L_le_12": {
            "n_Lambda_pos": len(lam_pos),
            "h_min": min(h_pos) if h_pos else None,
            "h_max": max(h_pos) if h_pos else None,
            "h_mean": sum(h_pos) / len(h_pos) if h_pos else None,
            "all_h_finite": all_finite,
            "n_Lambda_nonpos": len(lam_neg),
            "n_nonpos_h_infinite_scan": h_neg_inf,
            "shortest_h_samples": sorted(
                lam_pos, key=lambda d: (d["h"], d["L"])
            )[:20],
        },
        "conjecture": "Lambda>0 => h(omega)<infty (no stable parity lift)",
    }
    out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
