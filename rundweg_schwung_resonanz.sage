#!/usr/bin/env sage
# -*- coding: utf-8 -*-
# Quadratischer Schwung u = q̄ q' — Resonanz via N(u) mod 5, mod 7
# Start: ./run_schwung_resonanz.sh  oder  sage rundweg_schwung_resonanz.sage

import csv
import os

from sage.all import Matrix, QQ, ZZ, QuaternionAlgebra, is_prime, factor

# ---------------------------------------------------------------------------
# Quaternionen-Algebra H = QuaternionAlgebra(QQ, -1, -1)
# Basis: 1, i, j, k  —  E,A,B,C ↔ 1,i,j,k
# ---------------------------------------------------------------------------

H = QuaternionAlgebra(QQ, -1, -1)
_i, _j, _k = H.gens()


def q_from_eabc(e, a, b, c):
    """q = e·1 + a·i + b·j + c·k."""
    return H(e) + H(a) * _i + H(b) * _j + H(c) * _k


def eabc_from_q(q):
    t = tuple(list(q))
    return (t[0], t[1], t[2], t[3])


def q_norm(q):
    return q.reduced_norm()


def is_hurwitz(q):
    """Hurwitz-Maximalordnung: ganz/halbzahlig, Summe ∈ ZZ."""
    c = eabc_from_q(q)
    all_int = all(x in ZZ for x in c)
    all_half = all(2 * x in ZZ for x in c) and not all_int
    if not (all_int or all_half):
        return False
    return sum(c) in ZZ


def kollaps_operator_q():
    """Pi_Gamma = (I + Gamma + Gamma^2 + Gamma^3) / 4."""
    R = Matrix(QQ, [[0, 0, 0, 1],
                    [1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0]])
    K = Matrix(QQ, [[0, 0, 0, 1],
                    [0, 0, 1, 0],
                    [0, 1, 0, 0],
                    [1, 0, 0, 0]])
    Gamma = R * K
    I = Matrix.identity(QQ, 4)
    return (I + Gamma + Gamma**2 + Gamma**3) / 4


def kollaps_eabc(e, a, b, c):
    Pi = kollaps_operator_q()
    v = Matrix(QQ, [[e], [a], [b], [c]])
    w = Pi * v
    return tuple(w.list())


def primquaternionen():
    """π_5 = 2+i (Norm 5), π_7 = 2+i+j+k (Norm 7)."""
    pi5 = H(2) + _i
    pi7 = H(2) + _i + _j + _k
    return {5: pi5, 7: pi7}


def norm_notwendig(q, pi):
    nq = q_norm(q)
    npi = q_norm(pi)
    if nq == 0:
        return True
    return (nq % npi == 0) if npi != 0 else False


def in_linksideal(q, pi):
    """q ∈ H·π  ⇔  h = q·π^{-1} ist Hurwitz."""
    if q.is_zero():
        return True, H(0), "trivial"
    if not norm_notwendig(q, pi):
        return False, None, "Norm-Kriterium verletzt"
    h = q * pi**(-1)
    ok = is_hurwitz(h)
    grund = "h=q*pi^-1 ist Hurwitz" if ok else "h=q*pi^-1 nicht Hurwitz"
    return ok, h, grund


def quadratischer_schwung(q, q_prime):
    """u = q̄ · q' (Konjugation via q.conjugate())."""
    return q.conjugate() * q_prime


def norm_mod_p(q, p):
    """N(q) mod p als Integer in [0, p-1]."""
    n = q_norm(q)
    if n == 0:
        return 0
    return int(n % int(p))


def schwung_resonanz_modulo(u, primes=(5, 7)):
    """Resonanz: N(u) ≡ 0 (mod p) und optional u ∈ I_p."""
    pqs = primquaternionen()
    out = {}
    nu = q_norm(u)
    for p in primes:
        pi = pqs[p]
        ok_l, _, grund_l = in_linksideal(u, pi)
        out[p] = {
            "N_mod": norm_mod_p(u, p),
            "N_resonanz": (nu % p == 0),
            "I_p": ok_l,
            "grund": grund_l,
        }
    return out


def schwung_report(q, q_prime, label_q="q", label_qp="q'", *, verbose=True):
    """Einzelpaar: u = q̄ q', Normvergleich und mod-5/7-Resonanz."""
    u = quadratischer_schwung(q, q_prime)
    nq, nqp, nu = q_norm(q), q_norm(q_prime), q_norm(u)
    res = schwung_resonanz_modulo(u)

    if verbose:
        print(f"\n--- Paar {label_q} × {label_qp} ---")
        print(f"  q     = {eabc_from_q(q)},      N(q)  = {nq}")
        print(f"  q'    = {eabc_from_q(q_prime)}, N(q') = {nqp}")
        print(f"  u=q̄q' = {eabc_from_q(u)},      N(u)  = {nu}")
        print(f"  N(u)/N(q) = {nu/nq if nq else '—'},  N(u)/N(q') = {nu/nqp if nqp else '—'}")
        for p in (5, 7):
            r = res[p]
            print(
                f"  mod {p}: N(u)≡{r['N_mod']}  "
                f"N≡0: {'JA' if r['N_resonanz'] else 'NEIN'}  "
                f"I_{p}: {'JA' if r['I_p'] else 'NEIN'}  ({r['grund']})"
            )

    any_norm = any(res[p]["N_resonanz"] for p in (5, 7))
    any_ideal = any(res[p]["I_p"] for p in (5, 7))
    return {
        "u": eabc_from_q(u),
        "N_u": nu,
        "N_q": nq,
        "N_qp": nqp,
        "resonanz_mod": any_norm,
        "resonanz_ideal": any_ideal,
        "detail": res,
    }


# H32 / Zeitwürfel (wie rundweg_hurwitz_resonanz.sage)
H32_POS = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59]
H32 = sorted([-r for r in H32_POS] + H32_POS)


def eabc_class_mod12(n):
    r = int(n) % 12
    if r == 1:
        return "E"
    if r == 5:
        return "A"
    if r == 7:
        return "B"
    if r == 11:
        return "C"
    return None


def quaternion_from_flavor_counts(counts):
    return (
        H(counts.get("E", 0))
        + H(counts.get("A", 0)) * _i
        + H(counts.get("B", 0)) * _j
        + H(counts.get("C", 0)) * _k
    )


def quaternion_from_label(label):
    vec = {
        "E": (1, 0, 0, 0),
        "A": (0, 1, 0, 0),
        "B": (0, 0, 1, 0),
        "C": (0, 0, 0, 1),
    }.get(label)
    if vec is None:
        return None
    return q_from_eabc(*vec)


def _read_csv_rows(path):
    if not os.path.isfile(path):
        return None
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _float_or_none(x):
    if x is None or x == "":
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if v != v:
        return None
    return v


def scan_schwung_resonanz(M, *, data_dir=None, events_csv=None, trans_csv=None):
    """
    Alle relevanten (q_k, q_{k+1})-Paare aus Zeitwürfel-CSV:
    u = conj(q_k)·q_{k+1}, Resonanz via N(u) mod 5/7 und I_5/I_7.
    """
    if data_dir is None:
        data_dir = os.getcwd()
    ev_path = events_csv or os.path.join(data_dir, "bamberger_zeitwuerfel_active_events.csv")
    tr_path = trans_csv or os.path.join(data_dir, "bamberger_zeitwuerfel_active_transitions.csv")

    print(f"\n=== scan_schwung_resonanz(M={M}) ===")
    print(f"  events: {ev_path}")
    print(f"  transitions: {tr_path}")

    transitions = _read_csv_rows(tr_path)
    events = _read_csv_rows(ev_path)

    pairs = []
    if transitions:
        seen = set()
        for row in transitions:
            comps = tuple(
                int(round(_float_or_none(row.get(c)) or 0))
                for c in ("u_e", "u_a", "u_b", "u_c")
            )
            if comps == (0, 0, 0, 0) or comps in seen:
                continue
            seen.add(comps)
            qu = q_from_eabc(*comps)
            res = schwung_resonanz_modulo(qu)
            pairs.append({"u": comps, "res": res, "N": q_norm(qu)})

    # Label-Folge aus events → Einheitsquaternion-Paare
    label_pairs = []
    if events:
        labels = [row.get("label", "") for row in events if row.get("label") in ("E", "A", "B", "C")]
        for k in range(len(labels) - 1):
            qk = quaternion_from_label(labels[k])
            qkp1 = quaternion_from_label(labels[k + 1])
            if qk is None or qkp1 is None:
                continue
            u = quadratischer_schwung(qk, qkp1)
            res = schwung_resonanz_modulo(u)
            label_pairs.append({
                "k": k,
                "labels": (labels[k], labels[k + 1]),
                "u": eabc_from_q(u),
                "N": q_norm(u),
                "res": res,
            })

    print(f"  Distinkte Schwung-Typen (CSV u_e..u_c): {len(pairs)}")
    print(f"  Label-Übergänge (E/A/B/C-Einheiten): {len(label_pairs)}")

    n_norm5 = sum(1 for p in pairs if p["res"][5]["N_resonanz"])
    n_norm7 = sum(1 for p in pairs if p["res"][7]["N_resonanz"])
    n_i5 = sum(1 for p in pairs if p["res"][5]["I_p"])
    n_i7 = sum(1 for p in pairs if p["res"][7]["I_p"])

    print(f"  CSV-Schwünge mit N(u)≡0 (mod 5): {n_norm5}/{len(pairs)}")
    print(f"  CSV-Schwünge mit N(u)≡0 (mod 7): {n_norm7}/{len(pairs)}")
    print(f"  CSV-Schwünge in I_5: {n_i5}/{len(pairs)}  in I_7: {n_i7}/{len(pairs)}")

    print("\n--- Erste distinkte CSV-Schwünge ---")
    for p in pairs[:12]:
        r5, r7 = p["res"][5], p["res"][7]
        print(
            f"  u={p['u']}  N={p['N']}  "
            f"mod5={r5['N_mod']}({'JA' if r5['N_resonanz'] else 'NEIN'})  "
            f"mod7={r7['N_mod']}({'JA' if r7['N_resonanz'] else 'NEIN'})  "
            f"I_5={'JA' if r5['I_p'] else 'NEIN'}  I_7={'JA' if r7['I_p'] else 'NEIN'}"
        )
    if len(pairs) > 12:
        print(f"  ... ({len(pairs) - 12} weitere)")

    ln_norm5 = sum(1 for p in label_pairs if p["res"][5]["N_resonanz"])
    ln_norm7 = sum(1 for p in label_pairs if p["res"][7]["N_resonanz"])
    ln_i5 = sum(1 for p in label_pairs if p["res"][5]["I_p"])
    ln_i7 = sum(1 for p in label_pairs if p["res"][7]["I_p"])
    print(f"\n  Label-Paare mit N(u)≡0 (mod 5): {ln_norm5}/{len(label_pairs)}")
    print(f"  Label-Paare mit N(u)≡0 (mod 7): {ln_norm7}/{len(label_pairs)}")
    print(f"  Label-Paare in I_5: {ln_i5}/{len(label_pairs)}  in I_7: {ln_i7}/{len(label_pairs)}")

    any_res = (
        n_norm5 + n_norm7 + n_i5 + n_i7 + ln_norm5 + ln_norm7 + ln_i5 + ln_i7
    ) > 0
    return {
        "M": M,
        "n_csv_pairs": len(pairs),
        "n_label_pairs": len(label_pairs),
        "csv_norm5": n_norm5,
        "csv_norm7": n_norm7,
        "csv_I5": n_i5,
        "csv_I7": n_i7,
        "label_norm5": ln_norm5,
        "label_norm7": ln_norm7,
        "label_I5": ln_i5,
        "label_I7": ln_i7,
        "any_resonanz": any_res,
        "pairs": pairs,
        "label_pairs": label_pairs,
    }


def test_zustaende():
    """Referenztests: q_lokal, q_asym mit sinnvollen q'-Nachbarn."""
    q_lokal = q_from_eabc(5, 2, 0, 1)
    q_asym = q_from_eabc(0, 1, 0, -1)
    q_glatt = q_from_eabc(*kollaps_eabc(5, 2, 0, 1))
    q_e = q_from_eabc(1, 0, 0, 0)

    print("Referenzzustände:")
    print(f"  q_lokal = {eabc_from_q(q_lokal)},  Π_Γ → {kollaps_eabc(5, 2, 0, 1)}")
    print(f"  q_asym  = {eabc_from_q(q_asym)}")
    print(f"  q_glatt = {eabc_from_q(q_glatt)} (kollabiertes q_lokal)")

    results = []
    tests = [
        (q_lokal, q_glatt, "q_lokal", "q_glatt"),
        (q_lokal, q_e, "q_lokal", "E"),
        (q_asym, q_glatt, "q_asym", "q_glatt"),
        (q_asym, q_e, "q_asym", "E"),
        (q_asym, q_lokal, "q_asym", "q_lokal"),
    ]
    for q, qp, lq, lqp in tests:
        results.append(schwung_report(q, qp, lq, lqp))

    any_test = any(r["resonanz_mod"] or r["resonanz_ideal"] for r in results)
    print(f"\n>>> Test-Paare: Resonanz mod 5/7 oder I_5∪I_7 = {'ja' if any_test else 'nein'}")
    return results, any_test


def main():
    print("=== QUADRATISCHER SCHWUNG — RESONANZ MOD 5/7 ===")
    print("u = q̄ q'  (Konjugation: q.conjugate()); N = reduced_norm")
    print()

    pqs = primquaternionen()
    print("Primquaternionen:")
    for p, pi in sorted(pqs.items()):
        print(f"  π_{p} = {eabc_from_q(pi)},  N = {q_norm(pi)}")
    print()

    _, any_test = test_zustaende()

    print("\n" + "=" * 60)
    print("M=113160 — Zeitwürfel-Schwung-Scan")
    print("=" * 60)
    scan = scan_schwung_resonanz(113160)
    any_scan = scan["any_resonanz"]
    print(
        f"\n>>> M=113160 Schwung-Scan: Resonanz N(u) mod 5/7 oder I_5∪I_7 = "
        f"{'ja' if any_scan else 'nein'}"
    )
    print(f"    CSV: mod5={scan['csv_norm5']}, mod7={scan['csv_norm7']}, "
          f"I_5={scan['csv_I5']}, I_7={scan['csv_I7']} "
          f"(von {scan['n_csv_pairs']} distinkten Typen)")
    print(f"    Label: mod5={scan['label_norm5']}, mod7={scan['label_norm7']}, "
          f"I_5={scan['label_I5']}, I_7={scan['label_I7']} "
          f"(von {scan['n_label_pairs']} Übergängen)")

    print("\n=== Ende quadratischer Schwung-Scan ===")
    print("LaTeX: kapitel_5_8_schalen_kollaps.tex; linearer Vergleich: ./run_hurwitz_resonanz.sh")


if __name__ == "__main__":
    main()
