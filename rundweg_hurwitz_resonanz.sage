#!/usr/bin/env sage
# -*- coding: utf-8 -*-
# Pfad B: Hurwitz-Gitter — arithmetische Resonanz (Π_Γ, Ideale I_5, I_7)
# Start: ./run_hurwitz_resonanz.sh  oder  sage rundweg_hurwitz_resonanz.sage

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
    """
    Hurwitz-Maximalordnung: Koeffizienten ganz oder halbzahlig,
    Summe der Koeffizienten in ZZ.
    """
    c = eabc_from_q(q)
    all_int = all(x in ZZ for x in c)
    all_half = all(2 * x in ZZ for x in c) and not all_int
    if not (all_int or all_half):
        return False
    return sum(c) in ZZ


def init_kollaps_operators():
    R = Matrix(QQ, [[0, 0, 0, 1],
                    [1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0]])
    K = Matrix(QQ, [[0, 0, 0, 1],
                    [0, 0, 1, 0],
                    [0, 1, 0, 0],
                    [1, 0, 0, 0]])
    Gamma = R * K
    return R, K, Gamma


def kollaps_operator_q():
    """Pi_Gamma = (I + Gamma + Gamma^2 + Gamma^3) / 4."""
    _, _, Gamma = init_kollaps_operators()
    I = Matrix.identity(QQ, 4)
    return (I + Gamma + Gamma**2 + Gamma**3) / 4


def kollaps_eabc(e, a, b, c):
    """Π_Γ auf (e,a,b,c): Ausgabe (e, (a+c)/2, b, (a+c)/2)."""
    Pi = kollaps_operator_q()
    v = Matrix(QQ, [[e], [a], [b], [c]])
    w = Pi * v
    return tuple(w.list())


def primquaternionen():
    """
    π_5 = 2+i (Norm 5), π_7 = 2+i+j+k (Norm 7, ganzzahlig in H).
    (3/2)(1+i+j+k) hat Norm 9 — nur als Negativbeispiel.
    """
    pi5 = H(2) + _i
    pi7 = H(2) + _i + _j + _k
    pi7_falsch = (H(3) / 2) * (H(1) + _i + _j + _k)
    return {
        5: pi5,
        7: pi7,
        "pi7_halbzahlig_falsch": pi7_falsch,
    }


def norm_notwendig(q, pi):
    """Notwendig: N(pi) | N(q) für q in H·pi (Linksideal)."""
    nq = q_norm(q)
    npi = q_norm(pi)
    if nq == 0:
        return True
    return (nq % npi == 0) if npi != 0 else False


def in_linksideal(q, pi, *, verbose=False):
    """
    q ∈ H·π  (Linksideal)  ⇔  ∃ h ∈ H (Hurwitz): q = h·π.
    Operativ: h = q·π^{-1}, prüfe ob h Hurwitz ist (π ≠ 0).
    """
    if q.is_zero():
        return True, H(0), "trivial"
    if not norm_notwendig(q, pi):
        return False, None, "Norm-Kriterium verletzt"
    h = q * pi**(-1)
    ok = is_hurwitz(h)
    grund = "h=q*pi^-1 ist Hurwitz" if ok else "h=q*pi^-1 nicht Hurwitz"
    if verbose:
        print(f"    h = {eabc_from_q(h)}, Hurwitz={ok}")
    return ok, h, grund


def in_rechtsideal_pi_H(q, pi, *, verbose=False):
    """q ∈ π·H  ⇔  h = π^{-1}·q ist Hurwitz."""
    if q.is_zero():
        return True, H(0), "trivial"
    if not norm_notwendig(q, pi):
        return False, None, "Norm-Kriterium verletzt"
    h = pi**(-1) * q
    ok = is_hurwitz(h)
    grund = "h=pi^-1*q ist Hurwitz" if ok else "h=pi^-1*q nicht Hurwitz"
    if verbose:
        print(f"    h = {eabc_from_q(h)}, Hurwitz={ok}")
    return ok, h, grund


def resonanz_report(q, label, primes):
    print(f"\n--- Zustand {label}: q = {eabc_from_q(q)}, N(q) = {q_norm(q)} ---")
    for p in primes:
        pi = primquaternionen()[p]
        print(f"  Ideal I_{p}, π_p = {eabc_from_q(pi)}, N(π_p) = {q_norm(pi)}")
        ok_l, _, grund_l = in_linksideal(q, pi, verbose=True)
        ok_r, _, grund_r = in_rechtsideal_pi_H(q, pi, verbose=True)
        print(f"    Linksideal H·π_p (q=hπ):     {'JA' if ok_l else 'NEIN'}  ({grund_l})")
        print(f"    Rechtsideal π_p·H (q=πh'):   {'JA' if ok_r else 'NEIN'}  ({grund_r})")


# H32-Hülle (aus H32_eabc.py) für M-Scan
H32_POS = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59]
H32 = sorted([-r for r in H32_POS] + H32_POS)

FLAVOR_TO_IDX = {"E": 0, "A": 1, "B": 2, "C": 3}


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
    """Aggregat q = n_E·1 + n_A·i + n_B·j + n_C·k aus EABC-Zählern."""
    return (
        H(counts.get("E", 0))
        + H(counts.get("A", 0)) * _i
        + H(counts.get("B", 0)) * _j
        + H(counts.get("C", 0)) * _k
    )


def quaternion_from_label(label):
    """E/A/B/C → Einheits-Quaternion (wie Jitter Zeit.py / H32_eabc)."""
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
    if v != v:  # NaN
        return None
    return v


def _counts_from_labels(labels):
    counts = {"E": 0, "A": 0, "B": 0, "C": 0}
    for lab in labels:
        if lab in counts:
            counts[lab] += 1
    return counts


def zeitwuerfel_resonanz_fuer_M(M, *, data_dir=None, events_csv=None, trans_csv=None):
    """
    Zeitwürfel-Zustände (Jitter Zeit.py-CSV) an Hurwitz-Ideale I_5, I_7.

    Kodierung (identisch zu H32_eabc / Sage):
      EABC mod 12; E/A/B/C → 1,i,j,k; Schwung u = conj(q_k)·q_{k+1}.

    Phasenkopplung an Phasenzentrum M (M≡0 mod 30):
      n_k aus CSV liegt typisch ≪ M → Kopplung über r̂ = n_k mod 30
      gegen belegte H32-Slots r mit is_prime(M+r).
    """
    if data_dir is None:
        data_dir = os.getcwd()
    ev_path = events_csv or os.path.join(data_dir, "bamberger_zeitwuerfel_active_events.csv")
    tr_path = trans_csv or os.path.join(data_dir, "bamberger_zeitwuerfel_active_transitions.csv")

    print(f"\n=== zeitwuerfel_resonanz_fuer_M(M={M}) ===")
    print(f"  events: {ev_path}")
    print(f"  transitions: {tr_path}")

    occupied = [r for r in H32 if is_prime(M + r)]
    mod30_slots = {}
    for r in occupied:
        key = int(r) % 30
        mod30_slots.setdefault(key, []).append(r)

    events = _read_csv_rows(ev_path)
    transitions = _read_csv_rows(tr_path)
    if events is None:
        print("  Keine Zeitwürfel-CSV — bitte zuerst: python3 'Jitter Zeit.py'")
        return {"loaded": False}

    active_labels = []
    phase_labels = []
    phase_match = 0
    phase_total = 0

    for row in events:
        lab = row.get("label", "")
        if lab not in ("E", "A", "B", "C"):
            continue
        active_labels.append(lab)
        nk = row.get("n_k")
        if nk is None:
            continue
        nk = int(float(nk))
        rhat = nk % 30
        if rhat not in mod30_slots:
            continue
        phase_total += 1
        phase_labels.append(lab)
        for r in mod30_slots[rhat]:
            if eabc_class_mod12(M + r) == lab:
                phase_match += 1
                break

    counts_all = _counts_from_labels(active_labels)
    counts_phase = _counts_from_labels(phase_labels)
    q_all = quaternion_from_flavor_counts(counts_all)
    q_phase = quaternion_from_flavor_counts(counts_phase)

    print(f"  Aktive EABC-Ereignisse: {len(active_labels)}")
    print(f"  Mod-30-gekoppelt an H32-Slots von M: {phase_total} (Flavor-Treffer {phase_match}/{phase_total})")
    print(f"  Zähler gesamt: {counts_all}  →  q = {eabc_from_q(q_all)}, N = {q_norm(q_all)}")
    print(f"  Zähler phasengekoppelt: {counts_phase}  →  q = {eabc_from_q(q_phase)}, N = {q_norm(q_phase)}")

    primes = [5, 7]
    resonanz_report(q_all, f"Zeitwürfel-Aggregat (aktiv, M={M})", primes)
    if phase_total > 0:
        resonanz_report(q_phase, f"Zeitwürfel phasengekoppelt (n_k mod 30 ∈ H32|M)", primes)

    schwung_hits = []
    if transitions:
        seen_u = set()
        print("\n--- Zeitwürfel-Schwung u (aktive Übergänge, ganzzahlig) ---")
        for row in transitions:
            comps = tuple(
                int(round(_float_or_none(row.get(c)) or 0))
                for c in ("u_e", "u_a", "u_b", "u_c")
            )
            if comps == (0, 0, 0, 0) or comps in seen_u:
                continue
            seen_u.add(comps)
            qu = q_from_eabc(*comps)
            ok5 = in_linksideal(qu, primquaternionen()[5])[0]
            ok7 = in_linksideal(qu, primquaternionen()[7])[0]
            schwung_hits.append({"u": comps, "I5": ok5, "I7": ok7})
            print(f"  u={comps}  N={q_norm(qu)}  I_5={'JA' if ok5 else 'NEIN'}  I_7={'JA' if ok7 else 'NEIN'}")
        n5 = sum(1 for h in schwung_hits if h["I5"])
        n7 = sum(1 for h in schwung_hits if h["I7"])
        print(f"  Distinkte Schwünge: {len(schwung_hits)}  mit I_5: {n5}  mit I_7: {n7}")

    return {
        "loaded": True,
        "M": M,
        "n_active": len(active_labels),
        "n_phase": phase_total,
        "phase_flavor_match": phase_match,
        "counts_all": counts_all,
        "counts_phase": counts_phase,
        "q_all": eabc_from_q(q_all),
        "q_phase": eabc_from_q(q_phase) if phase_total else None,
        "I5_agg": in_linksideal(q_all, primquaternionen()[5])[0],
        "I7_agg": in_linksideal(q_all, primquaternionen()[7])[0],
        "I5_phase": in_linksideal(q_phase, primquaternionen()[5])[0] if phase_total else None,
        "I7_phase": in_linksideal(q_phase, primquaternionen()[7])[0] if phase_total else None,
        "schwung_hits": schwung_hits,
    }


def scan_resonanz(M, *, max_unit_norm=12, tabulate_units=False):
    """
    M=113160: H32-Phasenzentrum (H32_eabc.test_H32), zugleich Ptolemäus-Zentrum
    10^7 + 113160 im Skript Ptolo Norm.py — Bamberger-Zeitwürfel (Jitter Zeit.py)
    ist ein separates Modul ohne direkte Hurwitz-Ideal-Kopplung im Repo.

    Heuristik: belegte H32-Slots → EABC-Zähler → Aggregat-Quaternion;
    optional kleine Hurwitz-Einheiten mit N ≤ max_unit_norm.
    """
    print(f"\n=== scan_resonanz(M={M}) ===")
    print(f"Faktorisierung M: {factor(M)}")
    if M % 30 != 0:
        print("Hinweis: M ist nicht ≡ 0 (mod 30); H32-Interpretation ist nur heuristisch.")

    occupied = [r for r in H32 if is_prime(M + r)]
    missing = [r for r in H32 if not is_prime(M + r)]
    print(f"H32: {len(occupied)} belegte / {len(missing)} freie Slots (von {len(H32)})")

    counts = {"E": 0, "A": 0, "B": 0, "C": 0}
    slot_states = []
    for r in occupied:
        fl = eabc_class_mod12(M + r)
        if fl is None:
            continue
        counts[fl] += 1
        slot_states.append((r, M + r, fl))

    print(f"EABC-Zähler auf belegten Slots: {counts}")
    q_agg = quaternion_from_flavor_counts(counts)
    print(f"Aggregat-Quaternion (Slot-Summe): {eabc_from_q(q_agg)}, N = {q_norm(q_agg)}")

    primes = [5, 7]
    resonanz_report(q_agg, f"H32-Aggregat M={M}", primes)

    # Einzel-Slots als Einheitsrichtungen E/A/B/C
    print("\n--- Einzel-Slots (Einheitsrichtung je Flavor) ---")
    rows = []
    for r, p, fl in slot_states[:16]:
        vec = {"E": (1, 0, 0, 0), "A": (0, 1, 0, 0), "B": (0, 0, 1, 0), "C": (0, 0, 0, 1)}[fl]
        qs = q_from_eabc(*vec)
        row = {"r": r, "p": p, "fl": fl, "I5": in_linksideal(qs, primquaternionen()[5])[0],
               "I7": in_linksideal(qs, primquaternionen()[7])[0]}
        rows.append(row)
        print(f"  r={r:+3d}  p={p}  {fl}:  I_5={'JA' if row['I5'] else 'NEIN'}  I_7={'JA' if row['I7'] else 'NEIN'}")
    if len(slot_states) > 16:
        print(f"  ... ({len(slot_states) - 16} weitere Slots)")

    if tabulate_units:
        print("\n--- Kleine Hurwitz-Einheiten (N ≤ %s) — Resonanz mit I_5, I_7 ---" % max_unit_norm)
        units = enumerate_hurwitz_units(max_norm=max_unit_norm)
        hit5 = [u for u in units if in_linksideal(u, primquaternionen()[5])[0]]
        hit7 = [u for u in units if in_linksideal(u, primquaternionen()[7])[0]]
        print(f"  Gefunden: {len(units)} Hurwitz-Elemente mit 0 < N ≤ {max_unit_norm}")
        print(f"  In I_5: {len(hit5)}   In I_7: {len(hit7)}")
        for u in hit5[:5]:
            print(f"    I_5: {eabc_from_q(u)}")
        for u in hit7[:5]:
            print(f"    I_7: {eabc_from_q(u)}")

    return {
        "M": M,
        "occupied": len(occupied),
        "counts": counts,
        "q_agg": eabc_from_q(q_agg),
        "slot_rows": rows,
    }


def enumerate_hurwitz_units(max_norm=12):
    """Brute-Force-Hurwitz-Elemente mit 0 < N(q) ≤ max_norm (für Heuristik)."""
    found = []
    bound = max_norm + 2
    half = [ZZ(n) / 2 for n in range(-bound, bound + 1)]
    ints = list(range(-bound, bound + 1))
    for mode in ("int", "half"):
        coeffs_list = ints if mode == "int" else half
        for a in coeffs_list:
            for b in coeffs_list:
                for c in coeffs_list:
                    for d in coeffs_list:
                        q = H(a) + H(b) * _i + H(c) * _j + H(d) * _k
                        if not is_hurwitz(q) or q.is_zero():
                            continue
                        n = q_norm(q)
                        if 0 < n <= max_norm and n == 1:
                            found.append(q)
    # Deduplizieren
    seen = set()
    out = []
    for q in found:
        key = eabc_from_q(q)
        if key not in seen:
            seen.add(key)
            out.append(q)
    return out


def main():
    print("=== ARITHMETISCHE RESONANZ-ANALYSE ÜBER HURWITZ-GITTER ===")
    print("Pfad B (#Energiedoku): E,A,B,C ↔ 1,i,j,k; Π_Γ aus Rundweg.py")
    print()

    Pi = kollaps_operator_q()
    print("Pi_Gamma (Kollaps-Operator über Q^4):")
    print(Pi)
    print()

    pqs = primquaternionen()
    print("Primquaternionen (Normverifikation in Sage):")
    for key, pi in pqs.items():
        if isinstance(key, str):
            print(f"  {key}: {eabc_from_q(pi)}  N = {q_norm(pi)}  (nicht für I_7 verwenden)")
        else:
            print(f"  π_{key}: {eabc_from_q(pi)}  N = {q_norm(pi)}  Hurwitz={is_hurwitz(pi)}")
    print()

    q_lokal = q_from_eabc(5, 2, 0, 1)
    q_asym = q_from_eabc(0, 1, 0, -1)
    koll_lokal = kollaps_eabc(5, 2, 0, 1)
    koll_asym = kollaps_eabc(0, 1, 0, -1)

    print("Testzustand q_lokal = (5, 2, 0, 1):")
    print(f"  Π_Γ-Kollaps: {koll_lokal}")
    print("Asymmetrischer Eigenvektor V_{λ=0}: q_asym = (0, 1, 0, -1) = i - k")
    print(f"  Π_Γ(q_asym): {koll_asym}  (erwartet 0 im λ=0-Raum)")

    resonanz_report(q_lokal, "q_lokal (5,2,0,1)", [5, 7])
    resonanz_report(q_asym, "q_asym (antisymmetrisch)", [5, 7])

    # Schnitt V_{λ=0} mit Linksideal: Skalare Vielfache von q_asym
    print("\n--- Schnitt V_{λ=0} ∩ H·π_p (nur Skalare in Q) ---")
    for p in [5, 7]:
        pi = pqs[p]
        any_hit = False
        for t in [1, 2, -1, QQ(1) / 2]:
            qt = H(t) * q_asym
            ok, _, _ = in_linksideal(qt, pi)
            if ok and not qt.is_zero():
                any_hit = True
                print(f"  p={p}, t={t}: JA")
        if not any_hit:
            print(f"  p={p}: kein nichttriviales rationales Vielfach in H·π_p")

    print("\n" + "=" * 60)
    print("M=113160 — Bamberger/H32-Kontext")
    print("  • H32_eabc.py: test_H32(113160) — Mod-30-Phasenzentrum, 32-Slot-Hülle")
    print("  • Ptolo Norm.py: Zentrum 10^7 + 113160 (lokale Primzahl-Statistik)")
    print("  • Jitter Zeit.py: Bamberger-Zeitwürfel → zeitwuerfel_resonanz_fuer_M (CSV/Mod-30)")
    print("=" * 60)

    scan_resonanz(113160, tabulate_units=False)
    zw = zeitwuerfel_resonanz_fuer_M(113160)
    if zw.get("loaded"):
        any_res = zw["I5_agg"] or zw["I7_agg"] or zw.get("I5_phase") or zw.get("I7_phase")
        any_res = any_res or any(h["I5"] or h["I7"] for h in zw.get("schwung_hits", []))
        print(f"\n>>> Zeitwürfel↔Hurwitz M=113160: geladen=ja, arithm. Resonanz I_5∪I_7={'ja' if any_res else 'nein'}")

    print("\n=== Ende Hurwitz-Resonanzanalyse ===")
    print("Empfehlung systematischer Einheiten-Scan: siehe PAPER_HURWITZ_RESONANZ.md")


if __name__ == "__main__" and not globals().get("_skip_hurwitz_main"):
    main()
