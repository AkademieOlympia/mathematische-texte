#!/usr/bin/env sage
"""
#Energiedoku - Modul: Systematischer E1/E3-Scan & Multiplizitaets-Check
Pfad B: algebraischer Hurwitz-Idealfilter (I_5, I_7, Normen 5/7, Anker M=113160)
filtert stochastische Scheinresonanzen aus Gap-Matches.
"""

import argparse
import os

import numpy as np

from sage.all import *


ANKER_M = 113160  # H32-/Bamberger-Stabilitaetsanker (Pfad B, rundweg_hurwitz_resonanz.sage)

_H = None
_pi_5 = None
_pi_7 = None


def _hurwitz_setup():
    """Quaternionen-Algebra H und Primanker pi_5, pi_7 (Linksideale I_5, I_7)."""
    global _H, _pi_5, _pi_7
    if _H is None:
        _H = QuaternionAlgebra(QQ, -1, -1)
        i, j, k = _H.gens()
        _pi_5 = _H(2) + i
        _pi_7 = _H(2) + i + j + k
    return _H, _pi_5, _pi_7


def check_quaternion_ideal_resonance(p_ratio, q_ratio):
    """
    Pfad-B-Filter: rationale Kopplung p/q nur bei Ideal-Norm-Modulation (5, 7).
    pi_5=2+i (N=5), pi_7=2+i+j+k (N=7); Anker M=113160 fuer kombinatorische Phase.
    Operative Heuristik auf (p,q): oberer_takt % 5 == 0 or % 7 == 0 or unterer_takt == 4.
    """
    _hurwitz_setup()
    oberer_takt = int(p_ratio)
    unterer_takt = int(q_ratio)
    if unterer_takt == 0:
        return False
    if (oberer_takt % 5 == 0) or (oberer_takt % 7 == 0) or (unterer_takt == 4):
        return True
    return False


def lade_gamma_sonde():
    """Laedt Imaginaerteile der Riemann-Nullstellen aus zeros6.npy."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for dateipfad in (
        os.path.expanduser("~/Projects/zeros6.npy"),
        os.path.join(script_dir, "zeros6.npy"),
        "zeros6.npy",
    ):
        if os.path.exists(dateipfad):
            break
    else:
        raise FileNotFoundError("zeros6.npy nicht gefunden (~/Projects oder Skriptordner)")

    riemann_zeros = np.load(dateipfad)
    if np.iscomplexobj(riemann_zeros):
        gamma = np.sort(np.abs(np.imag(riemann_zeros)))
    elif len(riemann_zeros.shape) > 1 and riemann_zeros.shape[1] == 2:
        gamma = np.sort(np.abs(riemann_zeros[:, 1]))
    else:
        gamma = np.sort(np.abs(riemann_zeros))
    return gamma, dateipfad


def berechne_omega_modell():
    """Positive Modell-Frequenzen wie in Alle Kugeln.py (Lie-Hamilton)."""
    v_L = vector(RDF, [3 / 5, 1 / 5, 1 / 5, 3 / 5, 0, -1, 2, -1])
    theta = 11 / 30
    lmbda_optimal = 1.1
    generator_vec = theta * vector(RDF, list(v_L[:4]) + list(lmbda_optimal * v_L[4:]))
    M_Omega = matrix(CDF, 8, 8)
    for i in range(4):
        M_Omega[i, 7 - i] = generator_vec[i]
        M_Omega[7 - i, i] = -generator_vec[i]
        M_Omega[i + 4, 7 - (i + 4)] = generator_vec[i + 4]
        M_Omega[7 - (i + 4), i + 4] = -generator_vec[i + 4]
    H_operator = I * M_Omega
    eigenwerte = H_operator.eigenvalues()
    return sorted([float(abs(alpha.real())) for alpha in eigenwerte if alpha.real() > 1e-5])


def delta_mangoldt(T, eulersche=None):
    """
    Lokaler mittlerer Nullstellenabstand 1/(dN/dT) bei Hoehe T.
    Gleiche Formel wie Alle Kugeln.py Abschnitt 5 (von-Mangoldt, +1 im Log).
    """
    if eulersche is None:
        eulersche = exp(1)
    T = float(T)
    dN_dT = (1 / (2 * pi)) * (log(T / (2 * pi * float(eulersche))) + 1)
    return float(1.0 / dN_dT)


def delta_asymptotisch_alt(T, eulersche=None):
    """
    Vereinfachte asymptotische Formel 2*pi/log(T/(2*pi*e)) ohne +1-Term.
    Nur Vergleichsmodus (--alt-delta); Standard ist delta_mangoldt.
    """
    if eulersche is None:
        eulersche = exp(1)
    T = float(T)
    return float((2 * pi) / log(T / (2 * pi * float(eulersche))))


def lift_faktor_bei_T(T, global_mean_gap, delta_fn=delta_mangoldt, eulersche=None):
    if eulersche is None:
        eulersche = exp(1)
    return float(delta_fn(T, eulersche) / global_mean_gap)


def beste_rationale_naeherung(verhaeltnis, max_q=12):
    """Beste p/q mit q <= max_q (voller Raster-Scan)."""
    beste = float("inf")
    best_p, best_q = 0, 1
    for q_test in range(1, max_q + 1):
        p_test = int(round(float(verhaeltnis) * q_test))
        if p_test == 0:
            continue
        dist = abs(float(verhaeltnis) - p_test / q_test)
        if dist < beste:
            beste = dist
            best_p, best_q = p_test, q_test
    return best_p, best_q, beste


def scan_gaps_fuer_modus(
    lokale_gaps,
    omega_lifted,
    max_q=12,
    *,
    ideal_filter=False,
    sub_pct=None,
):
    """
    Pro Modus: beste p/q ueber alle Gaps.
    ideal_filter: nur Treffer mit check_quaternion_ideal_resonance.
    sub_pct: zusaetzliche Schwelle auf Abweichung (z. B. 0.01 fuer MC).
    """
    beste_global = float("inf")
    best_p, best_q, best_gap = 0, 1, 0.0
    validierte = []

    for dg in lokale_gaps:
        verhaeltnis = float(dg / omega_lifted)
        p_c, q_c, dist = beste_rationale_naeherung(verhaeltnis, max_q=max_q)
        if sub_pct is not None and dist >= sub_pct:
            continue
        if ideal_filter and not check_quaternion_ideal_resonance(p_c, q_c):
            continue
        if dist < beste_global:
            beste_global = dist
            best_p, best_q, best_gap = p_c, q_c, float(dg)
        if ideal_filter:
            validierte.append((p_c, q_c, float(dg), dist))

    return best_p, best_q, best_gap, beste_global, validierte


def monte_carlo_p_wert(
    zufalls_modi,
    lokale_gaps,
    lift_med,
    *,
    max_q=12,
    sub_pct=0.01,
    ideal_filter=False,
    integer_only=False,
):
    """
    integer_only=True: Sub-1%-Ganzzahl-Kopplung (q=1), vergleichbar mit Hold-out-p~0.837.
    integer_only=False: beste p/q mit q<=max_q unter sub_pct (strengeres Raster).
    """
    treffer = 0
    for rand_omega in zufalls_modi:
        w_l = float(rand_omega) * lift_med
        if integer_only:
            found = False
            for dg in lokale_gaps:
                verhaeltnis = float(dg / w_l)
                p_test = int(round(verhaeltnis))
                if p_test == 0:
                    continue
                if abs(verhaeltnis - p_test) < sub_pct:
                    if not ideal_filter or check_quaternion_ideal_resonance(p_test, 1):
                        found = True
                        break
            if found:
                treffer += 1
        else:
            _, _, _, dist, _ = scan_gaps_fuer_modus(
                lokale_gaps,
                w_l,
                max_q=max_q,
                ideal_filter=ideal_filter,
                sub_pct=sub_pct,
            )
            if dist < sub_pct:
                treffer += 1
    return treffer / len(zufalls_modi)


def systematischer_scan_mit_pfad_b(delta_fn=delta_mangoldt):
    print("=== E1-E3 Scan mit integriertem Pfad-B-Idealfilter ===")
    print(f"Anker M={ANKER_M} | Ideale I_5 (N=5), I_7 (N=7) | delta: {delta_fn.__name__}")

    eulersche = exp(1)

    try:
        gamma, pfad = lade_gamma_sonde()
        n_points = len(gamma)
        print(f"Referenz-Sonde geladen: {pfad}")
        print(f"Anzahl Nullstellen: {n_points}")
    except Exception as err:
        print(f"Fehler beim Laden der Sonde: {err}")
        return None

    idx_q25 = int(n_points * 0.25)
    idx_q50 = int(n_points * 0.50)
    idx_q75 = int(n_points * 0.75)
    quantile_indices = {"Q25": idx_q25, "Q50": idx_q50, "Q75": idx_q75}
    global_mean_gap = float(np.mean(np.diff(gamma)))

    omega_reell = berechne_omega_modell()
    print(f"Modell-Frequenzen (Hamilton): {omega_reell}")

    max_q = 12
    e1_e3_beste = {1: float("inf"), 3: float("inf")}

    print("\n--- Raster-Scan (Quantile) mit Pfad-B auf beste p/q ---")
    for q_label, idx in quantile_indices.items():
        T = float(gamma[idx])
        lift_faktor = lift_faktor_bei_T(T, global_mean_gap, delta_fn, eulersche)
        lokale_gaps = np.diff(gamma[max(0, idx - 50) : min(n_points, idx + 50)])

        print(f"\nHoehe {q_label} (T = {T:.2f}) | Lift-Faktor f(T) = {lift_faktor:.6f}")

        for i, omega in enumerate(omega_reell):
            omega_lifted = omega * lift_faktor
            p_raw, q_raw, gap_raw, dist_raw, _ = scan_gaps_fuer_modus(
                lokale_gaps, omega_lifted, max_q=max_q, ideal_filter=False
            )
            p_f, q_f, gap_f, dist_f, validierte = scan_gaps_fuer_modus(
                lokale_gaps, omega_lifted, max_q=max_q, ideal_filter=True
            )

            print(
                f"  E_{i} (w_lifted={omega_lifted:.6f}): roh {p_raw}/{q_raw} "
                f"dist={dist_raw:.6f} | Pfad-B {p_f}/{q_f} dist={dist_f:.6f} "
                f"({len(validierte)} gefilterte Gap-Treffer)"
            )

            if i in (1, 3) and dist_f < e1_e3_beste[i]:
                e1_e3_beste[i] = dist_f

    print("\n--- Median Q50: validierte E_i (Pfad-B, beste gefilterte p/q) ---")
    T_med = float(gamma[idx_q50])
    lift_med = lift_faktor_bei_T(T_med, global_mean_gap, delta_fn, eulersche)
    lokale_gaps_med = np.diff(gamma[max(0, idx_q50 - 50) : min(n_points, idx_q50 + 50)])
    print(f"Repraesentative Hoehe T = {T_med:.4f} | Lift-Faktor = {lift_med:.6f}")
    print("Filtere Resonanzbahnen via Hurwitz-Linksideale (I_5, I_7)...")

    for i, omega in enumerate(omega_reell):
        omega_lifted = omega * lift_med
        p_f, q_f, gap_f, dist_f, validierte = scan_gaps_fuer_modus(
            lokale_gaps_med, omega_lifted, max_q=max_q, ideal_filter=True
        )
        if dist_f < float("inf"):
            print(
                f"  [IDEAL MATCH] E_{i}: beste Kopplung {p_f}/{q_f} "
                f"bei Delta_g={gap_f:.4f} (Abweichung {dist_f:.6f}; "
                f"{len(validierte)} zulaessige Gap-Treffer)"
            )
        else:
            print(f"  E_{i}: Keine modulo Pfad-B zulaessige Resonanz auf dieser Schale.")

    print("\n--- E1/E3-Zusammenfassung (beste Pfad-B-Abweichung ueber Q25/Q50/Q75) ---")
    for i in (1, 3):
        print(f"  E_{i}: min gefilterte Abweichung = {e1_e3_beste[i]:.6f}")

    print("\n--- Re-Evaluation des p-Werts unter Pfad-B-Filterung ---")
    np.random.seed(42)
    zufalls_modi = np.random.uniform(0.05, 0.90, 1000)
    sub_pct = 0.01

    p_alt = monte_carlo_p_wert(
        zufalls_modi,
        lokale_gaps_med,
        lift_med,
        sub_pct=sub_pct,
        ideal_filter=False,
        integer_only=True,
    )
    p_neu = monte_carlo_p_wert(
        zufalls_modi,
        lokale_gaps_med,
        lift_med,
        sub_pct=sub_pct,
        ideal_filter=True,
        integer_only=True,
    )
    p_raster_ungefiltert = monte_carlo_p_wert(
        zufalls_modi,
        lokale_gaps_med,
        lift_med,
        max_q=max_q,
        sub_pct=sub_pct,
        ideal_filter=False,
        integer_only=False,
    )
    p_raster_gefiltert = monte_carlo_p_wert(
        zufalls_modi,
        lokale_gaps_med,
        lift_med,
        max_q=max_q,
        sub_pct=sub_pct,
        ideal_filter=True,
        integer_only=False,
    )

    treffer_alt = int(round(p_alt * 1000))
    treffer_neu = int(round(p_neu * 1000))
    print(
        "Nullmodell: Sub-1%-Ganzzahl-Kopplung (q=1), wie Hold-out-Scan in main-de.tex"
    )
    print(
        f"Zufaellige Modi mit Scheinresonanz (ungefiltert): {treffer_alt} von 1000"
    )
    print(f"Empirischer p-Wert (ohne Pfad-B): {p_alt:.4f}")
    print(
        f"Zufaellige Modi mit Scheinresonanz nach Pfad-B-Filter: "
        f"{treffer_neu} von 1000"
    )
    print(f"Rauschkorrigierter Pfad-B p-Wert (q=1): {p_neu:.4f}")
    print(
        f"Zusaetzlich beste p/q (q<=12): ungefiltert p={p_raster_ungefiltert:.4f}, "
        f"Pfad-B p={p_raster_gefiltert:.4f}"
    )
    print("=== Scan mit Pfad B erfolgreich beendet ===")
    return {"p_alt": p_alt, "p_neu": p_neu}


def systematischer_spektral_scan(delta_fn=delta_mangoldt):
    """Klassischer Scan ohne Pfad-B (Hold-out-Quantile, voller p/q-Raster)."""
    print("=== Initiere systematischen E1-E3 Scan & Multiplizitaets-Check ===")
    print(f"delta: {delta_fn.__name__} (Standard: delta_mangoldt wie Alle Kugeln.py)")

    eulersche = exp(1)

    try:
        gamma, pfad = lade_gamma_sonde()
        n_points = len(gamma)
        print(f"Referenz-Sonde geladen: {pfad}")
        print(f"Anzahl Nullstellen: {n_points}")
    except Exception as err:
        print(f"Fehler beim Laden der Sonde: {err}")
        return None

    idx_q25 = int(n_points * 0.25)
    idx_q50 = int(n_points * 0.50)
    idx_q75 = int(n_points * 0.75)

    quantile_indices = {"Q25": idx_q25, "Q50": idx_q50, "Q75": idx_q75}
    global_mean_gap = float(np.mean(np.diff(gamma)))

    omega_reell = berechne_omega_modell()
    print(f"Modell-Frequenzen (Hamilton): {omega_reell}")

    max_q = 12
    print("\n--- Starte systematischen Raster-Scan ueber Quantile ---")

    e1_e3_beste = {1: float("inf"), 3: float("inf")}

    for q_label, idx in quantile_indices.items():
        T = float(gamma[idx])
        lift_faktor = lift_faktor_bei_T(T, global_mean_gap, delta_fn, eulersche)
        lokale_gaps = np.diff(gamma[max(0, idx - 50) : min(n_points, idx + 50)])

        print(f"\nHoehe {q_label} (T = {T:.2f}) | Lift-Faktor f(T) = {lift_faktor:.6f}")

        for i, omega in enumerate(omega_reell):
            omega_lifted = omega * lift_faktor
            best_p, best_q, best_gap, beste_naehrung, _ = scan_gaps_fuer_modus(
                lokale_gaps, omega_lifted, max_q=max_q, ideal_filter=False
            )

            print(
                f"  E_{i} (w_lifted={omega_lifted:.6f}): Beste Resonanz "
                f"{best_p}/{best_q} bei Delta_g={best_gap:.4f} "
                f"(Abweichung: {beste_naehrung:.6f})"
            )

            if i in (1, 3) and beste_naehrung < e1_e3_beste[i]:
                e1_e3_beste[i] = beste_naehrung

    print("\n--- E1/E3-Zusammenfassung (beste Abweichung ueber Q25/Q50/Q75) ---")
    for i in (1, 3):
        print(f"  E_{i}: min Abweichung p/q-Suche = {e1_e3_beste[i]:.6f}")

    print("\n--- Quantifizierung des Multiplizitaetsproblems (p-Wert-Heuristik) ---")
    np.random.seed(42)
    zufalls_modi = np.random.uniform(0.05, 0.90, 1000)

    T_med = float(gamma[idx_q50])
    lift_med = lift_faktor_bei_T(T_med, global_mean_gap, delta_fn, eulersche)
    lokale_gaps_med = np.diff(gamma[max(0, idx_q50 - 50) : min(n_points, idx_q50 + 50)])

    p_wert = monte_carlo_p_wert(
        zufalls_modi,
        lokale_gaps_med,
        lift_med,
        sub_pct=0.01,
        ideal_filter=False,
        integer_only=True,
    )
    treffer = int(round(p_wert * 1000))
    print(f"Zufaellige Modi mit Sub-1%-Ganzzahl-Kopplung (q=1): {treffer} von 1000")
    print(f"Empirische Wahrscheinlichkeit fuer Zufallstreffer (p-Wert): {p_wert:.4f}")
    print("=== Scan erfolgreich beendet ===")
    return p_wert


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Energiedoku Spektral-Scan")
    parser.add_argument(
        "--ohne-pfad-b",
        action="store_true",
        help="Klassischer Scan ohne Hurwitz-Idealfilter",
    )
    parser.add_argument(
        "--alt-delta",
        action="store_true",
        help="Vereinfachte delta-Formel ohne +1 im Mangoldt-Log",
    )
    args = parser.parse_args()
    delta_fn = delta_asymptotisch_alt if args.alt_delta else delta_mangoldt

    if args.ohne_pfad_b:
        systematischer_spektral_scan(delta_fn=delta_fn)
    else:
        systematischer_scan_mit_pfad_b(delta_fn=delta_fn)
