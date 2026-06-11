#!/usr/bin/env sage
"""
#Energiedoku - Modul: Riemann-Nullstellen-Kopplung (Spektral-Test)
Direkter Abgleich der Eigenwerte des Hamilton-Operators mit zeros6.npy.
"""

import os
import numpy as np

try:
    from sage.all import *
except ImportError as e:
    raise SystemExit(
        'Dieses Skript benötigt SageMath. Bitte ausführen mit: sage "Alle Kugeln.py"'
    ) from e

def spektral_test_riemann_nullstellen():
    print("=== Starte Riemann-Resonanz-Test gegen zeros6.npy ===")
    
    # 1. Laden der empirischen Nullstellensonde (Projektpfad)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for dateipfad in (
        os.path.expanduser("~/Projects/zeros6.npy"),
        os.path.join(script_dir, "zeros6.npy"),
        "zeros6.npy",
    ):
        if os.path.exists(dateipfad):
            break
        
    try:
        riemann_zeros = np.load(dateipfad)
        print(f"Datei geladen: Shape = {riemann_zeros.shape}, Typ = {riemann_zeros.dtype}")
    except Exception as e:
        print(f"Fehler beim Laden der Nullstellendatei: {e}")
        return

    # Extrahieren der Imaginärteile (die klassischen Riemannschen Quanten-Frequenzen gamma)
    if np.iscomplexobj(riemann_zeros):
        gamma_empirisch = np.sort(np.abs(np.imag(riemann_zeros)))
    elif len(riemann_zeros.shape) > 1 and riemann_zeros.shape[1] == 2:
        gamma_empirisch = np.sort(np.abs(riemann_zeros[:, 1]))
    else:
        gamma_empirisch = np.sort(np.abs(riemann_zeros))
        
    print(f"Anzahl verf" + u"\u00fc" + f"gbarer Riemann-Frequenzen f" + u"\u00fc" + f"r den Test: {len(gamma_empirisch)}")
    print(f"Erste empirische Nullstellen-Frequenzen: {gamma_empirisch[:3]}")
    print("-" * 60)

    # 2. Berechnung des Lie-Hamilton-Spektrums (Aus Ihrem verifizierten Lauf)
    # Basis-Liftvektor v_L aus der E8-Doppelkugel-Torsion und der 11/30 Primvierling-Taktphase
    v_L = vector(RDF, [3/5, 1/5, 1/5, 3/5, 0, -1, 2, -1])
    theta = 11/30
    lmbda_optimal = 1.1  # Kalibrierter Wert aus der Parameterstudie
    
    generator_vec = theta * vector(RDF, list(v_L[:4]) + list(lmbda_optimal * v_L[4:]))
    
    # Aufbau der schiefsymmetrischen 8x8 Matrix M_Omega
    M_Omega = matrix(CDF, 8, 8)
    for i in range(4):
        M_Omega[i, 7-i] = generator_vec[i]
        M_Omega[7-i, i] = -generator_vec[i]
        M_Omega[i+4, 7-(i+4)] = generator_vec[i+4]
        M_Omega[7-(i+4), i+4] = -generator_vec[i+4]
        
    # Hermitescher Hamilton-Operator H = i * M_Omega nach dem Hilbert-P\u00f3lya-W\u00f6rterbuch
    H_operator = I * M_Omega
    eigenwerte = H_operator.eigenvalues()
    
    # Extraktion der positiven Energieniveaus (Resonanzfrequenzen des Modells)
    omega_modell = sorted([abs(alpha.real()) for alpha in eigenwerte if alpha.real() > 1e-5])
    print(f"Berechnete positive Modell-Frequenzen (\u03c9_k): {omega_modell}")
    print("-" * 60)

    # 4. Asymptotischer statistischer Spektral-Abgleich (Weyl-Mittelung der Gaps)
    print("=== Starte Asymptotischen Riemann-Resonanz-Test ===")
    gaps = np.diff(gamma_empirisch)
    mean_gap = float(np.mean(gaps))
    print(f"Mittlerer Abstand zwischen Nullstellen (mean_gap): {mean_gap:.6f}")
    for idx, omega in enumerate(omega_modell):
        skalen_faktor = mean_gap / omega
        print(f"  E_{idx} (\u03c9={omega:.6f}): skalen_faktor = {skalen_faktor:.6f}", end="")
        if abs(skalen_faktor - round(skalen_faktor)) < 0.1:
            knoten = int(round(skalen_faktor))
            print(f"  -> Resonanz-Knoten (nahe {knoten})")
        else:
            print()
    print("-" * 60)

    # 5. Logarithmischer Lift (lokale von-Mangoldt-Dichte)
    # Mikroskopische omega_k liegen in Modell-Einheiten; makroskopische Nullstellen-
    # abstände folgen dN/dT. Lift: omega_lifted = omega * (delta(T)/mean_gap), analog zu
    # skalen_faktor = mean_gap/omega in Abschnitt 4 — dort global, hier lokal bei Hoehe T.
    print("=== Logarithmischer Lift (asymptotische lokale Dichte) ===")

    Eulersche = exp(1)

    def N_von_mangoldt(T):
        T = float(T)
        return float((T / (2 * pi)) * log(T / (2 * pi * Eulersche)))

    def delta_asymptotisch(T):
        """Lokaler mittlerer Nullstellenabstand 1/(dN/dT) bei Hoehe T."""
        T = float(T)
        dN_dT = (1 / (2 * pi)) * (log(T / (2 * pi * Eulersche)) + 1)
        return float(1.0 / dN_dT)

    n = len(gamma_empirisch)
    T_median = float(np.median(gamma_empirisch))
    T_mitte = float(gamma_empirisch[n // 2])
    T_mittel = float(np.mean(gamma_empirisch))
    T_repr = T_mitte

    for label, T in (
        ("Median", T_median),
        ("Index n/2", T_mitte),
        ("Mittel", T_mittel),
    ):
        delta_T = delta_asymptotisch(T)
        delta_alt = float((2 * pi) / log(float(T) / (2 * pi)))
        lift_faktor = delta_T / mean_gap
        print(
            f"  T ({label}) = {T:.4f}: delta(T) = {delta_T:.6f}, "
            f"2*pi/log(T/2pi) = {delta_alt:.6f}, lift_faktor = {lift_faktor:.6f}"
        )

    delta_repr = delta_asymptotisch(T_repr)
    lift_faktor_repr = delta_repr / mean_gap
    N_repr = N_von_mangoldt(T_repr)
    print(f"Repraesentative Hoehe T = {T_repr:.4f} (Index n/2), N(T) ~ {N_repr:.2f}")
    print(f"Lift-Faktor f(T) = delta(T)/mean_gap = {lift_faktor_repr:.6f}")
    print("Geliftete Modell-Frequenzen (omega_lifted = omega * f(T)):")

    k_mitte = n // 2
    gap_lokal_1 = float(gamma_empirisch[k_mitte + 1] - gamma_empirisch[k_mitte])
    gap_lokal_2 = float(gamma_empirisch[k_mitte + 2] - gamma_empirisch[k_mitte + 1])
    gap_lokal_mittel = float(np.mean(np.diff(gamma_empirisch[k_mitte - 2 : k_mitte + 3])))

    toleranz_lift = 0.08
    for idx, omega in enumerate(omega_modell):
        omega_lifted = omega * lift_faktor_repr
        print(f"  E_{idx}: omega = {omega:.6f} -> omega_lifted = {omega_lifted:.6f}", end="")
        for name, gap_ref in (
            ("Delta_gamma lokal", gap_lokal_1),
            ("Delta_gamma+1", gap_lokal_2),
            ("mittl. 5 Nachbarn", gap_lokal_mittel),
            ("mean_gap global", mean_gap),
        ):
            verhaeltnis = gap_ref / omega_lifted
            if abs(verhaeltnis - round(verhaeltnis)) < toleranz_lift:
                print(
                    f"  | Resonanz mit {name} ({gap_ref:.4f}): "
                    f"~ {int(round(verhaeltnis))} (Verh. {verhaeltnis:.4f})",
                    end="",
                )
        skalen_lift = mean_gap / omega_lifted
        if abs(skalen_lift - round(skalen_lift)) < 0.1:
            print(f"  | skalen_faktor (global) ~ {int(round(skalen_lift))}", end="")
        print()
    print("-" * 60)

    # 3. Mathematischer Kreuzabgleich (Resonanz-Suche)
    # Wir testen, ob die Abstände zwischen den empirischen Riemann-Nullstellen (Gaps)
    # ganzzahlige oder rationale Resonanzen mit unseren Gruppen-Frequenzen bilden.
    print("Suche nach rationalen Symmetrie-Kopplungen (Kollaps-Kriterien):")
    
    for i in range(min(len(gamma_empirisch) - 1, 5)):
        gap_riemann = gamma_empirisch[i+1] - gamma_empirisch[i]
        
        for idx, omega in enumerate(omega_modell):
            verhaeltnis = gap_riemann / omega
            # Prüfe auf harmonische oder subharmonische Resonanz (nahe an einfachen Brüchen)
            if abs(verhaeltnis - round(verhaeltnis)) < 0.05:
                print(f"  -> MATCH: Riemann-Abstand \u0394_{i}({gap_riemann:.4f}) koppelt harmonisch mit Modell-Energie E_{idx}({omega:.4f})")
                print(f"     Verh" + u"\u00e4" + f"ltnis: {verhaeltnis:.4f} \u2248 {round(verhaeltnis)}")

    print("=== Spektral-Test erfolgreich abgeschlossen ===")

if __name__ == "__main__":
    spektral_test_riemann_nullstellen()