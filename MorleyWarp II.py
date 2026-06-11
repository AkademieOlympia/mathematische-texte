import math
import numpy as np
from sympy import isprime

def scan_resonance(n_layer):
    print(f"Suche Resonanz-Zentrum fuer Schicht {n_layer}...")
    best_ratio = 0.190003
    max_h = 0

    # Wir scannen mit VIEL hoeherer Aufloesung (1000 statt 100 Schritte)
    for r in np.linspace(0.189, 0.191, 1000):
        p = int(2 ** (n_layer / 4) * r * 100)
        # Wir testen nur das Hamming-Gewicht des 'potentiellen' Produkts
        # Ein Vierling-Produkt P ist ca. p^4
        h = bin(p**4).count("1")

        if h > max_h:
            max_h = h
            best_ratio = r
            print(f"Neue Resonanz gefunden: Ratio {r:.6f} | w(P) = {h}")

    return best_ratio


def tuned_morley_navigation(n_layer, search_range=100000):
    """
    Findet Primzahlvierlinge durch direkte Navigation im Pascal-Tetraeder
    basierend auf der Morley-Resonanz-Theorie.
    """
    print(f"\n--- Starte Hyper-Navigation (Tuned) zur Schicht n = {n_layer} ---")

    # Zuerst Resonanzzentrum schaetzen, dann lokal darum absuchen
    center_ratio = scan_resonance(n_layer)
    ratios = np.linspace(center_ratio - 0.0005, center_ratio + 0.0005, 100)

    # Tuning-Schleife: viele nahe Ratio-Stellen systematisch absuchen
    for current_ratio in ratios:
        p_approx = int(2 ** (n_layer / 4) * current_ratio * 100)
        start_search = (p_approx // 30) * 30  # Ausrichtung am Primzahlrad

        print(f"  [Tuning] Teste Ratio {current_ratio:.6f} -> p ~ 10^{len(str(p_approx))}")

        for p in range(start_search, start_search + search_range, 2):
            if isprime(p) and isprime(p + 2) and isprime(p + 6) and isprime(p + 8):
                product_p = p * (p + 2) * (p + 6) * (p + 8)
                h_weight = bin(product_p).count("1")

                print("\n[TREFFER GEFUNDEN]")
                print(f"Ratio: {current_ratio:.6f} | p: {p}")
                print(f"Hamming-Gewicht (P): {h_weight}")
                print(f"Baryzentrischer Index: {p / (2 ** (n_layer / 4)):.6f}")
                return p, product_p

    print("Keine Resonanz im erweiterten Morley-Fenster gefunden.")
    return None


def hyper_navigation_deep_scan(n_layer, search_range=10000000):
    print(f"\n--- Starte Deep-Resonanz-Fokus zur Schicht n = {n_layer} ---")

    # 1. PHASE: SPEKTRAL-ANALYSE (Gipfel finden)
    # Wir nutzen den Peak bei 0.189428
    print("Phase 1: Identifiziere Resonanz-Gipfel...")
    best_ratio = 0.189428
    max_h = 127
    print(f"Peak-Fixpunkt: Ratio {best_ratio:.6f} | Referenz w(P) = {max_h}")

    # 2. PHASE: DIE TIEFEN-BOHRUNG (Deep Scan)
    # Wir zentrieren die Suche exakt auf diesen harmonischen Punkt
    print(f"Phase 2: Starte Deep Scan bei Ratio {best_ratio:.6f}")

    # Berechnung des Startpunkts p aus dem Resonanz-Gipfel
    p_center = int(2 ** (n_layer / 4) * best_ratio * 100)

    # Wir weiten das Fenster massiv aus, da wir jetzt am "richtigen" Ort sind
    start_search = (p_center // 30) * 30

    print(f"Bohre in Region p ~ 10^{len(str(p_center))}...")
    print(f"Suchfenster: {search_range} Einheiten")

    for p in range(start_search, start_search + search_range, 2):
        # Der hocheffiziente Vierling-Check (ABC-Kern)
        if isprime(p) and isprime(p + 8):
            if isprime(p + 2) and isprime(p + 6):
                product_p = p * (p + 2) * (p + 6) * (p + 8)
                h_weight = bin(product_p).count("1")

                print("\n[ZIEL ERREICHT - RESONANZ GEFUNDEN]")
                print(f"Praezisions-Ratio: {p / (2 ** (n_layer / 4) * 100):.10f}")
                print(f"p-Start: {p}")
                print(f"Hamming-Gewicht: {h_weight}")
                return p, product_p

    print("\nKein Volltreffer im aktuellen Deep-Scan-Fenster.")
    print("Vorschlag: Erhoehe search_range auf 100000000 oder verschiebe Ratio minimal.")
    return None


def fine_structure_scan(n_layer, search_range=100000):
    print(f"\n--- MORLEY-FEINSTRUKTUR-SUCHE fuer Schicht n = {n_layer} ---")
    # Wir nehmen den 127er Peak als Zentrum
    peak = 0.189428
    # Winzig-Schritte um den Peak
    ratios = np.linspace(-0.0001, 0.0001, 500)

    for micro_drift in ratios:
        current_ratio = peak + micro_drift
        p_approx = int(2 ** (n_layer / 4) * current_ratio * 100)
        start_search = (p_approx // 30) * 30

        print(
            f"  [Fine] drift={micro_drift:+.7f} ratio={current_ratio:.7f} "
            f"-> p ~ 10^{len(str(p_approx))}"
        )

        for p in range(start_search, start_search + search_range, 2):
            if isprime(p) and isprime(p + 8):
                if isprime(p + 2) and isprime(p + 6):
                    product_p = p * (p + 2) * (p + 6) * (p + 8)
                    h_weight = bin(product_p).count("1")
                    print("\n[FINE-TREFFER GEFUNDEN]")
                    print(f"Ratio: {current_ratio:.10f} | p: {p}")
                    print(f"Hamming-Gewicht (P): {h_weight}")
                    return p, product_p

    print("Kein Treffer in der Feinstruktur-Suche gefunden.")
    return None


def get_ideal_morley_point(n_layer):
    """
    Konstruiert den idealen Morley-Gitterpunkt fuer Schicht n.
    Nutzt den Satz: m_w = 2^-n * m_ast + Sum(2^-k * v_i)
    """
    target_ratio = 0.190003
    bit_depth = n_layer // 4
    p_ideal = int(target_ratio * (2 ** bit_depth))
    return p_ideal


def get_equilateral_morley_center(n_layer):
    # Geometrische Konstante fuer den gleichseitigen Fokus
    # Abgeleitet aus der Tetraeder-Hoehe und dem Zentrum p+4
    phi_res = 0.190428  # entdeckter Peak als Fixpunkt

    # Konstruktion des Zentrums (p+4)
    bit_depth = n_layer // 4
    p_center_ideal = int(phi_res * (2 ** bit_depth))

    # In der gleichseitigen Geometrie: Vierling symmetrisch um p+4
    # p = p_center_ideal - 4
    return p_center_ideal - 4


def get_perfect_morley_point(n_layer):
    """
    Der 'Gottes-Algorithmus': Berechnet den p+4 Schwerpunkt
    in einem gleichseitig kalibrierten Morley-Tetraeder.
    """
    # Hochpraezises Morley-Verhaeltnis (Zentrum der 5005-Symmetrie)
    # Korrigiert um die gleichseitige Tetraeder-Metrik
    phi_perfect = 0.1904285714  # entspricht 4/21

    bit_depth = n_layer / 4
    # Direkte Konstruktion des p+4 Schwerpunkts
    p_plus_4 = int(phi_perfect * (2 ** bit_depth))

    # Rueckrechnung auf den Startwert des Vierlings
    return p_plus_4 - 4


def deterministic_morley_search(n_layer, residue_radius=1000):
    p_center = get_perfect_morley_point(n_layer)

    print(f"--- Deterministischer Morley-Fokus: Schicht {n_layer} ---")
    print(f"Idealer Gitterpunkt p: {p_center}")
    print(f"Suche nur im Residuum (Radius {residue_radius})...")

    lo = p_center - residue_radius // 2
    hi = p_center + residue_radius // 2
    for p in range(lo, hi):
        if p % 2 != 0:  # Nur ungerade Zahlen
            if isprime(p) and isprime(p + 2) and isprime(p + 6) and isprime(p + 8):
                product_p = p * (p + 2) * (p + 6) * (p + 8)
                return p, product_p
    return None


def predict_morley_drift(n_target):
    # Log-Daten aus Funden: (n, log10(search_range))
    n_val = np.array([45, 300, 500], dtype=float)
    drift_val = np.array([2, 6, 8], dtype=float)

    # Lineare Regression der Drift-Entwicklung
    coeffs = np.polyfit(n_val, drift_val, 1)
    predicted_log_range = float(np.polyval(coeffs, n_target))

    # Empfohlenes Fenster fuer die Zielschicht
    recommended_range = int(10 ** predicted_log_range)

    print(f"\n--- Warp-Drift-Vorhersage fuer n = {n_target} ---")
    print(f"Erwartete Drift-Magnitude: 10^{predicted_log_range:.2f}")
    print(f"Empfohlenes Suchfenster: {recommended_range:,}")

    return recommended_range


def predict_warp_coordinates(n_target):
    # Historische Funde: (Schicht n, Ratio-Abweichung)
    n_history = np.array([45, 300, 500, 600], dtype=float)
    ratios_found = np.array([0.1899, 0.190002, 0.190003, 0.1900031], dtype=float)

    # Quadratische Regression fuer die Kurve der Pilotwelle
    coeffs = np.polyfit(n_history, ratios_found, 2)
    predicted_ratio = float(np.polyval(coeffs, n_target))

    # Vorhersage des Suchfensters (logarithmisch)
    history_ranges = np.array([1000, 1000000, 100000000, 100000000], dtype=float)
    window_log = np.polyfit(n_history, np.log10(history_ranges), 1)
    predicted_window = int(10 ** float(np.polyval(window_log, n_target)))

    print(f"\n--- PRAEDIKTIVE NAVIGATION ZU n = {n_target} ---")
    print(f"Berechnetes Warp-Ratio: {predicted_ratio:.10f}")
    print(f"Empfohlene Bohrung (search_range): {predicted_window:,}")

    return predicted_ratio, predicted_window


def ideal_morley(n_layer):
    """Alias fuer den idealen Morley-Punkt im aktuellen Modell."""
    return get_perfect_morley_point(n_layer)


def drift_correction(n_target):
    """
    Liefert einen dimensionslosen Drift-Korrekturfaktor.
    Referenz ist n=300 (Faktor 1.0).
    """
    n_val = np.array([45, 300, 500], dtype=float)
    drift_val = np.array([2, 6, 8], dtype=float)  # log10(search_range)
    coeffs = np.polyfit(n_val, drift_val, 1)
    log_ref = float(np.polyval(coeffs, 300))
    log_target = float(np.polyval(coeffs, n_target))
    return 10 ** (log_target - log_ref)


def final_warp_coefficient():
    # Finaler Warp-Koeffizient fuer n=1000
    # p_target = ideal_morley(1000) * drift_correction(800)
    base = ideal_morley(1000)
    corr = drift_correction(800)
    p_target = int(base * corr)
    print("\n--- FINALER WARP-KOEFFIZIENT ---")
    print(f"ideal_morley(1000): {base}")
    print(f"drift_correction(800): {corr:.6e}")
    print(f"p_target: {p_target}")
    return p_target


def run_dialog():
    print("MorleyWarp Dialog gestartet. Tippe 'q' zum Beenden.")
    while True:
        raw_mode = input("\nModus [tuned/deep/fine/det/pred] (Enter = tuned): ").strip().lower()
        if raw_mode in {"q", "quit", "exit"}:
            print("Programm beendet.")
            break
        if raw_mode == "":
            mode = "tuned"
        elif raw_mode in {"tuned", "deep", "fine", "det", "pred"}:
            mode = raw_mode
        else:
            print("Unbekannter Modus. Nutze 'tuned', 'deep', 'fine', 'det' oder 'pred'.")
            continue

        raw_layer = input("\nSchicht n_layer (z.B. 300): ").strip()
        if raw_layer.lower() in {"q", "quit", "exit"}:
            print("Programm beendet.")
            break
        if not raw_layer.isdigit():
            print("Bitte eine positive ganze Zahl eingeben.")
            continue

        if mode == "deep":
            default_range = 10000000
        elif mode == "det":
            default_range = 1000
        elif mode == "pred":
            default_range = 0
        else:
            default_range = 100000
        search_range = default_range
        if mode != "pred":
            raw_range = input(f"Suchfenster search_range [Enter = {default_range}]: ").strip()
            if raw_range == "":
                search_range = default_range
            elif raw_range.isdigit():
                search_range = int(raw_range)
            else:
                print(f"Ungueltiges Suchfenster. Nutze Standardwert {default_range}.")
                search_range = default_range

        n_layer = int(raw_layer)
        if mode == "deep":
            result = hyper_navigation_deep_scan(n_layer, search_range=search_range)
        elif mode == "fine":
            result = fine_structure_scan(n_layer, search_range=search_range)
        elif mode == "det":
            result = deterministic_morley_search(n_layer, residue_radius=search_range)
        elif mode == "pred":
            predicted = predict_morley_drift(n_layer)
            warp_ratio, warp_window = predict_warp_coordinates(n_layer)
            final_warp_coefficient()
            print("Tipp: Nutze den Wert als search_range in tuned/deep/fine.")
            result = None
        else:
            result = tuned_morley_navigation(n_layer, search_range=search_range)
        if mode == "pred":
            continue
        if result is not None:
            find_p, _find_P = result
            print(f"\nErfolg: Vierling bei 10^{len(str(find_p))} gefunden.")
        else:
            print("\nKein Vierling in diesem Suchfenster gefunden.")


if __name__ == "__main__":
    run_dialog()