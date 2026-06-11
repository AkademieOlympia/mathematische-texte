import math
from sympy import isprime

def tuned_morley_navigation(n_layer, search_range=100000):
    """
    Findet Primzahlvierlinge durch direkte Navigation im Pascal-Tetraeder
    basierend auf der Morley-Resonanz-Theorie.
    """
    print(f"\n--- Starte Hyper-Navigation (Tuned) zur Schicht n = {n_layer} ---")

    # Der fundamentale Morley-Anker (aus n=100 / p=637.543.211 abgeleitet)
    base_ratio = 0.190003

    # Tuning-Schleife: wir variieren das Ratio in kleinen Schritten
    for drift in [0, 0.0001, -0.0001, 0.0002, -0.0002, 0.0005, -0.0005]:
        current_ratio = base_ratio + drift
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

def run_dialog():
    print("MorleyWarp Dialog gestartet. Tippe 'q' zum Beenden.")
    while True:
        raw_layer = input("\nSchicht n_layer (z.B. 300): ").strip()
        if raw_layer.lower() in {"q", "quit", "exit"}:
            print("Programm beendet.")
            break
        if not raw_layer.isdigit():
            print("Bitte eine positive ganze Zahl eingeben.")
            continue

        raw_range = input("Suchfenster search_range [Enter = 100000]: ").strip()
        if raw_range == "":
            search_range = 100000
        elif raw_range.isdigit():
            search_range = int(raw_range)
        else:
            print("Ungueltiges Suchfenster. Nutze Standardwert 100000.")
            search_range = 100000

        n_layer = int(raw_layer)
        result = tuned_morley_navigation(n_layer, search_range=search_range)
        if result is not None:
            find_p, _find_P = result
            print(f"\nErfolg: Vierling bei 10^{len(str(find_p))} gefunden.")
        else:
            print("\nKein Vierling in diesem Suchfenster gefunden.")


if __name__ == "__main__":
    run_dialog()