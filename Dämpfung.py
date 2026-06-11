import math


def fmt(x: float, fixed_decimals: int = 15, sci_threshold: float = 1e-6) -> str:
    """
    Formatiert Zahlen gemischt:
    - |x| < sci_threshold (und x != 0) -> wissenschaftlich
    - sonst -> feste Nachkommastellen
    """
    if x != 0.0 and abs(x) < sci_threshold:
        return f"{x:.6e}"
    return f"{x:.{fixed_decimals}f}"


def analysiere_daempfung():
    # Zielwerte im eabc-Modell
    phi = (1 + math.sqrt(5)) / 2
    log_phi = math.log(phi)
    log_bg = math.pi / 6

    print("--- ANATOMIE DER ARITHMETISCHEN DÄMPFUNG ---")
    print(f"Ziel-Energie ln(phi):           {fmt(log_phi)}")
    print(f"Urfeld-Energie pi/6:            {fmt(log_bg)}")
    print(f"Erforderliche Gesamtdämpfung:   {fmt(log_phi - log_bg)}\n")

    # Isoliere die Taylor-Komponenten des Tschebyscheff-Ramanujan-Bias
    bias_linear = math.exp(-math.pi)
    projektion_2 = 0.5 * math.exp(-2 * math.pi)
    fluktuation_5 = math.exp(-5 * math.pi)

    print("Diskrete Kraft-Komponenten (k=1):")
    print(f"  1. Tschebyscheff-Bias (e^-pi):       {fmt(bias_linear * 100, 8)} %")
    print(f"  2. Quadr. Projektion (1/2 * e^-2pi): {fmt(projektion_2 * 100, 8)} %")
    print(f"  3. Ikosaeder-Fluktuation (e^-5pi):   {fmt(fluktuation_5 * 100, 8)} %")

    # Konvergenzverlauf der Skalenkaskade
    aktueller_log = log_bg
    print("\nSkalenkaskade (Dissipativer Kollaps gegen die Struktur):")
    for k in range(1, 4):
        num = 1 + math.exp(-5 * (2 * k - 1) * math.pi)
        den = 1 + math.exp(-(2 * k - 1) * math.pi)
        delta_k = math.log(num / den)
        aktueller_log += delta_k
        fehler = log_phi - aktueller_log
        print(
            f"  Schnitt k={k}: "
            f"Aktuelle Energie = {fmt(aktueller_log)} | "
            f"Restfehler = {fmt(fehler)}"
        )


if __name__ == "__main__":
    analysiere_daempfung()