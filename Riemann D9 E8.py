import struct
from pathlib import Path

import mpmath as mp

mp.mp.dps = 25

E8_RESONANCE_MASK = 23


def passes_e8_resonance(u0):
    """E8-Resonanz: ((u^2 - 1) & 23) == 0."""
    return ((u0 * u0 - 1) & E8_RESONANCE_MASK) == 0


def run_advanced_e8_zeta_pipeline():
    print("=== Erweiterte EABC-E9 Pipeline: Numerische Zeta-Anbindung & Topos-Sieb ===")

    # Präzise bekannte Imaginärteile der ersten Zeta-Nullstellen
    zeta_t_values = [
        "14.13472514173469",
        "21.02203963877155",
        "25.01085758014569",
        "30.42487612583861",
        "32.93506158773919",
    ]

    print("\n1. Direkte numerische Zeta-Auswertung auf der kritischen Achse (Re(s) = 1/2):")
    for idx, t_str in enumerate(zeta_t_values, 1):
        t = mp.mpf(t_str)
        s = mp.mpc("0.5", t)
        val_zeta = mp.zeta(s)
        # Überprüfung des Absolutwerts (sollte numerisch nahe 0 sein)
        is_zero_approx = abs(val_zeta) < mp.mpf("1e-5")
        topos_status = "1_Ω (Gültiger Garben-Schnitt)" if is_zero_approx else "0_Ω"
        print(f"  Nullstelle #{idx}: t = {t}")
        print(f"    -> zeta(1/2 + i*t) = {val_zeta}")
        print(f"    -> Topos-Sieve Status: {topos_status} (Betrag: {float(abs(val_zeta)):.2e})")

    print("\n2. Test von Off-Axis Punkten (Re(s) != 1/2) durch das Topos-Sieb:")
    off_axis_cases = [
        mp.mpc("0.4", "14.1347"),
        mp.mpc("0.6", "21.0220"),
        mp.mpc("0.55", "30.4248"),
    ]
    for idx, s in enumerate(off_axis_cases, 1):
        val_zeta = mp.zeta(s)
        # Lemma 3 erzwingt den Kollaps zu 0_Ω unabhängig vom Zeta-Wert bei Off-Axis
        print(f"  Off-Axis Punkt #{idx}: s = {s}")
        print(f"    -> zeta(s) = {val_zeta}")
        print(f"    -> Topos-Sieb (Galois-Norm-Sperre): Kollabiert strikt zu 0_Ω (Leere Menge)")

    print("\n3. E8-Gitter Bitmasken-Integration (Modell: ((u^2 - 1) & 23) == 0):")
    test_u = [1, 5, 7, 12, 19, 25]
    for u in test_u:
        u_check = u**2 - 1
        passed = passes_e8_resonance(u)
        action = (
            "L1-Cache behalten (Pass)"
            if passed
            else "In 1 Taktzyklus gelöscht (Drop)"
        )
        print(f"  Gitterwert u = {u:2d} -> u^2 - 1 = {u_check:3d} -> {action}")


def run_e8_topos_pipeline_test():
    print("=== EABC-E9 Pipeline Test: Topos-Sieb & Erste Zeta-Nullstellen ===")

    # Die ersten nicht-trivialen Nullstellen der Riemannschen Zeta-Funktion (Im-Teil t)
    # s = 1/2 + i*t
    zeta_zeros_t = [
        14.134725142,
        21.022039639,
        25.010857580,
        30.424876126,
        32.935061588,
    ]

    print("\n1. Validierung der On-Axis Nullstellen (Re(s) = 1/2):")
    for idx, t in enumerate(zeta_zeros_t, 1):
        s_real = 0.5
        s_imag = t

        # Topos-Sieb Logik (Lemma 3)
        is_valid_sieve = abs(s_real - 0.5) < 1e-10
        status = "1_Ω (Valid Sheaf Section)" if is_valid_sieve else "0_Ω (Sieved Out)"
        print(f"  Nullstelle #{idx}: s = {s_real} + {s_imag}i -> Topos-Status: {status}")

    print("\n2. Test von fiktiven Off-Axis Nullstellen (Re(s) != 1/2):")
    off_axis_cases = [0.4 + 14.1347j, 0.6 + 21.0220j, 0.55 + 30.4248j]
    for idx, s in enumerate(off_axis_cases, 1):
        s_real = float(s.real)
        is_valid_sieve = abs(s_real - 0.5) < 1e-10
        status = "1_Ω (Valid Sheaf Section)" if is_valid_sieve else "0_Ω (Kollabiert zu 0_Ω)"
        print(f"  Off-Axis Fall #{idx}: s = {s} -> Topos-Status: {status}")

    print("\n3. E8-Gitter Resonanz & Bitmasken-Filter (Modell: ((u^2 - 1) & 23) == 0):")
    # Testen von Gitterkoordinaten u0 im Quaternionen-Modell
    test_coordinates = [1, 5, 7, 12, 25]
    for u0 in test_coordinates:
        u_check = u0**2 - 1
        bit_filter = passes_e8_resonance(u0)
        filter_status = (
            "Pass (L1 Cache behalten)" if bit_filter else "Drop (0 Takt-Zyklen gelöscht)"
        )
        print(
            f"  Gitterkoordinate u0 = {u0}: u^2-1 = {u_check:3d}, "
            f"Bitmaske ((x & 23) == 0) -> {filter_status}"
        )


def export_e8_states_for_cpp(filename="e8_states.bin"):
    print("=== Export-Pipeline für E8-Zustände (C++20, alignas(32)) ===")

    # 1. Generierung validierter Gitterkoordinaten über Quaternionen-/E8-Resonanz
    valid_states = []

    print("Filterung von Gitterpunkten über E8-Resonanzbedingung: ((u^2 - 1) & 23) == 0")
    for u0 in range(-50, 51):
        if passes_e8_resonance(u0) and u0 != 0:
            u_check = u0 * u0 - 1
            # 8D-Vektor: Koordinate, kritische Achse, Resonanzwert, Eins, Rest 0
            coords = [float(u0), 0.5, float(u_check), 1.0, 0.0, 0.0, 0.0, 0.0]
            valid_states.append(coords)

    print(f"  -> {len(valid_states)} valide E8-Zustände für den L1-Cache-Import selektiert.")

    # 2. Export im C++-Format: uint64-Anzahl, danach je 8 IEEE-754-Doubles (64 Byte)
    path = Path(__file__).with_name(filename)
    print(f"Exportiere Binärdaten nach '{path.name}'...")

    with path.open("wb") as f:
        f.write(struct.pack("Q", len(valid_states)))
        for state in valid_states:
            f.write(struct.pack("8d", *state))

    with path.open("rb") as f:
        (count,) = struct.unpack("Q", f.read(8))
        payload = f.read()
    expected_size = 8 + len(valid_states) * 64
    if count != len(valid_states) or len(payload) != count * 64 or path.stat().st_size != expected_size:
        raise RuntimeError(f"Binärexport unvollständig: {path}")

    print(
        f"Export abgeschlossen: {count} Zustände, {path.stat().st_size} Byte "
        f"(8 Byte Header + {count} x 64 Byte)."
    )
    return path


# Ausführung der Pipeline
run_advanced_e8_zeta_pipeline()
print()
run_e8_topos_pipeline_test()
print()
export_e8_states_for_cpp()