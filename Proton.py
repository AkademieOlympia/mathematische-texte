"""
Analyse: Kernladungszahlen im 24er-System (Plichta/Energiedoku).
Ohne Sage – reines Python.
"""

import math

def analyze_isotopes_and_24():
    """Stabile Kernladungszahlen (Z) im 24er-System."""
    stable_elements = [1, 3, 6, 8, 11, 13, 19, 20, 26, 28, 50, 82]

    element_names = {
        1: "Wasserstoff", 3: "Lithium", 6: "Kohlenstoff", 8: "Sauerstoff",
        11: "Natrium", 13: "Aluminium", 19: "Kalium", 20: "Calcium",
        26: "Eisen", 28: "Nickel", 50: "Zinn", 82: "Blei"
    }

    print("### Analyse: Kernladungszahlen im 24er-System ###\n")
    print(f"{'Z':>4} | {'Name':>12} | {'Z mod 24':>10} | {'Typ'}")
    print("-" * 45)

    for z in stable_elements:
        mod_val = z % 24
        on_ray = "STRAHL" if math.gcd(mod_val, 24) == 1 else "ACHSE"
        name = element_names.get(z, "Element")
        print(f"{z:4} | {name:12} | {mod_val:10} | {on_ray}")


if __name__ == "__main__":
    analyze_isotopes_and_24()
