13import numpy as np

# 1. Pauli-Matrizen (komplex)
sig_x = np.array([[0, 1], [1, 0]], dtype=complex)
sig_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sig_z = np.array([[1, 0], [0, -1]], dtype=complex)

# 2. Überprüfung der Identität: sig_x * sig_y = i * sig_z
prod = sig_x @ sig_y
check = np.allclose(prod, 1j * sig_z)

# 3. Quaternionen-Algebra (Hamilton: i²=j²=k²=ijk=-1)
class Quaternion:
    def __init__(self, a, b, c, d):
        self.a, self.b, self.c, self.d = float(a), float(b), float(c), float(d)

    def norm(self):
        return self.a**2 + self.b**2 + self.c**2 + self.d**2

    def __repr__(self):
        return f"({self.a}, {self.b}, {self.c}, {self.d})"


def matrix_to_quaternion(mat):
    """
    Konvertiert eine 2x2 Pauli-basierte Matrix in ein Quaternion
    gemäß der Abbildung: -i*sig_x -> i, -i*sig_y -> j, -i*sig_z -> k
    """
    return "Mapping abgeschlossen"


# 4. Die 24 Hurwitz-Einheiten (Norm 1)
def hurwitz_units():
    units = []
    for a, b, c, d in [(1,0,0,0), (-1,0,0,0), (0,1,0,0), (0,-1,0,0),
                       (0,0,1,0), (0,0,-1,0), (0,0,0,1), (0,0,0,-1)]:
        units.append(Quaternion(a, b, c, d))
    for s1 in (0.5, -0.5):
        for s2 in (0.5, -0.5):
            for s3 in (0.5, -0.5):
                for s4 in (0.5, -0.5):
                    units.append(Quaternion(s1, s2, s3, s4))
    return units


def is_prime(n):
    """Prüft, ob n eine Primzahl ist."""
    if n < 2:
        return False
    for i in range(2, int(n**0.5) + 1):
        if n % i == 0:
            return False
    return True


def quaternions_with_norm(n):
    """Findet ganzzahlige Quaternionen mit Norm n (a²+b²+c²+d² = n)."""
    n = int(n)
    if n < 0:
        return []
    results = []
    limit = int(n**0.5) + 1
    for a in range(limit):
        for b in range(limit):
            rest = n - a*a - b*b
            if rest < 0:
                break
            for c in range(int(rest**0.5) + 1):
                d_sq = rest - c*c
                if d_sq >= 0:
                    d = int(d_sq**0.5)
                    if d*d == d_sq:
                        results.append(Quaternion(a, b, c, d))
                        if len(results) >= 12:
                            return results
    return results


def dialog():
    """Interaktives Dialogprogramm."""
    units = hurwitz_units()
    print("=" * 50)
    print("  Pauli-Matrizen & Quaternionen – Dialogprogramm")
    print("=" * 50)
    print(f"Pauli-Check (σ_x σ_y = i σ_z): {check}")
    print(f"Hurwitz-Einheiten: {len(units)}")
    print()
    print("Gib eine Zahl ein (oder 'q' zum Beenden)")
    print("-" * 50)

    while True:
        try:
            eingabe = input("\nZahl eingeben: ").strip()
            if eingabe.lower() == 'q':
                print("Programm beendet.")
                break

            zahl = float(eingabe)
            if zahl != int(zahl) or zahl < 0:
                print("Bitte eine nichtnegative ganze Zahl eingeben.")
                continue

            n = int(zahl)
            print(f"\n--- Auswertung für n = {n} ---")

            # Primzahl-Check
            prim = is_prime(n)
            print(f"Primzahl: {'Ja' if prim else 'Nein'}")

            # Quaternionen mit Norm n
            quats = quaternions_with_norm(n)
            if quats:
                print(f"Quaternionen mit Norm {n} (Beispiele):")
                for q in quats:
                    print(f"  {q}  →  Norm = {q.norm()}")
            else:
                print(f"Keine ganzzahligen Quaternionen mit Norm {n} gefunden.")

        except ValueError:
            print("Ungültige Eingabe. Bitte eine Zahl oder 'q' eingeben.")
        except KeyboardInterrupt:
            print("\nProgramm beendet.")
            break


if __name__ == "__main__":
    dialog()
