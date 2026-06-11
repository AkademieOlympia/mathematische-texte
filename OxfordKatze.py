class Quaternion:
    def __init__(self, r, i, j, k):
        self.r = float(r)
        self.i = float(i)
        self.j = float(j)
        self.k = float(k)

    def __add__(self, other):
        return Quaternion(
            self.r + other.r,
            self.i + other.i,
            self.j + other.j,
            self.k + other.k
        )

    def __sub__(self, other):
        return Quaternion(
            self.r - other.r,
            self.i - other.i,
            self.j - other.j,
            self.k - other.k
        )

    def __mul__(self, other):
        # Nicht-kommutative Multiplikationsregel für Quaternionen
        r = self.r*other.r - self.i*other.i - self.j*other.j - self.k*other.k
        i = self.r*other.i + self.i*other.r + self.j*other.k - self.k*other.j
        j = self.r*other.j - self.i*other.k + self.j*other.r + self.k*other.i
        k = self.r*other.k + self.i*other.j - self.j*other.i + self.k*other.r
        return Quaternion(r, i, j, k)

    def norm(self):
        return self.r**2 + self.i**2 + self.j**2 + self.k**2

    def is_hurwitz(self):
        # Prüft, ob alle Koeffizienten ganzzahlig ODER alle halbzahlig (auf .5 endend) sind
        coeffs = [self.r, self.i, self.j, self.k]
        all_int = all(c.is_integer() for c in coeffs)
        all_half = all((c - 0.5).is_integer() for c in coeffs)
        return all_int or all_half

    def __str__(self):
        def format_term(val, unit):
            if val == 0: return ""
            sign = "+" if val > 0 else "-"
            abs_val = abs(val)
            val_str = f"{abs_val}" if abs_val != 1.0 or unit == "" else ""
            if val_str.endswith(".0"): val_str = val_str[:-2]
            return f" {sign} {val_str}{unit}"

        res = f"{self.r if not self.r.is_integer() else int(self.r)}"
        res += format_term(self.i, "i")
        res += format_term(self.j, "j")
        res += format_term(self.k, "k")
        return res.strip()


def bamberger_model_oxford_simulation():
    print("=== #Energiedoku: EABC Hurwitz-Quaternionen-Simulation (Reines Python) ===")
    print("Analon zum Oxforder Triskeezing-Prozess im Gitter\n")
    
    # 1. Definition der beiden nicht-kommutierenden Operatoren (Basiselemente)
    # Operator A: Linearer Operator (Norm 2)
    A = Quaternion(1, 1, 0, 0)
    # Operator B: Die fundamentale Hurwitz-Einheit (Norm 1)
    B = Quaternion(0.5, 0.5, 0.5, 0.5)
    
    print(f"Operator A (Kopf-Zustand): {A} (Norm: {A.norm()})")
    print(f"Operator B (Fallen-Zustand): {B} (Norm: {B.norm()})")
    print("-" * 60)
    
    # 2. Berechnung der Multiplikationspfade und des Kommutators
    AB = A * B
    BA = B * A
    commutator = AB - BA
    
    print(f"Pfad AB: {AB}")
    print(f"Pfad BA: {BA}")
    print(f"Arithmetischer Kommutator [A, B]: {commutator}")
    print(f"Kommutator-Norm (Scherungs-Energie): {commutator.norm()}")
    print("-" * 60)
    
    # 3. Simulation des Triskeezings (Wechselwirkung 3. Ordnung)
    triskeeze_links = (A * B) * A
    triskeeze_rechts = A * (B * A)
    
    print(f"Triskeezing (Linksassoziativ):  {triskeeze_links}")
    print(f"Triskeezing (Rechtsassoziativ): {triskeeze_rechts}")
    print(f"Resultierende End-Norm: {triskeeze_links.norm()}")
    print("-" * 60)
    
    # 4. EABC-Stabilitätsprüfung (Mid-Circuit Projektion)
    print("EABC-Kohärenz-Prüfung (Projektion):")
    print(f"Ist Endzustand im Hurwitz-Gitter stabil? -> {triskeeze_links.is_hurwitz()}")
    
    # 5. Phasenstruktur (Projektion in die reelle und i-Ebene)
    print("\nProjektion auf die komplexe Ebene (Symmetrie-Analyse):")
    print(f"Reale Komponente (Katzengewicht): {triskeeze_links.r}")
    print(f"Imaginäre i-Komponente (Interferenz): {triskeeze_links.i}")

if __name__ == "__main__":
    bamberger_model_oxford_simulation()