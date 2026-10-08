from sage.all import *

class OctonionicFourLevelSpace:
    """
    Implementierung der Vier-Ebenen-Architektur für oktonionische Experimente
    basierend auf dem eabc-renorm Preprint (PR #143).
    
    Ebene I  : Vektorraum-Einbettung R^8 / R^7 (Linear)
    Ebene II : G2-Kreuzprodukt & Jacobiator auf Im(O) (Malcev-Struktur)
    Ebene III: Vollmultiplikation, Assoziator, Ternutator (Alternative Algebra)
    Ebene IV : Analytische Schnittstelle (Zeta-Schalen / Form)
    """
    def __init__(self):
        # Basis-Algebra: Reelle Oktonionen O via Cayley-Dickson
        self.O = Octonions(RR)
        # Standard-Basis-Vektoren (e0 = 1, e1..e7 = Im(O))
        self.basis = list(self.O.basis())
        self.e0 = self.basis[0]
        self.e_im = self.basis[1:] # e1 bis e7

    # =========================================================================
    # EBENE I: Linear- & Normstrukturen
    # =========================================================================
    def scalar_product(self, x, y):
        """Standard-Skalarprodukt <x,y> auf R^8."""
        return RR(x.real() * y.real() + (x.im() * y.im()).real())

    def norm_sq(self, x):
        """Multiplikative Norm N(x) = ||x||^2 (Ebene III)."""
        return RR((x * x.conjugate()).real())

    def Q_form(self, x):
        """Gittertheoretische quadratische Form Q(x) = 1/2 * ||x||^2 (Ebene I)."""
        return self.norm_sq(x) / RR(2)

    # =========================================================================
    # EBENE II: G2-Kreuzprodukt & Malcev-Geometrie
    # =========================================================================
    def cross_product(self, x, y):
        """
        G2-Kreuzprodukt x × y = 1/2 * Im(x*y - y*x) auf Im(O).
        """
        comm = x * y - y * x
        # Extraktion des Imaginärteils als oktonionisches Element
        return (comm - comm.conjugate()) / RR(4)

    def jacobiator(self, x, y, z):
        """
        Kreuzprodukt-Jacobiator J_x(x,y,z) = (x×y)×z + (y×z)×x + (z×x)×y.
        """
        c_xy_z = self.cross_product(self.cross_product(x, y), z)
        c_yz_x = self.cross_product(self.cross_product(y, z), x)
        c_zx_y = self.cross_product(self.cross_product(z, x), y)
        return c_xy_z + c_yz_x + c_zx_y

    # =========================================================================
    # EBENE III: Alternativ-Algebra, Assoziator & Ternutator
    # =========================================================================
    def associator(self, x, y, z):
        """Oktonionenassoziator Assoc(x,y,z) = (x*y)*z - x*(y*z)."""
        return (x * y) * z - x * (y * z)

    def ternutator(self, x, y, z):
        """
        Antisymmetrisierter ternärer Kommutator T(x,y,z).
        T(x,y,z) = 1/6 * sum_{sigma in S3} sgn(sigma) * x_sigma(1)*(x_sigma(2)*x_sigma(3))
        """
        elems = [x, y, z]
        t_sum = self.O(0)
        
        # Iteration über alle Permutationen der Gruppe S3
        for perm in Permutations(3):
            sgn = perm.sign()
            i, j, k = perm[0]-1, perm[1]-1, perm[2]-1
            p_val = elems[i] * (elems[j] * elems[k])
            t_sum += sgn * p_val
            
        return t_sum / RR(6)

    # =========================================================================
    # EBENE IV: Analytische Zeta-Schicht Interface
    # =========================================================================
    def epstein_zeta_shell_weight(self, n):
        """
        Beispiel-Schnittstelle: Berechnet das E8-Schalen-Gewicht 240*sigma_3(n).
        (Demonstriert die funktionale Entkopplung der Zeta-Ebene)
        """
        if n <= 0:
            return 0
        sigma3 = sum(d**3 for d in divisors(n))
        return 240 * sigma3


# =============================================================================
# EXEKUTIVES AUDIT & EXPERIMENTE
# =============================================================================
def run_paper_experiments():
    space = OctonionicFourLevelSpace()
    print("=== 1. Verifikation der Jacobiator-Proportionalität: J_x = 3/2 * Assoc ===")
    
    # Zufällige rein imaginäre Oktonionen erzeugen
    import random
    def random_im():
        coeffs = [0] + [random.uniform(-2, 2) for _ in range(7)]
        return space.O(coeffs)

    x, y, z = random_im(), random_im(), random_im()
    
    J_val = space.jacobiator(x, y, z)
    Assoc_val = space.associator(x, y, z)
    diff = J_val - RR(1.5) * Assoc_val
    
    print(f"Random Im(O) Tripel Test: ||J_x - 1.5*Assoc|| = {space.norm_sq(diff):.2e}")
    assert space.norm_sq(diff) < 1e-12, "Jacobiator Identität fehlgeschlagen!"
    print("-> Identität J_x = 3/2 * Assoc [PASS]\n")

    print("=== 2. Gegenbeispiel: Ternutator-Non-Proportionalität (Prop 5.2) ===")
    e1, e2, e3 = space.e_im[0], space.e_im[1], space.e_im[2] # Quaternionische Basis
    
    assoc_q = space.associator(e1, e2, e3)
    tern_q = space.ternutator(e1, e2, e3)
    
    print(f"Quaternionen-Tripel (e1, e2, e3):")
    print(f"  Assoc(e1, e2, e3)  = {assoc_q}  (Norm = {space.norm_sq(assoc_q)})")
    print(f"  Ternutator(e1, e2, e3) = {tern_q}  (Norm = {space.norm_sq(tern_q)})")
    
    # Prüfung der Widerlegung von T = c * Assoc
    assert space.norm_sq(assoc_q) < 1e-12 and space.norm_sq(tern_q) > 0.99, "Ternutator Test fehlgeschlagen!"
    print("-> Ternutator T != c * Assoc erfolgreich am Quaternionen-Gegenbeispiel gezeigt [PASS]\n")

    print("=== 3. Ebene I vs. Ebene III Normdistinktion ===")
    q_x = space.Q_form(x)
    n_x = space.norm_sq(x)
    print(f"  Quadratische Form Q(x)   = {q_x:.4f}")
    print(f"  Multiplikative Norm N(x) = {n_x:.4f}")
    print(f"  Verhältnis N(x) / Q(x)   = {n_x / q_x:.1f} (Erwartet: 2.0)")

if __name__ == "__main__":
    run_paper_experiments()