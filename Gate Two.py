#!/usr/bin/env sage
from sage.all import *
import itertools

# ==============================================================================
# CLASS: SignedTripleSystem (E, χ)
# ==============================================================================
class SignedTripleSystem:
    def __init__(self, n, triples_dict):
        """
        n: Dimension |E|
        triples_dict: Dictionary {(a, b, c): val} für a < b < c mit val in {-1, 1}.
        """
        self.n = n
        self.chi_dict = {}
        # Vollständige alternierende Fortsetzung auf E^3
        for (a, b, c), val in triples_dict.items():
            for perm, sign in zip(itertools.permutations([a, b, c]), [1, -1, 1, -1, 1, -1]):
                # Vorzeichen der Permutation berechnen
                orig = (a, b, c)
                p_sign = 1 if self._perm_parity(orig, perm) == 0 else -1
                self.chi_dict[perm] = val * p_sign

    def _perm_parity(self, orig, perm):
        # Berechnet die Parität der Permutation
        inv_count = 0
        p = list(perm)
        for i in range(len(p)):
            for j in range(i + 1, len(p)):
                if orig.index(p[i]) > orig.index(p[j]):
                    inv_count += 1
        return inv_count % 2

    def chi(self, a, b, c):
        return self.chi_dict.get((a, b, c), 0)

    # --------------------------------------------------------------------------
    # GATE 1: Induzierte 3-Form & Metrik-Matrix (Nicht-Entartungs-Test)
    # --------------------------------------------------------------------------
    def gate1_metric_matrix(self):
        """
        Erzeugt die Gram-Matrix G der induzierten Metrik aus der 3-Form φ_χ.
        In d=7: G_ij = -1/6 * tr(A_i * A_j) mit A_i = (chi(i, a, b))_{a,b}
        """
        n = self.n
        G = matrix(QQ, n, n)
        
        # Konstruiere die Adjazenz/Kontraktionsmatrizen A_i
        A = []
        for i in range(n):
            Ai = matrix(QQ, n, n)
            for a in range(n):
                for b in range(n):
                    Ai[a, b] = self.chi(i, a, b)
            A.append(Ai)
        
        # Induzierte bilineare Spur-Form (Hitchin-Kontraktion)
        for i in range(n):
            for j in range(n):
                G[i, j] = (A[i] * A[j]).trace()
        
        return G

    def test_gate1(self):
        """
        Prüft Nicht-Entartung (Rang) und Positiv-Definitheit.
        """
        G = self.gate1_metric_matrix()
        det_G = G.det()
        is_full_rank = (G.rank() == self.n)
        is_positive = False
        if is_full_rank:
            eigenvals = G.eigenvalues()
            is_positive = all(ev > 0 for ev in eigenvals)
        return {
            "G_matrix": G,
            "det": det_G,
            "full_rank": is_full_rank,
            "positive_definite": is_positive,
            "pass_gate1": is_full_rank and is_positive
        }

    # --------------------------------------------------------------------------
    # KREUZPRODUKT-EVALUATION: x ×_φ y
    # --------------------------------------------------------------------------
    def cross_product(self, x, y):
        """
        Berechnet das Vektorkreuzprodukt z = x ×_φ y auf Basis der 3-Form χ.
        z_k = sum_{i,j} chi(i, j, k) * x_i * y_j
        """
        n = self.n
        z = [0] * n
        for k in range(n):
            val = 0
            for i in range(n):
                for j in range(n):
                    val += self.chi(i, j, k) * x[i] * y[j]
            z[k] = val
        return vector(z)

    # --------------------------------------------------------------------------
    # GATE 2a (EXPLORER): Empirischer & Symbolischer Kompositionsnorm-Test
    # --------------------------------------------------------------------------
    def test_gate2_explorer(self, num_samples=1000):
        """
        Testet die Kompositionsnorm ||x × y||² = ||x||² ||y||² - g(x,y)²
        empirisch an Zufallsvektoren (Falsifikator).
        """
        g1_res = self.test_gate1()
        if not g1_res["full_rank"]:
            return {"pass_gate2_explorer": False, "reason": "Gate 1 (Nicht-Entartung) fehlgeschlagen."}
        
        G = g1_res["G_matrix"]

        def norm_sq(v):
            return (v * G * v)

        def inner_prod(u, v):
            return (u * G * v)

        # Empirischer Falsifikationstest mit Zufallsvektoren
        for _ in range(num_samples):
            x = vector([QQ.random_element(min=-10, max=10) for _ in range(self.n)])
            y = vector([QQ.random_element(min=-10, max=10) for _ in range(self.n)])
            
            cross_xy = self.cross_product(x, y)
            
            left_side = norm_sq(cross_xy)
            right_side = norm_sq(x) * norm_sq(y) - (inner_prod(x, y))**2
            
            if left_side != right_side:
                return {
                    "pass_gate2_explorer": False,
                    "counterexample": (x, y),
                    "lhs": left_side,
                    "rhs": right_side,
                    "diff": left_side - right_side
                }

        return {"pass_gate2_explorer": True, "samples_tested": num_samples}

    def test_gate2_symbolic(self):
        """
        Prüft die Kompositionsnorm SYMBOLISCH über den Polynomring R[x_0..x_n, y_0..y_n].
        Liefert ein exaktes mathematisches Zertifikat (ALL-QUANTOR ∀x,y).
        """
        g1_res = self.test_gate1()
        if not g1_res["full_rank"]:
            return False, "Gate 1 fehlgeschlagen."

        G = g1_res["G_matrix"]
        
        # Symbolischer Polynomring
        R = PolynomialRing(QQ, 2 * self.n, 'var')
        vars_list = R.gens()
        x_sym = vector(vars_list[:self.n])
        y_sym = vector(vars_list[self.n:])

        def norm_sq_sym(v):
            return (v * G * v)

        def inner_sym(u, v):
            return (u * G * v)

        cross_xy = self.cross_product(x_sym, y_sym)
        
        left_side = norm_sq_sym(cross_xy)
        right_side = norm_sq_sym(x_sym) * norm_sq_sym(y_sym) - (inner_sym(x_sym, y_sym))**2
        
        diff_poly = left_side - right_side
        is_identically_zero = (diff_poly == 0)
        
        return is_identically_zero, diff_poly


# ==============================================================================
# BEISPIEL-TESTS & EXECUTION
# ==============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("GATE-2-EXPLORER: TEST DES OKTONIONEN-TRIPELSYSTEMS (d = 7)")
    print("=" * 70)

    # 1. Das klassische Fano-Ebenen Tripelsystem (7D Oktonionen-Basis)
    # Fano-Tripel: (1,2,3), (1,4,5), (1,7,6), (2,4,6), (2,5,7), (3,4,7), (3,6,5)
    fano_triples = {
        (0, 1, 2): 1,
        (0, 3, 4): 1,
        (0, 6, 5): 1,
        (1, 3, 5): 1,
        (1, 4, 6): 1,
        (2, 3, 6): 1,
        (2, 5, 4): 1
    }

    fano_system = SignedTripleSystem(n=7, triples_dict=fano_triples)

    # Gate 1 Test
    g1_out = fano_system.test_gate1()
    print(f"\n[GATE 1] Rang: {g1_out['G_matrix'].rank()}/{fano_system.n}")
    print(f"[GATE 1] Determinante: {g1_out['det']}")
    print(f"[GATE 1] Pass Gate 1: {g1_out['pass_gate1']}")

    # Gate 2 Explorer (Empirisch)
    g2_exp = fano_system.test_gate2_explorer(num_samples=500)
    print(f"\n[GATE 2 Explorer] Pass (Empirisch): {g2_exp['pass_gate2_explorer']}")

    # Gate 2 Symbolic (Exakter Nachweis ∀x,y)
    is_zero, poly_diff = fano_system.test_gate2_symbolic()
    print(f"[GATE 2 Proof] Symbolisch Identisch Null (∀x,y): {is_zero}")

    print("\n" + "=" * 70)
    print("GATE-2-EXPLORER: TEST EINES DEGENERORTEN / DEFEKTEN SYSTEMSS (d = 5)")
    print("=" * 70)

    # 2. Ein beliebiges 5D Tripelsystem (Soll an Gate 1 oder Gate 2 scheitern)
    defective_triples = {
        (0, 1, 2): 1,
        (1, 2, 3): 1,
        (2, 3, 4): 1
    }
    defective_system = SignedTripleSystem(n=5, triples_dict=defective_triples)

    g1_def = defective_system.test_gate1()
    print(f"\n[GATE 1] Rang: {g1_def['G_matrix'].rank()}/5")
    print(f"[GATE 1] Pass Gate 1: {g1_def['pass_gate1']}")

    g2_def = defective_system.test_gate2_explorer()
    print(f"\n[GATE 2 Explorer] Resultat: {g2_def}")