import numpy as np
from itertools import product

class HurwitzCosmosEngine:
    def __init__(self, precision=4):
        self.units = self._generate_hurwitz_units()
        self.dimension = 4

    def _generate_hurwitz_units(self):
        """Erzeugt das fundamentale Hurwitz-Alphabet (24 Einheiten)."""
        units = []
        # 8 Lipschitz-Einheiten (Achsen)
        for i in range(4):
            for sign in [-1, 1]:
                unit = [0] * 4
                unit[i] = sign
                units.append(unit)
        # 16 Hurwitz-Halbe (Hyperwürfel-Ecken)
        for signs in product([-0.5, 0.5], repeat=4):
            units.append(list(signs))
        return np.array(units)

    @staticmethod
    def q_multiply(q1, q2):
        """Hamilton-Multiplikation für 4D-Raumzeit-Interaktionen."""
        a1, b1, c1, d1 = q1
        a2, b2, c2, d2 = q2
        return np.array([
            a1*a2 - b1*b2 - c1*c2 - d1*d2,
            a1*b2 + b1*a2 + c1*d2 - d1*c2,
            a1*c2 - b1*d2 + c1*a2 + d1*b2,
            a1*d2 + b1*c2 - c1*b2 + d1*a2
        ])

    def get_littlewood_fluctuation(self, scale=0.01):
        """
        Erzeugt Fluktuationen basierend auf Littlewood-Koeffizienten (+/- 1).
        Diese simulieren die Quanten-Unschärfe im Gitter.
        """
        coeffs = np.random.choice([-1, 1], size=(len(self.units), 4))
        return coeffs * scale

    def evolve_metrics(self, time, mass_center=None, g_intensity=0.1):
        """
        Berechnet den Zustand des Universums zu einem Zeitpunkt t.
        Integriert Expansion (a_t) und Krümmung (Gravitation).
        """
        # Skalenfaktor (Expansion)
        a_t = np.sqrt(time + 0.1)
        
        # Basis-Expansion
        state = self.units * a_t
        
        # Gravitative Krümmung hinzufügen
        if mass_center is not None:
            for i, p in enumerate(state):
                dist_vec = p - mass_center
                dist = np.linalg.norm(dist_vec)
                # Verzerrung der Metrik
                warp = 1 - (g_intensity / (dist + 0.1))
                state[i] = mass_center + dist_vec * warp
        
        # Littlewood-Quantenrauschen hinzufügen
        state += self.get_littlewood_fluctuation(scale=0.05 / (time + 1))
        
        return state

# --- ANWENDUNGSBEISPIEL ---
# Initialisierung des Modells
cosmos = HurwitzCosmosEngine()

# Simulation: Zustand des Universums bei t=5 mit einer Masse im Zentrum
universe_now = cosmos.evolve_metrics(time=5.0, mass_center=np.array([0,0,0,0]))

print("Anzahl aktiver Hurwitz-Knoten:", len(universe_now))
print("Erster Knoten-Vektor (expandiert & gekrümmt):", universe_now[0])