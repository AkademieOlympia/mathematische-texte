import numpy as np

from hantel import PyrochloreGeometry


class EABCHamiltonian:
    def __init__(self, geometry, j_local=1.0, j_shear=0.5):
        """
        geometry: Eine Instanz der PyrochloreGeometry-Klasse
        j_local: Kopplungskonstante für das 'Katzengewicht' (Ising-analog)
        j_shear: Kopplungskonstante für die topologische Scherung
        """
        self.geo = geometry
        self.j_local = j_local
        self.j_shear = j_shear

    def calculate_bond_energy(self, sub_a, sub_b, data_a, data_b, r_vector):
        """
        Berechnet die Energie einer Bindung (Bond) zwischen zwei Untergittern.

        sub_a, sub_b: Indizes der Untergitter (0 bis 3)
        data_a, data_b: Dictionaries mit Hurwitz-Werten
                        z.B. {'katzen': -1.0, 'shear_norm': 2.0}
        r_vector: Distanzvektor zwischen den beiden Plätzen
        """
        r_len = np.linalg.norm(r_vector)
        if r_len < 1e-5:
            return 0.0

        E_local = self.j_local * (data_a['katzen'] * data_b['katzen'])

        # Projektion der lokalen x-Achsen (analog zur Dipol-Projektion, Eq. S18)
        axis_projection = np.dot(self.geo.x_hat[sub_a], self.geo.x_hat[sub_b])

        avg_shear_norm = (data_a['shear_norm'] + data_b['shear_norm']) / 2.0
        E_shear = self.j_shear * (axis_projection / (r_len**3)) * avg_shear_norm

        return E_local + E_shear


if __name__ == "__main__":
    geo = PyrochloreGeometry(r_nn=2.66)

    eabc_ham = EABCHamiltonian(geo, j_local=1.0, j_shear=0.8)

    state_plus = {'katzen': -1.0, 'shear_norm': 2.0}   # In
    state_minus = {'katzen': 1.0, 'shear_norm': 2.0}  # Out

    tetrahedron_states = {
        0: state_plus,
        1: state_plus,
        2: state_minus,
        3: state_minus,
    }

    # Eckpositionen eines Up-Tetraeders: r_nn * sqrt(3/8) entlang z_hat
    positions = geo.z_hat * (geo.r_nn / 2.0) * np.sqrt(1.5)

    total_tetrahedron_energy = 0.0

    for a in range(4):
        for b in range(a + 1, 4):
            r_vec = positions[b] - positions[a]
            E_bond = eabc_ham.calculate_bond_energy(
                sub_a=a,
                sub_b=b,
                data_a=tetrahedron_states[a],
                data_b=tetrahedron_states[b],
                r_vector=r_vec,
            )
            total_tetrahedron_energy += E_bond
            print(f"Bindung {a + 1}<->{b + 1} Energie: {E_bond:.4f}")

    print("-" * 60)
    print(f"Gesamt-Energie des oktonischen Tetraeders: {total_tetrahedron_energy:.4f}")
