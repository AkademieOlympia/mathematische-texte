#!/usr/bin/env python3
"""Oktonische Gitterkoppelung: Duo-Tetraeder mit gemeinsamem Knoten 0."""

import numpy as np

from eabc_hamiltonian import EABCHamiltonian
from hantel import PyrochloreGeometry


class PyrochloreDuoTetrahedron:
    def __init__(self, geometry, hamiltonian):
        self.geo = geometry
        self.ham = hamiltonian

    def simulate_deconfinement(self, state_plus, state_minus):
        states_ice_network = {
            0: state_plus,   # gemeinsamer Knoten
            1: state_plus,
            2: state_minus,
            3: state_minus,
            4: state_plus,
            5: state_minus,
            6: state_minus,
        }
        states_with_defect = states_ice_network.copy()
        states_with_defect[0] = state_minus

        E_ground = self._calculate_total_cluster_energy(states_ice_network)
        E_defect = self._calculate_total_cluster_energy(states_with_defect)

        print("=== Oktonische Gitterkoppelung (Duo-Tetraeder) ===")
        print(f"Energie im kollektiven Grundzustand: {E_ground:.4f}")
        print(f"Energie mit lokalem Defekt auf dem Gitterknoten: {E_defect:.4f}")
        print(f"Anregungsenergie (Delta E): {E_defect - E_ground:.4f}")
        print("-" * 60)

        if (E_defect - E_ground) < 4.0:
            print("PHYSIKALISCHE INTERPRETATION: Niedrige Defekt-Energie!")
            print("Deine oktonischen Hurwitz-Ladungen sind dekonfiniert.")
            print("Sie können sich wie die deconfined Spinonen im U(1) Spin-Liquid bewegen.")
        else:
            print("PHYSIKALISCHE INTERPRETATION: Hohe elastische Spannung.")
            print("Die Nicht-Assoziativität führt zu einem Confinement (Ladungs-Einschließung).")

        return E_ground, E_defect, E_defect - E_ground

    def _calculate_total_cluster_energy(self, states):
        total_E = 0.0
        pos = self.geo.z_hat * (self.geo.r_nn / 2.0) * np.sqrt(1.5)
        for a in [0, 1, 2, 3]:
            for b in [0, 1, 2, 3]:
                if a < b:
                    total_E += self.ham.calculate_bond_energy(
                        a, b, states[a], states[b], pos[b] - pos[a]
                    )
        for a in [0, 4, 5, 6]:
            for b in [0, 4, 5, 6]:
                if a < b:
                    sub_a = a if a < 4 else a - 3
                    sub_b = b if b < 4 else b - 3
                    total_E += self.ham.calculate_bond_energy(
                        sub_a,
                        sub_b,
                        states[a],
                        states[b],
                        pos[sub_b % 4] - pos[sub_a % 4],
                    )
        return total_E


if __name__ == "__main__":
    geo = PyrochloreGeometry(r_nn=2.66)
    ham = EABCHamiltonian(geo, j_local=1.0, j_shear=0.8)
    duo = PyrochloreDuoTetrahedron(geo, ham)
    state_plus = {"katzen": -1.0, "shear_norm": 2.0}
    state_minus = {"katzen": 1.0, "shear_norm": 2.0}
    duo.simulate_deconfinement(state_plus, state_minus)
