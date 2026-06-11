#!/usr/bin/env python3
"""Thermischer Confinement-Deconfinement-Übergang im EABC-Duo-Tetraeder-Modell."""

import numpy as np
import matplotlib.pyplot as plt

from eabc_hamiltonian import EABCHamiltonian
from hantel import PyrochloreGeometry
from pyrochlore_duo_tetrahedron import PyrochloreDuoTetrahedron


class ThermalEABCSimulator:
    def __init__(self, duo_system):
        self.duo = duo_system
        self.state_plus = {"katzen": -1.0, "shear_norm": 2.0}
        self.state_minus = {"katzen": 1.0, "shear_norm": 2.0}
        _, _, self.delta_E = self.duo.simulate_deconfinement(
            self.state_plus, self.state_minus
        )

    def run_temperature_scan(self, T_max=10.0, steps=100):
        """Boltzmann-Mittel der Defektdichte rho_m(T) mit k_B = 1."""
        temperatures = np.linspace(0.1, T_max, steps)
        defect_densities = []

        for T in temperatures:
            boltzmann_factor = np.exp(-self.delta_E / T)
            avg_defect_density = boltzmann_factor / (1.0 + boltzmann_factor)
            defect_densities.append(avg_defect_density)

        return temperatures, np.array(defect_densities)

    def estimate_critical_temperature(self, T_max=12.0, steps=120):
        T, densities = self.run_temperature_scan(T_max=T_max, steps=steps)
        susceptibility = np.gradient(densities, T)
        idx = np.argmax(susceptibility)
        return T[idx], T, densities, susceptibility


if __name__ == "__main__":
    geo = PyrochloreGeometry(r_nn=2.66)
    ham = EABCHamiltonian(geo, j_local=1.0, j_shear=0.8)
    duo = PyrochloreDuoTetrahedron(geo, ham)
    simulator = ThermalEABCSimulator(duo)

    T_c, T, densities, susceptibility = simulator.estimate_critical_temperature(
        T_max=12.0, steps=120
    )

    print(f"Anregungsenergie Delta E = {simulator.delta_E:.4f}")
    print(f"Kritische Temperatur T_c = {T_c:.4f}")
    print(f"Defektdichte bei T_c: rho_m(T_c) = {densities[np.argmax(susceptibility)]:.4f}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 6), sharex=True)

    ax1.plot(T, densities, "b-", linewidth=2)
    ax1.axvline(T_c, color="r", linestyle="--", label=rf"$T_c = {T_c:.2f}$")
    ax1.set_ylabel(r"$\rho_m$ (Defektdichte)")
    ax1.set_title("Thermischer Confinement-Deconfinement-Uebergang (EABC)")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.plot(T, susceptibility, "g-", linewidth=2)
    ax2.axvline(T_c, color="r", linestyle="--")
    ax2.set_xlabel(r"Temperatur $T$ (Einheiten mit $k_B = 1$)")
    ax2.set_ylabel(r"$d\rho_m/dT$")
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("thermal_transition_eabc.pdf", bbox_inches="tight")
    plt.savefig("thermal_transition_eabc.png", dpi=150, bbox_inches="tight")
    print("Plot gespeichert: thermal_transition_eabc.pdf / thermal_transition_eabc.png")
