import numpy as np
import matplotlib.pyplot as plt

class NeutrinoCoolingSystem:
    def __init__(self, gluon_resonance=0.0137, kappa=0.05):
        self.Lambda_G = gluon_resonance
        self.kappa = kappa # Kühl-Effizienz der schwachen Wechselwirkung
        self.xi_c = 1.0    # Kritische TMO-Schwelle

    def simulate_stability(self, time_steps, neutrino_density):
        temp = 0.5 # Start-Temperatur (stabil)
        history = []
        
        for t in time_steps:
            # Erhitzung durch die 67-Sigma Feinstruktur
            heating = 0.1 * self.Lambda_G * np.random.normal(1, 0.1)
            
            # Kühlung durch Neutrino-Tunnelung
            cooling = self.kappa * neutrino_density * (temp - 0.1)
            
            # Netto-Änderung
            temp += (heating - cooling)
            
            # Sicherheits-Check: Symmetriebruch?
            stability = "STABLE" if temp < self.xi_c else "BROKEN"
            history.append((temp, stability))
            
        return np.array(history)

# Simulation: Vergleich mit wenig vs. viel Neutrinos
steps = np.linspace(0, 100, 500)
system = NeutrinoCoolingSystem()

low_nu = system.simulate_stability(steps, neutrino_density=0.01)
high_nu = system.simulate_stability(steps, neutrino_density=0.15)

# Visualisierung
plt.figure(figsize=(10, 5))
plt.plot(steps, low_nu[:, 0].astype(float), label="Niedrige Neutrino-Dichte (Risiko)", color='orange')
plt.plot(steps, high_nu[:, 0].astype(float), label="Hohe Neutrino-Dichte (Stabil)", color='cyan')
plt.axhline(y=1.0, color='red', linestyle='--', label="Kritische Schwelle xi_c")
plt.title("Thermodynamik der Littlewood-Dynamik: Neutrinos als Raumzeit-Kühlmittel")
plt.xlabel("Kosmische Zeit (T)")
plt.ylabel("Gitter-Temperatur (Varianz)")
plt.legend()
plt.show()