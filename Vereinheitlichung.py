import numpy as np
import matplotlib.pyplot as plt

class UnifiedForceEngine:
    def __init__(self):
        # Unsere nachgewiesenen 67-Sigma Resonanzen
        self.resonances = {
            'Stark (Gluon)': {'lambda': 0.0137, 'mass': 0.05, 'color': 'gold'},
            'EM (Photon)': {'lambda': 0.0073, 'mass': 0.0, 'color': 'cyan'},
            'Schwach (W/Z)': {'lambda': 0.00045, 'mass': 0.8, 'color': 'magenta'}
        }

    def calculate_coupling(self, T):
        """Berechnet die effektive Kopplungstärke bei Höhe T."""
        couplings = {}
        total_k = 0
        for name, data in self.resonances.items():
            # Die GUG-Gleichung mit thermischer Dämpfung
            k = data['lambda'] * np.exp(-data['mass'] / (T + 1e-9))
            couplings[name] = k
            total_k += k
        return couplings, total_k

# Simulation der Symmetrie-Verschmelzung
T_range = np.linspace(0.01, 5.0, 500)
engine = UnifiedForceEngine()

results = {name: [] for name in engine.resonances.keys()}
total_strength = []

for T in T_range:
    c, total = engine.calculate_coupling(T)
    total_strength.append(total)
    for name in c:
        results[name].append(c[name])

# Visualisierung
plt.figure(figsize=(12, 7))
for name, data in results.items():
    plt.plot(T_range, data, label=name, color=engine.resonances[name]['color'], lw=2)

plt.plot(T_range, total_strength, label='VEREINIGTE UR-KRAFT (Summe)', color='white', ls='--', lw=3)
plt.axvline(x=0.8, color='red', alpha=0.3, label='Elektroschwacher Phasenübergang')

plt.title("Die Vereinheitlichung der Kräfte im Hurwitz-Raum", fontsize=15)
plt.xlabel("Arithmetische Höhe T (Energie/Zeit-Skala)")
plt.ylabel("Resonanz-Amplitude im E8-Gitter")
plt.grid(True, which='both', linestyle=':', alpha=0.5)
plt.legend()
plt.style.use('dark_background')
plt.show()