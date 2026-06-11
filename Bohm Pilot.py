import os
import sys

import numpy as np
import matplotlib

# Ohne GUI: Plot speichern und Skript beendet sofort (Prints laufen immer).
# Fenster: BOHM_PILOT_SHOW=1 python "Bohm Pilot.py"
_SHOW = os.environ.get("BOHM_PILOT_SHOW", os.environ.get("HEUREKA_SHOW", "0")) == "1"
if not _SHOW:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import zeta

# 1. Hilfsfunktionen für Primzahlen
def get_primes(n_max):
    primes = []
    is_prime = [True] * (n_max + 1)
    for p in range(2, n_max + 1):
        if is_prime[p]:
            primes.append(p)
            for i in range(p * p, n_max + 1, p):
                is_prime[i] = False
    return primes

primes = get_primes(8000) # Wir brauchen genug Spielraum für die ersten 1000 Primzahlen

# 2. Suche nach Konstellationen (Vierlinge und Fünflinge)
def find_constellations(primes):
    quads = []  # (p, p+2, p+6, p+8)
    quints = [] # (p, p+2, p+6, p+8, p+12) oder (p, p+4, p+6, p+10, p+12)
    
    for i in range(len(primes) - 5):
        p = primes[i]
        # Check für Vierlinge
        if primes[i+3] == p + 8:
            quads.append(primes[i:i+4])
        # Check für Fünflinge (engste Konstellationen)
        if primes[i+4] == p + 12:
            quints.append(primes[i:i+5])
            
    return quads, quints

quads, quints = find_constellations(primes)

# 3. Die Arithmetische Pilotwelle (Approximation über Riemann-Nullstellen)
# Wir nutzen die ersten 10 imaginären Teile der Nullstellen der Zeta-Funktion
gamma = [14.1347, 21.0220, 25.0109, 30.4249, 32.9351, 37.5862, 40.9187, 43.3271, 48.0052, 49.7738]

def pilot_wave(x):
    # Die "Dichte-Welle" nach Riemann/Bohm
    wave = np.zeros_like(x, dtype=float)
    for g in gamma:
        wave += np.sin(g * np.log(x)) / np.sqrt(x)
    return wave

# --- Visualisierung ---

x_vals = np.linspace(10, 1000, 2000)
wave_vals = pilot_wave(x_vals)

plt.figure(figsize=(15, 8))
plt.plot(x_vals, wave_vals, label='Arithmetische Pilotwelle (Riemann-Interferenz)', color='cyan', alpha=0.6)

# Markiere Vierlinge (ABC-Kerne)
for q in quads[:10]:
    plt.axvline(x=q[0], color='gold', linestyle='--', alpha=0.5, label='Vierling (ABC-Kern)' if q == quads[0] else "")

# Markiere Fünflinge (Instabile Partikel)
for qn in quints[:5]:
    plt.scatter(qn[0], 0, color='red', s=100, zorder=5, label='Fünfling (P5-Zustand)' if qn == quints[0] else "")

plt.title("EAB-Modell: Pilotwellen-Interferenz und Primzahl-Konstellationen")
plt.xlabel("Zahlenraum (x)")
plt.ylabel("Amplituden-Dichte")
plt.legend()
plt.grid(True, which='both', linestyle=':', alpha=0.3)

_out = os.path.join(os.path.dirname(os.path.abspath(__file__)) or ".", "bohm_pilot_wave.png")
plt.savefig(_out, dpi=150, bbox_inches="tight")
print(f"[Bohm Pilot] Plot gespeichert: {_out}", file=sys.stderr, flush=True)

# 4. Numerische Analyse der Zerfälle (P5 -> P4 + C)
print(f"Gefundene ABC-Kerne (Vierlinge): {len(quads)}")
print(f"Gefundene P5-Zustände (Fünflinge): {len(quints)}")
print("-" * 50)

for i in range(min(5, len(quints))):
    p5 = quints[i]
    # Suche den nächsten Vierling als Zerfallsprodukt
    # Wir nehmen den Vierling, der im P5 enthalten ist
    p4 = p5[:4]
    c_residue = p5[4] - sum(p4) # Dies ist eine symbolische Massendifferenz
    eigen_energy_c = sum(p5) - sum(p4) # Arithmetischer Erhaltungssatz
    
    print(f"Reaktion {i+1}: P5{p5} --> ABC{p4} + C-Feld")
    print(f"  Erhaltungssatz: Sum(P5) = {sum(p5)} | Sum(ABC) + C = {sum(p4) + eigen_energy_c}")
    print(f"  Eigenenergie C (Rydberg-Schale): {eigen_energy_c}")
    print("-" * 30)

if _SHOW:
    plt.show()
plt.close()