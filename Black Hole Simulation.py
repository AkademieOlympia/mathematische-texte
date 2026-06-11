import numpy as np

class EABCBlackHoleSimulation:
    def __init__(self, initial_mass, initial_kappa):
        """
        Initialisiert das System im EABC-Kontext.
        initial_mass: Startmasse des Systems
        initial_kappa: Start-Oberflächengravitation (analog zur Temperatur)
        """
        self.mass = float(initial_mass)
        self.kappa = float(initial_kappa)
        self.time_steps = 0
        self.history = []

    def absorb_gravitational_wave(self, wave_energy, eabc_resonance_factor):
        """
        Simuliert die Absorption einer Welle. Im EABC-Modell verändert die 
        Resonanz der Richtungen (e,a,b,c) die innere Struktur (Torsion/Spin),
        was die Oberflächengravitation kappa gegen 0 treibt.
        """
        self.time_steps += 1
        
        # Die Masse wächst durch die Energie der Welle
        self.mass += wave_energy
        
        # Deterministischer Abbau von Kappa durch die eabc-Resonanzkopplung
        # Je näher kappa an 0 kommt, desto stärker wirkt die arithmetische Ordnung
        decay_rate = wave_energy * eabc_resonance_factor * (1.0 + 1.0 / (self.kappa + 0.01))
        
        # Abzug von Kappa (darf nicht negativ werden)
        self.kappa = max(0.0, self.kappa - decay_rate)
        
        # Protokollierung für die Heuristik
        self.history.append((self.time_steps, self.mass, self.kappa))

    def run_simulation(self, wave_sequence, resonance_factor=0.05):
        """
        Läuft durch eine Sequenz von einlaufenden Gravitationswellen.
        """
        print(f"Start der EABC-Simulation...")
        print(f"Schritt 0: Masse = {self.mass:.4f}, Kappa (Temperatur) = {self.kappa:.4f}")
        print("-" * 60)
        
        for i, wave_energy in enumerate(wave_sequence):
            self.absorb_gravitational_wave(wave_energy, resonance_factor)
            
            # Ausgabe der Zwischenschritte
            if i % 2 == 0 or self.kappa == 0:
                print(f"Schritt {self.time_steps}: Welle={wave_energy:.2f} -> "
                      f"Masse={self.mass:.4f}, Kappa={self.kappa:.4f}")
            
            # Abbruchbedingung: Extremaler Zustand (Fixpunkt) in endlicher Zeit erreicht
            if self.kappa == 0.0:
                print("-" * 60)
                print(f"🎯 EXTREMALER ZUSTAND IN ENDLICHER ZEIT ERREICHT!")
                print(f"Fixpunkt bei Schritt {self.time_steps} (Masse: {self.mass:.4f})")
                break
        else:
            print("-" * 60)
            print("Simulation beendet. Extremaler Zustand wurde nicht exakt erreicht.")

# --- Parameter für die numerische Analyse ---

# Eine Sequenz von hochenergetischen Gravitationswellen (Vakuum-Energiefluktuationen)
# Im EABC-Modell entspricht dies einer oszillierenden Energiezufuhr
np.random.seed(42) # Für reproduzierbare Kalibrierung
eingangs_wellen = np.random.uniform(0.5, 1.5, size=20)

# Instanziierung: Masse = 10.0, Start-Kappa = 1.5
sim = EABCBlackHoleSimulation(initial_mass=10.0, initial_kappa=1.5)

# Start der Simulation mit einem definierten EABC-Resonanzfaktor
sim.run_simulation(eingangs_wellen, resonance_factor=0.08)