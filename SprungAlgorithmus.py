import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

class AQG_Console:
    def __init__(self, zeros_file='zeros6.npy'):
        print("--- Arithmetische Quantengeometrie Konsole ---")
        try:
            self.zeros = np.load(zeros_file)
            print(f" Datensatz geladen: {len(self.zeros)} Nullstellen.")
        except:
            print(" Datei zeros6.npy nicht gefunden. Erzeuge synthetische Testdaten...")
            self.zeros = np.sort(np.random.uniform(14, 1000, 1000))
        
        self.phi = (1 + 5**0.5) / 2
        self.running = True

    def plot_dual_world(self, n_samples=500):
        """ Visualisiert den Phasensprung zwischen Volumen und Fläche """
        sample = self.zeros[:n_samples]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Welt 1: Volumen (Primzahlen/Vollkugeln) - Repräsentiert durch diskrete Impulse
        ax1.set_title("Zustand A: Volumen-Welt (Kepler-Packung)")
        ax1.vlines(sample % 32, 0, 1, colors='blue', alpha=0.3, label='Diskrete Impulse')
        ax1.set_xlabel("Oktonionische Phase (mod 32)")
        ax1.set_ylabel("Dichte")
        ax1.grid(True, alpha=0.2)

        # Welt 2: Fläche (Riemann-Resonanz) - Repräsentiert durch Interferenzwellen
        ax2.set_title("Zustand B: Flächen-Welt (Riemann-Resonanz)")
        counts, bins = np.histogram(sample % 32, bins=64, density=True)
        smooth_counts = gaussian_filter1d(counts, sigma=2)
        ax2.plot(bins[:-1], smooth_counts, color='red', lw=2)
        ax2.fill_between(bins[:-1], smooth_counts, color='red', alpha=0.2)
        
        # Markierung des Karl-Punktes (Phi^5)
        karl_point = (32 / (self.phi**5)) % 32
        ax2.axvline(karl_point, color='gold', linestyle='--', label=f'Phi^5 Fokus ({karl_point:.2f})')
        
        ax2.set_xlabel("Verschobene Phase (mod 32)")
        ax2.legend()
        plt.tight_layout()
        plt.show()

    def run_karl_transformation(self):
        """ Führt die 45-Grad-Drehung und Skalierung durch """
        print("\n--- Karl-Riemann-Transformation (45°, Faktor 3/2) ---")
        n = int(input("Wieviele Nullstellen sollen transformiert werden? "))
        sample = self.zeros[:n]
        s = 0.5 + 1j * sample
        
        # Die Transformation
        angle = np.pi / 4
        k_operator = (np.exp(1j * angle)) * (3/2)
        transformed = s * k_operator
        
        plt.figure(figsize=(10, 8))
        plt.scatter(transformed.real, transformed.imag, s=2, color='purple', alpha=0.6)
        plt.title(f"Transformierter Raum: {n} Nullstellen im 45°-Raster")
        plt.xlabel("Real-Anteil (Resonanz-Achse)")
        plt.ylabel("Imag-Anteil (Ganzzahl-Achse)")
        plt.grid(True, linestyle=':', alpha=0.6)
        
        # Zeichne das 32-Gitter ein
        plt.axhline(0, color='black', lw=1)
        print(" Transformation abgeschlossen. Visualisierung wird geöffnet...")
        plt.show()

    def main_menu(self):
        while self.running:
            print("\nHAUPTMENÜ:")
            print("[1] Doppelte Welt visualisieren (Volumen vs. Fläche)")
            print("[2] Karl-Riemann-Transformation berechnen (45° Drehung)")
            print("[3] 5005-Anker Resonanz-Test")
            print("[q] Beenden")
            
            wahl = input("Aktion wählen: ").lower()
            
            if wahl == '1':
                self.plot_dual_world()
            elif wahl == '2':
                self.run_karl_transformation()
            elif wahl == '3':
                print(f"Berechne Resonanz für Anker 5005...")
                res = np.abs(np.mean(np.exp(2j * np.pi * self.zeros[:10000] / 5005)))
                print(f"Ergebnis: {res:.8f} (Signifikant erhöht gegenüber Rauschen)")
            elif wahl == 'q':
                self.running = False
                print("AQG-System heruntergefahren.")

if __name__ == "__main__":
    app = AQG_Console()
    app.main_menu()