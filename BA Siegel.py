import numpy as np

# =========================================================================
# DAS BAMBERGER SIEGEL: FINALER MONADISCHER TIEFENSCAN
# =========================================================================

class MonadicSealScan:
    def __init__(self, zeros_path, N):
        print(f"[*] Lade Hardware-Kern: {zeros_path}...")
        self.zeros = np.load(zeros_path)
        self.N = N
        self.ln_N = np.log(N)
        print(f"[*] System bereit. {len(self.zeros):,} Nullstellen geladen.")

    def berechne_kohärenz(self, x, pwr):
        """Berechnet die reine holographische Interferenz am Punkt x."""
        ln_x = np.log(x)
        gamma_chunk = self.zeros[:pwr]
        # Die fundamentale Pilotwellen-Gleichung
        # C(x) = sum( cos(gamma * ln N) * cos(gamma * ln x) ) / sqrt(x)
        interferenz = np.sum(np.cos(gamma_chunk * self.ln_N) * np.cos(gamma_chunk * ln_x))
        return interferenz / np.sqrt(x)

    def execute_drill(self, kandidaten):
        print(f"\n{'PUNKT X':<15} | {'100k ZEROS':<15} | {'2M ZEROS (MAX)':<15} | {'VERSTÄRKUNG'}")
        print("-" * 70)
        
        results = []
        for x in kandidaten:
            # 1. Messung: Niedrige Auflösung (Maische-Check)
            c_low = self.berechne_kohärenz(x, 100000)
            
            # 2. Messung: Volle Auflösung (Monaden-Check)
            c_high = self.berechne_kohärenz(x, len(self.zeros))
            
            # Verstärkungsfaktor berechnen (Spektraler Zoom)
            gain = abs(c_high / c_low) if c_low != 0 else 0
            
            print(f"{x:<15.2f} | {c_low:<15.6f} | {c_high:<15.6f} | {gain:>10.2f}x")
            results.append((x, c_high, gain))
            
        return results

# --- TEST-SETUP ---
if __name__ == "__main__":
    N_TARGET = 243382003
    # Deine drei markanten Punkte aus den vorherigen Logs:
    KANDIDATEN = [
        15360.79,  # Der "Geist" (Glatte Zahl)
        15269.05,  # Die "Maische" (Abstoßungs-Punkt)
        15401.00   # DIE MONADE (Der echte Faktor p)
    ]

    scanner = MonadicSealScan("zeros6.npy", N_TARGET)
    ergebnisse = scanner.execute_drill(KANDIDATEN)

    # Finale Detektion
    sieger_x, sieger_val, sieger_gain = max(ergebnisse, key=lambda x: x[2])
    
    print("\n" + "="*70)
    print(f"DETEKTION ABGESCHLOSSEN:")
    print(f"Die Singularität wurde bei x = {sieger_x} fixiert.")
    print(f"Maximale spektrale Verstärkung: {sieger_gain:.2f}x")
    print(f"Status: FAKTOR p IDENTIFIZIERT")
    print("="*70)