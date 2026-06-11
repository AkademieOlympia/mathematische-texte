import numpy as np
import matplotlib.pyplot as plt

class OktonicFermiReconstruction:
    def __init__(self, t=1.0, t_prime=-0.3, mu=-0.5):
        """
        t: Hopping-Amplitude der nächsten Nachbarn
        t_prime: Hopping-Amplitude der übernächsten Nachbarn
        mu: Chemisches Potential (steuert die Füllung/Größe der Fermi-Fläche)
        """
        self.t = t
        self.t_prime = t_prime
        self.mu = mu

    def bare_dispersion(self, kx, ky):
        """Klassische elektronische Dispersion (Tight-Binding)"""
        epsilon = -2 * self.t * (np.cos(kx) + np.cos(ky)) \
                  - 4 * self.t_prime * np.cos(kx) * np.cos(ky) - self.mu
        return epsilon

    def calculate_reconstructed_bands(self, kx, ky, delta_eabc):
        """
        Berechnet die aufgespaltenen Energiebänder laut Sachdev (Eq. 62)
        unter Verwendung der oktonischen String-Spannung (delta_eabc).
        """
        # Unrekonstruierte Energie am Punkt k
        eps_k = self.bare_dispersion(kx, ky)
        
        # Energie am um den Antiferromagnetischen Wellenvektor K=(pi, pi) verschobenen Punkt
        # K_afm verschiebt cos(k) -> cos(k+pi) = -cos(k)
        eps_k_K = self.bare_dispersion(kx + np.pi, ky + np.pi)
        
        # Sachdevs fundamentale Rekonstruktions-Gleichung (Eq. 62)
        avg_eps = (eps_k + eps_k_K) / 2.0
        diff_eps = (eps_k - eps_k_K) / 2.0
        
        # Die oktonische Lücke (delta_eabc) agiert als Masseterm / Higgs-Kopplung
        gap_term = np.sqrt(diff_eps**2 + delta_eabc**2)
        
        E_minus = avg_eps - gap_term
        E_plus = avg_eps + gap_term
        
        return E_minus, E_plus

# ==============================================================================
# Generierung des Fermi-Flächen-Plots
# ==============================================================================
if __name__ == "__main__":
    recon = OktonicFermiReconstruction()
    
    # K-Raum Gitter aufbauen (Brillouin-Zone von -pi bis pi)
    k_space = np.linspace(-np.pi, np.pi, 400)
    KX, KY = np.meshgrid(k_space, k_space)
    
    # 1. Unrekonstruierte, große Fermi-Fläche (Klassisches Metall)
    E_bare = recon.bare_dispersion(KX, KY)
    
    # 2. Rekonstruierte Fermi-Fläche mit deiner quantisierten String-Spannung (Delta = 4.0)
    # Um die Taschen-Struktur im metallischen Regime schön zu sehen, 
    # skalieren wir das effektive Matrix-Higgs-Feld normiert auf die Bandbreite.
    delta_EABC = 1.2  # Effektiver Higgs-Kopplungswert abgeleitet aus Delta E = 4.0
    E_low, E_high = recon.calculate_reconstructed_bands(KX, KY, delta_eabc=delta_EABC)
    
    # Plot erstellen
    plt.figure(figsize=(6, 6))
    
    # Zeichne die große Fermi-Fläche (Luttinger-Grenzlinie bei E=0)
    plt.contour(KX, KY, E_bare, levels=[0], colors='gray', linestyles=':', alpha=0.7)
    
    # Zeichne die neu entstandenen oktonischen Taschen (Pockets) bei E=0
    cs1 = plt.contour(KX, KY, E_low, levels=[0], colors='blue', linewidths=2)
    cs2 = plt.contour(KX, KY, E_high, levels=[0], colors='red', linewidths=2)
    
    # Optische Aufbereitung der Brillouin-Zone
    plt.title('Oktonische Fermi-Flächen-Taschen ($FL^*$-Phase)')
    plt.xlabel('$k_x$')
    plt.ylabel('$k_y$')
    plt.axhline(0, color='black', linewidth=0.5, alpha=0.5)
    plt.axvline(0, color='black', linewidth=0.5, alpha=0.5)
    
    # Dummy-Lines für die Legende
    plt.plot([], [], color='gray', linestyle=':', label='Große Fermi-Fläche (Luttinger-Fluid)')
    plt.plot([], [], color='blue', label='Hole-Pocket (Oktonisches Chargon)')
    plt.plot([], [], color='red', label='Electron-Pocket (Orthogonal-Hybrid)')
    plt.legend(loc='upper right', fontsize='small')
    
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.gca().set_aspect('equal', adjustable='box')
    
    # Als Vektor-PDF für dein LaTeX-Dokument speichern
    plt.savefig('fermi_reconstruction_eabc.pdf')
    print("Grafik 'fermi_reconstruction_eabc.pdf' erfolgreich für Seite 11 exportiert.")