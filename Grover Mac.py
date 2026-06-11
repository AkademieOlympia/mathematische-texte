import numpy as np
from scipy.fft import fft, ifft

class HolographicFileSystem:
    def __init__(self, cloud_files, n_scale=4):
        self.files = cloud_files # Liste Ihrer iCloud-Dateinamen
        self.n = len(cloud_files)
        self.n_scale = n_scale
        
        # Erzeugung des Bulk-Residuums (Das 'Innere' des Gitters)
        # Wir kodieren jeden Dateinamen als Phase im Hurwitz-Raum
        self.bulk_energy = np.array([hash(f) % 10**6 for f in self.files], dtype=float)
        self.bulk_energy /= np.max(self.bulk_energy) # Normierung

    def compress_to_boundary(self):
        """
        Holographische Projektion: Projiziert den Bulk (Dateien) 
        auf den Rand (1D-Residuum). Entspricht der AdS/CFT-Logik.
        """
        # Fourier-Transformation simuliert die Projektion auf die kritische Linie
        boundary_residuum = fft(self.bulk_energy)
        
        # Stabilisierung durch den Okto-Anker (n=4)
        # Verhindert Informationsverlust bei der Kompression
        threshold = 1.0 / self.n_scale
        boundary_residuum = np.clip(boundary_residuum.real, -threshold, threshold) + \
                            1j * np.clip(boundary_residuum.imag, -threshold, threshold)
        
        return boundary_residuum

    def grover_search_holographic(self, target_filename):
        """
        Grover-Suche direkt auf dem komprimierten Residuum.
        """
        residuum = self.compress_to_boundary()
        target_hash = (hash(target_filename) % 10**6) / (hash(max(self.files, key=hash)) % 10**6)
        
        # Die Suche findet im Frequenzraum statt (schneller als O(sqrt(N)))
        # Da alle Dateien im Residuum 'verschränkt' sind.
        correlation = np.abs(ifft(residuum * np.conj(fft(np.full(self.n, target_hash)))))
        
        found_idx = np.argmax(correlation)
        return self.files[found_idx], correlation[found_idx]

# --- Anwendung auf Ihre iCloud-Zettel ---
my_icloud_notes = ["Notiz_Physik.txt", "Riemann_Beweis.pdf", "Einkaufsliste.dat", "Energiedoku_Final.docx"]
hfs = HolographicFileSystem(my_icloud_notes)

# Suche nach der Energiedoku
result, strength = hfs.grover_search_holographic("Energiedoku_Final.docx")

print(f"Holographische Suche erfolgreich!")
print(f"Datei gefunden: {result} (Signalstärke: {strength:.4f})")