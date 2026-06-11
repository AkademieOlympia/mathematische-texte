import numpy as np

class PyrochloreGeometry:
    def __init__(self, r_nn=1.0):
        """
        Initialisiert die Geometrie eines Pyrochlor-Gitters basierend auf dem Paper.
        r_nn: Nächste-Nachbar-Distanz (Nearest-neighbor distance)
        """
        self.r_nn = r_nn
        
        # 1. Definition der lokalen z-Achsen für die 4 Untergitter (Eq. S2 / S3)
        # Diese zeigen von den Zentren der Up-Tetraeder nach außen zu den Ecken.
        self.z_hat = np.array([
            [1, 1, 1],
            [1, -1, -1],
            [-1, 1, -1],
            [-1, -1, 1]
        ], dtype=float) / np.sqrt(3)
        
        # 2. Definition der lokalen y-Achsen (Eq. S4)
        self.y_hat = np.array([
            [0, 1, -1],
            [-1, 0, -1],
            [-1, -1, 0],
            [-1, 1, 0]
        ], dtype=float) / np.sqrt(2)
        
        # 3. Berechnung der lokalen x-Achsen via Kreuzprodukt: x = y x z (Seite 9)
        self.x_hat = np.array([np.cross(self.y_hat[a], self.z_hat[a]) for a in range(4)])

    def verify_orthonormality(self):
        """Überprüft, ob die lokalen Koordinatensysteme rechtshändig und orthogonal sind."""
        print("=== Starte Orthonormalitätsprüfung der lokalen Achsen ===")
        for a in range(4):
            # Beträge prüfen (sollten 1 sein)
            len_x = np.linalg.norm(self.x_hat[a])
            len_y = np.linalg.norm(self.y_hat[a])
            len_z = np.linalg.norm(self.z_hat[a])
            
            # Skalarprodukte prüfen (sollten 0 sein)
            dot_xy = np.dot(self.x_hat[a], self.y_hat[a])
            dot_yz = np.dot(self.y_hat[a], self.z_hat[a])
            dot_zx = np.dot(self.z_hat[a], self.x_hat[a])
            
            print(f"Untergitter {a+1}:")
            print(f"  Längen (x, y, z): {len_x:.4f}, {len_y:.4f}, {len_z:.4f}")
            print(f"  Skalarprodukte (xy, yz, zx): {dot_xy:.4f}, {dot_yz:.4f}, {dot_zx:.4f}")
            
    def apply_mirror_plane(self, state_vectors, sublattice_idx):
        """
        Beispiel für eine Symmetrieoperation (Spiegelebene M).
        Im Paper (Eq. S8) bildet die Spiegelebene M: +x -> -x, +y -> +y, +z -> -z ab.
        """
        # Lokale Komponenten extrahieren
        v = state_vectors[sublattice_idx]
        
        # Transformation laut Symmetrie-Klassifikation im Paper
        v_transformed = np.array([
            -v[0],  # x-Komponente invertiert bei Spiegelung
             v[1],  # y-Komponente bleibt invariant
            -v[2]   # z-Komponente invertiert
        ])
        return v_transformed

# ==============================================================================
# TEST-AUSFÜHRUNG
# ==============================================================================
if __name__ == "__main__":
    # Geometrie initialisieren (z.B. mit r_nn = 2.66 Angström für Ce2Sn2O7)
    pyro = PyrochloreGeometry(r_nn=2.66)
    
    # 1. Geometrische Konsistenz des Gitters prüfen
    pyro.verify_orthonormality()
    
    print("\n=== Symmetrie-Test für lokale Pseudospins ===")
    # Ein hypothetischer lokaler Zustand auf Untergitter 0 (Index 0)
    # Format: [v_x, v_y, v_z] repräsentiert die Ausrichtung im lokalen System
    # z.B. reines dipolares Moment entlang z_i oder oktupolares Moment entlang y_i
    local_state = np.array([1.0, 1.0, 1.0]) 
    
    print(f"Ursprünglicher lokaler Vektor: {local_state}")
    
    # Spiegelung anwenden (Methode A aus der vorherigen Antwort)
    transformed_state = pyro.apply_mirror_plane(
        state_vectors={0: local_state}, 
        sublattice_idx=0
    )
    
    print(f"Nach Spiegelung M (Eq. S8):   {transformed_state}")
    print("\nHinweis: Wenn Sie Ihr oktonisches EABC-Modell einbetten, "
          "müssen Sie die Multiplikationstabelle Ihrer oktonischen algebraischen "
          "Elemente so erweitern, dass sie invariant unter diesen Gittertransformationen bleibt.")