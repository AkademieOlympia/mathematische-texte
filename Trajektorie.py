from sage.all import *

def generate_quat_trajectory(limit):
    """
    Erzeugt eine Liste von 3D-Koordinaten (i, j, k) 
    basierend auf der quaternionischen Primzahlsumme.
    """
    Q = QuaternionAlgebra(QQ, -1, -1)
    current_state = Q(0)
    path = [(0, 0, 0)]
    
    for p in primes(limit + 1):
        if p == 2: continue
        
        # Gewichtung (logarithmisch für bessere Skalierung)
        step_size = 1.0 / log(p)
        
        if p % 4 == 1:
            current_state += step_size * Q([0, 1, 0, 0]) # i-Stoß
        else:
            current_state += step_size * Q([0, 0, 1, 0]) # j-Stoß
            
        # Wir extrahieren die Koeffizienten für die 3D-Darstellung
        # k (Index 3) zeigt hier die 'induzierte' Drift
        coords = current_state.coefficient_tuple()
        path.append((coords[1], coords[2], coords[3]))
        
    return path

# Berechnung für die ersten 2000 Primzahlen
trajectory = generate_quat_trajectory(2000)
print(f"Trajektorie mit {len(trajectory)} Punkten berechnet.")
# In SageMath könnte man dies nun mit line3d(trajectory) plotten.