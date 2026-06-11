import numpy as np

class Quaternion:
    def __init__(self, w, x, y, z):
        self.w = float(w)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def norm(self):
        return np.sqrt(self.w**2 + self.x**2 + self.y**2 + self.z**2)

    def normalize(self):
        n = self.norm()
        if n == 0: return self
        return Quaternion(self.w/n, self.x/n, self.y/n, self.z/n)

    def to_rotation_matrix(self):
        """Konvertiert das Quaternion in eine 3x3 Rotationsmatrix (SO(3))"""
        # Nutzt die reinen imaginären Anteile für die räumliche Rotation
        w, x, y, z = self.w, self.x, self.y, self.z
        return np.array([
            [1 - 2*(y**2 + z**2),     2*(x*y - w*z),     2*(x*z + w*y)],
            [    2*(x*y + w*z), 1 - 2*(x**2 + z**2),     2*(y*z - w*x)],
            [    2*(x*z - w*y),     2*(y*z + w*x), 1 - 2*(x**2 + y**2)]
        ])

def slerp(q1, q2, t):
    """Sphärische lineare Interpolation (Scherungsfrei)"""
    dot = q1.w*q2.w + q1.x*q2.x + q1.y*q2.y + q1.z*q2.z
    
    # Kürzesten Weg auf der Hyper-Sphäre wählen
    if dot < 0.0:
        q2 = Quaternion(-q2.w, -q2.x, -q2.y, -q2.z)
        dot = -dot

    if dot > 0.9995:
        # Wenn zu nah beieinander, nutze lineare Interpolation als Grenzfall
        return Quaternion(q1.w + t*(q2.w - q1.w), q1.x + t*(q2.x - q1.x),
                          q1.y + t*(q2.y - q1.y), q1.z + t*(q2.z - q1.z)).normalize()

    theta_0 = np.arccos(dot)
    theta = theta_0 * t
    sin_theta_0 = np.sin(theta_0)
    sin_theta = np.sin(theta)

    s1 = np.sin(theta_0 - theta) / sin_theta_0
    s2 = sin_theta / sin_theta_0

    return Quaternion(
        s1*q1.w + s2*q2.w,
        s1*q1.x + s2*q2.x,
        s1*q1.y + s2*q2.y,
        s1*q1.z + s2*q2.z
    )

def analyze_seismic_shear(R):
    """
    Analysiert die numerische Matrix R auf Deformationen und Scherung.
    Berechnet den metrischen Schertensor und extrahiert die Spurkomponenten.
    """
    # G = R^T * R (Muss im idealen rigiden Fall die Einheitsmatrix sein)
    G = np.dot(R.T, R)
    I = np.eye(3)
    
    # 1. Metrischer Schertensor (Spurfrei)
    spur_G = np.trace(G)
    epsilon_shear = G - (1.0/3.0) * spur_G * I
    frob_norm_shear = np.linalg.norm(epsilon_shear, 'fro')
    
    # 2. Geophysikalische Analogie: Extraktion des Momententensors (M = 0.5 * (R + R^T) - I)
    # Misst die Abweichung von der reinen Starrheit als Deformationskomponente
    M = 0.5 * (R + R.T) - I
    spur_M = np.trace(M) # Isotrope Komponente (Explosion / Implosion)
    
    # Spurfreie Komponente (Scherung + CLVD)
    M_deviatoric = M - (1.0/3.0) * spur_M * I
    
    return frob_norm_shear, spur_M, M_deviatoric

# --- SIMULATION DES EABC-THERMOSTATS ---

if __name__ == "__main__":
    print("=== #Energiedoku: Numerische Scherungs-Analyse im quaternionischen Raum ===")
    
    # Startzustand auf einer Hurwitz-Schale (ideales Einheitsquaternion)
    q_start = Quaternion(0.5, 0.5, 0.5, 0.5) # Klassischer Hurwitz-Gitterpunkt
    
    # Zielzustand für die Transformation
    q_ziel = Quaternion(0.0, 1.0, 0.0, 0.0)
    
    # Wir induzieren künstlich einen numerischen Drift/Fehler (Akkumulations-Scherung)
    print("\n[Simuliere iterative Gitter-Transformationen mit numerischem Drift]")
    q_drift = Quaternion(0.5, 0.52, 0.48, 0.5) # Leicht deformierte Norm
    
    R_ideal = q_start.normalize().to_rotation_matrix()
    R_drift = q_drift.to_rotation_matrix() # Nicht mehr exakt orthogonal
    
    # Analyse der unkorrigierten Drift-Matrix
    shear_energy, isotrop_spur, dev_tensor = analyze_seismic_shear(R_drift)
    
    print(f"-> Akkumulierte Scherenergie (Frobenius-Norm): {shear_energy:.6f}")
    print(f"-> Seismische Spur (Isotrope Volumendehnung):  {isotrop_spur:.6f}")
    print("-> Deviatorischer Momententensor (Scherung/CLVD):")
    print(dev_tensor)
    
    # EABC-Thermostat Regulierungs-Grenzwert
    CRITICAL_LIMIT = 0.01
    print(f"\n[EABC-Thermostat aktiv | Grenzwert: {CRITICAL_LIMIT}]")
    
    if shear_energy > CRITICAL_LIMIT:
        print("🚨 CRITICAL DRIFT DETECTED: Symmetriebruch auf der Hurwitz-Schale!")
        print("-> Aktiviere starrheitserhaltende Renormierung...")
        q_korrigiert = q_drift.normalize()
        R_korrigiert = q_korrigiert.to_rotation_matrix()
        
        new_shear, _, _ = analyze_seismic_shear(R_korrigiert)
        print(f"-> Zustand stabilisiert. Restliche Scherenergie: {new_shear:.1e}")
    else:
        print("✅ Gitterstabilität innerhalb der Toleranz.")
        
    # Demonstration der scherungsfreien Interpolation (SLERP vs LERP Analogie)
    print("\n[Vergleich der Zustandsübergänge bei t = 0.5 (Mitte des Phasenpfads)]")
    t = 0.5
    
    # Scherungsbehaftetes LERP (Linear Blending Artifact / Candy-Wrapper)
    q_lerp_raw = Quaternion(
        (1-t)*q_start.w + t*q_ziel.w,
        (1-t)*q_start.x + t*q_ziel.x,
        (1-t)*q_start.y + t*q_ziel.y,
        (1-t)*q_start.z + t*q_ziel.z
    )
    print(f"-> LERP Norm (Volumenkollaps-Maß): {q_lerp_raw.norm():.4f} (Abweichung von 1.0 erzeugt mechanischen Stress)")
    
    # Scherungsfreies SLERP
    q_slerp = slerp(q_start, q_ziel, t)
    print(f"-> SLERP Norm (Perfekte Rigidität):  {q_slerp.norm():.4f} (Geodätischer Pfad gewahrt)")