import math

def kalibriere_kennzahlen():
    print("=== STARTE PHYSIKALISCHE RE/RA-ZUORDNUNG UND EICHUNG ===")
    
    M = 113160
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    
    # 1. RATIONALER SEKTOR (Auszug der markanten Schalen bis Radius 4)
    # Format: (e, a, b, c) -> Norm
    rationale_fixpunkte = {
        (1, 0, 0, 0): 1.0,
        (2, 0, 0, 0): 2.0,
        (2, 0, 2, 1): 3.0,
        (4, 0, 0, 0): 4.0,
        (4, 0, 3, 0): 5.0,
        (4, 0, 4, 2): 6.0
    }
    
    # Bestimmung des Impuls-Skalenfaktors Lambda
    # Ziel: Die am tiefsten verdichtete rationale Schale (Norm=6.0) korrespondiert 
    # mit der viskosen Grenze der laminaren Rohrströmung.
    # Re = Lambda * (M / Norm)
    empirisch_re_krit = 2300.0
    norm_basis_re = 6.0
    Lambda_impuls = empirisch_re_krit / (M / norm_basis_re)
    
    print(f"Kalibrierter Impuls-Skalenfaktor \u039b_impuls: {Lambda_impuls:.8f}")
    print("-" * 75)
    print(f"{'EABC-Rational':<15} | {'Norm':<6} | {'M / Norm':<10} | {'Effektive Re':<12} | {'Status'}")
    print("-" * 75)
    
    for koord, norm in sorted(rationale_fixpunkte.items(), key=lambda x: x[1]):
        re_effektiv = Lambda_impuls * (M / norm)
        status = "Laminarer Fixpunkt" if re_effektiv > empirisch_re_krit else "Kritische Grenze (Umschlag)"
        print(f"{str(koord):<15} | {norm:<6.1f} | {int(M/norm):<10} | {re_effektiv:<12.2f} | {status}")
        
    print("\n")
    
    # 2. \u03a6-SEKTOR (Kanonische Zeugen der Radien 2 und 3)
    # Format: (e, a, b, c) -> (Norm, algebraischer Koeffizient alpha von M²/N²)
    # Wir nutzen die exakten Fließkomma-Entsprechungen Ihrer Terminal-Ausgabe
    phi_zeugen = {
        (0, -1, 0, 0): (1.618034, 25610371200.0),   # Norm = \u03a6
        (-1, -1, 0, 0): (1.902113, 7683111360.0),
        (-2, -1, -1, 0): (2.760079, 2186251200.0),
        (0, -2, 0, 0): (3.236068, 6402592800.0),
        (-2, -2, 0, 0): (3.804226, 1920777840.0),
        (-3, -2, 0, 0): (4.412724, 1061893440.0)
    }
    
    # Bestimmung des thermischen Skalenfaktors Lambda
    # Ziel: Der infinitesimale Urzeuge (0, -1, 0, 0) markiert das Einsetzen 
    # der thermalen Konvektionswalzen (Ra_krit \u2248 1708)
    empirisch_ra_krit = 1708.0
    alpha_urzeuge = phi_zeugen[(0, -1, 0, 0)][1]
    Lambda_therm = empirisch_ra_krit / alpha_urzeuge
    
    print(f"Kalibrierter Thermischer Skalenfaktor \u039b_therm: {Lambda_therm:e}")
    print("-" * 75)
    print(f"{'EABC-\u03a6-Sektor':<15} | {'Norm':<8} | {'Algebraischer Kern \u03b1':<20} | {'Effektive Ra':<12}")
    print("-" * 75)
    
    for koord, (norm, alpha) in sorted(phi_zeugen.items(), key=lambda x: x[1][0]):
        ra_effektiv = Lambda_therm * alpha
        print(f"{str(koord):<15} | {norm:<8.5f} | {alpha:<20.1f} | {ra_effektiv:<12.2f}")
        
    print("-" * 75)
    print("=== KALIBRIERUNG ERFOLGREICH BEENDET ===")

if __name__ == "__main__":
    kalibriere_kennzahlen()