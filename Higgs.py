from sage.all import *

# Konstanten
v = 246.22
m_H = 125.10
m_W = 80.38
m_Z = 91.19

def analyze_ratios():
    """
    Untersucht die Verhältnisse der Bosonen-Massen im Hinblick 
    auf geometrische Faktoren (Pi, e, Wurzeln).
    """
    # Das Verhältnis von Higgs zu VEV
    ratio_hv = m_H / v
    # Das Verhältnis von W zu Z (Weinberg-Winkel Cosinus)
    cos_theta_w = m_W / m_Z
    
    print(f"Verhältnis m_H / v: {float(ratio_hv):.4f} (nahe 1/2)")
    print(f"Weinberg-Cosinus (m_W/m_Z): {float(cos_theta_w):.4f}")
    
    # Arithmetische Hypothese: m_H ~ v / sqrt(phi + 1) wobei phi der goldene Schnitt ist
    phi = (1 + sqrt(5)) / 2
    hyp_mH = v / sqrt(phi + 2)
    print(f"Hypothetische Masse aus Geometrie: {float(hyp_mH):.2f} GeV")

analyze_ratios()