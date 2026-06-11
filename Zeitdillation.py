# SageMath 10.5: Kompensation der Zeit-Dilatation
from sage.all import sqrt, bernoulli, RR

def bamberger_zeit_faktor(v_prozent_c):
    """
    Berechnet den Lorentz-Faktor und wendet die 
    arithmetische Heilung an.
    """
    v = v_prozent_c / 100.0
    # Klassischer Lorentz-Faktor
    gamma_klassisch = 1 / sqrt(1 - v**2)
    
    # Arithmetische Heilung: Wir nutzen den Zähler von B_4 (1/30)
    # um die Phasenverschiebung zu berechnen
    b4_korrektur = abs(bernoulli(4))
    # Im Bamberger Modell reduziert die harmonische Kopplung die Dilatation
    gamma_geheilt = gamma_klassisch * (1 - b4_korrektur * v**2)
    
    return gamma_klassisch.n(), gamma_geheilt.n()

v_test = 99.9 # 99.9% der Lichtgeschwindigkeit
klassik, geheilt = bamberger_zeit_faktor(v_test)

print(f"Klassische Zeit-Dilatation bei {v_test}% c: {klassik}")
print(f"Bamberger stabilisierte Eigenzeit: {geheilt}")