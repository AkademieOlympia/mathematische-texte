# SageMath-Skript — starten mit: ./run_rundweg.sh  oder  sage -python Rundweg.py
# Nicht mit pyenv/python3: dort fehlt das Modul "sage". Siehe RUN.md.
try:
    from sage.all import Matrix, QQ
except ModuleNotFoundError:
    import sys
    print(
        "Modul 'sage' fehlt. Bitte starten:\n"
        "  ./run_rundweg.sh\n"
        "oder:\n"
        "  sage -python Rundweg.py\n"
        "(Details in RUN.md)"
    )
    sys.exit(1)

def init_kollaps_operators():
    """
    Initialisiert die Operatoren für den Schalenübergang.
    Kombiniert R und K zu einem holonomen Rundweg Gamma.
    """
    # R (Rotation) und K (Komplement) wie zuvor
    R = Matrix(QQ, [[0, 0, 0, 1],
                    [1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0]])
    
    K = Matrix(QQ, [[0, 0, 0, 1],
                    [0, 0, 1, 0],
                    [0, 1, 0, 0],
                    [1, 0, 0, 0]])
    
    # Gamma als kombinierter Operator (Übergang durch die Schale)
    # Erhält die von Ihnen postulierte Symmetrie (E,C) und (A,B)
    Gamma = R * K
    
    return Gamma

def berechne_projektor(Gamma):
    """Berechnet den Mittelungsoperator Pi_Gamma = (I + Gamma + Gamma**2 + Gamma**3) / 4."""
    I = Matrix.identity(QQ, 4)
    G1 = Gamma
    G2 = Gamma**2
    G3 = Gamma**3
    
    # Mittelung über den geschlossenen Zyklus L=4 (nur ** und /, kein ^)
    Pi = (I + G1 + G2 + G3) / 4
    return Pi

def drucke_q_eigenraeume(Pi_Gamma):
    """Gibt rationale Eigenräume von Pi_Gamma über QQ aus."""
    print("=== ARITHMETISCHE STRUKTUR VON PI_GAMMA ÜBER Q ===")
    for eigenwert, raum in Pi_Gamma.eigenspaces_right():
        dim = raum.dimension()
        print(f"Eigenwert λ = {eigenwert}, Dimension = {dim}")
        for idx, basisvektor in enumerate(raum.basis(), start=1):
            koordinaten = tuple(basisvektor)
            print(f"  Basisvektor {idx}: {koordinaten}")
    print("-" * 40)

# Analyse des Kollapses
Gamma = init_kollaps_operators()
Pi_Gamma = berechne_projektor(Gamma)

print("Der Kollaps-Operator Pi_Gamma:")
print(Pi_Gamma)
print("-" * 40)

drucke_q_eigenraeume(Pi_Gamma)

# Test mit einem asymmetrischen Zustand aus dem Kleinschen Kern
# Beispiel: Starkes Ungleichgewicht (E=5, A=2, B=0, C=1)
v_lokal = Matrix(QQ, [[5], [2], [0], [1]])
print("Lokaler Zustand (Fluktuation im Kleinschen Kern):")
print(v_lokal.T)

v_glatt = Pi_Gamma * v_lokal
print("\nProjizierter Zustand nach dem Kollaps (S_glatt):")
print(v_glatt.T)
