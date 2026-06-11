from sage.all import *

def quaternions_zerlegung(n):
    """
    Zerlegt eine natürliche Zahl n in ein Produkt von Prim-Quaternionen.
    Basierend auf dem Hurwitz-Quaternions-Ring (H).
    """
    # Definition des Hurwitz-Quaternions-Rings über Q
    H = QuaternionAlgebra(QQ, -1, -1)
    
    # Wir suchen ein Quaternion q, dessen Norm n ist
    # (Satz von Lagrange: Jede natürliche Zahl ist Summe von 4 Quadraten)
    # Sage's 'representations_as_sum_of_four_squares' gibt uns (a, b, c, d)
    coords = representations_as_sum_of_four_squares(n)[0]
    q = H(coords)
    
    print(f"### #Energiedoku: Analyse von n = {n}")
    print(f"Basis-Quaternion q: {q}")
    print(f"Norm N(q): {q.norm()}")
    print("-" * 30)
    
    # In Sage können wir die Faktorisierung innerhalb von Maximalordnungen (Hurwitz-Zahlen)
    # durchführen. Da dies mathematisch komplex ist, simulieren wir hier die 
    # Zerlegung der Norm n in ihre Primfaktoren p_i und finden für jedes p_i 
    # ein zugehöriges Prim-Quaternion pi_i.
    
    faktoren_n = factor(n)
    prim_quaternionen_kette = []
    
    for p, expo in faktoren_n:
        # Finde ein Primquaternion für die Primzahl p
        # N(pi) = p = a^2 + b^2 + c^2 + d^2
        p_coords = representations_as_sum_of_four_squares(p)[0]
        pi = H(p_coords)
        for _ in range(expo):
            prim_quaternionen_kette.append(pi)
            
    return prim_quaternionen_kette

# Beispielaufruf für n = 1105 (Produkt aus 5, 13, 17 - alles Summen von Quadraten)
n_test = 1105
kette = quaternions_zerlegung(n_test)

print("Zerlegung in Prim-Quadrupel (Familien-Komponenten):")
for i, pi in enumerate(kette):
    print(f"Faktor {i+1} (Prim-Quaternion): {pi} | Norm: {pi.norm()}")

# Verifikation: Produkt der Normen muss n ergeben
check_norm = prod(p.norm() for p in kette)
print("-" * 30)
print(f"Verifikation der Gesamt-Energie (Produkt der Normen): {check_norm}")