from sage.all import *

def find_prime_quadruplets(count=100):
    quadruplets = []
    q = ZZ(5) # Quadrupletts beginnen oberhalb von 3
    
    while len(quadruplets) < count:
        if q.is_prime() and (q+2).is_prime() and (q+6).is_prime() and (q+8).is_prime():
            quadruplets.append((q, q+2, q+6, q+8))
        q = q + 2 # Nur ungerade Zahlen prüfen
    return quadruplets

# Berechne die ersten 100 Quadrupletts
quads = find_prime_quadruplets(100)

print(f"{'Nr.':<4} | {'q (Familie E)':<14} | {'Produkt P_ABC':<16} | {'3. Wurzel':<14} | {'Delta (q - Wurzel)':<18}")
print("-" * 75)

# Ausgabe der ersten 20 zur Visualisierung im Chat
for i, (q, q2, q6, q8) in enumerate(quads[:20], 1):
    P_ABC = ZZ(q2 * q6 * q8)
    
    # Numerische Berechnung der dritten Wurzel mit hoher Präzision
    root_3 = RealField(64)(P_ABC)^(1/3)
    delta = RealField(64)(q) - root_3
    
    print(f"{i:<4} | {int(q):<14} | {int(P_ABC):<16} | {float(root_3):<14.5f} | {float(delta):+.5f}")

# Statistische Zusammenfassung für die vollen 100 Quadrupletts
deltas = []
for q, q2, q6, q8 in quads:
    P_ABC = ZZ(q2 * q6 * q8)
    root_3 = RealField(64)(P_ABC)^(1/3)
    deltas.append(RealField(64)(q) - root_3)

mean_delta = mean(deltas)
print("\n" + "="*75)
print(f"Statistische Analyse für alle {len(quads)} Quadrupletts:")
print(f"Mittleres Delta (Asymmetrie E vs. ABC): {float(mean_delta):.5f}")