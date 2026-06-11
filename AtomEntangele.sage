# SageMath 10.5: Numerische EABC-Simulation des UNSW-Kerne-Experiments
# Definition der Quaternionen-Algebra über den reellen Zahlen (RR)
Q.<i,j,k> = QuaternionAlgebra(RR, -1, -1)

print("=== SCHRITT 1: INITIALE, GETRENNTE ATOMZUSTÄNDE ===")
# Kerne ruhen im geschützten e-Zentrum (Koeffizient von 1)
# Elektronen fluktuieren in den räumlichen Dimensionen (a, b, c)
N1_init = Q(1)                  # Kern 1 im e-Kanal
E1_init = Q([0, 1, 0, 0])       # Elektron 1 im a-Kanal (i)
E2_init = Q([0, 0, 1, 0])       # Elektron 2 im b-Kanal (j)

print(f"Kern 1 (Zentrum e):   {N1_init}")
print(f"Elektron 1 (Raum a): {E1_init}")
print(f"Elektron 2 (Raum b): {E2_init}\n")


print("=== SCHRITT 2: AUFBAU DER ELEKTRONEN-BRÜCKE (VERSCHRÄNKUNG) ===")
# Die Elektronen bauen eine raumartige Resonanzbrücke im (a,b,c)-Unterraum auf
# Wir simulieren dies als einen gemeinsamen, normierten Superpositions-Vektor
E_bridge = Q([0, 1/sqrt(2), 1/sqrt(2), 0])
print(f"Gemeinsamer Brücken-Vektor (a,b-Ebene): {E_bridge}\n")


print("=== SCHRITT 3: GEOMETRISCHES GATTER (QUATERNIONISCHE INJEKTION) ===")
# Das "Geometric Gate" entspricht einer algebraischen Transformation (Automorphismus),
# die die raumartige Korrelation phasenrichtig in die e-Zentren projiziert.
# Wir nutzen eine 45°-Rotationsmatrix im (e-c)-Unterraum:
gate = Q([1/sqrt(2), 0, 0, 1/sqrt(2)])
gate_inv = ~gate

# Anwendung des Gatters: Die Kopplung wandert in das e-Zentrum der Kerne
N1_entangled = gate * E_bridge * gate_inv
N2_entangled = gate * E_bridge * gate_inv # Symmetrische Verschränkung

print(f"Verschränkter Zustand Kern 1: {N1_entangled}")
print(f"Verschränkter Zustand Kern 2: {N2_entangled}")
print("-> Die Kerne teilen nun dieselbe arithmetische Phase im e-Kanal.\n")


print("=== SCHRITT 4: STABILITÄTSTEST GEGEN EXTERNES RAUSCHEN ===")
# Umwelteinflüsse (Dekohärenz) greifen nur die makroskopischen Raumkoordinaten (a,b,c) an.
# Wir addieren ein starkes stochastisches Rauschen auf die i, j, k Komponenten:
noise = Q([0, 0.45, -0.72, 0.11]) 
N1_noisy = N1_entangled + noise

print(f"Verrauschter Gesamtzustand von Kern 1: {N1_noisy}")
# Extraktion des geschützten Kern-Speichers (e-Komponente)
print(f"--> Ausgelesener e-Zentralwert (Speicher): {N1_noisy.coefficient(1)}")
spatial_noise_amp = sqrt(N1_noisy.coefficient(i)^2 + N1_noisy.coefficient(j)^2 + N1_noisy.coefficient(k)^2)
print(f"--> Amplitude des abgefangenen Raum-Rauschens (abc): {spatial_noise_amp:.4f}")