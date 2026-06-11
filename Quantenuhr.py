#!/usr/bin/env sage
import argparse
from sage.all import *

def simuliere_quanten_pendeluhr(schritte=100, stoerung_amplitude=0.15):
    print("======================================================================")
    # Das EABC-Modell nutzt die Hurwitz-Quaternionen als zugrundeliegendes Koordinatengitter
    print("Initialisiere EABC-Modell auf dem Hurwitz-Quaternionen-Gitter...")
    print("======================================================================")
    
    # Basis-Algebra definieren (Hamilton-Quaternionen über den rationalen Zahlen)
    Q = QuaternionAlgebra(QQ, -1, -1, names='i,j,k')
    i, j, k = Q.gens()
    
    # 1. Definition der fundamentalen EABC-Signaturen als Basis-Quaternionen
    E = Q(1)   # Zeitliche Basis / Skalarer Fluss
    A = i      # Raumdimension A
    B = j      # Raumdimension B
    C = k      # Raumdimension C
    
    # 2. Definition des Pendel-Zustands (Das ungestörte Atom)
    # Wir wählen ein stabiles Hurwitz-Element (z.B. (1 + i + j + k)/2) als Generator
    # des idealen, zyklischen Orbits.
    hurwitz_generator = (1 + i + j + k) / 2
    
    # Zustand des Atoms (aktuelles Pendel-Quaternion)
    q_atom = Q(1)
    
    # Globale Uhr-Variablen
    uhr_ticks = 0
    akkumulierte_entropie = 0.0
    
    print(f"Startzustand des Atoms (Pendel): {q_atom}")
    print(f"Simuliere {schritte} Zeitschritte mit Störungsamplitude {stoerung_amplitude}...\n")
    
    print(f"{'Schritt':<8} | {'Zustand (Atom)':<30} | {'Ptolem. Residuum':<18} | {'Ticks':<6}")
    print("-" * 75)
    
    for schritt in range(1, schritte + 1):
        # --- PHASE 1: Freie Evolution (Das Pendel schwingt arithmetisch) ---
        # Multiplikation mit dem Generator treibt das System auf dem Gitter voran
        q_atom = q_atom * hurwitz_generator
        
        # --- PHASE 2: Umwelt-Fluktuation / Dekohärenz ---
        # Eine kleine, kontinuierliche Störung verzerrt die ideale Gittergeometrie.
        # Dies simuliert das thermische Rauschen oder das Schrotrauschen der Photonen.
        stoerung_vektor = [
            RDF.random_element(-stoerung_amplitude, stoerung_amplitude)
            for _ in range(4)
        ]
        stoerung = Q(stoerung_vektor)
        q_atom_gestoert = q_atom + stoerung
        
        # --- PHASE 3: Die Ptolemäische Hemmung (Messung & Rückwirkung) ---
        # Die Hemmung prüft die geometrische Rigidiät des deformierten Quads.
        # Im EABC-Modell testen wir, wie stark das gestoerte Quaternion von der 
        # idealen algebraischen Norm (Ganzzahligkeit auf dem Hurwitz-Gitter) abweicht.
        
        # Extraktion der EABC-Komponenten des aktuellen Zustands
        e_val, a_val, b_val, c_val = q_atom_gestoert.coefficient_tuple()
        
        # Ptolemäische Balance-Bedingung: Abweichung von der algebraischen Identität
        # Hier nutzen wir das Residuum zur nächstgelegenen Hurwitz-Ganzzahl als Maß
        # für die Deformationsenergie des Gitters.
        norm_wert = q_atom_gestoert.reduced_norm()
        residuum = N(abs(norm_wert - round(norm_wert)))
        
        # Schwellenwert-Kriterium für die Hemmung (Wann rastet das Zahnrad ein?)
        # Wenn das Residuum unter einem kritischen Wert liegt, detektiert die Hemmung
        # eine erfolgreiche "Phasen-Koinzidenz" -> Die Uhr tickt.
        schwellenwert = 0.25
        
        if residuum < schwellenwert:
            uhr_ticks += 1
            # Quanten-Rückwirkung (Backaction): Die Hemmung zwingt das Atom 
            # zurück auf das ideale Gitter (Projektion / Kollaps der Wellenfunktion)
            q_atom = Q([round(x) if x.is_integral() else float(round(2*x))/2 for x in q_atom_gestoert.coefficient_tuple()])
            
            # Die quantisierte Korrektur erzeugt eine Entropie-Kostenkomponente
            akkumulierte_entropie += float(schwellenwert - residuum)
            ticker_ausgabe = "TICK!"
        else:
            # Das System verfehlt die Resonanz, die Zeit "fließt" ohne zu quantisieren
            q_atom = q_atom_gestoert
            ticker_ausgabe = ""
            
        # Ausgabe ausgewählter Schritte zur Analyse
        if schritt % 5 == 0 or ticker_ausgabe == "TICK!":
            zustand_str = f"({float(e_val):.2f})E + ({float(a_val):.2f})A"
            print(f"{schritt:<8} | {zustand_str:<30} | {residuum:<18.5f} | {uhr_ticks:<6} {ticker_ausgabe}")
            
    print("=" * 75)
    print(" SIMULATION BEENDET")
    print("=" * 75)
    print(f"Gemessene Quanten-Ticks       : {uhr_ticks}")
    print(f"Reale arithmetische Schritte  : {schritte}")
    print(f"Generierte Gitter-Entropie    : {akkumulierte_entropie:.4f}")
    
    # Berechnung der Uhrenpräzision im Sinne des thermodynamischen Limits von Brunelli
    if akkumulierte_entropie > 0:
        praezision = uhr_ticks / akkumulierte_entropie
        print(f"Berechnete Uhrenpräzision (η) : {praezision:.4f} Ticks/Entropie")
    else:
        print("Uhrenpräzision (η)            : Nicht definierbar (keine Entropie erzeugt)")
    print("======================================================================")

# Ausführen der Simulation
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simulation der Quantenuhr auf dem Hurwitz-Quaternionen-Gitter."
    )
    parser.add_argument(
        "--N",
        type=int,
        default=40,
        help="Anzahl der arithmetischen Zeitschritte (Default: 40).",
    )
    parser.add_argument(
        "--stoerung",
        type=float,
        default=0.18,
        help="Stoerungsamplitude fuer Dekohärenz (Default: 0.18).",
    )
    args = parser.parse_args()

    simuliere_quanten_pendeluhr(args.N, args.stoerung)