#!/usr/bin/env sage
"""
#Energiedoku - Modul: Nullstellen-Sonden-Test (zeros6.npy)
Validierung der spektralen Daten gegen die kritische Linie und E8-Schalen.
"""

import numpy as np
from sage.all import *

def teste_nullstellen_datei(dateipfad):
    print(f"=== Starte numerischen Test der Datei: {dateipfad} ===")
    
    # 1. Laden des Arrays aus dem lokalen Projektspeicher
    try:
        data = np.load(dateipfad)
    except Exception as e:
        print(f"Fehler beim Laden der Datei: {e}")
        return

    print(f"Array-Struktur: Shape = {data.shape}, Datentyp = {data.dtype}")
    print("-" * 60)
    
    # 2. Mathematische Fallunterscheidung nach Datencharakter
    if np.iscomplexobj(data) or (len(data.shape) > 1 and data.shape[1] == 2):
        # Falls die Nullstellen als komplexe Zahlen oder (Real, Imag)-Paare vorliegen
        if not np.iscomplexobj(data):
            real_parts = data[:, 0]
            imag_parts = data[:, 1]
        else:
            real_parts = np.real(data)
            imag_parts = np.imag(data)
            
        print(f"Anzahl gefundener Nullstellen: {len(data)}")
        print(f"Mittelwert des Realteils: {np.mean(real_parts):.6f}")
        print(f"Standardabweichung des Realteils: {np.std(real_parts):.6f}")
        
        # Riemann-Check: Wie nah liegen die empirischen Nullstellen an Re(s) = 0.5?
        nahe_kritische_linie = np.isclose(real_parts, 0.5, atol=1e-4)
        anzahl_auf_linie = np.sum(nahe_kritische_linie)
        prozent_auf_linie = (anzahl_auf_linie / len(data)) * 100
        
        print(f"\n-> Resonanz-Check (Kritische Linie Re(s) = 1/2):")
        print(f"   {anzahl_auf_linie} von {len(data)} Werten liegen auf der kritischen Linie (Toleranz 1e-4).")
        print(f"   Das entspricht einer Übereinstimmung von {prozent_auf_linie:.2f}%.")
        
        # Ausgabe der ersten spektralen Frequenzen (Imaginärteile)
        print(f"\nDie ersten 5 Resonanzknoten (Imaginärteile sortiert):")
        sorted_imags = np.sort(np.abs(imag_parts))
        for idx, val in enumerate(sorted_imags[:5]):
            print(f"  Resonanz_Knoten_{idx}: {val:.6f}")

    else:
        # Falls die Datei reine reelle Eigenwerte (Energiespektrum des Operators H) enthält
        print(f"Anzahl extrahierter reeller Eigenwerte: {len(data)}")
        print(f"Minimaler Energiewert (E_min): {np.min(data):.6f}")
        print(f"Maximaler Energiewert (E_max): {np.max(data):.6f}")
        
        # Symmetrie-Check für stehende Wellen (+E vs -E)
        symmetrie_fehler = np.abs(np.sum(data))
        print(f"\n-> Symmetrie-Check des Hamilton-Spektrums:")
        print(f"   Gesamtsumme des Spektrums (Soll nahe 0 sein): {symmetrie_fehler:.10e}")
        
    print("=== Test der Datei erfolgreich beendet ===")

# Ausführen auf Ihrem Pfad
if __name__ == "__main__":
    import os
    pfad = os.path.expanduser("~/Projects/zeros6.npy")
    if os.path.exists(pfad):
        teste_nullstellen_datei(pfad)
    else:
        # Alternativ im aktuellen Verzeichnis prüfen
        teste_nullstellen_datei("zeros6.npy")