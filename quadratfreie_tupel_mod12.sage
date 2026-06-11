#!/usr/bin/env sage
"""
Tupel Modulo 12 Generator
==========================
Dieses Skript erzeugt eine CSV-Datei mit allen Zahlen von 1 bis N,
die genau 4 verschiedene Primfaktoren haben, wobei jeder Primfaktor
aus einer der Restklassen {1, 5, 7, 11} modulo 12 stammt und alle
vier Restklassen vertreten sind. Jede gefundene Zahl wird mit ihren 4 Restklassen in die CSV geschrieben.
"""

from sage.all import *
import csv
import sys

# Erlaubte Restklassen modulo 12
ALLOWED_RESIDUES = [1, 5, 7, 11]

def has_exactly_four_different_residues(n):
    """
    Prüft, ob n das Produkt von genau 4 Primzahlen ist,
    wobei jede Primzahl aus {1, 5, 7, 11} mod 12 stammt
    und alle vier Restklassen {1, 5, 7, 11} vertreten sind.
    
    Das bedeutet: n = p1 * p2 * p3 * p4, wobei:
    - p1, p2, p3, p4 sind Primzahlen (können gleich sein)
    - p1 mod 12, p2 mod 12, p3 mod 12, p4 mod 12 sind in {1, 5, 7, 11}
    - Alle vier Restklassen {1, 5, 7, 11} müssen vertreten sein
    """
    if n == 1:
        return False
    
    fac = factor(n)
    
    # Prüfe, ob genau 4 verschiedene Primfaktoren vorhanden sind
    # (d.h. das Produkt von genau 4 Primzahlen, alle verschieden)
    if len(fac) != 4:
        return False
    
    # Sammle die Restklassen der Primfaktoren
    residues = set()
    for p, e in fac:
        residue = p % 12
        # Prüfe, ob Primfaktor in erlaubten Restklassen ist
        if residue not in ALLOWED_RESIDUES:
            return False
        residues.add(residue)
    
    # Prüfe, ob alle vier Restklassen {1, 5, 7, 11} vertreten sind
    return residues == set(ALLOWED_RESIDUES)

def get_factorization_string(n):
    """
    Erzeugt eine String-Darstellung der Faktorisierung.
    Beispiel: 5^3*17*29^2
    """
    fac = factor(n)
    parts = []
    for p, e in sorted(fac):
        if e == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{e}")
    return "*".join(parts)

def get_residue_sequence(n):
    """
    Gibt die Restklassen der Primfaktoren in der Reihenfolge ihres Erscheinens zurück.
    Da wir nur Zahlen mit genau 4 verschiedenen Primfaktoren betrachten,
    gibt es genau 4 Restklassen.
    Beispiel: n = 5*13*7*11 -> [5, 1, 7, 11] (da 5 mod 12 = 5, 13 mod 12 = 1, 7 mod 12 = 7, 11 mod 12 = 11)
    """
    if n == 1:
        return []
    
    fac = factor(n)
    residues = []
    for p, e in sorted(fac):
        # Jeder Primfaktor kommt genau einmal vor (quadratfrei, genau 4 verschiedene)
        residue = p % 12
        residues.append(residue)
    return residues

def generate_quadratfreie_tupel(n, output_file="quadratfreie_tupel_mod12.csv"):
    """
    Generiert alle Zahlen von 1 bis n mit genau 4 verschiedenen Primfaktoren,
    wobei jeder Primfaktor aus {1, 5, 7, 11} mod 12 stammt und alle vier
    Restklassen vertreten sind.
    
    Args:
        n: Obergrenze (natürliche Zahl)
        output_file: Name der Ausgabe-CSV-Datei
    """
    print(f"=== Tupel Modulo 12 Generator ===")
    print(f"Obergrenze: n = {n}")
    print(f"SageMath Version: {version()}\n")
    
    results = []
    count = 0
    
    print("Durchsuche Zahlen von 1 bis {}...".format(n))
    print("Suche nach Zahlen mit genau 4 verschiedenen Primfaktoren,")
    print("wobei jeder Primfaktor aus {{1, 5, 7, 11}} mod 12 stammt")
    print("und alle vier Restklassen vertreten sind.\n")
    
    # Debug: Zeige ein Beispiel
    test_num = 5005  # 5 * 13 * 7 * 11
    if test_num <= n:
        print(f"Test mit Beispiel: {test_num} = 5 * 13 * 7 * 11")
        print(f"  Restklassen: 5 mod 12 = {5 % 12}, 13 mod 12 = {13 % 12}, 7 mod 12 = {7 % 12}, 11 mod 12 = {11 % 12}")
        print(f"  Prüfung: {has_exactly_four_different_residues(test_num)}\n")
    
    for num in range(1, n + 1):
        # Prüfe, ob num das Produkt von genau 4 Primzahlen ist,
        # wobei alle Primzahlen aus {1, 5, 7, 11} mod 12 stammen
        # und alle vier Restklassen {1, 5, 7, 11} vertreten sind
        if not has_exactly_four_different_residues(num):
            continue
        
        # Zahl gefunden - füge zur Liste hinzu
        residues = get_residue_sequence(num)
        fac_str = get_factorization_string(num)
        
        # Erstelle Dictionary mit natürlicher Zahl und Restklassen
        record = {"n": num}
        # Füge Restklassen in der Reihenfolge ihres Erscheinens hinzu
        for i, r in enumerate(residues, 1):
            record[f"r{i}"] = r
        
        results.append(record)
        count += 1
        
        # Fortschrittsanzeige alle 10000 Zahlen
        if num % 10000 == 0:
            print(f"  Fortschritt: {num}/{n} geprüft, {count} gefunden...")
    
    print(f"\nGefunden: {count} Zahlen mit genau 4 verschiedenen Primfaktoren")
    print(f"         (jeweils einer aus jeder Restklasse {{1, 5, 7, 11}} mod 12)")
    
    if count == 0:
        print("\n⚠️  Warnung: Keine Zahlen gefunden!")
        print("   Mögliche Gründe:")
        print("   - N ist zu klein (benötigt mindestens 5005 = 5*13*7*11)")
        print("   - Prüfen Sie, ob die Bedingungen korrekt sind")
    
    # Immer genau 4 Restklassen (r1, r2, r3, r4)
    fieldnames = ["n", "r1", "r2", "r3", "r4"]
    
    # Schreibe CSV-Datei (auch wenn leer, damit die Datei existiert)
    print(f"\nSchreibe Ergebnisse nach '{output_file}'...")
    try:
        with open(output_file, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            if results:
                writer.writerows(results)
        print(f"✓ Fertig! {count} Einträge in '{output_file}' geschrieben.")
        if count > 0:
            print(f"CSV-Format: Erste Spalte 'n' (natürliche Zahl), dann Restklassen r1, r2, r3, r4")
    except Exception as e:
        print(f"✗ Fehler beim Schreiben der CSV-Datei: {e}")
    
    # Zeige erste 10 Einträge als Beispiel
    if results:
        print("\nErste 10 Einträge:")
        print("-" * 80)
        for i, rec in enumerate(results[:10], 1):
            n_val = rec["n"]
            residues_str = ", ".join([str(rec.get(f"r{j}", "")) for j in range(1, 5)])
            fac_str = get_factorization_string(n_val)
            print(f"{i:3d}. n={n_val:6d} | Restklassen: [{residues_str}] | Faktorisierung: {fac_str}")
        if len(results) > 10:
            print(f"... und {len(results) - 10} weitere Einträge")

# Hauptausführung - funktioniert sowohl von Kommandozeile als auch in SageMath-Konsole
def main():
    """Hauptfunktion für die Ausführung"""
    # Versuche N als Kommandozeilenargument
    n_val = None
    if len(sys.argv) > 1:
        try:
            n_val = int(sys.argv[1])
            print(f"N wurde als Kommandozeilenargument übergeben: {n_val}")
        except ValueError:
            pass
    
    # Falls kein Argument, frage interaktiv
    if n_val is None:
        print("=" * 60)
        print("Tupel Modulo 12 Generator")
        print("=" * 60)
        print("\nDieses Programm findet alle Zahlen von 1 bis N,")
        print("die genau 4 verschiedene Primfaktoren haben,")
        print("wobei jeder Primfaktor aus einer der Restklassen {1, 5, 7, 11} modulo 12 stammt")
        print("und alle vier Restklassen vertreten sind.\n")
        try:
            n_input = input("Bitte geben Sie die Obergrenze N ein: ")
            n_val = int(n_input)
        except (ValueError, EOFError, KeyboardInterrupt):
            print("\nFehler: Ungültige Eingabe oder Abbruch.")
            print("Verwenden Sie: sage quadratfreie_tupel_mod12.sage <N>")
            print("Oder in SageMath-Konsole: generate_quadratfreie_tupel(1000)")
            return
    
    if n_val < 1:
        print("Fehler: N muss >= 1 sein.")
        return
    
    # Standard-Ausgabedatei
    output_file = "quadratfreie_tupel_mod12.csv"
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    
    # Führe Berechnung aus
    generate_quadratfreie_tupel(n_val, output_file)

# Führe main() aus, wenn das Skript direkt aufgerufen wird
if __name__ == "__main__":
    main()
