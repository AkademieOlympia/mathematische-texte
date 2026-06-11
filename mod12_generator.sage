#!/usr/bin/env sage
"""
Mod-12 Einzelmodul-Generator (SageMath Native Version)
======================================================
Dieses Skript erzeugt alle Zahlen ≤ n, deren Primfaktoren ausschließlich 
aus genau einer der mod-12 Klassen {1, 5, 7, 11} stammen.

Es nutzt parallele Verarbeitung und optimierte SageMath-Funktionen.
"""

from sage.all import *
import csv
import sys
import time

# --- Konfiguration ---
GROUP_ORDER = [1, 5, 7, 11]
DEFAULT_N = 1000
OUT_PREFIX = "mod12_single_module"

def get_factorization_string(n):
    """Nutzt SageMath's optimierte Faktorisierung für die String-Darstellung."""
    fac = factor(n)
    return str(fac).replace(" ", "")

@parallel(ncpus=4)
def generate_for_residue(n, primes, residue):
    """
    Parallele Funktion zur Generierung der Zahlen für eine Restklasse.
    """
    records = []
    # Pre-cache für schnellere Index-Abfrage
    prime_to_idx = {p: i + 1 for i, p in enumerate(primes)}
    
    def dfs(start_i, current_val, current_indices):
        if current_val > 1:
            # Wir berechnen die Faktorisierung erst am Ende oder hier
            # Da wir die Faktoren schon kennen, nutzen wir sie direkt
            fac = factor(current_val)
            
            records.append({
                "group_mod12": residue,
                "value": current_val,
                "factorization": str(fac).replace(" ", ""),
                "indices_str": " ".join(map(str, sorted(current_indices))),
                "len": len(current_indices),
                "min_index": min(current_indices),
                "max_index": max(current_indices),
            })
            
        for i in range(start_i, len(primes)):
            p = primes[i]
            if current_val * p > n:
                break
            
            current_indices.append(i + 1)
            dfs(i, current_val * p, current_indices)
            current_indices.pop()

    dfs(0, 1, [])
    return records

def run_generator(n=DEFAULT_N, out_prefix=OUT_PREFIX):
    print(f"=== SageMath Mod-12 Generator ===")
    print(f"Obergrenze: n = {n}")
    print(f"SageMath Version: {version()}")
    
    start_time = time.time()
    
    # 1. Primzahlen generieren
    print(f"Erzeuge Primzahlen bis {n}...")
    primes_up_to_n = prime_range(n + 1)
    
    groups = {r: [] for r in GROUP_ORDER}
    for p in primes_up_to_n:
        r = p % 12
        if r in groups:
            groups[r].append(p)
            
    for r in GROUP_ORDER:
        print(f"  Klasse {r:2d}: {len(groups[r])} Primzahlen")

    # 2. Parallel Zahlen generieren
    print("\nStarte parallele Generierung...")
    # Wir bereiten die Argumente für @parallel vor
    # Format: [(args, kwargs), ...]
    tasks = [((n, groups[r], r), {}) for r in GROUP_ORDER if groups[r]]
    
    all_records = []
    for result in generate_for_residue(tasks):
        # result ist ((n, ps, r), records)
        residue_records = result[1]
        all_records.extend(residue_records)
        print(f"  Klasse {result[0][0][2]:2d}: {len(residue_records)} Zahlen gefunden")

    print(f"\nGesamtanzahl: {len(all_records)} Zahlen erzeugt.")
    
    # 3. Sortieren
    print("Sortiere Ergebnisse...")
    # Sortierung nach Wert
    by_value = sorted(all_records, key=lambda x: (x["value"], x["group_mod12"]))
    
    # 4. CSV Export
    fields = ["group_mod12", "value", "factorization", "indices_str", "len", "min_index", "max_index"]
    
    out_file = f"{out_prefix}_by_value.csv"
    with open(out_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(by_value)
        
    end_time = time.time()
    print(f"✓ Fertig! CSV geschrieben nach: {out_file}")
    print(f"Dauer: {end_time - start_time:.2f} Sekunden")

if __name__ == "__main__":
    # Falls n als Argument übergeben wird
    n_val = DEFAULT_N
    if len(sys.argv) > 1:
        try:
            n_val = int(sys.argv[1])
        except ValueError:
            pass
    run_generator(n_val)
