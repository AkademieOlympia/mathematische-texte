#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mod-12 Einzelmodul-Generator (SageMath-Python Version)
Optimiert für die Hoffbauer-Antiprimes Forschung.
"""

import sys
import csv
import time
from sage import *

# --- Konfiguration ---
GROUP_ORDER = [1, 5, 7, 11]

def get_sage_factorization(n):
    """Gibt die Faktorisierung als kompakten String zurück, z.B. 5^2*11."""
    return str(factor(n)).replace(" ", "")

def generate_group_numbers(n, primes, residue):
    """
    Generiert alle Zahlen <= n, deren Primfaktoren nur aus 'primes' stammen.
    Nutzt DFS für effiziente Pfadsuche.
    """
    records = []
    
    def dfs(start_i, current_val, indices):
        if current_val > 1:
            records.append({
                "group_mod12": residue,
                "value": int(current_val),
                "factorization": get_sage_factorization(current_val),
                "indices_str": " ".join(map(str, sorted(indices))),
                "len": len(indices),
                "min_index": min(indices),
                "max_index": max(indices),
            })
            
        for i in range(start_i, len(primes)):
            p = primes[i]
            if current_val * p > n:
                break
            
            indices.append(i + 1)
            dfs(i, current_val * p, indices)
            indices.pop()

    dfs(0, 1, [])
    return records

def main():
    # Parameter n aus Argumenten oder Default
    n = 1000
    if len(sys.argv) > 1:
        try:
            n = int(sys.argv[1])
        except ValueError:
            pass
            
    print(f"=== SageMath Mod-12 Generator ===")
    print(f"SageMath Version: {version()}")
    print(f"Ziel: n = {n}")
    sys.stdout.flush()
    
    start_time = time.time()
    
    # 1. Primzahlen mit Sage's prime_range (sehr schnell)
    all_primes = prime_range(n + 1)
    print(f"Primes up to {n}: {len(all_primes)}")
    
    # 2. Gruppieren nach Restklassen
    groups = {r: [p for p in all_primes if p % 12 == r] for r in GROUP_ORDER}
    
    for r in GROUP_ORDER:
        print(f"  Klasse {r:2d}: {len(groups[r])} Primzahlen")
    
    # 3. Generiere Zahlen
    all_records = []
    for r in GROUP_ORDER:
        if not groups[r]: continue
        print(f"Generiere für Klasse {r}...")
        records = generate_group_numbers(n, groups[r], r)
        all_records.extend(records)
        print(f"  -> {len(records)} Zahlen gefunden")
        
    # 4. Sortieren nach Wert
    print("Sortiere Ergebnisse...")
    all_records.sort(key=lambda x: (x["value"], x["group_mod12"]))
    
    # 5. CSV Export
    out_file = "mod12_single_module_by_value.csv"
    fields = ["group_mod12", "value", "factorization", "indices_str", "len", "min_index", "max_index"]
    
    with open(out_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_records)
        
    elapsed = time.time() - start_time
    print(f"\n✓ Fertig! {len(all_records)} Datensätze in {elapsed:.2f}s erzeugt.")
    print(f"Ausgabe: {out_file}")

if __name__ == "__main__":
    main()
