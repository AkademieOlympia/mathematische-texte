#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ZweiEABC (SageMath Version)
Erzeugt 4-Tupel von Primzahlen aus den Restklassen mod 12.
"""

import sys
import csv
import time
from sage.all import *

GROUP_ORDER = [1, 5, 7, 11]

def get_sage_factorization(n):
    return str(factor(n)).replace(" ", "")

def main():
    # Parameter n aus Argumenten oder Default
    n = 1000
    if len(sys.argv) > 1:
        try:
            n = int(sys.argv[1])
        except ValueError:
            pass
            
    print(f"=== SageMath ZweiEABC (4-Tupel) ===")
    print(f"Goal: n = {n}")
    
    start_time = time.time()
    
    # 1. Primzahlen generieren
    all_primes = prime_range(n + 1)
    
    # 2. Gruppieren nach Restklassen (E=1, A=5, B=7, C=11)
    primes_by_group = {r: [p for p in all_primes if p % 12 == r] for r in GROUP_ORDER}
    
    for r in GROUP_ORDER:
        print(f"  Klasse {r:2d}: {len(primes_by_group[r])} Primzahlen")
    
    # 3. 4-Tupel erzeugen
    print("\nErzeuge 4-Tupel...")
    tuples = []
    
    # Da wir in Sage sind, können wir die Iteration effizienter gestalten
    # Aber für n=1000 ist das verschachtelte Loop okay.
    # Für sehr großes n wäre ein generator besser.
    
    for e in primes_by_group[1]:
        for a in primes_by_group[5]:
            for b in primes_by_group[7]:
                for c in primes_by_group[11]:
                    tuples.append({
                        "E": int(e),
                        "A": int(a),
                        "B": int(b),
                        "C": int(c),
                        "E_fac": get_sage_factorization(e),
                        "A_fac": get_sage_factorization(a),
                        "B_fac": get_sage_factorization(b),
                        "C_fac": get_sage_factorization(c),
                        "sum": int(e + a + b + c),
                        "product": int(e * a * b * c),
                    })
                    
    # 4. Sortieren
    print(f"Sortiere {len(tuples)} Tupel...")
    tuples.sort(key=lambda t: (t["sum"], t["product"], t["E"]))
    
    # 5. CSV Export
    out_prefix = "mod12_single_module"
    out_file = f"{out_prefix}_4tuples.csv"
    fields = ["E", "A", "B", "C", "E_fac", "A_fac", "B_fac", "C_fac", "sum", "product"]
    
    with open(out_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(tuples)
        
    elapsed = time.time() - start_time
    print(f"\n✓ Fertig! {len(tuples)} 4-Tupel in {elapsed:.2f}s erzeugt.")
    print(f"Ausgabe: {out_file}")

if __name__ == "__main__":
    main()