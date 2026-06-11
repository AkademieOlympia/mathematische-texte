import time
import statistics
from collections import defaultdict
from sage.all import *

def run_comparative_benchmark(cracker, bit_depths=[40, 60, 80, 96], samples=50):
    print(f"\n==========================================================")
    print(f"   RANDOMISIERTER LAUFZEIT-VERGLEICH: E8 vs. MILLER-RABIN")
    print(f"==========================================================")
    print(f"Samples pro Bit-Stufe: {samples}")
    print(f"Zusammensetzung: Zufallsmix aus E, A, B, C Klassen")
    
    # Speicher für Ergebnisse
    results = {}

    for bits in bit_depths:
        print(f"\n>>> Teste {bits}-Bit Zahlen...")
        
        times_mr = []
        times_e8 = []
        success_e8 = 0
        c_wall_hits = 0
        
        # Fortschrittsbalken-Simulation
        for i in range(samples):
            # 1. Erzeuge zufällige Semiprimzahl (p*q)
            # Wir sorgen für echte Zufälligkeit bei den Faktoren
            p = next_prime(ZZ.random_element(2**(bits//2-1), 2**(bits//2)))
            q = next_prime(ZZ.random_element(2**(bits//2-1), 2**(bits//2)))
            n = p * q
            
            # 2. Miller-Rabin Zeitmessung (Referenz)
            t0 = time.perf_counter()
            _ = is_prime(n) # Sage nutzt MR + PSW
            t_mr = time.perf_counter() - t0
            times_mr.append(t_mr)
            
            # 3. E8-Cracker Zeitmessung
            t1 = time.perf_counter()
            factor, meta = cracker.crack(n)
            t_e8 = time.perf_counter() - t1
            times_e8.append(t_e8)
            
            if factor is not None:
                success_e8 += 1
            
            if "Shadow" in str(meta.get("method", "")):
                c_wall_hits += 1

        # Statistik berechnen
        avg_mr = statistics.mean(times_mr) * 1000 # in ms
        avg_e8 = statistics.mean(times_e8) * 1000 # in ms
        
        # Speed-Faktor (Wie viel mal langsamer ist E8?)
        # Schutz vor Division durch Null bei extrem schnellem MR
        slowdown = avg_e8 / avg_mr if avg_mr > 0 else 0
        
        results[bits] = {
            "avg_mr_ms": avg_mr,
            "avg_e8_ms": avg_e8,
            "slowdown": slowdown,
            "success_rate": (success_e8 / samples) * 100,
            "c_wall_ratio": (c_wall_hits / samples) * 100
        }
        
        # Live-Ausgabe pro Stufe
        print(f"   MR (Ø): {avg_mr:.4f} ms")
        print(f"   E8 (Ø): {avg_e8:.4f} ms")
        print(f"   Faktor: {slowdown:.1f}x langsamer als MR")
        print(f"   Trefferquote E8: {results[bits]['success_rate']:.1f}%")
        print(f"   Davon C-Wand (Schatten): {results[bits]['c_wall_ratio']:.1f}% der Fälle")

    # Finale Tabelle
    print(f"\n\n{'Bits':<6} | {'MR Ø (ms)':<12} | {'E8 Ø (ms)':<12} | {'Faktor':<8} | {'Erfolg %':<8}")
    print("-" * 60)
    for bits in bit_depths:
        r = results[bits]
        print(f"{bits:<6} | {r['avg_mr_ms']:<12.4f} | {r['avg_e8_ms']:<12.4f} | {r['slowdown']:<8.1f}x | {r['success_rate']:<8.1f}")

    return results

# --- Ausführung ---
if __name__ == "__main__":
    # Stelle sicher, dass E8PrimeCracker definiert ist
    if 'cracker' not in locals():
        cracker = E8PrimeCracker()
        
    # Wir testen realistische Bereiche für Prototypen
    # 40 Bit = Sofort
    # 80 Bit = Kryptografisch interessant (alte Sicherheit)
    # 96 Bit = Härtetest für Python/Sage Skripte
    data = run_comparative_benchmark(cracker, bit_depths=[40, 64, 80, 96], samples=50)