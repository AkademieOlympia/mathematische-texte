# pyright: reportMissingImports=false, reportUndefinedVariable=false
import os
import subprocess
import sys

SAGE_BIN = "/Applications/SageMath-10-8.app/Contents/Frameworks/Sage.framework/Versions/Current/local/bin/sage"

try:
    from sage.all import *
except ModuleNotFoundError:
    if __name__ == "__main__" and os.environ.get("DIFFUSION_REEXEC_UNDER_SAGE") != "1":
        env = os.environ.copy()
        env["DIFFUSION_REEXEC_UNDER_SAGE"] = "1"
        completed = subprocess.run([SAGE_BIN, os.path.abspath(__file__), *sys.argv[1:]], env=env)
        raise SystemExit(completed.returncode)
    raise RuntimeError(
        "Dieses Skript benoetigt SageMath. Starte es mit 'sage Diffusion.py' "
        "oder direkt mit dem SageMath-10.8-Binaer."
    )

def simulate_quaternion_diffusion(steps=10):
    # Definieren des Quaternionen-Algebren-Raums über den rationalen Zahlen
    H = QuaternionAlgebra(QQ, -1, -1)
    i, j, k = H.gens()
    
    # Beispielhafte diskrete Punkte (Gitter) im Raum
    # Wir betrachten Hurwitz-Quaternionen: (a+bi+cj+dk) wobei a,b,c,d alle Ganz oder alle Halbganz sind
    points = [a + b*i + c*j + d*k for a in range(-2, 3) 
                                    for b in range(-2, 3) 
                                    for c in range(-2, 3) 
                                    for d in range(-2, 3)]
    
    def is_quaternion_prime(q):
        # Ein Quaternion ist prim, wenn seine Norm eine Primzahl in Z ist
        n = q.reduced_norm()
        return is_prime(ZZ(n))

    # Initialisierung der Konzentration C (z.B. Gauß-Verteilung im Zentrum)
    concentration = {q: float(exp(-q.reduced_norm())) for q in points}
    
    # Diffusions-Schleife (vereinfachtes Ficksches Gesetz)
    # J = -D * grad(C)
    new_concentration = concentration.copy()
    
    for _ in range(steps):
        for q in points:
            # D ist höher, wenn das Ziel-Quaternion eine Primzahlstruktur hat
            D = 0.5 if is_quaternion_prime(q) else 0.1
            
            # Numerischer Gradient (Differenz zu Nachbarn - vereinfacht)
            # In einem echten Modell würden wir die Dirac-Gleichung auf dem Gitter nutzen
            neighbor_sum = 0
            for dq in [i, -i, j, -j, k, -k]:
                neighbor = q + dq
                if neighbor in concentration:
                    neighbor_sum += concentration[neighbor]
            
            # Update-Regel (Diskretisierte Diffusionsgleichung)
            new_concentration[q] += D * (neighbor_sum / 6.0 - concentration[q])
            
        concentration = new_concentration.copy()

    return concentration

# Start der Simulation
result = simulate_quaternion_diffusion(steps=5)
print(f"Simulation abgeschlossen. Anzahl der berechneten Knoten: {len(result)}")
# Beispielwert eines Primzahl-Knotens anzeigen
H = QuaternionAlgebra(QQ, -1, -1)
i, j, k = H.gens()
example_q = 1 + 1*i + 1*j + 0*k # Norm ist 3 (Primzahl)
if example_q in result:
    value = result[example_q]
    print(f"Konzentration bei Hurwitz-Primzahl {example_q}: {value}")
    print(
        "Interpretation: Es wurden 625 Gitterpunkte im Quaternionenraum simuliert. "
        "Der ausgegebene Wert misst die verbleibende lokale Konzentration nach 5 "
        "Diffusionsschritten am primen Knoten 1 + i + j."
    )
    print(
        "Ein positiver Wert wie hier bedeutet, dass an diesem Knoten noch Substanz "
        "gebunden ist; je groesser der Wert, desto staerker die lokale Anreicherung "
        "im Vergleich zu staerker ausgeduennten Bereichen."
    )