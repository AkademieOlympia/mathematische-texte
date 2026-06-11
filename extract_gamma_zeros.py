"""
Liest die ersten 100.000 Zeilen aus gamma_zeros_100.txt (oder zeros_gamma.txt)
und speichert sie als gamma_zeros_100k.txt (eine Zahl pro Zeile).
"""
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) or "."
# Quelle: gamma_zeros_100.txt (2M) oder zeros_gamma.txt
for name in ("gamma_zeros_100.txt", "zeros_gamma.txt"):
    candidate = os.path.join(SCRIPT_DIR, name)
    if os.path.isfile(candidate):
        SOURCE = candidate
        break
else:
    SOURCE = os.path.join(SCRIPT_DIR, "gamma_zeros_100.txt")
OUTPUT = os.path.join(SCRIPT_DIR, "gamma_zeros_100k.txt")
N_LINES = 100_000

if not os.path.isfile(SOURCE):
    print(f"Datei nicht gefunden. Bitte gamma_zeros_100.txt oder zeros_gamma.txt in den Ordner legen:")
    print(f"  {SCRIPT_DIR}")
    exit(1)

written = 0
with open(SOURCE, "r") as f_in, open(OUTPUT, "w") as f_out:
    for i, line in enumerate(f_in):
        if i >= N_LINES:
            break
        f_out.write(line)
        written += 1

print(f"Erste {written:,} Zeilen aus {SOURCE}")
print(f"→ gespeichert als {OUTPUT}")
