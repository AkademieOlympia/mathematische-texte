import csv
import gzip
import os
import sys
import time
from pathlib import Path

N = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000_000
here = Path(__file__).resolve().parent
out_path = here / f"eabc_mod12_mod240_{N}.csv.gz"

family = {1: "e", 5: "a", 7: "b", 11: "c"}

count = 0
t0 = time.time()
with gzip.open(out_path, "wt", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["n", "mod12", "familie", "mod240", "is_101_mod240", "is_191_mod240"])
    for n in range(1, N + 1):
        r12 = n % 12
        if r12 in family:
            r240 = n % 240
            w.writerow([n, r12, family[r12], r240, int(r240 == 101), int(r240 == 191)])
            count += 1

size_mb = os.path.getsize(out_path) / (1024 * 1024)
print(f"Datei erstellt: {out_path}")
print(f"Zeilen ohne Header: {count:,}")
print(f"Dateigröße: {size_mb:.2f} MB")
print(f"Laufzeit: {time.time() - t0:.1f} s")