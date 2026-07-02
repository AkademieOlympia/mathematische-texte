# SageMath generator for chunked CSV files of twin primes and prime quadruplets.
# Usage in Sage: load("generate_twin_and_quadruplet_chunks.sage"); export_chunks(10^8, 10^6, "out")
# WARNING: x=10^14 is a major distributed/HPC-scale job and yields ~1.3578e12 twin rows.

import os, csv

def export_chunks(limit, chunk_size, outdir):
    os.makedirs(outdir, exist_ok=True)
    for a in range(2, limit + 1, chunk_size):
        b = min(a + chunk_size - 1, limit)
        # Need primes up to b+8 to detect twins p,p+2 and quadruplets p,p+2,p+6,p+8.
        P = set(prime_range(max(2, a - 8), b + 9))
        twin_path = os.path.join(outdir, f"twins_{a}_{b}.csv")
        quad_path = os.path.join(outdir, f"quadruplets_{a}_{b}.csv")
        with open(twin_path, "w", newline="") as ft:
            w = csv.writer(ft)
            w.writerow(["p", "p_plus_2"])
            for p in sorted(P):
                if a <= p <= b and p + 2 in P:
                    w.writerow([p, p + 2])
        with open(quad_path, "w", newline="") as fq:
            w = csv.writer(fq)
            w.writerow(["p", "p_plus_2", "p_plus_6", "p_plus_8"])
            for p in sorted(P):
                if a <= p <= b and all((p + d) in P for d in (2, 6, 8)):
                    w.writerow([p, p + 2, p + 6, p + 8])
        print(f"done {a}..{b}")
