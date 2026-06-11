"""
#Energiedoku: Performance-Tests für QMBM-Modelle
Miller-Rabin, Trial Division, Makro-Quantencomputer-Resonanz.

Lauf:
  python test_qmbm_performance.py
  python test_qmbm_performance.py --quick   # Weniger Läufe, schneller
"""

import time
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import numpy as np


# --- Miller-Rabin (probabilistischer Primzahltest) ---
def _miller_rabin_witness(a, n, r, d):
    """Ein Zeuge für zusammengesetztes n."""
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return False
    for _ in range(r - 1):
        x = pow(x, 2, n)
        if x == n - 1:
            return False
    return True


def miller_rabin(n, k=10):
    """
    Miller-Rabin Primzahltest (probabilistisch).
    k Runden: Fehlerwahrscheinlichkeit ≤ 4^(-k).
    """
    n = int(n)
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
    # Kleine Primzahlen als Basis (deterministisch für n < 2^64 mit k=12)
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37][:k]
    for a in bases:
        if a >= n:
            continue
        if _miller_rabin_witness(a, n, r, d):
            return False
    return True


def trial_division(n):
    """Klassische Trial Division (für Vergleich)."""
    n = int(n)
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(n**0.5) + 1, 2):
        if n % i == 0:
            return False
    return True


# --- Benchmark-Hilfsfunktionen ---
BENCHMARK_N_RUNS = 100


def _benchmark(name, func, *args, n_runs=None):
    """Führt func(*args) n_runs mal aus und gibt mittlere Laufzeit in µs zurück."""
    n_runs = n_runs or BENCHMARK_N_RUNS
    times = []
    for _ in range(n_runs):
        t0 = time.perf_counter()
        func(*args)
        times.append((time.perf_counter() - t0) * 1e6)
    return np.mean(times), np.std(times)


def run_primality_benchmarks():
    """Vergleicht Miller-Rabin, Trial Division und Cache-Lookup."""
    from importlib.util import spec_from_file_location, module_from_spec
    spec = spec_from_file_location("qmbm", SCRIPT_DIR / "Second QMBM Modell.py")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    _load_primes_2m = mod._load_primes_2m
    _is_prime = mod._is_prime

    _load_primes_2m()

    # Testzahlen: klein (Cache), mittel (Trial), groß (Miller-Rabin Vorteil)
    test_cases = [
        (61, "klein (Cache)"),
        (32452843, "2M-te Primzahl"),
        (32452844, "2M-te Primzahl + 1 (zusammengesetzt)"),
        (10**9 + 7, "große Primzahl"),
        (10**9 + 9, "große zusammengesetzte"),
    ]

    n_runs = __import__("test_qmbm_performance").BENCHMARK_N_RUNS
    print("=" * 70)
    print(f"  PRIMZAHLTEST-PERFORMANCE (Mittelwert ± Std über {n_runs} Läufe, µs)")
    print("=" * 70)

    for n, desc in test_cases:
        mr_mean, mr_std = _benchmark("MR", miller_rabin, n)
        td_mean, td_std = _benchmark("TD", trial_division, n)
        ip_mean, ip_std = _benchmark("_is_prime", _is_prime, n)

        mr_ok = miller_rabin(n)
        td_ok = trial_division(n)
        ip_ok = _is_prime(n)
        ok = "OK" if (mr_ok == td_ok == ip_ok) else "DIFF"

        print(f"\n  n = {n} ({desc}) [{ok}]")
        print(f"    Miller-Rabin:   {mr_mean:8.2f} ± {mr_std:6.2f} µs")
        print(f"    Trial Division: {td_mean:8.2f} ± {td_std:6.2f} µs")
        print(f"    _is_prime:      {ip_mean:8.2f} ± {ip_std:6.2f} µs")

    # Batch: 10.000 is_prime-Checks (10K Qubit-Programm-Simulation)
    print("\n" + "-" * 70)
    print("  BATCH: 10.000 is_prime-Checks (simuliert 10K Qubit-Programm)")
    n_batch = 10_000
    t0 = time.perf_counter()
    for i in range(n_batch):
        _is_prime(i + 1)
    elapsed = (time.perf_counter() - t0) * 1e6
    print(f"    Gesamt: {elapsed:.0f} µs | pro Check: {elapsed/n_batch:.2f} µs")


def run_quantum_benchmarks():
    """Benchmark für Makro-Quantencomputer-Operationen."""
    from importlib.util import spec_from_file_location, module_from_spec
    spec = spec_from_file_location("qmbm", SCRIPT_DIR / "Second QMBM Modell.py")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    MacroQuantumComputer = mod.MacroQuantumComputer
    EnergiedokuWaveModel = mod.EnergiedokuWaveModel

    print("\n" + "=" * 70)
    print("  MAKRO-QUANTENCOMPUTER PERFORMANCE")
    print("=" * 70)

    mqc = MacroQuantumComputer(auto_load_persistent=True)
    if mqc.zeros is None:
        print("  [SKIP] Riemann-Nullstellen nicht verfügbar (zeros6.npz/npy)")
        return

    # resonance_check
    n_runs = 50
    t0 = time.perf_counter()
    for _ in range(n_runs):
        mqc.resonance_check(61, subset_size=100_000)
    elapsed = (time.perf_counter() - t0) * 1000
    print(f"\n  resonance_check(61): {elapsed/n_runs:.2f} ms pro Aufruf ({n_runs} Läufe)")

    # quantum_spectrum (512 Punkte)
    t_vals = np.linspace(0.01, 15.0, 512)
    t0 = time.perf_counter()
    amps = mqc.quantum_spectrum(t_vals)
    elapsed = (time.perf_counter() - t0) * 1000
    print(f"  quantum_spectrum(512 Punkte): {elapsed:.1f} ms")

    # run_10k_quantum_program
    model = EnergiedokuWaveModel(use_quantum_backend=True)
    model.load_from_quantum_computer()
    model.load_primes_2m()
    t0 = time.perf_counter()
    model.run_10k_quantum_program(num_qubits=10_000)
    elapsed = (time.perf_counter() - t0) * 1000
    print(f"  run_10k_quantum_program(10K): {elapsed:.0f} ms")


def run_correctness_tests():
    """Korrektheit: Miller-Rabin vs. 2M-Primzahlen-Cache."""
    from importlib.util import spec_from_file_location, module_from_spec
    spec = spec_from_file_location("qmbm", SCRIPT_DIR / "Second QMBM Modell.py")
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    _load_primes_2m = mod._load_primes_2m
    _is_prime = mod._is_prime

    primes = _load_primes_2m()
    if primes is None:
        print("[SKIP] Primzahlen nicht verfügbar")
        return

    # Stichprobe: erste 1000, letzte 1000, zufällige 1000
    sample = list(primes[:1000]) + list(primes[-1000:])
    rng = np.random.default_rng(42)
    idx = rng.choice(len(primes), size=min(1000, len(primes)), replace=False)
    sample.extend(primes[idx].tolist())

    errors = 0
    for p in sample:
        if not miller_rabin(p):
            errors += 1
            if errors <= 3:
                print(f"  FEHLER: Miller-Rabin sagt {p} nicht prim (Cache: prim)")

    # Zusammengesetzte: Nachbarn von Primzahlen
    for p in primes[:100]:
        c = p + 1
        if c % 2 != 0 and miller_rabin(c) and c not in set(primes.tolist()):
            errors += 1

    print("\n" + "=" * 70)
    print("  KORREKTHEIT: Miller-Rabin vs. 2M-Primzahlen-Cache")
    print("=" * 70)
    print(f"  Stichprobe: {len(sample)} Primzahlen")
    print(f"  Fehler: {errors}")
    print("  Status:", "PASS" if errors == 0 else "FAIL")


if __name__ == "__main__":
    import test_qmbm_performance as _m
    _m.BENCHMARK_N_RUNS = 20 if "--quick" in sys.argv else 100

    print("\n#Energiedoku QMBM Performance-Tests\n")
    run_correctness_tests()
    run_primality_benchmarks()
    run_quantum_benchmarks()
    print("\n" + "=" * 70)
    print("  Fertig.")
    print("=" * 70 + "\n")
