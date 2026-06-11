"""
Asymmetrie-Skala für Primzahlphasen (oktonisch motivierter Bias, ohne SageMath).
"""

from __future__ import annotations

import math
import sys
from typing import List


def primes_range_lt(n: int) -> List[int]:
    """
    Primzahlen p mit p < n (analog zu Sage ``primes_range(n)``).
    """
    if n <= 2:
        return []
    is_prime = bytearray(b"\x01") * n
    is_prime[0:2] = b"\x00\x00"
    for i in range(2, int(math.isqrt(n - 1)) + 1):
        if is_prime[i]:
            for j in range(i * i, n, i):
                is_prime[j] = 0
    return [i for i in range(2, n) if is_prime[i]]


def riemann_octonionic_bias(n_range: int = 1000, bias: float = 0.01) -> float:
    """
    Mittelwert von cos(phase) mit phase = (p * bias) mod 2π über alle Primen p < n_range.
    """
    primes = primes_range_lt(n_range)
    if not primes:
        return float("nan")
    s = 0.0
    two_pi = 2 * math.pi
    for p in primes:
        phase = (p * bias) % two_pi
        s += math.cos(phase)
    return s / len(primes)


def main() -> None:
    n = 50_000
    if len(sys.argv) > 1:
        try:
            n = int(sys.argv[1])
        except ValueError:
            print("Usage: python3 Nullstellenverschiebung.py [n_range]", file=sys.stderr)
            sys.exit(1)
    v = riemann_octonionic_bias(n, bias=0.01)
    print(f"Symmetrie-Abweichung der Primzahl-Phasen (Mittel cos(phase), p < {n}): {v}")


# Direkt- und Notebook-/runpy-Start
if __name__ in ("__main__", "<run_path>"):
    main()
