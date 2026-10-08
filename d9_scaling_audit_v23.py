import math
import time
import random
import sys
import csv
import os
from typing import List, Dict, Any, Tuple
import numpy as np
import sympy

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from d9_quantum_master_app_v23 import D9QuantumMasterApp

# Global App Instance
_GLOBAL_D9_APP = D9QuantumMasterApp()

def factor_v23(N: int) -> Dict[str, Any]:
    """
    Direct Python API Call for v23.0 Factorization Routine.
    """
    return _GLOBAL_D9_APP.factorize_full_recursive(N)

def generate_balanced_prime(bits: int) -> int:
    while True:
        candidate = random.getrandbits(bits) | (1 << (bits - 1)) | 1
        if sympy.isprime(candidate):
            return candidate

def generate_semiprime_pair(total_bits: int, ratio_class: str = "A") -> Tuple[int, int, int, float]:
    half_bits = total_bits // 2
    while True:
        p = generate_balanced_prime(half_bits)
        if ratio_class == "A":
            min_q = int(p * 1.00)
            max_q = int(p * 1.05)
        else:
            min_q = int(p * 1.30)
            max_q = int(p * 1.50)
            
        q_cand = sympy.nextprime(min_q)
        if q_cand <= max_q:
            N = p * q_cand
            if N.bit_length() == total_bits:
                ratio = q_cand / float(p)
                return N, min(p, q_cand), max(p, q_cand), ratio

if __name__ == "__main__":
    print("Testing factor_v23 wrapper...")
    N, p, q, ratio = generate_semiprime_pair(64, "A")
    res = factor_v23(N)
    print(f"64-bit Test: N={N}, Factors={res['prime_factors']}, Time={res['time_ms']:.2f}ms")
