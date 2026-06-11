"""
Hardy–Littlewood / Bateman–Horn: lokale Faktoren fürs Zwillingsmuster (0,2) in Restklassen mod 60.

- Globale Zwillingskonstante: C2 = ∏_{p>2} (1 - 1/(p-1)²)  (twin prime constant).
- Unter der üblichen Bateman–Horn-Form: die sechs zulässigen Klassen
  {11, 17, 29, 41, 47, 59} teilen die gleiche asymptotische Dichte (Kante 1/6 je Kante).

Empirische Anpassung: Abweichung der beobachteten Klassenhäufigkeit von 1/6
(Singular-Series-Proxy auf endlichem N), Skalierung der Kanal-Deltas d_r.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

# Twin-Kanäle (kleinerer Term des Zwillings) mod 60, die für p>3 tatsächlich vorkommen.
TWIN_MOD60: Final[tuple[int, ...]] = (11, 17, 29, 41, 47, 59)
MOD: Final[int] = 60


@dataclass(frozen=True, slots=True)
class TwinC2:
    c2: float
    two_c2: float
    prime_upper: int
    n_primes: int


def _primes_to(n: int) -> list[int]:
    if n < 2:
        return []
    s = bytearray(b"\1") * (n + 1)
    s[0:2] = b"\0\x00"
    r = int(math.isqrt(n))
    for p in range(2, r + 1):
        if s[p]:
            step = p
            start = p * p
            s[start : n + 1 : step] = b"\0" * ((n - start) // step + 1)
    return [i for i in range(2, n + 1) if s[i]]


def compute_twin_prime_constant_C2(upper: int = 2_000_000) -> TwinC2:
    """
    C2 = ∏_{p>2, p≤upper} (1 - 1/(p-1)²); konvergiert schnell; upper groß wählen für 15+ Stellen.
    """
    primes = _primes_to(upper)
    prod = 1.0
    n_used = 0
    for p in primes:
        if p <= 2:
            continue
        prod *= 1.0 - 1.0 / (p - 1) ** 2
        n_used += 1
    return TwinC2(c2=prod, two_c2=2.0 * prod, prime_upper=upper, n_primes=n_used)


def bateman_horn_theoretical_pmf_mod60() -> dict[int, float]:
    """
    Asymptotische PMF P(Start ≡ a) für a ∈ TWIN_MOD60: uniform 1/6.
    """
    w = 1.0 / len(TWIN_MOD60)
    return {a: w for a in TWIN_MOD60}


def hl_correction_vs_uniform(count_by_channel: dict[int, int], total_twins: int) -> dict[int, float]:
    """
    Empirische HL-Abweichung: κ_a = f_a / (1/6) = 6 n_a / T, wobei f_a = n_a / T.
    """
    if total_twins <= 0:
        raise ValueError("total_twins muss > 0 sein")
    w_hl = 1.0 / 6.0
    kappa: dict[int, float] = {}
    for a in TWIN_MOD60:
        n_a = int(count_by_channel.get(a, 0))
        f_a = n_a / float(total_twins)
        kappa[a] = f_a / w_hl
    return kappa


def delta_normalized_resid_emp(
    d_r_11: float,
    d_r_47: float,
    kappa_11: float,
    kappa_47: float,
) -> float:
    """
    P_11,47^resid(emp) = d_47/κ_47 - d_11/κ_11
    (Heuristik: beobachtete Klassenanteile gegen HL-Uniform 1/6 rechnen).
    """
    if kappa_11 == 0.0 or kappa_47 == 0.0:
        return float("nan")
    return (d_r_47 / kappa_47) - (d_r_11 / kappa_11)
