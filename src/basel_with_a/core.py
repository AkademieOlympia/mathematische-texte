"""Core number-theoretic routines for basel-with-a."""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt, pi, sqrt


@dataclass(frozen=True, slots=True)
class BaselSummary:
    prime_count: int
    nth_prime: int
    prime_square_sum: int
    zeta_2: float
    difference: float
    index_value_a: float


def _is_prime(candidate: int) -> bool:
    if candidate < 2:
        return False
    if candidate == 2:
        return True
    if candidate % 2 == 0:
        return False

    limit = isqrt(candidate)
    divisor = 3
    while divisor <= limit:
        if candidate % divisor == 0:
            return False
        divisor += 2
    return True


def first_n_primes(count: int) -> list[int]:
    if count < 1:
        raise ValueError("count muss mindestens 1 sein")

    primes: list[int] = []
    candidate = 2
    while len(primes) < count:
        if _is_prime(candidate):
            primes.append(candidate)
        candidate = 3 if candidate == 2 else candidate + 2
    return primes


def nth_prime(count: int) -> int:
    return first_n_primes(count)[-1]


def analyze_first_primes(count: int = 100) -> BaselSummary:
    primes = first_n_primes(count)
    prime_square_sum = sum(prime * prime for prime in primes)
    zeta_2 = (pi**2) / 6.0

    return BaselSummary(
        prime_count=count,
        nth_prime=primes[-1],
        prime_square_sum=prime_square_sum,
        zeta_2=zeta_2,
        difference=prime_square_sum - zeta_2,
        index_value_a=sqrt(prime_square_sum),
    )
