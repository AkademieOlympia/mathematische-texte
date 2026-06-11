import math
import random
import numpy as np


def classify(n):
    r = n % 12
    if r == 1:
        return "E"
    if r == 5:
        return "A"
    if r == 7:
        return "B"
    if r == 11:
        return "C"
    return "X"


def energy_bin(n):
    x = math.log(n)
    if x < 4:
        return "L1"
    if x < 6:
        return "L2"
    if x < 8:
        return "L3"
    if x < 10:
        return "L4"
    return "L5"


def is_prime(n):
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0:
        return False
    if n < 10_000:
        for i in range(3, int(n**0.5) + 1, 2):
            if n % i == 0:
                return False
        return True

    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2

    for a in (2, 3, 5, 7, 11):
        if a >= n:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True


def type_of(n):
    if is_prime(n):
        return "P"
    if any(n % i == 0 for i in range(2, 50)):
        return "M"
    return "S"


def gap_bin(prev, curr):
    g = abs(curr - prev)
    if g <= 2:
        return "g2"
    if g <= 4:
        return "g4"
    if g <= 10:
        return "g10"
    return "g_big"


def state_of(prev, n):
    return (
        classify(n),
        energy_bin(n),
        type_of(n),
        gap_bin(prev, n)
    )


def simulate(start=101, steps=20000):
    n = start
    prev = n
    states = []

    for _ in range(steps):
        states.append(state_of(prev, n))

        prev = n

        if random.random() < 0.4:
            n = 2 * n
        elif n % 2 == 0:
            n = n // 2 + 3
        else:
            n = n + random.randint(-5, 5)

        if n < 2:
            n = 2

    return states


def build_markov(states):
    unique = sorted(set(states))
    idx = {s: i for i, s in enumerate(unique)}

    P = np.zeros((len(unique), len(unique)))

    for i in range(len(states) - 1):
        a = states[i]
        b = states[i + 1]
        P[idx[a], idx[b]] += 1

    row_sums = P.sum(axis=1, keepdims=True)
    P = np.divide(P, row_sums, out=np.zeros_like(P), where=row_sums != 0)

    return P, unique


def spectrum(P):
    eigvals = np.linalg.eigvals(P)
    return sorted(eigvals, key=lambda x: -abs(x))


if __name__ == "__main__":
    random.seed(42)

    states = simulate(start=101, steps=20000)

    print("ANZAHL STATES:", len(states))
    print("UNIQUE STATES:", len(set(states)))

    P, unique = build_markov(states)
    eigvals = spectrum(P)

    print("\nstates", len(unique))
    print("top", eigvals[:5])

    modes = [(l, abs(l), np.angle(l)) for l in eigvals[:5]]
    print("modes", modes)