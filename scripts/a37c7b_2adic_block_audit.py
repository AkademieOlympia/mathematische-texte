from math import log2
from collections import Counter

LIMIT = 100_000
MAX_STEPS = 500
MAX_BLOCK = 40


def v2(n: int) -> int:
    if n <= 0:
        raise ValueError("v2 requires n > 0")
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def U(n: int):
    if n <= 0 or n % 2 == 0:
        raise ValueError("U is defined here on positive odd integers")
    m = 3 * n + 1
    a = v2(m)
    return m >> a, a


def odd_trajectory(start: int, max_steps=MAX_STEPS):
    if start <= 0 or start % 2 == 0:
        raise ValueError("odd positive start required")

    n = start
    states = [n]
    vals = []

    for _ in range(max_steps):
        if n == 1:
            break
        n, a = U(n)
        states.append(n)
        vals.append(a)

    return states, vals


def first_descent_block(start: int, max_block=MAX_BLOCK):
    if start <= 1 or start % 2 == 0:
        return None

    n = start
    A = 0

    for r in range(1, max_block + 1):
        n, a = U(n)
        A += a

        E = A - r * log2(3.0)

        if n < start:
            return {
                "start": start,
                "end": n,
                "length": r,
                "valuation_sum": A,
                "excess": E,
                "ratio": n / start,
            }

    return None


def valuation_block(start: int, r: int):
    n = start
    vals = []
    A = 0

    for _ in range(r):
        n, a = U(n)
        vals.append(a)
        A += a

    return {
        "start": start,
        "end": n,
        "valuations": tuple(vals),
        "valuation_sum": A,
        "excess": A - r * log2(3.0),
        "ratio": n / start,
    }


def audit(limit=LIMIT):
    starts = range(3, limit + 1, 2)

    first_v2 = Counter()
    descent_lengths = Counter()
    unresolved = []
    a1_cases = 0
    a1_descents = []

    max_descent_length = 0
    max_descent_start = None

    for start in starts:
        _, a0 = U(start)
        first_v2[a0] += 1

        if a0 == 1:
            a1_cases += 1

        block = first_descent_block(start)

        if block is None:
            unresolved.append(start)
            continue

        r = block["length"]
        descent_lengths[r] += 1

        if r > max_descent_length:
            max_descent_length = r
            max_descent_start = start

        if a0 == 1:
            a1_descents.append(block)

    print("A3.7c7b — 2-Adic Block Audit")
    print(f"odd starts audited : {(limit - 1) // 2}")
    print(f"range              : 3 .. {limit}")
    print(f"max block          : {MAX_BLOCK}")
    print()
    print("Initial v2 distribution")
    for a in sorted(first_v2):
        print(f"v2={a:2d} : {first_v2[a]}")
    print()
    print(f"a0 = 1 cases       : {a1_cases}")
    print(f"unresolved <= {MAX_BLOCK}: {len(unresolved)}")
    print(f"max first-descent length : {max_descent_length}")
    print(f"witness start            : {max_descent_start}")

    if max_descent_start is not None:
        print()
        print("Longest first-descent witness")
        print(first_descent_block(max_descent_start))

    print()
    print("First-descent block lengths")
    for r in sorted(descent_lengths):
        print(f"r={r:2d} : {descent_lengths[r]}")

    if unresolved:
        print()
        print("First unresolved starts")
        print(unresolved[:50])

    print()
    print("Selected classical trajectories")
    for start in (27, 871, 6171):
        states, vals = odd_trajectory(start)
        print(
            f"{start:5d}: odd_steps={len(vals):3d} "
            f"sum_v2={sum(vals):4d} "
            f"reached_1={states[-1] == 1}"
        )


if __name__ == "__main__":
    audit()
