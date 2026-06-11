import time
from math import isqrt


def simple_sieve(limit: int):
    """Kleines Basissieb bis sqrt(X)."""
    if limit < 2:
        return []

    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start:limit + 1:p] = b"\x00" * (((limit - start) // p) + 1)

    return [p for p in range(2, limit + 1) if is_prime[p]]


def segmented_sieve(limit: int, segment_size: int = 1_000_000):
    """
    Segmentiertes Sieb bis limit.
    Rückgabe:
      - is_prime_global: bool-Feld für O(1)-Nachschlagen bis limit
      - primes: Primzahlliste bis limit
    """
    if limit < 2:
        return bytearray(limit + 1), []

    base_primes = simple_sieve(isqrt(limit))
    is_prime_global = bytearray(b"\x00") * (limit + 1)
    primes = []

    low = 2
    while low <= limit:
        high = min(low + segment_size - 1, limit)
        segment = bytearray(b"\x01") * (high - low + 1)

        for p in base_primes:
            start = max(p * p, ((low + p - 1) // p) * p)
            if start > high:
                continue
            segment[start - low: high - low + 1: p] = b"\x00" * (((high - start) // p) + 1)

        if low == 2:
            pass
        else:
            # gerade Zahlen > 2 sind sowieso nicht prim; optional könnte man odd-only optimieren
            pass

        for i, flag in enumerate(segment):
            n = low + i
            if flag and n >= 2:
                is_prime_global[n] = 1
                primes.append(n)

        low = high + 1

    return is_prime_global, primes


def familie(p: int) -> str:
    r = p % 12
    if r == 1:
        return "E"
    if r == 5:
        return "A"
    if r == 7:
        return "B"
    if r == 11:
        return "C"
    raise ValueError(f"Primzahl {p} ist nicht in EABC klassifizierbar.")


def build_family_lists(primes: list[int]):
    fam = {"E": [], "A": [], "B": [], "C": []}
    for p in primes:
        if p > 3:
            fam[familie(p)].append(p)
    return fam


def benchmark_direct(is_prime: bytearray, primes: list[int]):
    starts = [p for p in primes if p > 3]
    max_n = len(is_prime) - 1

    t0 = time.perf_counter()

    twin_candidates = 0
    drill_i_candidates = 0
    drill_ii_candidates = 0
    quad_candidates = 0

    twins = 0
    drill_i = 0
    drill_ii = 0
    quads = 0

    for p in starts:
        if p + 2 <= max_n:
            twin_candidates += 1
            if is_prime[p + 2]:
                twins += 1

        if p + 6 <= max_n:
            drill_i_candidates += 1
            if is_prime[p + 2] and is_prime[p + 6]:
                drill_i += 1

            drill_ii_candidates += 1
            if is_prime[p + 4] and is_prime[p + 6]:
                drill_ii += 1

        if p + 8 <= max_n:
            quad_candidates += 1
            if is_prime[p + 2] and is_prime[p + 6] and is_prime[p + 8]:
                quads += 1

    t1 = time.perf_counter()

    return {
        "search_time": t1 - t0,
        "twin_candidates": twin_candidates,
        "drill_i_candidates": drill_i_candidates,
        "drill_ii_candidates": drill_ii_candidates,
        "quad_candidates": quad_candidates,
        "twins": twins,
        "drill_i": drill_i,
        "drill_ii": drill_ii,
        "quads": quads,
    }


def benchmark_eabc_filtered(is_prime: bytearray, fam: dict[str, list[int]]):
    starts_ac = fam["A"] + fam["C"]
    starts_eb = fam["E"] + fam["B"]
    max_n = len(is_prime) - 1

    t0 = time.perf_counter()

    twin_candidates = 0
    drill_i_candidates = 0
    drill_ii_candidates = 0
    quad_candidates = 0

    twins = 0
    drill_i = 0
    drill_ii = 0
    quads = 0

    for p in starts_ac:
        if p + 2 <= max_n:
            twin_candidates += 1
            if is_prime[p + 2]:
                twins += 1

        if p + 6 <= max_n:
            drill_i_candidates += 1
            if is_prime[p + 2] and is_prime[p + 6]:
                drill_i += 1

        if p + 8 <= max_n:
            quad_candidates += 1
            if is_prime[p + 2] and is_prime[p + 6] and is_prime[p + 8]:
                quads += 1

    for p in starts_eb:
        if p + 6 <= max_n:
            drill_ii_candidates += 1
            if is_prime[p + 4] and is_prime[p + 6]:
                drill_ii += 1

    t1 = time.perf_counter()

    return {
        "search_time": t1 - t0,
        "twin_candidates": twin_candidates,
        "drill_i_candidates": drill_i_candidates,
        "drill_ii_candidates": drill_ii_candidates,
        "quad_candidates": quad_candidates,
        "twins": twins,
        "drill_i": drill_i,
        "drill_ii": drill_ii,
        "quads": quads,
    }


def pct(hit: int, cand: int) -> float:
    return 100.0 * hit / cand if cand else 0.0


def print_report(limit: int, segment_size: int = 1_000_000):
    t0 = time.perf_counter()
    is_prime, primes = segmented_sieve(limit, segment_size=segment_size)
    t1 = time.perf_counter()

    t2 = time.perf_counter()
    fam = build_family_lists(primes)
    t3 = time.perf_counter()

    direct = benchmark_direct(is_prime, primes)
    filt = benchmark_eabc_filtered(is_prime, fam)
    t4 = time.perf_counter()

    print("=" * 78)
    print(f"Segmentierter EABC-Benchmark bis {limit:,}".replace(",", "."))
    print("=" * 78)
    print(f"Anzahl Primzahlen:      {len(primes)}")
    print(f"Segmentiertes Sieb:     {t1 - t0:.6f} s")
    print(f"Familienaufbau:         {t3 - t2:.6f} s")
    print(f"Direkte Suchphase:      {direct['search_time']:.6f} s")
    print(f"EABC-Suchphase:         {filt['search_time']:.6f} s")
    print(f"Gesamt inkl. Aufbau:    {t4 - t0:.6f} s")
    print()

    print("Familiengrößen (ohne 2,3):")
    for k in ["E", "A", "B", "C"]:
        print(f"  {k}: {len(fam[k])}")
    print()

    rows = [
        ("reguläre Primzwillinge",
         direct["twin_candidates"], direct["twins"],
         filt["twin_candidates"], filt["twins"]),
        ("primitive Drillinge Typ I",
         direct["drill_i_candidates"], direct["drill_i"],
         filt["drill_i_candidates"], filt["drill_i"]),
        ("primitive Drillinge Typ II",
         direct["drill_ii_candidates"], direct["drill_ii"],
         filt["drill_ii_candidates"], filt["drill_ii"]),
        ("echte Primvierlinge",
         direct["quad_candidates"], direct["quads"],
         filt["quad_candidates"], filt["quads"]),
    ]

    print("Kandidaten- und Treffervergleich")
    print("-" * 78)
    for name, dc, dh, fc, fh in rows:
        print(name)
        print(f"  direkt: Kandidaten={dc}, Treffer={dh}, Quote={pct(dh, dc):.4f}%")
        print(f"  EABC:   Kandidaten={fc}, Treffer={fh}, Quote={pct(fh, fc):.4f}%")
        print()

    print("Geschwindigkeitsfaktoren der reinen Suchphase")
    print("-" * 78)
    if filt["search_time"] > 0:
        print(f"direkt / EABC = {direct['search_time'] / filt['search_time']:.3f}")
    else:
        print("direkt / EABC = n/a")


if __name__ == "__main__":
    # Beispiel:
    #print_report(5_000_000, segment_size=500_000)

    # Für größere Grenzen:
    print_report(10_000_000, segment_size=1_000_000)
    print_report(50_000_000, segment_size=2_000_000)