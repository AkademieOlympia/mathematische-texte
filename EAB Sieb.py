import time
from math import isqrt


def odd_index(n: int) -> int:
    """Index der ungeraden Zahl n >= 3 im odd-only-Feld."""
    return (n - 3) // 2


def simple_sieve_odd(limit: int):
    """
    Odd-only Basissieb bis limit.
    Rückgabe: Liste aller Primzahlen <= limit.
    """
    if limit < 2:
        return []
    if limit == 2:
        return [2]

    size = odd_index(limit) + 1
    is_prime_odd = bytearray(b"\x01") * size

    max_p = isqrt(limit)
    for i in range(odd_index(3), odd_index(max_p) + 1):
        if is_prime_odd[i]:
            p = 2 * i + 3
            start = p * p
            step = 2 * p
            for n in range(start, limit + 1, step):
                is_prime_odd[odd_index(n)] = 0

    primes = [2]
    primes.extend(2 * i + 3 for i, flag in enumerate(is_prime_odd) if flag)
    return primes


class OddPrimeTable:
    """
    Odd-only Primtabelle für O(1)-Nachschlagen:
    - 2 wird separat behandelt
    - ungerade Zahlen >= 3 liegen im Bitfeld
    """
    def __init__(self, limit: int):
        self.limit = limit
        size = odd_index(limit) + 1 if limit >= 3 else 0
        self.bits = bytearray(b"\x00") * size

    def set_prime(self, n: int):
        if n == 2:
            return
        if n >= 3 and n % 2 == 1:
            self.bits[odd_index(n)] = 1

    def is_prime(self, n: int) -> bool:
        if n == 2:
            return True
        if n < 2 or n > self.limit or n % 2 == 0:
            return False
        return bool(self.bits[odd_index(n)])


def segmented_sieve_odd(limit: int, segment_size: int = 1_000_000):
    """
    Odd-only segmentiertes Sieb bis limit.
    Rückgabe:
      - OddPrimeTable für O(1)-Abfragen
      - Primzahlliste bis limit
    """
    if limit < 2:
        return OddPrimeTable(limit), []

    base_primes = simple_sieve_odd(isqrt(limit))
    table = OddPrimeTable(limit)
    primes = [2] if limit >= 2 else []

    # nur ungerade Segmente
    low = 3
    if low % 2 == 0:
        low += 1

    while low <= limit:
        high = min(low + segment_size - 1, limit)
        if high % 2 == 0:
            high -= 1
        if high < low:
            break

        seg_len = ((high - low) // 2) + 1
        segment = bytearray(b"\x01") * seg_len

        for p in base_primes:
            if p == 2:
                continue
            p2 = p * p
            if p2 > high:
                break

            start = max(p2, ((low + p - 1) // p) * p)
            if start % 2 == 0:
                start += p

            step = 2 * p
            for n in range(start, high + 1, step):
                segment[(n - low) // 2] = 0

        for i, flag in enumerate(segment):
            if flag:
                n = low + 2 * i
                table.set_prime(n)
                primes.append(n)

        low = high + 2

    return table, primes


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


def benchmark_direct(table: OddPrimeTable, primes: list[int]):
    """
    Direkte Suche auf allen Primstarts p > 3.
    """
    starts = [p for p in primes if p > 3]
    max_n = table.limit

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
            if table.is_prime(p + 2):
                twins += 1

        if p + 6 <= max_n:
            drill_i_candidates += 1
            if table.is_prime(p + 2) and table.is_prime(p + 6):
                drill_i += 1

            drill_ii_candidates += 1
            if table.is_prime(p + 4) and table.is_prime(p + 6):
                drill_ii += 1

        if p + 8 <= max_n:
            quad_candidates += 1
            if table.is_prime(p + 2) and table.is_prime(p + 6) and table.is_prime(p + 8):
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


def benchmark_eabc_filtered(table: OddPrimeTable, fam: dict[str, list[int]]):
    """
    EABC-gefilterte Suche:
    - Zwillinge und Vierlinge nur von A,C
    - Drillinge I nur von A,C
    - Drillinge II nur von E,B
    """
    starts_ac = fam["A"] + fam["C"]
    starts_eb = fam["E"] + fam["B"]
    max_n = table.limit

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
            if table.is_prime(p + 2):
                twins += 1

        if p + 6 <= max_n:
            drill_i_candidates += 1
            if table.is_prime(p + 2) and table.is_prime(p + 6):
                drill_i += 1

        if p + 8 <= max_n:
            quad_candidates += 1
            if table.is_prime(p + 2) and table.is_prime(p + 6) and table.is_prime(p + 8):
                quads += 1

    for p in starts_eb:
        if p + 6 <= max_n:
            drill_ii_candidates += 1
            if table.is_prime(p + 4) and table.is_prime(p + 6):
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
    table, primes = segmented_sieve_odd(limit, segment_size=segment_size)
    t1 = time.perf_counter()

    t2 = time.perf_counter()
    fam = build_family_lists(primes)
    t3 = time.perf_counter()

    direct = benchmark_direct(table, primes)
    filt = benchmark_eabc_filtered(table, fam)
    t4 = time.perf_counter()

    print("=" * 82)
    print(f"Odd-only segmentierter EABC-Benchmark bis {limit:,}".replace(",", "."))
    print("=" * 82)
    print(f"Anzahl Primzahlen:      {len(primes)}")
    print(f"Odd-only Segmentsieb:   {t1 - t0:.6f} s")
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
    print("-" * 82)
    for name, dc, dh, fc, fh in rows:
        print(name)
        print(f"  direkt: Kandidaten={dc}, Treffer={dh}, Quote={pct(dh, dc):.4f}%")
        print(f"  EABC:   Kandidaten={fc}, Treffer={fh}, Quote={pct(fh, fc):.4f}%")
        print()

    print("Geschwindigkeitsfaktor der reinen Suchphase")
    print("-" * 82)
    if filt["search_time"] > 0:
        print(f"direkt / EABC = {direct['search_time'] / filt['search_time']:.3f}")
    else:
        print("direkt / EABC = n/a")


if __name__ == "__main__":
    #print_report(5_000_000, segment_size=500_000)

    # Weitere sinnvolle Tests:
    # print_report(10_000_000, segment_size=1_000_000)
    print_report(50_000_000, segment_size=2_000_000)
    print_report(100_000_000, segment_size=5_000_000)