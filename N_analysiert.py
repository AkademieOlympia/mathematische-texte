import math
from collections import Counter

# ------------------------------------------------------------
# Hilfsfunktionen
# ------------------------------------------------------------

def factorint(n: int) -> Counter:
    """Einfache Primfaktorzerlegung per Trial Division."""
    if n < 1:
        raise ValueError("n muss >= 1 sein")
    factors = Counter()
    x = n

    while x % 2 == 0:
        factors[2] += 1
        x //= 2

    p = 3
    while p * p <= x:
        while x % p == 0:
            factors[p] += 1
            x //= p
        p += 2

    if x > 1:
        factors[x] += 1

    return factors


def product_from_factors(factors: Counter) -> int:
    out = 1
    for p, e in factors.items():
        out *= p ** e
    return out


def smooth_part(n: int, basis: set[int]) -> int:
    fac = factorint(n)
    out = 1
    for p in basis:
        out *= p ** fac.get(p, 0)
    return out


def kernel_part(n: int, basis: set[int]) -> int:
    return n // smooth_part(n, basis)


def small_part(n: int, lam: int) -> int:
    fac = factorint(n)
    out = 1
    for p, e in fac.items():
        if p <= lam:
            out *= p ** e
    return out


def family_label_mod12(p: int) -> str | None:
    """
    Für Primzahlen > 3:
      E: 1 mod 12
      A: 5 mod 12
      B: 7 mod 12
      C: 11 mod 12
    Für 2,3 geben wir None zurück.
    """
    if p in (2, 3):
        return None
    r = p % 12
    if r == 1:
        return "E"
    if r == 5:
        return "A"
    if r == 7:
        return "B"
    if r == 11:
        return "C"
    raise ValueError(f"Unerwartete Restklasse für Primzahl {p}")


def family_parts_of_kernel(m: int) -> dict[str, int]:
    """Zerlegt den Restkern in E/A/B/C-Produkte."""
    fac = factorint(m)
    out = {"E": 1, "A": 1, "B": 1, "C": 1}
    for p, e in fac.items():
        fam = family_label_mod12(p)
        if fam is not None:
            out[fam] *= p ** e
    return out


def log_safe(x: int | float) -> float:
    if x <= 0:
        raise ValueError("log nur für x > 0")
    return math.log(x)


def family_weights(m: int) -> dict[str, float]:
    """Logarithmische Familiengewichte rho_X."""
    if m <= 1:
        return {"E": 0.0, "A": 0.0, "B": 0.0, "C": 0.0}
    fam_parts = family_parts_of_kernel(m)
    lm = log_safe(m)
    out = {}
    for fam in ("E", "A", "B", "C"):
        val = fam_parts[fam]
        out[fam] = 0.0 if val == 1 else log_safe(val) / lm
    return out


def family_entropy(weights: dict[str, float]) -> float:
    """Shannon-artige Entropie -sum rho log rho."""
    s = 0.0
    for v in weights.values():
        if v > 0:
            s -= v * math.log(v)
    return s


def epsilon_B(n: int, basis: set[int]) -> float:
    """Energetischer Quotient epsilon_B."""
    if n <= 1:
        return 0.0
    m = kernel_part(n, basis)
    return 0.0 if m == 1 else log_safe(m) / log_safe(n)


def curvature_K(n: int, basis: set[int]) -> float:
    """K_B(n) = epsilon_B(n) * C_fam(n)."""
    m = kernel_part(n, basis)
    if m == 1:
        return 0.0
    weights = family_weights(m)
    return epsilon_B(n, basis) * family_entropy(weights)


def family_vector(n: int, basis: set[int]) -> tuple[float, float, float, float]:
    """Gerichteter Familienvektor (rho_E, rho_A, rho_B, rho_C)."""
    m = kernel_part(n, basis)
    w = family_weights(m)
    return (w["E"], w["A"], w["B"], w["C"])


def euclidean_distance(v1, v2) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(v1, v2)))


def structural_distance(n: int, m: int, basis: set[int], lam: float = 1.0) -> float:
    """d_str^B(n,m)."""
    e1 = epsilon_B(n, basis)
    e2 = epsilon_B(m, basis)
    v1 = family_vector(n, basis)
    v2 = family_vector(m, basis)
    return math.sqrt((e1 - e2) ** 2 + lam * euclidean_distance(v1, v2) ** 2)


def format_family_support(m: int) -> str:
    """Gibt z.B. 'A·B' oder 'E·A·B' zurück."""
    fam_parts = family_parts_of_kernel(m)
    labels = [fam for fam in ("E", "A", "B", "C") if fam_parts[fam] > 1]
    return "1" if not labels else "·".join(labels)


# ------------------------------------------------------------
# Strukturdatensatz für eine Zahl
# ------------------------------------------------------------

def structure_record(n: int, basis: set[int], lam_small: int = 13) -> dict:
    fac = factorint(n)
    sB = smooth_part(n, basis)
    mB = n // sB
    kL = small_part(n, lam_small)

    fam_parts = family_parts_of_kernel(mB)
    weights = family_weights(mB)
    eps = epsilon_B(n, basis)
    c_fam = family_entropy(weights)
    K = eps * c_fam

    return {
        "n": n,
        "factorization": dict(fac),
        "sB": sB,
        "mB": mB,
        "kL": kL,
        "fam_parts": fam_parts,
        "fam_support": format_family_support(mB),
        "rho_E": weights["E"],
        "rho_A": weights["A"],
        "rho_B": weights["B"],
        "rho_C": weights["C"],
        "epsilon": eps,
        "C_fam": c_fam,
        "K": K,
        "family_vector": (weights["E"], weights["A"], weights["B"], weights["C"]),
    }


# ------------------------------------------------------------
# Ausgabehilfen
# ------------------------------------------------------------

def factorization_to_string(fac: dict[int, int]) -> str:
    parts = []
    for p in sorted(fac):
        e = fac[p]
        parts.append(f"{p}" if e == 1 else f"{p}^{e}")
    return " * ".join(parts)


def print_records(records: list[dict], max_rows: int | None = None) -> None:
    rows = records if max_rows is None else records[:max_rows]

    header = (
        f"{'n':>6} | {'Faktorisierung':<20} | {'sB':>6} | {'mB':>6} | "
        f"{'Fam':<8} | {'eps':>7} | {'C_fam':>7} | {'K':>7} | "
        f"{'rho=(E,A,B,C)':<30}"
    )
    print(header)
    print("-" * len(header))

    for r in rows:
        rho_str = (
            f"({r['rho_E']:.3f}, {r['rho_A']:.3f}, "
            f"{r['rho_B']:.3f}, {r['rho_C']:.3f})"
        )
        print(
            f"{r['n']:6d} | "
            f"{factorization_to_string(r['factorization']):<20} | "
            f"{r['sB']:6d} | {r['mB']:6d} | "
            f"{r['fam_support']:<8} | "
            f"{r['epsilon']:7.3f} | {r['C_fam']:7.3f} | {r['K']:7.3f} | "
            f"{rho_str:<30}"
        )


# ------------------------------------------------------------
# Exploration bis N
# ------------------------------------------------------------

def records_up_to(N: int, basis: set[int], lam_small: int = 13, start: int = 2) -> list[dict]:
    return [structure_record(n, basis, lam_small=lam_small) for n in range(start, N + 1)]


def top_by(records: list[dict], key: str, k: int = 20, reverse: bool = True) -> list[dict]:
    return sorted(records, key=lambda r: r[key], reverse=reverse)[:k]


def nearest_neighbors(target_n: int, records: list[dict], basis: set[int], lam: float = 1.0, k: int = 10):
    out = []
    for r in records:
        n = r["n"]
        if n == target_n:
            continue
        d = structural_distance(target_n, n, basis=basis, lam=lam)
        out.append((d, n))
    out.sort()
    return out[:k]


# ------------------------------------------------------------
# Demo
# ------------------------------------------------------------

if __name__ == "__main__":
    basis = {2, 3}

    demo_numbers = [360, 420, 441, 455, 540, 693, 945]
    demo_records = [structure_record(n, basis) for n in demo_numbers]

    print("\nBeispielzahlen:\n")
    print_records(demo_records)

    print("\nTop nach Strukturkrümmung K bis 500:\n")
    recs_500 = records_up_to(500, basis)
    print_records(top_by(recs_500, "K", k=15))

    print("\nTop nach Familienentropie C_fam bis 500:\n")
    print_records(top_by(recs_500, "C_fam", k=15))

    print("\nTop nach energetischem Quotienten epsilon bis 500:\n")
    print_records(top_by(recs_500, "epsilon", k=15))

    target = 420
    print(f"\nNächste strukturelle Nachbarn von {target} bis 500:\n")
    neigh = nearest_neighbors(target, recs_500, basis, lam=1.0, k=12)
    for d, n in neigh:
        print(f"n={n:4d}, d_str={d:.4f}")