import math
import csv
from collections import Counter, defaultdict
from pathlib import Path

# ============================================================
# Grundeinstellungen
# ============================================================

BASIS = {2, 3}
LAMBDA_SMALL = 13
DEFAULT_MAX_N = 1000

# Gewichte für die verfeinerte Strukturdistanz
W_EPS = 1.0
W_KAPPA = 1.0
W_RHO = 1.0

# Ausgabeordner
OUTPUT_DIR = Path(".")


# ============================================================
# Arithmetische Grundfunktionen
# ============================================================

def factorint(n: int) -> Counter:
    """Primfaktorzerlegung von n >= 1 per Trial Division."""
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


def factorization_to_string(fac: dict[int, int]) -> str:
    if not fac:
        return "1"
    parts = []
    for p in sorted(fac):
        e = fac[p]
        parts.append(f"{p}" if e == 1 else f"{p}^{e}")
    return " * ".join(parts)


def log_safe(x: int | float) -> float:
    if x <= 0:
        raise ValueError("log ist nur für x > 0 definiert")
    return math.log(x)


# ============================================================
# Glatter Anteil, Restkern, kleiner Anteil
# ============================================================

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


def sigma_B(n: int, basis: set[int]) -> float:
    """Anteil des glatten Trägers relativ zur Gesamtgröße."""
    if n <= 1:
        return 0.0
    sB = smooth_part(n, basis)
    return log_safe(sB) / log_safe(n)


def epsilon_B(n: int, basis: set[int]) -> float:
    """Anteil des Restkerns relativ zur Gesamtgröße."""
    if n <= 1:
        return 0.0
    mB = kernel_part(n, basis)
    return 0.0 if mB == 1 else log_safe(mB) / log_safe(n)


def kappa_lambda(n: int, lam: int) -> float:
    """Anteil kleiner Primträger relativ zur Gesamtgröße."""
    if n <= 1:
        return 0.0
    kL = small_part(n, lam)
    return log_safe(kL) / log_safe(n)


# ============================================================
# Familien modulo 12
# ============================================================

def family_label_mod12(p: int) -> str | None:
    """
    Für Primzahlen > 3:
      E: 1 mod 12
      A: 5 mod 12
      B: 7 mod 12
      C: 11 mod 12
    Für 2, 3 -> None
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
    raise ValueError(f"Unerwartete Restklasse {r} für Primzahl {p}")


def family_parts_of_kernel(m: int) -> dict[str, int]:
    """
    Zerlegt einen {2,3}-freien Restkern in E/A/B/C-Produkte.

    Die Familienklassifikation modulo 12 ist hier nur für Primfaktoren > 3
    sinnvoll. Falls 2 oder 3 noch im Restkern vorkommen, wäre die spätere
    Normierung der Familiengewichte mathematisch inkonsistent. Deshalb brechen
    wir in diesem Fall bewusst ab, statt stillschweigend falsche rho-Werte zu
    erzeugen.
    """
    fac = factorint(m)
    if 2 in fac or 3 in fac:
        raise ValueError(
            "family_parts_of_kernel erwartet einen {2,3}-freien Restkern; "
            "2 oder 3 sind noch in m enthalten."
        )
    out = {"E": 1, "A": 1, "B": 1, "C": 1}
    for p, e in fac.items():
        fam = family_label_mod12(p)
        if fam is not None:
            out[fam] *= p ** e
    return out


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
    s = 0.0
    for v in weights.values():
        if v > 0:
            s -= v * math.log(v)
    return s


def family_support_string(m: int) -> str:
    fam_parts = family_parts_of_kernel(m)
    labels = [fam for fam in ("E", "A", "B", "C") if fam_parts[fam] > 1]
    return "1" if not labels else "·".join(labels)


def family_vector(n: int, basis: set[int]) -> tuple[float, float, float, float]:
    mB = kernel_part(n, basis)
    w = family_weights(mB)
    return (w["E"], w["A"], w["B"], w["C"])


# ============================================================
# Strukturgrößen
# ============================================================

def curvature_K(n: int, basis: set[int]) -> float:
    """Bisherige Strukturkrümmung K = epsilon * C_fam."""
    mB = kernel_part(n, basis)
    if mB == 1:
        return 0.0
    w = family_weights(mB)
    return epsilon_B(n, basis) * family_entropy(w)


def curvature_K2(n: int, basis: set[int]) -> float:
    """
    Alternative Strukturkrümmung:
    absolute Restkernmasse mal Entropie, normiert über log(1+sB)
    """
    if n <= 1:
        return 0.0
    sB = smooth_part(n, basis)
    mB = kernel_part(n, basis)
    if mB == 1:
        return 0.0
    c_fam = family_entropy(family_weights(mB))
    return (log_safe(mB) / log_safe(1 + sB)) * c_fam


def curvature_K3(n: int, basis: set[int], alpha: float = 1.0) -> float:
    """
    Alternative Strukturkrümmung:
    max(log mB - alpha log sB, 0) * C_fam
    """
    if n <= 1:
        return 0.0
    sB = smooth_part(n, basis)
    mB = kernel_part(n, basis)
    if mB == 1:
        return 0.0
    c_fam = family_entropy(family_weights(mB))
    val = max(log_safe(mB) - alpha * log_safe(max(sB, 1)), 0.0)
    return val * c_fam


# ============================================================
# Strukturdistanzen
# ============================================================

def euclidean_distance(v1, v2) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(v1, v2)))


def family_distance(n: int, m: int, basis: set[int]) -> float:
    return euclidean_distance(family_vector(n, basis), family_vector(m, basis))


def structural_distance_v1(
    n: int,
    m: int,
    basis: set[int],
    lambda_rho: float = 1.0
) -> float:
    """
    Alte Distanz:
    sqrt((eps_n-eps_m)^2 + lambda * ||rho_n-rho_m||^2)
    """
    e1 = epsilon_B(n, basis)
    e2 = epsilon_B(m, basis)
    r1 = family_vector(n, basis)
    r2 = family_vector(m, basis)
    return math.sqrt((e1 - e2) ** 2 + lambda_rho * euclidean_distance(r1, r2) ** 2)


def structural_distance_v2(
    n: int,
    m: int,
    basis: set[int],
    lam_small: int,
    w_eps: float = 1.0,
    w_kappa: float = 1.0,
    w_rho: float = 1.0
) -> float:
    """
    Neue Distanz:
    sqrt(w_eps*(eps diff)^2 + w_kappa*(kappa diff)^2 + w_rho*||rho diff||^2)
    """
    e1 = epsilon_B(n, basis)
    e2 = epsilon_B(m, basis)

    k1 = kappa_lambda(n, lam_small)
    k2 = kappa_lambda(m, lam_small)

    r1 = family_vector(n, basis)
    r2 = family_vector(m, basis)

    return math.sqrt(
        w_eps * (e1 - e2) ** 2
        + w_kappa * (k1 - k2) ** 2
        + w_rho * euclidean_distance(r1, r2) ** 2
    )


# ============================================================
# Strukturdatensatz
# ============================================================

def structure_record(n: int, basis: set[int], lam_small: int) -> dict:
    fac = factorint(n)
    sB = smooth_part(n, basis)
    mB = kernel_part(n, basis)
    kL = small_part(n, lam_small)

    fam_parts = family_parts_of_kernel(mB)
    rho = family_weights(mB)

    eps = epsilon_B(n, basis)
    sig = sigma_B(n, basis)
    kap = kappa_lambda(n, lam_small)
    c_fam = family_entropy(rho)
    K = curvature_K(n, basis)
    K2 = curvature_K2(n, basis)
    K3 = curvature_K3(n, basis, alpha=1.0)

    return {
        "n": n,
        "factorization": dict(fac),
        "factorization_str": factorization_to_string(dict(fac)),
        "sB": sB,
        "mB": mB,
        "kL": kL,
        "sigma": sig,
        "epsilon": eps,
        "kappa": kap,
        "fam_E_part": fam_parts["E"],
        "fam_A_part": fam_parts["A"],
        "fam_B_part": fam_parts["B"],
        "fam_C_part": fam_parts["C"],
        "fam_support": family_support_string(mB),
        "rho_E": rho["E"],
        "rho_A": rho["A"],
        "rho_B": rho["B"],
        "rho_C": rho["C"],
        "C_fam": c_fam,
        "K": K,
        "K2": K2,
        "K3": K3,
        "family_vector": (rho["E"], rho["A"], rho["B"], rho["C"]),
    }


def records_up_to(N: int, basis: set[int], lam_small: int, start: int = 2) -> list[dict]:
    return [structure_record(n, basis, lam_small) for n in range(start, N + 1)]


# ============================================================
# Sortierung / Anzeige
# ============================================================

def top_by(records: list[dict], key: str, k: int = 20, reverse: bool = True) -> list[dict]:
    return sorted(records, key=lambda r: r[key], reverse=reverse)[:k]


def print_records(records: list[dict], max_rows: int | None = None) -> None:
    rows = records if max_rows is None else records[:max_rows]

    header = (
        f"{'n':>6} | {'Faktorisierung':<22} | {'sB':>6} | {'mB':>6} | {'kL':>6} | "
        f"{'Fam':<8} | {'sig':>6} | {'eps':>6} | {'kap':>6} | "
        f"{'C_fam':>7} | {'K':>7} | {'K2':>7}"
    )
    print(header)
    print("-" * len(header))

    for r in rows:
        print(
            f"{r['n']:6d} | "
            f"{r['factorization_str']:<22} | "
            f"{r['sB']:6d} | {r['mB']:6d} | {r['kL']:6d} | "
            f"{r['fam_support']:<8} | "
            f"{r['sigma']:6.3f} | {r['epsilon']:6.3f} | {r['kappa']:6.3f} | "
            f"{r['C_fam']:7.3f} | {r['K']:7.3f} | {r['K2']:7.3f}"
        )


# ============================================================
# Nachbarschaften
# ============================================================

def nearest_neighbors_v2(
    target_n: int,
    records: list[dict],
    basis: set[int],
    lam_small: int,
    k: int = 10,
    w_eps: float = 1.0,
    w_kappa: float = 1.0,
    w_rho: float = 1.0
):
    out = []
    for r in records:
        n = r["n"]
        if n == target_n:
            continue
        d = structural_distance_v2(
            target_n, n,
            basis=basis,
            lam_small=lam_small,
            w_eps=w_eps,
            w_kappa=w_kappa,
            w_rho=w_rho
        )
        out.append((d, n))
    out.sort(key=lambda x: x[0])
    return out[:k]


# ============================================================
# Restkern-Fasern
# ============================================================

def fiber_by_kernel(records: list[dict], kernel_value: int) -> list[dict]:
    return [r for r in records if r["mB"] == kernel_value]


def fiber_map(records: list[dict]) -> dict[int, list[int]]:
    out = defaultdict(list)
    for r in records:
        out[r["mB"]].append(r["n"])
    return dict(out)


def print_fiber(records: list[dict], kernel_value: int, max_rows: int | None = None) -> None:
    fib = fiber_by_kernel(records, kernel_value)
    fib_sorted = sorted(fib, key=lambda r: r["n"])
    print(f"\nRestkern-Faser mB = {kernel_value}:")
    print_records(fib_sorted, max_rows=max_rows)


def top_fibers(records: list[dict], min_size: int = 5, top_k: int = 20) -> list[tuple[int, int]]:
    fmap = fiber_map(records)
    pairs = [(mB, len(vals)) for mB, vals in fmap.items() if len(vals) >= min_size]
    pairs.sort(key=lambda x: (-x[1], x[0]))
    return pairs[:top_k]


# ============================================================
# Cluster nach Familiensupport
# ============================================================

def cluster_by_family_support(records: list[dict]) -> dict[str, list[int]]:
    out = defaultdict(list)
    for r in records:
        out[r["fam_support"]].append(r["n"])
    return dict(out)


# ============================================================
# CSV-Export
# ============================================================

def export_records_to_csv(records: list[dict], filepath: str | Path) -> None:
    filepath = Path(filepath)

    fieldnames = [
        "n",
        "factorization_str",
        "sB",
        "mB",
        "kL",
        "sigma",
        "epsilon",
        "kappa",
        "fam_support",
        "fam_E_part",
        "fam_A_part",
        "fam_B_part",
        "fam_C_part",
        "rho_E",
        "rho_A",
        "rho_B",
        "rho_C",
        "C_fam",
        "K",
        "K2",
        "K3",
    ]

    with filepath.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in records:
            row = {key: r[key] for key in fieldnames}
            writer.writerow(row)


def export_plot_data(records: list[dict], filepath: str | Path) -> None:
    """
    CSV für spätere Plots, z.B.
    x = epsilon, y = C_fam, Farbe = kappa, Größe = K
    """
    filepath = Path(filepath)
    fieldnames = [
        "n",
        "epsilon",
        "kappa",
        "C_fam",
        "K",
        "K2",
        "K3",
        "rho_E",
        "rho_A",
        "rho_B",
        "rho_C",
        "mB",
        "sB",
        "fam_support",
    ]
    with filepath.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in records:
            row = {key: r[key] for key in fieldnames}
            writer.writerow(row)


# ============================================================
# Berichtsfunktionen
# ============================================================

def print_neighbors_report(
    target: int,
    records: list[dict],
    basis: set[int],
    lam_small: int,
    k: int = 12
) -> None:
    neigh = nearest_neighbors_v2(
        target, records, basis, lam_small,
        k=k,
        w_eps=W_EPS, w_kappa=W_KAPPA, w_rho=W_RHO
    )
    print(f"\nNächste strukturelle Nachbarn v2 von {target}:")
    for d, n in neigh:
        print(f"n={n:5d}, d_str2={d:.5f}")


def print_family_clusters(records: list[dict], min_size: int = 5) -> None:
    clusters = cluster_by_family_support(records)
    print("\nCluster nach Familiensupport:")
    for key in sorted(clusters.keys(), key=lambda x: (len(x), x)):
        vals = clusters[key]
        if len(vals) >= min_size:
            print(f"{key:>5}: {len(vals):4d} Zahlen")


# ============================================================
# Hauptprogramm
# ============================================================

def main():
    N = DEFAULT_MAX_N
    basis = BASIS
    lam_small = LAMBDA_SMALL

    print(f"Analysiere Zahlen bis N={N}, Basis={basis}, Lambda={lam_small}\n")

    # Beispielzahlen
    demo_numbers = [360, 420, 441, 455, 540, 693, 945]
    demo_records = [structure_record(n, basis, lam_small) for n in demo_numbers]

    print("Beispielzahlen:\n")
    print_records(demo_records)

    # Gesamtdaten
    records = records_up_to(N, basis, lam_small)

    print("\nTop nach K:")
    print_records(top_by(records, "K", k=15))

    print("\nTop nach K2:")
    print_records(top_by(records, "K2", k=15))

    print("\nTop nach C_fam:")
    print_records(top_by(records, "C_fam", k=15))

    print("\nTop nach epsilon:")
    print_records(top_by(records, "epsilon", k=15))

    print("\nTop nach kappa:")
    print_records(top_by(records, "kappa", k=15))

    # Faserbeispiele
    print_fiber(records, kernel_value=35, max_rows=30)
    print_fiber(records, kernel_value=77, max_rows=30)
    print_fiber(records, kernel_value=143, max_rows=30)

    # Größte Fasern
    print("\nGrößte Restkern-Fasern:")
    for mB, size in top_fibers(records, min_size=5, top_k=20):
        print(f"mB={mB:5d} -> Größe {size}")

    # Cluster
    print_family_clusters(records, min_size=10)

    # Nachbarn
    print_neighbors_report(420, records, basis, lam_small, k=12)
    print_neighbors_report(455, records, basis, lam_small, k=12)

    # CSV-Export
    export_records_to_csv(records, OUTPUT_DIR / "arith_struktur_bis_1000.csv")
    export_plot_data(records, OUTPUT_DIR / "arith_plotdaten_bis_1000.csv")
    print("\nCSV-Dateien geschrieben:")
    print(" - arith_struktur_bis_1000.csv")
    print(" - arith_plotdaten_bis_1000.csv")


if __name__ == "__main__":
    main()