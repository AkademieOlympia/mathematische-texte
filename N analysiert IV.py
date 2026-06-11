import math
import csv
from collections import defaultdict
from pathlib import Path

# ============================================================
# Einstellungen
# ============================================================

N = 100000
BASIS = {2, 3}
KAPPA_SCALES = [5, 7, 11, 13]

W_EPS = 1.0
W_RHO = 1.0
W_KAPPA = 1.0

OUTPUT_DIR = Path(".")


# ============================================================
# SPF-Sieb: kleinster Primfaktor
# ============================================================

def build_spf(limit: int) -> list[int]:
    spf = list(range(limit + 1))
    spf[0] = 0
    if limit >= 1:
        spf[1] = 1

    for i in range(2, int(limit ** 0.5) + 1):
        if spf[i] == i:
            step = i
            start = i * i
            for j in range(start, limit + 1, step):
                if spf[j] == j:
                    spf[j] = i
    return spf


SPF = build_spf(N)


# ============================================================
# Grundfunktionen
# ============================================================

def factorint_from_spf(n: int, spf: list[int]) -> dict[int, int]:
    fac: dict[int, int] = {}
    x = n
    while x > 1:
        p = spf[x]
        fac[p] = fac.get(p, 0) + 1
        x //= p
    return fac


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
        raise ValueError("log nur für x > 0")
    return math.log(x)


def euclidean_distance(v1, v2) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(v1, v2)))


# ============================================================
# Mod-12-Familien
# ============================================================

def family_label_mod12(p: int) -> str | None:
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


# ============================================================
# Strukturrecord für eine Zahl
# ============================================================

def structure_record(n: int) -> dict:
    fac = factorint_from_spf(n, SPF)

    sB = 1
    mB = 1
    kappa_parts = {lam: 1 for lam in KAPPA_SCALES}

    fam_parts = {"E": 1, "A": 1, "B": 1, "C": 1}

    for p, e in fac.items():
        pe = p ** e

        if p in BASIS:
            sB *= pe
        else:
            mB *= pe
            fam = family_label_mod12(p)
            if fam is not None:
                fam_parts[fam] *= pe

        for lam in KAPPA_SCALES:
            if p <= lam:
                kappa_parts[lam] *= pe

    ln = log_safe(n)
    lsB = log_safe(sB) if sB > 1 else 0.0
    lmB = log_safe(mB) if mB > 1 else 0.0

    sigma = lsB / ln if n > 1 else 0.0
    epsilon = lmB / ln if mB > 1 else 0.0

    kappas = {}
    for lam in KAPPA_SCALES:
        val = kappa_parts[lam]
        kappas[lam] = log_safe(val) / ln if val > 1 else 0.0

    # Familiengewichte
    if mB > 1:
        rho = {}
        for fam in ("E", "A", "B", "C"):
            val = fam_parts[fam]
            rho[fam] = (log_safe(val) / lmB) if val > 1 else 0.0
    else:
        rho = {"E": 0.0, "A": 0.0, "B": 0.0, "C": 0.0}

    # Entropie
    C_fam = 0.0
    for v in rho.values():
        if v > 0:
            C_fam -= v * math.log(v)

    # Krümmungen
    K_bal = epsilon * C_fam if mB > 1 else 0.0
    K_mass = (lmB / log_safe(1 + sB)) * C_fam if mB > 1 else 0.0
    K_small = K_bal * kappas[13] if mB > 1 else 0.0

    fam_support = "·".join([fam for fam in ("E", "A", "B", "C") if fam_parts[fam] > 1])
    if not fam_support:
        fam_support = "1"

    return {
        "n": n,
        "factorization_str": factorization_to_string(fac),
        "sB": sB,
        "mB": mB,
        "sigma": sigma,
        "epsilon": epsilon,
        "kappa_5": kappas[5],
        "kappa_7": kappas[7],
        "kappa_11": kappas[11],
        "kappa_13": kappas[13],
        "fam_support": fam_support,
        "rho_E": rho["E"],
        "rho_A": rho["A"],
        "rho_B": rho["B"],
        "rho_C": rho["C"],
        "C_fam": C_fam,
        "K_bal": K_bal,
        "K_mass": K_mass,
        "K_small": K_small,
        "family_vector": (rho["E"], rho["A"], rho["B"], rho["C"]),
        "kappa_vector": (kappas[5], kappas[7], kappas[11], kappas[13]),
    }


# ============================================================
# Distanz v3
# ============================================================

def structural_distance_v3(rec1: dict, rec2: dict) -> float:
    e1, e2 = rec1["epsilon"], rec2["epsilon"]
    r1, r2 = rec1["family_vector"], rec2["family_vector"]
    k1, k2 = rec1["kappa_vector"], rec2["kappa_vector"]

    return math.sqrt(
        W_EPS * (e1 - e2) ** 2
        + W_RHO * euclidean_distance(r1, r2) ** 2
        + W_KAPPA * euclidean_distance(k1, k2) ** 2
    )


# ============================================================
# Gesamtdaten
# ============================================================

def build_records(limit: int) -> list[dict]:
    return [structure_record(n) for n in range(2, limit + 1)]


def top_by(records: list[dict], key: str, k: int = 20) -> list[dict]:
    return sorted(records, key=lambda r: r[key], reverse=True)[:k]


def print_records(records: list[dict], max_rows: int | None = None) -> None:
    rows = records if max_rows is None else records[:max_rows]
    header = (
        f"{'n':>8} | {'Faktorisierung':<24} | {'sB':>8} | {'mB':>8} | {'Fam':<8} | "
        f"{'eps':>6} | {'k5':>5} | {'k7':>5} | {'k11':>5} | {'k13':>5} | "
        f"{'C_fam':>7} | {'K_bal':>7} | {'K_mass':>8} | {'K_small':>8}"
    )
    print(header)
    print("-" * len(header))
    for r in rows:
        print(
            f"{r['n']:8d} | "
            f"{r['factorization_str']:<24} | "
            f"{r['sB']:8d} | {r['mB']:8d} | {r['fam_support']:<8} | "
            f"{r['epsilon']:6.3f} | "
            f"{r['kappa_5']:5.3f} | {r['kappa_7']:5.3f} | {r['kappa_11']:5.3f} | {r['kappa_13']:5.3f} | "
            f"{r['C_fam']:7.3f} | {r['K_bal']:7.3f} | {r['K_mass']:8.3f} | {r['K_small']:8.3f}"
        )


# ============================================================
# Fasern / Cluster / Nachbarn
# ============================================================

def fiber_by_kernel(records: list[dict], kernel_value: int) -> list[dict]:
    return [r for r in records if r["mB"] == kernel_value]


def print_fiber(records: list[dict], kernel_value: int, max_rows: int | None = None) -> None:
    fib = sorted(fiber_by_kernel(records, kernel_value), key=lambda r: r["n"])
    print(f"\nRestkern-Faser mB = {kernel_value}:")
    print_records(fib, max_rows=max_rows)


def top_fibers(records: list[dict], min_size: int = 5, top_k: int = 20):
    fmap = defaultdict(list)
    for r in records:
        fmap[r["mB"]].append(r["n"])
    pairs = [(mB, len(vals)) for mB, vals in fmap.items() if len(vals) >= min_size]
    pairs.sort(key=lambda x: (-x[1], x[0]))
    return pairs[:top_k]


def cluster_by_family_support(records: list[dict]) -> dict[str, int]:
    counts = defaultdict(int)
    for r in records:
        counts[r["fam_support"]] += 1
    return dict(sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])))


def nearest_neighbors(records: list[dict], target_n: int, k: int = 12):
    target = next(r for r in records if r["n"] == target_n)
    out = []
    for r in records:
        if r["n"] == target_n:
            continue
        d = structural_distance_v3(target, r)
        out.append((d, r["n"]))
    out.sort(key=lambda x: x[0])
    return out[:k]


# ============================================================
# Export
# ============================================================

def export_csv(records: list[dict], path: str | Path) -> None:
    path = Path(path)
    fieldnames = [
        "n", "factorization_str", "sB", "mB", "sigma", "epsilon",
        "kappa_5", "kappa_7", "kappa_11", "kappa_13",
        "fam_support", "rho_E", "rho_A", "rho_B", "rho_C",
        "C_fam", "K_bal", "K_mass", "K_small"
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in records:
            writer.writerow({k: r[k] for k in fieldnames})


# ============================================================
# Main
# ============================================================

def main():
    print(f"Analysiere Zahlen bis N={N}, Basis={BASIS}, Kappa-Skalen={KAPPA_SCALES}\n")

    records = build_records(N)

    demo_numbers = [360, 420, 441, 455, 540, 693, 945]
    demo = [next(r for r in records if r["n"] == n) for n in demo_numbers]

    print("Beispielzahlen:\n")
    print_records(demo)

    print("\nTop nach K_bal:\n")
    print_records(top_by(records, "K_bal", 20))

    print("\nTop nach K_mass:\n")
    print_records(top_by(records, "K_mass", 20))

    print("\nTop nach K_small:\n")
    print_records(top_by(records, "K_small", 20))

    print_fiber(records, 35, 20)
    print_fiber(records, 77, 20)
    print_fiber(records, 143, 20)

    print("\nGrößte Restkern-Fasern:")
    for mB, size in top_fibers(records, min_size=8, top_k=25):
        print(f"mB={mB:8d} -> Größe {size}")

    print("\nCluster nach Familiensupport:")
    clusters = cluster_by_family_support(records)
    for key, count in list(clusters.items())[:20]:
        print(f"{key:>8}: {count:6d}")

    print("\nNächste strukturelle Nachbarn von 420:")
    for d, n in nearest_neighbors(records, 420, 15):
        print(f"n={n:8d}, d_str3={d:.6f}")

    print("\nNächste strukturelle Nachbarn von 455:")
    for d, n in nearest_neighbors(records, 455, 15):
        print(f"n={n:8d}, d_str3={d:.6f}")

    export_csv(records, OUTPUT_DIR / "arith_struktur_v3_bis_100000.csv")
    print("\nCSV geschrieben: arith_struktur_v3_bis_100000.csv")


if __name__ == "__main__":
    main()