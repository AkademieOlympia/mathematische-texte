# pyright: reportMissingImports=false
# SageMath 10.5: Carnot-Ptolemaeus-Validierung im eabc-Modell
import csv
from collections import deque
from functools import lru_cache
import math
import os
import random
import subprocess
import sys
import time
from datetime import datetime

SAGE_BIN = "/Applications/SageMath-10-8.app/Contents/Frameworks/Sage.framework/Versions/Current/local/bin/sage"

try:
    from sage.all import (
        acos,
        graphics_array,
        I,
        OctonionAlgebra,
        QQ,
        QuaternionAlgebra,
        RR,
        ZZ,
        bernoulli,
        four_squares,
        is_prime,
        log,
        matrix,
        matrix_plot,
        next_prime,
        pi,
        prime_divisors,
        sqrt,
        text,
        zeta,
    )
except ModuleNotFoundError:
    if __name__ == "__main__" and os.environ.get("CARNOT_REEXEC_UNDER_SAGE") != "1":
        env = os.environ.copy()
        env["CARNOT_REEXEC_UNDER_SAGE"] = "1"
        completed = subprocess.run([SAGE_BIN, os.path.abspath(__file__), *sys.argv[1:]], env=env)
        raise SystemExit(completed.returncode)
    raise RuntimeError(
        "Dieses Skript benoetigt SageMath. Starte es mit 'sage Carnot.py' "
        "oder direkt mit dem SageMath-10.8-Binaer."
    )

# Quaternionen-Algebra definieren
Q = QuaternionAlgebra(QQ, -1, -1)
QR = QuaternionAlgebra(RR, -1, -1)
O = OctonionAlgebra(QQ, -1, -1, -1)

# Vier symmetrische Testzustände des Carnot-Ptolemaeus-Zirkels
ZUSTAENDE = (
    Q([1, 0, 0, 0]),   # A
    Q([0, 1, 0, 0]),   # B
    Q([-1, 0, 0, 0]),  # C
    Q([0, -1, 0, 0]),  # D
)

def ptolemy_check(A, B, C, D):
    """Prüft die Sehnenviereck-Bedingung für vier quaternionische Zustände."""
    # Seitenlängen (Distanzen zwischen den Zuständen)
    a = sqrt((B - A).reduced_norm()) # Seite AB
    b = sqrt((C - B).reduced_norm()) # Seite BC
    c = sqrt((D - C).reduced_norm()) # Seite CD
    d = sqrt((A - D).reduced_norm()) # Seite DA
    
    # Diagonalen
    e = sqrt((C - A).reduced_norm()) # Diagonale AC
    f = sqrt((D - B).reduced_norm()) # Diagonale BD
    
    ptolemy_sum = a*c + b*d
    diagonal_prod = e*f
    
    diff = abs(ptolemy_sum - diagonal_prod)
    return ptolemy_sum, diagonal_prod, diff

def berechne_system_resonanz(rotation_offset, zustaende=ZUSTAENDE):
    """Summiert |zeta(1/2 + i*t)| ueber alle gegebenen Zustaende."""
    gesamt_resonanz = RR(0)
    for q in zustaende:
        # Fuer Quaternionen-Elemente ist die reduzierte Norm die stabile Groesse.
        norm_wert = RR(q.reduced_norm()) * RR(rotation_offset)
        s = RR(1) / 2 + I * norm_wert
        gesamt_resonanz += abs(zeta(s))

    return gesamt_resonanz


def scan_energieniveaus(start=14.0, ende=30.0, schritte=100, zustaende=ZUSTAENDE):
    """Sucht im vorgegebenen Intervall das Resonanzminimum."""
    if schritte < 2:
        raise ValueError("schritte muss mindestens 2 sein.")

    delta = (RR(ende) - RR(start)) / (schritte - 1)
    energieniveaus = [RR(start) + n * delta for n in range(schritte)]
    ergebnisse = [(e, berechne_system_resonanz(e, zustaende)) for e in energieniveaus]
    optimal_e, min_res = min(ergebnisse, key=lambda x: x[1])
    return optimal_e, min_res, ergebnisse


def find_prime_clusters(limit_triplets=100, limit_quads=50):
    """Findet Primzahldrillinge und Primzahlvierlinge in Standardmustern."""
    triplets = []
    quads = []
    p = 2

    while len(triplets) < limit_triplets or len(quads) < limit_quads:
        p = next_prime(p)

        # Klassisches Vierlingsmuster: (p, p+2, p+6, p+8)
        if is_prime(p + 2) and is_prime(p + 6) and is_prime(p + 8):
            if len(quads) < limit_quads:
                quads.append((p, p + 2, p + 6, p + 8))

        # Zulaessiges Drillingsmuster: (p, p+2, p+6)
        if is_prime(p + 2) and is_prime(p + 6):
            if len(triplets) < limit_triplets:
                triplets.append((p, p + 2, p + 6))

    return triplets, quads


def get_zeta_signature_weights(M=500):
    """
    Schaetzt EABC-Gewichte aus Primteilern der Bernoulli-Nenner.
    Primzahlen >= 5 werden ueber ihre Restklasse mod 12 klassifiziert.
    """
    if M < 1:
        raise ValueError("M muss mindestens 1 sein.")

    counts = {"A": 0, "B": 0, "C": 0, "E": 0}
    for m in range(1, M + 1):
        denominator = bernoulli(2 * m).denom()
        for prime in prime_divisors(denominator):
            if prime >= 5:
                residue = int(prime % 12)
                if residue == 5:
                    counts["A"] += 1
                elif residue == 7:
                    counts["B"] += 1
                elif residue == 11:
                    counts["C"] += 1
                elif residue == 1:
                    counts["E"] += 1

    total = sum(counts.values())
    if total == 0:
        return {key: 0.0 for key in counts}
    return {key: counts[key] / total for key in counts}


def residue_signature_mod12(n):
    """Ordnet eine Zahl ihrer mod-12-Signatur zu."""
    residue = int(n % 12)
    if residue == 5:
        return "A"
    if residue == 7:
        return "B"
    if residue == 11:
        return "C"
    if residue == 1:
        return "E"
    return None


def get_signature_open_probabilities(weights, base_open=0.20, spread=0.50):
    """
    Leitet aus den Zeta-Signaturen Oeffnungswahrscheinlichkeiten ab.
    Die Gesamtbelegung bleibt moderat, damit Perkolation nicht trivial wird.
    """
    max_weight = max(weights.values()) if weights and max(weights.values()) > 0 else 1.0
    probabilities = {}
    for signature in ("A", "B", "C", "E"):
        normalized = weights.get(signature, 0.0) / max_weight
        probabilities[signature] = min(0.95, base_open + spread * normalized)
    probabilities["default"] = base_open
    return probabilities


def build_signature_percolation_grid(height, tower_width=48, weights=None, weighted=True, rng=None):
    """
    Baut ein probabilistisches Gitter fuer die mod-12-Signaturen.
    Beim gewichteten Modell haengt die Oeffnung von der Zeta-Signatur ab,
    beim Zufallsmodell nur von derselben globalen Dichte.
    """
    if height < 1 or tower_width < 1:
        raise ValueError("height und tower_width muessen positiv sein.")

    rng = rng or random.Random()
    weights = weights or get_zeta_signature_weights()
    probabilities = get_signature_open_probabilities(weights)
    mean_open = (sum(probabilities[sig] for sig in ("A", "B", "C", "E")) + 8 * probabilities["default"]) / 12.0

    grid = []
    for row in range(height):
        grid_row = []
        for col in range(tower_width):
            n = row * tower_width + col + 1
            signature = residue_signature_mod12(n)
            if weighted:
                open_probability = probabilities.get(signature, probabilities["default"])
            else:
                open_probability = mean_open
            grid_row.append(rng.random() < open_probability)
        grid.append(grid_row)
    return grid


def coarse_grain_percolation_grid(open_grid, block_size):
    """
    Fasst ein Gitter zu Bloecken zusammen.
    Ein Block bleibt nur dann offen, wenn alle Zellen darin offen sind.
    """
    if block_size < 1:
        raise ValueError("block_size muss mindestens 1 sein.")

    nrows = len(open_grid)
    ncols = len(open_grid[0]) if open_grid else 0
    coarse = []
    for row_start in range(0, nrows, block_size):
        coarse_row = []
        for col_start in range(0, ncols, block_size):
            block_is_open = True
            for row in range(row_start, min(row_start + block_size, nrows)):
                for col in range(col_start, min(col_start + block_size, ncols)):
                    if not open_grid[row][col]:
                        block_is_open = False
                        break
                if not block_is_open:
                    break
            coarse_row.append(block_is_open)
        coarse.append(coarse_row)
    return coarse


def has_top_bottom_path(open_grid):
    """Prueft per BFS, ob ein Weg von oben nach unten existiert."""
    if not open_grid or not open_grid[0]:
        return False

    nrows = len(open_grid)
    ncols = len(open_grid[0])
    queue = deque()
    seen = set()

    for col in range(ncols):
        if open_grid[0][col]:
            queue.append((0, col))
            seen.add((0, col))

    while queue:
        row, col = queue.popleft()
        if row == nrows - 1:
            return True
        for drow, dcol in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            next_row = row + drow
            next_col = col + dcol
            if 0 <= next_row < nrows and 0 <= next_col < ncols:
                if open_grid[next_row][next_col] and (next_row, next_col) not in seen:
                    seen.add((next_row, next_col))
                    queue.append((next_row, next_col))
    return False


def find_min_b_percolation(z, tower_width=48, weights=None, weighted=True, rng=None):
    """
    Findet die kleinste Blockskala b, auf der keine Top-Bottom-Perkolation
    mehr moeglich ist.
    """
    open_grid = build_signature_percolation_grid(
        z,
        tower_width=tower_width,
        weights=weights,
        weighted=weighted,
        rng=rng,
    )
    max_block = min(z, tower_width)
    for block_size in range(1, max_block + 1):
        coarse_grid = coarse_grain_percolation_grid(open_grid, block_size)
        if not has_top_bottom_path(coarse_grid):
            return block_size
    return None


def simulate_zeta_impact(z, N_samples=50, tower_width=48, M=500, seed=0):
    """
    Vergleicht ein zeta-gewichtetes Perkolationsmodell mit einem gleich dichten
    Zufallsmodell. Dies ist ein erster Proxy fuer die Phase-4-Idee.
    """
    if N_samples < 1:
        raise ValueError("N_samples muss mindestens 1 sein.")

    weights = get_zeta_signature_weights(M)
    weighted_b_values = []
    random_b_values = []

    for sample in range(N_samples):
        weighted_rng = random.Random(seed + sample)
        random_rng = random.Random(seed + 10_000 + sample)
        weighted_b = find_min_b_percolation(
            z,
            tower_width=tower_width,
            weights=weights,
            weighted=True,
            rng=weighted_rng,
        )
        random_b = find_min_b_percolation(
            z,
            tower_width=tower_width,
            weights=weights,
            weighted=False,
            rng=random_rng,
        )
        weighted_b_values.append(weighted_b if weighted_b is not None else min(z, tower_width) + 1)
        random_b_values.append(random_b if random_b is not None else min(z, tower_width) + 1)

    weighted_mean = sum(weighted_b_values) / len(weighted_b_values)
    random_mean = sum(random_b_values) / len(random_b_values)
    weighted_std = (
        sum((value - weighted_mean) ** 2 for value in weighted_b_values) / len(weighted_b_values)
    ) ** 0.5
    random_std = (
        sum((value - random_mean) ** 2 for value in random_b_values) / len(random_b_values)
    ) ** 0.5

    return {
        "weights": weights,
        "weighted_b_values": weighted_b_values,
        "random_b_values": random_b_values,
        "weighted_mean_b": weighted_mean,
        "random_mean_b": random_mean,
        "weighted_std_b": weighted_std,
        "random_std_b": random_std,
        "delta_mean_b": weighted_mean - random_mean,
    }


def run_full_carnot_model(W, z, N, H=None, mode=None):
    """
    Verdichtet die vorhandenen Carnot-Bausteine zu einem einzelnen Modelllauf
    fuer Invarianz- und Skalierungstests.
    """
    if H is None:
        H = 10 * int(z)
    if W < 1 or z < 1 or N < 1 or H < 1:
        raise ValueError("W, z, N und H muessen positiv sein.")

    started_at = time.time()
    triplet_limit = min(max(int(W + z + 50), 200), 4000)
    triplets, _ = find_prime_clusters(limit_triplets=triplet_limit, limit_quads=0)
    cluster = triplets[int(N) % len(triplets)]

    rotation_offset = RR(H) / RR(max(1, z))
    states = tuple(QR(state.coefficient_tuple()) for state in cluster_to_quaternions(cluster))
    zeta_resonance_sum = sum(zeta_resonance(state, rotation_offset) for state in states)

    iterations = min(120, max(25, int(W // 2)))
    virtual_result = find_virtual_complement(cluster, iterations=iterations)

    weights = get_zeta_signature_weights(min(500, max(50, int(z))))
    bmax = find_min_b_percolation(
        z=int(H),
        tower_width=int(W),
        weights=weights,
        weighted=True,
        rng=random.Random(int(N) % (2 ** 31)),
    )
    if bmax is None:
        bmax = min(int(H), int(W)) + 1

    healed_states = states + (virtual_result["D_virt"],)
    phase = float(calculate_phase_shift(healed_states))

    return {
        "cluster": cluster,
        "torsion": float(calculate_octonion_torsion_active(cluster)),
        "bmax": int(bmax),
        "ptolemy_defect": float(virtual_result["start_error"]),
        "ptolemy_diff": float(virtual_result["start_error"]),
        "zeta_resonance": float(zeta_resonance_sum),
        "start_error": float(virtual_result["start_error"]),
        "end_error": float(virtual_result["end_error"]),
        "rotation_offset": float(rotation_offset),
        "phase": phase,
        "mode": mode or "default",
        "execution_time": time.time() - started_at,
    }


def _group_experiment_rows(rows):
    grouped = {}
    for row in rows:
        key = (row["W"], row["z"], row["R"])
        grouped.setdefault(key, []).append(row)
    return grouped


def summarize_tauB_invariance(rows):
    """Aggregiert die Rohdaten des tauB-Experiments je (W, z, R)."""
    summary_rows = []
    for (W, z, R), group in sorted(_group_experiment_rows(rows).items()):
        torsion_values = [row["torsion"] for row in group]
        bmax_values = [row["bmax"] for row in group]
        end_error_values = [row["end_error"] for row in group]

        torsion_mean = sum(torsion_values) / len(torsion_values)
        torsion_variance = sum((value - torsion_mean) ** 2 for value in torsion_values) / len(torsion_values)
        torsion_sorted = sorted(torsion_values)
        mid = len(torsion_sorted) // 2
        if len(torsion_sorted) % 2 == 1:
            torsion_median = torsion_sorted[mid]
        else:
            torsion_median = (torsion_sorted[mid - 1] + torsion_sorted[mid]) / 2

        summary_rows.append(
            {
                "W": W,
                "z": z,
                "R": R,
                "torsion_mean": torsion_mean,
                "torsion_std": torsion_variance ** 0.5,
                "torsion_median": torsion_median,
                "bmax_mean": sum(bmax_values) / len(bmax_values),
                "bmax_max": max(bmax_values),
                "end_error_mean": sum(end_error_values) / len(end_error_values),
            }
        )
    return summary_rows


def summarize_tauB_by_W(rows):
    """Aggregiert die tauB-Rohdaten nur nach W fuer den Invarianz-Check."""
    grouped = {}
    for row in rows:
        grouped.setdefault(row["W"], []).append(row)

    summary_rows = []
    for W, group in sorted(grouped.items()):
        torsion_values = [row["torsion"] for row in group]
        bmax_values = [row["bmax"] for row in group]
        end_error_values = [row["end_error"] for row in group]

        torsion_mean = sum(torsion_values) / len(torsion_values)
        torsion_variance = sum((value - torsion_mean) ** 2 for value in torsion_values) / len(torsion_values)
        torsion_sorted = sorted(torsion_values)
        mid = len(torsion_sorted) // 2
        if len(torsion_sorted) % 2 == 1:
            torsion_median = torsion_sorted[mid]
        else:
            torsion_median = (torsion_sorted[mid - 1] + torsion_sorted[mid]) / 2

        healing_success_rate = sum(1 for value in end_error_values if value < 1e-3) / len(end_error_values)
        summary_rows.append(
            {
                "W": int(W),
                "mu_torsion": torsion_mean,
                "sigma_torsion": torsion_variance ** 0.5,
                "median_torsion": torsion_median,
                "bmax_worst": max(bmax_values),
                "bmax_ratio": (sum(bmax_values) / len(bmax_values)) / W,
                "healing_success_rate": healing_success_rate,
                "tau_hat_W": torsion_mean / W,
                "tau_log_W": torsion_mean / float(log(RR(W))),
            }
        )
    return summary_rows


def _format_summary_table(rows, columns):
    """Formatiert eine kleine Tabelle fuer Konsolenausgaben ohne Zusatzbibliotheken."""
    widths = {}
    for key, header in columns:
        values = [header] + [str(row[key]) for row in rows]
        widths[key] = max(len(value) for value in values)

    lines = []
    header_line = " | ".join(header.ljust(widths[key]) for key, header in columns)
    separator_line = "-+-".join("-" * widths[key] for key, _ in columns)
    lines.append(header_line)
    lines.append(separator_line)
    for row in rows:
        lines.append(" | ".join(str(row[key]).ljust(widths[key]) for key, _ in columns))
    return "\n".join(lines)


def write_csv_table(path, rows, fieldnames):
    """Schreibt eine Liste homogener Dictionaries als CSV-Datei."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def run_tauB_invariance_experiment(
    W_values=(48, 96, 192),
    z_values=(307, 503, 701, 1009),
    R_values=(2, 3, 4),
    samples=30,
    seed=0,
    raw_csv_path="tauB_invariance_test.csv",
    summary_csv_path="tauB_invariance_summary.csv",
):
    """
    Fuehrt den parameterisierten tauB-Invarianztest aus und speichert
    Rohdaten sowie eine gruppierte Zusammenfassung als CSV.
    """
    if samples < 1:
        raise ValueError("samples muss mindestens 1 sein.")

    rng = random.Random(seed)
    rows = []

    for W in W_values:
        for z in z_values:
            for R in R_values:
                for _ in range(samples):
                    N = rng.randint(int(z ** R), int(10 * (z ** R)))
                    H = 10 * int(z)
                    result = run_full_carnot_model(W=W, z=z, N=N, H=H)
                    rows.append(
                        {
                            "W": int(W),
                            "z": int(z),
                            "R": int(R),
                            "N": int(N),
                            "H": int(H),
                            "torsion": result["torsion"],
                            "bmax": result["bmax"],
                            "ptolemy_defect": result["ptolemy_defect"],
                            "zeta_resonance": result["zeta_resonance"],
                            "start_error": result["start_error"],
                            "end_error": result["end_error"],
                        }
                    )

    summary_rows = summarize_tauB_invariance(rows)
    write_csv_table(
        raw_csv_path,
        rows,
        [
            "W",
            "z",
            "R",
            "N",
            "H",
            "torsion",
            "bmax",
            "ptolemy_defect",
            "zeta_resonance",
            "start_error",
            "end_error",
        ],
    )
    write_csv_table(
        summary_csv_path,
        summary_rows,
        [
            "W",
            "z",
            "R",
            "torsion_mean",
            "torsion_std",
            "torsion_median",
            "bmax_mean",
            "bmax_max",
            "end_error_mean",
        ],
    )

    return {
        "rows": rows,
        "summary": summary_rows,
        "raw_csv_path": raw_csv_path,
        "summary_csv_path": summary_csv_path,
    }


def run_invariance_analysis():
    """
    Fuehrt den Bamberger Invarianz-Test ueber die bereits eingebettete
    tauB-Experimentlogik aus und schreibt zwei CSV-Berichte.
    """
    W_values = (48, 96, 192)
    z_values = (307, 503, 701, 1009)
    R_values = (2, 3, 4)
    samples_per_config = 30
    total_runs = len(W_values) * len(z_values) * len(R_values) * samples_per_config

    print(f"Starte Phase 5: Invarianz-Eichung ({total_runs} Einzellaeufe)...")
    experiment = run_tauB_invariance_experiment(
        W_values=W_values,
        z_values=z_values,
        R_values=R_values,
        samples=samples_per_config,
        seed=42,
        raw_csv_path="tauB_invariance_raw_data.csv",
        summary_csv_path="tauB_invariance_grouped_by_W_z_R.csv",
    )

    summary_rows = summarize_tauB_by_W(experiment["rows"])
    write_csv_table(
        "tauB_summary_report.csv",
        summary_rows,
        [
            "W",
            "mu_torsion",
            "sigma_torsion",
            "median_torsion",
            "bmax_worst",
            "bmax_ratio",
            "healing_success_rate",
            "tau_hat_W",
            "tau_log_W",
        ],
    )

    printable_rows = []
    for row in summary_rows:
        printable_rows.append(
            {
                "W": row["W"],
                "mu_torsion": f"{row['mu_torsion']:.6f}",
                "sigma_torsion": f"{row['sigma_torsion']:.6f}",
                "median_torsion": f"{row['median_torsion']:.6f}",
                "bmax_worst": row["bmax_worst"],
                "bmax_ratio": f"{row['bmax_ratio']:.6f}",
                "healing_success_rate": f"{row['healing_success_rate']:.4f}",
                "tau_hat_W": f"{row['tau_hat_W']:.6f}",
                "tau_log_W": f"{row['tau_log_W']:.6f}",
            }
        )

    print()
    print("### Bamberger Invarianz-Zusammenfassung ###")
    print(
        _format_summary_table(
            printable_rows,
            [
                ("W", "W"),
                ("mu_torsion", "mu_torsion"),
                ("sigma_torsion", "sigma_torsion"),
                ("median_torsion", "median_torsion"),
                ("bmax_worst", "bmax_worst"),
                ("bmax_ratio", "bmax_ratio"),
                ("healing_success_rate", "healing_success_rate"),
                ("tau_hat_W", "tau_hat_W"),
                ("tau_log_W", "tau_log_W"),
            ],
        )
    )

    mu_values = [row["mu_torsion"] for row in summary_rows]
    tau_hat_values = [row["tau_hat_W"] for row in summary_rows]
    mu_mean = sum(mu_values) / len(mu_values)
    tau_hat_mean = sum(tau_hat_values) / len(tau_hat_values)
    mu_spread = (sum((value - mu_mean) ** 2 for value in mu_values) / len(mu_values)) ** 0.5
    tau_hat_spread = (
        sum((value - tau_hat_mean) ** 2 for value in tau_hat_values) / len(tau_hat_values)
    ) ** 0.5

    print()
    if mu_spread < 0.5:
        print("BEFUND: Fall A (Invarianz bestaetigt). Delta_torsion ist stabil.")
    elif tau_hat_spread < mu_spread:
        print("BEFUND: Fall B (Skalierung detektiert). Delta_torsion ist extensiv.")
    else:
        print("BEFUND: Fall C (Drift/Fluktuation). Kein universeller Parameter.")

    return {
        "raw_rows": experiment["rows"],
        "grouped_summary": experiment["summary"],
        "W_summary": summary_rows,
        "raw_csv_path": experiment["raw_csv_path"],
        "grouped_csv_path": experiment["summary_csv_path"],
        "summary_csv_path": "tauB_summary_report.csv",
    }


def _circular_mean_and_resultant(phases):
    """Berechnet zirkulären Mittelwert und Resultantenlänge für Phasen."""
    if not phases:
        raise ValueError("phases darf nicht leer sein.")

    mean_sin = sum(math.sin(value) for value in phases) / len(phases)
    mean_cos = sum(math.cos(value) for value in phases) / len(phases)
    mean_phi = math.atan2(mean_sin, mean_cos)
    resultant_length = math.sqrt(mean_cos ** 2 + mean_sin ** 2)
    return mean_phi, resultant_length


def run_real_test_matrix(
    W_values=(48, 96, 192),
    z_values=(307, 503, 701, 1009),
    R_values=(2, 3, 4),
    samples_per_config=30,
    seed=42,
):
    """
    Fuehrt den echten Testlauf fuer die Invarianz von Torsion und Phase ueber W aus.
    Schreibt Rohdaten und eine nach W aggregierte Zusammenfassung als CSV.
    """
    if samples_per_config < 1:
        raise ValueError("samples_per_config muss mindestens 1 sein.")

    rng = random.Random(seed)
    total_runs = len(W_values) * len(z_values) * len(R_values) * samples_per_config
    rows = []
    start_total = time.time()

    print(f"--- START ECHTER TEST: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ---")
    print(f"Geplante Einzellaeufe: {total_runs}")

    for W in W_values:
        print(f"Pruefe Aufloesung W = {W}...")
        for z in z_values:
            for R in R_values:
                for _ in range(samples_per_config):
                    N = rng.randint(int(z ** R), int(10 * (z ** R)))
                    H = 10 * int(z)
                    data = run_full_carnot_model(W=W, z=z, N=N, H=H, mode="Full-Tors")
                    rows.append(
                        {
                            "W": int(W),
                            "z": int(z),
                            "R": int(R),
                            "N": int(N),
                            "H": int(H),
                            "torsion": data["torsion"],
                            "phase": data["phase"],
                            "bmax": data["bmax"],
                            "ptolemy_defect": data["ptolemy_defect"],
                            "zeta_resonance": data["zeta_resonance"],
                            "start_error": data["start_error"],
                            "end_error": data["end_error"],
                            "duration": data["execution_time"],
                            "success": data["end_error"] < 1e-4,
                        }
                    )

    raw_filename = f"BM_RealTest_Raw_{int(time.time())}.csv"
    write_csv_table(
        raw_filename,
        rows,
        [
            "W",
            "z",
            "R",
            "N",
            "H",
            "torsion",
            "phase",
            "bmax",
            "ptolemy_defect",
            "zeta_resonance",
            "start_error",
            "end_error",
            "duration",
            "success",
        ],
    )

    grouped = {}
    for row in rows:
        grouped.setdefault(row["W"], []).append(row)

    summary_rows = []
    for W, group in sorted(grouped.items()):
        torsion_values = [row["torsion"] for row in group]
        phase_values = [row["phase"] for row in group]
        bmax_values = [row["bmax"] for row in group]
        success_values = [1.0 if row["success"] else 0.0 for row in group]

        mean_phi, phase_R = _circular_mean_and_resultant(phase_values)
        torsion_mean = sum(torsion_values) / len(torsion_values)
        torsion_variance = sum((value - torsion_mean) ** 2 for value in torsion_values) / len(torsion_values)

        summary_rows.append(
            {
                "W": int(W),
                "mu_torsion": torsion_mean,
                "sigma_torsion": torsion_variance ** 0.5,
                "phi_circ_mean": mean_phi,
                "phase_R": phase_R,
                "theta_porosity": (sum(bmax_values) / len(bmax_values)) / W,
                "healing_success_rate": sum(success_values) / len(success_values),
            }
        )

    summary_filename = "Bamberger_Invarianz_Final.csv"
    write_csv_table(
        summary_filename,
        summary_rows,
        [
            "W",
            "mu_torsion",
            "sigma_torsion",
            "phi_circ_mean",
            "phase_R",
            "theta_porosity",
            "healing_success_rate",
        ],
    )

    printable_rows = []
    for row in summary_rows:
        printable_rows.append(
            {
                "W": row["W"],
                "mu_torsion": f"{row['mu_torsion']:.6f}",
                "sigma_torsion": f"{row['sigma_torsion']:.6f}",
                "phi_circ_mean": f"{row['phi_circ_mean']:.6f}",
                "phase_R": f"{row['phase_R']:.6f}",
                "theta_porosity": f"{row['theta_porosity']:.6f}",
                "healing_success_rate": f"{row['healing_success_rate']:.4f}",
            }
        )

    print()
    print("### ECHTER TEST ERGEBNISSE ###")
    print(
        _format_summary_table(
            printable_rows,
            [
                ("W", "W"),
                ("mu_torsion", "mu_torsion"),
                ("sigma_torsion", "sigma_torsion"),
                ("phi_circ_mean", "phi_circ_mean"),
                ("phase_R", "phase_R"),
                ("theta_porosity", "theta_porosity"),
                ("healing_success_rate", "healing_success_rate"),
            ],
        )
    )

    duration_total = time.time() - start_total
    print()
    print(f"Daten in {raw_filename} und {summary_filename} gesichert. Dauer: {duration_total:.2f}s")

    return {
        "rows": rows,
        "summary": summary_rows,
        "raw_csv_path": raw_filename,
        "summary_csv_path": summary_filename,
        "duration_seconds": duration_total,
    }


@lru_cache(maxsize=None)
def quat_representation(n):
    """Gibt eine Darstellung von n als Summe von vier Quadraten zurueck."""
    return tuple(int(x) for x in four_squares(int(n)))


def integer_to_quaternion(n):
    """Bettet eine Zahl ueber ihre Vier-Quadrate-Darstellung in Q ein."""
    return Q(quat_representation(n))


def prime_octonion(p_quad):
    """Mappt einen Primzahlvierling auf ein Oktonion in den ersten vier Komponenten."""
    if len(p_quad) != 4:
        raise ValueError("p_quad muss genau vier Eintraege enthalten.")
    return O([p_quad[0], p_quad[1], p_quad[2], p_quad[3], 0, 0, 0, 0])


def octonion_associator(a, b, c):
    """Misst die Nicht-Assoziativitaet dreier Oktonionen."""
    return (a * b) * c - a * (b * c)


def cluster_to_octonion(cluster):
    """Bettet ein Cluster in den oktonionischen Raum ein und fuellt mit Nullen auf."""
    if len(cluster) > 8:
        raise ValueError("cluster darf hoechstens 8 Eintraege enthalten.")
    return O(list(cluster) + [0] * (8 - len(cluster)))


def get_octonion_torsion(cluster_a, cluster_b, cluster_c):
    """Misst die oktonionische Torsion als Norm des Assoziators dreier Cluster."""
    o1 = cluster_to_octonion(cluster_a)
    o2 = cluster_to_octonion(cluster_b)
    o3 = cluster_to_octonion(cluster_c)
    return RR(octonion_associator(o1, o2, o3).norm())


def calculate_octonion_torsion(cluster):
    """
    Berechnet eine lokale oktonionische Torsion aus einem Drilling
    ueber drei verschobene 8-dimensionale Einbettungen.
    """
    if len(cluster) != 3:
        raise ValueError("cluster muss genau 3 Eintraege enthalten.")

    p0, p1, p2 = (RR(value) for value in cluster)
    o1 = O([p0, p1, p2, 0, 0, 0, 0, 0])
    o2 = O([0, p0, p1, p2, 0, 0, 0, 0])
    o3 = O([0, 0, p0, p1, p2, 0, 0, 0])

    assoziator = octonion_associator(o1, o2, o3)
    scale = (p0 * p1 * p2) ** (RR(1) / 3)
    if scale == 0:
        raise ValueError("Normierung darf nicht 0 sein.")
    return RR(assoziator.norm() / scale)


def calculate_octonion_torsion_active(cluster):
    """
    Erzwingt den Sprung aus der assoziativen Unteralgebra
    ueber eine nichttriviale Fano-Triade.
    """
    if len(cluster) != 3:
        raise ValueError("cluster muss genau 3 Eintraege enthalten.")

    p0, p1, p2 = (RR(value) for value in cluster)
    o1 = O([0, p0, 0, 0, 0, 0, 0, 0])  # e1
    o2 = O([0, 0, p1, 0, 0, 0, 0, 0])  # e2
    o3 = O([0, 0, 0, 0, p2, 0, 0, 0])  # e4

    assoziator = octonion_associator(o1, o2, o3)
    scale = (p0 * p1 * p2) ** (RR(1) / 3)
    if scale == 0:
        raise ValueError("Normierung darf nicht 0 sein.")
    return RR(assoziator.norm() / scale)


def get_imaginary_vector(n):
    """Extrahiert den imaginaeren Anteil der Vier-Quadrate-Darstellung."""
    _, b, c, d = quat_representation(n)
    return Q([0, b, c, d])


def quaternion_imaginary_part(q):
    """Projiziert eine Quaternion auf ihren rein imaginaeren Anteil."""
    _, b, c, d = q.coefficient_tuple()
    return q.parent()([0, b, c, d])


def cluster_to_quaternions(cluster):
    """Uebersetzt ein Primzahl-Cluster in Quaternionen-Zustaende."""
    return tuple(integer_to_quaternion(n) for n in cluster)


def imag_vector_dot(u, v):
    """Skalarprodukt zweier rein imaginaerer Quaternionen."""
    _, ux, uy, uz = u.coefficient_tuple()
    _, vx, vy, vz = v.coefficient_tuple()
    return ux * vx + uy * vy + uz * vz


def cross_product_area(u, v):
    """Betrag des Kreuzprodukts der imaginaeren Anteile."""
    area_sq = u.reduced_norm() * v.reduced_norm() - imag_vector_dot(u, v) ** 2
    if area_sq < 0:
        area_sq = 0
    return sqrt(area_sq)


def calculate_cycle_work(quadruplet):
    """Summiert die Kreuzprodukt-Flaechen eines geschlossenen Viererzyklus."""
    vectors = [get_imaginary_vector(p) for p in quadruplet]
    work = RR(0)
    for index, current_vector in enumerate(vectors):
        next_vector = vectors[(index + 1) % len(vectors)]
        work += cross_product_area(current_vector, next_vector)
    return work


def calculate_cycle_work_for_states(states):
    """Summiert die Kreuzprodukt-Flaechen eines Zyklus aus Quaternionen-Zustaenden."""
    vectors = [quaternion_imaginary_part(state) for state in states]
    work = RR(0)
    for index, current_vector in enumerate(vectors):
        next_vector = vectors[(index + 1) % len(vectors)]
        work += cross_product_area(current_vector, next_vector)
    return work


def calculate_phase_shift(quadruplet_quats):
    """Berechnet die kumulative Phasenabweichung relativ zu 2*pi."""
    total_phase = RR(0)
    vectors = [quaternion_imaginary_part(q) for q in quadruplet_quats]

    for index, v1 in enumerate(vectors):
        v2 = vectors[(index + 1) % len(vectors)]
        norm_product = v1.reduced_norm() * v2.reduced_norm()
        if norm_product == 0:
            continue

        cos_theta = RR(imag_vector_dot(v1, v2) / sqrt(norm_product))
        cos_theta = min(RR(1), max(RR(-1), cos_theta))
        total_phase += acos(cos_theta)

    return total_phase - 2 * pi


def zeta_resonance(q, rotation_offset=1.0):
    """Misst die Zeta-Resonanz einer Quaternion ueber ihre reduzierte Norm."""
    t = RR(q.reduced_norm()) * RR(rotation_offset)
    s = RR(1) / 2 + I * t
    return abs(zeta(s))


def kanonische_energie_minimierung(starr_8, frei_4, rotation_offset=1.0, epsilon=1e-12):
    """
    Bewertet eine 8+4-Konfiguration.

    `starr_8` definiert den Carnot-Rahmen, `frei_4` die beweglichen Resonanzzustaende.
    """
    if len(starr_8) != 8:
        raise ValueError("starr_8 muss genau 8 Quaternionen enthalten.")
    if len(frei_4) != 4:
        raise ValueError("frei_4 muss genau 4 Quaternionen enthalten.")

    rahmen_energie = calculate_cycle_work_for_states(starr_8)
    resonanz = sum(zeta_resonance(q, rotation_offset) for q in frei_4)
    resonanz = max(RR(epsilon), RR(resonanz))
    return rahmen_energie / resonanz


def analysiere_8plus4_konfiguration(starr_8, frei_4, rotation_offset=1.0, epsilon=1e-12):
    """Gibt alle Kennzahlen der 8+4-Energiedoku-Konfiguration zurueck."""
    rahmen_energie = calculate_cycle_work_for_states(starr_8)
    resonanz = sum(zeta_resonance(q, rotation_offset) for q in frei_4)
    resonanz = max(RR(epsilon), RR(resonanz))
    return {
        "rahmen_energie": rahmen_energie,
        "resonanz": resonanz,
        "effizienz_index": kanonische_energie_minimierung(
            starr_8,
            frei_4,
            rotation_offset=rotation_offset,
            epsilon=epsilon,
        ),
    }


def validate_alpha_coupling(index, mu_target):
    """
    Berechnet eine heuristische Alpha^-1-Vorhersage aus dem eabc-Defekt.
    """
    defekt = mu_target - index
    if defekt == 0:
        raise ValueError("index und mu_target duerfen nicht identisch sein.")
    return (index / defekt) * (RR(8) / RR(12))


def local_eabc_index_for_triplet(triplet):
    """Ein lokaler eabc-Index aus den imaginaeren Vier-Quadrate-Anteilen."""
    q_sum = sum(get_imaginary_vector(x).reduced_norm() for x in triplet)
    return RR(q_sum) / RR(12)


def batch_process_kanon(limit=1000, observable=None, mu_target=1836.15267):
    """
    Fuehrt einen statistischen Kanon-Check fuer viele Primzahldrillinge durch.
    """
    if limit < 1:
        raise ValueError("limit muss mindestens 1 sein.")

    observable = RR(observable if observable is not None else "0.13387")
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    triplets, _ = find_prime_clusters(limit_triplets=limit, limit_quads=0)
    beta_distribution = []

    for triplet in triplets:
        xi_local = local_eabc_index_for_triplet(triplet)
        beta_local = (RR(mu_target) - xi_local) / observable
        beta_distribution.append(float(beta_local))

    mean_beta = sum(beta_distribution) / len(beta_distribution)
    variance = sum((value - mean_beta) ** 2 for value in beta_distribution) / len(beta_distribution)
    std_beta = variance ** 0.5

    return {
        "count": len(beta_distribution),
        "mean_beta": mean_beta,
        "std_beta": std_beta,
        "min_beta": min(beta_distribution),
        "max_beta": max(beta_distribution),
        "beta_distribution": beta_distribution,
        "clusters": triplets,
    }


@lru_cache(maxsize=None)
def cached_virtual_complement_energy(triplet, iterations=20, step_size=0.15, gradient_eps=0.05):
    """Cached die virtuelle Komplementenergie fuer wiederholte Analysen desselben Drillings."""
    result = find_virtual_complement(
        triplet,
        iterations=iterations,
        step_size=step_size,
        gradient_eps=gradient_eps,
    )
    return RR(result["virtual_energy"])


def calculate_damped_beta(xi_local, mu_target, observable, cluster):
    """
    Berechnet ein stabilisiertes Beta mit glatter Stress-Daempfung
    und einer gecachten virtuellen Heilungsenergie.
    """
    observable = RR(observable)
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    stress_matrix = calculate_stress_matrix(cluster)
    det_val = RR(abs(stress_matrix.det()))
    damping = RR(1) + log(det_val + 1, 10)

    healing_energy = cached_virtual_complement_energy(tuple(cluster), iterations=20)
    baseline = RR(1) + RR(xi_local)
    healing_tension = healing_energy / baseline

    effective_obs = observable * damping * (RR(1) + healing_tension)
    return (RR(mu_target) - RR(xi_local)) / effective_obs


def batch_process_stabilized(limit=1000, observable=None, mu_target=1836.15267):
    """Fuehrt den Kanon-Check mit stabilisierter Beta-Formel aus."""
    if limit < 1:
        raise ValueError("limit muss mindestens 1 sein.")

    observable = RR(observable if observable is not None else "0.13387")
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    triplets, _ = find_prime_clusters(limit_triplets=limit, limit_quads=0)
    stable_betas = []

    for triplet in triplets:
        xi_local = local_eabc_index_for_triplet(triplet)
        beta_value = calculate_damped_beta(xi_local, mu_target, observable, triplet)
        stable_betas.append(float(beta_value))

    mean_beta = sum(stable_betas) / len(stable_betas)
    variance = sum((value - mean_beta) ** 2 for value in stable_betas) / len(stable_betas)

    return {
        "count": len(stable_betas),
        "new_mean": mean_beta,
        "new_std": variance ** 0.5,
        "min_beta": min(stable_betas),
        "max_beta": max(stable_betas),
        "beta_distribution": stable_betas,
        "clusters": triplets,
    }


def calculate_damped_beta_prime_scaled(xi_local, mu_target, observable, cluster):
    """
    Alternative Beta-Formel mit Primzahl-Skalierung der Heilungsspannung.
    Diese Variante behaelt die vom Benutzer vorgeschlagene Normierung bei.
    """
    observable = RR(observable)
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    stress_matrix = calculate_stress_matrix(cluster)
    det_val = RR(abs(stress_matrix.det()))
    damping = max(RR(1), log(det_val + 1, 10))

    first_prime = RR(cluster[0])
    if first_prime == 0:
        raise ValueError("cluster[0] darf nicht 0 sein.")

    healing_energy = cached_virtual_complement_energy(tuple(cluster), iterations=20)
    healing_tension = abs(healing_energy - first_prime) / first_prime

    effective_obs = observable * damping * (RR(1) + healing_tension)
    return (RR(mu_target) - RR(xi_local)) / effective_obs


def batch_process_stabilized_prime_scaled(limit=1000, observable=None, mu_target=1836.15267):
    """Vergleichsvariante mit primzahlskalierter Heilungsspannung."""
    if limit < 1:
        raise ValueError("limit muss mindestens 1 sein.")

    observable = RR(observable if observable is not None else "0.13387")
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    triplets, _ = find_prime_clusters(limit_triplets=limit, limit_quads=0)
    stable_betas = []

    for triplet in triplets:
        xi_local = local_eabc_index_for_triplet(triplet)
        beta_value = calculate_damped_beta_prime_scaled(
            xi_local,
            mu_target,
            observable,
            triplet,
        )
        stable_betas.append(float(beta_value))

    mean_beta = sum(stable_betas) / len(stable_betas)
    variance = sum((value - mean_beta) ** 2 for value in stable_betas) / len(stable_betas)

    return {
        "count": len(stable_betas),
        "new_mean": mean_beta,
        "new_std": variance ** 0.5,
        "min_beta": min(stable_betas),
        "max_beta": max(stable_betas),
        "beta_distribution": stable_betas,
        "clusters": triplets,
    }


def run_stabilized_kanon_check(limit=1000, observable=None, mu_target=1836.15267):
    """Kompakter Wrapper fuer den stabilisierten statistischen Kanon-Check."""
    stats = batch_process_stabilized(limit=limit, observable=observable, mu_target=mu_target)
    return {
        "count": stats["count"],
        "mean": stats["new_mean"],
        "std": stats["new_std"],
        "min": stats["min_beta"],
        "max": stats["max_beta"],
        "beta_distribution": stats["beta_distribution"],
        "clusters": stats["clusters"],
    }


def run_prime_scaled_kanon_check(limit=1000, observable=None, mu_target=1836.15267):
    """Kompakter Wrapper fuer die primzahlskalierte Vergleichsvariante."""
    stats = batch_process_stabilized_prime_scaled(
        limit=limit,
        observable=observable,
        mu_target=mu_target,
    )
    return {
        "count": stats["count"],
        "mean": stats["new_mean"],
        "std": stats["new_std"],
        "min": stats["min_beta"],
        "max": stats["max_beta"],
        "beta_distribution": stats["beta_distribution"],
        "clusters": stats["clusters"],
    }


def get_chebyshev_weight(p, strength=RR("0.15"), cap=RR("0.25")):
    """
    Milde Chebyshev-Korrektur nach Restklasse mod 4.
    Die Abweichung von 1 wird bewusst begrenzt, damit der Bias nicht dominiert.
    """
    p = RR(p)
    if p <= 2:
        return RR(1)

    delta = min(cap, strength / log(p))
    if ZZ(p) % 4 == 3:
        return RR(1) - delta
    return RR(1) + delta


def calculate_eabc_final_bias_beta(
    xi_local,
    mu_target,
    observable,
    cluster,
    strength=RR("0.15"),
    cap=RR("0.25"),
):
    """Stabilisierte Beta-Formel mit milder Chebyshev-Bias-Korrektur."""
    observable = RR(observable)
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    first_prime = RR(cluster[0])
    if first_prime == 0:
        raise ValueError("cluster[0] darf nicht 0 sein.")

    stress_matrix = calculate_stress_matrix(cluster)
    det_val = RR(abs(stress_matrix.det()))
    damping = RR(1) + log(det_val + 1, 10)

    healing_energy = cached_virtual_complement_energy(tuple(cluster), iterations=30)
    spannung = abs(healing_energy - first_prime) / first_prime

    bias_correction = sum(
        get_chebyshev_weight(p, strength=RR(strength), cap=RR(cap)) for p in cluster
    ) / RR(len(cluster))
    effective_obs = observable * bias_correction * (RR(1) + spannung) * damping
    return (RR(mu_target) - RR(xi_local)) / effective_obs


def run_eabc_final_bias_check(
    limit=1000,
    observable=None,
    mu_target=1836.15267,
    strength=RR("0.15"),
    cap=RR("0.25"),
):
    """Wrapper fuer die finale Chebyshev-Bias-Korrektur."""
    if limit < 1:
        raise ValueError("limit muss mindestens 1 sein.")

    observable = RR(observable if observable is not None else "0.13387")
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    triplets, _ = find_prime_clusters(limit_triplets=limit, limit_quads=0)
    bias_betas = []

    for triplet in triplets:
        xi_local = local_eabc_index_for_triplet(triplet)
        beta_value = calculate_eabc_final_bias_beta(
            xi_local,
            mu_target,
            observable,
            triplet,
            strength=RR(strength),
            cap=RR(cap),
        )
        bias_betas.append(float(beta_value))

    mean_beta = sum(bias_betas) / len(bias_betas)
    variance = sum((value - mean_beta) ** 2 for value in bias_betas) / len(bias_betas)

    return {
        "count": len(bias_betas),
        "mean": mean_beta,
        "std": variance ** 0.5,
        "min": min(bias_betas),
        "max": max(bias_betas),
        "beta_distribution": bias_betas,
        "clusters": triplets,
    }


def compare_chebyshev_profiles(limit=1000, observable=None, mu_target=1836.15267):
    """Vergleicht milde, mittlere und starke Chebyshev-Bias-Profile."""
    profiles = {
        "mild": {"strength": RR("0.10"), "cap": RR("0.20")},
        "medium": {"strength": RR("0.15"), "cap": RR("0.25")},
        "strong": {"strength": RR("0.20"), "cap": RR("0.30")},
    }
    return {
        name: run_eabc_final_bias_check(
            limit=limit,
            observable=observable,
            mu_target=mu_target,
            strength=params["strength"],
            cap=params["cap"],
        )
        for name, params in profiles.items()
    }


def calculate_mid_stable_beta(xi_local, mu_target, observable, cluster):
    """
    Mittelstarke Beta-Berechnung mit Log-Daempfung zur Potenz 1.5.
    Sie bildet den Kompromiss zwischen glatter und ultra-starker Stabilisierung.
    """
    observable = RR(observable)
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    first_prime = RR(cluster[0])
    if first_prime == 0:
        raise ValueError("cluster[0] darf nicht 0 sein.")

    stress_matrix = calculate_stress_matrix(cluster)
    det_val = RR(abs(stress_matrix.det()))
    damping = max(RR(1), log(det_val + 1, 10) ** RR("1.5"))

    healing_energy = cached_virtual_complement_energy(tuple(cluster), iterations=30)
    spannung = abs(healing_energy - first_prime) / first_prime

    effektive_obs = observable * damping * (RR(1) + spannung)
    return (RR(mu_target) - RR(xi_local)) / effektive_obs


def run_mid_stable_kanon_check(limit=1000, observable=None, mu_target=1836.15267):
    """Wrapper fuer die mittlere Beta-Stabilisierung mit Potenz 1.5."""
    if limit < 1:
        raise ValueError("limit muss mindestens 1 sein.")

    observable = RR(observable if observable is not None else "0.13387")
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    triplets, _ = find_prime_clusters(limit_triplets=limit, limit_quads=0)
    mid_betas = []

    for triplet in triplets:
        xi_local = local_eabc_index_for_triplet(triplet)
        beta_value = calculate_mid_stable_beta(xi_local, mu_target, observable, triplet)
        mid_betas.append(float(beta_value))

    mean_beta = sum(mid_betas) / len(mid_betas)
    variance = sum((value - mean_beta) ** 2 for value in mid_betas) / len(mid_betas)

    return {
        "count": len(mid_betas),
        "mean": mean_beta,
        "std": variance ** 0.5,
        "min": min(mid_betas),
        "max": max(mid_betas),
        "beta_distribution": mid_betas,
        "clusters": triplets,
    }


def calculate_resonance_eich_beta(
    xi_local,
    mu_target,
    observable,
    cluster,
    rotation_offset=14.1616,
):
    """
    Finalisierte Eichung: Der Daempfungsexponent ist eine Funktion
    der lokalen Zeta-Amplitude.
    """
    observable = RR(observable)
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    first_prime = RR(cluster[0])
    if first_prime == 0:
        raise ValueError("cluster[0] darf nicht 0 sein.")

    stress_matrix = calculate_stress_matrix(cluster)
    det_val = RR(abs(stress_matrix.det()))

    t_local = RR(integer_to_quaternion(cluster[1]).reduced_norm()) * RR(rotation_offset)
    zeta_amp = RR(abs(zeta(RR(1) / 2 + I * t_local)))

    # Harmonische Cluster mit kleiner Zeta-Amplitude werden weniger hart gedaempft.
    gamma = RR(1) + zeta_amp
    damping = max(RR(1), log(det_val + 1, 10) ** gamma)

    healing_energy = cached_virtual_complement_energy(tuple(cluster), iterations=40)
    spannung = abs(healing_energy - first_prime) / first_prime

    effektive_obs = observable * damping * (RR(1) + spannung)
    return (RR(mu_target) - RR(xi_local)) / effektive_obs


def run_resonance_eich_kanon_check(
    limit=1000,
    observable=None,
    mu_target=1836.15267,
    rotation_offset=14.1616,
):
    """Wrapper fuer die resonanzgesteuerte Zeta-Eichung."""
    if limit < 1:
        raise ValueError("limit muss mindestens 1 sein.")

    observable = RR(observable if observable is not None else "0.13387")
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    triplets, _ = find_prime_clusters(limit_triplets=limit, limit_quads=0)
    resonance_betas = []

    for triplet in triplets:
        xi_local = local_eabc_index_for_triplet(triplet)
        beta_value = calculate_resonance_eich_beta(
            xi_local,
            mu_target,
            observable,
            triplet,
            rotation_offset=rotation_offset,
        )
        resonance_betas.append(float(beta_value))

    mean_beta = sum(resonance_betas) / len(resonance_betas)
    variance = sum((value - mean_beta) ** 2 for value in resonance_betas) / len(resonance_betas)

    return {
        "count": len(resonance_betas),
        "mean": mean_beta,
        "std": variance ** 0.5,
        "min": min(resonance_betas),
        "max": max(resonance_betas),
        "beta_distribution": resonance_betas,
        "clusters": triplets,
    }


def calculate_ultra_stable_beta(xi_local, mu_target, observable, cluster):
    """
    Ultra-stabile Beta-Berechnung mit quadratischer Log-Daempfung.
    Isoliert Singularitaeten im Primzahlraum besonders stark.
    """
    observable = RR(observable)
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    first_prime = RR(cluster[0])
    if first_prime == 0:
        raise ValueError("cluster[0] darf nicht 0 sein.")

    stress_matrix = calculate_stress_matrix(cluster)
    det_val = RR(abs(stress_matrix.det()))

    # Quadratische Daempfung: 10^12 Stress fuehrt naeherungsweise auf Faktor 144.
    damping = max(RR(1), log(det_val + 1, 10) ** 2)

    healing_energy = cached_virtual_complement_energy(tuple(cluster), iterations=40)
    spannung = abs(healing_energy - first_prime) / first_prime

    effektive_obs = observable * damping * (RR(1) + spannung)
    return (RR(mu_target) - RR(xi_local)) / effektive_obs


def run_ultra_stabilized_kanon_check(limit=1000, observable=None, mu_target=1836.15267):
    """Wrapper fuer die ultra-stabilisierte Beta-Auswertung mit Quadrat-Daempfung."""
    if limit < 1:
        raise ValueError("limit muss mindestens 1 sein.")

    observable = RR(observable if observable is not None else "0.13387")
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    triplets, _ = find_prime_clusters(limit_triplets=limit, limit_quads=0)
    ultra_betas = []

    for triplet in triplets:
        xi_local = local_eabc_index_for_triplet(triplet)
        beta_value = calculate_ultra_stable_beta(xi_local, mu_target, observable, triplet)
        ultra_betas.append(float(beta_value))

    mean_beta = sum(ultra_betas) / len(ultra_betas)
    variance = sum((value - mean_beta) ** 2 for value in ultra_betas) / len(ultra_betas)

    return {
        "count": len(ultra_betas),
        "mean": mean_beta,
        "std": variance ** 0.5,
        "min": min(ultra_betas),
        "max": max(ultra_betas),
        "beta_distribution": ultra_betas,
        "clusters": triplets,
    }


def calculate_ultra_chebyshev_beta(xi_local, mu_target, observable, cluster):
    """
    Koppelt quadratische Log-Daempfung mit einer milden Chebyshev-Modulation
    der Heilungsspannung.
    """
    observable = RR(observable)
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    first_prime = RR(cluster[0])
    if first_prime == 0:
        raise ValueError("cluster[0] darf nicht 0 sein.")

    stress_matrix = calculate_stress_matrix(cluster)
    det_val = RR(abs(stress_matrix.det()))
    damping = max(RR(1), log(det_val + 1, 10) ** 2)

    rho = sum(get_chebyshev_weight(p) for p in cluster) / RR(len(cluster))

    healing_energy = cached_virtual_complement_energy(tuple(cluster), iterations=50)
    spannung = (abs(healing_energy - first_prime) / first_prime) * rho

    effektive_obs = observable * damping * (RR(1) + spannung)
    return (RR(mu_target) - RR(xi_local)) / effektive_obs


def run_ultra_chebyshev_check(limit=1000, observable=None, mu_target=1836.15267):
    """Wrapper fuer die ultra-chebyshevische Hybridkopplung."""
    if limit < 1:
        raise ValueError("limit muss mindestens 1 sein.")

    observable = RR(observable if observable is not None else "0.13387")
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    triplets, _ = find_prime_clusters(limit_triplets=limit, limit_quads=0)
    hybrid_betas = []

    for triplet in triplets:
        xi_local = local_eabc_index_for_triplet(triplet)
        beta_value = calculate_ultra_chebyshev_beta(xi_local, mu_target, observable, triplet)
        hybrid_betas.append(float(beta_value))

    mean_beta = sum(hybrid_betas) / len(hybrid_betas)
    variance = sum((value - mean_beta) ** 2 for value in hybrid_betas) / len(hybrid_betas)

    return {
        "count": len(hybrid_betas),
        "mean": mean_beta,
        "std": variance ** 0.5,
        "min": min(hybrid_betas),
        "max": max(hybrid_betas),
        "beta_distribution": hybrid_betas,
        "clusters": triplets,
    }


def calculate_final_torsion_beta(xi_local, mu_target, observable, cluster):
    """
    Die vollstaendige eabc-Gleichung inkl. oktonionischer Torsion.
    """
    observable = RR(observable)
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    first_prime = RR(cluster[0])
    if first_prime == 0:
        raise ValueError("cluster[0] darf nicht 0 sein.")

    stress_matrix = calculate_stress_matrix(cluster)
    det_val = RR(abs(stress_matrix.det()))
    damping = max(RR(1), log(det_val + 1, 10) ** 2)

    tau = calculate_octonion_torsion(cluster)
    rho = sum(get_chebyshev_weight(p) for p in cluster) / RR(len(cluster))

    healing_energy = cached_virtual_complement_energy(tuple(cluster), iterations=50)
    spannung = abs(healing_energy - first_prime) / first_prime

    effektive_obs = observable * damping * (RR(1) + spannung * rho + tau)
    return (RR(mu_target) - RR(xi_local)) / effektive_obs


def calculate_full_torsion_beta(xi_local, mu_target, observable, cluster):
    """Rueckwaertskompatibler Alias fuer die finale Volltorsions-Gleichung."""
    return calculate_final_torsion_beta(xi_local, mu_target, observable, cluster)


def run_full_torsion_check(limit=1000, observable=None, mu_target=1836.15267):
    """Wrapper fuer das voll gekoppelte Torsionsmodell."""
    if limit < 1:
        raise ValueError("limit muss mindestens 1 sein.")

    observable = RR(observable if observable is not None else "0.13387")
    if observable == 0:
        raise ValueError("observable darf nicht 0 sein.")

    triplets, _ = find_prime_clusters(limit_triplets=limit, limit_quads=0)
    torsion_betas = []

    for triplet in triplets:
        xi_local = local_eabc_index_for_triplet(triplet)
        beta_value = calculate_final_torsion_beta(xi_local, mu_target, observable, triplet)
        torsion_betas.append(float(beta_value))

    mean_beta = sum(torsion_betas) / len(torsion_betas)
    variance = sum((value - mean_beta) ** 2 for value in torsion_betas) / len(torsion_betas)

    return {
        "count": len(torsion_betas),
        "mean": mean_beta,
        "std": variance ** 0.5,
        "min": min(torsion_betas),
        "max": max(torsion_betas),
        "beta_distribution": torsion_betas,
        "clusters": triplets,
    }


def compare_beta_models(limit=1000, observable=None, mu_target=1836.15267):
    """Vergleicht die stabilisierten Beta-Modelle anhand derselben Cluster."""
    return {
        "stabilized": run_stabilized_kanon_check(
            limit=limit,
            observable=observable,
            mu_target=mu_target,
        ),
        "prime_scaled": run_prime_scaled_kanon_check(
            limit=limit,
            observable=observable,
            mu_target=mu_target,
        ),
        "chebyshev_bias": run_eabc_final_bias_check(
            limit=limit,
            observable=observable,
            mu_target=mu_target,
        ),
        "mid": run_mid_stable_kanon_check(
            limit=limit,
            observable=observable,
            mu_target=mu_target,
        ),
        "resonance_eich": run_resonance_eich_kanon_check(
            limit=limit,
            observable=observable,
            mu_target=mu_target,
        ),
        "ultra": run_ultra_stabilized_kanon_check(
            limit=limit,
            observable=observable,
            mu_target=mu_target,
        ),
        "ultra_chebyshev": run_ultra_chebyshev_check(
            limit=limit,
            observable=observable,
            mu_target=mu_target,
        ),
        "full_torsion": run_full_torsion_check(
            limit=limit,
            observable=observable,
            mu_target=mu_target,
        ),
    }


def analyze_outlier(distribution, clusters, target_beta=10.9446):
    """Findet den Drilling mit der groessten Abweichung vom Ziel-Beta."""
    if not distribution or not clusters or len(distribution) != len(clusters):
        raise ValueError("distribution und clusters muessen gleich lang und nicht leer sein.")

    diffs = [abs(beta - target_beta) for beta in distribution]
    outlier_idx = diffs.index(max(diffs))
    return {
        "cluster": clusters[outlier_idx],
        "beta": distribution[outlier_idx],
        "delta": diffs[outlier_idx],
        "index": outlier_idx,
    }


def calculate_stress_matrix(cluster):
    """Berechnet eine 3x3-Stress-Matrix fuer einen Primzahldrilling."""
    pts = [QR(state.coefficient_tuple()) for state in cluster_to_quaternions(cluster)]
    stress_matrix = matrix(RR, 3, 3)

    for row in range(3):
        for col in range(3):
            product_imag = quaternion_imaginary_part(pts[row] * pts[col])
            stress_matrix[row, col] = sqrt(product_imag.reduced_norm())

    return stress_matrix


def build_annotated_matrix_plot(data_matrix, cmap, title):
    """Erzeugt eine beschriftete Heatmap mit Matrixwerten und Achsenlabels."""
    nrows = data_matrix.nrows()
    ncols = data_matrix.ncols()
    plot = matrix_plot(data_matrix, cmap=cmap, colorbar=True)

    for row in range(nrows):
        for col in range(ncols):
            value = float(data_matrix[row, col])
            x = col + RR(0.5)
            y = nrows - row - RR(0.5)
            plot += text(f"{value:.1f}", (x, y), color="white", fontsize=10)

    for col in range(ncols):
        x = col + RR(0.5)
        plot += text(f"C{col + 1}", (x, nrows + RR(0.18)), color="black", fontsize=10)

    for row in range(nrows):
        y = nrows - row - RR(0.5)
        plot += text(f"R{row + 1}", (RR(-0.18), y), color="black", fontsize=10)

    plot += text(title, (ncols / RR(2), nrows + RR(0.42)), color="black", fontsize=12)
    return plot


def plot_stress_comparison(kanon_matrix, outlier_matrix, output_path):
    """Speichert zwei beschriftete Stress-Matrizen als nebeneinander liegende Heatmaps."""
    p1 = build_annotated_matrix_plot(kanon_matrix, "viridis", "Idealer Kanon")
    p2 = build_annotated_matrix_plot(outlier_matrix, "inferno", "Instabiler Ausreisser")
    combined = graphics_array([p1, p2])
    combined.save(
        output_path,
        figsize=[12, 5],
        axes=False,
        title="Stress-Matrix Vergleich: Kanon vs. Ausreisser",
    )
    return output_path


def ptolemy_defect(points):
    """Gibt nur den Ptolemaeus-Defekt eines Quaternionen-Vierlings zurueck."""
    return ptolemy_check(*points)[2]


def simulate_resilience(
    quadruplet,
    iterations=50,
    step_size=0.1,
    perturbation=None,
    gradient_eps=0.05,
):
    """Simuliert die Rueckkehr zum Ptolemaeus-Gleichgewicht nach einer Stoerung."""
    states = [QR(state.coefficient_tuple()) for state in cluster_to_quaternions(quadruplet)]
    perturbation = perturbation or QR([0, RR(1) / 2, RR(1) / 2, RR(1) / 2])
    states[2] += perturbation

    direction_x = QR([0, 1, 0, 0])
    direction_y = QR([0, 0, 1, 0])
    direction_z = QR([0, 0, 0, 1])
    current_pts = list(states)
    current_step = RR(step_size)
    eps = RR(gradient_eps)
    history = []

    for _ in range(iterations):
        err = ptolemy_defect(current_pts)
        history.append(float(err))

        gradient_components = []
        for direction in (direction_x, direction_y, direction_z):
            plus_pts = list(current_pts)
            minus_pts = list(current_pts)
            delta = QR([0, *(eps * coeff for coeff in direction.coefficient_tuple()[1:])])
            plus_pts[2] += delta
            minus_pts[2] -= delta
            plus_err = ptolemy_defect(plus_pts)
            minus_err = ptolemy_defect(minus_pts)
            gradient_components.append((plus_err - minus_err) / (2 * eps))

        candidate_pts = list(current_pts)
        candidate_pts[2] -= QR([0, *(current_step * component for component in gradient_components)])
        candidate_err = ptolemy_defect(candidate_pts)

        if candidate_err <= err:
            current_pts = candidate_pts
            current_step = min(RR(step_size), current_step * RR(1.05))
        else:
            current_step *= RR(0.5)

    final_error = ptolemy_defect(current_pts)
    return {
        "history": history,
        "start_error": history[0],
        "end_error": float(final_error),
        "final_states": tuple(current_pts),
    }


def find_virtual_complement(
    triplet,
    iterations=80,
    step_size=0.15,
    gradient_eps=0.05,
):
    """
    Findet zu einem Primzahldrilling ein virtuelles viertes Quaternion D,
    das den Ptolemaeus-Defekt moeglichst klein macht.
    """
    A, B, C = [QR(state.coefficient_tuple()) for state in cluster_to_quaternions(triplet)]
    current_D = A - B + C

    directions = (
        QR([1, 0, 0, 0]),
        QR([0, 1, 0, 0]),
        QR([0, 0, 1, 0]),
        QR([0, 0, 0, 1]),
    )
    current_step = RR(step_size)
    eps = RR(gradient_eps)
    history = []

    for _ in range(iterations):
        current_error = ptolemy_defect((A, B, C, current_D))
        history.append(float(current_error))

        gradient_components = []
        for direction in directions:
            plus_D = current_D + eps * direction
            minus_D = current_D - eps * direction
            plus_error = ptolemy_defect((A, B, C, plus_D))
            minus_error = ptolemy_defect((A, B, C, minus_D))
            gradient_components.append((plus_error - minus_error) / (2 * eps))

        candidate_D = current_D - QR([current_step * component for component in gradient_components])
        candidate_error = ptolemy_defect((A, B, C, candidate_D))

        if candidate_error <= current_error:
            current_D = candidate_D
            current_step = min(RR(step_size), current_step * RR(1.05))
        else:
            current_step *= RR(0.5)

    final_error = ptolemy_defect((A, B, C, current_D))
    return {
        "triplet": triplet,
        "A": A,
        "B": B,
        "C": C,
        "D_virt": current_D,
        "virtual_energy": current_D.reduced_norm(),
        "start_error": history[0],
        "end_error": float(final_error),
        "history": history,
    }


def analysiere_vierling(vierling, start=14.0, ende=30.0, schritte=100):
    """Berechnet Ptolemaeus-Defekt und Zeta-Minimum fuer einen Vierling."""
    zustaende = cluster_to_quaternions(vierling)
    p_sum, d_prod, diff = ptolemy_check(*zustaende)
    optimal_e, min_res, _ = scan_energieniveaus(start, ende, schritte, zustaende)
    return {
        "vierling": vierling,
        "zustaende": zustaende,
        "ptolemy_sum": p_sum,
        "diagonal_prod": d_prod,
        "ptolemy_diff": diff,
        "optimal_e": optimal_e,
        "min_res": min_res,
        "cycle_work": calculate_cycle_work(vierling),
    }


def ranglisten_fuer_vierlinge(vierlinge, top_n=3):
    """Erzeugt Ranglisten nach Ptolemaeus-Defekt und Zeta-Resonanz."""
    analysen = [analysiere_vierling(vierling) for vierling in vierlinge]
    nach_ptolemaeus = sorted(analysen, key=lambda item: (item["ptolemy_diff"], item["vierling"]))
    nach_resonanz = sorted(analysen, key=lambda item: (item["min_res"], item["vierling"]))
    return nach_ptolemaeus[:top_n], nach_resonanz[:top_n]


def drucke_rangliste(titel, eintraege, schluessel):
    """Gibt eine kompakte Rangliste fuer Vierlinge aus."""
    print(titel)
    for index, eintrag in enumerate(eintraege, start=1):
        print(
            f"{index}. {eintrag['vierling']} | "
            f"{schluessel}: {eintrag[schluessel].n()} | "
            f"E*: {eintrag['optimal_e'].n()}"
        )


def drucke_abschnitt(titel):
    """Formatiert einen gut sichtbaren Abschnittstitel fuer die Konsolenausgabe."""
    print(f"### {titel} ###")


def main():
    A, B, C, D = ZUSTAENDE
    p_sum, d_prod, error = ptolemy_check(A, B, C, D)
    optimal_e, min_res, _ = scan_energieniveaus()
    triplets, quads = find_prime_clusters()
    erster_drilling = triplets[0]
    erster_vierling = quads[0]
    beste_ptolemaeus, beste_resonanz = ranglisten_fuer_vierlinge(quads, top_n=3)
    fav_ptolemy = beste_ptolemaeus[0]
    fav_zeta = beste_resonanz[0]
    resilienz_ptolemy = simulate_resilience(fav_ptolemy["vierling"])
    resilienz_zeta = simulate_resilience(fav_zeta["vierling"])
    starr_8 = beste_ptolemaeus[0]["zustaende"] + beste_ptolemaeus[1]["zustaende"]
    frei_4 = fav_zeta["zustaende"]
    konfiguration_8plus4 = analysiere_8plus4_konfiguration(starr_8, frei_4)
    virtuelles_komplement = find_virtual_complement(triplets[-1])
    natuerlicher_vierling_phase = calculate_phase_shift(fav_ptolemy["zustaende"])
    geheilter_drilling_phase = calculate_phase_shift(
        tuple(QR(state.coefficient_tuple()) for state in cluster_to_quaternions(triplets[-1]))
        + (virtuelles_komplement["D_virt"],)
    )

    drucke_abschnitt("Grundlagen des Zirkels")
    print(f"Ptolemaeus-Summe (ac + bd):   {p_sum.n()}")
    print(f"Diagonal-Produkt (ef):        {d_prod.n()}")
    print(f"Abweichung (Entropie-Mass):   {error.n()}")
    print(f"Optimales Energieniveau:      {optimal_e.n()}")
    print(f"Minimale Zeta-Amplitude:      {min_res.n()}")
    print()

    drucke_abschnitt("Kanon-Konfiguration")
    print(f"Erster Drillings-Kanon: {erster_drilling} -> eabc-Resonanz aktiv")
    print(f"100. Drillings-Kanon:   {triplets[-1]}")
    print(
        f"Vier-Quadrate-Darstellung von {erster_drilling[0]}: "
        f"{quat_representation(erster_drilling[0])}"
    )
    print("-" * 40)
    print(f"Erster Vierlings-Kanon: {erster_vierling} -> Ptolemaeus-stabil")
    print(f"50. Vierlings-Kanon:    {quads[-1]}")
    print(
        f"Vier-Quadrate-Darstellung von {erster_vierling[0]}: "
        f"{quat_representation(erster_vierling[0])}"
    )
    print()

    drucke_abschnitt("Ranglisten und Arbeitsflaechen")
    drucke_rangliste("Top-3 Vierlinge nach kleinstem Ptolemaeus-Defekt:", beste_ptolemaeus, "ptolemy_diff")
    print()
    drucke_rangliste("Top-3 Vierlinge nach staerkster Zeta-Resonanz:", beste_resonanz, "min_res")
    print()
    print(f"Arbeitsflaeche (Ptolemaeus-Favorit): {fav_ptolemy['cycle_work'].n()}")
    print(f"Arbeitsflaeche (Zeta-Favorit):     {fav_zeta['cycle_work'].n()}")
    print()

    drucke_abschnitt("Resilienz und Heilung")
    print(
        "Resilienz-Check (Ptolemaeus-Favorit): "
        f"Start-Fehler {resilienz_ptolemy['start_error']:.4f}, "
        f"End-Fehler {resilienz_ptolemy['end_error']:.4f}"
    )
    print(
        "Resilienz-Check (Zeta-Favorit):      "
        f"Start-Fehler {resilienz_zeta['start_error']:.4f}, "
        f"End-Fehler {resilienz_zeta['end_error']:.4f}"
    )
    print()
    print("8+4 Energiedoku-Konfiguration:")
    print(f"Rahmen-Energie: {konfiguration_8plus4['rahmen_energie'].n()}")
    print(f"Resonanzsumme:  {konfiguration_8plus4['resonanz'].n()}")
    print(f"Effizienz-Index: {konfiguration_8plus4['effizienz_index'].n()}")
    print()
    print(f"Virtuelles Komplement fuer Drilling {virtuelles_komplement['triplet']}:")
    print(f"Virtuelle Quaternione: {virtuelles_komplement['D_virt']}")
    print(
        "Benötigtes Energieniveau (virtuelle Norm): "
        f"{float(virtuelles_komplement['virtual_energy']):.2f}"
    )
    print(
        "Ptolemaeus-Defekt: "
        f"{virtuelles_komplement['start_error']:.4f} -> "
        f"{virtuelles_komplement['end_error']:.4f}"
    )
    print(f"Phasenverschiebung Natürlicher Vierling: {float(natuerlicher_vierling_phase):.6f} rad")
    print(f"Phasenverschiebung Geheilter Drilling:   {float(geheilter_drilling_phase):.6f} rad")
    print()
    phys_mu = RR("1836.15267")
    eabc_index = konfiguration_8plus4["effizienz_index"]
    diff = phys_mu - eabc_index
    phase_corr = abs(geheilter_drilling_phase) / (2 * pi)
    zeta_corr = RR(1) / konfiguration_8plus4["resonanz"]
    multiplikativ_index = eabc_index * (1 + phase_corr / 100)
    additiv_phase_index = eabc_index + phase_corr

    if phase_corr == 0:
        alpha_phase = RR(0)
    else:
        alpha_phase = diff / phase_corr
    kalibrierter_phase_index = eabc_index + alpha_phase * phase_corr

    combined_corr = phase_corr * zeta_corr
    if combined_corr == 0:
        beta_combined = RR(0)
    else:
        beta_combined = diff / combined_corr
    kalibrierter_combined_index = eabc_index + beta_combined * combined_corr

    drucke_abschnitt("Korrekturmodelle der eabc-Effizienz")
    print(f"Original Index: {eabc_index.n()}")
    print(f"Zielwert (Proton/Elektron): {phys_mu.n()}")
    print(f"Urspruengliche Differenz: {diff.n()}")
    print(f"Multiplikativer Phasenterm: {multiplikativ_index.n()}")
    print(f"Restdifferenz multiplikativ: {(phys_mu - multiplikativ_index).n()}")
    print(f"Additiver Phasenterm (unkalibriert): {additiv_phase_index.n()}")
    print(f"Restdifferenz additiv: {(phys_mu - additiv_phase_index).n()}")
    print(f"Kalibrierter Phasenkoeffizient alpha: {alpha_phase.n()}")
    print(f"Kalibriertes Phasenmodell: {kalibrierter_phase_index.n()}")
    print(f"Kombinierte Observable Phase/Zeta: {combined_corr.n()}")
    print(f"Kalibrierter Kombinationskoeffizient beta: {beta_combined.n()}")
    print(f"Kalibriertes Kombinationsmodell: {kalibrierter_combined_index.n()}")
    print()
    batch_stats = batch_process_kanon(limit=1000, observable=combined_corr, mu_target=phys_mu)
    drucke_abschnitt("Statistischer Kanon-Check")
    print(f"Anzahl Cluster: {batch_stats['count']}")
    print(f"Mittelwert Beta: {batch_stats['mean_beta']:.4f}")
    print(f"Standardabweichung: {batch_stats['std_beta']:.4f}")
    print(f"Minimum Beta: {batch_stats['min_beta']:.4f}")
    print(f"Maximum Beta: {batch_stats['max_beta']:.4f}")
    print()
    outlier_info = analyze_outlier(batch_stats["beta_distribution"], batch_stats["clusters"])
    kanon_cluster = fav_ptolemy["vierling"][:3]
    kanon_stress_matrix = calculate_stress_matrix(kanon_cluster)
    stress_matrix = calculate_stress_matrix(outlier_info["cluster"])
    stress_plot_path = plot_stress_comparison(
        kanon_stress_matrix,
        stress_matrix,
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "stress_comparison.png"),
    )
    drucke_abschnitt("Kollaps-Analyse des instabilsten Drillings")
    print(f"Idealer Kanon (Stress-Referenz): {kanon_cluster}")
    print(f"Cluster: {outlier_info['cluster']}")
    print(
        f"Abweichendes Beta: {outlier_info['beta']:.4f} "
        f"(Delta: {outlier_info['delta']:.4f})"
    )
    print(f"Determinante der Stress-Matrix: {stress_matrix.det().n()}")
    print(f"Heatmap gespeichert unter: {stress_plot_path}")
    print()
    oktonion_kanon = prime_octonion(fav_ptolemy["vierling"])
    oktonion_zeta = prime_octonion(fav_zeta["vierling"])
    oktonion_basis = prime_octonion(erster_vierling)
    oktonion_asso = octonion_associator(oktonion_kanon, oktonion_zeta, oktonion_basis)
    drucke_abschnitt("Oktonionische Erweiterung")
    print(f"Kanon-Oktonion: {oktonion_kanon}")
    print(f"Zeta-Oktonion:  {oktonion_zeta}")
    print(f"Assoziator-Norm: {oktonion_asso.norm().n()}")
    torsion_kanon_zeta_basis = get_octonion_torsion(
        fav_ptolemy["vierling"],
        fav_zeta["vierling"],
        erster_vierling,
    )
    torsion_kanon_heilung_hybrid = get_octonion_torsion(
        kanon_cluster,
        virtuelles_komplement["triplet"],
        outlier_info["cluster"],
    )
    print("Oktonionische Torsionen:")
    print(f"Kanon/Zeta/Basis:     {torsion_kanon_zeta_basis.n()}")
    print(f"Kanon/Heilung/Hybrid: {torsion_kanon_heilung_hybrid.n()}")
    print()
    alpha_inv = RR("137.035999")
    kopplung = diff / alpha_inv

    drucke_abschnitt("Alpha-Kopplung und Vorhersagen")
    print("Alpha-Kopplung im eabc-Modell:")
    print(f"Inverses Feinstrukturmass: {alpha_inv.n()}")
    print(f"Benoetigter Korrekturfaktor pro Alpha-Einheit: {kopplung.n()}")
    print()
    alpha_inv_pred = (eabc_index / diff) * (log(phys_mu) / phys_mu) * RR(137)
    alpha_inv_geo = (eabc_index / RR(12)) - (diff * pi)
    print("Alpha-Vorhersage im eabc-Modell:")
    print(f"Bamberger Index (Struktur): {eabc_index.n()}")
    print(f"Vorhergesagtes Alpha^-1 (kalibriert): {alpha_inv_pred.n()}")
    print(f"Vorhergesagtes Alpha^-1 (geometrisch): {alpha_inv_geo.n()}")
    print(f"Referenzwert (CODATA): {alpha_inv.n()}")
    print(f"Abweichung geometrisch: {(alpha_inv - alpha_inv_geo).n()}")
    print()
    alpha_calc = validate_alpha_coupling(eabc_index, phys_mu)
    print("Bamberger Eich-Modul:")
    print(f"Theoretisches Alpha^-1: {alpha_calc.n()}")
    print(f"Abweichung zum CODATA-Wert: {(alpha_inv - alpha_calc).n()}")
    print()

    drucke_abschnitt("Stabilisierter statistischer Kanon-Check")
    stats = run_stabilized_kanon_check(limit=1000, observable=combined_corr, mu_target=phys_mu)
    print(f"Anzahl Cluster: {stats['count']}")
    print(f"Stabilisierter Mittelwert Beta: {stats['mean']:.4f}")
    print(f"Stabilisierte Standardabweichung: {stats['std']:.4f}")
    print(f"Beta-Bereich: [{stats['min']:.4f} bis {stats['max']:.4f}]")
    print()

    phase_healing_factor = RR(12)
    alpha_inv_final = alpha_inv_geo + (abs(geheilter_drilling_phase) / pi) * phase_healing_factor
    print("Optimierte Alpha-Vorhersage (Phase + Geometrie):")
    print(f"Alpha^-1 (eabc-stabilisiert): {RR(alpha_inv_final).n()}")
    print(f"Differenz zum CODATA-Wert:    {RR(alpha_inv - alpha_inv_final).n()}")
    print()

    drucke_abschnitt("Chebyshev-Bias-Korrektur")
    chebyshev_stats = run_eabc_final_bias_check(
        limit=1000,
        observable=combined_corr,
        mu_target=phys_mu,
    )
    print(f"Anzahl Cluster: {chebyshev_stats['count']}")
    print(f"Mittelwert Beta (Chebyshev): {chebyshev_stats['mean']:.4f}")
    print(f"Standardabweichung:          {chebyshev_stats['std']:.4f}")
    print(f"Beta-Bereich:                [{chebyshev_stats['min']:.4f} bis {chebyshev_stats['max']:.4f}]")
    print()

    drucke_abschnitt("Chebyshev-Profile im Vergleich")
    chebyshev_profiles = compare_chebyshev_profiles(
        limit=1000,
        observable=combined_corr,
        mu_target=phys_mu,
    )
    for profile_name, profile_stats in (
        ("Mild", chebyshev_profiles["mild"]),
        ("Medium", chebyshev_profiles["medium"]),
        ("Strong", chebyshev_profiles["strong"]),
    ):
        print(
            f"{profile_name:8s} | "
            f"Mittel {profile_stats['mean']:.4f} | "
            f"Std {profile_stats['std']:.4f} | "
            f"Bereich [{profile_stats['min']:.4f}, {profile_stats['max']:.4f}]"
        )
    best_chebyshev_profile = min(
        (
            ("Mild", chebyshev_profiles["mild"]),
            ("Medium", chebyshev_profiles["medium"]),
            ("Strong", chebyshev_profiles["strong"]),
        ),
        key=lambda entry: (
            entry[1]["std"],
            entry[1]["max"] - entry[1]["min"],
            abs(entry[1]["mean"]),
        ),
    )
    best_profile_width = best_chebyshev_profile[1]["max"] - best_chebyshev_profile[1]["min"]
    print("-" * 40)
    print(
        "Bestes Chebyshev-Profil: "
        f"{best_chebyshev_profile[0]} "
        f"(Std {best_chebyshev_profile[1]['std']:.4f}, "
        f"Breite {best_profile_width:.4f})"
    )
    print()

    drucke_abschnitt("Mittlere Beta-Stabilisierung")
    mid_stats = run_mid_stable_kanon_check(
        limit=1000,
        observable=combined_corr,
        mu_target=phys_mu,
    )
    print(f"Anzahl Cluster: {mid_stats['count']}")
    print(f"Mittelwert Beta (Potenz 1.5): {mid_stats['mean']:.4f}")
    print(f"Standardabweichung:           {mid_stats['std']:.4f}")
    print(f"Beta-Bereich:                 [{mid_stats['min']:.4f} bis {mid_stats['max']:.4f}]")
    print()

    drucke_abschnitt("Resonanz-Eichung mit Zeta-Katalysator")
    resonance_stats = run_resonance_eich_kanon_check(
        limit=1000,
        observable=combined_corr,
        mu_target=phys_mu,
    )
    print(f"Anzahl Cluster: {resonance_stats['count']}")
    print(f"Mittelwert Beta (resonant):   {resonance_stats['mean']:.4f}")
    print(f"Standardabweichung:           {resonance_stats['std']:.4f}")
    print(f"Beta-Bereich:                 [{resonance_stats['min']:.4f} bis {resonance_stats['max']:.4f}]")
    print()

    drucke_abschnitt("Finale Riemann-Resonanz-Eichung")
    r_mean = resonance_stats["mean"]
    r_std = resonance_stats["std"]
    print(f"Mittelwert Beta (resonant): {r_mean:.4f}")
    print(f"Standardabweichung:         {r_std:.4f}")

    okto_torsion_const = RR("34.518167")
    resonance_beta_reference = RR("10.9446")
    codata_alpha_reference = RR("137.036")
    alpha_inv_raw = alpha_inv_geo + (abs(geheilter_drilling_phase) / pi) * phase_healing_factor
    alpha_inv_resonance = alpha_inv_raw + (
        RR(r_mean) / resonance_beta_reference
    ) * (okto_torsion_const / codata_alpha_reference)

    coherence = max(0.0, (1 - r_std / abs(r_mean)) * 100) if r_mean != 0 else 0.0

    print("-" * 40)
    print(f"Alpha^-1 (eabc-Resonanz):   {RR(alpha_inv_resonance).n()}")
    print(f"Praezision (Delta zu CODATA): {RR(alpha_inv - alpha_inv_resonance).n()}")
    print(f"Systemische Kohaerenz:       {coherence:.2f} %")
    print()

    drucke_abschnitt("Ultra-stabilisierter Kanon-Check")
    ultra_stats = run_ultra_stabilized_kanon_check(
        limit=1000,
        observable=combined_corr,
        mu_target=phys_mu,
    )
    print(f"Anzahl Cluster: {ultra_stats['count']}")
    print(f"Mittelwert Beta (quadratisch): {ultra_stats['mean']:.4f}")
    print(f"Standardabweichung:            {ultra_stats['std']:.4f}")
    print(f"Beta-Bereich:                  [{ultra_stats['min']:.4f} bis {ultra_stats['max']:.4f}]")
    print()

    ultra_phase_factor = RR(12)
    ultra_beta_reference = RR("10.9446")
    codata_alpha_reference = RR("137.036")
    okto_torsion = RR("34.518")
    alpha_inv_raw = alpha_inv_geo + (abs(geheilter_drilling_phase) / pi) * ultra_phase_factor
    alpha_inv_ultra = alpha_inv_raw + (
        RR(ultra_stats["mean"]) / ultra_beta_reference
    ) * (okto_torsion / codata_alpha_reference)

    print("Alpha-Vorhersage inkl. oktonionischer Torsion:")
    print(f"Alpha^-1 (eabc-ultra):      {RR(alpha_inv_ultra).n()}")
    print(f"Rest-Differenz zu CODATA:   {RR(alpha_inv - alpha_inv_ultra).n()}")
    print()

    drucke_abschnitt("Bamberger Vollendung: Ultra-Chebyshev-Kopplung")
    hybrid_stats = run_ultra_chebyshev_check(
        limit=1000,
        observable=combined_corr,
        mu_target=phys_mu,
    )
    h_mean = hybrid_stats["mean"]
    h_std = hybrid_stats["std"]
    hybrid_beta_reference = RR("10.9446")
    okto_torsion_const = RR("34.518167")
    codata_alpha_reference = RR("137.036")
    alpha_inv_eabc = alpha_inv_raw + (
        RR(h_mean) / hybrid_beta_reference
    ) * (okto_torsion_const / codata_alpha_reference)
    print(f"Hybrid-Mittelwert (Beta):    {h_mean:.6f}")
    print(f"Hybrid-Standardabweichung:   {h_std:.4f}")
    print("-" * 45)
    print(f"Alpha^-1 (eabc-Hybrid):      {RR(alpha_inv_eabc).n()}")
    print(f"CODATA-Referenz:             {alpha_inv.n()}")
    print(f"Finale Praezision (Delta):   {RR(alpha_inv - alpha_inv_eabc).n()}")
    if abs(RR(alpha_inv - alpha_inv_eabc)) < RR("0.1"):
        print("STATUS: #Energiedoku - Die Feinstrukturkonstante ist geometrisch geeicht.")
    print()

    drucke_abschnitt("Vollmodell: Torsion, Bias und Heilung")
    full_torsion_stats = run_full_torsion_check(
        limit=1000,
        observable=combined_corr,
        mu_target=phys_mu,
    )
    full_mean = full_torsion_stats["mean"]
    full_std = full_torsion_stats["std"]
    full_beta_reference = RR("10.9446")
    full_torsion_ref = RR("34.518167")
    codata_alpha_reference = RR("137.036")
    alpha_inv_full = alpha_inv_raw + (
        RR(full_mean) / full_beta_reference
    ) * (full_torsion_ref / codata_alpha_reference)
    print(f"Vollmodell-Mittelwert (Beta): {full_mean:.6f}")
    print(f"Vollmodell-Standardabw.:     {full_std:.4f}")
    print("-" * 45)
    print(f"Alpha^-1 (eabc-Vollmodell):  {RR(alpha_inv_full).n()}")
    print(f"CODATA-Referenzwert:         {alpha_inv.n()}")
    print(f"Finale Praezision (Delta):   {RR(alpha_inv - alpha_inv_full).n()}")
    print()

    drucke_abschnitt("Vergleich der Alpha-Endwerte")
    alpha_comparison = (
        ("Stabilisiert", RR(alpha_inv_final), RR(alpha_inv - alpha_inv_final)),
        ("Resonant", RR(alpha_inv_resonance), RR(alpha_inv - alpha_inv_resonance)),
        ("Ultra", RR(alpha_inv_ultra), RR(alpha_inv - alpha_inv_ultra)),
        ("Hybrid", RR(alpha_inv_eabc), RR(alpha_inv - alpha_inv_eabc)),
        ("Vollmodell", RR(alpha_inv_full), RR(alpha_inv - alpha_inv_full)),
    )
    for model_name, alpha_value, alpha_delta in alpha_comparison:
        print(
            f"{model_name:12s} | "
            f"Alpha^-1 {alpha_value.n()} | "
            f"Delta {alpha_delta.n()}"
        )
    best_alpha_model = min(alpha_comparison, key=lambda entry: abs(entry[2]))
    print("-" * 40)
    print(
        "Bestes Alpha-Modell: "
        f"{best_alpha_model[0]} "
        f"(Delta {RR(best_alpha_model[2]).n()})"
    )
    print()

    drucke_abschnitt("Modellvergleich der Beta-Stabilisierung")
    model_comparison = compare_beta_models(
        limit=1000,
        observable=combined_corr,
        mu_target=phys_mu,
    )
    for model_name, model_stats in (
        ("Stabilized", model_comparison["stabilized"]),
        ("Prime-Scaled", model_comparison["prime_scaled"]),
        ("Chebyshev", model_comparison["chebyshev_bias"]),
        ("Mid-1.5", model_comparison["mid"]),
        ("Resonance", model_comparison["resonance_eich"]),
        ("Ultra", model_comparison["ultra"]),
        ("Ultra-Cheb", model_comparison["ultra_chebyshev"]),
        ("Full-Tors.", model_comparison["full_torsion"]),
    ):
        print(
            f"{model_name:12s} | "
            f"Mittel {model_stats['mean']:.4f} | "
            f"Std {model_stats['std']:.4f} | "
            f"Bereich [{model_stats['min']:.4f}, {model_stats['max']:.4f}]"
        )
    print()

    drucke_abschnitt("Rangfolge der Beta-Modelle")
    ranked_models = sorted(
        (
            (
                model_name,
                model_stats,
                model_stats["max"] - model_stats["min"],
            )
            for model_name, model_stats in (
                ("Stabilized", model_comparison["stabilized"]),
                ("Prime-Scaled", model_comparison["prime_scaled"]),
                ("Chebyshev", model_comparison["chebyshev_bias"]),
                ("Mid-1.5", model_comparison["mid"]),
                ("Resonance", model_comparison["resonance_eich"]),
                ("Ultra", model_comparison["ultra"]),
                ("Ultra-Cheb", model_comparison["ultra_chebyshev"]),
                ("Full-Tors.", model_comparison["full_torsion"]),
            )
        ),
        key=lambda entry: (entry[1]["std"], entry[2], abs(entry[1]["mean"])),
    )
    for rank, (model_name, model_stats, width) in enumerate(ranked_models, start=1):
        print(
            f"{rank}. {model_name:12s} | "
            f"Std {model_stats['std']:.4f} | "
            f"Breite {width:.4f} | "
            f"Mittel {model_stats['mean']:.4f}"
        )
    best_beta_model = ranked_models[0]
    print("-" * 40)
    print(
        "Bestes Beta-Modell: "
        f"{best_beta_model[0]} "
        f"(Std {best_beta_model[1]['std']:.4f}, "
        f"Breite {best_beta_model[2]:.4f})"
    )
    print()

    drucke_abschnitt("Zeta-Perkolations-Proxy")
    zeta_impact = simulate_zeta_impact(
        z=96,
        N_samples=40,
        tower_width=48,
        M=500,
        seed=42,
    )
    print(f"Zeta-Gewichtung (EABC): {zeta_impact['weights']}")
    print(f"Gewichtetes Mittel b_min: {zeta_impact['weighted_mean_b']:.4f}")
    print(f"Gewichtete Streuung:       {zeta_impact['weighted_std_b']:.4f}")
    print(f"Zufalls-Mittel b_min:      {zeta_impact['random_mean_b']:.4f}")
    print(f"Zufalls-Streuung:          {zeta_impact['random_std_b']:.4f}")
    print(f"Delta Mittelwert:          {zeta_impact['delta_mean_b']:.4f}")
    if zeta_impact["delta_mean_b"] < 0:
        print("Interpretation: Die Zeta-Gewichtung senkt im Mittel die kritische Blockskala.")
    elif zeta_impact["delta_mean_b"] > 0:
        print("Interpretation: Die Zeta-Gewichtung erhoeht im Mittel die kritische Blockskala.")
    else:
        print("Interpretation: Gewichtetes Modell und Zufallsmodell liegen gleichauf.")


if __name__ == "__main__":
    main()