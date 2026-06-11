from decimal import Decimal, getcontext
from pathlib import Path

import numpy as np


getcontext().prec = 50
SCRIPT_VERSION = "Alpha_EABC.py :: Zusatztests-v2 :: 2026-06-07"

M_Pl = Decimal("1.2209e19")
alpha_G_exact = Decimal(1) / Decimal(20)
alpha_codata = Decimal(1) / Decimal("137.035999177")
epsilon_CODATA = Decimal("137.035999177") - Decimal(137)
DEFAULT_CENTER = 138
DEFAULT_WINDOW = 20

MODEL_16 = 16
MODEL_4 = 4
MODEL_8 = 8
MODEL_20 = 20
MODEL_137 = 137
MODEL_320 = MODEL_16 * MODEL_20
APRIORI_CENTER = MODEL_137 - (MODEL_16 + MODEL_4) - MODEL_4
APRIORI_WINDOW = MODEL_4
APRIORI_PRIME_INDEX = 2 * (MODEL_16 - 1)

base_dir = Path(__file__).resolve().parent


def load_norm_gaps(path: Path) -> np.ndarray:
    zeros = np.load(path)
    gaps = np.diff(zeros)
    mid = zeros[:-1]
    return gaps * np.log(mid / (2 * np.pi)) / (2 * np.pi)


def zeta_witness(norm_gaps: np.ndarray, center: int, window: int):
    local = norm_gaps[center - window : center + window]
    median_gap = np.median(local)
    epsilon_zeta = 1 - median_gap
    alpha_inv_zeta = 137 + epsilon_zeta
    alpha_zeta = Decimal(str(1 / alpha_inv_zeta))
    return epsilon_zeta, alpha_zeta, alpha_inv_zeta


def delta_epsilon(epsilon_zeta: float) -> Decimal:
    return abs(Decimal(str(epsilon_zeta)) - epsilon_CODATA)


# Monopolwerte: Dirac- und BPS-Monopolwerte für einen gegebenen alpha-Kandidaten.
def monopole_observables(alpha: Decimal, alpha_inv, M_X: Decimal, n: int = 1):
    """Dirac- und BPS-Monopolwerte für einen gegebenen alpha-Kandidaten.

    Natürliche Einheiten: alpha = e^2/(4*pi), Dirac: eg/(4*pi)=n/2.
    Daraus folgen alpha_m = g^2/(4*pi) = n^2/(4 alpha)
    und g/e = n/(2 alpha).
    """
    alpha_inv_dec = Decimal(str(alpha_inv))
    n_dec = Decimal(n)
    alpha_m = (n_dec * n_dec) / (Decimal(4) * alpha)
    g_over_e = n_dec / (Decimal(2) * alpha)
    bps_mass_backtest = M_X / alpha
    bps_relative_error = abs((bps_mass_backtest - M_M) / M_M)
    return {
        "n": n,
        "alpha_inv": alpha_inv_dec,
        "alpha_m": alpha_m,
        "alpha_m_from_alpha_inv": (n_dec * n_dec * alpha_inv_dec) / Decimal(4),
        "g_over_e": g_over_e,
        "M_X_over_M_M": M_X / M_M,
        "M_M_from_BPS": bps_mass_backtest,
        "BPS_relative_error": bps_relative_error,
    }


def print_monopole_block(label: str, alpha: Decimal, alpha_inv, M_X: Decimal):
    obs = monopole_observables(alpha, alpha_inv, M_X)
    print(f"--- Monopolwerte ({label}) ---")
    print("Dirac-Zahl n =", obs["n"])
    print("alpha_m = 1/(4 alpha) =", obs["alpha_m"])
    print("alpha_m aus alpha_inv/4 =", obs["alpha_m_from_alpha_inv"])
    print("g/e = 1/(2 alpha) =", obs["g_over_e"])
    print("M_X / M_M =", obs["M_X_over_M_M"])
    print("BPS-Rücktest M_X / alpha =", obs["M_M_from_BPS"], "GeV")
    print("BPS-relativer Fehler =", obs["BPS_relative_error"])


def scan_top5(norm_gaps: np.ndarray):
    n = len(norm_gaps)
    results = []
    for center in range(100, 201):
        for window in range(3, 31):
            if center - window < 0 or center + window > n:
                continue
            eps, _, alpha_inv = zeta_witness(norm_gaps, center, window)
            results.append((delta_epsilon(eps), center, window, eps, alpha_inv))
    results.sort(key=lambda r: r[0])
    return results[:5]


# --- Zusatzfunktionen für A-PRIORI-Statistik und Tests ---

def trimmed_mean(values: np.ndarray, trim_fraction: float = 0.1) -> float:
    if len(values) == 0:
        return float("nan")
    sorted_values = np.sort(values)
    k = int(len(sorted_values) * trim_fraction)
    if k == 0 or 2 * k >= len(sorted_values):
        return float(np.mean(sorted_values))
    return float(np.mean(sorted_values[k:-k]))


def epsilon_statistics(norm_gaps: np.ndarray, center: int, window: int):
    local = norm_gaps[center - window : center + window]
    if len(local) == 0:
        raise ValueError("Leeres Fenster in epsilon_statistics")
    stats = {
        "median": float(np.median(local)),
        "mean": float(np.mean(local)),
        "trimmed_mean_10pct": trimmed_mean(local, 0.1),
        "q25": float(np.quantile(local, 0.25)),
        "q75": float(np.quantile(local, 0.75)),
        "std": float(np.std(local)),
        "min": float(np.min(local)),
        "max": float(np.max(local)),
    }
    stats["epsilon_median"] = 1.0 - stats["median"]
    stats["epsilon_mean"] = 1.0 - stats["mean"]
    stats["epsilon_trimmed_mean_10pct"] = 1.0 - stats["trimmed_mean_10pct"]
    stats["iqr"] = stats["q75"] - stats["q25"]
    return stats


def epsilon_eabc_corrected(norm_gaps: np.ndarray, center: int, window: int):
    """EABC-ε-Korrektur über Quartile im lokalen Fenster (nur Modellzahlen)."""
    stats = epsilon_statistics(norm_gaps, center, window)
    epsilon_0 = stats["epsilon_median"]
    q25 = stats["q25"]
    q75 = stats["q75"]
    sigma_loc = stats["std"]

    a_quartile = (q75 - 1.0) - (1.0 - q25)
    quartile_correction = (q25 + q75 - 2.0) / MODEL_320
    epsilon_eabc = epsilon_0 + quartile_correction

    corr_sigma_16_137 = sigma_loc / (MODEL_16 * MODEL_137)
    corr_sigma_4_8_137 = (MODEL_4 * sigma_loc) / (MODEL_8 * MODEL_137)
    epsilon_corr_sigma_16_137 = epsilon_0 - corr_sigma_16_137
    epsilon_corr_sigma_4_8_137 = epsilon_0 - corr_sigma_4_8_137

    alpha_inv_eabc = Decimal(MODEL_137) + Decimal(str(epsilon_eabc))
    alpha_eabc = Decimal(1) / alpha_inv_eabc
    delta_eabc = abs(Decimal(str(epsilon_eabc)) - epsilon_CODATA)
    delta_epsilon_0 = abs(Decimal(str(epsilon_0)) - epsilon_CODATA)

    return {
        "epsilon_0": epsilon_0,
        "q25": q25,
        "q75": q75,
        "sigma_loc": sigma_loc,
        "a_quartile": a_quartile,
        "quartile_correction": quartile_correction,
        "epsilon_eabc": epsilon_eabc,
        "corr_sigma_16_137": corr_sigma_16_137,
        "corr_sigma_4_8_137": corr_sigma_4_8_137,
        "epsilon_corr_sigma_16_137": epsilon_corr_sigma_16_137,
        "epsilon_corr_sigma_4_8_137": epsilon_corr_sigma_4_8_137,
        "alpha_inv_eabc": alpha_inv_eabc,
        "alpha_eabc": alpha_eabc,
        "delta_eabc": delta_eabc,
        "delta_epsilon_0": delta_epsilon_0,
    }


def print_epsilon_eabc_block(norm_gaps: np.ndarray, center: int, window: int):
    corr = epsilon_eabc_corrected(norm_gaps, center, window)
    m_x_eabc = corr["alpha_eabc"] * M_M
    improvement = (
        float(corr["delta_epsilon_0"] / corr["delta_eabc"])
        if corr["delta_eabc"] > 0
        else float("inf")
    )

    print("=== EABC ε-Korrektur (Divisionsalgebren) ===")
    print(
        f"Formel: ε_EABC = ε_0 + (q_25 + q_75 - 2) / ({MODEL_16}·{MODEL_20}) "
        f"= ε_0 + (q_25 + q_75 - 2) / {MODEL_320}"
    )
    print(f"Fenster: center = {center}, window = {window}")
    print("q_25 =", corr["q25"], "q_75 =", corr["q75"])
    print(
        "A = (q_75 - 1) - (1 - q_25) = q_25 + q_75 - 2 =",
        corr["a_quartile"],
    )
    print("σ_loc (lokale Std) =", corr["sigma_loc"])
    print()
    print("--- Korrekturterm ---")
    print(
        f"Quartil-Korrektur = (q_25 + q_75 - 2) / {MODEL_320} =",
        corr["quartile_correction"],
    )
    print()
    print("--- Modellinterne Alternativen (nur Vergleich, kein Fit) ---")
    print(
        f"1) ε_corr = ε_0 - σ_loc / ({MODEL_16}·{MODEL_137}) =",
        corr["epsilon_corr_sigma_16_137"],
        f"(Term = {corr['corr_sigma_16_137']})",
    )
    print(
        f"2) ε_corr = ε_0 - {MODEL_4}·σ_loc / ({MODEL_8}·{MODEL_137}) =",
        corr["epsilon_corr_sigma_4_8_137"],
        f"(Term = {corr['corr_sigma_4_8_137']}, zu stark bei naiver Anwendung)",
    )
    print(
        f"3) ε_EABC (Quartil, primär) = ε_0 + (q_25 + q_75 - 2) / {MODEL_320} =",
        corr["epsilon_eabc"],
    )
    print()
    print("--- ε-Vergleich zu CODATA ---")
    print("ε_CODATA =", epsilon_CODATA)
    print("ε_0 (Median) =", corr["epsilon_0"], "|delta| =", corr["delta_epsilon_0"])
    print("ε_EABC (korrigiert) =", corr["epsilon_eabc"], "|delta| =", corr["delta_eabc"])
    print("Verbesserungsfaktor |delta_0| / |delta_EABC| ≈", improvement)
    print()
    print("--- Physikalische Ableitungen (ε_EABC) ---")
    print("1/alpha_EABC =", corr["alpha_inv_eabc"])
    print("alpha_EABC =", corr["alpha_eabc"])
    print("M_X_EABC =", m_x_eabc, "GeV")
    print_monopole_block("EABC ε-korrigiert", corr["alpha_eabc"], corr["alpha_inv_eabc"], m_x_eabc)


def robustness_by_window(norm_gaps: np.ndarray, center: int, max_window: int = 20):
    rows = []
    for window in range(1, max_window + 1):
        if center - window < 0 or center + window > len(norm_gaps):
            continue
        eps, _, alpha_inv = zeta_witness(norm_gaps, center, window)
        rows.append((delta_epsilon(eps), window, eps, alpha_inv))
    rows.sort(key=lambda row: row[0])
    return rows


def neighborhood_by_center(norm_gaps: np.ndarray, center_min: int = 100, center_max: int = 130, window: int = 4):
    rows = []
    for center in range(center_min, center_max + 1):
        if center - window < 0 or center + window > len(norm_gaps):
            continue
        eps, _, alpha_inv = zeta_witness(norm_gaps, center, window)
        rows.append((delta_epsilon(eps), center, eps, alpha_inv))
    rows.sort(key=lambda row: row[0])
    return rows


def block_recurrence_test(norm_gaps: np.ndarray, base_center: int, window: int, block_size: int = 1000, max_blocks: int = 10):
    rows = []
    n = len(norm_gaps)
    for block_index in range(max_blocks):
        offset = block_index * block_size
        center = offset + base_center
        if center - window < 0 or center + window > n:
            break
        eps, _, alpha_inv = zeta_witness(norm_gaps, center, window)
        rows.append((block_index, center, eps, alpha_inv, delta_epsilon(eps)))
    return rows


def monte_carlo_permutation_test(norm_gaps: np.ndarray, center: int, window: int, trials: int = 20000, seed: int = 137):
    """Zufallstest mit zufällig gezogenen Gap-Fenstern gleicher Länge aus denselben normierten Gaps."""
    rng = np.random.default_rng(seed)
    local_len = 2 * window
    target_delta = float(delta_epsilon(zeta_witness(norm_gaps, center, window)[0]))
    n = len(norm_gaps)
    hits = 0
    best_delta = float("inf")
    best_epsilon = None
    best_start = None
    for _ in range(trials):
        start = int(rng.integers(0, n - local_len + 1))
        sample = norm_gaps[start : start + local_len]
        eps = 1.0 - float(np.median(sample))
        d = abs(eps - float(epsilon_CODATA))
        if d <= target_delta:
            hits += 1
        if d < best_delta:
            best_delta = d
            best_epsilon = eps
            best_start = start
    return {
        "trials": trials,
        "target_delta": target_delta,
        "hits": hits,
        "p_value_estimate": hits / trials,
        "best_delta": best_delta,
        "best_epsilon": best_epsilon,
        "best_start": best_start,
    }


def gue_gap_surrogate_test(norm_gaps: np.ndarray, center: int, window: int, trials: int = 5000, matrix_size: int = 160, seed: int = 113):
    """Kleiner GUE-Surrogattest.

    Erzeugt komplex-hermitesche Zufallsmatrizen, nimmt zentrale Eigenwertabstände,
    normiert sie auf mittleren Gap 1 und prüft denselben Median-Defekt.
    Das ist kein vollständiger RMT-Beweis, aber ein nützlicher Kontrolltest.
    """
    rng = np.random.default_rng(seed)
    target_delta = float(delta_epsilon(zeta_witness(norm_gaps, center, window)[0]))
    local_len = 2 * window
    hits = 0
    best_delta = float("inf")
    best_epsilon = None
    for _ in range(trials):
        real = rng.normal(size=(matrix_size, matrix_size))
        imag = rng.normal(size=(matrix_size, matrix_size))
        hermitian = (real + real.T) / 2.0 + 1j * (imag - imag.T) / 2.0
        eigenvalues = np.linalg.eigvalsh(hermitian)
        gaps = np.diff(eigenvalues)
        center_index = len(gaps) // 2
        bulk = gaps[max(0, center_index - 40) : min(len(gaps), center_index + 40)]
        mean_gap = float(np.mean(bulk))
        if mean_gap <= 0:
            continue
        norm_bulk = bulk / mean_gap
        if len(norm_bulk) < local_len:
            continue
        start = len(norm_bulk) // 2 - window
        sample = norm_bulk[start : start + local_len]
        eps = 1.0 - float(np.median(sample))
        d = abs(eps - float(epsilon_CODATA))
        if d <= target_delta:
            hits += 1
        if d < best_delta:
            best_delta = d
            best_epsilon = eps
    return {
        "trials": trials,
        "target_delta": target_delta,
        "hits": hits,
        "p_value_estimate": hits / trials,
        "best_delta": best_delta,
        "best_epsilon": best_epsilon,
    }


def signature_metrics(norm_gaps: np.ndarray, center: int, window: int):
    """Mehrdimensionale EABC-Signatur für ein Fenster (center, window)."""
    eps, alpha, alpha_inv = zeta_witness(norm_gaps, center, window)
    M_X = alpha * M_M

    obs = monopole_observables(alpha, alpha_inv, M_X)
    obs_cod = monopole_observables(alpha_codata, Decimal(1) / alpha_codata, M_X_codata)

    s1_alpha = float(delta_epsilon(eps) / epsilon_CODATA)
    s2_monopole = float(abs((obs["alpha_m"] - obs_cod["alpha_m"]) / obs_cod["alpha_m"]))
    s3_mx = float(abs((M_X - M_X_codata) / M_X_codata))

    candidate_delta = float(delta_epsilon(eps))
    neighbor_deltas = []
    for dc, dw in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        c2 = center + dc
        w2 = window + dw
        if w2 < 1 or c2 - w2 < 0 or c2 + w2 > len(norm_gaps):
            continue
        eps2, _, _ = zeta_witness(norm_gaps, c2, w2)
        neighbor_deltas.append(float(delta_epsilon(eps2)))

    if neighbor_deltas and candidate_delta > 0:
        s4_local_contrast = float(np.mean(neighbor_deltas) / candidate_delta)
    else:
        s4_local_contrast = 0.0

    eps_windows = []
    for w2 in [window - 1, window, window + 1]:
        if w2 < 1 or center - w2 < 0 or center + w2 > len(norm_gaps):
            continue
        eps2, _, _ = zeta_witness(norm_gaps, center, w2)
        eps_windows.append(float(eps2))

    s5_window_std = float(np.std(eps_windows)) if len(eps_windows) >= 2 else float("inf")

    eps_floor = 1e-30
    quality_Q = (
        1.0 / max(s1_alpha, eps_floor)
        * 1.0 / max(s2_monopole, eps_floor)
        * 1.0 / max(s3_mx, eps_floor)
        * max(s4_local_contrast, eps_floor)
        * 1.0 / max(s5_window_std, eps_floor)
    )

    return {
        "center": center,
        "window": window,
        "epsilon": eps,
        "alpha_inv": alpha_inv,
        "alpha": alpha,
        "M_X": M_X,
        "alpha_m": obs["alpha_m"],
        "g_over_e": obs["g_over_e"],
        "s1_alpha": s1_alpha,
        "s2_monopole": s2_monopole,
        "s3_mx": s3_mx,
        "s4_local_contrast": s4_local_contrast,
        "s5_window_std": s5_window_std,
        "Q": quality_Q,
    }


def print_signature_block(title: str, metrics: dict):
    print(title)
    print("center =", metrics["center"], "window =", metrics["window"])
    print("epsilon =", metrics["epsilon"])
    print("1/alpha =", metrics["alpha_inv"])
    print("M_X =", metrics["M_X"], "GeV")
    print("alpha_m =", metrics["alpha_m"])
    print("g/e =", metrics["g_over_e"])
    print("s1_alpha_rel =", metrics["s1_alpha"])
    print("s2_monopole_rel =", metrics["s2_monopole"])
    print("s3_MX_rel =", metrics["s3_mx"])
    print("s4_local_contrast =", metrics["s4_local_contrast"])
    print("s5_window_std =", metrics["s5_window_std"])
    print("Q =", metrics["Q"])


def signature_monte_carlo_test(
    norm_gaps: np.ndarray,
    center: int,
    window: int,
    trials: int = 50000,
    seed: int = 137113,
    center_min: int = 100,
    center_max: int = 5000,
    window_min: int = 2,
    window_max: int = 30,
):
    """Monte-Carlo-Test der gesamten Signatur Q.

    Anders als der einfache epsilon-Test fordert dieser Test gleichzeitig:
    Alpha-Nähe, Monopolkopplung, M_X-Nähe, lokalen Kontrast und
    Fensterrobustheit.
    """
    rng = np.random.default_rng(seed)
    q_target = signature_metrics(norm_gaps, center, window)["Q"]

    n = len(norm_gaps)
    center_max_eff = min(center_max, n - window_max - 2)
    hits = 0
    best_q = -float("inf")
    best = None

    for _ in range(trials):
        c = int(rng.integers(center_min, center_max_eff + 1))
        w = int(rng.integers(window_min, window_max + 1))
        if c - w < 0 or c + w > n:
            continue
        m = signature_metrics(norm_gaps, c, w)
        q = m["Q"]
        if q >= q_target:
            hits += 1
        if q > best_q:
            best_q = q
            best = m

    return {
        "trials": trials,
        "Q_target": q_target,
        "hits": hits,
        "p_estimate": hits / trials,
        "best": best,
    }


def _delta_eabc(norm_gaps: np.ndarray, center: int, window: int) -> float:
    corr = epsilon_eabc_corrected(norm_gaps, center, window)
    return float(corr["delta_eabc"])


def hard_signature_pass(norm_gaps, center, window, target_delta=None):
    if target_delta is None:
        target_delta = _delta_eabc(norm_gaps, APRIORI_CENTER, APRIORI_WINDOW)
    if window not in (3, 4):
        return False
    d0 = _delta_eabc(norm_gaps, center, window)
    if d0 > target_delta:
        return False
    # lokaler Rang im Bereich center +/- 15
    local = []
    for c2 in range(max(window, center - 15), min(len(norm_gaps) - window, center + 15) + 1):
        local.append((_delta_eabc(norm_gaps, c2, window), c2))
    local.sort()
    if local[0][1] != center:
        return False
    # lokaler Kontrast gegen direkte Nachbarn
    neighbor_deltas = []
    for dc, dw in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        c2 = center + dc
        w2 = window + dw
        if w2 >= 1 and c2 - w2 >= 0 and c2 + w2 <= len(norm_gaps):
            neighbor_deltas.append(_delta_eabc(norm_gaps, c2, w2))
    if not neighbor_deltas:
        return False
    contrast = float(np.mean(neighbor_deltas)) / d0 if d0 > 0 else float("inf")
    return contrast >= 100


def hard_signature_scan(norm_gaps, center_min=100, center_max=5000, window_min=2, window_max=30):
    rows = []
    total = 0
    for c in range(center_min, min(center_max, len(norm_gaps) - window_max - 1) + 1):
        for w in range(window_min, window_max + 1):
            if c - w < 0 or c + w > len(norm_gaps):
                continue
            total += 1
            if hard_signature_pass(norm_gaps, c, w):
                corr = epsilon_eabc_corrected(norm_gaps, c, w)
                rows.append(
                    (
                        c,
                        w,
                        corr["epsilon_eabc"],
                        corr["alpha_inv_eabc"],
                        float(corr["delta_eabc"]),
                    )
                )
    return total, rows


def factor_int(n: int):
    factors = []
    d = 2
    while d * d <= n:
        count = 0
        while n % d == 0:
            n //= d
            count += 1
        if count:
            factors.append((d, count))
        d += 1 if d == 2 else 2
    if n > 1:
        factors.append((n, 1))
    return factors


def is_prime_int(n: int) -> bool:
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def digit_sum(n: int) -> int:
    return sum(int(ch) for ch in str(abs(n)))


def prime_pair_structure(n: int):
    return {
        "n": n,
        "n+1": n + 1,
        "2n+1": 2 * n + 1,
        "n_prime": is_prime_int(n),
        "2n+1_prime": is_prime_int(2 * n + 1),
        "n_mod16": n % 16,
        "n_mod4": n % 4,
        "pair_sum": n + (n + 1),
        "n_times_n_minus_1": n * (n - 1),
        "n_times_n_plus_1": n * (n + 1),
        "factors_n": factor_int(n),
        "factors_2n_plus_1": factor_int(2 * n + 1),
    }


def eabc_center_diagnostics(norm_gaps: np.ndarray, centers=(113, 2494, 4125), windows=(3, 4)):
    rows = []
    for c in centers:
        for w in windows:
            if c - w < 0 or c + w > len(norm_gaps):
                continue
            eps, alpha, alpha_inv = zeta_witness(norm_gaps, c, w)
            d = float(delta_epsilon(eps))
            metrics = signature_metrics(norm_gaps, c, w)
            rows.append({
                "center": c,
                "window": w,
                "epsilon": eps,
                "alpha_inv": alpha_inv,
                "delta": d,
                "is_prime": is_prime_int(c),
                "factors": factor_int(c),
                "digit_sum": digit_sum(c),
                "mod4": c % 4,
                "mod16": c % 16,
                "mod20": c % 20,
                "mod137": c % 137,
                "distance_to_137k": min(c % 137, 137 - (c % 137)),
                "s4_local_contrast": metrics["s4_local_contrast"],
                "s5_window_std": metrics["s5_window_std"],
                "M_X": metrics["M_X"],
                "alpha_m": metrics["alpha_m"],
                "g_over_e": metrics["g_over_e"],
            })
    return rows


def blindtest_other_zero_files(zero_files, primary_label: str, center: int, window: int):
    """Blindtest: dieselbe eingefrorene Regel (center, window) auf weiteren Dateien."""
    rows = []
    for label, path in zero_files:
        if label == primary_label:
            continue
        norm_gaps = load_norm_gaps(path)
        rows.append((label, signature_metrics(norm_gaps, center, window)))
    return rows


M_M = alpha_G_exact.sqrt() * M_Pl
M_X_codata = alpha_codata * M_M

print("=== SCRIPT-CHECK ===")
print("Version =", SCRIPT_VERSION)
print("Datei =", Path(__file__).resolve())
print()

print("=== CODATA ===")
print("alpha =", alpha_codata)
print("1/alpha =", Decimal(1) / alpha_codata)
print("epsilon_CODATA =", epsilon_CODATA)
print("M_X =", M_X_codata, "GeV")
print_monopole_block("CODATA", alpha_codata, Decimal(1) / alpha_codata, M_X_codata)

zero_files = [("zeros6.npy", base_dir / "zeros6.npy")]
optional_path = base_dir / "zeros6(9).npy"
if optional_path.exists():
    zero_files.append(("zeros6(9).npy", optional_path))

primary_norm_gaps = load_norm_gaps(zero_files[0][1])
corr_ap = epsilon_eabc_corrected(primary_norm_gaps, APRIORI_CENTER, APRIORI_WINDOW)
eps_ap_0 = corr_ap["epsilon_0"]
eps_ap = corr_ap["epsilon_eabc"]
alpha_ap = corr_ap["alpha_eabc"]
alpha_inv_ap = corr_ap["alpha_inv_eabc"]
d_eps_ap = corr_ap["delta_eabc"]
d_eps_ap_0 = corr_ap["delta_epsilon_0"]
M_X_ap = alpha_ap * M_M

print()
print("=== A-PRIORI ===")
print(
    f"Regel: c = {MODEL_137} - ({MODEL_16}+{MODEL_4}) - {MODEL_4} = {APRIORI_CENTER}, "
    f"w = {APRIORI_WINDOW}"
)
print("Modellzahlen:", {MODEL_16, MODEL_4, MODEL_20, MODEL_137})
print(
    f"Primstützung: {APRIORI_CENTER} = p_{APRIORI_PRIME_INDEX}, "
    f"{APRIORI_PRIME_INDEX} = 2({MODEL_16}-1)"
)
print("center =", APRIORI_CENTER, "window =", APRIORI_WINDOW)
print("ε_0 (Median, Referenz) =", eps_ap_0, "|delta| =", d_eps_ap_0)
print("ε_EABC (primär) =", eps_ap)
print("1/alpha =", alpha_inv_ap)
print("|epsilon_EABC - epsilon_CODATA| =", d_eps_ap)
print(
    f"|delta| / epsilon_CODATA ≈ {float(d_eps_ap / epsilon_CODATA * 100):.4f}% vs CODATA"
)
print("M_X =", M_X_ap, "GeV")
print_monopole_block("A-PRIORI (ε_EABC)", alpha_ap, alpha_inv_ap, M_X_ap)

print()
print_epsilon_eabc_block(primary_norm_gaps, APRIORI_CENTER, APRIORI_WINDOW)

print()
print("=== n(n±1)-Kern-Rand-Kopplung ===")
n_core = APRIORI_CENTER
c_core, w_core = n_core, 4
c_edge, w_edge = n_core + 1, 3

eps_core, alpha_core, alpha_inv_core = zeta_witness(primary_norm_gaps, c_core, w_core)
eps_edge, alpha_edge, alpha_inv_edge = zeta_witness(primary_norm_gaps, c_edge, w_edge)

weight_core = Decimal(n_core) / Decimal(2 * n_core + 1)
weight_edge = Decimal(n_core + 1) / Decimal(2 * n_core + 1)

eps_pair = weight_core * Decimal(str(eps_core)) + weight_edge * Decimal(str(eps_edge))
alpha_inv_pair = Decimal(137) + eps_pair
alpha_pair = Decimal(1) / alpha_inv_pair
M_X_pair = alpha_pair * M_M

print("n =", n_core)
print("n(n-1) =", n_core * (n_core - 1))
print("n(n+1) =", n_core * (n_core + 1))
print("2n+1 =", 2 * n_core + 1)
print("core: center =", c_core, "window =", w_core, "epsilon =", eps_core)
print("edge: center =", c_edge, "window =", w_edge, "epsilon =", eps_edge)
print("weight_core =", weight_core)
print("weight_edge =", weight_edge)
print("epsilon_pair =", eps_pair)
print("epsilon_CODATA =", epsilon_CODATA)
print("|epsilon_pair - epsilon_CODATA| =", abs(eps_pair - epsilon_CODATA))
print("1/alpha_pair =", alpha_inv_pair)
print("alpha_pair =", alpha_pair)
print("M_X_pair =", M_X_pair, "GeV")
print_monopole_block("n(n+1)-Kopplung", alpha_pair, alpha_inv_pair, M_X_pair)

corr_core_eabc = epsilon_eabc_corrected(primary_norm_gaps, c_core, w_core)
corr_edge_eabc = epsilon_eabc_corrected(primary_norm_gaps, c_edge, w_edge)
eps_core_eabc = corr_core_eabc["epsilon_eabc"]
eps_edge_eabc = corr_edge_eabc["epsilon_eabc"]

eps_pair_eabc = (
    weight_core * Decimal(str(eps_core_eabc))
    + weight_edge * Decimal(str(eps_edge_eabc))
)
alpha_inv_pair_eabc = Decimal(137) + eps_pair_eabc
alpha_pair_eabc = Decimal(1) / alpha_inv_pair_eabc
M_X_pair_eabc = alpha_pair_eabc * M_M
delta_pair = abs(eps_pair - epsilon_CODATA)
delta_pair_eabc = abs(eps_pair_eabc - epsilon_CODATA)
delta_core_eabc = corr_core_eabc["delta_eabc"]

print()
print("--- ε_EABC-Kopplung (primär) ---")
print("core: center =", c_core, "window =", w_core, "ε_EABC =", eps_core_eabc)
print("edge: center =", c_edge, "window =", w_edge, "ε_EABC =", eps_edge_eabc)
print("ε_EABC pair (primär) =", eps_pair_eabc)
print("|ε_EABC pair - ε_CODATA| =", delta_pair_eabc)
print("1/alpha_pair_EABC =", alpha_inv_pair_eabc)
print("alpha_pair_EABC =", alpha_pair_eabc)
print("M_X_pair_EABC =", M_X_pair_eabc, "GeV")
print_monopole_block("n(n+1)-Kopplung ε_EABC", alpha_pair_eabc, alpha_inv_pair_eabc, M_X_pair_eabc)

print()
print("--- Vergleich Median-Paar vs ε_EABC-Paar vs einzelner Kern ε_EABC ---")
print("ε_CODATA =", epsilon_CODATA)
print("|delta| Median-Paar =", delta_pair)
print("|delta| ε_EABC-Paar (primär) =", delta_pair_eabc)
print("|delta| einzelner Kern ε_EABC (c=113, w=4) =", delta_core_eabc)
print("Median-Paar vs Kern ε_EABC: Faktor =", float(delta_pair / delta_core_eabc) if delta_core_eabc > 0 else float("inf"))
print(
    "ε_EABC-Paar vs Kern ε_EABC: Faktor =",
    float(delta_pair_eabc / delta_core_eabc) if delta_core_eabc > 0 else float("inf"),
)
print(
    "ε_EABC-Paar vs Median-Paar: Faktor =",
    float(delta_pair_eabc / delta_pair) if delta_pair > 0 else float("inf"),
)

print()
print("=== Zweite EABC-Randkorrektur ===")
sigma_core = epsilon_statistics(primary_norm_gaps, 113, 4)["std"]
sigma_edge = epsilon_statistics(primary_norm_gaps, 114, 3)["std"]
sigma_bar = (Decimal(str(sigma_core)) + Decimal(str(sigma_edge))) / Decimal(2)
second_term = sigma_bar / Decimal(16 * 20 * 137)
eps_final = eps_pair_eabc - second_term
alpha_inv_final = Decimal(137) + eps_final
alpha_final = Decimal(1) / alpha_inv_final
M_X_final = alpha_final * M_M

print("sigma_core =", sigma_core)
print("sigma_edge =", sigma_edge)
print("sigma_bar =", sigma_bar)
print("second_term = sigma_bar / (16*20*137) =", second_term)
print("epsilon_pair_EABC =", eps_pair_eabc)
print("epsilon_final =", eps_final)
print("epsilon_CODATA =", epsilon_CODATA)
print("|epsilon_final - epsilon_CODATA| =", abs(eps_final - epsilon_CODATA))
print("1/alpha_final =", alpha_inv_final)
print("alpha_final =", alpha_final)
print("M_X_final =", M_X_final, "GeV")
print_monopole_block("zweite EABC-Randkorrektur", alpha_final, alpha_inv_final, M_X_final)

print()
print("=== n(n±1)-Strukturtest der harten Treffer ===")
hard_centers = [113, 114, 3876]
for c in hard_centers:
    n = c if c % 2 == 1 else c - 1
    s = prime_pair_structure(n)
    print()
    print("Ausgangspunkt center =", c)
    print("gewähltes n =", s["n"])
    print("n+1 =", s["n+1"])
    print("2n+1 =", s["2n+1"])
    print("n prim =", s["n_prime"])
    print("2n+1 prim =", s["2n+1_prime"])
    print("n mod 4 =", s["n_mod4"])
    print("n mod 16 =", s["n_mod16"])
    print("n(n-1) =", s["n_times_n_minus_1"])
    print("n(n+1) =", s["n_times_n_plus_1"])
    print("Faktoren n =", s["factors_n"])
    print("Faktoren 2n+1 =", s["factors_2n_plus_1"])

print()
print("=== Zusatztests A-PRIORI-Regel ===")

stats_ap = epsilon_statistics(primary_norm_gaps, APRIORI_CENTER, APRIORI_WINDOW)
print("Statistiken im A-priori-Fenster:")
print("median =", stats_ap["median"], "epsilon_median =", stats_ap["epsilon_median"])
print("mean =", stats_ap["mean"], "epsilon_mean =", stats_ap["epsilon_mean"])
print("trimmed_mean_10pct =", stats_ap["trimmed_mean_10pct"], "epsilon_trimmed =", stats_ap["epsilon_trimmed_mean_10pct"])
print("q25 =", stats_ap["q25"], "q75 =", stats_ap["q75"], "iqr =", stats_ap["iqr"])
print("std =", stats_ap["std"], "min =", stats_ap["min"], "max =", stats_ap["max"])

print()
print("Fensterrobustheit bei center=113 — Top 10 Fensterbreiten:")
for rank, (d, window, eps_w, alpha_inv_w) in enumerate(
    robustness_by_window(primary_norm_gaps, APRIORI_CENTER, 20)[:10], 1
):
    print(
        f"{rank}. window={window:2d} epsilon={eps_w:.12f} "
        f"1/alpha={alpha_inv_w:.12f} |delta|={d:.12e}"
    )

print()
print("Nachbarschaftstest bei window=4 — Top 10 Zentren:")
for rank, (d, center, eps_c, alpha_inv_c) in enumerate(
    neighborhood_by_center(primary_norm_gaps, 100, 130, APRIORI_WINDOW)[:10], 1
):
    print(
        f"{rank}. center={center:3d} epsilon={eps_c:.12f} "
        f"1/alpha={alpha_inv_c:.12f} |delta|={d:.12e}"
    )

print()
print("Block-Wiederkehrtest: center = block_offset + 113, window=4")
for block_index, center, eps_b, alpha_inv_b, d_b in block_recurrence_test(
    primary_norm_gaps, APRIORI_CENTER, APRIORI_WINDOW, block_size=1000, max_blocks=10
):
    print(
        f"block={block_index:2d} center={center:6d} epsilon={eps_b:.12f} "
        f"1/alpha={alpha_inv_b:.12f} |delta|={d_b:.12e}"
    )

print()
print("Monte-Carlo-Fenstertest aus denselben Riemann-Gaps:")
mc = monte_carlo_permutation_test(primary_norm_gaps, APRIORI_CENTER, APRIORI_WINDOW)
print("trials =", mc["trials"])
print("target_delta =", mc["target_delta"])
print("hits =", mc["hits"])
print("p_estimate =", mc["p_value_estimate"])
print("best_delta =", mc["best_delta"], "best_epsilon =", mc["best_epsilon"], "best_start =", mc["best_start"])

print()
print("GUE-Surrogattest:")
gue = gue_gap_surrogate_test(primary_norm_gaps, APRIORI_CENTER, APRIORI_WINDOW)
print("trials =", gue["trials"])
print("target_delta =", gue["target_delta"])
print("hits =", gue["hits"])
print("p_estimate =", gue["p_value_estimate"])
print("best_delta =", gue["best_delta"], "best_epsilon =", gue["best_epsilon"])

print()
print("=== Mehrdimensionale EABC-Signatur ===")
sig_ap = signature_metrics(primary_norm_gaps, APRIORI_CENTER, APRIORI_WINDOW)
print_signature_block("A-priori-Signatur:", sig_ap)

print()
print("Signatur-Monte-Carlo-Test:")
sig_mc = signature_monte_carlo_test(primary_norm_gaps, APRIORI_CENTER, APRIORI_WINDOW)
print("trials =", sig_mc["trials"])
print("Q_target =", sig_mc["Q_target"])
print("hits =", sig_mc["hits"])
print("p_estimate =", sig_mc["p_estimate"])
if sig_mc["best"] is not None:
    print_signature_block("Beste zufällige Signatur:", sig_mc["best"])

print()
print("=== Harte EABC-Signatur ===")
total_hard, hard_rows = hard_signature_scan(primary_norm_gaps)
print("Kandidaten gesamt =", total_hard)
print("harte Treffer =", len(hard_rows))
print("p_hard ≈", len(hard_rows) / total_hard if total_hard else None)
for i, (c, w, eps_h, alpha_inv_h, d_h) in enumerate(hard_rows[:20], 1):
    print(
        f"{i}. center={c} window={w} "
        f"epsilon={eps_h:.12f} 1/alpha={alpha_inv_h:.12f} |delta|={d_h:.12e}"
    )

print()
print("=== Diagnose der harten Trefferzentren ===")
diag_rows = eabc_center_diagnostics(primary_norm_gaps)
for row in diag_rows:
    print()
    print(f"center={row['center']} window={row['window']}")
    print("epsilon =", row["epsilon"])
    print("1/alpha =", row["alpha_inv"])
    print("|delta epsilon| =", row["delta"])
    print("is_prime =", row["is_prime"])
    print("factors =", row["factors"])
    print("digit_sum =", row["digit_sum"])
    print("mod4 =", row["mod4"], "mod16 =", row["mod16"], "mod20 =", row["mod20"], "mod137 =", row["mod137"])
    print("distance_to_137k =", row["distance_to_137k"])
    print("local_contrast =", row["s4_local_contrast"])
    print("window_std =", row["s5_window_std"])
    print("M_X =", row["M_X"], "GeV")
    print("alpha_m =", row["alpha_m"])
    print("g/e =", row["g_over_e"])

print()
print("Blindtest auf weiteren Nullstellen-Dateien:")
blind_rows = blindtest_other_zero_files(zero_files, "zeros6.npy", APRIORI_CENTER, APRIORI_WINDOW)
if not blind_rows:
    print("Keine zweite Nullstellen-Datei gefunden. Erwartet optional: zeros6(9).npy")
else:
    for label, metrics in blind_rows:
        print_signature_block(f"Blindtest {label}:", metrics)

for label, path in zero_files:
    norm_gaps = load_norm_gaps(path)
    eps, alpha_zeta, _ = zeta_witness(norm_gaps, DEFAULT_CENTER, DEFAULT_WINDOW)
    d_eps = delta_epsilon(eps)
    M_X_zeta = alpha_zeta * M_M

    print()
    print(f"=== Zeta ({label}) ===")
    print("center =", DEFAULT_CENTER, "window =", DEFAULT_WINDOW)
    print("epsilon_zeta =", eps)
    print("|epsilon - epsilon_CODATA| =", d_eps)
    print("alpha =", alpha_zeta)
    print("1/alpha =", Decimal(1) / alpha_zeta)
    print("M_X =", M_X_zeta, "GeV")
    print_monopole_block(f"Zeta {label}", alpha_zeta, Decimal(1) / alpha_zeta, M_X_zeta)

    print()
    print(f"=== Kontrollscan 137-Zone ({label}) — Top 5 ===")
    for rank, (d, center, window, eps_s, alpha_inv_s) in enumerate(scan_top5(norm_gaps), 1):
        alpha_s = Decimal(str(1 / alpha_inv_s))
        M_X_s = alpha_s * M_M
        obs_s = monopole_observables(alpha_s, alpha_inv_s, M_X_s)
        print(
            f"{rank}. center={center} window={window} "
            f"epsilon={eps_s:.12f} 1/alpha={alpha_inv_s:.12f} |delta|={d:.12e} "
            f"M_X={M_X_s:.6e} GeV alpha_m={obs_s['alpha_m']:.12f} g/e={obs_s['g_over_e']:.6f}"
        )

print()
print("=== Gemeinsam ===")
print("alpha_G =", alpha_G_exact)
print("M_M / M_Pl =", M_M / M_Pl)
print("M_M =", M_M, "GeV")
print("1/alpha_G =", Decimal(1) / alpha_G_exact)


def modulo_profile_deep():
    numbers = [113, 114, 227, 228]
    moduli = [4, 8, 16, 20, 24, 32, 40, 64, 80, 96, 112, 113, 114, 137, 160, 227, 228, 240, 320]
    print()
    print("=== Tiefe Modulo-Struktur: 113, 114, 227, 228 ===")
    for x in numbers:
        print()
        print(f"x = {x}")
        print("prime =", is_prime_int(x), "factors =", factor_int(x), "digit_sum =", digit_sum(x))
        for m in moduli:
            r = x % m
            flags = []
            if r == 0:
                flags.append("0")
            if r == 1:
                flags.append("+1")
            if r == m - 1:
                flags.append("-1")
            if m % 2 == 0 and r == m // 2:
                flags.append("half")
            if m % 4 == 0 and r in (m // 4, 3 * m // 4):
                flags.append("quarter")
            print(f"mod {m:3d}: r={r:3d}, flags={','.join(flags) if flags else '-'}")
    print()
    print("=== EABC-Sequenztest ===")
    n = 113
    seq = [n, n + 1, 2 * n + 1, 2 * n + 2]
    print("n =", n)
    print("[n, n+1, 2n+1, 2n+2] =", seq)
    print("mod 16 =", [x % 16 for x in seq])
    print("mod 20 =", [x % 20 for x in seq])
    print("mod 137 =", [x % 137 for x in seq])
    print("prime flags =", [is_prime_int(x) for x in seq])
    print()
    print("Kernrelationen:")
    print("113 + 114 =", 113 + 114)
    print("227 + 1 =", 227 + 1)
    print("2 * 114 =", 2 * 114)
    print("227 = 2*113 + 1 ->", 2 * 113 + 1)
    print("228 = 2*(113+1) ->", 2 * (113 + 1))


modulo_profile_deep()
