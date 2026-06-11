import math
from fractions import Fraction
import numpy as np
from pathlib import Path


# ============================================================
# EABC-RaQM Toy Model
# Diskreter Hilbertraum, rationalisiertes Hadamard,
# Bell/CHSH-Test, Kohärenz-Skalierung
# ============================================================


# ------------------------------------------------------------
# 1. EABC-Struktur
# ------------------------------------------------------------

EABC = {
    "E": 1,
    "A": 5,
    "B": 7,
    "C": 11,
}

# Klein-Vierergruppe modulo 12
def eabc_mul(x, y):
    return (EABC[x] * EABC[y]) % 12


def eabc_name(value):
    for k, v in EABC.items():
        if v == value:
            return k
    return None


def eabc_product(x, y):
    return eabc_name(eabc_mul(x, y))


print("EABC-Multiplikation:")
for x in EABC:
    for y in EABC:
        print(f"{x} * {y} = {eabc_product(x, y)}")
    print()


# ------------------------------------------------------------
# 2. Rationalisierung von 1/sqrt(2)
# ------------------------------------------------------------

def best_rational_approx(x, max_denominator):
    return Fraction(x).limit_denominator(max_denominator)


target = 1 / math.sqrt(2)

print("\nRationalisierte Hadamard-Faktoren:")
for qmax in [10, 100, 1000, 10_000, 100_000]:
    f = best_rational_approx(target, qmax)
    err = abs(float(f) - target)
    print(f"max_den={qmax:>6}  approx={f}  wert={float(f):.12f}  fehler={err:.3e}")


# ------------------------------------------------------------
# 3. Standard-Hadamard und EABC-rationalisiertes Hadamard
# ------------------------------------------------------------

def standard_hadamard():
    return (1 / math.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=float)


def rational_hadamard(max_denominator=1000):
    f = best_rational_approx(1 / math.sqrt(2), max_denominator)
    return float(f) * np.array([[1, 1], [1, -1]], dtype=float), f


H_std = standard_hadamard()
H_rat, frac = rational_hadamard(1000)

print("\nStandard-Hadamard:")
print(H_std)

print("\nRationalisiertes EABC-Hadamard:")
print(frac)
print(H_rat)

print("\nUnitaritätsfehler H_rat^T H_rat - I:")
unit_error = H_rat.T @ H_rat - np.eye(2)
print(unit_error)
print("Frobenius-Fehler:", np.linalg.norm(unit_error))


# ------------------------------------------------------------
# 4. Fehlerakkumulation unter wiederholten Hadamard-Anwendungen
# ------------------------------------------------------------

def hadamard_error_growth(max_steps=50, max_denominator=1000):
    Hs = standard_hadamard()
    Hr, f = rational_hadamard(max_denominator)

    psi_std = np.array([1.0, 0.0])
    psi_rat = np.array([1.0, 0.0])

    errors = []

    for step in range(1, max_steps + 1):
        psi_std = Hs @ psi_std
        psi_rat = Hr @ psi_rat

        # Renormierung optional:
        # psi_rat = psi_rat / np.linalg.norm(psi_rat)

        err = np.linalg.norm(psi_std - psi_rat)
        norm_rat = np.linalg.norm(psi_rat)
        errors.append((step, err, norm_rat))

    return errors


print("\nFehlerwachstum nach Hadamard-Kaskaden:")
errors = hadamard_error_growth(max_steps=20, max_denominator=1000)
for step, err, norm_rat in errors:
    print(f"step={step:>3}  zustandsfehler={err:.3e}  norm_rat={norm_rat:.12f}")


# ------------------------------------------------------------
# 5. EABC-diskretisierte Phasen
# ------------------------------------------------------------

def discrete_phase(k, N):
    return 2 * math.pi * k / N


def closest_discrete_phase(theta, N):
    k = round(theta * N / (2 * math.pi)) % N
    theta_q = discrete_phase(k, N)
    return k, theta_q, abs(theta - theta_q)


print("\nDiskrete Phasenapproximation:")
for N in [12, 24, 60, 120, 360, 1000]:
    theta = math.pi / 7
    k, theta_q, err = closest_discrete_phase(theta, N)
    print(f"N={N:>4}  k={k:>4}  theta_q={theta_q:.8f}  fehler={err:.3e}")


# ------------------------------------------------------------
# 6. Bell/CHSH-Korrelation
# Standard-QM: E(a,b) = -cos(a-b)
# EABC: Winkel werden auf diskretes Phasengitter projiziert
# ------------------------------------------------------------

def corr_standard(a, b):
    return -math.cos(a - b)


def corr_eabc(a, b, N):
    _, aq, _ = closest_discrete_phase(a, N)
    _, bq, _ = closest_discrete_phase(b, N)
    return -math.cos(aq - bq)


def chsh_value(N=None):
    # Standard-Wahl für maximale Quantenverletzung
    a = 0
    ap = math.pi / 2
    b = math.pi / 4
    bp = -math.pi / 4

    if N is None:
        Eab = corr_standard(a, b)
        Eabp = corr_standard(a, bp)
        Eapb = corr_standard(ap, b)
        Eapbp = corr_standard(ap, bp)
    else:
        Eab = corr_eabc(a, b, N)
        Eabp = corr_eabc(a, bp, N)
        Eapb = corr_eabc(ap, b, N)
        Eapbp = corr_eabc(ap, bp, N)

    # CHSH-Konvention
    S = abs(Eab + Eabp + Eapb - Eapbp)
    return S


print("\nCHSH-Test:")
print("Standard-QM S =", chsh_value())

for N in [8, 12, 16, 24, 32, 64, 128, 256, 512, 1024]:
    print(f"EABC N={N:>4}  S={chsh_value(N):.8f}")


# ------------------------------------------------------------
# 7. Kohärenzmodell gegen Qubit-Zahl
# Palmer-artige Skalierung
# ------------------------------------------------------------

def coherence_exponential(n, epsilon):
    """
    Toy-Annahme:
    globale Fehlerzahl wächst ~ 2^n.
    Kohärenz K = exp(-epsilon * 2^n)
    """
    return math.exp(-epsilon * (2 ** n))


def find_nmax(epsilon, threshold=0.5, n_limit=2000):
    for n in range(n_limit + 1):
        if coherence_exponential(n, epsilon) < threshold:
            return n
    return None


print("\nKohärenzgrenze K(n)=exp(-epsilon*2^n):")

for exponent in [20, 40, 80, 120, 160]:
    epsilon = 10 ** (-exponent)
    nmax = find_nmax(epsilon, threshold=0.5, n_limit=1000)
    print(f"epsilon=1e-{exponent:<3}  n_max≈{nmax}")


# ------------------------------------------------------------
# 8. Ziel: epsilon so wählen, dass n_max ~ 400
# ------------------------------------------------------------

target_n = 400
threshold = 0.5

# exp(-epsilon*2^n)=threshold
# epsilon = -ln(threshold)/2^n

epsilon_for_400 = -math.log(threshold) / (2 ** target_n)

print("\nEpsilon für n_max≈400:")
print("epsilon ≈", epsilon_for_400)
print("log10(epsilon) ≈", math.log10(epsilon_for_400))


# ------------------------------------------------------------
# 9. EABC-Holonomie-Spielzeugmodell
# ------------------------------------------------------------

def eabc_phase_factor(channel, k, N):
    """
    Kanal + diskrete Phase.
    Für das Toy-Modell nehmen wir komplexe Phasen.
    Quaternionische Nichtkommutativität wird hier noch nicht voll eingebaut.
    """
    theta = discrete_phase(k, N)

    channel_phase = {
        "E": 0,
        "A": math.pi / 2,
        "B": math.pi,
        "C": 3 * math.pi / 2,
    }[channel]

    return np.exp(1j * (theta + channel_phase))


def channel_quarter_turns(channel):
    """
    Phase in Viertelumdrehungen:
    E=0, A=1, B=2, C=3  (entspricht 0, pi/2, pi, 3pi/2)
    """
    return {"E": 0, "A": 1, "B": 2, "C": 3}[channel]


def loop_group_product(loop):
    p = "E"
    for channel, _ in loop:
        p = eabc_product(p, channel)
    return p


def loop_phase_units(loop, N):
    """
    Gesamtphase in Einheiten von 2*pi/N.
    Exakte Schließung: loop_phase_units(loop, N) % N == 0.
    Voraussetzung: N durch 4 teilbar.
    """
    if N % 4 != 0:
        raise ValueError("N muss durch 4 teilbar sein für exakte pi/2-Phasen.")
    quarter_unit = N // 4
    return sum(k for _, k in loop) + quarter_unit * sum(channel_quarter_turns(c) for c, _ in loop)


def is_compatible_loop(loop, N):
    if N % 4 != 0:
        return False
    group_ok = loop_group_product(loop) == "E"
    phase_ok = (loop_phase_units(loop, N) % N) == 0
    return group_ok and phase_ok


def loop_holonomy(loop, N):
    """
    loop = Liste von Tupeln: (channel, k)
    """
    z = 1 + 0j
    for channel, k in loop:
        z *= eabc_phase_factor(channel, k, N)
    return z


def random_loop(length=6, N=120):
    return random_loop_with_rng(np.random.default_rng(), length=length, N=N)


def random_loop_with_rng(rng, length=6, N=120):
    channels = list(EABC.keys())
    return [(rng.choice(channels), int(rng.integers(0, N))) for _ in range(length)]


def generate_compatible_loop(rng, length=6, N=120):
    """
    Erzeugt einen Loop mit
    - Gruppenprodukt = E
    - Phasensumme = 0 mod 2*pi
    """
    if length < 2:
        raise ValueError("Loop-Länge muss >= 2 sein.")
    if N % 4 != 0:
        raise ValueError("Für kompatible Loops muss N durch 4 teilbar sein.")

    channels = list(EABC.keys())
    prefix_channels = [rng.choice(channels) for _ in range(length - 1)]
    prefix_k = [int(rng.integers(0, N)) for _ in range(length - 1)]

    # Schließe Gruppenprodukt auf E:
    # (prod_prefix) * last = E  ->  last = prod_prefix (im Klein-Vierer).
    prod_prefix = "E"
    for c in prefix_channels:
        prod_prefix = eabc_product(prod_prefix, c)
    last_channel = prod_prefix

    quarter_unit = N // 4
    prefix_phase_units = sum(prefix_k) + quarter_unit * sum(channel_quarter_turns(c) for c in prefix_channels)
    last_phase_units = quarter_unit * channel_quarter_turns(last_channel)

    last_k = int((-prefix_phase_units - last_phase_units) % N)
    loop = list(zip(prefix_channels, prefix_k)) + [(last_channel, last_k)]
    return loop


def generate_frustrated_loop(rng, length=6, N=120):
    """
    Startet von kompatiblem Loop und stört gezielt die Phasenschließung.
    Gruppenprodukt bleibt E, Holonomie wird i.A. != 1.
    """
    loop = generate_compatible_loop(rng, length=length, N=N)
    delta_k = int(rng.integers(1, N))  # garantiert != 0
    channel, k = loop[-1]
    loop[-1] = (channel, (k + delta_k) % N)
    return loop


def holonomy_coherence(num_loops=10_000, loop_length=6, N=120, rng=None):
    if rng is None:
        rng = np.random.default_rng()

    vals = []
    for _ in range(num_loops):
        loop = random_loop_with_rng(rng, loop_length, N)
        vals.append(loop_holonomy(loop, N))

    vals = np.array(vals)
    K = abs(np.mean(vals))
    R = np.mean(abs(1 - vals))
    return K, R


def loop_class_coherence_scan(
    N_values,
    num_loops=5000,
    loop_length=6,
    seed=777,
):
    """
    Vergleicht kompatible und frustrierte EABC-Loops.
    """
    rng = np.random.default_rng(seed)
    rows = []

    for N in N_values:
        if N % 4 != 0:
            continue

        compatible_vals = []
        frustrated_vals = []
        for _ in range(num_loops):
            loop_c = generate_compatible_loop(rng, length=loop_length, N=N)
            loop_f = generate_frustrated_loop(rng, length=loop_length, N=N)
            compatible_vals.append(loop_holonomy(loop_c, N))
            frustrated_vals.append(loop_holonomy(loop_f, N))

        compatible_vals = np.array(compatible_vals)
        frustrated_vals = np.array(frustrated_vals)

        c_K = abs(np.mean(compatible_vals))
        c_R = np.mean(abs(1 - compatible_vals))
        c_ok = np.mean([is_compatible_loop(generate_compatible_loop(rng, loop_length, N), N) for _ in range(200)])

        f_K = abs(np.mean(frustrated_vals))
        f_R = np.mean(abs(1 - frustrated_vals))
        f_omega_neq_1 = np.mean(np.abs(1 - frustrated_vals) > 1e-12)

        rows.append(
            {
                "N": N,
                "compatible_K": float(c_K),
                "compatible_R": float(c_R),
                "compatible_valid_rate": float(c_ok),
                "frustrated_K": float(f_K),
                "frustrated_R": float(f_R),
                "frustrated_omega_neq_1_rate": float(f_omega_neq_1),
            }
        )

    return rows


def save_loop_class_plot(rows, filename="eabc_loop_classes.png"):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("\nPlot (Loop-Klassen) übersprungen: matplotlib nicht verfügbar.", exc)
        return None

    if not rows:
        return None

    n_values = np.array([r["N"] for r in rows], dtype=int)
    c_k = np.array([r["compatible_K"] for r in rows], dtype=float)
    f_k = np.array([r["frustrated_K"] for r in rows], dtype=float)
    c_r = np.array([r["compatible_R"] for r in rows], dtype=float)
    f_r = np.array([r["frustrated_R"] for r in rows], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    axes[0].plot(n_values, c_k, "o-", label="kompatibel")
    axes[0].plot(n_values, f_k, "o-", label="frustriert")
    axes[0].set_title("Kohärenz K nach Loop-Klasse")
    axes[0].set_xlabel("N")
    axes[0].set_ylabel("K = |mean(Omega)|")
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    axes[1].plot(n_values, c_r, "o-", label="kompatibel")
    axes[1].plot(n_values, f_r, "o-", label="frustriert")
    axes[1].set_title("Krümmungsproxy R nach Loop-Klasse")
    axes[1].set_xlabel("N")
    axes[1].set_ylabel("R = mean(|1-Omega|)")
    axes[1].grid(alpha=0.3)
    axes[1].legend()

    out_path = Path(filename).resolve()
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


def fit_frustrated_scaling(rows):
    """
    Fit für K_frust(N):
    1) Potenzgesetz: K = A * N^{-alpha}
    2) Exponential: K = B * exp(-c N)
    Rückgabe inkl. R^2 im jeweiligen linearen Log-Raum.
    """
    n_values = np.array([r["N"] for r in rows], dtype=float)
    k_values = np.array([r["frustrated_K"] for r in rows], dtype=float)

    mask = k_values > 0
    n = n_values[mask]
    k = k_values[mask]

    log_n = np.log(n)
    log_k = np.log(k)

    # Potenzfit: log(K) = b0 + b1 log(N), alpha = -b1, A = exp(b0)
    p_slope, p_intercept = np.polyfit(log_n, log_k, 1)
    alpha = -float(p_slope)
    A = float(math.exp(p_intercept))
    log_k_pred_power = p_intercept + p_slope * log_n

    # Expfit: log(K) = e0 + e1 N, c = -e1, B = exp(e0)
    e_slope, e_intercept = np.polyfit(n, log_k, 1)
    c = -float(e_slope)
    B = float(math.exp(e_intercept))
    log_k_pred_exp = e_intercept + e_slope * n

    def r2(y, y_pred):
        ss_res = float(np.sum((y - y_pred) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        return 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    r2_power = r2(log_k, log_k_pred_power)
    r2_exp = r2(log_k, log_k_pred_exp)

    return {
        "n": n,
        "k": k,
        "power": {
            "A": A,
            "alpha": alpha,
            "r2_log": r2_power,
            "k_pred": np.exp(log_k_pred_power),
        },
        "exp": {
            "B": B,
            "c": c,
            "r2_log": r2_exp,
            "k_pred": np.exp(log_k_pred_exp),
        },
    }


def bootstrap_alpha_ci(rows, n_boot=1000, seed=424242):
    """
    Bootstrap für den Power-Law-Exponent alpha.
    Resampling mit Zurücklegen auf den (N, K_frust)-Paaren.
    """
    fit = fit_frustrated_scaling(rows)
    n = fit["n"]
    k = fit["k"]
    m = len(n)
    if m < 3:
        return None

    rng = np.random.default_rng(seed)
    alphas = []
    r2_values = []

    for _ in range(n_boot):
        idx = rng.integers(0, m, size=m)
        n_b = n[idx]
        k_b = k[idx]

        log_n = np.log(n_b)
        log_k = np.log(k_b)

        slope, intercept = np.polyfit(log_n, log_k, 1)
        alpha_b = -float(slope)
        log_pred = intercept + slope * log_n

        ss_res = float(np.sum((log_k - log_pred) ** 2))
        ss_tot = float(np.sum((log_k - np.mean(log_k)) ** 2))
        r2_b = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

        alphas.append(alpha_b)
        r2_values.append(r2_b)

    alphas = np.array(alphas, dtype=float)
    r2_values = np.array(r2_values, dtype=float)

    ci_low, ci_high = np.percentile(alphas, [2.5, 97.5])
    return {
        "n_boot": n_boot,
        "alpha_mean": float(np.mean(alphas)),
        "alpha_median": float(np.median(alphas)),
        "alpha_std": float(np.std(alphas, ddof=1)),
        "alpha_ci95_low": float(ci_low),
        "alpha_ci95_high": float(ci_high),
        "alpha_min": float(np.min(alphas)),
        "alpha_max": float(np.max(alphas)),
        "r2_mean": float(np.mean(r2_values)),
        "r2_median": float(np.median(r2_values)),
    }


def save_frustrated_fit_plot(fit, filename="eabc_frustrated_scaling_fit.png"):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("\nPlot (Skalierungsfit) übersprungen: matplotlib nicht verfügbar.", exc)
        return None

    n = fit["n"]
    k = fit["k"]
    k_pred_power = fit["power"]["k_pred"]
    k_pred_exp = fit["exp"]["k_pred"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    # Semilogy: gut für Exponentialtrend
    axes[0].semilogy(n, k, "o", label="Daten K_frust")
    axes[0].semilogy(n, k_pred_exp, "-", label="Exp-Fit")
    axes[0].semilogy(n, k_pred_power, "--", label="Power-Fit")
    axes[0].set_title("K_frust(N) (semilogy)")
    axes[0].set_xlabel("N")
    axes[0].set_ylabel("K_frust")
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    # Log-Log: gut für Potenztrend
    axes[1].loglog(n, k, "o", label="Daten K_frust")
    axes[1].loglog(n, k_pred_power, "-", label="Power-Fit")
    axes[1].loglog(n, k_pred_exp, "--", label="Exp-Fit")
    axes[1].set_title("K_frust(N) (log-log)")
    axes[1].set_xlabel("N")
    axes[1].set_ylabel("K_frust")
    axes[1].grid(alpha=0.3)
    axes[1].legend()

    out_path = Path(filename).resolve()
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


def collect_frustrated_nl_data(
    N_values,
    L_values,
    num_loops=4000,
    seed=999,
):
    """
    Erzeugt Datensatz K_frust(N, L) und R_frust(N, L).
    """
    rng = np.random.default_rng(seed)
    rows = []
    for L in L_values:
        for N in N_values:
            if N % 4 != 0:
                continue
            vals = []
            for _ in range(num_loops):
                loop_f = generate_frustrated_loop(rng, length=L, N=N)
                vals.append(loop_holonomy(loop_f, N))
            vals = np.array(vals)
            kf = abs(np.mean(vals))
            rf = np.mean(abs(1 - vals))
            rows.append(
                {
                    "N": int(N),
                    "L": int(L),
                    "NL": int(N * L),
                    "K_frust": float(kf),
                    "R_frust": float(rf),
                }
            )
    return rows


def first_primes(count, start_at=5):
    primes = []
    n = max(2, start_at)
    while len(primes) < count:
        is_prime = True
        limit = int(math.sqrt(n))
        for d in range(2, limit + 1):
            if n % d == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(n)
        n += 1
    return primes


def residue_to_channel_mod12(residue):
    mapping = {1: "E", 5: "A", 7: "B", 11: "C"}
    return mapping.get(residue % 12, "E")


def structured_channel_sequence(family, length):
    if family == "cycle":
        base = ["E", "A", "C", "B"]
        return [base[i % len(base)] for i in range(length)]
    if family == "counter_cycle":
        base = ["E", "B", "C", "A"]
        return [base[i % len(base)] for i in range(length)]
    if family == "prime_residues":
        ps = first_primes(length + 10, start_at=5)
        channels = [residue_to_channel_mod12(p % 12) for p in ps]
        return channels[:length]
    if family == "prime_plus4_centers":
        ps = first_primes(length + 10, start_at=5)
        channels = []
        for p in ps:
            channels.append(residue_to_channel_mod12(p % 12))
            if len(channels) >= length:
                break
            channels.append(residue_to_channel_mod12((p + 4) % 12))
            if len(channels) >= length:
                break
        return channels[:length]
    raise ValueError(f"Unbekannte Loop-Familie: {family}")


def structured_frustrated_loop(rng, family, length=6, N=120):
    """
    Strukturierter Channel-Pfad + gezielte Frustration.
    Die letzte Phase wird so gesetzt, dass die Gesamtholonomie nichttrivial bleibt.
    """
    if length < 2:
        raise ValueError("Loop-Länge muss >= 2 sein.")
    if N % 4 != 0:
        raise ValueError("N muss durch 4 teilbar sein.")

    channels = structured_channel_sequence(family, length)
    k_values = [int(rng.integers(0, N)) for _ in range(length)]

    quarter_unit = N // 4
    base_units = sum(k_values[:-1]) + quarter_unit * sum(channel_quarter_turns(c) for c in channels[:-1])
    last_units = quarter_unit * channel_quarter_turns(channels[-1])

    # Wähle last_k nahe Schließung, dann störender Shift delta_k != 0.
    close_k = int((-base_units - last_units) % N)
    delta_k = int(rng.integers(1, N))
    k_values[-1] = (close_k + delta_k) % N

    return list(zip(channels, k_values))


def collect_structured_frustrated_nl_data(
    family,
    N_values,
    L_values,
    num_loops=3000,
    seed=1001,
):
    rng = np.random.default_rng(seed)
    rows = []
    for L in L_values:
        for N in N_values:
            if N % 4 != 0:
                continue
            vals = []
            for _ in range(num_loops):
                loop_f = structured_frustrated_loop(rng, family=family, length=L, N=N)
                vals.append(loop_holonomy(loop_f, N))
            vals = np.array(vals)
            rows.append(
                {
                    "family": family,
                    "N": int(N),
                    "L": int(L),
                    "NL": int(N * L),
                    "K_frust": float(abs(np.mean(vals))),
                    "R_frust": float(np.mean(abs(1 - vals))),
                }
            )
    return rows


def fit_power_single_variable(x_values, y_values):
    """
    Fit y = C * x^{-p} über log(y)=a+b log(x), p=-b.
    """
    x = np.array(x_values, dtype=float)
    y = np.array(y_values, dtype=float)
    mask = (x > 0) & (y > 0)
    x = x[mask]
    y = y[mask]

    lx = np.log(x)
    ly = np.log(y)
    slope, intercept = np.polyfit(lx, ly, 1)
    p = -float(slope)
    C = float(math.exp(intercept))
    ly_pred = intercept + slope * lx
    y_pred = np.exp(ly_pred)
    ss_res = float(np.sum((ly - ly_pred) ** 2))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    r2_log = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return {
        "C": C,
        "p": p,
        "r2_log": r2_log,
        "x": x,
        "y": y,
        "y_pred": y_pred,
    }


def fit_frustrated_nl_models(rows):
    """
    Vergleicht drei 1D-Skalierungsachsen:
    K(N,L) gegen N, gegen L, gegen NL.
    """
    n_vals = [r["N"] for r in rows]
    l_vals = [r["L"] for r in rows]
    nl_vals = [r["NL"] for r in rows]
    k_vals = [r["K_frust"] for r in rows]

    fit_n = fit_power_single_variable(n_vals, k_vals)
    fit_l = fit_power_single_variable(l_vals, k_vals)
    fit_nl = fit_power_single_variable(nl_vals, k_vals)

    return {"N": fit_n, "L": fit_l, "NL": fit_nl}


def save_nl_scaling_plot(rows, fits, filename="eabc_frustrated_nl_scaling.png"):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("\nPlot (N/L/NL) übersprungen: matplotlib nicht verfügbar.", exc)
        return None

    n = np.array([r["N"] for r in rows], dtype=float)
    l = np.array([r["L"] for r in rows], dtype=float)
    nl = np.array([r["NL"] for r in rows], dtype=float)
    k = np.array([r["K_frust"] for r in rows], dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))

    # K vs N
    order = np.argsort(n)
    axes[0].loglog(n, k, "o", alpha=0.7, label="Daten")
    axes[0].loglog(n[order], fits["N"]["y_pred"][np.argsort(fits["N"]["x"])], "-", label="Fit")
    axes[0].set_title("K_frust vs N")
    axes[0].set_xlabel("N")
    axes[0].set_ylabel("K_frust")
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    # K vs L
    axes[1].loglog(l, k, "o", alpha=0.7, label="Daten")
    axes[1].loglog(
        fits["L"]["x"][np.argsort(fits["L"]["x"])],
        fits["L"]["y_pred"][np.argsort(fits["L"]["x"])],
        "-",
        label="Fit",
    )
    axes[1].set_title("K_frust vs L")
    axes[1].set_xlabel("L")
    axes[1].set_ylabel("K_frust")
    axes[1].grid(alpha=0.3)
    axes[1].legend()

    # K vs NL
    axes[2].loglog(nl, k, "o", alpha=0.7, label="Daten")
    axes[2].loglog(
        fits["NL"]["x"][np.argsort(fits["NL"]["x"])],
        fits["NL"]["y_pred"][np.argsort(fits["NL"]["x"])],
        "-",
        label="Fit",
    )
    axes[2].set_title("K_frust vs N*L")
    axes[2].set_xlabel("N*L")
    axes[2].set_ylabel("K_frust")
    axes[2].grid(alpha=0.3)
    axes[2].legend()

    out_path = Path(filename).resolve()
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


def save_structured_family_scaling_plot(all_rows, family_fits, filename="eabc_structured_family_scaling.png"):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("\nPlot (strukturierte Familien) übersprungen: matplotlib nicht verfügbar.", exc)
        return None

    families = list(family_fits.keys())
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes = axes.flatten()

    for i, family in enumerate(families[:4]):
        ax = axes[i]
        rows = [r for r in all_rows if r["family"] == family]
        x = np.array([r["N"] for r in rows], dtype=float)
        y = np.array([r["K_frust"] for r in rows], dtype=float)
        fit = family_fits[family]["N"]
        order = np.argsort(fit["x"])

        ax.loglog(x, y, "o", alpha=0.6, label="Daten")
        ax.loglog(fit["x"][order], fit["y_pred"][order], "-", label="Fit vs N")
        ax.set_title(f"{family} (R2_N={fit['r2_log']:.3f})")
        ax.set_xlabel("N")
        ax.set_ylabel("K_frust")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    out_path = Path(filename).resolve()
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


def generate_open_string(family, length, N, rng):
    channels = []
    if family == "random":
        base = list(EABC.keys())
        channels = [rng.choice(base) for _ in range(length)]
    else:
        channels = structured_channel_sequence(family, length)
    k_values = [int(rng.integers(0, N)) for _ in range(length)]
    return list(zip(channels, k_values))


def generate_closed_string(length, N, rng):
    """
    Erzeugt einen string-artigen Pfad mit exakter Endschließung:
    Gruppenprodukt E und Gesamtphase 0 mod 2*pi.
    """
    if length <= 0:
        return []
    if N % 4 != 0:
        raise ValueError("Für geschlossene Strings muss N durch 4 teilbar sein.")
    if length == 1:
        return [("E", 0)]

    channels = list(EABC.keys())
    prefix_channels = [rng.choice(channels) for _ in range(length - 1)]
    prefix_k = [int(rng.integers(0, N)) for _ in range(length - 1)]

    prod_prefix = "E"
    for c in prefix_channels:
        prod_prefix = eabc_product(prod_prefix, c)
    last_channel = prod_prefix

    quarter_unit = N // 4
    prefix_phase_units = sum(prefix_k) + quarter_unit * sum(channel_quarter_turns(c) for c in prefix_channels)
    last_phase_units = quarter_unit * channel_quarter_turns(last_channel)
    last_k = int((-prefix_phase_units - last_phase_units) % N)

    return list(zip(prefix_channels, prefix_k)) + [(last_channel, last_k)]


def endpoint_and_cumulative_energy(path, N):
    """
    Endpunktenergie: |Omega_d - 1|^2
    Kumulative Energie: Summe ueber Prefix-Defekte.
    """
    omega = 1 + 0j
    cumulative = 0.0
    for channel, k in path:
        omega *= eabc_phase_factor(channel, k, N)
        cumulative += abs(omega - 1) ** 2
    endpoint = abs(omega - 1) ** 2
    return float(endpoint), float(cumulative)


def prefix_energies(path, N):
    omega = 1 + 0j
    endpoint_vals = np.zeros(len(path), dtype=float)
    cumulative_vals = np.zeros(len(path), dtype=float)
    run = 0.0
    for i, (channel, k) in enumerate(path):
        omega *= eabc_phase_factor(channel, k, N)
        e = float(abs(omega - 1) ** 2)
        run += e
        endpoint_vals[i] = e
        cumulative_vals[i] = run
    return endpoint_vals, cumulative_vals


def measure_open_strings_by_prefix(family, max_d, N, num_samples=400, seed=2027):
    rng = np.random.default_rng(seed)
    sum_endpoint = np.zeros(max_d, dtype=float)
    sum_endpoint_sq = np.zeros(max_d, dtype=float)
    sum_cum = np.zeros(max_d, dtype=float)
    sum_cum_sq = np.zeros(max_d, dtype=float)

    for _ in range(num_samples):
        path = generate_open_string(family, max_d, N, rng)
        endpoint_vals, cumulative_vals = prefix_energies(path, N)
        sum_endpoint += endpoint_vals
        sum_endpoint_sq += endpoint_vals ** 2
        sum_cum += cumulative_vals
        sum_cum_sq += cumulative_vals ** 2

    mean_endpoint = sum_endpoint / num_samples
    var_endpoint = np.maximum(0.0, sum_endpoint_sq / num_samples - mean_endpoint**2)
    mean_cum = sum_cum / num_samples
    var_cum = np.maximum(0.0, sum_cum_sq / num_samples - mean_cum**2)

    return {
        "d": np.arange(1, max_d + 1, dtype=int),
        "mean_endpoint": mean_endpoint,
        "var_endpoint": var_endpoint,
        "mean_cumulative": mean_cum,
        "var_cumulative": var_cum,
        "num_samples": num_samples,
        "family": family,
    }


def measure_closed_strings(max_d, N, num_samples=60, seed=2031):
    rng = np.random.default_rng(seed)
    mean_endpoint = np.zeros(max_d, dtype=float)
    var_endpoint = np.zeros(max_d, dtype=float)
    mean_cum = np.zeros(max_d, dtype=float)
    var_cum = np.zeros(max_d, dtype=float)

    for d in range(1, max_d + 1):
        endpoints = []
        cumulatives = []
        for _ in range(num_samples):
            path = generate_closed_string(d, N, rng)
            e_end, e_cum = endpoint_and_cumulative_energy(path, N)
            endpoints.append(e_end)
            cumulatives.append(e_cum)
        endpoints = np.array(endpoints, dtype=float)
        cumulatives = np.array(cumulatives, dtype=float)
        mean_endpoint[d - 1] = float(np.mean(endpoints))
        var_endpoint[d - 1] = float(np.var(endpoints))
        mean_cum[d - 1] = float(np.mean(cumulatives))
        var_cum[d - 1] = float(np.var(cumulatives))

    return {
        "d": np.arange(1, max_d + 1, dtype=int),
        "mean_endpoint": mean_endpoint,
        "var_endpoint": var_endpoint,
        "mean_cumulative": mean_cum,
        "var_cumulative": var_cum,
        "num_samples": num_samples,
        "family": "closed",
    }


def linear_fit_with_r2(x_values, y_values):
    x = np.array(x_values, dtype=float)
    y = np.array(y_values, dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": float(r2),
        "x": x,
        "y_pred": y_pred,
    }


def save_confinement_plot(results, filename="eabc_confinement_scan.png"):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("\nPlot (Confinement) übersprungen: matplotlib nicht verfügbar.", exc)
        return None

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # 1) Mittlere kumulative Energie
    ax = axes[0, 0]
    for label, data in results.items():
        ax.plot(data["d"], data["mean_cumulative"], label=label, linewidth=1.1)
    ax.set_title("Mittlere kumulative Frustrationsenergie E(d)")
    ax.set_xlabel("d")
    ax.set_ylabel("E_cum(d)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

    # 2) Endpunktenergie
    ax = axes[0, 1]
    for label, data in results.items():
        ax.plot(data["d"], data["mean_endpoint"], label=label, linewidth=1.1)
    ax.set_title("Mittlere Endpunktenergie |Omega_d-1|^2")
    ax.set_xlabel("d")
    ax.set_ylabel("E_end(d)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

    # 3) Varianz kumulative Energie
    ax = axes[1, 0]
    for label, data in results.items():
        ax.plot(data["d"], data["var_cumulative"], label=label, linewidth=1.1)
    ax.set_title("Varianz von E_cum(d)")
    ax.set_xlabel("d")
    ax.set_ylabel("Var[E_cum(d)]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

    # 4) Relative Trennung offen vs geschlossen (random)
    ax = axes[1, 1]
    if "open_random" in results and "closed" in results:
        d = results["open_random"]["d"]
        ratio = results["open_random"]["mean_cumulative"] / np.maximum(
            1e-12, results["closed"]["mean_cumulative"]
        )
        ax.plot(d, ratio, color="tab:red", linewidth=1.3)
        ax.set_title("Verhältnis E_open_random / E_closed")
        ax.set_xlabel("d")
        ax.set_ylabel("Ratio")
        ax.grid(alpha=0.3)
    else:
        ax.text(0.15, 0.5, "Ratio-Panel nicht verfügbar", transform=ax.transAxes)
        ax.axis("off")

    out_path = Path(filename).resolve()
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


def save_confinement_delta_plot(delta_results, filename="eabc_confinement_delta.png"):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("\nPlot (Confinement-Delta) übersprungen: matplotlib nicht verfügbar.", exc)
        return None

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    # Delta kumulative Energie
    ax = axes[0]
    for label, data in delta_results.items():
        ax.plot(data["d"], data["delta_cumulative"], linewidth=1.2, label=label)
    ax.set_title("Delta kumulative Energie")
    ax.set_xlabel("d")
    ax.set_ylabel("E_open_cum(d) - E_closed_cum(d)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    # Delta Endpunktenergie
    ax = axes[1]
    for label, data in delta_results.items():
        ax.plot(data["d"], data["delta_endpoint"], linewidth=1.2, label=label)
    ax.set_title("Delta Endpunktenergie")
    ax.set_xlabel("d")
    ax.set_ylabel("E_open_end(d) - E_closed_end(d)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    out_path = Path(filename).resolve()
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


def measure_open_string_endpoint_stats(family, max_d, N, num_samples=300, seed=2041, tol=1e-12):
    """
    Für offene Strings:
    - mean_endpoint(d) = E[|Omega_d-1|^2]
    - p_nontrivial(d) = P(|Omega_d-1| > tol)
    """
    rng = np.random.default_rng(seed)
    sum_endpoint = np.zeros(max_d, dtype=float)
    nontrivial_count = np.zeros(max_d, dtype=float)

    for _ in range(num_samples):
        path = generate_open_string(family, max_d, N, rng)
        endpoint_vals, _ = prefix_energies(path, N)
        sum_endpoint += endpoint_vals
        nontrivial_count += (endpoint_vals > tol).astype(float)

    mean_endpoint = sum_endpoint / num_samples
    p_nontrivial = nontrivial_count / num_samples
    return {
        "d": np.arange(1, max_d + 1, dtype=int),
        "mean_endpoint": mean_endpoint,
        "p_nontrivial": p_nontrivial,
        "num_samples": num_samples,
        "family": family,
    }


def measure_closed_string_endpoint_stats(max_d, N, num_samples=80, seed=2042, tol=1e-12):
    rng = np.random.default_rng(seed)
    mean_endpoint = np.zeros(max_d, dtype=float)
    p_nontrivial = np.zeros(max_d, dtype=float)

    for d in range(1, max_d + 1):
        endpoint_vals = []
        nontrivial = 0
        for _ in range(num_samples):
            path = generate_closed_string(d, N, rng)
            e_end, _ = endpoint_and_cumulative_energy(path, N)
            endpoint_vals.append(e_end)
            if e_end > tol:
                nontrivial += 1
        endpoint_vals = np.array(endpoint_vals, dtype=float)
        mean_endpoint[d - 1] = float(np.mean(endpoint_vals))
        p_nontrivial[d - 1] = nontrivial / num_samples

    return {
        "d": np.arange(1, max_d + 1, dtype=int),
        "mean_endpoint": mean_endpoint,
        "p_nontrivial": p_nontrivial,
        "num_samples": num_samples,
        "family": "closed",
    }


def save_confinement_string_energy_plot(results, filename="eabc_confinement_string_energy.png"):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("\nPlot (String-Energie) übersprungen: matplotlib nicht verfügbar.", exc)
        return None

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    # Hard string energy: d * 1_{Omega!=1}
    ax = axes[0]
    for label, data in results.items():
        d = data["d"]
        e_hard = d * data["p_nontrivial"]
        ax.plot(d, e_hard, linewidth=1.2, label=label)
    ax.set_title("Hard String-Energie d * P(Omega!=1)")
    ax.set_xlabel("d")
    ax.set_ylabel("E_hard(d)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    # Soft string energy: d * E[|Omega-1|^2]
    ax = axes[1]
    for label, data in results.items():
        d = data["d"]
        e_soft = d * data["mean_endpoint"]
        ax.plot(d, e_soft, linewidth=1.2, label=label)
    ax.set_title("Soft String-Energie d * E[|Omega-1|^2]")
    ax.set_xlabel("d")
    ax.set_ylabel("E_soft(d)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    out_path = Path(filename).resolve()
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


def holonomy_coherence_stats(
    N,
    num_batches=25,
    loops_per_batch=4000,
    loop_length=6,
    seed=12345,
):
    """
    Robustere Schätzung über mehrere unabhängige Batches.
    Liefert Mittelwert, Standardabweichung, Standardfehler und 95%-CI.
    """
    rng = np.random.default_rng(seed)
    k_values = []
    r_values = []

    for _ in range(num_batches):
        K, R = holonomy_coherence(
            num_loops=loops_per_batch,
            loop_length=loop_length,
            N=N,
            rng=rng,
        )
        k_values.append(K)
        r_values.append(R)

    k_values = np.array(k_values, dtype=float)
    r_values = np.array(r_values, dtype=float)

    k_mean = float(np.mean(k_values))
    k_std = float(np.std(k_values, ddof=1))
    k_se = k_std / math.sqrt(num_batches)
    k_ci95 = 1.96 * k_se

    r_mean = float(np.mean(r_values))
    r_std = float(np.std(r_values, ddof=1))
    r_se = r_std / math.sqrt(num_batches)
    r_ci95 = 1.96 * r_se

    return {
        "K_mean": k_mean,
        "K_std": k_std,
        "K_se": k_se,
        "K_ci95": k_ci95,
        "R_mean": r_mean,
        "R_std": r_std,
        "R_se": r_se,
        "R_ci95": r_ci95,
        "num_batches": num_batches,
        "loops_per_batch": loops_per_batch,
        "total_loops": num_batches * loops_per_batch,
    }


def is_eabc_compatible_n(N):
    """
    Für die CHSH-Standardwinkel (0, pi/2, +/-pi/4) ist N % 8 == 0
    die natürliche Kompatibilitätsbedingung.
    """
    return N % 8 == 0


def chsh_scan(n_min=8, n_max=512):
    s_qm = chsh_value()
    rows = []
    for N in range(n_min, n_max + 1):
        s_n = chsh_value(N)
        rows.append(
            {
                "N": N,
                "S": s_n,
                "delta": abs(s_qm - s_n),
                "compatible": is_eabc_compatible_n(N),
            }
        )
    return rows, s_qm


def save_chsh_deviation_plot(rows, s_qm, filename="chsh_deviation_vs_N.png"):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("\nPlot übersprungen: matplotlib nicht verfügbar.", exc)
        return None

    n_values = np.array([r["N"] for r in rows], dtype=int)
    deviations = np.array([r["delta"] for r in rows], dtype=float)
    compatible_mask = np.array([r["compatible"] for r in rows], dtype=bool)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(n_values, deviations, linewidth=1.5, color="tab:blue", label="|S(N)-S_QM|")
    ax.scatter(
        n_values[compatible_mask],
        deviations[compatible_mask],
        s=14,
        color="tab:green",
        label="EABC-kompatibel (N % 8 == 0)",
        zorder=3,
    )
    ax.set_title("CHSH-Abweichung vom QM-Wert")
    ax.set_xlabel("N (Phasengittergröße)")
    ax.set_ylabel("|S(N) - 2*sqrt(2)|")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right")

    out_path = Path(filename).resolve()
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


def summarize_incompatible_mod_classes(rows):
    """
    Nur nicht-kompatible N (N % 8 != 0):
    gruppiert nach Restklasse mod 8 und liefert robuste Kennzahlen.
    """
    groups = {r: [] for r in range(1, 8)}
    for row in rows:
        if row["compatible"]:
            continue
        residue = row["N"] % 8
        groups[residue].append(row)

    summary = []
    for residue in range(1, 8):
        vals = groups[residue]
        deltas = np.array([v["delta"] for v in vals], dtype=float)
        if len(deltas) == 0:
            continue
        max_idx = int(np.argmax(deltas))
        max_item = vals[max_idx]
        summary.append(
            {
                "residue": residue,
                "count": len(vals),
                "mean_delta": float(np.mean(deltas)),
                "median_delta": float(np.median(deltas)),
                "std_delta": float(np.std(deltas, ddof=1)) if len(deltas) > 1 else 0.0,
                "max_delta": float(deltas[max_idx]),
                "N_at_max": int(max_item["N"]),
            }
        )
    summary.sort(key=lambda x: x["mean_delta"], reverse=True)
    return summary


def save_incompatible_mod8_plot(rows, filename="chsh_incompatible_mod8.png"):
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print("\nPlot (mod8) übersprungen: matplotlib nicht verfügbar.", exc)
        return None

    incompatible = [r for r in rows if not r["compatible"]]
    residues = sorted(set(r["N"] % 8 for r in incompatible))

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

    # Links: Streuung delta vs N, nach mod-8-Farbe.
    for residue in residues:
        sub = [r for r in incompatible if (r["N"] % 8) == residue]
        n_values = np.array([r["N"] for r in sub], dtype=int)
        deltas = np.array([r["delta"] for r in sub], dtype=float)
        axes[0].scatter(n_values, deltas, s=12, alpha=0.75, label=f"N mod 8 = {residue}")
    axes[0].set_title("Nicht-kompatible N: |S(N)-S_QM|")
    axes[0].set_xlabel("N")
    axes[0].set_ylabel("|S(N)-2*sqrt(2)|")
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=8, ncol=2)

    # Rechts: mittlere Abweichung pro Restklasse.
    means = []
    for residue in residues:
        sub = [r["delta"] for r in incompatible if (r["N"] % 8) == residue]
        means.append(float(np.mean(sub)))
    axes[1].bar([str(r) for r in residues], means, color="tab:orange")
    axes[1].set_title("Mittelwert |delta| pro Restklasse")
    axes[1].set_xlabel("N mod 8")
    axes[1].set_ylabel("mean(|S(N)-2*sqrt(2)|)")
    axes[1].grid(alpha=0.3, axis="y")

    out_path = Path(filename).resolve()
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


print("\nRandomisierte EABC-Holonomie:")
for N in [12, 24, 60, 120, 360]:
    stats = holonomy_coherence_stats(
        N=N,
        num_batches=25,
        loops_per_batch=4000,
        loop_length=6,
        seed=2026 + N,
    )
    print(
        f"N={N:>4}  K={stats['K_mean']:.5f}±{stats['K_ci95']:.5f} (95%-CI)  "
        f"R={stats['R_mean']:.5f}±{stats['R_ci95']:.5f}  "
        f"loops={stats['total_loops']}"
    )


print("\nCHSH-Scan über viele N:")
scan_rows, s_qm = chsh_scan(n_min=8, n_max=512)
compatible_count = sum(1 for row in scan_rows if row["compatible"])

print(f"Referenz S_QM = {s_qm:.8f}")
print(f"Gescannt: {len(scan_rows)} N-Werte, davon kompatibel: {compatible_count}")
print("Ausgewählte Punkte (N, S(N), |delta|, kompatibel):")
for row in scan_rows:
    if row["N"] in [8, 12, 16, 24, 32, 64, 128, 256, 400, 512]:
        compat = "ja" if row["compatible"] else "nein"
        print(
            f"N={row['N']:>4}  S={row['S']:.8f}  "
            f"|delta|={row['delta']:.3e}  kompatibel={compat}"
        )

plot_path = save_chsh_deviation_plot(scan_rows, s_qm)
if plot_path is not None:
    print(f"CHSH-Abweichungsplot gespeichert unter: {plot_path}")

print("\nNicht-kompatible N (N mod 8 != 0):")
mod_summary = summarize_incompatible_mod_classes(scan_rows)
for item in mod_summary:
    print(
        f"mod8={item['residue']}  count={item['count']:>3}  "
        f"mean={item['mean_delta']:.5e}  median={item['median_delta']:.5e}  "
        f"std={item['std_delta']:.5e}  max={item['max_delta']:.5e} @N={item['N_at_max']}"
    )

mod_plot_path = save_incompatible_mod8_plot(scan_rows)
if mod_plot_path is not None:
    print(f"Modulo-Analyseplot gespeichert unter: {mod_plot_path}")

print("\nEABC-Loop-Klassen: kompatibel vs. frustriert")
loop_rows = loop_class_coherence_scan(
    N_values=[8, 12, 16, 20, 24, 28, 32, 40, 48, 64, 80, 96, 120, 160, 192, 256],
    num_loops=5000,
    loop_length=6,
    seed=20260528,
)
for row in loop_rows:
    print(
        f"N={row['N']:>4}  "
        f"K_komp={row['compatible_K']:.5f}  R_komp={row['compatible_R']:.5e}  "
        f"K_frust={row['frustrated_K']:.5f}  R_frust={row['frustrated_R']:.5f}  "
        f"valid_rate={row['compatible_valid_rate']:.3f}  "
        f"Omega!=1_rate={row['frustrated_omega_neq_1_rate']:.3f}"
    )

loop_plot_path = save_loop_class_plot(loop_rows)
if loop_plot_path is not None:
    print(f"Loop-Klassenplot gespeichert unter: {loop_plot_path}")

print("\nSkalierungsfit für K_frust(N):")
fit = fit_frustrated_scaling(loop_rows)
alpha = fit["power"]["alpha"]
A = fit["power"]["A"]
r2_power = fit["power"]["r2_log"]
c = fit["exp"]["c"]
B = fit["exp"]["B"]
r2_exp = fit["exp"]["r2_log"]

print(f"Power-Law:     K ~ A * N^(-alpha),  A={A:.6f}, alpha={alpha:.6f}, R2_log={r2_power:.6f}")
print(f"Exponential:   K ~ B * exp(-c*N),   B={B:.6f}, c={c:.6f}, R2_log={r2_exp:.6f}")
better = "Power-Law" if r2_power > r2_exp else "Exponential"
print(f"Besserer Fit (nach R2_log): {better}")

fit_plot_path = save_frustrated_fit_plot(fit)
if fit_plot_path is not None:
    print(f"Skalierungsfit-Plot gespeichert unter: {fit_plot_path}")

print("\nBootstrap-Stabilität für alpha (Power-Law):")
alpha_boot = bootstrap_alpha_ci(loop_rows, n_boot=1000, seed=20260528)
if alpha_boot is None:
    print("Zu wenige Datenpunkte für Bootstrap.")
else:
    print(
        f"n_boot={alpha_boot['n_boot']}  "
        f"alpha_mean={alpha_boot['alpha_mean']:.6f}  "
        f"alpha_median={alpha_boot['alpha_median']:.6f}  "
        f"alpha_std={alpha_boot['alpha_std']:.6f}"
    )
    print(
        f"95%-KI(alpha)=[{alpha_boot['alpha_ci95_low']:.6f}, {alpha_boot['alpha_ci95_high']:.6f}]  "
        f"alpha_min={alpha_boot['alpha_min']:.6f}  alpha_max={alpha_boot['alpha_max']:.6f}"
    )
    print(
        f"Bootstrap-R2_log: mean={alpha_boot['r2_mean']:.6f}  "
        f"median={alpha_boot['r2_median']:.6f}"
    )

print("\n2D-Skalierungstest K_frust(N, L):")
nl_num_loops = 8000
nl_rows = collect_frustrated_nl_data(
    N_values=[8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256],
    L_values=[4, 6, 8, 10, 12, 16],
    num_loops=nl_num_loops,
    seed=20260601,
)
fits_nl = fit_frustrated_nl_models(nl_rows)

print(
    f"Fit nur N:   K~C*N^-p,    C={fits_nl['N']['C']:.6f}, "
    f"p={fits_nl['N']['p']:.6f}, R2_log={fits_nl['N']['r2_log']:.6f}"
)
print(
    f"Fit nur L:   K~C*L^-p,    C={fits_nl['L']['C']:.6f}, "
    f"p={fits_nl['L']['p']:.6f}, R2_log={fits_nl['L']['r2_log']:.6f}"
)
print(
    f"Fit auf NL:  K~C*(NL)^-p, C={fits_nl['NL']['C']:.6f}, "
    f"p={fits_nl['NL']['p']:.6f}, R2_log={fits_nl['NL']['r2_log']:.6f}"
)

best_axis = max(["N", "L", "NL"], key=lambda key: fits_nl[key]["r2_log"])
print(f"Beste 1D-Skalierungsachse (R2_log): {best_axis}")

nl_plot_path = save_nl_scaling_plot(nl_rows, fits_nl)
if nl_plot_path is not None:
    print(f"N/L/NL-Skalierungsplot gespeichert unter: {nl_plot_path}")

print("\nStrukturierte Loops: Familienvergleich K_frust(N,L)")
families = ["cycle", "counter_cycle", "prime_residues", "prime_plus4_centers"]
structured_num_loops = 8000
all_structured_rows = []
family_fits = {}
for i, family in enumerate(families):
    rows_f = collect_structured_frustrated_nl_data(
        family=family,
        N_values=[8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256],
        L_values=[4, 6, 8, 10, 12, 16],
        num_loops=structured_num_loops,
        seed=20260610 + i,
    )
    all_structured_rows.extend(rows_f)
    fits_f = fit_frustrated_nl_models(rows_f)
    family_fits[family] = fits_f
    best_axis_f = max(["N", "L", "NL"], key=lambda key: fits_f[key]["r2_log"])
    print(
        f"{family:>20}: "
        f"R2_N={fits_f['N']['r2_log']:.3f}, "
        f"R2_L={fits_f['L']['r2_log']:.3f}, "
        f"R2_NL={fits_f['NL']['r2_log']:.3f}, "
        f"best={best_axis_f}"
    )

structured_plot_path = save_structured_family_scaling_plot(all_structured_rows, family_fits)
if structured_plot_path is not None:
    print(f"Strukturierte-Familien-Plot gespeichert unter: {structured_plot_path}")

print("\nHigh-Stat-Nachtest: counter_cycle")
counter_highstat_loops = 20000
counter_highstat_rows = collect_structured_frustrated_nl_data(
    family="counter_cycle",
    N_values=[8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256],
    L_values=[4, 6, 8, 10, 12, 16, 20, 24],
    num_loops=counter_highstat_loops,
    seed=20260699,
)
counter_highstat_fits = fit_frustrated_nl_models(counter_highstat_rows)
counter_best_axis = max(["N", "L", "NL"], key=lambda key: counter_highstat_fits[key]["r2_log"])
print(
    f"counter_cycle (high-stat): "
    f"R2_N={counter_highstat_fits['N']['r2_log']:.3f}, "
    f"R2_L={counter_highstat_fits['L']['r2_log']:.3f}, "
    f"R2_NL={counter_highstat_fits['NL']['r2_log']:.3f}, "
    f"best={counter_best_axis}, loops={counter_highstat_loops}"
)

print("\nEABC-Confinement-Modell: Stringscan d=1..1000")
max_d = 1000
confinement_N = 120
open_samples = 400
closed_samples = 60

confinement_results = {
    "open_random": measure_open_strings_by_prefix(
        family="random", max_d=max_d, N=confinement_N, num_samples=open_samples, seed=20261001
    ),
    "open_cycle": measure_open_strings_by_prefix(
        family="cycle", max_d=max_d, N=confinement_N, num_samples=open_samples, seed=20261002
    ),
    "open_counter_cycle": measure_open_strings_by_prefix(
        family="counter_cycle", max_d=max_d, N=confinement_N, num_samples=open_samples, seed=20261003
    ),
    "open_prime_residues": measure_open_strings_by_prefix(
        family="prime_residues", max_d=max_d, N=confinement_N, num_samples=open_samples, seed=20261004
    ),
    "open_prime_plus4": measure_open_strings_by_prefix(
        family="prime_plus4_centers", max_d=max_d, N=confinement_N, num_samples=open_samples, seed=20261005
    ),
    "closed": measure_closed_strings(
        max_d=max_d, N=confinement_N, num_samples=closed_samples, seed=20261010
    ),
}

print(
    f"Samples: open={open_samples} pro Familie, closed={closed_samples}, "
    f"N={confinement_N}, max_d={max_d}"
)

print("Linearfits E_cum(d)=a*d+b (ab d>=50):")
for label, data in confinement_results.items():
    d = data["d"]
    y = data["mean_cumulative"]
    mask = d >= 50
    fit = linear_fit_with_r2(d[mask], y[mask])
    print(
        f"{label:>18}: slope={fit['slope']:.6f}, intercept={fit['intercept']:.3f}, "
        f"R2={fit['r2']:.4f}, Var@d=1000={data['var_cumulative'][-1]:.3f}"
    )

confinement_plot_path = save_confinement_plot(confinement_results)
if confinement_plot_path is not None:
    print(f"Confinement-Plot gespeichert unter: {confinement_plot_path}")

print("\nConfinement-Differenztest: DeltaE(d)=E_open(d)-E_closed(d)")
closed_data = confinement_results["closed"]
delta_results = {}
for label, data in confinement_results.items():
    if label == "closed":
        continue
    delta_results[label] = {
        "d": data["d"],
        "delta_cumulative": data["mean_cumulative"] - closed_data["mean_cumulative"],
        "delta_endpoint": data["mean_endpoint"] - closed_data["mean_endpoint"],
    }

print("Linearfits fuer Delta kumulative Energie (ab d>=50):")
for label, data in delta_results.items():
    mask = data["d"] >= 50
    fit = linear_fit_with_r2(data["d"][mask], data["delta_cumulative"][mask])
    print(
        f"{label:>18}: slope={fit['slope']:.6f}, intercept={fit['intercept']:.3f}, R2={fit['r2']:.4f}"
    )

print("Linearfits fuer Delta Endpunktenergie (ab d>=50):")
for label, data in delta_results.items():
    mask = data["d"] >= 50
    fit = linear_fit_with_r2(data["d"][mask], data["delta_endpoint"][mask])
    print(
        f"{label:>18}: slope={fit['slope']:.6f}, intercept={fit['intercept']:.3f}, R2={fit['r2']:.4f}"
    )

delta_plot_path = save_confinement_delta_plot(delta_results)
if delta_plot_path is not None:
    print(f"Confinement-Delta-Plot gespeichert unter: {delta_plot_path}")

print("\nAlternative String-Energie (offene Randladung):")
string_max_d = 1000
string_N = 120
string_open_samples = 300
string_closed_samples = 80

string_results = {
    "open_random": measure_open_string_endpoint_stats(
        "random", max_d=string_max_d, N=string_N, num_samples=string_open_samples, seed=20261101
    ),
    "open_cycle": measure_open_string_endpoint_stats(
        "cycle", max_d=string_max_d, N=string_N, num_samples=string_open_samples, seed=20261102
    ),
    "open_counter_cycle": measure_open_string_endpoint_stats(
        "counter_cycle", max_d=string_max_d, N=string_N, num_samples=string_open_samples, seed=20261103
    ),
    "open_prime_residues": measure_open_string_endpoint_stats(
        "prime_residues", max_d=string_max_d, N=string_N, num_samples=string_open_samples, seed=20261104
    ),
    "open_prime_plus4": measure_open_string_endpoint_stats(
        "prime_plus4_centers", max_d=string_max_d, N=string_N, num_samples=string_open_samples, seed=20261105
    ),
    "closed": measure_closed_string_endpoint_stats(
        max_d=string_max_d, N=string_N, num_samples=string_closed_samples, seed=20261110
    ),
}

print(
    f"Samples (String-Energie): open={string_open_samples} pro Familie, "
    f"closed={string_closed_samples}, N={string_N}, max_d={string_max_d}"
)
print("Linearfits ab d>=50 fuer E_hard(d)=d*P(Omega!=1) und E_soft(d)=d*E[|Omega-1|^2]:")
for label, data in string_results.items():
    d = data["d"]
    mask = d >= 50
    e_hard = d * data["p_nontrivial"]
    e_soft = d * data["mean_endpoint"]
    fit_hard = linear_fit_with_r2(d[mask], e_hard[mask])
    fit_soft = linear_fit_with_r2(d[mask], e_soft[mask])
    print(
        f"{label:>18}: "
        f"hard_slope={fit_hard['slope']:.6f}, hard_R2={fit_hard['r2']:.4f}; "
        f"soft_slope={fit_soft['slope']:.6f}, soft_R2={fit_soft['r2']:.4f}"
    )

string_plot_path = save_confinement_string_energy_plot(string_results)
if string_plot_path is not None:
    print(f"String-Energie-Plot gespeichert unter: {string_plot_path}")


print("\nFertig.")