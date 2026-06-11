import numpy as np
from dataclasses import dataclass
from typing import List, Dict, Tuple
import random

# ------------------------------------------------------------
# Zustände
# ------------------------------------------------------------

STATES = [
    "EABC",                  # S1
    "ABCE",                  # S2
    "ABBA",                  # S3
    "BAAB",                  # S4
    "AABB",                  # S5
    "ABAB",                  # S6
    "C+Delta3(EAB)",         # S7
    "E+Delta3(ABC)",         # S8
    "Delta2+Delta2",         # S9
    "C"                      # S10 absorbierend
]

STATE_INDEX = {s: i for i, s in enumerate(STATES)}
ABSORBING_INDEX = STATE_INDEX["C"]


# ------------------------------------------------------------
# Hilfsfunktionen
# ------------------------------------------------------------

def sigmoid(x: float) -> float:
    return 1.0 / (1.0 + np.exp(-x))


@dataclass
class EnergyParams:
    theta_p: float = 3.0   # Schwelle für Pinski-Kanäle
    kappa_p: float = 2.0   # Steilheit
    theta_d: float = 5.0   # Schwelle für direkte Detektorkopplung
    kappa_d: float = 1.5   # Steilheit


def alpha_pinski(sqrt_s: float, params: EnergyParams) -> float:
    return sigmoid(params.kappa_p * (sqrt_s - params.theta_p))


def gamma_detector(sqrt_s: float, params: EnergyParams) -> float:
    return sigmoid(params.kappa_d * (sqrt_s - params.theta_d))


def beta_neutral(sqrt_s: float, params: EnergyParams) -> float:
    return 1.0 - alpha_pinski(sqrt_s, params)


def normalize_row(row: np.ndarray) -> np.ndarray:
    total = row.sum()
    if total <= 0:
        raise ValueError("Zeile kann nicht normiert werden: Summe <= 0.")
    return row / total


# ------------------------------------------------------------
# Energieabhängige Matrix
# ------------------------------------------------------------

def build_transition_matrix(sqrt_s: float, params: EnergyParams) -> np.ndarray:
    a = alpha_pinski(sqrt_s, params)   # Pinski-Aktivierung
    b = beta_neutral(sqrt_s, params)   # neutrale Bindungsneigung
    g = gamma_detector(sqrt_s, params) # direkte Detektorkopplung

    P = np.zeros((10, 10), dtype=float)

    # --------------------------------------------------------
    # S1 = EABC
    # geladener Vierer: Pinski-Kanäle öffnen mit Energie
    # --------------------------------------------------------
    row = np.array([
        0.00,        # EABC
        0.10 * b,    # ABCE
        0.00,        # ABBA
        0.00,        # BAAB
        0.10 * b,    # AABB
        0.10 * b,    # ABAB
        0.45 * a,    # C+Delta3(EAB)
        0.25 * a,    # E+Delta3(ABC)
        0.10 * b,    # Delta2+Delta2
        0.10 * g     # C
    ])
    P[0] = normalize_row(row)

    # --------------------------------------------------------
    # S2 = ABCE
    # ähnlich wie S1, aber etwas stärker auf E+Delta3
    # --------------------------------------------------------
    row = np.array([
        0.10 * b,    # EABC
        0.00,
        0.00,
        0.00,
        0.05 * b,
        0.10 * b,
        0.25 * a,
        0.35 * a,
        0.05 * b,
        0.10 * g
    ])
    P[1] = normalize_row(row)

    # --------------------------------------------------------
    # S3 = ABBA
    # neutraler Spiegelzustand: bei kleiner Energie bevorzugt 2+2
    # --------------------------------------------------------
    row = np.array([
        0.00,
        0.00,
        0.10 * b,
        0.15 * b,
        0.10 * b,
        0.00,
        0.10 * a,
        0.00,
        0.55 * b,
        0.10 * g
    ])
    P[2] = normalize_row(row)

    # --------------------------------------------------------
    # S4 = BAAB
    # neutraler Spiegelzustand
    # --------------------------------------------------------
    row = np.array([
        0.00,
        0.00,
        0.15 * b,
        0.10 * b,
        0.10 * b,
        0.00,
        0.10 * a,
        0.00,
        0.55 * b,
        0.10 * g
    ])
    P[3] = normalize_row(row)

    # --------------------------------------------------------
    # S5 = AABB
    # Übergangszustand / Verzweigungspunkt
    # --------------------------------------------------------
    row = np.array([
        0.05 * a,    # EABC
        0.00,
        0.30 * b,    # ABBA
        0.10 * b,    # BAAB
        0.05 * b,    # AABB
        0.20,        # ABAB (immer relevant)
        0.20 * a,    # C+Delta3(EAB)
        0.05 * a,    # E+Delta3(ABC)
        0.20 * b,    # Delta2+Delta2
        0.05 * g     # C
    ])
    P[4] = normalize_row(row)

    # --------------------------------------------------------
    # S6 = ABAB
    # intermediär: teils gerichtet, teils neutral
    # --------------------------------------------------------
    row = np.array([
        0.05 * a,    # EABC
        0.10 * a,    # ABCE
        0.05 * b,    # ABBA
        0.00,
        0.15 * b,    # AABB
        0.05,
        0.30 * a,    # C+Delta3(EAB)
        0.10 * a,    # E+Delta3(ABC)
        0.05 * b,    # Delta2+Delta2
        0.15 * g     # C
    ])
    P[5] = normalize_row(row)

    # --------------------------------------------------------
    # S7 = C+Delta3(EAB)
    # relaxiert stark in C
    # --------------------------------------------------------
    row = np.array([
        0.00,
        0.00,
        0.00,
        0.00,
        0.00,
        0.00,
        0.15 * (1.0 - g),  # metastabil
        0.05 * (1.0 - g),  # Wechsel in anderen 1+3-Zustand
        0.10,
        0.70 + 0.15 * g    # C
    ])
    P[6] = normalize_row(row)

    # --------------------------------------------------------
    # S8 = E+Delta3(ABC)
    # relaxiert ebenfalls stark in C
    # --------------------------------------------------------
    row = np.array([
        0.00,
        0.00,
        0.00,
        0.00,
        0.00,
        0.00,
        0.05 * (1.0 - g),
        0.15 * (1.0 - g),
        0.10,
        0.70 + 0.15 * g
    ])
    P[7] = normalize_row(row)

    # --------------------------------------------------------
    # S9 = Delta2+Delta2
    # neutraler 2+2-Zerfallszustand
    # --------------------------------------------------------
    row = np.array([
        0.00,
        0.00,
        0.10 * b,
        0.10 * b,
        0.00,
        0.00,
        0.00,
        0.00,
        0.40 * b + 0.10,  # metastabil neutral
        0.40 * g + 0.10   # C
    ])
    P[8] = normalize_row(row)

    # --------------------------------------------------------
    # S10 = C absorbierend
    # --------------------------------------------------------
    P[9, 9] = 1.0

    return P


# ------------------------------------------------------------
# Ausgabe-Funktionen
# ------------------------------------------------------------

def pretty_print_matrix(P: np.ndarray) -> None:
    header = " " * 16 + " ".join(f"{s[:12]:>12}" for s in STATES)
    print(header)
    for i, row in enumerate(P):
        print(f"{STATES[i][:15]:>15} " + " ".join(f"{x:12.3f}" for x in row))


def simulate_path(
    P: np.ndarray,
    start_state: str,
    max_steps: int = 50,
    rng_seed: int = None
) -> List[str]:
    if rng_seed is not None:
        random.seed(rng_seed)

    current = STATE_INDEX[start_state]
    path = [STATES[current]]

    for _ in range(max_steps):
        if current == ABSORBING_INDEX:
            break
        probs = P[current]
        current = random.choices(range(len(STATES)), weights=probs, k=1)[0]
        path.append(STATES[current])

    return path


def estimate_absorption_stats(
    P: np.ndarray,
    start_state: str,
    n_runs: int = 5000,
    max_steps: int = 200,
    rng_seed: int = 123
) -> Dict[str, float]:
    random.seed(rng_seed)
    start_idx = STATE_INDEX[start_state]

    absorbed = 0
    steps_to_absorb = []

    for _ in range(n_runs):
        current = start_idx
        for step in range(max_steps):
            if current == ABSORBING_INDEX:
                absorbed += 1
                steps_to_absorb.append(step)
                break
            current = random.choices(range(len(STATES)), weights=P[current], k=1)[0]

    absorption_prob = absorbed / n_runs
    mean_steps = float(np.mean(steps_to_absorb)) if steps_to_absorb else float("nan")

    return {
        "absorption_probability": absorption_prob,
        "mean_steps_to_absorption": mean_steps,
        "n_absorbed": absorbed,
        "n_runs": n_runs
    }


def transient_distribution(
    P: np.ndarray,
    start_state: str,
    n_steps: int = 10
) -> np.ndarray:
    v = np.zeros(len(STATES))
    v[STATE_INDEX[start_state]] = 1.0

    history = [v.copy()]
    for _ in range(n_steps):
        v = v @ P
        history.append(v.copy())

    return np.array(history)


# ------------------------------------------------------------
# Beispiel-Nutzung
# ------------------------------------------------------------

if __name__ == "__main__":
    params = EnergyParams(theta_p=3.0, kappa_p=2.0, theta_d=5.0, kappa_d=1.5)

    for sqrt_s in [1.5, 3.5, 6.0]:
        print("\n" + "=" * 80)
        print(f"Energie sqrt(s) = {sqrt_s}")
        P = build_transition_matrix(sqrt_s, params)

        print("\nÜbergangsmatrix:")
        pretty_print_matrix(P)

        print("\nBeispielpfad ab EABC:")
        path = simulate_path(P, start_state="EABC", max_steps=20, rng_seed=42)
        print(" -> ".join(path))

        print("\nAbsorptionsstatistik ab EABC:")
        stats = estimate_absorption_stats(P, start_state="EABC", n_runs=3000, max_steps=200)
        for k, v in stats.items():
            print(f"{k}: {v}")

        print("\nVerteilung über die ersten 6 Schritte ab ABBA:")
        hist = transient_distribution(P, start_state="ABBA", n_steps=6)
        for t, dist in enumerate(hist):
            summary = ", ".join(f"{STATES[i]}={dist[i]:.3f}" for i in range(len(STATES)) if dist[i] > 0.02)
            print(f"t={t}: {summary}")