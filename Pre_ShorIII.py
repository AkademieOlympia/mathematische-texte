import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import product
from typing import Optional, Dict, Tuple, Any

import matplotlib.pyplot as plt
import numpy as np

try:
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import squareform

    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False

# ============================================================
# Alphabet und absorbierende Kanäle
# ============================================================

ALPHABET = ("E", "A", "B", "C")
ABSORBING = ["A", "B", "C", "E", "0"]

WORDS = ["".join(p) for p in product(ALPHABET, repeat=4)]
ALL_STATES = WORDS + ABSORBING

WORD_INDEX = {w: i for i, w in enumerate(WORDS)}
STATE_INDEX = {s: i for i, s in enumerate(ALL_STATES)}
N_TOTAL = len(ALL_STATES)


def clip(x, lo=0.0, hi=1.0):
    return max(lo, min(hi, x))


def sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def profile_of_word(w: str):
    return tuple(sorted(Counter(w).values(), reverse=True))


def reverse_word(w: str):
    return w[::-1]


def swap_adjacent(w: str, i: int):
    s = list(w)
    s[i], s[i + 1] = s[i + 1], s[i]
    return "".join(s)


def all_local_neighbors(w: str):
    """Lokale Nachbarn via Swaps, Reverse und zyklische Rotation."""
    nbrs = set()
    for i in range(3):
        nbrs.add(swap_adjacent(w, i))
    nbrs.add(reverse_word(w))
    nbrs.add(w[1:] + w[0])
    nbrs.add(w[-1] + w[:-1])
    nbrs.discard(w)
    return sorted(nbrs)


CYCLIC_PLUS = {"EABC", "ABCE", "BCEA", "CEAB"}
CYCLIC_MINUS = {"EACB", "ACBE", "CBEA", "BEAC"}
LETTER_ORDER = {"E": 0, "A": 1, "B": 2, "C": 3}


def classify_22_word(w: str):
    if profile_of_word(w) != (2, 2):
        raise ValueError(f"{w} ist nicht vom Typ (2,2)")

    symbols = []
    seen = set()
    for ch in w:
        if ch not in seen:
            symbols.append(ch)
            seen.add(ch)
    x, y = symbols[0], symbols[1]

    block = {x + x + y + y, y + y + x + x}
    alt = {x + y + x + y, y + x + y + x}
    mirror = {x + y + y + x, y + x + x + y}

    if w in mirror:
        return "N"
    if w in alt:
        return "A"
    if w in block:
        return "B"
    raise RuntimeError(f"Unerwartete (2,2)-Struktur fuer {w}")


def macro_class(w: str) -> str:
    prof = profile_of_word(w)

    if prof == (4,):
        return "U"
    if prof == (3, 1):
        return "D"
    if prof == (2, 2):
        return classify_22_word(w)
    if prof == (2, 1, 1):
        return "M"
    if prof == (1, 1, 1, 1):
        if w in CYCLIC_PLUS:
            return "Z+"
        if w in CYCLIC_MINUS:
            return "Z-"
        return "Z0"

    raise RuntimeError(f"Unbekanntes Profil fuer {w}")


def morley_M(w: str) -> float:
    cls = macro_class(w)
    return {
        "N": 0.00,
        "B": 0.30,
        "A": 0.70,
        "M": 0.55,
        "D": 0.25,
        "U": 0.10,
        "Z+": 1.00,
        "Z-": 1.00,
        "Z0": 0.85,
    }[cls]


def chirality_chi(w: str) -> int:
    cls = macro_class(w)
    if cls == "Z+":
        return +1
    if cls == "Z-":
        return -1
    if cls == "A":
        return +1 if w < w[::-1] else -1
    return 0


def walter_T(w: str) -> float:
    cls = macro_class(w)
    return {
        "N": 0.80,
        "B": 0.50,
        "A": 0.35,
        "M": 0.45,
        "D": 0.20,
        "U": 0.15,
        "Z+": 0.20,
        "Z-": 0.20,
        "Z0": 0.20,
    }[cls]


@dataclass
class QModelParams:
    lambda_M: float = 0.5
    lambda_chi: float = 0.3
    lambda_T: float = 0.4
    q_threshold_neutral: float = 0.18


def activity_Q(w: str, params: QModelParams = QModelParams()) -> float:
    q = (
        params.lambda_M * morley_M(w)
        + params.lambda_chi * abs(chirality_chi(w))
        - params.lambda_T * walter_T(w)
    )
    return clip(q, 0.0, 1.0)


def sigma_projection(w: str, params: QModelParams = QModelParams()) -> str:
    q = activity_Q(w, params)
    cls = macro_class(w)
    cnt = Counter(w)

    if q < params.q_threshold_neutral:
        return "0"

    if cls == "Z+":
        mapping = {
            "EABC": "C",
            "ABCE": "E",
            "BCEA": "A",
            "CEAB": "B",
        }
        return mapping.get(w, "mix")

    if cls == "Z-":
        mapping = {
            "EACB": "B",
            "ACBE": "E",
            "CBEA": "A",
            "BEAC": "C",
        }
        return mapping.get(w, "mix")

    if cls == "A":
        syms = sorted(cnt.keys())
        return syms[0] if syms[0] in "EABC" else "mix"

    if cls == "B":
        syms = sorted(cnt.keys())
        return "0" if activity_Q(w, params) < 0.25 else syms[0]

    if cls == "N":
        return "0"

    if cls == "M":
        return [k for k, v in cnt.items() if v == 2][0]

    if cls == "D":
        return [k for k, v in cnt.items() if v == 3][0]

    if cls == "U":
        return "0" if w[0] == "E" else w[0]

    if cls == "Z0":
        return w[-1]

    return "0"


@dataclass
class EnergyParams:
    theta_p: float = 3.0
    kappa_p: float = 2.0
    theta_d: float = 5.0
    kappa_d: float = 1.3


def alpha_pinski(sqrt_s: float, p: EnergyParams) -> float:
    return sigmoid(p.kappa_p * (sqrt_s - p.theta_p))


def gamma_detector(sqrt_s: float, p: EnergyParams) -> float:
    return sigmoid(p.kappa_d * (sqrt_s - p.theta_d))


def beta_neutral(sqrt_s: float, p: EnergyParams) -> float:
    return 1.0 - alpha_pinski(sqrt_s, p)


def word_to_word_weight(w: str, w2: str, sqrt_s: float, qparams: QModelParams, eparams: EnergyParams) -> float:
    """Gewicht fuer interne Uebergaenge W -> W'."""
    a = alpha_pinski(sqrt_s, eparams)
    b = beta_neutral(sqrt_s, eparams)

    m1, m2 = morley_M(w), morley_M(w2)
    t1, t2 = walter_T(w), walter_T(w2)
    q1, q2 = activity_Q(w, qparams), activity_Q(w2, qparams)
    chi1, chi2 = chirality_chi(w), chirality_chi(w2)
    c1, c2 = macro_class(w), macro_class(w2)

    weight = 0.0

    if c2 in ("N", "B"):
        weight += b * (0.20 + 0.25 * t1 - 0.10 * m1)
    if c2 in ("A", "Z+", "Z-", "Z0", "M"):
        weight += a * (0.15 + 0.20 * m1 + 0.15 * q1 - 0.10 * t1)
    if c2 in ("D", "U"):
        weight += 0.05 + 0.05 * (1.0 - t1)
    if chi1 != 0 and chi2 != 0 and chi1 == chi2:
        weight += 0.08
    if c1 == c2:
        weight += 0.05

    s1 = sigma_projection(w, qparams)
    s2 = sigma_projection(w2, qparams)
    if s1 == s2 and s1 != "0":
        weight += 0.08

    weight += 0.05 * (t1 + t2)
    _ = m2, q2
    return max(weight, 0.0)


def word_to_absorbing_weight(
    w: str, absorb_state: str, sqrt_s: float, qparams: QModelParams, eparams: EnergyParams
) -> float:
    """Gewicht fuer W -> A/B/C/E/0."""
    g = gamma_detector(sqrt_s, eparams)
    a = alpha_pinski(sqrt_s, eparams)

    m = morley_M(w)
    t = walter_T(w)
    q = activity_Q(w, qparams)
    sig = sigma_projection(w, qparams)

    if absorb_state == "0":
        return max((1.0 - g) * (0.10 + 0.35 * t - 0.15 * m - 0.10 * q) + 0.05 * (sig == "0"), 0.0)

    base = g * (0.03 + 0.30 * q + 0.20 * m - 0.12 * t) + a * (0.02 + 0.08 * m)

    if sig == absorb_state:
        base += 0.25
    elif sig != "0":
        base += 0.04

    return max(base, 0.0)


def build_full_transition_matrix(
    sqrt_s: float, qparams: QModelParams = QModelParams(), eparams: EnergyParams = EnergyParams()
) -> np.ndarray:
    P = np.zeros((N_TOTAL, N_TOTAL), dtype=float)

    for w in WORDS:
        i = STATE_INDEX[w]
        row = np.zeros(N_TOTAL, dtype=float)

        for w2 in all_local_neighbors(w):
            j = STATE_INDEX[w2]
            row[j] += word_to_word_weight(w, w2, sqrt_s, qparams, eparams)

        row[i] += 0.03 + 0.04 * walter_T(w)

        for abs_state in ABSORBING:
            j = STATE_INDEX[abs_state]
            row[j] += word_to_absorbing_weight(w, abs_state, sqrt_s, qparams, eparams)

        s = row.sum()
        if s <= 0:
            row[i] = 1.0
        else:
            row = row / s

        P[i] = row

    for astate in ABSORBING:
        i = STATE_INDEX[astate]
        P[i, i] = 1.0

    return P


def transient_indices():
    return [STATE_INDEX[w] for w in WORDS]


def absorbing_indices():
    return [STATE_INDEX[s] for s in ABSORBING]


def fundamental_matrix(P: np.ndarray):
    t_idx = transient_indices()
    a_idx = absorbing_indices()
    Q = P[np.ix_(t_idx, t_idx)]
    R = P[np.ix_(t_idx, a_idx)]
    N = np.linalg.inv(np.eye(Q.shape[0]) - Q)
    B = N @ R
    return Q, R, N, B


def exact_absorption_probabilities(P: np.ndarray, start_word: str):
    _, _, _, B = fundamental_matrix(P)
    i = WORD_INDEX[start_word]
    probs = B[i]
    return {ABSORBING[k]: probs[k] for k in range(len(ABSORBING))}


def exact_mean_absorption_steps(P: np.ndarray, start_word: str):
    _, _, N, _ = fundamental_matrix(P)
    t = N @ np.ones((N.shape[0], 1))
    i = WORD_INDEX[start_word]
    return float(t[i, 0])


def quasi_stationary_distribution(P: np.ndarray):
    Q = transient_submatrix(P)
    eigvals, eigvecs = np.linalg.eig(Q.T)
    idx = np.argmax(np.real(eigvals))
    vec = np.real(eigvecs[:, idx])
    vec = np.abs(vec)
    vec = vec / vec.sum()
    return {WORDS[i]: vec[i] for i in range(len(WORDS))}


def transient_submatrix(P):
    """Gibt die transiente Teilmatrix Q der absorbierenden Kette zurueck."""
    Q, _, _, _ = fundamental_matrix(P)
    return Q


# ============================================================
# Zusätzliche Wortinvarianten und verfeinerte Dynamik
# ============================================================

def hamming_distance(w1: str, w2: str) -> int:
    return sum(a != b for a, b in zip(w1, w2))


def inversion_count(w: str) -> int:
    vals = [LETTER_ORDER[ch] for ch in w]
    inv = 0
    for i in range(len(vals)):
        for j in range(i + 1, len(vals)):
            if vals[i] > vals[j]:
                inv += 1
    return inv


def palindromic_score(w: str) -> int:
    score = 0
    if w[0] == w[3]:
        score += 1
    if w[1] == w[2]:
        score += 1
    return score


def cyclic_reference_distance(w: str, ref: str) -> int:
    rots = [ref[i:] + ref[:i] for i in range(4)]
    return min(hamming_distance(w, r) for r in rots)


def cyclic_plus_score(w: str) -> float:
    d = cyclic_reference_distance(w, "EABC")
    return 1.0 - d / 4.0


def cyclic_minus_score(w: str) -> float:
    d = cyclic_reference_distance(w, "EACB")
    return 1.0 - d / 4.0


def tiny_hash_bias(w: str) -> float:
    """Winziger deterministischer Bias gegen Spektralentartung."""
    s = sum((i + 1) * (ord(ch) - 64) for i, ch in enumerate(w))
    return 1e-4 * (s % 17)


def mixed_subclass(w: str) -> str | None:
    if profile_of_word(w) != (2, 1, 1):
        return None

    cnt = Counter(w)
    double_symbol = [k for k, v in cnt.items() if v == 2][0]
    pos = [i for i, ch in enumerate(w) if ch == double_symbol]

    if pos[1] == pos[0] + 1:
        return "M1"
    if pos in ([0, 2], [1, 3]):
        return "M3"
    return "M2"


def refined_macro_class(w: str) -> str:
    cls = macro_class(w)
    if cls != "M":
        return cls
    return mixed_subclass(w)


def class_is_core_like(cls: str) -> bool:
    return cls in ("N", "B", "M1")


def class_is_active_like(cls: str) -> bool:
    return cls in ("A", "M3", "Z+", "Z-", "Z0")


def class_is_mixed_internal(cls: str) -> bool:
    return cls in ("M1", "M2", "M3")


def refined_morley_M(w: str) -> float:
    cls = refined_macro_class(w)

    if cls == "N":
        return 0.00
    if cls == "B":
        return 0.28
    if cls == "A":
        return 0.72
    if cls == "M1":
        return 0.42
    if cls == "M2":
        return 0.56
    if cls == "M3":
        return 0.68
    if cls == "D":
        return 0.22
    if cls == "U":
        return 0.08
    if cls == "Z+":
        return 1.00
    if cls == "Z-":
        return 1.00
    if cls == "Z0":
        return 0.65 + 0.20 * max(cyclic_plus_score(w), cyclic_minus_score(w))
    return 0.3


def refined_walter_T(w: str) -> float:
    cls = refined_macro_class(w)

    if cls == "N":
        return 0.82
    if cls == "B":
        return 0.56
    if cls == "A":
        return 0.30
    if cls == "M1":
        return 0.58
    if cls == "M2":
        return 0.44
    if cls == "M3":
        return 0.32
    if cls == "D":
        return 0.18
    if cls == "U":
        return 0.12
    if cls in ("Z+", "Z-"):
        return 0.18
    if cls == "Z0":
        return 0.22
    return 0.35


def refined_activity_Q(
    w: str, lambda_M=0.52, lambda_chi=0.28, lambda_T=0.38, lambda_pal=0.10
) -> float:
    M = refined_morley_M(w)
    chi = abs(chirality_chi(w))
    T = refined_walter_T(w)
    pal = palindromic_score(w) / 2.0
    q = lambda_M * M + lambda_chi * chi - lambda_T * T - lambda_pal * pal
    return max(0.0, min(1.0, q))


@dataclass
class RefinedIIIParams:
    # interne Dynamik stärken
    internal_scale: float = 2.2

    # Absorption schwächen
    absorption_scale: float = 0.58

    # Selbstübergänge erhöhen
    self_loop_base: float = 0.045
    self_loop_walter: float = 0.11

    # spezielle Kopplungen der M-Unterklassen
    m1_to_core_boost: float = 0.22
    m2_internal_boost: float = 0.18
    m3_to_active_boost: float = 0.24

    # arithmetische interne Verstärkung
    n_internal_bias: float = 0.18


@dataclass
class HurwitzPressureParams:
    H_inf: float = 20.0
    z_boost: float = 0.18
    core_boost: float = 0.14
    absorb_bias: float = 0.08


def hurwitz_pressure_factor_for_class(cls: str, hparams: HurwitzPressureParams) -> float:
    H = hparams.H_inf

    if cls in ("N", "B", "M1"):
        return 1.0 + hparams.core_boost * H

    if cls in ("M2",):
        return 1.0 + 0.5 * hparams.core_boost * H

    if cls in ("A", "M3", "Z+", "Z-", "Z0"):
        return 1.0 + hparams.z_boost * H

    return 1.0 + 0.04 * H


def structured_neighbors(w: str):
    nbrs = set(all_local_neighbors(w))

    if len(w) == 4:
        nbrs.add(w[0] + w[1] + w[1] + w[0])
        nbrs.add(w[0] + w[2] + w[2] + w[0])

    nbrs.add(w[1:] + w[0])
    nbrs.add(w[-1] + w[:-1])

    cnt = Counter(w)
    doubles = [k for k, v in cnt.items() if v >= 2]
    if doubles:
        d = doubles[0]
        rest = "".join(ch for ch in w if ch != d)
        candidate = d + d + rest[:2]
        if len(candidate) == 4:
            nbrs.add(candidate)

    nbrs.discard(w)
    nbrs = {u for u in nbrs if len(u) == 4 and all(ch in ALPHABET for ch in u)}
    return sorted(nbrs)


def refined_word_to_word_weight(
    w: str,
    w2: str,
    sqrt_s: float,
    qparams: QModelParams,
    eparams: EnergyParams,
    rparams: RefinedIIIParams = RefinedIIIParams(),
) -> float:
    a = alpha_pinski(sqrt_s, eparams)
    b = beta_neutral(sqrt_s, eparams)
    _ = gamma_detector(sqrt_s, eparams)

    M1, M2 = refined_morley_M(w), refined_morley_M(w2)
    T1, T2 = refined_walter_T(w), refined_walter_T(w2)
    Q1, Q2 = refined_activity_Q(w), refined_activity_Q(w2)
    chi1, chi2 = chirality_chi(w), chirality_chi(w2)
    cls1 = refined_macro_class(w)
    cls2 = refined_macro_class(w2)

    dH = hamming_distance(w, w2)
    dI = abs(inversion_count(w) - inversion_count(w2))
    pal2 = palindromic_score(w2)

    sig1 = sigma_projection(w, qparams)
    sig2 = sigma_projection(w2, qparams)

    locality = np.exp(-0.55 * dH - 0.20 * dI)
    store = b * (0.20 + 0.30 * T2 + 0.10 * pal2 - 0.10 * M2)
    active = a * (0.18 + 0.28 * M2 + 0.18 * Q2 - 0.15 * T2)
    chir = 0.05 if (chi1 != 0 and chi1 == chi2) else 0.0
    sig_bonus = 0.06 if (sig1 == sig2 and sig1 != "0") else 0.0
    pal_bonus = 0.08 * (sig2 == "0") * pal2

    subclass_bonus = 0.0
    if cls1 == "M1" and cls2 in ("N", "B"):
        subclass_bonus += rparams.m1_to_core_boost
    if cls1 == "M2" and cls2 in ("M2", "M3"):
        subclass_bonus += rparams.m2_internal_boost
    if cls1 == "M3" and cls2 in ("A", "Z+", "Z-", "Z0"):
        subclass_bonus += rparams.m3_to_active_boost
    if cls1 in ("N", "B", "M1", "M2", "M3") and cls2 not in ABSORBING:
        subclass_bonus += rparams.n_internal_bias

    weight = rparams.internal_scale * locality * (store + active + chir + sig_bonus + pal_bonus + subclass_bonus)
    weight += tiny_hash_bias(w2)
    _ = M1, T1, Q1
    return max(weight, 0.0)


def refined_word_to_word_weight_v3(
    w: str,
    w2: str,
    sqrt_s: float,
    qparams: QModelParams,
    eparams: EnergyParams,
    rparams: RefinedIIIParams,
    hparams: HurwitzPressureParams,
    N: int | None = None,
) -> float:
    a = alpha_pinski(sqrt_s, eparams)
    b = beta_neutral(sqrt_s, eparams)

    c1 = refined_macro_class(w)
    c2 = refined_macro_class(w2)

    M1, M2 = refined_morley_M(w), refined_morley_M(w2)
    T1, T2 = refined_walter_T(w), refined_walter_T(w2)
    Q1, Q2 = refined_activity_Q(w), refined_activity_Q(w2)
    chi1, chi2 = chirality_chi(w), chirality_chi(w2)

    dH = hamming_distance(w, w2)
    dI = abs(inversion_count(w) - inversion_count(w2))
    pal2 = palindromic_score(w2)

    sig1 = sigma_projection(w, qparams)
    sig2 = sigma_projection(w2, qparams)

    locality = np.exp(-0.45 * dH - 0.14 * dI)
    store = b * (0.18 + 0.28 * T2 + 0.12 * pal2 - 0.08 * M2)
    active = a * (0.14 + 0.24 * M2 + 0.16 * Q2 - 0.10 * T2)
    base = rparams.internal_scale * locality * (store + active)

    if chi1 != 0 and chi1 == chi2:
        base += 0.05

    if sig1 == sig2 and sig1 != "0":
        base += 0.05

    if c1 == "M1" and c2 in ("N", "B", "M1"):
        base += rparams.m1_to_core_boost

    if c1 == "M2" and c2 in ("M1", "M2", "M3"):
        base += rparams.m2_internal_boost

    if c1 == "M3" and c2 in ("A", "Z+", "Z-", "Z0", "M3"):
        base += rparams.m3_to_active_boost

    if class_is_core_like(c1) and class_is_core_like(c2):
        base += 0.08

    if class_is_active_like(c1) and class_is_active_like(c2):
        base += 0.08

    if sig2 == "0":
        base += 0.05 * pal2 + 0.05 * T2

    if N is not None:
        if N % 3 == 0 and c2 in ("B", "M2"):
            base += rparams.n_internal_bias

        if N % 5 == 0 and c2 in ("A", "M3"):
            base += 0.16

        if N % 7 == 0 and c2 in ("Z0", "Z-", "M3"):
            base += 0.16

        if N % 2 == 0 and c2 in ("N", "B", "M1"):
            base += 0.12

    hp1 = hurwitz_pressure_factor_for_class(c1, hparams)
    hp2 = hurwitz_pressure_factor_for_class(c2, hparams)
    base *= (0.35 * hp1 + 0.65 * hp2)

    base += tiny_hash_bias(w2)
    _ = M1, T1, Q1
    return max(base, 0.0)


def refined_word_to_absorbing_weight(
    w: str,
    absorb_state: str,
    sqrt_s: float,
    qparams: QModelParams,
    eparams: EnergyParams,
    rparams: RefinedIIIParams = RefinedIIIParams(),
) -> float:
    g = gamma_detector(sqrt_s, eparams)
    a = alpha_pinski(sqrt_s, eparams)

    M = refined_morley_M(w)
    T = refined_walter_T(w)
    Q = refined_activity_Q(w)
    sig = sigma_projection(w, qparams)
    pal = palindromic_score(w) / 2.0

    if absorb_state == "0":
        val = (1.0 - g) * (0.08 + 0.42 * T + 0.18 * pal - 0.18 * M - 0.12 * Q)
        if sig == "0":
            val += 0.08
        return max(val, 0.0)

    val = g * (0.02 + 0.34 * Q + 0.24 * M - 0.14 * T) + a * (0.03 + 0.10 * M)

    if sig == absorb_state:
        val += 0.22
    elif sig != "0":
        val += 0.03

    val = rparams.absorption_scale * val + tiny_hash_bias(w)
    return max(val, 0.0)


def refined_word_to_absorbing_weight_v3(
    w: str,
    absorb_state: str,
    sqrt_s: float,
    qparams: QModelParams,
    eparams: EnergyParams,
    rparams: RefinedIIIParams,
    hparams: HurwitzPressureParams,
    N: int | None = None,
) -> float:
    g = gamma_detector(sqrt_s, eparams)
    a = alpha_pinski(sqrt_s, eparams)

    M = refined_morley_M(w)
    T = refined_walter_T(w)
    Q = refined_activity_Q(w)
    sig = sigma_projection(w, qparams)
    pal = palindromic_score(w) / 2.0
    cls = refined_macro_class(w)

    if absorb_state == "0":
        val = (1.0 - g) * (0.08 + 0.46 * T + 0.16 * pal - 0.14 * M - 0.08 * Q)
        if sig == "0":
            val += 0.08
        if cls in ("N", "B", "M1"):
            val += 0.05
        H = hparams.H_inf
        val *= (1.0 + 0.06 * H * (cls in ("N", "B", "M1")))
        return max(rparams.absorption_scale * val, 0.0)

    val = g * (0.01 + 0.22 * Q + 0.16 * M - 0.12 * T) + a * (0.01 + 0.05 * M)

    if sig == absorb_state:
        val += 0.12
    elif sig != "0":
        val += 0.015

    if N is not None:
        if N % 5 == 0 and absorb_state == "A":
            val += 0.04
        if N % 7 == 0 and absorb_state == "C":
            val += 0.04
        if N % 3 == 0 and absorb_state == "B":
            val += 0.03

    H = hparams.H_inf
    val *= (1.0 + hparams.absorb_bias * H * (cls in ("A", "M3", "Z+", "Z-", "Z0")))

    val += tiny_hash_bias(w)
    return max(rparams.absorption_scale * val, 0.0)


def build_full_transition_matrix_refined(
    sqrt_s: float,
    qparams: QModelParams = QModelParams(),
    eparams: EnergyParams = EnergyParams(),
    rparams: RefinedIIIParams = RefinedIIIParams(),
) -> np.ndarray:
    P = np.zeros((N_TOTAL, N_TOTAL), dtype=float)

    for w in WORDS:
        i = STATE_INDEX[w]
        row = np.zeros(N_TOTAL, dtype=float)

        for w2 in structured_neighbors(w):
            j = STATE_INDEX[w2]
            r = refined_word_to_word_weight(w, w2, sqrt_s, qparams, eparams, rparams)
            row[j] += r

        row[i] += rparams.self_loop_base + rparams.self_loop_walter * refined_walter_T(w) + tiny_hash_bias(w)

        for astate in ABSORBING:
            j = STATE_INDEX[astate]
            row[j] += refined_word_to_absorbing_weight(w, astate, sqrt_s, qparams, eparams, rparams)

        s = row.sum()
        if s <= 0:
            row[i] = 1.0
        else:
            row /= s
        P[i] = row

    for astate in ABSORBING:
        i = STATE_INDEX[astate]
        P[i, i] = 1.0

    return P


def build_full_transition_matrix_refined_v3(sqrt_s: float,
                                            qparams: QModelParams = QModelParams(),
                                            eparams: EnergyParams = EnergyParams(),
                                            rparams: RefinedIIIParams = RefinedIIIParams(),
                                            hparams: HurwitzPressureParams = HurwitzPressureParams(),
                                            N: int | None = None) -> np.ndarray:
    P = np.zeros((N_TOTAL, N_TOTAL), dtype=float)

    for w in WORDS:
        i = STATE_INDEX[w]
        row = np.zeros(N_TOTAL, dtype=float)

        nbrs = structured_neighbors(w)

        for w2 in nbrs:
            j = STATE_INDEX[w2]
            row[j] += refined_word_to_word_weight_v3(
                w, w2, sqrt_s, qparams, eparams, rparams, hparams, N=N
            )

        row[i] += (
            rparams.self_loop_base
            + rparams.self_loop_walter * refined_walter_T(w)
            + tiny_hash_bias(w)
        )

        for astate in ABSORBING:
            j = STATE_INDEX[astate]
            row[j] += refined_word_to_absorbing_weight_v3(
                w, astate, sqrt_s, qparams, eparams, rparams, hparams, N=N
            )

        s = row.sum()
        if s <= 0:
            row[i] = 1.0
        else:
            row /= s
        P[i] = row

    for astate in ABSORBING:
        i = STATE_INDEX[astate]
        P[i, i] = 1.0

    return P

# ============================================================
# N-abhängige Symbolgewichte
# ============================================================

def symbol_weights_from_N(N: int):
    """
    Erweitertes arithmetisches Encoding N -> Gewichte für E,A,B,C
    mit zusätzlicher Sensitivität auf 11 und 13.
    """
    weights = {}

    weights["E"] = 1.0
    weights["A"] = 1.0
    weights["B"] = 1.0
    weights["C"] = 1.0

    # grobe Teilbarkeit
    if N % 2 == 0:
        weights["E"] += 0.8
    if N % 3 == 0:
        weights["E"] += 0.5
        weights["B"] += 0.4
    if N % 5 == 0:
        weights["A"] += 0.5
    if N % 7 == 0:
        weights["C"] += 0.5

    # Version V: 11/13 etwas schwächer, damit sie nicht alles überdominieren
    if N % 11 == 0:
        weights["E"] += 0.14
        weights["B"] += 0.11
        weights["C"] += 0.08
    if N % 13 == 0:
        weights["E"] += 0.14
        weights["A"] += 0.11
        weights["C"] += 0.08

    # Restklassen
    if N % 4 == 1:
        weights["A"] += 0.3
    if N % 3 == 1:
        weights["B"] += 0.3
    if N % 12 in (1, 5, 7, 11):
        weights["C"] += 0.3

    # Erweiterte Restklassenstruktur
    if N % 11 in (1, 10):
        weights["B"] += 0.06
        weights["E"] += 0.05
    if N % 13 in (1, 12):
        weights["A"] += 0.06
        weights["E"] += 0.05

    return weights


# ============================================================
# N-abhängiger Wortbonus
# ============================================================

def word_bonus_from_N(w: str, N: int):
    """
    Zusätzlicher Strukturbonus für Wörter abhängig von N.
    Noch heuristisch.
    """
    bonus = 1.0
    cls = refined_macro_class(w)
    sig = sigma_projection(w)

    # gerade / durch 3 teilbare Zahlen stärker E-/Neutralitätsnah
    if N % 2 == 0 and cls in ("N", "B", "M1", "M2", "M3"):
        bonus += 0.2
    if N % 3 == 0 and "B" in w:
        bonus += 0.15

    # mod-12 Struktur
    if N % 12 in (1, 5, 7, 11) and cls in ("Z+", "Z-", "A"):
        bonus += 0.2

    # falls Projektion direkt zu kleinen Restklassensignalen passt
    if N % 5 == 0 and sig == "A":
        bonus += 0.15
    if N % 7 == 0 and sig == "C":
        bonus += 0.15

    # neue 11-/13-Sensitivität (Version V: reduziert)
    if N % 11 == 0:
        if cls in ("B", "M1", "M2", "M3", "Z0"):
            bonus += 0.09
        if sig in ("B", "E"):
            bonus += 0.05

    if N % 13 == 0:
        if cls in ("A", "M1", "M2", "M3", "Z0"):
            bonus += 0.09
        if sig in ("A", "E"):
            bonus += 0.05

    # Version V: mixed explizit stärken für Mehrfachfaktoren aus kleinen Achsen
    small_hits = sum([
        int(N % 3 == 0),
        int(N % 5 == 0),
        int(N % 7 == 0),
        int(N % 11 == 0),
        int(N % 13 == 0),
    ])

    if small_hits >= 2 and cls in ("M1", "M2", "M3"):
        bonus += 0.18
    if small_hits >= 3 and cls in ("M1", "M2", "M3"):
        bonus += 0.12

    # Version V: family-7 leicht stabilisieren
    if N % 7 == 0 and cls in ("Z+", "Z-", "Z0", "A"):
        bonus += 0.08

    # Version V.1.1: lokaler mixed-Boost nur für 9x7-Kopplung
    # Ziel: 63 eher in mixed statt rein in family-3, aber 21 nicht mitziehen
    if (N % 9 == 0) and (N % 7 == 0):
        if cls == "M":
            bonus += 0.10
        if cls in ("A", "B"):
            bonus += 0.04

    return bonus


# ============================================================
# Startverteilung pi_N über Wörter
# ============================================================

def initial_distribution_from_N(N: int):
    """
    Gibt eine normierte Startverteilung pi_N auf den 256 Wörtern zurück.
    """
    sw = symbol_weights_from_N(N)
    raw = np.zeros(len(WORDS), dtype=float)

    for i, w in enumerate(WORDS):
        prod_weight = 1.0
        for ch in w:
            prod_weight *= sw[ch]
        prod_weight *= word_bonus_from_N(w, N)
        raw[i] = prod_weight

    raw /= raw.sum()
    return raw


# ============================================================
# N-abhängige Modifikation des Operators
# ============================================================

@dataclass
class NOperatorParams:
    spectral_bias_strength: float = 0.25
    absorption_bias_strength: float = 0.35
    use_refined_base: bool = True


def build_N_dependent_transition_matrix_v3(N: int,
                                           sqrt_s: float,
                                           qparams=QModelParams(),
                                           eparams=EnergyParams(),
                                           rparams=RefinedIIIParams(),
                                           hparams=HurwitzPressureParams()):
    return build_full_transition_matrix_refined_v3(
        sqrt_s=sqrt_s,
        qparams=qparams,
        eparams=eparams,
        rparams=rparams,
        hparams=hparams,
        N=N
    )


def build_N_dependent_transition_matrix(N: int,
                                        sqrt_s: float,
                                        qparams=QModelParams(),
                                        eparams=EnergyParams(),
                                        nparams=NOperatorParams(),
                                        hparams=HurwitzPressureParams(),
                                        rparams=RefinedIIIParams()):
    """
    Baut P_N aus der Basis-Matrix, modifiziert durch N-abhängige Faktoren.
    """
    base_builder = build_full_transition_matrix_refined if nparams.use_refined_base else build_full_transition_matrix
    if nparams.use_refined_base:
        P = base_builder(sqrt_s, qparams, eparams, rparams).copy()
    else:
        P = base_builder(sqrt_s, qparams, eparams).copy()

    # Wortzustände anpassen
    for w in WORDS:
        i = STATE_INDEX[w]
        row = P[i].copy()

        cls = refined_macro_class(w)
        sig = sigma_projection(w)
        hp = hurwitz_pressure_factor_for_class(cls, hparams)

        # Faktor für interne Speicherung vs. Projektion
        internal_factor = 1.0
        absorb_factor = 1.0

        # arithmetische Heuristik
        if N % 2 == 0 and cls in ("N", "B", "M1", "M2", "M3"):
            internal_factor += nparams.spectral_bias_strength

        if N % 3 == 0 and "B" in w:
            internal_factor += 0.15

        if N % 5 == 0 and sig == "A":
            absorb_factor += nparams.absorption_bias_strength

        if N % 7 == 0 and sig == "C":
            absorb_factor += nparams.absorption_bias_strength

        if N % 12 in (1, 5, 7, 11) and cls in ("Z+", "Z-", "A"):
            internal_factor += 0.10

        # neue interne Achsen für 11 und 13 (reduziert)
        if N % 11 == 0 and cls in ("B", "M2", "Z0"):
            internal_factor += 0.08
        if N % 13 == 0 and cls in ("A", "M3", "Z0"):
            internal_factor += 0.08

        # Version V: family-7 stabilisieren
        if N % 7 == 0 and cls in ("Z+", "Z-", "Z0", "A"):
            internal_factor += 0.10

        # Version V: mixed-Bereich für Mehrfachhits explizit verstärken
        small_hits = sum([
            int(N % 3 == 0),
            int(N % 5 == 0),
            int(N % 7 == 0),
            int(N % 11 == 0),
            int(N % 13 == 0),
        ])

        if small_hits >= 2 and cls in ("M1", "M2", "M3", "A", "B"):
            internal_factor += 0.10
        if small_hits >= 3 and cls in ("M1", "M2", "M3", "A", "B", "Z0"):
            internal_factor += 0.08

        # Version V.1.1: zusätzlicher 9x7-Mixed-Drift
        if (N % 9 == 0) and (N % 7 == 0):
            if cls in ("M1", "M2", "M3", "A", "B"):
                internal_factor += 0.08

        # Hurwitzdruck als Sektorverstärker
        internal_factor *= hp

        # interne Wortzustände
        for w2 in WORDS:
            j = STATE_INDEX[w2]
            row[j] *= internal_factor

        # absorbierende Kanäle
        for a in ABSORBING:
            j = STATE_INDEX[a]
            if a == sig:
                row[j] *= absorb_factor * (1.0 + hparams.absorb_bias * hparams.H_inf)
            elif a == "0" and cls == "N":
                row[j] *= 1.0 + 0.15 + 0.02 * hparams.H_inf
            else:
                row[j] *= 1.0

        # Version V: mixed-Zahlen weniger stark in reine Familien absorbieren lassen
        if small_hits >= 2:
            for a in ("A", "B", "C", "E"):
                j = STATE_INDEX[a]
                row[j] *= 0.92
            j0 = STATE_INDEX["0"]
            row[j0] *= 1.03

        # Version V.1.1: bei 9x7-Kopplung reine Absorption weiter leicht dämpfen
        if (N % 9 == 0) and (N % 7 == 0):
            for a in ("A", "B", "C", "E"):
                j = STATE_INDEX[a]
                row[j] *= 0.97
            j0 = STATE_INDEX["0"]
            row[j0] *= 1.01

        # normieren
        s = row.sum()
        row = row / s
        P[i] = row

    return P


# ============================================================
# Gemischte Startanalyse mit pi_N
# ============================================================

def absorption_profile_from_initial_distribution(P: np.ndarray, pi_words: np.ndarray):
    """
    Mischt die exakten Absorptionswahrscheinlichkeiten über eine Startverteilung pi_words.
    """
    profile = {a: 0.0 for a in ABSORBING}
    mean_steps = 0.0

    for i, w in enumerate(WORDS):
        p = pi_words[i]
        probs = exact_absorption_probabilities(P, w)
        steps = exact_mean_absorption_steps(P, w)

        for a in ABSORBING:
            profile[a] += p * probs[a]
        mean_steps += p * steps

    profile["mean_steps"] = mean_steps
    return profile


# ============================================================
# Spektralprofil eines N
# ============================================================

def leading_eigenvalues_of_transient(P: np.ndarray, k: int = 12):
    Q = transient_submatrix(P)
    vals = np.linalg.eigvals(Q)
    vals = np.array(sorted(vals, key=lambda z: -abs(z)))
    return vals[:k]

def qsd_macro_profile(P: np.ndarray):
    qsd = quasi_stationary_distribution(P)
    agg = defaultdict(float)
    for w, p in qsd.items():
        agg[refined_macro_class(w)] += p
    return dict(agg)

def N_spectral_profile(N: int,
                       sqrt_s: float,
                       qparams=QModelParams(),
                       eparams=EnergyParams(),
                       rparams=RefinedIIIParams(),
                       nparams=NOperatorParams(),
                       hparams=HurwitzPressureParams(),
                       k_eigs: int = 8):
    piN = initial_distribution_from_N(N)

    PN = build_N_dependent_transition_matrix(
        N, sqrt_s,
        qparams=qparams,
        eparams=eparams,
        nparams=nparams,
        hparams=hparams,
        rparams=rparams
    )

    absorption = absorption_profile_from_initial_distribution(PN, piN)
    eigs = leading_eigenvalues_of_transient(PN, k=k_eigs)
    qsd_macro = qsd_macro_profile(PN)

    e_feats = stabilized_e_quadruple_features(N)
    s2_feats = stabilized_s2_features(N)

    return {
        "N": N,
        "sqrt_s": sqrt_s,
        "initial_distribution": piN,
        "absorption": absorption,
        "eigenvalues": eigs,
        "qsd_macro": qsd_macro,
        "e_features": e_feats,
        "s2_features": s2_feats,
    }


# ============================================================
# Vergleich mehrerer Zahlen
# ============================================================

def compare_numbers(N_list, sqrt_s: float,
                    qparams=QModelParams(),
                    eparams=EnergyParams(),
                    nparams=NOperatorParams(),
                    rparams=RefinedIIIParams(),
                    hparams=HurwitzPressureParams()):
    profiles = []
    for N in N_list:
        profiles.append(
            N_spectral_profile(
                N=N,
                sqrt_s=sqrt_s,
                qparams=qparams,
                eparams=eparams,
                rparams=rparams,
                nparams=nparams,
                hparams=hparams,
                k_eigs=8
            )
        )
    return profiles


def print_number_profiles(profiles):
    for prof in profiles:
        print("=" * 80)
        print(f"N = {prof['N']}, sqrt(s) = {prof['sqrt_s']}")
        print("Absorptionsprofil:")
        for k, v in prof["absorption"].items():
            if k == "mean_steps":
                print(f"  {k}: {v:.4f}")
            else:
                print(f"  {k}: {v:.4f}")

        print("Führende Eigenwerte:")
        for z in prof["eigenvalues"]:
            print(f"  {z:.6f}")

        print("QSD-Makroprofil:")
        for cls, v in sorted(prof["qsd_macro"].items()):
            print(f"  {cls}: {v:.4f}")


# ============================================================
# Distanz zwischen Zahlenprofilen
# ============================================================

def profile_vector(profile):
    """
    Baut einen numerischen Vektor aus Absorption + QSD-Makro + führenden Eigenwertbeträgen.
    """
    abs_part = [
        profile["absorption"]["A"],
        profile["absorption"]["B"],
        profile["absorption"]["C"],
        profile["absorption"]["E"],
        profile["absorption"]["0"],
        profile["absorption"]["mean_steps"],
    ]

    macros = ["N", "A", "B", "M1", "M2", "M3", "D", "U", "Z+", "Z-", "Z0"]
    qsd_part = [profile["qsd_macro"].get(m, 0.0) for m in macros]

    eig_part = [abs(z) for z in profile["eigenvalues"]]

    return np.array(abs_part + qsd_part + eig_part, dtype=float)

def pairwise_profile_distances(profiles):
    """
    Einfache euklidische Distanzen der Profilvektoren.
    """
    vecs = {p["N"]: profile_vector(p) for p in profiles}
    Ns = list(vecs.keys())
    D = np.zeros((len(Ns), len(Ns)))

    for i, Ni in enumerate(Ns):
        for j, Nj in enumerate(Ns):
            D[i, j] = np.linalg.norm(vecs[Ni] - vecs[Nj])

    return Ns, D


# ============================================================
# E-Quattropel-Features
# ============================================================

def divisors_of_n(n: int):
    ds = []
    r = int(math.isqrt(n))
    for d in range(1, r + 1):
        if n % d == 0:
            ds.append(d)
            if d * d != n:
                ds.append(n // d)
    return sorted(ds)


def residue_family_12(n: int):
    r = n % 12
    if r == 1:
        return "e"
    if r == 5:
        return "a"
    if r == 7:
        return "b"
    if r == 11:
        return "c"
    return "other"


def abc_triplet_factorizations_of_e_number(q: int):
    """
    Sucht Faktorisierungen q = x*y*z mit
    x ≡ 5 mod 12, y ≡ 7 mod 12, z ≡ 11 mod 12.
    Sehr brute-force, daher nur für kleine q geeignet.
    """
    ds = divisors_of_n(q)
    sols = []

    for x in ds:
        if residue_family_12(x) != "a":
            continue
        rem1 = q // x
        if x * rem1 != q:
            continue

        ds2 = divisors_of_n(rem1)
        for y in ds2:
            if residue_family_12(y) != "b":
                continue
            z = rem1 // y
            if y * z != rem1:
                continue
            if residue_family_12(z) != "c":
                continue
            sols.append((x, y, z))

    return sols


def e_factor_pairs(n: int):
    """
    Alle Faktorpaarungen n = u*v mit u <= v und u,v ≡ 1 mod 12.
    """
    pairs = []
    for d in divisors_of_n(n):
        q = n // d
        if d > q:
            continue
        if d % 12 == 1 and q % 12 == 1:
            pairs.append((d, q))
    return pairs


def best_e_balance_score(n: int, pairs):
    """
    Liefert:
    - best_e_balance: je kleiner |log u - log v|, desto symmetrischer
    - sqrt_e_distance: normierter Abstand von sqrt(n)
    """
    if not pairs:
        return {
            "has_e_pair": 0.0,
            "best_e_balance": 0.0,
            "sqrt_e_distance": 0.0,
            "best_pair": None,
        }

    s = math.sqrt(n)
    best_pair = None
    best_balance = None
    best_dist = None

    for u, v in pairs:
        balance = abs(math.log(u) - math.log(v))
        dist = abs(v - s) / max(1.0, s)

        if best_balance is None or balance < best_balance:
            best_balance = balance
            best_dist = dist
            best_pair = (u, v)

    return {
        "has_e_pair": 1.0,
        "best_e_balance": float(best_balance),
        "sqrt_e_distance": float(best_dist),
        "best_pair": best_pair,
    }


def e_quadruple_witnesses_from_pairs(n: int, pairs):
    """
    Sehr einfache Quattropel-Zeugen:
    Wenn (u1,v1) und (u2,v2) zwei E-Faktorpaarungen von n sind, dann
    u1*v1 = u2*v2 = n.
    """
    if len(pairs) < 2:
        return {
            "has_quadruple_witness": 0.0,
            "num_quadruple_witnesses": 0.0,
            "quadruple_symmetry_score": 0.0,
        }

    witnesses = 0
    best_score = None

    for i in range(len(pairs)):
        for j in range(i + 1, len(pairs)):
            (u1, v1) = pairs[i]
            (u2, v2) = pairs[j]

            witnesses += 1

            score = (
                abs(math.log(u1) - math.log(u2))
                + abs(math.log(v1) - math.log(v2))
            )

            if best_score is None or score < best_score:
                best_score = score

    return {
        "has_quadruple_witness": 1.0,
        "num_quadruple_witnesses": float(witnesses),
        "quadruple_symmetry_score": float(best_score if best_score is not None else 0.0),
    }


def e_quadruple_features(n: int):
    pairs = e_factor_pairs(n)
    bal = best_e_balance_score(n, pairs)
    wit = e_quadruple_witnesses_from_pairs(n, pairs)

    return {
        "has_e_pair": bal["has_e_pair"],
        "num_e_factor_pairs": float(len(pairs)),
        "best_e_balance": bal["best_e_balance"],
        "sqrt_e_distance": bal["sqrt_e_distance"],
        "has_quadruple_witness": wit["has_quadruple_witness"],
        "num_quadruple_witnesses": wit["num_quadruple_witnesses"],
        "quadruple_symmetry_score": wit["quadruple_symmetry_score"],
    }


# ============================================================
# Summe-zweier-Quadrate-Features
# ============================================================

def prime_factorization(n: int):
    """
    Sehr einfache Primfaktorzerlegung für moderate n.
    Rückgabe: dict p -> Exponent
    """
    fac = {}
    m = int(n)
    d = 2
    while d * d <= m:
        while m % d == 0:
            fac[d] = fac.get(d, 0) + 1
            m //= d
        d = 3 if d == 2 else d + 2
    if m > 1:
        fac[m] = fac.get(m, 0) + 1
    return fac


def is_sum_of_two_squares(n: int) -> bool:
    """
    Fermat-Euler:
    n = x^2 + y^2 genau dann, wenn jeder Primfaktor p ≡ 3 mod 4
    mit geradem Exponenten vorkommt.
    """
    if n < 0:
        return False
    if n == 0:
        return True
    fac = prime_factorization(n)
    for p, e in fac.items():
        if p % 4 == 3 and (e % 2 == 1):
            return False
    return True


def num_sum_of_two_squares_repr(n: int):
    """
    Zählt Darstellungen n = x^2 + y^2 mit 0 <= x <= y.
    Für moderate n völlig ausreichend.
    """
    if n < 0:
        return 0
    count = 0
    r = int(math.isqrt(n))
    for x in range(r + 1):
        y2 = n - x * x
        y = int(math.isqrt(y2))
        if y >= x and y * y == y2:
            count += 1
    return count


def s2_features(n: int):
    """
    Basale Quadratsummen-Features.
    """
    flag = 1.0 if is_sum_of_two_squares(n) else 0.0
    count = float(num_sum_of_two_squares_repr(n)) if flag > 0.5 else 0.0
    return {
        "is_sum_of_two_squares": flag,
        "num_sum_of_two_squares_repr": count,
    }


def s2_e_pair_features(n: int):
    """
    Prüft die E-Faktorpaarungen n = u*v mit u,v ≡ 1 mod 12 zusätzlich
    auf Quadratsummen-Eigenschaft.
    """
    pairs = e_factor_pairs(n)
    if not pairs:
        return {
            "num_s2_e_factor_pairs": 0.0,
            "best_e_pair_s2_flag": 0.0,
        }

    num_s2 = 0
    best_pair = best_e_balance_score(n, pairs)["best_pair"]
    best_flag = 0.0

    for u, v in pairs:
        if is_sum_of_two_squares(u) and is_sum_of_two_squares(v):
            num_s2 += 1

    if best_pair is not None:
        u, v = best_pair
        if is_sum_of_two_squares(u) and is_sum_of_two_squares(v):
            best_flag = 1.0

    return {
        "num_s2_e_factor_pairs": float(num_s2),
        "best_e_pair_s2_flag": float(best_flag),
    }


# ============================================================
# Erweiterte Profilanalyse
# ============================================================

def stabilized_e_quadruple_features(n: int):
    """
    Numerisch stabilisierte Version:
    - Flags separat
    - Zählwerte log-skaliert
    - keine extremen Sentinel-Werte
    """
    eq = e_quadruple_features(n)

    def safe_log1p(x: float) -> float:
        return math.log1p(max(0.0, x))

    def soft_clip(x: float, scale: float = 4.0) -> float:
        # Sättigt große Werte weich gegen 1
        return x / (scale + abs(x))

    return {
        "has_e_pair": eq["has_e_pair"],
        "log_num_e_factor_pairs": safe_log1p(eq["num_e_factor_pairs"]),
        "best_e_balance": soft_clip(eq["best_e_balance"], scale=4.0),
        "sqrt_e_distance": soft_clip(eq["sqrt_e_distance"], scale=6.0),
        "has_quadruple_witness": eq["has_quadruple_witness"],
        "log_num_quadruple_witnesses": safe_log1p(eq["num_quadruple_witnesses"]),
        "quadruple_symmetry_score": soft_clip(eq["quadruple_symmetry_score"], scale=5.0),
    }


def stabilized_s2_features(n: int):
    """
    Numerisch stabilisierte S2-Features.
    """
    s2 = s2_features(n)
    s2e = s2_e_pair_features(n)
    return {
        "is_sum_of_two_squares": s2["is_sum_of_two_squares"],
        "log_num_sum_of_two_squares_repr": math.log1p(s2["num_sum_of_two_squares_repr"]),
        "log_num_s2_e_factor_pairs": math.log1p(s2e["num_s2_e_factor_pairs"]),
        "best_e_pair_s2_flag": s2e["best_e_pair_s2_flag"],
    }


# ---------------------------------------------------------------------
# Prime/Composite-Gate
# ---------------------------------------------------------------------

NEUTRAL_PRIME_LIKE_LABEL = "neutral-prime-like"


def _safe_get(d: Dict[str, Any], key: str, default: float = 0.0) -> float:
    v = d.get(key, default)
    try:
        return float(v)
    except Exception:
        return default


def flatten_profile(profile: Dict[str, Any]) -> Dict[str, float]:
    """Flattens the nested profile dict for gate consumption."""
    out = {}
    
    # Absorption
    for k, v in profile["absorption"].items():
        if k == "mean_steps":
            out["mean_steps"] = float(v)
        else:
            out[f"abs_{k}"] = float(v)
            
    # QSD
    for k, v in profile["qsd_macro"].items():
        out[f"qsd_{k}"] = float(v)
        
    # Eigenvalues
    for i, val in enumerate(profile["eigenvalues"]):
        out[f"eig_{i+1}"] = float(abs(val))
        
    # E-Features
    if "e_features" in profile:
        out.update(profile["e_features"])
        
    # S2-Features
    if "s2_features" in profile:
        out.update(profile["s2_features"])
        
    return out


def compute_prime_composite_gate(features: Dict[str, float]) -> Dict[str, float]:
    """
    Liefert weiche Gate-Scores:

      - prime_score:    Evidenz für neutral-prime-like
      - composite_score:Evidenz für zusammengesetzte Struktur

    Die Idee ist bewusst heuristisch und soll NICHT die eigentliche
    Geometrie ersetzen, sondern nur grobe Fehlentscheidungen abfedern.
    """
    abs0 = _safe_get(features, "abs_0")
    absE = _safe_get(features, "abs_E")
    absA = _safe_get(features, "abs_A")
    absB = _safe_get(features, "abs_B")
    absC = _safe_get(features, "abs_C")
    mean_steps = _safe_get(features, "mean_steps")

    # E-Quattropel-Merkmale
    has_e_pair = _safe_get(features, "has_e_pair")
    log_num_e_factor_pairs = _safe_get(features, "log_num_e_factor_pairs")
    best_e_balance = _safe_get(features, "best_e_balance")
    sqrt_e_distance = _safe_get(features, "sqrt_e_distance")
    has_quadruple_witness = _safe_get(features, "has_quadruple_witness")
    log_num_quadruple_witnesses = _safe_get(features, "log_num_quadruple_witnesses")
    quadruple_symmetry_score = _safe_get(features, "quadruple_symmetry_score")

    # Quadratsummen-Merkmale
    is_sum_of_two_squares = _safe_get(features, "is_sum_of_two_squares")
    log_num_sum_of_two_squares_repr = _safe_get(features, "log_num_sum_of_two_squares_repr")
    log_num_s2_e_factor_pairs = _safe_get(features, "log_num_s2_e_factor_pairs")
    best_e_pair_s2_flag = _safe_get(features, "best_e_pair_s2_flag")

    # QSD / Spektral grob
    qsdA = _safe_get(features, "qsd_A")
    qsdB = _safe_get(features, "qsd_B")
    qsdN = _safe_get(features, "qsd_N")
    qsdU = _safe_get(features, "qsd_U")
    eig1 = _safe_get(features, "eig_1")
    eig2 = _safe_get(features, "eig_2")
    eig7 = _safe_get(features, "eig_7")
    eig8 = _safe_get(features, "eig_8")

    # -----------------------------
    # Prime-like Evidenz
    # -----------------------------
    #
    # Neutral-prime-like scheint sich in deinen Läufen vor allem dadurch
    # auszuzeichnen, dass:
    #   - wenig klare E-/Quattropel-Zeugen existieren,
    #   - keine starke Quadratsummen-/e-Paar-Struktur vorliegt,
    #   - die Absorption vergleichsweise "neutral" bleibt.
    #
    prime_score = 0.0
    prime_score += 1.10 * abs0
    prime_score += 0.55 * qsdN
    prime_score += 0.35 * qsdU
    prime_score += 0.25 * (eig1 + eig2 + eig7 + eig8)

    # Strafterm für klare zusammengesetzte Struktur
    prime_score -= 1.30 * has_e_pair
    prime_score -= 0.70 * has_quadruple_witness
    prime_score -= 0.60 * is_sum_of_two_squares * best_e_pair_s2_flag
    prime_score -= 0.45 * log_num_e_factor_pairs
    prime_score -= 0.35 * log_num_quadruple_witnesses
    prime_score -= 0.30 * log_num_s2_e_factor_pairs

    # Kleine Belohnung, wenn wirklich "nichts" an e-Struktur da ist
    if has_e_pair < 0.5 and has_quadruple_witness < 0.5:
        prime_score += 0.45

    # -----------------------------
    # Composite-like Evidenz
    # -----------------------------
    composite_score = 0.0
    composite_score += 1.30 * has_e_pair
    composite_score += 1.00 * has_quadruple_witness
    composite_score += 0.75 * log_num_e_factor_pairs
    composite_score += 0.60 * log_num_quadruple_witnesses
    composite_score += 0.50 * is_sum_of_two_squares * best_e_pair_s2_flag
    composite_score += 0.35 * log_num_s2_e_factor_pairs

    # "gute" e-Paarung -> stärkere composite-Evidenz
    if best_e_balance > 0.0:
        composite_score += 0.35 / (1.0 + best_e_balance)
    if sqrt_e_distance > 0.0:
        composite_score += 0.35 / (1.0 + sqrt_e_distance)
    if quadruple_symmetry_score > 0.0:
        composite_score += 0.20 / (1.0 + quadruple_symmetry_score)

    # Leichte allgemeine Strukturbeiträge
    composite_score += 0.20 * absE
    composite_score += 0.10 * (absA + absB + absC)
    composite_score += 0.10 * mean_steps

    return {
        "prime_score": float(prime_score),
        "composite_score": float(composite_score),
        "gate_margin": float(prime_score - composite_score),
    }


def gate_adjusted_distances(
    raw_distances: Dict[str, float],
    features: Dict[str, float],
) -> Dict[str, float]:
    """
    Modifiziert die Klassendistanzen sanft per Prime/Composite-Gate.
    """
    d = dict(raw_distances)
    gate = compute_prime_composite_gate(features)
    margin = gate["gate_margin"]

    # Tuning-Konstanten:
    PRIME_PUSH = 0.85
    PRIME_PULL = 0.55
    COMP_PUSH = 0.95
    COMP_PULL = 0.60
    MIXED_DAMP = 0.20

    family_labels = [
        "family-3",
        "family-5",
        "family-7",
        "family-11",
        "family-13",
        "family-11-13",
    ]

    if margin > 0.0:
        # prime-like bevorzugen
        if NEUTRAL_PRIME_LIKE_LABEL in d:
            d[NEUTRAL_PRIME_LIKE_LABEL] = max(
                0.0,
                d[NEUTRAL_PRIME_LIKE_LABEL] - PRIME_PUSH * margin
            )
        for lbl in family_labels:
            if lbl in d:
                d[lbl] += PRIME_PULL * margin
        if "mixed" in d:
            d["mixed"] += MIXED_DAMP * margin
    else:
        # composite-like bevorzugen
        comp_margin = -margin
        if NEUTRAL_PRIME_LIKE_LABEL in d:
            d[NEUTRAL_PRIME_LIKE_LABEL] += COMP_PUSH * comp_margin
        for lbl in family_labels:
            if lbl in d:
                d[lbl] = max(0.0, d[lbl] - COMP_PULL * comp_margin)
        if "mixed" in d:
            d["mixed"] = max(0.0, d["mixed"] - 0.15 * comp_margin)

    return d


MACRO_KEYS_EXTENDED = ["A", "B", "D", "M1", "M2", "M3", "N", "U", "Z+", "Z-", "Z0"]


def profile_vector_v2(profile: dict) -> np.ndarray:
    """
    Baut einen numerischen Profilvektor aus:
    - Absorptionsprofil
    - mean_steps
    - QSD-Makroprofil
    - führenden Eigenwerten (Beträge)
    """
    abs_part = [
        profile["absorption"]["A"],
        profile["absorption"]["B"],
        profile["absorption"]["C"],
        profile["absorption"]["E"],
        profile["absorption"]["0"],
        profile["absorption"]["mean_steps"],
    ]

    qsd_part = [profile["qsd_macro"].get(k, 0.0) for k in MACRO_KEYS_EXTENDED]
    eig_part = [abs(z) for z in profile["eigenvalues"]]
    return np.array(abs_part + qsd_part + eig_part, dtype=float)


def profile_vector_v3(profile: dict, H_inf: float = 20.0) -> np.ndarray:
    abs_part = [
        profile["absorption"]["A"],
        profile["absorption"]["B"],
        profile["absorption"]["C"],
        profile["absorption"]["E"],
        profile["absorption"]["0"],
        profile["absorption"]["mean_steps"],
    ]

    qsd_part = [profile["qsd_macro"].get(k, 0.0) for k in MACRO_KEYS_EXTENDED]
    eig_part = [abs(z) for z in profile["eigenvalues"]]

    hz = H_inf * (
        profile["qsd_macro"].get("Z+", 0.0)
        + profile["qsd_macro"].get("Z-", 0.0)
        + profile["qsd_macro"].get("Z0", 0.0)
    )
    h0 = H_inf * profile["absorption"]["0"]
    hs = H_inf * profile["absorption"]["mean_steps"]

    # Stabilisierte E-Quattropel-Features
    eq = stabilized_e_quadruple_features(profile["N"])
    e_part_raw = [
        eq["has_e_pair"],
        eq["log_num_e_factor_pairs"],
        eq["best_e_balance"],
        eq["sqrt_e_distance"],
        eq["has_quadruple_witness"],
        eq["log_num_quadruple_witnesses"],
        eq["quadruple_symmetry_score"],
    ]

    # Stabilisierte S2-Features
    s2 = stabilized_s2_features(profile["N"])
    s2_part_raw = [
        s2["is_sum_of_two_squares"],
        s2["log_num_sum_of_two_squares_repr"],
        s2["log_num_s2_e_factor_pairs"],
        s2["best_e_pair_s2_flag"],
    ]

    # Sanfte Blockgewichtung
    # Stark zurückgenommen: nur noch Zusatzsignal, nicht Leitgeometrie
    E_BLOCK_WEIGHT = 0.10
    S2_BLOCK_WEIGHT = 0.05

    e_part = [E_BLOCK_WEIGHT * x for x in e_part_raw]
    s2_part = [S2_BLOCK_WEIGHT * x for x in s2_part_raw]

    return np.array(
        abs_part + qsd_part + eig_part + [hz, h0, hs] + e_part + s2_part,
        dtype=float
    )


def build_profile_matrix(profiles: list, H_inf: float = 20.0):
    Ns = [p["N"] for p in profiles]
    X = np.vstack([profile_vector_v3(p, H_inf=H_inf) for p in profiles])
    return Ns, X


def build_labeled_profile_matrix(profiles, labels_map, H_inf=20.0):
    rows = []
    ys = []
    ns = []

    for p in profiles:
        N = p["N"]
        if N in labels_map:
            rows.append(profile_vector_v3(p, H_inf=H_inf))
            ys.append(labels_map[N])
            ns.append(N)

    X = np.vstack(rows)
    return ns, X, ys


def compute_class_centroids(Xs, y):
    centroids = {}
    classes = set(y)
    for c in classes:
        indices = [i for i, label in enumerate(y) if label == c]
        if indices:
            centroids[c] = np.mean(Xs[indices], axis=0)
    return centroids


def classify_by_nearest_centroid(x, centroids):
    best_c = None
    min_d = float("inf")
    for c, center in centroids.items():
        d = np.linalg.norm(x - center)
        if d < min_d:
            min_d = d
            best_c = c
    return best_c, min_d


def train_and_report_classifier(profiles, labels_map, H_inf=20.0):
    Ns, X, y = build_labeled_profile_matrix(profiles, labels_map, H_inf=H_inf)
    Xs, mu, sigma = standardize_matrix(X)
    centroids = compute_class_centroids(Xs, y)

    print("\nKlassifikator-Report")
    print("-" * 60)

    # Dictionary for fast profile access by N
    prof_map = {p["N"]: p for p in profiles}

    for i, N in enumerate(Ns):
        # Raw Classification
        raw_dists = {}
        for label, center in centroids.items():
            raw_dists[label] = float(np.linalg.norm(Xs[i] - center))

        # Gate Logic
        if N in prof_map:
            feats = flatten_profile(prof_map[N])
            gated = gate_adjusted_distances(raw_dists, feats)
            pred = min(gated, key=gated.get)
            dist = gated[pred]
        else:
            pred, dist = classify_by_nearest_centroid(Xs[i], centroids)

        ok = "OK" if pred == y[i] else "MISS"
        print(f"N={N:>4d}  true={y[i]:>12s}  pred={pred:>12s}  dist={dist:.4f}  {ok}")

    return {
        "Ns": Ns,
        "X": X,
        "Xs": Xs,
        "y": y,
        "centroids": centroids,
        "mu": mu,
        "sigma": sigma,
    }


def predict_number_family(N: int, sqrt_s: float,
                          trained: dict,
                          qparams=QModelParams(),
                          eparams=EnergyParams(),
                          rparams=RefinedIIIParams(),
                          nparams=NOperatorParams(),
                          hparams=HurwitzPressureParams()):
    profile = N_spectral_profile(
        N=N,
        sqrt_s=sqrt_s,
        qparams=qparams,
        eparams=eparams,
        rparams=rparams,
        nparams=nparams,
        hparams=hparams,
        k_eigs=8,
    )

    x = profile_vector_v3(profile, H_inf=hparams.H_inf).reshape(1, -1)
    x_scaled = (x - trained["mu"]) / trained["sigma"]

    # Raw Distances berechnen
    raw_dists = {}
    for label, center in trained["centroids"].items():
        raw_dists[label] = float(np.linalg.norm(x_scaled[0] - center))

    # Gate anwenden
    feats = flatten_profile(profile)
    gated = gate_adjusted_distances(raw_dists, feats)

    pred = min(gated, key=gated.get)
    return profile, pred, gated[pred]


def standardize_matrix(X: np.ndarray):
    mu = X.mean(axis=0)
    sigma = X.std(axis=0)
    sigma[sigma == 0] = 1.0
    Xs = (X - mu) / sigma
    return Xs, mu, sigma


def pca_2d(X: np.ndarray):
    """
    PCA über SVD, Rückgabe:
    - Y: 2D-Projektion
    - singular values
    - explained variance ratio
    """
    Xc = X - X.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    _ = U

    Y = Xc @ Vt[:2].T
    variances = (S ** 2) / max(1, (X.shape[0] - 1))
    total_var = variances.sum()
    evr = variances / total_var if total_var > 0 else variances
    return Y, S, evr, Vt


def plot_pca(Ns, X, title="PCA der Zahlenprofile"):
    Xs, _, _ = standardize_matrix(X)
    Y, S, evr, Vt = pca_2d(Xs)
    _ = S, Vt

    plt.figure(figsize=(9, 7))
    plt.scatter(Y[:, 0], Y[:, 1])

    for i_idx, N in enumerate(Ns):
        plt.text(Y[i_idx, 0], Y[i_idx, 1], str(N), fontsize=10)

    plt.xlabel(f"PC1 ({evr[0]*100:.1f}% Varianz)")
    plt.ylabel(f"PC2 ({evr[1]*100:.1f}% Varianz)")
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return Y, evr


def pairwise_distances_from_matrix(X: np.ndarray):
    n = X.shape[0]
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            D[i, j] = np.linalg.norm(X[i] - X[j])
    return D


def plot_distance_heatmap(Ns, D, title="Distanz-Heatmap"):
    plt.figure(figsize=(8, 7))
    plt.imshow(D, aspect="auto")
    plt.colorbar(label="euklidische Distanz")
    plt.xticks(range(len(Ns)), Ns, rotation=45)
    plt.yticks(range(len(Ns)), Ns)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def plot_dendrogram_from_distances(Ns, D, title="Hierarchisches Clustering"):
    if not SCIPY_AVAILABLE:
        print("SciPy nicht verfügbar: Dendrogramm wird übersprungen.")
        return

    condensed = squareform(D, checks=False)
    Z = linkage(condensed, method="average")

    plt.figure(figsize=(10, 5))
    dendrogram(Z, labels=[str(N) for N in Ns])
    plt.title(title)
    plt.ylabel("Cluster-Distanz")
    plt.tight_layout()
    plt.show()


def number_label(N: int) -> str:
    """
    Rein heuristische Beschriftung für die Plot-Ausgabe.
    """
    flags = []
    if N % 2 == 0:
        flags.append("2")
    if N % 3 == 0:
        flags.append("3")
    if N % 5 == 0:
        flags.append("5")
    if N % 7 == 0:
        flags.append("7")
    if not flags:
        return f"{N} [sonst]"
    return f"{N} [{' '.join(flags)}]"


def plot_pca_with_labels(Ns, X, title="PCA der Zahlenprofile"):
    Xs, _, _ = standardize_matrix(X)
    Y, S, evr, Vt = pca_2d(Xs)
    _ = S, Vt

    plt.figure(figsize=(10, 7))
    plt.scatter(Y[:, 0], Y[:, 1])

    for i_idx, N in enumerate(Ns):
        plt.text(Y[i_idx, 0], Y[i_idx, 1], number_label(N), fontsize=9)

    plt.xlabel(f"PC1 ({evr[0]*100:.1f}% Varianz)")
    plt.ylabel(f"PC2 ({evr[1]*100:.1f}% Varianz)")
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return Y, evr


def analyze_number_profiles(profiles, title_prefix=""):
    Ns, X = build_profile_matrix(profiles)

    print("\nProfilmatrix:")
    print("Zahlen:", Ns)
    print("Shape:", X.shape)

    Y, evr = plot_pca_with_labels(Ns, X, title=f"{title_prefix} PCA")
    print("\nErklärte Varianz (erste Komponenten):")
    for i, val in enumerate(evr[:5]):
        print(f"PC{i+1}: {val:.4f}")

    Xs, _, _ = standardize_matrix(X)
    D = pairwise_distances_from_matrix(Xs)
    plot_distance_heatmap(Ns, D, title=f"{title_prefix} Distanz-Heatmap")
    plot_dendrogram_from_distances(Ns, D, title=f"{title_prefix} Dendrogramm")

    return Ns, X, D, Y, evr


# ============================================================
# Feature-Namen für PCA-Loadings
# ============================================================

def feature_names_for_profiles(X_shape_1: int):
    base = [
        "abs_A", "abs_B", "abs_C", "abs_E", "abs_0", "mean_steps",
        "qsd_A", "qsd_B", "qsd_D", "qsd_M1", "qsd_M2", "qsd_M3",
        "qsd_N", "qsd_U", "qsd_Z+", "qsd_Z-", "qsd_Z0",
    ]
    # 3 Hurwitz + 7 E-Quattropel + 4 S2 = 14 Merkmale
    n_special = 14
    n_eigs = X_shape_1 - len(base) - n_special
    eigs = [f"eig_{idx+1}" for idx in range(n_eigs)]
    hurwitz = ["hurwitz_z", "hurwitz_0", "hurwitz_steps"]
    equad = [
        "has_e_pair", "log_num_e_pairs", "best_e_bal", "sqrt_e_dist",
        "has_quad_wit", "log_num_quad_wit", "quad_sym_score"
    ]
    s2_names = [
        "is_s2", "log_num_s2_repr", "log_num_s2_e_pairs", "best_e_pair_s2"
    ]
    return base + eigs + hurwitz + equad + s2_names


def print_pca_loadings(X, Vt, n_components=3, top_k=12):
    names = feature_names_for_profiles(X.shape[1])

    for comp in range(n_components):
        print(f"\nLoadings PC{comp+1}:")
        pairs = list(zip(names, Vt[comp]))
        pairs = sorted(pairs, key=lambda kv: -abs(kv[1]))
        for name, val in pairs[:top_k]:
            print(f"{name:>12s}: {val: .4f}")


def print_pca_coordinates_3d(Ns, X):
    Xs, _, _ = standardize_matrix(X)
    Y, S, evr, Vt = pca_2d(Xs)
    _ = Y, S

    Xc = Xs - Xs.mean(axis=0, keepdims=True)
    if Vt.shape[0] < 3:
        U, Sfull, Vtfull = np.linalg.svd(Xc, full_matrices=False)
        _ = U
        Vt_use = Vtfull
        Y3 = Xc @ Vt_use[:3].T
        evr_full = (Sfull**2) / max(1, (Xc.shape[0] - 1))
        evr_full = evr_full / evr_full.sum()
    else:
        Vt_use = Vt
        Y3 = Xc @ Vt_use[:3].T
        evr_full = evr

    print("\nPCA-Koordinaten:")
    for idx, number in enumerate(Ns):
        print(f"{number:>4d}: PC1={Y3[idx,0]: .4f}, PC2={Y3[idx,1]: .4f}, PC3={Y3[idx,2]: .4f}")

    print("\nErklärte Varianz (erste Komponenten):")
    for idx, val in enumerate(evr_full[:5]):
        print(f"PC{idx+1}: {val:.4f}")

    print_pca_loadings(Xs, Vt_use, n_components=3, top_k=14)

    return Y3, evr_full, Vt_use


def family_mean_vectors(profiles, families: dict, H_inf: float = 1.0):
    prof_map = {p["N"]: p for p in profiles}
    result = {}

    for name, nums in families.items():
        vecs = [profile_vector_v3(prof_map[n], H_inf=H_inf) for n in nums if n in prof_map]
        if vecs:
            result[name] = np.mean(np.vstack(vecs), axis=0)

    return result


def print_family_centers(profiles):
    families = {
        "prime": [13, 17, 19, 23, 29, 31, 37, 41, 43],
        "3-family": [15, 21, 33, 39, 51, 57, 69],
        "5-family": [35, 55, 65, 85, 95],
        "7-family": [77, 91, 119, 133],
    }

    prof_map = {p["N"]: p for p in profiles}
    print("\nFamilien-Mittelwerte (nur Absorption + mean_steps):")
    for fam, nums in families.items():
        avail = [prof_map[n] for n in nums if n in prof_map]
        if not avail:
            continue
        absA = np.mean([p["absorption"]["A"] for p in avail])
        absB = np.mean([p["absorption"]["B"] for p in avail])
        absC = np.mean([p["absorption"]["C"] for p in avail])
        absE = np.mean([p["absorption"]["E"] for p in avail])
        abs0 = np.mean([p["absorption"]["0"] for p in avail])
        steps = np.mean([p["absorption"]["mean_steps"] for p in avail])
        print(f"{fam:>8s}: A={absA:.4f}, B={absB:.4f}, C={absC:.4f}, E={absE:.4f}, 0={abs0:.4f}, steps={steps:.4f}")


# ============================================================
# Demo
# ============================================================

TRAIN_LABELS = {
    # Version V.1: 13 taxonomisch als family-13 auffassen
    13: "family-13",
    17: "neutral-prime-like", 19: "neutral-prime-like", 23: "neutral-prime-like",
    29: "neutral-prime-like", 31: "neutral-prime-like", 37: "neutral-prime-like", 41: "neutral-prime-like", 43: "neutral-prime-like",

    21: "family-3", 33: "family-3", 39: "family-3", 51: "family-3", 57: "family-3", 69: "family-3",

    55: "family-5", 65: "family-5", 85: "family-5", 95: "family-5",

    77: "family-7", 91: "family-7", 119: "family-7", 133: "family-7",

    # neue 11/13-Familien
    121: "family-11",
    169: "family-13",
    143: "family-11-13",
    187: "family-11",
    209: "family-11",
    221: "family-13",

    # Version V: mixed robuster trainieren
    15: "mixed", 35: "mixed", 45: "mixed", 63: "mixed", 75: "mixed", 105: "mixed",
}

if __name__ == "__main__":
    hparams = HurwitzPressureParams(H_inf=20.0)

    test_numbers = [
        13, 17, 19, 23, 29, 31, 37, 41, 43,
        15, 21, 33, 39, 45, 51, 57, 63, 69, 75, 105,
        35, 55, 65, 85, 95,
        77, 91, 119, 133,
        121, 143, 169, 187, 209, 221
    ]

    profiles = compare_numbers(
        test_numbers,
        sqrt_s=3.5,
        hparams=hparams
    )
    print_number_profiles(profiles)

    trained = train_and_report_classifier(
        profiles,
        TRAIN_LABELS,
        H_inf=hparams.H_inf
    )

    for N_test in [27, 45, 63, 75, 105, 143, 121, 169, 187, 209, 221]:
        prof, pred, dist = predict_number_family(
            N_test,
            sqrt_s=3.5,
            trained=trained,
            hparams=hparams
        )
        print(f"\nTestzahl {N_test}: pred={pred}, dist={dist:.4f}")

    # S2-Kurztest
    test_vals = [121, 143, 169, 385, 1001, 105]
    print("\nE-Quattropel-Kurztest")
    print("-" * 60)
    for n in test_vals:
        feat = stabilized_e_quadruple_features(n)
        print(f"N={n}")
        for k, v in feat.items():
            print(f"  {k}: {v}")
        print()

    print("\nS2-Kurztest")
    print("-" * 60)
    for n in test_vals:
        feat = stabilized_s2_features(n)
        print(f"N={n}")
        for k, v in feat.items():
            print(f"  {k}: {v}")
        print()
