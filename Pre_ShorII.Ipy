import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import product

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
            row[j] += refined_word_to_word_weight(w, w2, sqrt_s, qparams, eparams, rparams)

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


def build_full_transition_matrix_refined_v3(
    sqrt_s: float,
    qparams: QModelParams = QModelParams(),
    eparams: EnergyParams = EnergyParams(),
    rparams: RefinedIIIParams = RefinedIIIParams(),
    N: int | None = None,
) -> np.ndarray:
    P = np.zeros((N_TOTAL, N_TOTAL), dtype=float)

    for w in WORDS:
        i = STATE_INDEX[w]
        row = np.zeros(N_TOTAL, dtype=float)

        nbrs = structured_neighbors(w)

        for w2 in nbrs:
            j = STATE_INDEX[w2]
            row[j] += refined_word_to_word_weight_v3(
                w, w2, sqrt_s, qparams, eparams, rparams, N=N
            )

        row[i] += (
            rparams.self_loop_base
            + rparams.self_loop_walter * refined_walter_T(w)
            + tiny_hash_bias(w)
        )

        for astate in ABSORBING:
            j = STATE_INDEX[astate]
            row[j] += refined_word_to_absorbing_weight_v3(
                w, astate, sqrt_s, qparams, eparams, rparams, N=N
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
    Erstes arithmetisches Encoding N -> Gewichte für E,A,B,C
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

    # Restklassen
    if N % 4 == 1:
        weights["A"] += 0.3
    if N % 3 == 1:
        weights["B"] += 0.3
    if N % 12 in (1, 5, 7, 11):
        weights["C"] += 0.3

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
                                           rparams=RefinedIIIParams()):
    return build_full_transition_matrix_refined_v3(
        sqrt_s=sqrt_s,
        qparams=qparams,
        eparams=eparams,
        rparams=rparams,
        N=N
    )


def build_N_dependent_transition_matrix(N: int,
                                        sqrt_s: float,
                                        qparams=QModelParams(),
                                        eparams=EnergyParams(),
                                        nparams=NOperatorParams(),
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

        # interne Wortzustände
        for w2 in WORDS:
            j = STATE_INDEX[w2]
            row[j] *= internal_factor

        # absorbierende Kanäle
        for a in ABSORBING:
            j = STATE_INDEX[a]
            if a == sig:
                row[j] *= absorb_factor
            elif a == "0" and cls == "N":
                row[j] *= 1.0 + 0.15
            else:
                row[j] *= 1.0

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
                       nparams=NOperatorParams(),
                       rparams=RefinedIIIParams(),
                       k_eigs: int = 8):
    """
    Komplettes Profil einer Zahl N.
    """
    piN = initial_distribution_from_N(N)
    PN = build_N_dependent_transition_matrix_v3(N, sqrt_s, qparams, eparams, rparams)

    absorption = absorption_profile_from_initial_distribution(PN, piN)
    eigs = leading_eigenvalues_of_transient(PN, k=k_eigs)
    qsd_macro = qsd_macro_profile(PN)

    return {
        "N": N,
        "sqrt_s": sqrt_s,
        "initial_distribution": piN,
        "absorption": absorption,
        "eigenvalues": eigs,
        "qsd_macro": qsd_macro,
    }


# ============================================================
# Vergleich mehrerer Zahlen
# ============================================================

def compare_numbers(N_list, sqrt_s: float):
    profiles = []
    for N in N_list:
        profiles.append(N_spectral_profile(N, sqrt_s))
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
# Erweiterte Profilanalyse
# ============================================================

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


def build_profile_matrix(profiles: list):
    Ns = [p["N"] for p in profiles]
    X = np.vstack([profile_vector_v2(p) for p in profiles])
    return Ns, X


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

    for i, N in enumerate(Ns):
        plt.text(Y[i, 0], Y[i, 1], str(N), fontsize=10)

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

    for i, N in enumerate(Ns):
        plt.text(Y[i, 0], Y[i, 1], number_label(N), fontsize=9)

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
    n_eigs = X_shape_1 - len(base)
    eigs = [f"eig_{idx+1}" for idx in range(n_eigs)]
    return base + eigs


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


def family_mean_vectors(profiles, families: dict):
    prof_map = {p["N"]: p for p in profiles}
    result = {}

    for name, nums in families.items():
        vecs = [profile_vector_v2(prof_map[n]) for n in nums if n in prof_map]
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

if __name__ == "__main__":
    extended_numbers = [
        13, 17, 19, 23, 29, 31, 37, 41, 43,
        15, 21, 33, 39, 51, 57, 69,
        35, 55, 65, 85, 95,
        77, 91, 119, 133
    ]

    extended_profiles = compare_numbers(extended_numbers, sqrt_s=3.5)

    print_number_profiles(extended_profiles)

    Ns_ext, X_ext = build_profile_matrix(extended_profiles)
    Y3, evr_ext, Vt_ext = print_pca_coordinates_3d(Ns_ext, X_ext)
    _ = Y3, evr_ext, Vt_ext

    Xs_ext, _, _ = standardize_matrix(X_ext)
    D_ext = pairwise_distances_from_matrix(Xs_ext)

    print("\nPaarweise Profildistanzen (großer Testkorpus):")
    print("      " + " ".join(f"{N:>6d}" for N in Ns_ext))
    for idx, number in enumerate(Ns_ext):
        row = " ".join(f"{D_ext[idx,j]:6.3f}" for j in range(len(Ns_ext)))
        print(f"{number:>5d} {row}")

    plot_pca_with_labels(Ns_ext, X_ext, title="Erweiterter Zahlentest: PCA")
    plot_distance_heatmap(Ns_ext, D_ext, title="Erweiterter Zahlentest: Distanz-Heatmap")
    plot_dendrogram_from_distances(Ns_ext, D_ext, title="Erweiterter Zahlentest: Dendrogramm")