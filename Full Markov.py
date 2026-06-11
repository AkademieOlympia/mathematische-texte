from itertools import product
from collections import Counter, defaultdict
from dataclasses import dataclass
import numpy as np
import math
import random

# ============================================================
# Alphabet und absorbierende Kanäle
# ============================================================

ALPHABET = ("E", "A", "B", "C")
ABSORBING = ["A", "B", "C", "E", "0"]

WORDS = [''.join(p) for p in product(ALPHABET, repeat=4)]
ALL_STATES = WORDS + ABSORBING

WORD_INDEX = {w: i for i, w in enumerate(WORDS)}
STATE_INDEX = {s: i for i, s in enumerate(ALL_STATES)}

N_WORDS = len(WORDS)              # 256
N_ABS = len(ABSORBING)            # 5
N_TOTAL = len(ALL_STATES)         # 261


# ============================================================
# Grundfunktionen
# ============================================================

def clip(x, lo=0.0, hi=1.0):
    return max(lo, min(hi, x))

def sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))

def profile_of_word(w: str):
    return tuple(sorted(Counter(w).values(), reverse=True))

def content_key(w: str):
    cnt = Counter(w)
    return tuple(sorted(cnt.items(), key=lambda kv: (kv[0], kv[1])))

def cyclic_rotations(w: str):
    return [w[i:] + w[:i] for i in range(len(w))]

def reverse_word(w: str):
    return w[::-1]

def swap_adjacent(w: str, i: int):
    s = list(w)
    s[i], s[i+1] = s[i+1], s[i]
    return ''.join(s)

def all_local_neighbors(w: str):
    """
    Lokale Nachbarn durch:
    - drei benachbarte Swaps
    - reverse
    - rotate left
    - rotate right
    """
    nbrs = set()
    for i in range(3):
        nbrs.add(swap_adjacent(w, i))
    nbrs.add(reverse_word(w))
    nbrs.add(w[1:] + w[0])
    nbrs.add(w[-1] + w[:-1])
    nbrs.discard(w)
    return sorted(nbrs)


# ============================================================
# Ausgezeichnete zyklische Klassen
# ============================================================

CYCLIC_PLUS = {"EABC", "ABCE", "BCEA", "CEAB"}
CYCLIC_MINUS = {"EACB", "ACBE", "CBEA", "BEAC"}


# ============================================================
# Makroklassifikation
# ============================================================

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

    block = {x+x+y+y, y+y+x+x}
    alt = {x+y+x+y, y+x+y+x}
    mirror = {x+y+y+x, y+x+x+y}

    if w in mirror:
        return "N"   # neutral / Spiegel
    elif w in alt:
        return "A"   # Alternation
    elif w in block:
        return "B"   # Block
    else:
        raise RuntimeError(f"Unerwartete (2,2)-Struktur für {w}")

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
        elif w in CYCLIC_MINUS:
            return "Z-"
        else:
            return "Z0"

    raise RuntimeError(f"Unbekanntes Profil für {w}")


# ============================================================
# Morley / Walter / Aktivität / Signatur
# ============================================================

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
        # einfache Projektion auf erstes Symbol der Paarung
        return syms[0] if syms[0] in "EABC" else "mix"

    if cls == "B":
        syms = sorted(cnt.keys())
        return "0" if activity_Q(w, params) < 0.25 else syms[0]

    if cls == "N":
        return "0"

    if cls == "M":
        # doppelt auftretendes Symbol als Hauptprojektion
        double_symbol = [k for k, v in cnt.items() if v == 2][0]
        return double_symbol

    if cls == "D":
        triple_symbol = [k for k, v in cnt.items() if v == 3][0]
        return triple_symbol

    if cls == "U":
        return "0" if w[0] == "E" else w[0]

    if cls == "Z0":
        # einfache Heuristik: letztes Zeichen
        return w[-1]

    return "0"


# ============================================================
# Energieparameter
# ============================================================

@dataclass
class EnergyParams:
    theta_p: float = 3.0   # Pinski-Schwelle
    kappa_p: float = 2.0
    theta_d: float = 5.0   # direkte Projektion
    kappa_d: float = 1.3

def alpha_pinski(sqrt_s: float, p: EnergyParams) -> float:
    return sigmoid(p.kappa_p * (sqrt_s - p.theta_p))

def gamma_detector(sqrt_s: float, p: EnergyParams) -> float:
    return sigmoid(p.kappa_d * (sqrt_s - p.theta_d))

def beta_neutral(sqrt_s: float, p: EnergyParams) -> float:
    return 1.0 - alpha_pinski(sqrt_s, p)


# ============================================================
# Übergangsgewichte Wort -> Wort
# ============================================================

def word_to_word_weight(w: str, w2: str, sqrt_s: float,
                        qparams: QModelParams,
                        eparams: EnergyParams) -> float:
    """
    Gewicht für interne Übergänge W -> W'
    """
    a = alpha_pinski(sqrt_s, eparams)
    b = beta_neutral(sqrt_s, eparams)

    M1, M2 = morley_M(w), morley_M(w2)
    T1, T2 = walter_T(w), walter_T(w2)
    Q1, Q2 = activity_Q(w, qparams), activity_Q(w2, qparams)
    chi1, chi2 = chirality_chi(w), chirality_chi(w2)

    c1, c2 = macro_class(w), macro_class(w2)

    # Basisterme
    weight = 0.0

    # neutrale Speicher-/Kernpfade
    if c2 in ("N", "B"):
        weight += b * (0.20 + 0.25 * T1 - 0.10 * M1)

    # gerichtete / aktive Pfade
    if c2 in ("A", "Z+", "Z-", "Z0", "M"):
        weight += a * (0.15 + 0.20 * M1 + 0.15 * Q1 - 0.10 * T1)

    # Defekt / uniforme Sektoren als schwache Streuung
    if c2 in ("D", "U"):
        weight += 0.05 + 0.05 * (1.0 - T1)

    # ähnliche Chiralität koppelt besser
    if chi1 != 0 and chi2 != 0 and chi1 == chi2:
        weight += 0.08

    # ähnliche Makroklasse koppelt etwas besser
    if c1 == c2:
        weight += 0.05

    # ähnliche Signatur koppelt besser
    s1 = sigma_projection(w, qparams)
    s2 = sigma_projection(w2, qparams)
    if s1 == s2 and s1 != "0":
        weight += 0.08

    # Walter-Kern bevorzugt interne Speicherung
    weight += 0.05 * (T1 + T2)

    return max(weight, 0.0)


# ============================================================
# Übergangsgewichte Wort -> absorbierender Kanal
# ============================================================

def word_to_absorbing_weight(w: str, absorb_state: str,
                             sqrt_s: float,
                             qparams: QModelParams,
                             eparams: EnergyParams) -> float:
    """
    Gewicht für W -> A/B/C/E/0
    """
    g = gamma_detector(sqrt_s, eparams)
    a = alpha_pinski(sqrt_s, eparams)

    M = morley_M(w)
    T = walter_T(w)
    Q = activity_Q(w, qparams)
    sig = sigma_projection(w, qparams)

    # neutraler Kanal
    if absorb_state == "0":
        return max(
            (1.0 - g) * (0.10 + 0.35 * T - 0.15 * M - 0.10 * Q)
            + 0.05 * (sig == "0"),
            0.0
        )

    # signaturtragende Kanäle
    base = g * (0.03 + 0.30 * Q + 0.20 * M - 0.12 * T) + a * (0.02 + 0.08 * M)

    if sig == absorb_state:
        base += 0.25
    elif sig != "0":
        base += 0.04

    return max(base, 0.0)


# ============================================================
# Gesamte Übergangsmatrix 261x261
# ============================================================

def build_full_transition_matrix(sqrt_s: float,
                                 qparams: QModelParams = QModelParams(),
                                 eparams: EnergyParams = EnergyParams()) -> np.ndarray:
    P = np.zeros((N_TOTAL, N_TOTAL), dtype=float)

    # Wortzustände
    for w in WORDS:
        i = STATE_INDEX[w]
        row = np.zeros(N_TOTAL, dtype=float)

        # lokale Wortnachbarn
        neighbors = all_local_neighbors(w)
        for w2 in neighbors:
            j = STATE_INDEX[w2]
            row[j] += word_to_word_weight(w, w2, sqrt_s, qparams, eparams)

        # schwacher Selbstübergang
        row[i] += 0.03 + 0.04 * walter_T(w)

        # absorbierende Kanäle
        for abs_state in ABSORBING:
            j = STATE_INDEX[abs_state]
            row[j] += word_to_absorbing_weight(w, abs_state, sqrt_s, qparams, eparams)

        # Normierung
        s = row.sum()
        if s <= 0:
            row[i] = 1.0
        else:
            row = row / s

        P[i] = row

    # absorbierende Zustände
    for astate in ABSORBING:
        i = STATE_INDEX[astate]
        P[i, i] = 1.0

    return P


# ============================================================
# Analysefunktionen
# ============================================================

def analyze_word(w: str,
                 qparams: QModelParams = QModelParams()):
    return {
        "word": w,
        "profile": profile_of_word(w),
        "class": macro_class(w),
        "M": morley_M(w),
        "chi": chirality_chi(w),
        "T": walter_T(w),
        "Q": round(activity_Q(w, qparams), 6),
        "Sigma": sigma_projection(w, qparams),
    }

def transient_indices():
    return [STATE_INDEX[w] for w in WORDS]

def absorbing_indices():
    return [STATE_INDEX[s] for s in ABSORBING]

def fundamental_matrix(P: np.ndarray):
    T_idx = transient_indices()
    A_idx = absorbing_indices()
    Q = P[np.ix_(T_idx, T_idx)]
    R = P[np.ix_(T_idx, A_idx)]
    N = np.linalg.inv(np.eye(Q.shape[0]) - Q)
    B = N @ R
    return Q, R, N, B

def exact_absorption_probabilities(P: np.ndarray, start_word: str):
    Qm, Rm, Nm, Bm = fundamental_matrix(P)
    i = WORD_INDEX[start_word]
    probs = Bm[i]
    return {ABSORBING[k]: probs[k] for k in range(len(ABSORBING))}

def exact_mean_absorption_steps(P: np.ndarray, start_word: str):
    Qm, Rm, Nm, Bm = fundamental_matrix(P)
    t = Nm @ np.ones((Nm.shape[0], 1))
    i = WORD_INDEX[start_word]
    return float(t[i, 0])

def quasi_stationary_distribution(P: np.ndarray):
    T_idx = transient_indices()
    Qm = P[np.ix_(T_idx, T_idx)]
    eigvals, eigvecs = np.linalg.eig(Qm.T)
    idx = np.argmax(np.real(eigvals))
    vec = np.real(eigvecs[:, idx])
    vec = np.abs(vec)
    vec = vec / vec.sum()
    return {WORDS[i]: vec[i] for i in range(len(WORDS))}

def simulate_path(P: np.ndarray, start_word: str, max_steps: int = 100, seed: int = 123):
    random.seed(seed)
    current = STATE_INDEX[start_word]
    path = [ALL_STATES[current]]

    for _ in range(max_steps):
        if ALL_STATES[current] in ABSORBING:
            break
        probs = P[current]
        current = random.choices(range(N_TOTAL), weights=probs, k=1)[0]
        path.append(ALL_STATES[current])

    return path


# ============================================================
# Aggregation nach Makroklassen
# ============================================================

def aggregate_qsd_by_macro(qsd_word: dict):
    agg = defaultdict(float)
    for w, p in qsd_word.items():
        agg[macro_class(w)] += p
    return dict(sorted(agg.items(), key=lambda kv: kv[0]))

def summarize_words():
    rows = [analyze_word(w) for w in WORDS]
    by_class = defaultdict(list)
    for r in rows:
        by_class[r["class"]].append(r)
    return rows, by_class


# ============================================================
# Demo
# ============================================================

if __name__ == "__main__":
    qparams = QModelParams()
    eparams = EnergyParams()

    # Beispielwörter
    examples = ["ABBA", "ABAB", "AABB", "EABC", "ABCE", "EACB", "AABC", "AAAE", "EEEE"]
    print("Einzelanalysen")
    print("-" * 70)
    for w in examples:
        print(analyze_word(w, qparams))

    # Klassenstatistik
    rows, by_class = summarize_words()
    print("\nMakroklassen")
    print("-" * 70)
    for cls in sorted(by_class.keys()):
        rs = by_class[cls]
        mean_M = np.mean([r["M"] for r in rs])
        mean_T = np.mean([r["T"] for r in rs])
        mean_Q = np.mean([r["Q"] for r in rs])
        print(f"{cls:>3s}: {len(rs):>3d} Wörter | <M>={mean_M:.3f}, <T>={mean_T:.3f}, <Q>={mean_Q:.3f}")

    # Energie-Slices
    for sqrt_s in [1.5, 3.5, 6.0]:
        print("\n" + "=" * 70)
        print(f"sqrt(s) = {sqrt_s}")

        P = build_full_transition_matrix(sqrt_s, qparams, eparams)

        for start in ["ABBA", "ABAB", "EABC"]:
            probs = exact_absorption_probabilities(P, start)
            mean_steps = exact_mean_absorption_steps(P, start)
            print(f"\nStartwort: {start}")
            print("Absorptionswahrscheinlichkeiten:", {k: round(v, 4) for k, v in probs.items()})
            print("Mittlere Absorptionsschritte:", round(mean_steps, 4))

        print("\nBeispielpfad ab EABC:")
        path = simulate_path(P, "EABC", max_steps=30, seed=42)
        print(" -> ".join(path))

        qsd = quasi_stationary_distribution(P)
        agg = aggregate_qsd_by_macro(qsd)
        print("\nQuasi-stationäre Masse nach Makroklassen:")
        for cls, val in agg.items():
            print(f"{cls:>3s}: {val:.4f}")