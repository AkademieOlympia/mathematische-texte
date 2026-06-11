from itertools import product
from collections import Counter, defaultdict
from dataclasses import dataclass
import math
import numpy as np

ALPHABET = ("E", "A", "B", "C")


# ============================================================
# Hilfsfunktionen
# ============================================================

def clip(x, lo=0.0, hi=1.0):
    return max(lo, min(hi, x))


def word_reverse(w: str) -> str:
    return w[::-1]


def cyclic_rotations(w: str):
    return [w[i:] + w[:i] for i in range(len(w))]


def profile_of_word(w: str):
    """
    Gibt das Besetzungsprofil als sortiertes Tupel zurück, z.B.
    EABC -> (1,1,1,1)
    ABBA -> (2,2)
    AABC -> (2,1,1)
    """
    counts = sorted(Counter(w).values(), reverse=True)
    return tuple(counts)


def content_key(w: str):
    """
    Inhaltsklasse als sortierte Paare, z.B.
    ABBA -> (('A',2),('B',2))
    EABC -> (('A',1),('B',1),('C',1),('E',1))
    """
    cnt = Counter(w)
    return tuple(sorted(cnt.items(), key=lambda kv: (kv[0], kv[1])))


def is_palindromic(w: str) -> bool:
    return w == w[::-1]


def unique_symbols(w: str):
    return sorted(set(w))


# ============================================================
# Zyklische / antizyklische Vollwörter
# ============================================================

CYCLIC_PLUS = {"EABC", "ABCE", "BCEA", "CEAB"}
CYCLIC_MINUS = {"EACB", "ACBE", "CBEA", "BEAC"}

# Optional: weitere Rotationsfamilien könnten aufgenommen werden,
# falls du E nicht ausgezeichnet behandeln willst. Für den Moment
# bleibt E bewusst privilegiert.


# ============================================================
# Klassifikation im Typ (2,2)
# ============================================================

def classify_22_word(w: str):
    """
    Nur aufrufen, wenn Profil (2,2) ist.
    Liefert:
    - 'N' für Spiegeltyp XYYX / YXXY
    - 'A' für Alternation XYXY / YXYX
    - 'B' für Blocktyp XXYY / YYXX
    """
    if profile_of_word(w) != (2, 2):
        raise ValueError(f"{w} ist nicht vom Typ (2,2)")

    # beide Symbole extrahieren
    symbols = []
    seen = set()
    for ch in w:
        if ch not in seen:
            symbols.append(ch)
            seen.add(ch)
    if len(symbols) != 2:
        raise ValueError(f"{w} sollte genau 2 verschiedene Symbole enthalten")

    x, y = symbols[0], symbols[1]

    block = {x+x+y+y, y+y+x+x}
    alt = {x+y+x+y, y+x+y+x}
    mirror = {x+y+y+x, y+x+x+y}

    if w in mirror:
        return "N"
    elif w in alt:
        return "A"
    elif w in block:
        return "B"
    else:
        raise RuntimeError(f"Unerwartete (2,2)-Struktur für {w}")


# ============================================================
# Makroklassen
# ============================================================

def macro_class(w: str) -> str:
    prof = profile_of_word(w)

    if prof == (4,):
        return "U"  # uniform

    if prof == (3, 1):
        return "D"  # Defekt

    if prof == (2, 2):
        return classify_22_word(w)  # N, A, B

    if prof == (2, 1, 1):
        return "M"  # gemischt

    if prof == (1, 1, 1, 1):
        if w in CYCLIC_PLUS:
            return "Z+"
        elif w in CYCLIC_MINUS:
            return "Z-"
        else:
            return "Z0"  # vollheterogen, aber nicht in der ausgezeichneten Zyklusliste

    raise RuntimeError(f"Unbekanntes Profil für {w}")


# ============================================================
# Morley-Kompatibilität M(W)
# ============================================================

def morley_M(w: str) -> float:
    cls = macro_class(w)

    if cls == "N":
        return 0.0
    if cls == "B":
        return 0.30
    if cls == "A":
        return 0.70
    if cls == "M":
        return 0.55
    if cls == "D":
        return 0.25
    if cls == "U":
        return 0.10
    if cls in ("Z+", "Z-"):
        return 1.0
    if cls == "Z0":
        return 0.85

    return 0.0


# ============================================================
# Chiralität chi(W)
# ============================================================

def chirality_chi(w: str) -> int:
    cls = macro_class(w)

    if cls == "Z+":
        return +1
    if cls == "Z-":
        return -1

    # Alternation kann man ebenfalls orientiert lesen.
    # Hier einfache Konvention:
    if cls == "A":
        # orientiere anhand lexikographischer Kleinheit gegenüber Reverse
        return +1 if w < w[::-1] else -1

    return 0


# ============================================================
# Walter-Kernfunktion T(W)
# ============================================================

def walter_T(w: str) -> float:
    cls = macro_class(w)

    if cls == "N":
        return 0.80
    if cls == "B":
        return 0.50
    if cls == "A":
        return 0.35
    if cls == "M":
        return 0.45
    if cls == "D":
        return 0.20
    if cls == "U":
        return 0.15
    if cls in ("Z+", "Z-", "Z0"):
        return 0.20

    return 0.30


# ============================================================
# Äußere Aktivität Q(W)
# ============================================================

@dataclass
class QModelParams:
    lambda_M: float = 0.5
    lambda_chi: float = 0.3
    lambda_T: float = 0.4
    q_threshold_neutral: float = 0.18


def activity_Q(w: str, params: QModelParams = QModelParams()) -> float:
    M = morley_M(w)
    chi = abs(chirality_chi(w))
    T = walter_T(w)
    q = params.lambda_M * M + params.lambda_chi * chi - params.lambda_T * T
    return clip(q, 0.0, 1.0)


# ============================================================
# Äußere Signatur Sigma(W)
# ============================================================

def sigma_projection(w: str, params: QModelParams = QModelParams()) -> str:
    q = activity_Q(w, params)
    cls = macro_class(w)
    cnt = Counter(w)

    if q < params.q_threshold_neutral:
        return "0"

    # ausgezeichnete zyklische Vollwörter
    if cls == "Z+":
        mapping = {
            "EABC": "C",
            "ABCE": "E",
            "BCEA": "A",
            "CEAB": "B",
        }
        return mapping.get(w, "Zplus")

    if cls == "Z-":
        mapping = {
            "EACB": "B-",
            "ACBE": "E-",
            "CBEA": "A-",
            "BEAC": "C-",
        }
        return mapping.get(w, "Zminus")

    # alternierende Doppelwörter
    if cls == "A":
        syms = sorted(cnt.keys())
        return f"sigma_{syms[0]}{syms[1]}"

    # Blockwörter
    if cls == "B":
        syms = sorted(cnt.keys())
        return f"tau_{syms[0]}{syms[1]}"

    # Spiegeltyp
    if cls == "N":
        return "0"

    # Gemischt
    if cls == "M":
        # doppelt auftretendes Symbol hervorheben
        double_symbol = [k for k, v in cnt.items() if v == 2][0]
        singles = sorted([k for k, v in cnt.items() if v == 1])
        return f"mix_{double_symbol}_{''.join(singles)}"

    # Defekt
    if cls == "D":
        triple_symbol = [k for k, v in cnt.items() if v == 3][0]
        single_symbol = [k for k, v in cnt.items() if v == 1][0]
        return f"def_{triple_symbol}|{single_symbol}"

    # Uniform
    if cls == "U":
        return f"uni_{w[0]}"

    # vollheterogen, aber nicht ausgezeichnet
    if cls == "Z0":
        return "full_mix"

    return "unknown"


# ============================================================
# Wortanalyse
# ============================================================

def analyze_word(w: str, params: QModelParams = QModelParams()):
    return {
        "word": w,
        "profile": profile_of_word(w),
        "content": content_key(w),
        "class": macro_class(w),
        "M": morley_M(w),
        "chi": chirality_chi(w),
        "T": walter_T(w),
        "Q": round(activity_Q(w, params), 6),
        "Sigma": sigma_projection(w, params),
    }


def generate_all_words():
    return [''.join(p) for p in product(ALPHABET, repeat=4)]


# ============================================================
# Statistik über alle 256 Wörter
# ============================================================

def summarize_all_words(params: QModelParams = QModelParams()):
    words = generate_all_words()
    analyzed = [analyze_word(w, params) for w in words]

    by_class = defaultdict(list)
    by_profile = defaultdict(list)
    by_sigma = defaultdict(list)

    for row in analyzed:
        by_class[row["class"]].append(row)
        by_profile[row["profile"]].append(row)
        by_sigma[row["Sigma"]].append(row)

    return analyzed, by_class, by_profile, by_sigma


# ============================================================
# Energieabhängige Makroklassen-Dynamik
# ============================================================

def sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


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


MACRO_STATES = ["N", "A", "B", "M", "D", "U", "Z+", "Z-", "Z0", "END"]


def macro_features(cls: str):
    """
    Repräsentative mittlere M/chi/T/Q-Werte pro Makroklasse.
    """
    reps = {
        "N":  {"M": 0.00, "chi": 0,  "T": 0.80, "Q": 0.00},
        "A":  {"M": 0.70, "chi": 1,  "T": 0.35, "Q": 0.50},
        "B":  {"M": 0.30, "chi": 0,  "T": 0.50, "Q": 0.10},
        "M":  {"M": 0.55, "chi": 0,  "T": 0.45, "Q": 0.20},
        "D":  {"M": 0.25, "chi": 0,  "T": 0.20, "Q": 0.05},
        "U":  {"M": 0.10, "chi": 0,  "T": 0.15, "Q": 0.00},
        "Z+": {"M": 1.00, "chi": 1,  "T": 0.20, "Q": 0.72},
        "Z-": {"M": 1.00, "chi": -1, "T": 0.20, "Q": 0.72},
        "Z0": {"M": 0.85, "chi": 0,  "T": 0.20, "Q": 0.35},
        "END":{"M": 0.00, "chi": 0,  "T": 0.00, "Q": 1.00},
    }
    return reps[cls]


def build_macro_transition_matrix(sqrt_s: float,
                                  eparams: EnergyParams = EnergyParams()):
    """
    Reduzierte Übergangsmatrix auf Makroklassenebene.
    'END' ist absorbierend.
    """
    a = alpha_pinski(sqrt_s, eparams)
    b = beta_neutral(sqrt_s, eparams)
    g = gamma_detector(sqrt_s, eparams)

    n = len(MACRO_STATES)
    idx = {s: i for i, s in enumerate(MACRO_STATES)}
    P = np.zeros((n, n), dtype=float)

    for cls in MACRO_STATES:
        i = idx[cls]
        if cls == "END":
            P[i, i] = 1.0
            continue

        feat = macro_features(cls)
        M = feat["M"]
        T = feat["T"]
        Q = feat["Q"]

        row = np.zeros(n)

        # neutrale Schleifen
        row[idx["N"]] += b * (0.25 + 0.35 * T - 0.15 * M)
        row[idx["B"]] += b * (0.10 + 0.15 * T)
        row[idx["M"]] += b * (0.10 + 0.10 * (1.0 - T))

        # gerichtete / Pinski-Kanäle
        row[idx["A"]]  += a * (0.10 + 0.20 * M)
        row[idx["Z+"]] += a * (0.10 + 0.25 * max(0.0, feat["chi"]))
        row[idx["Z-"]] += a * (0.10 + 0.25 * max(0.0, -feat["chi"]))
        row[idx["Z0"]] += a * (0.10 + 0.10 * M)

        # Defekt/Uniform-Sektoren
        row[idx["D"]] += 0.05 + 0.10 * (1.0 - T)
        row[idx["U"]] += 0.02 + 0.05 * (1.0 - M)

        # direkte Absorption
        row[idx["END"]] += g * (0.08 + 0.35 * Q + 0.15 * M - 0.10 * T)

        # Selbstübergang schwach erlauben
        row[i] += 0.05 + 0.05 * T

        # Negativschutz und Normierung
        row = np.maximum(row, 0.0)
        s = row.sum()
        if s == 0:
            row[i] = 1.0
        else:
            row = row / s

        P[i] = row

    return P, idx


# ============================================================
# Hilfsroutinen für Ausgabe
# ============================================================

def print_class_summary(by_class):
    print("Makroklassen-Zusammenfassung")
    print("-" * 60)
    for cls in sorted(by_class.keys()):
        rows = by_class[cls]
        mean_M = sum(r["M"] for r in rows) / len(rows)
        mean_T = sum(r["T"] for r in rows) / len(rows)
        mean_Q = sum(r["Q"] for r in rows) / len(rows)
        print(f"{cls:>3s} : {len(rows):>3d} Wörter | "
              f"<M>={mean_M:.3f}, <T>={mean_T:.3f}, <Q>={mean_Q:.3f}")


def print_sample_words(by_class, n=8):
    print("\nBeispielwörter pro Klasse")
    print("-" * 60)
    for cls in sorted(by_class.keys()):
        words = [r["word"] for r in by_class[cls][:n]]
        print(f"{cls:>3s}: {', '.join(words)}")


def print_macro_matrix(P, idx):
    order = MACRO_STATES
    print("\nReduzierte Übergangsmatrix")
    print("-" * 60)
    header = "       " + "".join(f"{c:>8s}" for c in order)
    print(header)
    for c in order:
        i = idx[c]
        row = "".join(f"{P[i, idx[d]]:8.3f}" for d in order)
        print(f"{c:>6s} {row}")


# ============================================================
# Hauptprogramm
# ============================================================

if __name__ == "__main__":
    qparams = QModelParams()
    analyzed, by_class, by_profile, by_sigma = summarize_all_words(qparams)

    print(f"Anzahl aller Wörter: {len(analyzed)}")
    print_class_summary(by_class)
    print_sample_words(by_class)

    # Beispiele
    examples = ["ABBA", "BAAB", "ABAB", "AABB", "EABC", "ABCE", "EACB", "AABC", "AAAE", "EEEE"]
    print("\nEinzelanalysen")
    print("-" * 60)
    for w in examples:
        print(analyze_word(w, qparams))

    # Energieabhängige Makro-Dynamik
    for sqrt_s in [1.5, 3.5, 6.0]:
        print(f"\n=== sqrt(s) = {sqrt_s} ===")
        P, idx = build_macro_transition_matrix(sqrt_s)
        print_macro_matrix(P, idx)