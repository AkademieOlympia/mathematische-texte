import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter, defaultdict

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import LeaveOneOut

try:
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import squareform
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False

from Pre_ShorIII import (
    compare_numbers,
    TRAIN_LABELS,
    HurwitzPressureParams,
    QModelParams,
    EnergyParams,
    NOperatorParams,
)

# ---------------------------------------------------------------------
# Strukturkorrektor
# ---------------------------------------------------------------------

# Rollback-Schalter:
# True  = strukturelle Features (E-Quattropel / S2) nur diagnostisch verwenden
# False = alte experimentelle Vermischung mit Klassifikator erlauben
ROLLBACK_TO_STABLE_CLASSIFIER = True

# Falls man die Berichte weiterhin sehen möchte:
PRINT_STRUCTURAL_DIAGNOSTICS = True

# Diese Features sollen im Rollback NICHT in die Distanz-/PCA-/LOO-Klassifikation eingehen.
STRUCTURAL_ONLY_FEATURES = {
    "has_e_pair",
    "log_num_e_factor_pairs",
    "best_e_balance",
    "sqrt_e_distance",
    "has_quadruple_witness",
    "log_num_quadruple_witnesses",
    "quadruple_symmetry_score",
    "is_sum_of_two_squares",
    "log_num_sum_of_two_squares_repr",
    "log_num_s2_e_factor_pairs",
    "best_e_pair_s2_flag",
}

# Sentinel-Werte, die in Diagnoseblöcken erlaubt sind, aber nicht in die Klassifikation gehören
STRUCTURAL_ZERO_FILL = 0.0

STRUCTURE_CORRECTOR_ENABLED = True
STRUCTURE_MARGIN_THRESHOLD = 1.25
STRUCTURE_MAX_BONUS = 0.90
STRUCTURE_MAX_PENALTY = 0.55

STRUCTURE_TARGET_CLASSES = {
    "family-11", "family-13", "family-11-13",
    "neutral-prime-like", "mixed",
    "structured-composite",
    "family-3", "family-5", "family-7"
}

def _safe_get(dct, key, default=0.0):
    try:
        val = dct.get(key, default)
        if val is None:
            return default
        return float(val)
    except Exception:
        return default

def _clip(x, lo, hi):
    return max(lo, min(hi, x))

def safe_log1p(x):
    return math.log1p(max(0.0, x))

def sanitize_structural_feature_value(x):
    """
    Strukturelle Diagnosewerte robust machen:
    - keine 999-Sentinel mehr
    - keine negativen Ersatzwerte
    - alles auf >= 0 clampen
    """
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return STRUCTURAL_ZERO_FILL
    return max(0.0, float(x))

def structure_signature(features):
    sig = {}
    for k in [
        "has_e_pair", "log_num_e_factor_pairs", "best_e_balance", "sqrt_e_distance",
        "has_quadruple_witness", "log_num_quadruple_witnesses", "quadruple_symmetry_score",
        "is_sum_of_two_squares", "log_num_sum_of_two_squares_repr",
        "log_num_s2_e_factor_pairs", "best_e_pair_s2_flag"
    ]:
        sig[k] = _safe_get(features, k, 0.0)
    return sig

def structure_strength(sig):
    e_balance_good = math.exp(-2.0 * sig["best_e_balance"]) if sig["best_e_balance"] > 0 else 0.0
    sqrt_good = math.exp(-2.0 * sig["sqrt_e_distance"]) if sig["sqrt_e_distance"] > 0 else 0.0
    quad_good = sig["has_quadruple_witness"] * (1.0 + 0.35 * sig["log_num_quadruple_witnesses"])
    s2_good = sig["is_sum_of_two_squares"] * (0.5 + 0.25 * sig["best_e_pair_s2_flag"] + 0.15 * sig["log_num_s2_e_factor_pairs"])
    epair_good = sig["has_e_pair"] * (0.4 + 0.15 * sig["log_num_e_factor_pairs"])

    raw = (
        0.32 * epair_good
        + 0.33 * quad_good
        + 0.18 * s2_good
        + 0.10 * e_balance_good
        + 0.07 * sqrt_good
    )
    return _clip(raw, 0.0, 1.0)

def class_structure_adjustments(sig):
    strength = structure_strength(sig)
    has_quad = sig["has_quadruple_witness"] > 0.5
    has_e_pair = sig["has_e_pair"] > 0.5
    has_s2 = sig["is_sum_of_two_squares"] > 0.5

    adj = defaultdict(float)

    # Families
    if has_e_pair:
        adj["family-11"] -= 0.32 * strength
        adj["family-13"] -= 0.32 * strength
    if has_s2:
        adj["family-11"] -= 0.18 * strength
        adj["family-13"] -= 0.18 * strength

    # Family 11-13
    if has_quad:
        adj["family-11-13"] -= 0.55 * strength
    elif has_e_pair and has_s2:
        adj["family-11-13"] -= 0.12 * strength

    # Neutral Prime Like (Malus bei Struktur)
    if has_quad or has_e_pair:
        adj["neutral-prime-like"] += 0.42 * strength
    if has_s2 and has_e_pair:
        adj["neutral-prime-like"] += 0.16 * strength

    # Mixed
    if has_e_pair and not has_quad:
        adj["mixed"] -= 0.08 * strength
    if has_quad:
        adj["mixed"] += 0.12 * strength # Malus

    # Structured Composite (Coarse)
    if has_e_pair or has_quad or has_s2:
        adj["structured-composite"] -= 0.30 * strength
    if has_quad:
        adj["structured-composite"] -= 0.30 * strength

    for k in list(adj.keys()):
        adj[k] = _clip(adj[k], -STRUCTURE_MAX_BONUS, STRUCTURE_MAX_PENALTY)
    return dict(adj)

def apply_structure_corrector(base_distances, features):
    if ROLLBACK_TO_STABLE_CLASSIFIER:
        return dict(base_distances), {"used": False, "reason": "rollback"}

    if not STRUCTURE_CORRECTOR_ENABLED:
        return dict(base_distances), {"used": False}

    if len(base_distances) < 2:
        return dict(base_distances), {"used": False}

    items = sorted(base_distances.items(), key=lambda kv: kv[1])
    best_cls, best_dist = items[0]
    second_cls, second_dist = items[1]
    margin = second_dist - best_dist

    if margin > STRUCTURE_MARGIN_THRESHOLD:
        return dict(base_distances), {"used": False, "reason": "large_margin"}

    # Prüfen ob relevante Klassen involviert
    relevant = False
    for cls in (best_cls, second_cls):
        if cls in STRUCTURE_TARGET_CLASSES:
            relevant = True
            break
    if not relevant:
         return dict(base_distances), {"used": False, "reason": "irrelevant_classes"}

    sig = structure_signature(features)
    adjustments = class_structure_adjustments(sig)

    corrected = dict(base_distances)
    for cls, delta in adjustments.items():
        if cls in corrected:
            corrected[cls] += delta

    return corrected, {
        "used": True,
        "margin_before": margin,
        "adjustments": adjustments
    }

# --------------------------------------------------------------------
# Datensatz
# --------------------------------------------------------------------

DEFAULT_NUMBERS = [
    13, 17, 19, 23, 29, 31, 37, 41, 43,
    15, 21, 33, 39, 45, 51, 57, 63, 69, 75, 105,
    35, 55, 65, 85, 95,
    77, 91, 119, 133,
    121, 143, 169, 187, 209, 221
]

# --------------------------------------------------------------------
# Neue zweistufige Klassifikation
# --------------------------------------------------------------------

FINE_FAMILIES = [
    "family-3",
    "family-5",
    "family-7",
    "family-11",
    "family-13",
    "family-11-13",
]

COARSE_LABELS = [
    "neutral-prime-like",
    "mixed",
    "structured-composite",
]


def coarse_label_from_true_label(label: str) -> str:
    if label == "neutral-prime-like":
        return "neutral-prime-like"
    if label == "mixed":
        return "mixed"
    return "structured-composite"


def _feature_vector(row, features):
    return np.array([float(row.get(f, 0.0)) for f in features], dtype=float)


def _centroids_from_rows(rows, label_key, features, allowed_labels=None):
    buckets = defaultdict(list)
    for r in rows:
        lab = r[label_key]
        if allowed_labels is not None and lab not in allowed_labels:
            continue
        buckets[lab].append(_feature_vector(r, features))

    centroids = {}
    for lab, arrs in buckets.items():
        if len(arrs) == 0:
            continue
        centroids[lab] = np.mean(np.vstack(arrs), axis=0)
    return centroids


def _nearest_centroid_predict(vec, centroids):
    best_lab = None
    best_dist = float("inf")
    all_dists = {}
    for lab, c in centroids.items():
        d = float(np.linalg.norm(vec - c))
        all_dists[lab] = d
        if d < best_dist:
            best_dist = d
            best_lab = lab
    return best_lab, best_dist, all_dists


def add_two_stage_targets(rows):
    for r in rows:
        if "true" in r:
            r["coarse_label"] = coarse_label_from_true_label(r["true"])


def build_two_stage_feature_sets():
    # Stufe 1: grobe Trennung
    coarse_features = [
        # bestehende starke Signale
        "abs_A", "abs_B", "abs_C", "abs_E", "abs_0",
        "mean_steps",
        "eig_1", "eig_2", "eig_3", "eig_4",
        "eig_5", "eig_6", "eig_7", "eig_8",
        "qsd_A", "qsd_B", "qsd_N", "qsd_U",
        # E-Quattropel / S2 nur als Strukturmarker
        "has_e_pair",
        "log_num_e_factor_pairs",
        "best_e_balance",
        "sqrt_e_distance",
        "has_quadruple_witness",
        "log_num_quadruple_witnesses",
        "quadruple_symmetry_score",
        "is_sum_of_two_squares",
        "log_num_sum_of_two_squares_repr",
        "log_num_s2_e_factor_pairs",
        "best_e_pair_s2_flag",
    ]

    # Stufe 2: nur strukturierte Familien
    fine_features = [
        "abs_A", "abs_B", "abs_C", "abs_E", "abs_0",
        "mean_steps",
        "eig_1", "eig_2", "eig_3", "eig_4",
        "eig_5", "eig_6", "eig_7", "eig_8",
        "qsd_A", "qsd_B", "qsd_N", "qsd_U",
        "has_e_pair",
        "log_num_e_factor_pairs",
        "best_e_balance",
        "sqrt_e_distance",
        "has_quadruple_witness",
        "log_num_quadruple_witnesses",
        "quadruple_symmetry_score",
        "is_sum_of_two_squares",
        "log_num_sum_of_two_squares_repr",
        "log_num_s2_e_factor_pairs",
        "best_e_pair_s2_flag",
    ]

    if ROLLBACK_TO_STABLE_CLASSIFIER:
         coarse_features = [f for f in coarse_features if f not in STRUCTURAL_ONLY_FEATURES]
         fine_features = [f for f in fine_features if f not in STRUCTURAL_ONLY_FEATURES]

    return coarse_features, fine_features


def normalize_rows_inplace(rows, features):
    if not rows:
        return None, None
    X = np.vstack([_feature_vector(r, features) for r in rows])
    mu = X.mean(axis=0)
    sigma = X.std(axis=0)
    sigma[sigma == 0.0] = 1.0
    for r in rows:
        v = _feature_vector(r, features)
        vn = (v - mu) / sigma
        for i, f in enumerate(features):
            r[f"__norm__{f}"] = float(vn[i])
    return mu, sigma


def _norm_feature_vector(row, features):
    return np.array([float(row[f"__norm__{f}"]) for f in features], dtype=float)


def _centroids_from_normalized_rows(rows, label_key, features, allowed_labels=None):
    buckets = defaultdict(list)
    for r in rows:
        lab = r[label_key]
        if allowed_labels is not None and lab not in allowed_labels:
            continue
        buckets[lab].append(_norm_feature_vector(r, features))

    centroids = {}
    for lab, arrs in buckets.items():
        centroids[lab] = np.mean(np.vstack(arrs), axis=0)
    return centroids


def two_stage_predict(train_rows, test_row, coarse_features, fine_features):
    # getrennte Normalisierung pro Stufe
    # Wir kopieren, um inplace-Normalisierung nicht auf Originaldaten wirken zu lassen
    work_train = [dict(r) for r in train_rows]
    work_test = dict(test_row)

    normalize_rows_inplace(work_train + [work_test], coarse_features)
    coarse_centroids = _centroids_from_normalized_rows(
        work_train,
        "coarse_label",
        coarse_features,
        allowed_labels=COARSE_LABELS,
    )
    coarse_vec = _norm_feature_vector(work_test, coarse_features)
    coarse_pred, coarse_dist, coarse_dists = _nearest_centroid_predict(coarse_vec, coarse_centroids)

    # +++ STRUKTURKORREKTUR COARSE +++
    corrected_coarse_dists, debug_coarse = apply_structure_corrector(coarse_dists, test_row)
    if debug_coarse.get("used", False):
         best_coarse = min(corrected_coarse_dists, key=corrected_coarse_dists.get)
         coarse_pred = best_coarse
         coarse_dist = corrected_coarse_dists[best_coarse]
         coarse_dists = corrected_coarse_dists
    # +++ ENDE KORREKTUR +++

    if coarse_pred in ("neutral-prime-like", "mixed"):
        return {
            "pred": coarse_pred,
            "coarse_pred": coarse_pred,
            "coarse_dist": coarse_dist,
            "coarse_dists": coarse_dists,
            "fine_pred": None,
            "fine_dist": None,
            "fine_dists": {},
            "debug_coarse": debug_coarse,
            "debug_fine": {}
        }

    # Stufe 2: nur strukturierte Komposita
    fine_train = [r for r in work_train if r["coarse_label"] == "structured-composite"]
    
    if not fine_train:
         return {
            "pred": coarse_pred, # Fallback
            "coarse_pred": coarse_pred,
            "coarse_dist": coarse_dist,
            "coarse_dists": coarse_dists,
            "fine_pred": None,
            "fine_dist": None,
            "fine_dists": {},
            "debug_coarse": debug_coarse,
            "debug_fine": {}
        }

    normalize_rows_inplace(fine_train + [work_test], fine_features)
    fine_centroids = _centroids_from_normalized_rows(
        fine_train,
        "true",
        fine_features,
        allowed_labels=FINE_FAMILIES,
    )
    fine_vec = _norm_feature_vector(work_test, fine_features)
    fine_pred, fine_dist, fine_dists = _nearest_centroid_predict(fine_vec, fine_centroids)

    # +++ STRUKTURKORREKTUR FINE +++
    corrected_fine_dists, debug_fine = apply_structure_corrector(fine_dists, test_row)
    if debug_fine.get("used", False):
         best_fine = min(corrected_fine_dists, key=corrected_fine_dists.get)
         fine_pred = best_fine
         fine_dist = corrected_fine_dists[best_fine]
         fine_dists = corrected_fine_dists
    # +++ ENDE KORREKTUR +++

    return {
        "pred": fine_pred,
        "coarse_pred": coarse_pred,
        "coarse_dist": coarse_dist,
        "coarse_dists": coarse_dists,
        "fine_pred": fine_pred,
        "fine_dist": fine_dist,
        "fine_dists": fine_dists,
        "debug_coarse": debug_coarse,
        "debug_fine": debug_fine
    }


def run_two_stage_loo(rows):
    coarse_features, fine_features = build_two_stage_feature_sets()
    results = []
    # LeaveOneOut von sklearn
    loo = LeaveOneOut()
    
    # Konvertieren zu Liste für Indexing
    rows_list = list(rows)
    idxs = np.arange(len(rows_list))

    print(f"Starte 2-Stufen LOO auf {len(rows_list)} Datensätzen...")

    for train_idx, test_idx in loo.split(idxs):
        train_rows = [rows_list[i] for i in train_idx]
        test_row = rows_list[int(test_idx[0])]
        
        pred_info = two_stage_predict(train_rows, test_row, coarse_features, fine_features)
        
        results.append({
            "N": test_row["N"],
            "true": test_row["true"],
            "pred": pred_info["pred"],
            "dist": pred_info["fine_dist"] if pred_info["fine_dist"] is not None else pred_info["coarse_dist"],
            "coarse_pred": pred_info["coarse_pred"],
        })
    return results


def print_two_stage_report(results):
    print("\nZweistufiger LOO-Report")
    print("-" * 72)
    ok = 0
    for r in results:
        hit = (r["true"] == r["pred"])
        ok += int(hit)
        stat = "OK" if hit else "MISS"
        print(f"N={r['N']:4d}  true={r['true']:>18s}  pred={r['pred']:>18s}  dist={r['dist']:.4f}  {stat}")
    print()
    print(f"LOO-Genauigkeit: {ok}/{len(results)} = {ok/len(results):.4f}")


# --------------------------------------------------------------------
# Feature-Engineering
# --------------------------------------------------------------------

def enrich_row_with_base_features(row):
    # Dummy-Funktion, falls wir später noch was brauchen
    # Basisfeatures kommen schon aus flatten_profile_to_row
    return row


def add_e_quattropel_features(row):
    # Erwartet N in row["N"]
    N = int(row["N"])

    # --- bestehende Suchlogik beibehalten ---
    # Hilfsdefinition: E = 1 mod 12
    def is_E(x):
        return x > 0 and (x % 12 == 1)

    e_pairs = []
    lim = int(math.isqrt(N))
    for d in range(1, lim + 1):
        if N % d == 0:
            q = N // d
            if is_E(d) and is_E(q):
                e_pairs.append((d, q))

    raw_num_e_factor_pairs = len(e_pairs)
    raw_best_e_balance = None
    raw_sqrt_e_distance = None
    
    if e_pairs:
        rootN = math.sqrt(N)
        balances = []
        dists = []
        for x, y in e_pairs:
            lx = math.log(x)
            ly = math.log(y)
            balances.append(abs(lx - ly))
            dists.append(min(abs(x - rootN), abs(y - rootN)))

        raw_best_e_balance = min(balances)
        raw_sqrt_e_distance = min(dists) / max(1.0, rootN)

    raw_num_quadruple_witnesses = 0
    raw_quadruple_symmetry_score = None
    
    if len(e_pairs) >= 2:
        witnesses = 0
        sym_scores = []
        logs = [(math.log(x), math.log(y)) for x, y in e_pairs]
        for i in range(len(logs)):
            for j in range(i + 1, len(logs)):
                b1 = abs(logs[i][0] - logs[i][1])
                b2 = abs(logs[j][0] - logs[j][1])
                delta = abs(b1 - b2)
                if delta < 0.75:
                    witnesses += 1
                    sym_scores.append(delta)
        
        raw_num_quadruple_witnesses = witnesses
        if sym_scores:
            raw_quadruple_symmetry_score = min(sym_scores)

    # --- Mapping auf saubere Diagnosewerte ---
    row["has_e_pair"] = 1.0 if raw_num_e_factor_pairs > 0 else 0.0
    row["log_num_e_factor_pairs"] = safe_log1p(raw_num_e_factor_pairs)
    row["best_e_balance"] = sanitize_structural_feature_value(raw_best_e_balance)
    row["sqrt_e_distance"] = sanitize_structural_feature_value(raw_sqrt_e_distance)
    row["has_quadruple_witness"] = 1.0 if raw_num_quadruple_witnesses > 0 else 0.0
    row["log_num_quadruple_witnesses"] = safe_log1p(raw_num_quadruple_witnesses)
    row["quadruple_symmetry_score"] = sanitize_structural_feature_value(raw_quadruple_symmetry_score)

    return row


def add_s2_features(row):
    N = int(row["N"])

    raw_is_s2 = False
    raw_num_repr = 0
    
    reps = []
    lim = int(math.isqrt(N))
    for a in range(lim + 1):
        b2 = N - a * a
        if b2 < 0:
            break
        b = int(math.isqrt(b2))
        if a * a + b * b == N:
            if a <= b:
                reps.append((a, b))

    if reps:
        raw_is_s2 = True
        raw_num_repr = len(reps)

    # Kopplung an E-Paar-Struktur (Heuristik wie im Patch)
    raw_num_s2_e_pairs = 0
    raw_best_flag = False
    
    if row.get("has_e_pair", 0.0) > 0.5 and reps:
         # Hier approximieren wir Anzahl s2-E-Paare
         # Hinweis: Hier nutzen wir den bereits berechneten log-Wert, daher expm1
         s2_e_pairs_approx = int(round(math.expm1(row.get("log_num_e_factor_pairs", 0.0))))
         raw_num_s2_e_pairs = s2_e_pairs_approx
         raw_best_flag = True

    # --- Mapping auf saubere Diagnosewerte ---
    row["is_sum_of_two_squares"] = 1.0 if raw_is_s2 else 0.0
    row["log_num_sum_of_two_squares_repr"] = safe_log1p(raw_num_repr)
    row["log_num_s2_e_factor_pairs"] = safe_log1p(raw_num_s2_e_pairs)
    row["best_e_pair_s2_flag"] = 1.0 if raw_best_flag else 0.0

    return row


def enrich_row_with_structural_features(row):
    row = add_e_quattropel_features(row)
    row = add_s2_features(row)
    return row


def add_old_gate_features(row):
    # alte Gate-Features bleiben optional erhalten,
    # werden aber nicht mehr für die Hauptklassifikation benutzt
    return row


def flatten_profile_to_row(N, profile, label):
    row = {"N": N, "true": label}
    # Absorption
    for k, v in profile["absorption"].items():
        if k == "mean_steps":
            row["mean_steps"] = float(v)
        else:
            row[f"abs_{k}"] = float(v)
    # QSD
    for k, v in profile["qsd_macro"].items():
        row[f"qsd_{k}"] = float(v)
    # Eigenvalues
    for i, val in enumerate(profile["eigenvalues"]):
        row[f"eig_{i+1}"] = float(abs(val))
    return row


def build_dataset(profiles, labels_map):
    rows = []
    for p in profiles:
        N = p["N"]
        if N in labels_map:
            row = flatten_profile_to_row(N, p, labels_map[N])
            row = enrich_row_with_base_features(row)
            row = enrich_row_with_structural_features(row)
            row = add_old_gate_features(row)
            rows.append(row)
    add_two_stage_targets(rows)
    return rows


# --------------------------------------------------------------------
# Berichte
# --------------------------------------------------------------------

def print_structural_diagnostics(test_numbers):
    if not PRINT_STRUCTURAL_DIAGNOSTICS:
        return

    print("\nE-Quattropel-Kurztest")
    print("------------------------------------------------------------")
    # Wir berechnen Features ad-hoc, indem wir Dummys erstellen
    for N in test_numbers:
        row = {"N": N}
        row = add_e_quattropel_features(row)
        print(f"N={N}")
        for k in [
            "has_e_pair",
            "log_num_e_factor_pairs",
            "best_e_balance",
            "sqrt_e_distance",
            "has_quadruple_witness",
            "log_num_quadruple_witnesses",
            "quadruple_symmetry_score",
        ]:
            if k in row:
                print(f"  {k}: {row[k]}")
        print()

    print("\nS2-Kurztest")
    print("------------------------------------------------------------")
    for N in test_numbers:
        row = {"N": N}
        # e_features vorher nötig für s2? Ja für Koppelung
        row = add_e_quattropel_features(row)
        row = add_s2_features(row)
        print(f"N={N}")
        for k in [
            "is_sum_of_two_squares",
            "log_num_sum_of_two_squares_repr",
            "log_num_s2_e_factor_pairs",
            "best_e_pair_s2_flag",
        ]:
            if k in row:
                print(f"  {k}: {row[k]}")
        print()


def print_pca_report(df, feature_cols):
    X = df[feature_cols].values
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)

    pca = PCA(n_components=2)
    Y = pca.fit_transform(Xs)
    evr = pca.explained_variance_ratio_

    mode_label = "stable-core" if ROLLBACK_TO_STABLE_CLASSIFIER else "full-experimental"
    print(f"\nPCA Report (Features: {len(feature_cols)}, mode={mode_label})")
    print(f"Explained Variance: PC1={evr[0]:.4f}, PC2={evr[1]:.4f}")
    
    plt.figure(figsize=(10, 7))
    plt.scatter(Y[:, 0], Y[:, 1])
    
    Ns = df["N"].values
    for i, n in enumerate(Ns):
        plt.text(Y[i, 0], Y[i, 1], str(n), fontsize=8)
        
    plt.xlabel(f"PC1 ({evr[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({evr[1]*100:.1f}%)")
    plt.title(f"Zweistufige Merkmale PCA ({mode_label})")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def print_confusion_for_two_stage(results):
    labels = [
        "family-11",
        "family-11-13",
        "family-13",
        "family-3",
        "family-5",
        "family-7",
        "mixed",
        "neutral-prime-like",
    ]
    # Filter labels to only those present if needed, or keep fixed list
    # Fixed list is better for comparison
    
    y_true = [r["true"] for r in results]
    y_pred = [r["pred"] for r in results]
    
    # Check if labels appear in data
    unique_labels = sorted(list(set(y_true) | set(y_pred)))
    # Ensure all our interest labels are there if possible
    
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    print("\nKonfusionsmatrix")
    print("-" * 72)
    print("true\\pred".ljust(25) + "".join(f"{lab:>18s}" for lab in labels))
    for i, lab in enumerate(labels):
        line = f"{lab:<25s}"
        for j in range(len(labels)):
            line += f"{cm[i, j]:>18d}"
        print(line)


# --------------------------------------------------------------------
# Main
# --------------------------------------------------------------------

def main():
    hparams = HurwitzPressureParams(H_inf=20.0)
    
    # Wir nehmen die Standardzahlen aus Pre_ShorIII.TRAIN_LABELS + Testzahlen
    # Aber Pre_ShorIII hat kein "numbers" list exportiert.
    # Wir nehmen DEFAULT_NUMBERS oben definiert.
    
    # Berechnung der Profile
    print("Berechne Profile...")
    profiles = compare_numbers(
        DEFAULT_NUMBERS,
        sqrt_s=3.5,
        hparams=hparams
    )
    
    # Dataset bauen
    rows = build_dataset(profiles, TRAIN_LABELS)
    df = pd.DataFrame(rows)
    
    # PCA Report
    pca_feature_cols = [
        c for c in df.columns
        if c not in ("N", "true", "coarse_label")
        and not c.startswith("__norm__")
    ]
    print_pca_report(df, pca_feature_cols)
    
    # LOO
    two_stage_results = run_two_stage_loo(rows)
    print_two_stage_report(two_stage_results)
    print_confusion_for_two_stage(two_stage_results)
    
    # Testzahlen Vorhersage
    print("\nBeispiel-Testzahlen")
    print("-" * 72)
    test_candidates = [27, 45, 63, 75, 105, 121, 143, 169, 187, 209, 221]
    
    # Wir müssen Profile für diese Zahlen berechnen
    test_profiles = compare_numbers(
        test_candidates,
        sqrt_s=3.5,
        hparams=hparams
    )
    
    # Features für Testzahlen
    coarse_feats, fine_feats = build_two_stage_feature_sets()
    
    for p in test_profiles:
        N = p["N"]
        # Dummy Label "unknown" für Test
        row = flatten_profile_to_row(N, p, "unknown")
        row = enrich_row_with_structural_features(row)
        
        # Train ist alle rows
        pred = two_stage_predict(rows, row, coarse_feats, fine_feats)
        dist = pred["fine_dist"] if pred["fine_dist"] is not None else pred["coarse_dist"]
        
        print(f"Testzahl {N}: pred={pred['pred']:>18s}, dist={dist:.4f}")

        # Debug-Output für Strukturkorrektur
        if pred.get("debug_coarse", {}).get("used"):
             d = pred["debug_coarse"]
             print(f"   [Coarse-Correction] margin={d.get('margin_before',0):.4f}")
             for k,v in d.get("adjustments", {}).items():
                 if abs(v) > 0.001: print(f"     {k}: {v:+.4f}")

        if pred.get("debug_fine", {}).get("used"):
             d = pred["debug_fine"]
             print(f"   [Fine-Correction]   margin={d.get('margin_before',0):.4f}")
             for k,v in d.get("adjustments", {}).items():
                 if abs(v) > 0.001: print(f"     {k}: {v:+.4f}")

    print_structural_diagnostics([121, 143, 169, 385, 1001, 105])


if __name__ == "__main__":
    main()
