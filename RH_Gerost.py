import json, os
import numpy as np
from math import isfinite

# ---------- utilities ----------
def ipr(vec: np.ndarray) -> float:
    p = np.abs(vec)**2
    p = p / (p.sum() + 1e-18)
    return float((p**2).sum())

def maxprob(vec: np.ndarray) -> float:
    p = np.abs(vec)**2
    p = p / (p.sum() + 1e-18)
    return float(p.max())

def top_nodes(vec: np.ndarray, k=10):
    p = np.abs(vec)**2
    p = p / (p.sum() + 1e-18)
    idx = np.argsort(-p)[:k]
    return idx.tolist(), p[idx].tolist()

def greedy_match(eigs, zetas, window=0.2):
    used = np.zeros(len(eigs), dtype=bool)
    matches = []
    for k, z in enumerate(zetas, start=1):
        diffs = np.abs(eigs - z)
        diffs[used] = np.inf
        j = int(np.argmin(diffs))
        if diffs[j] <= window:
            used[j] = True
            matches.append((k, float(z), int(j), float(eigs[j]), float(diffs[j])))
    return matches

# ---------- placeholders you must implement ----------
def build_operator(params):
    # TODO: build your BM operator / matrix here
    raise NotImplementedError

def compute_spectrum(A, num_modes):
    # TODO: eigensolver; return eigvals (sorted) and eigenvectors
    raise NotImplementedError

def get_zeta_targets(kind="first_20"):
    # TODO: load your stored gamma_k list; use your existing file (e.g. Rieman100000.tex or extracted list)
    # return numpy array of ordinates gamma_k
    raise NotImplementedError

# ---------- runner ----------
def run_experiment(params_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    with open(params_path, "r") as f:
        params = json.load(f)

    # seed
    seed = params["rng"]["seed"]
    np.random.seed(seed)

    # build operator
    A = build_operator(params)

    # spectrum
    num_modes = params["measurement"]["num_modes"]
    eigvals, eigvecs = compute_spectrum(A, num_modes)
    eigvals = np.array(eigvals, dtype=float)

    # zeta matching
    zetas = get_zeta_targets(params["measurement"].get("zeta_targets", "first_20"))
    window = params["measurement"].get("resonance_window", 0.2)
    matches = greedy_match(eigvals, zetas, window=window)

    # particle classification near matched modes
    ipr_thr = params["measurement"].get("particle_ipr_threshold", 0.15)
    loc_reports = []
    for (k, z, j, f, err) in matches:
        vec = eigvecs[:, j]
        I = ipr(vec)
        M = maxprob(vec)
        nodes, probs = top_nodes(vec, k=10)
        cls = "PARTICLE" if I <= ipr_thr else "EXTENDED"
        loc_reports.append({
            "zeta_k": k, "gamma": z, "mode_index": j, "eigenvalue": f, "abs_error": err,
            "ipr": I, "max_probability": M, "top_nodes": nodes, "top_probs": probs,
            "classification": cls
        })

    # write outputs
    with open(os.path.join(out_dir, "params.json"), "w") as f:
        json.dump(params, f, indent=2)

    # spectrum.csv
    csv_path = os.path.join(out_dir, "spectrum.csv")
    with open(csv_path, "w") as f:
        f.write("mode_index,eigenvalue\n")
        for i, v in enumerate(eigvals):
            if isfinite(v):
                f.write(f"{i},{v:.10f}\n")

    with open(os.path.join(out_dir, "localization.json"), "w") as f:
        json.dump(loc_reports, f, indent=2)

    return matches, loc_reports