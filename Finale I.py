from __future__ import annotations

import math

from dataclasses import dataclass

from pathlib import Path

from typing import Dict, List, Sequence, Tuple

import numpy as np

import pandas as pd

import matplotlib.pyplot as plt

# ------------------------------------------------------------

# Grundlegende Zahlentheorie

# ------------------------------------------------------------

def is_prime(n: int) -> bool:

    if n < 2:

        return False

    if n in (2, 3):

        return True

    if n % 2 == 0 or n % 3 == 0:

        return False

    r = int(math.isqrt(n))

    f = 5

    while f <= r:

        if n % f == 0 or n % (f + 2) == 0:

            return False

        f += 6

    return True

def primes_up_to(limit: int) -> List[int]:

    if limit < 2:

        return []

    sieve = np.ones(limit + 1, dtype=bool)

    sieve[:2] = False

    for p in range(2, int(limit**0.5) + 1):

        if sieve[p]:

            sieve[p * p : limit + 1 : p] = False

    return np.flatnonzero(sieve).tolist()

def prime_index_map(limit: int) -> Dict[int, int]:

    ps = primes_up_to(limit)

    return {p: i + 1 for i, p in enumerate(ps)}

def v_p(n: int, p: int) -> int:

    c = 0

    while n > 0 and n % p == 0:

        n //= p

        c += 1

    return c

def factorize(n: int) -> Dict[int, int]:

    if n < 1:

        raise ValueError("n must be >= 1")

    if n == 1:

        return {}

    rem = n

    fac: Dict[int, int] = {}

    for p in (2, 3):

        c = 0

        while rem % p == 0:

            rem //= p

            c += 1

        if c:

            fac[p] = c

    f = 5

    while f * f <= rem:

        for cand in (f, f + 2):

            c = 0

            while rem % cand == 0:

                rem //= cand

                c += 1

            if c:

                fac[cand] = c

        f += 6

    if rem > 1:

        fac[rem] = fac.get(rem, 0) + 1

    return fac

def radical(n: int) -> int:

    out = 1

    for p in factorize(n):

        out *= p

    return out

# ------------------------------------------------------------

# ERABC / Quaternionische Zustände

# ------------------------------------------------------------

def erabc_coords(n: int) -> Tuple[int, int, int, int]:

    fac = factorize(n)

    E = A = B = C = 0

    for p, e in fac.items():

        if p <= 3:

            continue

        r = p % 12

        if r == 1:

            E += e

        elif r == 5:

            A += e

        elif r == 7:

            B += e

        elif r == 11:

            C += e

    return E, A, B, C

def h_vec(n: int) -> np.ndarray:

    return np.array(erabc_coords(n), dtype=float)

def aggregate_state(values: Sequence[int]) -> np.ndarray:

    return np.sum([h_vec(n) for n in values], axis=0)

def direction_defect(values: Sequence[int]) -> float:

    H = aggregate_state(values)

    norm = np.linalg.norm(H)

    uhat = 0.5 * np.array([1.0, 1.0, 1.0, 1.0])

    return float(np.linalg.norm(H / norm - uhat))

def variance_of_values(values: Sequence[int]) -> float:

    arr = np.asarray(values, dtype=float)

    return float(np.mean((arr - arr.mean()) ** 2))

# ------------------------------------------------------------

# Mehrlinge

# ------------------------------------------------------------

@dataclass(frozen=True)

class Twin:

    p: int

    @property

    def values(self) -> Tuple[int, int]:

        return (self.p, self.p + 2)

@dataclass(frozen=True)

class TripletI:

    p: int

    @property

    def values(self) -> Tuple[int, int, int]:

        return (self.p, self.p + 2, self.p + 6)

@dataclass(frozen=True)

class TripletII:

    p: int

    @property

    def values(self) -> Tuple[int, int, int]:

        return (self.p, self.p + 4, self.p + 6)

@dataclass(frozen=True)

class Quadruplet:

    p: int

    @property

    def values(self) -> Tuple[int, int, int, int]:

        return (self.p, self.p + 2, self.p + 6, self.p + 8)

@dataclass(frozen=True)

class Quintuplet:

    p: int

    @property

    def values(self) -> Tuple[int, int, int, int, int]:

        return (self.p, self.p + 2, self.p + 6, self.p + 8, self.p + 12)

def first_twins(N: int, search_limit: int = 10_000_000) -> List[Twin]:

    out = []

    for p in range(3, search_limit):

        if is_prime(p) and is_prime(p + 2):

            out.append(Twin(p))

            if len(out) >= N:

                return out

    raise RuntimeError("Search limit too small for twins")

def first_triplets_I(N: int, search_limit: int = 20_000_000) -> List[TripletI]:

    out = []

    for p in range(3, search_limit):

        if is_prime(p) and is_prime(p + 2) and is_prime(p + 6):

            out.append(TripletI(p))

            if len(out) >= N:

                return out

    raise RuntimeError("Search limit too small for triplets-I")

def first_triplets_II(N: int, search_limit: int = 20_000_000) -> List[TripletII]:

    out = []

    for p in range(3, search_limit):

        if is_prime(p) and is_prime(p + 4) and is_prime(p + 6):

            out.append(TripletII(p))

            if len(out) >= N:

                return out

    raise RuntimeError("Search limit too small for triplets-II")

def first_quadruplets(N: int, search_limit: int = 20_000_000) -> List[Quadruplet]:

    out = []

    for p in range(5, search_limit):

        if is_prime(p) and is_prime(p + 2) and is_prime(p + 6) and is_prime(p + 8):

            out.append(Quadruplet(p))

            if len(out) >= N:

                return out

    raise RuntimeError("Search limit too small for quadruplets")

def first_quintuplets(N: int, search_limit: int = 200_000_000) -> List[Quintuplet]:

    out = []

    if N <= 0:

        return out

    if search_limit > 17 and all(is_prime(n) for n in (5, 7, 11, 13, 17)):

        out.append(Quintuplet(5))

        if len(out) >= N:

            return out

    for p in range(11, search_limit, 30):

        if (

            is_prime(p)

            and is_prime(p + 2)

            and is_prime(p + 6)

            and is_prime(p + 8)

            and is_prime(p + 12)

        ):

            out.append(Quintuplet(p))

            if len(out) >= N:

                return out

    raise RuntimeError("Search limit too small for quintuplets")

# ------------------------------------------------------------

# Einheitliche Features für alle Mehrlingstypen

# ------------------------------------------------------------

def interval_radical(values: Sequence[int]) -> int:

    a, b = min(values), max(values)

    prod = 1

    for n in range(a, b + 1):

        prod *= n

    return radical(prod)

def common_feature_vector(values: Sequence[int], idx_map: Dict[int, int]) -> np.ndarray:

    values = tuple(values)

    idxs = np.array([idx_map[n] for n in values], dtype=float)

    span = float(max(values) - min(values))

    geo_mean = float(np.prod(np.asarray(values, dtype=float)) ** (1.0 / len(values)))

    return np.array(

        [

            math.log(values[0]),                  # log-Start

            float(idxs.mean()),                  # mittlerer Primindex

            float(idxs[-1] - idxs[0]),           # Indexspanne

            math.log(interval_radical(values)),  # log Intervallradikal

            float(sum(v_p(n, 5) for n in values)), # 5-Tiefe

            direction_defect(values),            # Richtungsdefekt

            span / geo_mean,                     # relative Dehnung

            variance_of_values(values),          # Varianz

            float(np.dot(aggregate_state(values), aggregate_state(values))),  # Normquadrat

        ],

        dtype=float,

    )

# ------------------------------------------------------------

# Spektral-Pipeline

# ------------------------------------------------------------

def standardize_features(X: np.ndarray) -> np.ndarray:

    mu = X.mean(axis=0)

    sd = X.std(axis=0)

    sd = np.where(sd == 0.0, 1.0, sd)

    return (X - mu) / sd

def pairwise_euclidean(X: np.ndarray) -> np.ndarray:

    sq = np.sum(X**2, axis=1, keepdims=True)

    D2 = sq + sq.T - 2 * X @ X.T

    D2 = np.maximum(D2, 0.0)

    return np.sqrt(D2)

def gaussian_kernel_from_dist(D: np.ndarray, sigma: float) -> np.ndarray:

    W = np.exp(-(D**2) / (2.0 * sigma**2))

    np.fill_diagonal(W, 0.0)

    return W

def laplacian(W: np.ndarray) -> np.ndarray:

    D = np.diag(W.sum(axis=1))

    return D - W

def spectrum(M: np.ndarray) -> np.ndarray:

    vals = np.linalg.eigvalsh(M)

    vals.sort()

    return vals

def zscore(x: np.ndarray) -> np.ndarray:

    x = np.asarray(x, dtype=float)

    s = x.std()

    return (x - x.mean()) / (s if s != 0 else 1.0)

def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:

    x = np.asarray(x, dtype=float)

    y = np.asarray(y, dtype=float)

    x0 = x - x.mean()

    y0 = y - y.mean()

    denom = np.linalg.norm(x0) * np.linalg.norm(y0)

    if denom == 0:

        return float("nan")

    return float(np.dot(x0, y0) / denom)

def rmse(x: np.ndarray, y: np.ndarray) -> float:

    return float(np.sqrt(np.mean((x - y) ** 2)))

def mae(x: np.ndarray, y: np.ndarray) -> float:

    return float(np.mean(np.abs(x - y)))

def spectral_experiment(X: np.ndarray, sigma: float) -> np.ndarray:

    Xs = standardize_features(X)

    D = pairwise_euclidean(Xs)

    W = gaussian_kernel_from_dist(D, sigma)

    L = laplacian(W)

    return spectrum(L)

def sigma_sweep(X: np.ndarray, sigmas: Sequence[float], m: int, ref: np.ndarray) -> pd.DataFrame:

    rows = []

    for sigma in sigmas:

        spec = spectral_experiment(X, sigma)[:m]

        rows.append({

            "sigma": sigma,

            "corr": pearson_corr(spec, ref[:m]),

            "rmse": rmse(spec, ref[:m]),

            "mae": mae(spec, ref[:m]),

        })

    return pd.DataFrame(rows)

def feature_matrix(patterns: Sequence[object], idx_map: Dict[int, int]) -> np.ndarray:

    return np.vstack([common_feature_vector(pattern.values, idx_map) for pattern in patterns])

def summary_row(name: str, patterns: Sequence[object], sweep: pd.DataFrame) -> Dict[str, object]:

    best_idx = sweep["corr"].idxmax()

    return {

        "Raum": name,

        "Anzahl": len(patterns),

        "kleinstes Muster": str(patterns[0].values),

        "groesstes Muster": str(patterns[-1].values),

        "beste Korrelation": float(sweep["corr"].max()),

        "bestes sigma": float(sweep.loc[best_idx, "sigma"]),

    }

# ------------------------------------------------------------

# Gemeinsame Studie

# ------------------------------------------------------------

N = 80

sigmas = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0]

m = 20

ref_path = Path("zeros6.npy")

ref = zscore(np.load(ref_path).astype(float).ravel()[:N])

space_specs = [

    ("Zwillinge", "sigma_sweep_twins_zeta.csv", first_twins),

    ("Drillinge Typ I", "sigma_sweep_triplets_I_zeta.csv", first_triplets_I),

    ("Drillinge Typ II", "sigma_sweep_triplets_II_zeta.csv", first_triplets_II),

    ("Vierlinge", "sigma_sweep_quadruplets_zeta.csv", first_quadruplets),

    ("Fuenflinge", "sigma_sweep_quintuplets_zeta.csv", first_quintuplets),

]

spaces = {name: builder(N) for name, _, builder in space_specs}

max_val = max(

    patterns[-1].values[-1] for patterns in spaces.values()

) + 20

idx_map = prime_index_map(max_val)

sweeps = {

    name: sigma_sweep(feature_matrix(spaces[name], idx_map), sigmas, m, ref)

    for name, _, _ in space_specs

}

summary = pd.DataFrame(

    [summary_row(name, spaces[name], sweeps[name]) for name, _, _ in space_specs]

)

outdir = Path("bm_results_parallel")

outdir.mkdir(parents=True, exist_ok=True)

summary_path = outdir / "summary_parallel_zeta.csv"

plot_path = outdir / "parallel_comparison_zeta.png"

summary.to_csv(summary_path, index=False)

for name, filename, _ in space_specs:

    sweeps[name].to_csv(outdir / filename, index=False)

# Vergleichsplot

plt.figure(figsize=(8, 5))

for name, _, _ in space_specs:

    plt.plot(sweeps[name]["sigma"], sweeps[name]["corr"], marker="o", label=name)

plt.xlabel("sigma")

plt.ylabel("Korrelation mit Zeta-Referenz")

plt.title("Parallelvergleich Mehrlingsräume vs. zeros6.npy")

plt.legend()

plt.tight_layout()

plt.savefig(plot_path, dpi=200, bbox_inches="tight")

plt.show()

print("Dateien gespeichert in:", outdir)

print("\nZusammenfassung:")

print(summary.to_string(index=False))

print("\nSigma-Sweep Vierlinge:")

print(sweeps["Vierlinge"].to_string(index=False))

print("\nSigma-Sweep Fuenflinge:")

print(sweeps["Fuenflinge"].to_string(index=False))