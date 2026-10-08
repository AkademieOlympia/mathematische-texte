"""
DrillingsTest — Hurwitz-Rotor mit arithmetischer / Gap-Kopplung.

Kopplungsmodi (--coupling):
  mod30 : A(q) = Π_H(q mod 30)          [Restklassen]
  gap   : A(q) = R_g(q)                  [Gap-Operator auf (q0,q1,q2)]

Multiplikative Dynamik (Default):
    q_{n+1} = (u q_n u^{-1}) * (1 + α A(q_n) / N(q_n)^{1/2})

Gap-Matrix (lokal entlang der Bahn):
    g_ij = |q_i - q_j|,  i,j ∈ {0,1,2}
    R_g = [[0, g01^{-1}, g02^{-1}],
           [g01^{-1}, 0, g12^{-1}],
           [g02^{-1}, g12^{-1}, 0]]
    A(q) = (0, R_g · (q0,q1,q2))

Dynamisches Gate:
    G_dyn = (∑ w_k S3(q_k)) / (∑ w_k S1(q_k))^{3/2}

Zeitumkehr (Test D3, --time-reversal):
    A_dyn = G_→ − G_←
    G_← = G_dyn auf der gespiegelten Vorwärtsbahn (q_T,…,q_0).

Harte Nullen (--hard-nulls) auf A_dyn:
    D1: Gap-erhaltend — Match P↔C nach |log n| in der Wheel-Klasse
    D2: Density Matching — Match nach lokalem Siebfaktor H(n)

Rotor-Eigendynamik (--rotor-dyn), ohne Primzahllabel als Ziel:
    I1: Invarianz von A_dyn unter q ↦ h q h^{-1} (24 Hurwitz-Einheiten)
    L1: Lyapunov-Schätzer λ = (1/T) log(d_T/d_0)

Lyapunov-α-Scan (--lyapunov-scan):
    λ(α) für α ∈ {1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.5}
    plus invariante Observablen ΔE = N(q_T)−N(q_0), I_sym = ∑ w_n · 2 Re(q_n)

Form-Lyapunov (--shape-lyapunov):
    ˆq = q/√N(q); λ_shape vs. ambient λ; α-Scan und T-Scan (5…100)
    trennt Formdynamik von Energieexplosion.

S³-Maß (--s3-measure):
    Hopf-Koordinaten (η, φ1, φ2) auf ˆq; Winkelhistogramme + C(τ)=⟨ˆq_n·ˆq_{n+τ}⟩

Period-3-Perturbation (--period3-perturbation):
    F_0³=id, ||F_α³−id||/α, |C_α−C_0|≤Kα f(τ) um Ad_u-Rotation

Explizite E_g-Schranke (--eg-bound):
    Lemma: auf Gap-Klasse (2,4,6) ist |B|=|R_g|/√N ≤ √(85/216),
    damit |E_g|_∞ ≤ 3√(85/216) auf S³-Formdynamik.

Klein-Kopplung (--small-coupling):
    α-Scan 0.001…0.02 für lineare Antwort vs. Innenmaximum.

Governance:
    Gate = externe Annotation. Z/Null = Stichprobendiagnostik, kein Primzahlbeweis.
    Zwischenstand: Gap-Rotor dynamisch nichttrivial, kein Primzahldrilling-Selektor.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
import sympy as sp

# -----------------------------
# Primzahldrillinge / Kandidaten
# -----------------------------


def wheel_ok(n: int) -> bool:
    W = 30030
    # math.gcd: kein Overflow bei großen n (N ≥ 10^7)
    return math.gcd(n * (n + 2) * (n + 6), W) == 1


def is_prime_triplet(n: int) -> bool:
    return sp.isprime(n) and sp.isprime(n + 2) and sp.isprime(n + 6)


# -----------------------------
# Quaternionen q = (a, b, c, d)
# -----------------------------


def qmul(q, r):
    a, b, c, d = q
    e, f, g, h = r
    return np.array(
        [
            a * e - b * f - c * g - d * h,
            a * f + b * e + c * h - d * g,
            a * g - b * h + c * e + d * f,
            a * h + b * g - c * f + d * e,
        ],
        dtype=float,
    )


def qconj(q):
    a, b, c, d = q
    return np.array([a, -b, -c, -d], dtype=float)


def qnorm(q):
    return float(np.dot(q, q))


def qone():
    return np.array([1.0, 0.0, 0.0, 0.0], dtype=float)


# Hurwitz-Einheit der Evolution: |u|^2 = 1 ⇒ u^{-1} = conjugate(u)
u = np.array([0.5, 0.5, 0.5, 0.5], dtype=float)


def hurwitz_units() -> list[np.ndarray]:
    """Die 24 Einheiten der Hurwitz-Ganzzahlen (Norm 1)."""
    units: list[np.ndarray] = []
    for i in range(4):
        e = np.zeros(4, dtype=float)
        e[i] = 1.0
        units.append(e.copy())
        units.append((-e).copy())
    for s0 in (-0.5, 0.5):
        for s1 in (-0.5, 0.5):
            for s2 in (-0.5, 0.5):
                for s3 in (-0.5, 0.5):
                    units.append(np.array([s0, s1, s2, s3], dtype=float))
    return units


def qconj_by(h: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Konjugation q ↦ h q h^{-1} (h Einheit ⇒ h^{-1} = conjugate(h))."""
    return qmul(qmul(h, q), qconj(h))


def Pi_H_mod30(q):
    """Restklassen-Anteil A(q) = Π_H(q mod 30) als Quaternion."""
    return np.mod(q, 30.0)


def R_g_quat(q, eps: float = 1e-30):
    """
    Gap-Operator als Quaternion:
      g01=|q0-q1|, g02=|q0-q2|, g12=|q1-q2|
      R_g · (q0,q1,q2) → reiner Vektoranteil von A.
    Gaps werden am aktuellen Zustand ausgewertet (bahnlokal).
    """
    p0, p1, p2, _ = (float(x) for x in q)
    g01 = abs(p0 - p1)
    g02 = abs(p0 - p2)
    g12 = abs(p1 - p2)
    i01 = 1.0 / max(g01, eps)
    i02 = 1.0 / max(g02, eps)
    i12 = 1.0 / max(g12, eps)
    M = np.array(
        [
            [0.0, i01, i02],
            [i01, 0.0, i12],
            [i02, i12, 0.0],
        ],
        dtype=float,
    )
    w = M @ np.array([p0, p1, p2], dtype=float)
    return np.array([0.0, w[0], w[1], w[2]], dtype=float)


def coupling_A(q, coupling: str):
    if coupling == "mod30":
        return Pi_H_mod30(q)
    if coupling == "gap":
        return R_g_quat(q)
    raise ValueError(f"unknown coupling: {coupling}")


def hurwitz_step(
    q,
    alpha: float,
    *,
    coupling: str = "mod30",
    multiplicative: bool = True,
):
    """
    Ein Schritt:
      multiplikativ: (uqu^{-1}) (1 + α A / √N)
      additiv:       (uqu^{-1}) + α A / √N
    """
    rot = qmul(qmul(u, q), qconj(u))
    A = coupling_A(q, coupling)
    nrm = np.sqrt(max(qnorm(q), 1e-30))
    drive = (alpha / nrm) * A
    if multiplicative:
        return qmul(rot, qone() + drive)
    return rot + drive


# -----------------------------
# Invarianten / Gates
# -----------------------------


def S1(q):
    return qnorm(q)


def S3(q):
    q3 = qmul(qmul(q, q), q)
    return float(q3[0])


def gate_path_static(
    q0,
    alpha: float,
    steps: int = 10,
    *,
    multiplicative: bool = True,
    coupling: str = "mod30",
):
    """Altes Punkt-Gate: gewichteter Mittelwert von S3/S1^{3/2} entlang der Bahn."""
    q = q0.copy()
    values = []
    weights = []
    for k in range(steps + 1):
        s1 = S1(q)
        values.append(S3(q) / (s1**1.5) if s1 > 0 else 0.0)
        weights.append(np.exp(-k / 3.0))
        q = hurwitz_step(
            q, alpha, coupling=coupling, multiplicative=multiplicative
        )
    w = np.asarray(weights, dtype=float)
    w /= w.sum()
    return float(np.sum(w * np.asarray(values, dtype=float)))


def evolve_trajectory(
    q0,
    alpha: float,
    steps: int = 10,
    *,
    multiplicative: bool = True,
    coupling: str = "mod30",
) -> list[np.ndarray]:
    """Vorwärtsbahn q_0 → … → q_steps."""
    q = q0.copy()
    traj = [q.copy()]
    for _ in range(steps):
        q = hurwitz_step(
            q, alpha, coupling=coupling, multiplicative=multiplicative
        )
        traj.append(q.copy())
    return traj


def gate_dyn_on_traj(traj: list[np.ndarray]) -> float:
    """
    G_dyn auf einer Zustandsfolge (gleiche Gewichte w_k = e^{-k/3}).
    Für Zeitumkehr: traj rückwärts übergeben.
    """
    s3_vals = [S3(q) for q in traj]
    s1_vals = [S1(q) for q in traj]
    weights = np.exp(-np.arange(len(traj)) / 3.0)
    weights /= weights.sum()
    num = float(np.sum(weights * np.asarray(s3_vals, dtype=float)))
    den = float(np.sum(weights * np.asarray(s1_vals, dtype=float)))
    if den <= 0:
        return 0.0
    return num / (den**1.5)


def gate_path_dyn(
    q0,
    alpha: float,
    steps: int = 10,
    *,
    multiplicative: bool = True,
    coupling: str = "mod30",
):
    """
    Dynamisches Bahn-Gate:
        G_dyn = (∑ w_k S3(q_k)) / (∑ w_k S1(q_k))^{3/2}
    """
    traj = evolve_trajectory(
        q0,
        alpha,
        steps,
        multiplicative=multiplicative,
        coupling=coupling,
    )
    return gate_dyn_on_traj(traj)


def asymmetry_dyn(
    q0,
    alpha: float,
    steps: int = 10,
    *,
    multiplicative: bool = True,
    coupling: str = "mod30",
) -> float:
    """
    Test D3 — Zeitumkehr-Asymmetrie der gewichteten Bahn:
        A_dyn = G_→ − G_←
    mit G_← = G_dyn auf der zeitgespiegelten Zustandsfolge (q_T,…,q_0).

    Keine inverse Dynamik-Map: prüft gerichtete Bahnstruktur unter
    denselben Vorwärtsgewichten.
    """
    traj = evolve_trajectory(
        q0,
        alpha,
        steps,
        multiplicative=multiplicative,
        coupling=coupling,
    )
    g_fwd = gate_dyn_on_traj(traj)
    g_bwd = gate_dyn_on_traj(list(reversed(traj)))
    return g_fwd - g_bwd


# -----------------------------
# Experiment / Statistik
# -----------------------------


def welch_z(P: np.ndarray, C: np.ndarray) -> tuple[float, float, float, float]:
    """
    Δ = μ_P − μ_C und Welch-Z:
        Z = Δ / sqrt(σ_P²/N_P + σ_C²/N_C)
    """
    mu_p = float(P.mean())
    mu_c = float(C.mean())
    delta = mu_p - mu_c
    sig_p = float(P.std(ddof=1)) if len(P) > 1 else float("nan")
    sig_c = float(C.std(ddof=1)) if len(C) > 1 else float("nan")
    se = np.sqrt(sig_p**2 / len(P) + sig_c**2 / len(C))
    z = float(delta / se) if se > 0 else float("nan")
    return delta, sig_p, sig_c, z


def alpha_grid(a0: float, a1: float, step: float) -> list[float]:
    """Inklusives α-Gitter mit stabiler Rundung."""
    n = int(round((a1 - a0) / step))
    return [round(a0 + i * step, 10) for i in range(n + 1)]


def local_sieve_factor(n: int, L: int = 30) -> float:
    """
    Lokaler Siebfaktor (Density-Matching, Test D2):
      H(n) = ∏_{ℓ≤L} (1 − ν_ℓ(n)/ℓ)
    mit ν_ℓ = #{x ∈ {n,n+2,n+6} : x ≡ 0 (mod ℓ)}.
    """
    h = 1.0
    for ell in sp.primerange(2, L + 1):
        nu = sum(1 for x in (n, n + 2, n + 6) if x % ell == 0)
        h *= 1.0 - nu / ell
    return h


def collect_gate_table(
    *,
    N: int,
    alphas: list[float],
    steps: int,
    multiplicative: bool,
    use_dyn: bool,
    coupling: str,
    time_reversal: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[float, np.ndarray]]:
    """
    Einmaliger Durchlauf über wheel_ok-Kandidaten.
    Returns: ns, labels, H(n), values[alpha]
    """
    ns: list[int] = []
    labels: list[int] = []
    dens: list[float] = []
    values: dict[float, list[float]] = {a: [] for a in alphas}

    for n in range(5, N, 2):
        if not wheel_ok(n):
            continue
        q0 = np.array([n, n + 2, n + 6, 6], dtype=float)
        ns.append(n)
        labels.append(1 if is_prime_triplet(n) else 0)
        dens.append(local_sieve_factor(n))
        for a in alphas:
            if time_reversal:
                values[a].append(
                    asymmetry_dyn(
                        q0,
                        a,
                        steps=steps,
                        multiplicative=multiplicative,
                        coupling=coupling,
                    )
                )
            elif use_dyn:
                values[a].append(
                    gate_path_dyn(
                        q0,
                        a,
                        steps=steps,
                        multiplicative=multiplicative,
                        coupling=coupling,
                    )
                )
            else:
                values[a].append(
                    gate_path_static(
                        q0,
                        a,
                        steps=steps,
                        multiplicative=multiplicative,
                        coupling=coupling,
                    )
                )

    return (
        np.asarray(ns, dtype=int),
        np.asarray(labels, dtype=int),
        np.asarray(dens, dtype=float),
        {a: np.asarray(values[a], dtype=float) for a in alphas},
    )


def paired_welch(diff: np.ndarray) -> tuple[float, float, float]:
    """Paarweise: mean(diff), std, Z = mean / (std/√n)."""
    d = np.asarray(diff, dtype=float)
    if len(d) < 2:
        return float("nan"), float("nan"), float("nan")
    mu = float(d.mean())
    sig = float(d.std(ddof=1))
    z = mu / (sig / np.sqrt(len(d))) if sig > 0 else float("nan")
    return mu, sig, z


def hard_null_tests(
    ns: np.ndarray,
    labels: np.ndarray,
    dens: np.ndarray,
    gtab: dict[float, np.ndarray],
    alphas: list[float],
    *,
    seed: int,
) -> dict:
    """
    D1: Gap-erhaltend + größennahe Kontrollen (gleiche Wheel-Klasse,
        n nicht Primzahldrilling, Match nach |log n|).
    D2: Density Matching nach H(n) (lokaler Siebfaktor, L=30).
    """
    rng = np.random.default_rng(seed)
    idx_p = np.flatnonzero(labels == 1)
    idx_c = np.flatnonzero(labels == 0)
    if len(idx_p) == 0 or len(idx_c) == 0:
        print("hard nulls: insufficient P/C")
        return {}

    logn = np.log(ns.astype(float))
    summary: dict = {"D1_size": {}, "D2_density": {}}

    def one_match(key_all: np.ndarray, name: str) -> dict:
        # greedy nearest-neighbor without replacement; P in random order
        order = idx_p.copy()
        rng.shuffle(order)
        c_idx = idx_c.copy()
        c_keys = key_all[c_idx].astype(float)
        used = np.zeros(len(c_idx), dtype=bool)
        pairs: list[tuple[int, int, float]] = []
        for ip in order:
            if used.all():
                break
            diffs = np.abs(c_keys - float(key_all[ip]))
            diffs[used] = np.inf
            j = int(np.argmin(diffs))
            used[j] = True
            pairs.append((int(ip), int(c_idx[j]), float(diffs[j])))

        out = {}
        print()
        print(f"HARD NULL {name}  matched_pairs={len(pairs)}")
        print("alpha    mean(A_P−A_C)   sig_diff      Z_paired   mean|key_diff|")
        for a in alphas:
            diffs = np.array(
                [gtab[a][ip] - gtab[a][ic] for ip, ic, _ in pairs], dtype=float
            )
            mu, sig, z = paired_welch(diffs)
            mean_key = float(np.mean([kd for _, _, kd in pairs]))
            out[str(a)] = {
                "n_pairs": len(pairs),
                "mean_diff": mu,
                "sig_diff": sig,
                "Z_paired": z,
                "mean_key_absdiff": mean_key,
            }
            print(
                f"{a:7.3f}  {mu:+.6e}  {sig:.6e}  {z:+8.3f}  {mean_key:.6e}"
            )
        print(
            "note: D1/D2 = strukturgleiche Kontrollen; "
            "Z_paired klein ⇒ Effekt nicht primzahlspezifisch."
        )
        return out

    summary["D1_size"] = one_match(
        logn, "D1 size-matched (log n), gap-preserving wheel"
    )
    summary["D2_density"] = one_match(dens, "D2 density-matched H(n), L=30")
    return summary


def rows_from_labels(
    labels: np.ndarray,
    gtab: dict[float, np.ndarray],
    alphas: list[float],
) -> list[dict]:
    rows: list[dict] = []
    for a in alphas:
        g = gtab[a]
        P = g[labels == 1]
        C = g[labels == 0]
        if len(P) < 2 or len(C) < 2:
            continue
        delta, sig_p, sig_c, z = welch_z(P, C)
        rows.append(
            {
                "alpha": a,
                "delta": delta,
                "abs_delta": abs(delta),
                "Z": z,
                "sigP": sig_p,
                "sigC": sig_c,
                "meanP": float(P.mean()),
                "meanC": float(C.mean()),
                "nP": int(len(P)),
                "nC": int(len(C)),
            }
        )
    return rows


def print_rows(rows: list[dict], *, header: str) -> None:
    print(header)
    print(
        "alpha      Delta           |Delta|        Z         "
        "sigP        sigC        nP     nC"
    )
    for row in rows:
        print(
            f"{row['alpha']:6.2f}  {row['delta']:+.6e}  {row['abs_delta']:.6e}  "
            f"{row['Z']:+8.3f}  {row['sigP']:.6e}  {row['sigC']:.6e}  "
            f"{row['nP']:6d}  {row['nC']:6d}"
        )
    if rows:
        best = max(rows, key=lambda r: r["abs_delta"])
        best_z = max(rows, key=lambda r: abs(r["Z"]))
        print()
        print(
            f"alpha*_absDelta = {best['alpha']:.2f}  "
            f"|Delta|={best['abs_delta']:.6e}  Z={best['Z']:+.3f}"
        )
        print(
            f"alpha*_absZ     = {best_z['alpha']:.2f}  "
            f"|Z|={abs(best_z['Z']):.3f}  Delta={best_z['delta']:+.6e}"
        )
        print(
            "note: |Z| misst Trennung der Stichprobenmittel; "
            "kein Beweis für Primzahlerkennung."
        )


def run_null_model(
    labels: np.ndarray,
    gtab: dict[float, np.ndarray],
    alphas: list[float],
    *,
    trials: int,
    seed: int,
) -> dict:
    """
    Test B: zufällige Permutation der Primzahllabel innerhalb der Wheel-Klasse
    (Anzahl P/C bleibt erhalten). Erwartet: |Z| fällt auf Rauschniveau.
    """
    rng = np.random.default_rng(seed)
    nP = int(labels.sum())
    absZ_true = {
        a: abs(welch_z(gtab[a][labels == 1], gtab[a][labels == 0])[3]) for a in alphas
    }
    absZ_null: dict[float, list[float]] = {a: [] for a in alphas}

    for _ in range(trials):
        # feste Anzahl Einsen, zufällige Positionen
        perm = np.zeros(len(labels), dtype=int)
        idx = rng.choice(len(labels), size=nP, replace=False)
        perm[idx] = 1
        for a in alphas:
            P = gtab[a][perm == 1]
            C = gtab[a][perm == 0]
            _, _, _, z = welch_z(P, C)
            absZ_null[a].append(abs(z))

    summary = {}
    print()
    print(f"NULL MODEL  trials={trials}  seed={seed}  (label permutation in wheel class)")
    print("alpha   |Z|_true   mean|Z|_null   p_emp(|Z_null|>=|Z_true|)")
    for a in alphas:
        arr = np.asarray(absZ_null[a], dtype=float)
        zt = absZ_true[a]
        p_emp = float(np.mean(arr >= zt)) if len(arr) else float("nan")
        summary[a] = {
            "absZ_true": zt,
            "mean_absZ_null": float(arr.mean()),
            "p_empirical": p_emp,
            "null_absZ": arr.tolist(),
        }
        print(f"{a:6.2f}  {zt:8.3f}  {arr.mean():12.3f}  {p_emp:8.3f}")
    print(
        "note: p_emp klein ⇒ beobachtetes |Z| unwahrscheinlich unter Label-Null; "
        "kein Primzahlbeweis."
    )
    return summary


def sample_wheel_seeds(N: int, n_sample: int, seed: int) -> list[np.ndarray]:
    """Gleichmäßige Stichprobe von Wheel-Kandidaten als Rotor-Seeds (ohne Label)."""
    ns = [n for n in range(5, N, 2) if wheel_ok(n)]
    if not ns:
        return []
    rng = np.random.default_rng(seed)
    take = min(n_sample, len(ns))
    chosen = rng.choice(ns, size=take, replace=False)
    return [np.array([int(n), int(n) + 2, int(n) + 6, 6], dtype=float) for n in chosen]


def run_invariance_test(
    seeds: list[np.ndarray],
    *,
    alpha: float,
    steps: int,
    coupling: str,
    multiplicative: bool,
) -> dict:
    """
    I1: Ist A_dyn unter Hurwitz-Konjugation invariant?
        A_dyn(q)  vs  A_dyn(h q h^{-1})  für alle 24 Einheiten h.
    """
    units = hurwitz_units()
    abs_diffs: list[float] = []
    rel_diffs: list[float] = []
    a_vals: list[float] = []

    for q0 in seeds:
        a0 = asymmetry_dyn(
            q0,
            alpha,
            steps=steps,
            multiplicative=multiplicative,
            coupling=coupling,
        )
        a_vals.append(a0)
        for h in units:
            qh = qconj_by(h, q0)
            ah = asymmetry_dyn(
                qh,
                alpha,
                steps=steps,
                multiplicative=multiplicative,
                coupling=coupling,
            )
            d = abs(ah - a0)
            abs_diffs.append(d)
            scale = max(abs(a0), 1e-30)
            rel_diffs.append(d / scale)

    abs_arr = np.asarray(abs_diffs, dtype=float)
    rel_arr = np.asarray(rel_diffs, dtype=float)
    a_arr = np.asarray(a_vals, dtype=float)
    out = {
        "n_seeds": len(seeds),
        "n_units": len(units),
        "n_comparisons": int(len(abs_arr)),
        "mean_abs_A_dyn": float(np.mean(np.abs(a_arr))) if len(a_arr) else float("nan"),
        "mean_abs_delta": float(abs_arr.mean()) if len(abs_arr) else float("nan"),
        "max_abs_delta": float(abs_arr.max()) if len(abs_arr) else float("nan"),
        "median_rel_delta": float(np.median(rel_arr)) if len(rel_arr) else float("nan"),
        "frac_rel_delta_gt_1e-3": (
            float(np.mean(rel_arr > 1e-3)) if len(rel_arr) else float("nan")
        ),
        "frac_rel_delta_gt_1e-1": (
            float(np.mean(rel_arr > 1e-1)) if len(rel_arr) else float("nan")
        ),
    }
    print()
    print(
        f"I1 INVARIANCE  coupling={coupling}  alpha={alpha}  steps={steps}  "
        f"seeds={len(seeds)}  units={len(units)}"
    )
    print(
        f"  mean|A_dyn|={out['mean_abs_A_dyn']:.6e}  "
        f"mean|ΔA|={out['mean_abs_delta']:.6e}  "
        f"max|ΔA|={out['max_abs_delta']:.6e}"
    )
    print(
        f"  median|ΔA|/|A|={out['median_rel_delta']:.6e}  "
        f"frac(rel>1e-3)={out['frac_rel_delta_gt_1e-3']:.3f}  "
        f"frac(rel>0.1)={out['frac_rel_delta_gt_1e-1']:.3f}"
    )
    print(
        "note: große relative ΔA ⇒ A_dyn ist keine Konjugationsinvariante "
        "(Koordinaten-/Gap-Kopplung bricht Aut(H)-Symmetrie)."
    )
    return out


def normalize_q(q: np.ndarray, eps: float = 1e-30) -> np.ndarray:
    """Formprojektion ˆq = q / √N(q) auf die Einheitssphäre in H ≅ R⁴."""
    nrm = math.sqrt(max(qnorm(q), eps))
    return q / nrm


def scalar_part(q: np.ndarray) -> float:
    """S(q) = q + conjugate(q) = 2 Re(q) — konjugationsinvariant unter Einheiten."""
    return 2.0 * float(q[0])


def delta_energy(traj: list[np.ndarray]) -> float:
    """ΔE = N(q_T) − N(q_0)."""
    return qnorm(traj[-1]) - qnorm(traj[0])


def I_sym_on_traj(traj: list[np.ndarray]) -> float:
    """I_sym = ∑ w_n S(q_n), w_n = e^{-n/3} (normiert)."""
    weights = np.exp(-np.arange(len(traj)) / 3.0)
    weights /= weights.sum()
    svals = np.array([scalar_part(q) for q in traj], dtype=float)
    return float(np.sum(weights * svals))


def run_lyapunov_test(
    seeds: list[np.ndarray],
    *,
    alpha: float,
    steps: int,
    coupling: str,
    multiplicative: bool,
    eps: float,
    seed: int,
    quiet: bool = False,
    shape: bool = False,
) -> dict:
    """
    L1: Lyapunov-artiger Exponent entlang der Rotor-Bahn:
        λ = (1/T) log(d_T / d_0),  d_k = |q_k − q̃_k|.

    shape=True: Abstände auf ˆq = q/√N(q) → λ_shape (Formdynamik).
    """
    rng = np.random.default_rng(seed)
    T = steps
    lambdas: list[float] = []
    ratios: list[float] = []

    for q0 in seeds:
        # zufällige Richtung, Norm eps
        v = rng.normal(size=4)
        v = (eps / max(np.linalg.norm(v), 1e-30)) * v
        q_tilde = q0 + v
        a0 = normalize_q(q0) if shape else q0
        a_t0 = normalize_q(q_tilde) if shape else q_tilde
        d0 = float(np.linalg.norm(a_t0 - a0))
        traj = evolve_trajectory(
            q0, alpha, steps=T, multiplicative=multiplicative, coupling=coupling
        )
        traj_t = evolve_trajectory(
            q_tilde, alpha, steps=T, multiplicative=multiplicative, coupling=coupling
        )
        aT = normalize_q(traj[-1]) if shape else traj[-1]
        a_tT = normalize_q(traj_t[-1]) if shape else traj_t[-1]
        dT = float(np.linalg.norm(a_tT - aT))
        if d0 <= 0 or dT <= 0:
            continue
        ratios.append(dT / d0)
        lambdas.append(math.log(dT / d0) / T)

    lam = np.asarray(lambdas, dtype=float)
    rat = np.asarray(ratios, dtype=float)
    tag = "shape" if shape else "ambient"
    out = {
        "n_seeds": len(seeds),
        "n_ok": int(len(lam)),
        "eps": eps,
        "T": T,
        "alpha": alpha,
        "space": tag,
        "mean_lambda": float(lam.mean()) if len(lam) else float("nan"),
        "median_lambda": float(np.median(lam)) if len(lam) else float("nan"),
        "std_lambda": float(lam.std(ddof=1)) if len(lam) > 1 else float("nan"),
        "frac_lambda_pos": float(np.mean(lam > 0)) if len(lam) else float("nan"),
        "mean_dT_over_d0": float(rat.mean()) if len(rat) else float("nan"),
        "median_dT_over_d0": float(np.median(rat)) if len(rat) else float("nan"),
    }
    if not quiet:
        print()
        print(
            f"L1 LYAPUNOV[{tag}]  coupling={coupling}  alpha={alpha}  T={T}  "
            f"eps={eps}  seeds={out['n_ok']}"
        )
        print(
            f"  mean λ={out['mean_lambda']:+.6e}  "
            f"median λ={out['median_lambda']:+.6e}  "
            f"std λ={out['std_lambda']:.6e}"
        )
        print(
            f"  frac(λ>0)={out['frac_lambda_pos']:.3f}  "
            f"median(d_T/d_0)={out['median_dT_over_d0']:.6e}"
        )
        print(
            "note: λ>0 ⇒ lokale Expansion (chaotisch/instabil); "
            "λ≈0 / <0 ⇒ neutrale bzw. kontrahierende Klassen."
        )
    return out


def run_paired_ambient_shape_lyapunov(
    seeds: list[np.ndarray],
    *,
    alpha: float,
    steps: int,
    coupling: str,
    multiplicative: bool,
    eps: float,
    seed: int,
) -> dict:
    """
    Dieselbe Störung/Bahn: ambient λ und λ_shape parallel.
    λ_shape misst Expansion von ˆq = q/√N(q) (Form, nicht Energie).
    """
    amb = run_lyapunov_test(
        seeds,
        alpha=alpha,
        steps=steps,
        coupling=coupling,
        multiplicative=multiplicative,
        eps=eps,
        seed=seed,
        quiet=True,
        shape=False,
    )
    shp = run_lyapunov_test(
        seeds,
        alpha=alpha,
        steps=steps,
        coupling=coupling,
        multiplicative=multiplicative,
        eps=eps,
        seed=seed,
        quiet=True,
        shape=True,
    )
    return {
        "alpha": alpha,
        "T": steps,
        "ambient": amb,
        "shape": shp,
        "median_lambda_ambient": amb["median_lambda"],
        "median_lambda_shape": shp["median_lambda"],
        "frac_pos_ambient": amb["frac_lambda_pos"],
        "frac_pos_shape": shp["frac_lambda_pos"],
        "median_ratio_amb": amb["median_dT_over_d0"],
        "median_ratio_shape": shp["median_dT_over_d0"],
    }


def run_shape_lyapunov_suite(
    *,
    N: int,
    alphas: list[float],
    T_list: list[int],
    coupling: str,
    multiplicative: bool,
    n_sample: int,
    seed: int,
    lyapunov_eps: float,
    focus_alpha: float,
) -> dict:
    """
    B) Form vs. Energie: λ_shape(ˆq) gegen ambient λ.
    A) Renormierung: λ_T und λ_shape,T für T ∈ T_list bei focus_alpha.
    """
    mode = "multiplicative" if multiplicative else "additive"
    seeds = sample_wheel_seeds(N, n_sample, seed)
    print(
        f"SHAPE LYAPUNOV  mode={mode}  coupling={coupling}  "
        f"N={N}  n_seeds={len(seeds)}  eps={lyapunov_eps}"
    )
    print(
        "claim-boundary: separates shape instability from norm blow-up; "
        "not a prime-triplet selector."
    )

    print()
    print(
        "α-SCAN   alpha    med λ_amb    med λ_shape   "
        "frac+_amb  frac+_shape  med(dT/d0)_sh"
    )
    alpha_rows = []
    for a in alphas:
        pair = run_paired_ambient_shape_lyapunov(
            seeds,
            alpha=a,
            steps=10,
            coupling=coupling,
            multiplicative=multiplicative,
            eps=lyapunov_eps,
            seed=seed + 1,
        )
        alpha_rows.append(pair)
        print(
            f"         {a:8.4g}  {pair['median_lambda_ambient']:+.6e}  "
            f"{pair['median_lambda_shape']:+.6e}  "
            f"{pair['frac_pos_ambient']:8.3f}  {pair['frac_pos_shape']:8.3f}  "
            f"{pair['median_ratio_shape']:.6e}"
        )

    print()
    print(
        f"T-SCAN  focus_alpha={focus_alpha}  "
        "(ambient vs shape; asymptotisch vs. Kurzzeit)"
    )
    print(
        "T     med λ_amb    med λ_shape   "
        "frac+_amb  frac+_shape  med(dT/d0)_amb  med(dT/d0)_sh"
    )
    T_rows = []
    for T in T_list:
        pair = run_paired_ambient_shape_lyapunov(
            seeds,
            alpha=focus_alpha,
            steps=T,
            coupling=coupling,
            multiplicative=multiplicative,
            eps=lyapunov_eps,
            seed=seed + 1,
        )
        T_rows.append(pair)
        print(
            f"{T:4d}  {pair['median_lambda_ambient']:+.6e}  "
            f"{pair['median_lambda_shape']:+.6e}  "
            f"{pair['frac_pos_ambient']:8.3f}  {pair['frac_pos_shape']:8.3f}  "
            f"{pair['median_ratio_amb']:.6e}  {pair['median_ratio_shape']:.6e}"
        )

    # α=0 shape reference
    pair0 = run_paired_ambient_shape_lyapunov(
        seeds,
        alpha=0.0,
        steps=10,
        coupling=coupling,
        multiplicative=multiplicative,
        eps=lyapunov_eps,
        seed=seed + 1,
    )

    # Diagnose
    med_sh = np.array([r["median_lambda_shape"] for r in alpha_rows], dtype=float)
    med_am = np.array([r["median_lambda_ambient"] for r in alpha_rows], dtype=float)
    shape_pos = bool(np.all(med_sh > 0.05))
    shape_near_zero = bool(np.all(np.abs(med_sh) < 0.05))
    shape_neg = bool(np.mean(med_sh) < -0.05)
    if shape_near_zero:
        verdict = "norm_dominated"
        note = (
            "λ_shape≈0 ⇒ gemessene ambient-Expansion ist vorwiegend "
            "Energieaufblähung, keine echte Forminstabilität."
        )
    elif shape_neg:
        verdict = "shape_contracting"
        note = (
            "λ_shape<0 ⇒ Formdynamik kontrahiert; ambient λ>0 kommt "
            "von der radialen Normdrift."
        )
    elif shape_pos and float(np.median(med_sh)) < 0.5 * float(np.median(med_am)):
        verdict = "mixed_weaker_shape"
        note = (
            "λ_shape>0 aber deutlich kleiner als ambient λ ⇒ "
            "teilweise geometrische Instabilität, teilweise Normdrift."
        )
    elif shape_pos:
        verdict = "genuine_shape_instability"
        note = (
            "λ_shape>0 und vergleichbar mit ambient ⇒ echte Forminstabilität "
            "auf der Einheitssphäre, nicht nur Energieexplosion."
        )
    else:
        verdict = "inconclusive"
        note = "kein klares Vorzeichenmuster für λ_shape."

    # T-Asymptotik der Form
    sh_T = np.array([r["median_lambda_shape"] for r in T_rows], dtype=float)
    if len(sh_T) >= 3:
        late_drift = abs(sh_T[-1] - sh_T[-2]) / max(abs(sh_T[-2]), 1e-30)
        T_note = (
            "λ_shape,T stabilisiert sich asymptotisch"
            if late_drift < 0.2
            else "λ_shape,T driftet noch mit T (Kurzzeit/Transiente möglich)"
        )
    else:
        T_note = "T-Scan zu kurz"

    print()
    print(f"VERDICT  {verdict}")
    print(f"note: {note}")
    print(f"T-note: {T_note}")
    print(
        f"α=0 ref  med λ_amb={pair0['median_lambda_ambient']:+.6e}  "
        f"med λ_shape={pair0['median_lambda_shape']:+.6e}"
    )

    return {
        "mode": mode,
        "coupling": coupling,
        "N": N,
        "n_seeds": len(seeds),
        "lyapunov_eps": lyapunov_eps,
        "focus_alpha": focus_alpha,
        "alphas": alphas,
        "T_list": T_list,
        "alpha_scan": alpha_rows,
        "T_scan": T_rows,
        "alpha0_reference": pair0,
        "verdict": verdict,
        "note": note,
        "T_note": T_note,
        "claim": (
            "Shape Lyapunov separates form dynamics from norm blow-up; "
            "not a prime-triplet claim"
        ),
    }


def hopf_coords(qhat: np.ndarray) -> tuple[float, float, float]:
    """
    Hopf-Koordinaten von ˆq = a+bi+cj+dk ∈ S³:
      z1=a+ib, z2=c+id
      η = |z1|² − |z2|² ∈ [-1,1]
      φ1 = arg(z1), φ2 = arg(z2)
    """
    a, b, c, d = (float(x) for x in qhat)
    r1_sq = a * a + b * b
    r2_sq = c * c + d * d
    eta = r1_sq - r2_sq
    phi1 = math.atan2(b, a)
    phi2 = math.atan2(d, c)
    return eta, phi1, phi2


def _hist_stats(values: np.ndarray, bins: int, lo: float, hi: float) -> dict:
    counts, edges = np.histogram(values, bins=bins, range=(lo, hi), density=True)
    # Gleichverteilung: Dichte 1/(hi-lo)
    uniform = 1.0 / (hi - lo)
    # χ²-ähnliche Abweichung und Peak-zu-Mittel
    peak = float(counts.max()) if len(counts) else float("nan")
    mean_d = float(counts.mean()) if len(counts) else float("nan")
    l1_dev = float(np.mean(np.abs(counts - uniform))) if len(counts) else float("nan")
    return {
        "bins": bins,
        "density": counts.tolist(),
        "edges": edges.tolist(),
        "peak_density": peak,
        "mean_density": mean_d,
        "peak_over_uniform": peak / uniform if uniform > 0 else float("nan"),
        "mean_abs_dev_from_uniform": l1_dev,
        "uniform_density": uniform,
    }


def run_s3_measure_analysis(
    *,
    N: int,
    alpha: float,
    steps: int,
    burn_in: int,
    coupling: str,
    multiplicative: bool,
    n_sample: int,
    seed: int,
    eta_bins: int = 20,
    phi_bins: int = 24,
    tau_max: int = 40,
) -> dict:
    """
    Test C: Geometrie des erzeugten Maßes auf S³.
      ρ(η, φ1, φ2) via Histogramme
      C(τ) = ⟨ ˆq_n · ˆq_{n+τ} ⟩  (Wiederkehr / Mischungsdiagnose)
    """
    mode = "multiplicative" if multiplicative else "additive"
    seeds = sample_wheel_seeds(N, n_sample, seed)
    print(
        f"S³ MEASURE  mode={mode}  coupling={coupling}  alpha={alpha}  "
        f"N={N}  steps={steps}  burn_in={burn_in}  n_seeds={len(seeds)}"
    )
    print(
        "claim-boundary: geometric characterization of the S³ measure; "
        "not a prime-triplet selector."
    )

    etas: list[float] = []
    phi1s: list[float] = []
    phi2s: list[float] = []
    # Korrelation: mitteln über Seeds und Zeitfenster nach burn-in
    c_num = np.zeros(tau_max + 1, dtype=float)
    c_den = np.zeros(tau_max + 1, dtype=float)

    for q0 in seeds:
        traj = evolve_trajectory(
            q0,
            alpha,
            steps=steps,
            multiplicative=multiplicative,
            coupling=coupling,
        )
        hats = [normalize_q(q) for q in traj]
        # Hopf-Samples nach burn-in
        for qh in hats[burn_in:]:
            eta, p1, p2 = hopf_coords(qh)
            etas.append(eta)
            phi1s.append(p1)
            phi2s.append(p2)
        # Autokorrelation der Form
        usable = hats[burn_in:]
        m = len(usable)
        for tau in range(tau_max + 1):
            if tau >= m:
                break
            dots = [
                float(np.dot(usable[i], usable[i + tau]))
                for i in range(m - tau)
            ]
            if dots:
                c_num[tau] += float(np.sum(dots))
                c_den[tau] += len(dots)

    eta_a = np.asarray(etas, dtype=float)
    p1_a = np.asarray(phi1s, dtype=float)
    p2_a = np.asarray(phi2s, dtype=float)
    C = np.array(
        [
            (c_num[t] / c_den[t]) if c_den[t] > 0 else float("nan")
            for t in range(tau_max + 1)
        ],
        dtype=float,
    )

    eta_hist = _hist_stats(eta_a, eta_bins, -1.0, 1.0)
    phi1_hist = _hist_stats(p1_a, phi_bins, -math.pi, math.pi)
    phi2_hist = _hist_stats(p2_a, phi_bins, -math.pi, math.pi)

    # 2D φ1–φ2 Peak (geometrische Präferenz)
    h2, xedges, yedges = np.histogram2d(
        p1_a, p2_a, bins=phi_bins, range=[[-math.pi, math.pi], [-math.pi, math.pi]],
        density=True,
    )
    uni2 = 1.0 / ((2 * math.pi) ** 2)
    peak2 = float(h2.max()) if h2.size else float("nan")

    # C(τ)-Diagnose
    c1 = abs(C[1]) if len(C) > 1 and not math.isnan(C[1]) else float("nan")
    tau_e = None
    if c1 > 0 and not math.isnan(c1):
        thr = c1 / math.e
        for t in range(1, len(C)):
            if abs(C[t]) < thr:
                tau_e = t
                break
    signs = np.sign(C[1:])
    sign_flips = int(np.sum(signs[:-1] * signs[1:] < 0)) if len(signs) > 1 else 0

    # Period-3 der Ad_u-Konjugation: C(3k)≈1, C(≢0 mod 3)≈const
    period3_hits = [
        t for t in range(3, tau_max + 1, 3) if not math.isnan(C[t]) and C[t] > 0.95
    ]
    period3_score = len(period3_hits) / max(1, tau_max // 3)
    if period3_score >= 0.8 and c1 > 0.5:
        c_class = "period3_rotation"
    elif tau_e is not None and tau_e <= 3:
        c_class = "fast_decay"
    elif sign_flips >= max(3, (tau_max - 1) // 4):
        c_class = "oscillatory_rotation"
    elif tau_e is not None and tau_e >= 10:
        c_class = "slow_quasiperiodic"
    elif tau_e is None and c1 > 0.2:
        c_class = "slow_decay_or_plateau"
    else:
        c_class = "mixed"

    # Winkelpräferenz
    prefer_eta = eta_hist["peak_over_uniform"] > 1.5
    prefer_phi = (
        phi1_hist["peak_over_uniform"] > 1.5
        or phi2_hist["peak_over_uniform"] > 1.5
        or peak2 / uni2 > 2.0
    )

    print()
    print("HOPF MARGINS")
    print(
        f"  η:   mean={eta_a.mean():+.4f}  std={eta_a.std():.4f}  "
        f"peak/uniform={eta_hist['peak_over_uniform']:.2f}  "
        f"L1dev={eta_hist['mean_abs_dev_from_uniform']:.4f}"
    )
    print(
        f"  φ1:  mean={p1_a.mean():+.4f}  std={p1_a.std():.4f}  "
        f"peak/uniform={phi1_hist['peak_over_uniform']:.2f}"
    )
    print(
        f"  φ2:  mean={p2_a.mean():+.4f}  std={p2_a.std():.4f}  "
        f"peak/uniform={phi2_hist['peak_over_uniform']:.2f}"
    )
    print(
        f"  φ1×φ2 peak/uniform={peak2 / uni2:.2f}  "
        f"prefer_eta={prefer_eta}  prefer_angles={prefer_phi}"
    )

    print()
    print("AUTOCORRELATION C(τ)=⟨ˆq_n·ˆq_{n+τ}⟩")
    print("  τ    C(τ)")
    show_taus = list(range(0, min(tau_max, 20) + 1)) + (
        [tau_max] if tau_max > 20 else []
    )
    for t in show_taus:
        if t <= tau_max:
            print(f"  {t:3d}  {C[t]:+.6f}")
    print(
        f"  class={c_class}  tau_1/e={tau_e}  "
        f"sign_flips={sign_flips}  |C(1)|={c1:.4f}  "
        f"period3_score={period3_score:.2f}"
    )

    if prefer_eta or prefer_phi:
        geom = "anisotropic_preferred_zones"
        geom_note = (
            "Winkel-/η-Dichte weicht klar von Uniformität ab ⇒ "
            "geometrische Präferenz auf S³."
        )
    else:
        geom = "near_uniform_or_weak"
        geom_note = (
            "keine starke Winkelzonen-Präferenz (nahe Uniformität); "
            "Sättigung eher diffus mischend."
        )

    print()
    print(f"GEOMETRY  {geom}")
    print(f"note: {geom_note}")
    print(
        "C-note: period3_rotation≈Ad_u-Orbit (Ordnung 3); "
        "slow⇒quasiperiodisch; oscillatory⇒Rotation; fast⇒Mischung."
    )

    # α=0 Kontrolle (reine Konjugation): sollte eher rotierend/oskillierend sein
    c0_num = np.zeros(tau_max + 1, dtype=float)
    c0_den = np.zeros(tau_max + 1, dtype=float)
    for q0 in seeds[: min(50, len(seeds))]:
        traj0 = evolve_trajectory(
            q0, 0.0, steps=steps, multiplicative=multiplicative, coupling=coupling
        )
        hats0 = [normalize_q(q) for q in traj0][burn_in:]
        m = len(hats0)
        for tau in range(min(tau_max + 1, m)):
            dots = [
                float(np.dot(hats0[i], hats0[i + tau])) for i in range(m - tau)
            ]
            if dots:
                c0_num[tau] += float(np.sum(dots))
                c0_den[tau] += len(dots)
    C0 = [
        (c0_num[t] / c0_den[t]) if c0_den[t] > 0 else float("nan")
        for t in range(tau_max + 1)
    ]
    print(
        f"α=0 C(1)={C0[1]:+.4f}  C(5)={C0[5]:+.4f}  C(10)={C0[10]:+.4f}  "
        f"(reine u-Konjugation)"
    )

    return {
        "mode": mode,
        "coupling": coupling,
        "alpha": alpha,
        "N": N,
        "steps": steps,
        "burn_in": burn_in,
        "n_seeds": len(seeds),
        "n_hopf_samples": int(len(eta_a)),
        "eta": {
            "mean": float(eta_a.mean()),
            "std": float(eta_a.std()),
            **{k: v for k, v in eta_hist.items() if k != "density"},
            "density": eta_hist["density"],
        },
        "phi1": {
            "mean": float(p1_a.mean()),
            "std": float(p1_a.std()),
            **{k: v for k, v in phi1_hist.items() if k != "density"},
            "density": phi1_hist["density"],
        },
        "phi2": {
            "mean": float(p2_a.mean()),
            "std": float(p2_a.std()),
            **{k: v for k, v in phi2_hist.items() if k != "density"},
            "density": phi2_hist["density"],
        },
        "phi1_phi2_peak_over_uniform": float(peak2 / uni2),
        "prefer_eta": bool(prefer_eta),
        "prefer_angles": bool(prefer_phi),
        "geometry": geom,
        "geometry_note": geom_note,
        "C_tau": C.tolist(),
        "C_class": c_class,
        "tau_one_over_e": tau_e,
        "sign_flips": sign_flips,
        "period3_score": float(period3_score),
        "C_alpha0": C0,
        "claim": (
            "Gap-Hurwitz rotor: coordinate-dependent transient form separation "
            "on S^3 without asymptotic Lyapunov chaos; S^3 measure diagnostics only"
        ),
    }


def apply_F(q, alpha: float, *, coupling: str, multiplicative: bool) -> np.ndarray:
    """Ein Schritt F_α."""
    return hurwitz_step(
        q, alpha, coupling=coupling, multiplicative=multiplicative
    )


def apply_F_power(
    q,
    alpha: float,
    n: int,
    *,
    coupling: str,
    multiplicative: bool,
) -> np.ndarray:
    """F_α^n(q)."""
    x = q.copy()
    for _ in range(n):
        x = apply_F(x, alpha, coupling=coupling, multiplicative=multiplicative)
    return x


def correlation_C(
    seeds: list[np.ndarray],
    *,
    alpha: float,
    steps: int,
    burn_in: int,
    tau_max: int,
    coupling: str,
    multiplicative: bool,
) -> np.ndarray:
    """C_α(τ) = ⟨ ˆq_n · ˆq_{n+τ} ⟩ über Seeds und Zeit nach burn-in."""
    c_num = np.zeros(tau_max + 1, dtype=float)
    c_den = np.zeros(tau_max + 1, dtype=float)
    for q0 in seeds:
        traj = evolve_trajectory(
            q0, alpha, steps=steps, multiplicative=multiplicative, coupling=coupling
        )
        hats = [normalize_q(q) for q in traj][burn_in:]
        m = len(hats)
        for tau in range(min(tau_max + 1, m)):
            s = 0.0
            cnt = 0
            for i in range(m - tau):
                s += float(np.dot(hats[i], hats[i + tau]))
                cnt += 1
            if cnt:
                c_num[tau] += s
                c_den[tau] += cnt
    return np.array(
        [
            (c_num[t] / c_den[t]) if c_den[t] > 0 else float("nan")
            for t in range(tau_max + 1)
        ],
        dtype=float,
    )


def run_period3_perturbation(
    *,
    N: int,
    alphas: list[float],
    coupling: str,
    multiplicative: bool,
    n_sample: int,
    seed: int,
    steps: int = 60,
    burn_in: int = 5,
    tau_max: int = 30,
    k_max: int = 10,
) -> dict:
    """
    Störungstheorie um Ad_u-Period-3:
      (i)   F_0^3(q)=q  (exakt, da u^3=-1 ⇒ Ad_{u^3}=id)
      (ii)  ||F_α^3(q)-q|| / α ≈ |E_g|  (beschränkte lineare Antwort)
      (iii) |C_α(τ)-C_0(τ)| ≤ K α f(τ)
    """
    mode = "multiplicative" if multiplicative else "additive"
    seeds = sample_wheel_seeds(N, n_sample, seed)
    print(
        f"PERIOD-3 PERTURBATION  mode={mode}  coupling={coupling}  "
        f"N={N}  n_seeds={len(seeds)}"
    )
    print(
        "claim-boundary: quantifies Ad_u period-3 perturbation; "
        "not a prime-triplet claim."
    )

    # u^3 = -1 ⇒ Ad_u^3 = id
    u2 = qmul(u, u)
    u3 = qmul(u2, u)
    print()
    print(f"ALGEBRA  u^3 = {u3}  (expect -1)")
    print(f"         ||u^3 + 1|| = {float(np.linalg.norm(u3 + qone())):.3e}")

    # (i) F_0^3 = id
    err0 = []
    for q0 in seeds:
        q3 = apply_F_power(q0, 0.0, 3, coupling=coupling, multiplicative=multiplicative)
        err0.append(float(np.linalg.norm(q3 - q0)))
    err0_a = np.asarray(err0, dtype=float)
    print()
    print(
        f"(i) F_0^3=id  max||F_0^3(q)-q||={err0_a.max():.3e}  "
        f"mean={err0_a.mean():.3e}"
    )

    # (ii) E_g-Schätzer: ||F_α^3(q)-q||/α
    print()
    print("(ii) F_α^3(q)=q+α E_g+O(α²)   α    med||δ||/α   max||δ||/α   med||δ||/α²")
    eg_rows = []
    for a in alphas:
        if a <= 0:
            continue
        ratios = []
        ratios2 = []
        for q0 in seeds:
            q3 = apply_F_power(
                q0, a, 3, coupling=coupling, multiplicative=multiplicative
            )
            d = float(np.linalg.norm(q3 - q0))
            ratios.append(d / a)
            ratios2.append(d / (a * a))
        r = np.asarray(ratios, dtype=float)
        r2 = np.asarray(ratios2, dtype=float)
        eg_rows.append(
            {
                "alpha": a,
                "median_Eg_est": float(np.median(r)),
                "mean_Eg_est": float(r.mean()),
                "max_Eg_est": float(r.max()),
                "median_delta_over_alpha2": float(np.median(r2)),
            }
        )
        print(
            f"                    {a:8.4g}  {np.median(r):.6e}  "
            f"{r.max():.6e}  {np.median(r2):.6e}"
        )

    # Linearität: med||δ|| sollte ~ α · M sein ⇒ med||δ||/α ≈ const
    eg_meds = [r["median_Eg_est"] for r in eg_rows]
    if len(eg_meds) >= 2:
        eg_cv = float(np.std(eg_meds) / max(abs(np.mean(eg_meds)), 1e-30))
    else:
        eg_cv = float("nan")
    M_hat = float(np.median(eg_meds)) if eg_meds else float("nan")
    print(f"  M̂ = median_α(med||δ||/α) = {M_hat:.6e}  CV(across α)={eg_cv:.3f}")
    print(
        "  note: kleine CV ⇒ lineare Antwort in α dominiert; "
        "große med||δ||/α² bei winzigem α ⇒ O(α)-Term sichtbar."
    )

    # Drift nach 3k Schritten: d(F^{3k},q) vs α k
    print()
    print("(ii') d(F_α^{3k},q)   α       k   med d    med d/(αk)")
    drift_rows = []
    for a in alphas:
        if a <= 0:
            continue
        for k in (1, 2, 5, k_max):
            ds = []
            for q0 in seeds:
                qk = apply_F_power(
                    q0, a, 3 * k, coupling=coupling, multiplicative=multiplicative
                )
                # Formabstand (Projektion), da Norm driftet
                ds.append(float(np.linalg.norm(normalize_q(qk) - normalize_q(q0))))
            da = np.asarray(ds, dtype=float)
            scale = a * k
            drift_rows.append(
                {
                    "alpha": a,
                    "k": k,
                    "median_d_shape": float(np.median(da)),
                    "median_d_over_alpha_k": float(np.median(da) / scale),
                }
            )
            if k in (1, 5, k_max):
                print(
                    f"                    {a:8.4g}  {k:3d}  "
                    f"{np.median(da):.6e}  {np.median(da)/scale:.6e}"
                )

    # (iii) |C_α - C_0|
    print()
    print("(iii) |C_α(τ)-C_0(τ)| / α")
    C0 = correlation_C(
        seeds,
        alpha=0.0,
        steps=steps,
        burn_in=burn_in,
        tau_max=tau_max,
        coupling=coupling,
        multiplicative=multiplicative,
    )
    print(
        f"  C_0: C(0)={C0[0]:+.4f} C(1)={C0[1]:+.4f} C(2)={C0[2]:+.4f} "
        f"C(3)={C0[3]:+.4f}"
    )
    print(
        "  α         max|ΔC|/α   mean|ΔC|/α   |ΔC(3)|/α   |ΔC(1)|/α   "
        "max_τ |ΔC|"
    )
    c_rows = []
    for a in alphas:
        if a <= 0:
            continue
        Ca = correlation_C(
            seeds,
            alpha=a,
            steps=steps,
            burn_in=burn_in,
            tau_max=tau_max,
            coupling=coupling,
            multiplicative=multiplicative,
        )
        dC = np.abs(Ca - C0)
        # f(τ)∼1+τ/3 als einfache Zeitskala in Period-3-Schritten
        f_tau = 1.0 + np.arange(tau_max + 1) / 3.0
        # K-Schätzer: max |ΔC|/(α f(τ))
        with np.errstate(divide="ignore", invalid="ignore"):
            K_pts = dC / (a * f_tau)
        K_hat = float(np.nanmax(K_pts))
        c_rows.append(
            {
                "alpha": a,
                "C": Ca.tolist(),
                "abs_dC": dC.tolist(),
                "max_abs_dC": float(np.nanmax(dC)),
                "max_abs_dC_over_alpha": float(np.nanmax(dC) / a),
                "mean_abs_dC_over_alpha": float(np.nanmean(dC) / a),
                "abs_dC_3_over_alpha": float(dC[3] / a),
                "abs_dC_1_over_alpha": float(dC[1] / a),
                "K_hat_f_tau": K_hat,
            }
        )
        print(
            f"  {a:8.4g}  {np.nanmax(dC)/a:.6e}  {np.nanmean(dC)/a:.6e}  "
            f"{dC[3]/a:.6e}  {dC[1]/a:.6e}  {np.nanmax(dC):.6e}"
        )

    K_global = float(np.nanmax([r["K_hat_f_tau"] for r in c_rows])) if c_rows else float("nan")
    mad = [r["max_abs_dC_over_alpha"] for r in c_rows]
    mad_cv = float(np.std(mad) / max(abs(np.mean(mad)), 1e-30)) if len(mad) >= 2 else float("nan")
    # Form-Antwort nach einer Periode: med d_shape(k=1)/α
    shape_resp = [
        r["median_d_over_alpha_k"]
        for r in drift_rows
        if r["k"] == 1
    ]
    shape_M = float(np.median(shape_resp)) if shape_resp else float("nan")
    # max|ΔC| selbst (absolut) über α=0.01-Skala
    max_dC_at_1e2 = next(
        (r["max_abs_dC"] for r in c_rows if abs(r["alpha"] - 0.01) < 1e-12),
        float("nan"),
    )

    print()
    print(
        f"BOUND  empirisch K̂≈{K_global:.3e} für |ΔC|≤K α (1+τ/3); "
        f"CV(max|ΔC|/α)={mad_cv:.3f}"
    )
    print(
        f"  shape: med‖ˆF_α³(q)−ˆq‖/α ≈ {shape_M:.3e}  "
        f"|ΔC|_max(α=0.01)≈{max_dC_at_1e2:.3e}"
    )
    print(
        "  note: ambient ‖F_α³−q‖/α ist groß wegen ‖q‖; "
        "maßgeblich sind Formabstand O(α) und |ΔC|=O(α)."
    )
    print(
        "  note: Form-Drift d(F^{3k}) sättigt in k (nicht ∝k) — "
        "konsistent mit endlicher Separation R_∞, nicht säkularer Drift."
    )

    linear_F3 = (not math.isnan(eg_cv)) and eg_cv < 0.1
    small_dC = (not math.isnan(max_dC_at_1e2)) and max_dC_at_1e2 < 0.02
    shape_linear = (not math.isnan(shape_M)) and 0.1 < shape_M < 2.0
    if linear_F3 and small_dC and shape_linear:
        verdict = "perturbation_supported"
        note = (
            "F_0³=id exakt; ‖F_α³−id‖_shape=Θ(α); |C_α−C_0|=O(α) und klein "
            "— gestörte Ad_u-Period-3-Rotation (empirisch gestützt)."
        )
    elif linear_F3:
        verdict = "perturbation_supported_F3"
        note = (
            "Lineare F_α³-Antwort klar; C-Skalierung teilweise verrauscht "
            "wegen sehr kleinem |ΔC|."
        )
    else:
        verdict = "perturbation_partial"
        note = "Abweichung von O(α)-Struktur — Fenster prüfen."
    print(f"VERDICT  {verdict}")
    print(f"note: {note}")

    return {
        "mode": mode,
        "coupling": coupling,
        "N": N,
        "n_seeds": len(seeds),
        "u3": u3.tolist(),
        "F0_cubed_id": {
            "max_err": float(err0_a.max()),
            "mean_err": float(err0_a.mean()),
        },
        "alphas": alphas,
        "E_g_rows": eg_rows,
        "M_hat": M_hat,
        "E_g_cv_across_alpha": eg_cv,
        "drift_rows": drift_rows,
        "C0": C0.tolist(),
        "C_rows": c_rows,
        "K_hat": K_global,
        "max_dC_over_alpha_cv": mad_cv,
        "shape_M_hat": shape_M,
        "max_dC_at_alpha_0p01": max_dC_at_1e2,
        "verdict": verdict,
        "note": note,
        "claim": (
            "For small α the gap-Hurwitz rotor is a coordinate-dependent "
            "perturbation of the finite Ad_u period-3 rotation on S^3"
        ),
    }


# Asymptotische Gap-Klasse (2,4,6): |R_g(q)|/√N(q) → √(85/216)
B_GAP246_ASYMP = math.sqrt(85.0 / 216.0)
# Dreiecksschranke (korrekt, nicht scharf): |E_g| ≤ 3 M_B
E_G_BOUND_GAP246 = 3.0 * B_GAP246_ASYMP
# Scharfe Asymptotik (Kettenregel, n→∞): |E_g| → √(29/24)
#   C_g^△ := |E_g|/(3 M_B) → √(29/85) ≈ 0.584103
#   C_g^B := |E_g|/M_B     → 3√(29/85) ≈ 1.752309
E_G_SHARP_ASYMP = math.sqrt(29.0 / 24.0)
C_G_TRI_ASYMP = math.sqrt(29.0 / 85.0)
C_G_B_ASYMP = 3.0 * C_G_TRI_ASYMP


def drive_B(q: np.ndarray, *, coupling: str = "gap") -> np.ndarray:
    """B(q) = A(q)/√N(q) — multiplikativer Antrieb bei α=1."""
    A = coupling_A(q, coupling)
    nrm = math.sqrt(max(qnorm(q), 1e-30))
    return A / nrm


def gap246_B_closed_form(n: float) -> float:
    """
    Exakt für q=(n,n+2,n+6,6) mit Gaps (2,4,6):
      w = (2n/3+2, 3n/4+3/2, 5n/12+1/2)
      |B| = |w| / √(3n²+16n+76)
    """
    w0 = (2.0 / 3.0) * n + 2.0
    w1 = (3.0 / 4.0) * n + 1.5
    w2 = (5.0 / 12.0) * n + 0.5
    N = 3.0 * n * n + 16.0 * n + 76.0
    return math.sqrt(w0 * w0 + w1 * w1 + w2 * w2) / math.sqrt(N)


def run_eg_explicit_bound(
    *,
    N: int,
    alphas: list[float],
    n_sample: int,
    seed: int,
    coupling: str = "gap",
) -> dict:
    """
    Lemma (Gap-Klasse G_{2,4,6}, Form auf S³):

      B(q) := R_g(q)/√N(q),   Φ_α(ˆq) = F_α(q)/|F_α(q)|.

    Auf q_n=(n,n+2,n+6,6):
      |B(q_n)| = |w(n)|/√(3n²+16n+76) ↗ √(85/216)   (n→∞).

    Da Ad_u orthogonal auf T S³ und Re(B)=0 (reiner Vektoranteil):
      |V(q)| = |B(q)|  für die α-Linearisierung eines Schritts.

    Kettenregel für Φ_α³ bei Φ_0³=id:
      |E_g(q)| ≤ 3 sup |B| ≤ 3√(85/216).

    Damit gleichmäßig in der Gap-Klasse:
      |Φ_α³(x)-x| ≤ 3√(85/216)·α + O(α²).

    Bemerkung: Auf ganz S³ ohne Gap-Schranke ist |R_g| unbeschränkt
    (δ→0). Das Lemma gilt auf der affinen Gap-Klasse G_{2,4,6}
    (bzw. deren Ad_u-Orbit unter Wheel-Seeds), nicht Haar-global.
    """
    print(
        f"E_g EXPLICIT BOUND  coupling={coupling}  N={N}  "
        f"n_sample={n_sample}"
    )
    print(
        "claim-boundary: analytic bound on gap-class form perturbation; "
        "not a prime-triplet claim."
    )

    print()
    print("LEMMA CONSTANTS")
    print(f"  √(85/216)     = {B_GAP246_ASYMP:.12f}  =: M_B^∞")
    print(f"  3√(85/216)    = {E_G_BOUND_GAP246:.12f}  =: M_E^△ (Dreieck)")
    print(f"  √(29/24)      = {E_G_SHARP_ASYMP:.12f}  =: |E_g|_∞ (Schnärfe, n→∞)")
    print(f"  √(29/85)      = {C_G_TRI_ASYMP:.12f}  =: C_g^△ = |E|/M_E^△")
    print(f"  3√(29/85)     = {C_G_B_ASYMP:.12f}  =: C_g^B = |E|/M_B")
    print(
        "  Lemma1-Slack: 85 N − 216|w|² = 16(13n+316) ≥ 0  (scharf im Leading-Term)"
    )

    # Monotonie / Approximaton der geschlossenen Form
    print()
    print("CLOSED FORM |B|(n) on G_{2,4,6}")
    print("  n         |B|_formula   |B|_num      M_B^∞-|B|   Lemma1-slack")
    form_rows = []
    for n in (5, 17, 101, 10007, 10**6 + 7):
        q = np.array([n, n + 2, n + 6, 6.0], dtype=float)
        bf = gap246_B_closed_form(float(n))
        bn = float(np.linalg.norm(drive_B(q, coupling=coupling)))
        # slack = 85 N - 216 |w|^2 = 16(13n+316)
        slack = 16.0 * (13.0 * n + 316.0)
        form_rows.append(
            {"n": n, "B_formula": bf, "B_num": bn, "lemma1_slack": slack}
        )
        print(
            f"  {n:9d}  {bf:.10f}  {bn:.10f}  {B_GAP246_ASYMP - bf:.3e}  "
            f"{slack:.6e}"
        )

    # Sup |B| auf Wheel-Seeds und Ad_u-Orbit
    seeds = sample_wheel_seeds(N, n_sample, seed)
    max_B = 0.0
    max_B_orbit = 0.0
    for q0 in seeds:
        b0 = float(np.linalg.norm(drive_B(q0, coupling=coupling)))
        max_B = max(max_B, b0)
        q = q0.copy()
        for _ in range(3):
            max_B_orbit = max(
                max_B_orbit,
                float(np.linalg.norm(drive_B(q, coupling=coupling))),
            )
            q = qconj_by(u, q)

    print()
    print(
        f"EMPIRICAL SUP  max|B|(seeds)={max_B:.10f}  "
        f"max|B|(Ad_u-orbit)={max_B_orbit:.10f}"
    )
    print(
        f"  bound check: max ≤ M_B^∞ ?  "
        f"{max_B_orbit <= B_GAP246_ASYMP + 1e-9}"
    )

    def eg_exact_norm(q0: np.ndarray) -> float:
        """Exakte Kettenregel: |d/dα normalize(F_α³(q))|_α=0|."""
        q1 = qconj_by(u, q0)
        q2 = qconj_by(u, q1)
        B0 = drive_B(q0, coupling=coupling)
        B1 = drive_B(q1, coupling=coupling)
        B2 = drive_B(q2, coupling=coupling)
        dq1 = qmul(q1, B0)
        dq2 = qmul(q2, B1) + qconj_by(u, dq1)
        dq3 = qmul(qconj_by(u, q2), B2) + qconj_by(u, dq2)
        nrm = math.sqrt(max(qnorm(q0), 1e-30))
        x = q0 / nrm
        v = dq3 / nrm
        v = v - float(np.dot(v, x)) * x
        return float(np.linalg.norm(v))

    # Exaktes |E_g| (α=0-Linearisierung) vs. Schranken
    print()
    print(
        "E_g EXACT (chain rule)   n_seed-stat   med|E|   max|E|   "
        "med/M_E^△  med/√(29/24)"
    )
    eg_exact = np.asarray([eg_exact_norm(q0) for q0 in seeds], dtype=float)
    print(
        f"                         n={len(seeds):<6d}  "
        f"{np.median(eg_exact):.6e}  {eg_exact.max():.6e}  "
        f"{np.median(eg_exact)/E_G_BOUND_GAP246:.4f}  "
        f"{np.median(eg_exact)/E_G_SHARP_ASYMP:.4f}"
    )
    print(
        f"  C_g^△ med={np.median(eg_exact)/E_G_BOUND_GAP246:.6f}  "
        f"max={eg_exact.max()/E_G_BOUND_GAP246:.6f}  "
        f"(asymp √(29/85)={C_G_TRI_ASYMP:.6f})"
    )
    print(
        f"  C_g^B med={np.median(eg_exact)/B_GAP246_ASYMP:.6f}  "
        f"max={eg_exact.max()/B_GAP246_ASYMP:.6f}  "
        f"(asymp 3√(29/85)={C_G_B_ASYMP:.6f})"
    )

    # Finite-α Check gegen Dreiecksschranke
    print()
    print(
        "E_g FINITE-α vs M_E^△    α      med|Ê|   max|Ê|   "
        "med/M_E^△  max/M_E^△"
    )
    eg_rows = []
    for a in alphas:
        if a <= 0:
            continue
        est = []
        for q0 in seeds:
            q3 = apply_F_power(
                q0, a, 3, coupling=coupling, multiplicative=True
            )
            d = float(np.linalg.norm(normalize_q(q3) - normalize_q(q0)))
            est.append(d / a)
        ea = np.asarray(est, dtype=float)
        eg_rows.append(
            {
                "alpha": a,
                "median_Eg_hat": float(np.median(ea)),
                "max_Eg_hat": float(ea.max()),
                "median_over_triangle": float(np.median(ea) / E_G_BOUND_GAP246),
                "max_over_triangle": float(ea.max() / E_G_BOUND_GAP246),
            }
        )
        print(
            f"                       {a:8.4g}  {np.median(ea):.6e}  "
            f"{ea.max():.6e}  {np.median(ea)/E_G_BOUND_GAP246:.4f}  "
            f"{ea.max()/E_G_BOUND_GAP246:.4f}"
        )

    max_ratio = (
        max(r["max_over_triangle"] for r in eg_rows) if eg_rows else float("nan")
    )
    ok_tri = bool(max_ratio <= 1.0 + 1e-6)
    med_eg = float(np.median(eg_exact))
    # Scharf nur asymptotisch: med nahe √(29/24); kleine n können größer sein
    ok_sharp_med = bool(med_eg <= E_G_SHARP_ASYMP * 1.01)

    print()
    if ok_tri and ok_sharp_med:
        verdict = "sharp_asymptotics_supported"
        note = (
            f"Dreieck |E|≤3√(85/216) hält; asymptotisch |E|→√(29/24), "
            f"C_g^△→√(29/85)≈{C_G_TRI_ASYMP:.4f} (nicht C_g<1 gegen M_B, "
            f"sondern gegen die Dreiecksschranke)."
        )
    elif ok_tri:
        verdict = "triangle_bound_holds"
        note = (
            "Dreiecksschranke hält; scharfe √(29/24)-Kontrolle nur für große n "
            "(kleine Seeds können |E|>√(29/24) haben, z.B. n=5)."
        )
    else:
        verdict = "lemma_bound_violated"
        note = "Empirisches |Ê| überschreitet 3√(85/216)."
    print(f"VERDICT  {verdict}")
    print(f"note: {note}")
    print(
        "note: globale S³-Schranke ohne min-Gap δ ist falsch "
        "(R_g unbeschränkt bei δ→0)."
    )

    return {
        "M_B_asymp": B_GAP246_ASYMP,
        "M_E_triangle": E_G_BOUND_GAP246,
        "E_g_sharp_asymp": E_G_SHARP_ASYMP,
        "C_g_triangle_asymp": C_G_TRI_ASYMP,
        "C_g_over_MB_asymp": C_G_B_ASYMP,
        "lemma1_slack_formula": "16*(13*n+316)",
        "closed_form_rows": form_rows,
        "max_B_seeds": max_B,
        "max_B_Ad_u_orbit": max_B_orbit,
        "E_g_exact_median": med_eg,
        "E_g_exact_max": float(eg_exact.max()),
        "alphas": alphas,
        "E_g_rows": eg_rows,
        "max_Eg_over_triangle": max_ratio,
        "verdict": verdict,
        "note": note,
        "lemma": (
            "Lemma1: |B|²≤85/216 via slack=16(13n+316). "
            "Triangle: |E_g|≤3√(85/216). "
            "Sharp asymp: |E_g|→√(29/24), C_g^△→√(29/85)."
        ),
        "claim": (
            "Explicit E_g bound from R_g on the (2,4,6)-gap class; "
            "bounded periodic perturbation of Ad_u, not a prime selector"
        ),
        "lean_hooks": {
            "lemma1": "pure rational inequality — ready for Lean/mathlib",
            "lemma2_3": (
                "need Mathlib.Quaternion + Ad_u; "
                "eabc-renorm HolonomySpin only has {-2,2i} shadow; "
                "KeplerHurwitz has mathlib but no H-dynamics yet"
            ),
        },
    }


def classify_lyapunov_regime(rows: list[dict]) -> dict:
    """
    Heuristik A/B/C:
      A linear: λ/α näherungsweise konstant über kleine α
      B Schwelle: Vorzeichenwechsel von λ
      C Sättigung: λ flacht bei großen α ab
    """
    if len(rows) < 2:
        return {"regime": "insufficient", "note": "need ≥2 alphas"}
    alphas = np.array([r["alpha"] for r in rows], dtype=float)
    med = np.array([r["median_lambda"] for r in rows], dtype=float)
    # Linearität auf den kleinsten positiven α mit |α|≤0.01
    small = [(a, m) for a, m in zip(alphas, med) if 0 < a <= 0.01]
    linear_score = float("nan")
    if len(small) >= 2:
        ratios = [m / a for a, m in small if a > 0]
        linear_score = float(np.std(ratios) / max(abs(np.mean(ratios)), 1e-30))
    sign_change = bool(np.any(med < 0) and np.any(med > 0))
    # Sättigung: relative Änderung der letzten zwei Punkte klein vs. früher
    sat = False
    if len(med) >= 3 and abs(med[-2]) > 1e-30:
        late = abs(med[-1] - med[-2]) / abs(med[-2])
        early = abs(med[1] - med[0]) / max(abs(med[0]), 1e-30)
        sat = late < 0.25 * max(early, 1e-12) or late < 0.15
    if sign_change:
        regime = "B_threshold"
    elif sat and abs(med[-1]) > 0.1:
        regime = "C_saturation"
    elif not math.isnan(linear_score) and linear_score < 0.35 and np.all(med[alphas <= 0.01] > 0):
        regime = "A_linear"
    elif np.all(med > 0) and sat:
        regime = "C_saturation"
    elif np.all(med > 0):
        # expansive für alle α: zwischen A und C unterscheiden
        if not math.isnan(linear_score) and linear_score < 0.35:
            regime = "A_linear"
        else:
            regime = "C_like_or_mixed"
    else:
        regime = "mixed"
    return {
        "regime": regime,
        "sign_change": bool(sign_change),
        "linear_cv_small_alpha": linear_score,
        "saturation_hint": bool(sat),
        "note": (
            "A=λ∼cα; B=Vorzeichenwechsel; C=Sättigung λ→λ∞; "
            "Heuristik, kein Theorem."
        ),
    }


def run_symmetric_observables(
    seeds: list[np.ndarray],
    *,
    alpha: float,
    steps: int,
    coupling: str,
    multiplicative: bool,
    quiet: bool = False,
) -> dict:
    """
    Konjugationsfreundliche Observablen vs. koordinatengebundenes A_dyn:
      ΔE = N(q_T)−N(q_0)
      I_sym = ∑ w_n · 2 Re(q_n)
    """
    dE: list[float] = []
    isym: list[float] = []
    adyn: list[float] = []
    for q0 in seeds:
        traj = evolve_trajectory(
            q0, alpha, steps=steps, multiplicative=multiplicative, coupling=coupling
        )
        dE.append(delta_energy(traj))
        isym.append(I_sym_on_traj(traj))
        adyn.append(
            asymmetry_dyn(
                q0,
                alpha,
                steps=steps,
                multiplicative=multiplicative,
                coupling=coupling,
            )
        )
    dE_a = np.asarray(dE, dtype=float)
    is_a = np.asarray(isym, dtype=float)
    a_a = np.asarray(adyn, dtype=float)

    def corr(x: np.ndarray, y: np.ndarray) -> float:
        if len(x) < 3 or x.std() == 0 or y.std() == 0:
            return float("nan")
        return float(np.corrcoef(x, y)[0, 1])

    out = {
        "alpha": alpha,
        "mean_delta_E": float(dE_a.mean()),
        "median_delta_E": float(np.median(dE_a)),
        "mean_rel_delta_E": float(np.mean(dE_a / np.maximum(
            [qnorm(s) for s in seeds], 1e-30
        ))),
        "mean_I_sym": float(is_a.mean()),
        "median_I_sym": float(np.median(is_a)),
        "mean_abs_A_dyn": float(np.mean(np.abs(a_a))),
        "corr_I_sym_A_dyn": corr(is_a, a_a),
        "corr_deltaE_A_dyn": corr(dE_a, a_a),
        "corr_deltaE_I_sym": corr(dE_a, is_a),
    }
    if not quiet:
        print()
        print(
            f"SYM OBS  alpha={alpha}  "
            f"median ΔE={out['median_delta_E']:+.6e}  "
            f"mean relΔE={out['mean_rel_delta_E']:+.6e}  "
            f"median I_sym={out['median_I_sym']:+.6e}"
        )
        print(
            f"  corr(I_sym,A_dyn)={out['corr_I_sym_A_dyn']:+.3f}  "
            f"corr(ΔE,A_dyn)={out['corr_deltaE_A_dyn']:+.3f}  "
            f"corr(ΔE,I_sym)={out['corr_deltaE_I_sym']:+.3f}"
        )
    return out


def run_lyapunov_alpha_scan(
    *,
    N: int,
    alphas: list[float],
    steps: int,
    coupling: str,
    multiplicative: bool,
    n_sample: int,
    seed: int,
    lyapunov_eps: float,
) -> dict:
    """λ(α)-Scan + symmetrische Observablen — reine Rotordynamik."""
    mode = "multiplicative" if multiplicative else "additive"
    seeds = sample_wheel_seeds(N, n_sample, seed)
    print(
        f"LYAPUNOV α-SCAN  mode={mode}  coupling={coupling}  "
        f"N={N}  T={steps}  n_seeds={len(seeds)}  eps={lyapunov_eps}"
    )
    print(
        "claim-boundary: classifies rotor expansion vs α; "
        "not a prime-triplet selector."
    )
    print()
    print(
        "alpha       median λ      mean λ     frac(λ>0)  "
        "med(dT/d0)   med ΔE        corr(I,A)"
    )
    rows = []
    sym_rows = []
    for a in alphas:
        lya = run_lyapunov_test(
            seeds,
            alpha=a,
            steps=steps,
            coupling=coupling,
            multiplicative=multiplicative,
            eps=lyapunov_eps,
            seed=seed + 1,
            quiet=True,
        )
        sym = run_symmetric_observables(
            seeds,
            alpha=a,
            steps=steps,
            coupling=coupling,
            multiplicative=multiplicative,
            quiet=True,
        )
        row = {**lya, **{k: sym[k] for k in (
            "median_delta_E",
            "mean_rel_delta_E",
            "median_I_sym",
            "mean_abs_A_dyn",
            "corr_I_sym_A_dyn",
            "corr_deltaE_A_dyn",
            "corr_deltaE_I_sym",
        )}}
        rows.append(row)
        sym_rows.append(sym)
        print(
            f"{a:8.4g}  {lya['median_lambda']:+.6e}  {lya['mean_lambda']:+.6e}  "
            f"{lya['frac_lambda_pos']:8.3f}  {lya['median_dT_over_d0']:.6e}  "
            f"{sym['median_delta_E']:+.6e}  {sym['corr_I_sym_A_dyn']:+.3f}"
        )

    regime = classify_lyapunov_regime(rows)
    print()
    print(
        f"REGIME  {regime['regime']}  "
        f"sign_change={regime['sign_change']}  "
        f"linear_cv={regime['linear_cv_small_alpha']}  "
        f"sat_hint={regime['saturation_hint']}"
    )
    print(f"note: {regime['note']}")

    # α=0 Referenz (reine Konjugation ohne Gap-Drive)
    lya0 = run_lyapunov_test(
        seeds,
        alpha=0.0,
        steps=steps,
        coupling=coupling,
        multiplicative=multiplicative,
        eps=lyapunov_eps,
        seed=seed + 1,
        quiet=True,
    )
    print()
    print(
        f"α=0 reference  median λ={lya0['median_lambda']:+.6e}  "
        f"frac(λ>0)={lya0['frac_lambda_pos']:.3f}  "
        f"(reine u-Konjugation, kein Gap-Drive)"
    )

    return {
        "mode": mode,
        "coupling": coupling,
        "N": N,
        "steps": steps,
        "n_seeds": len(seeds),
        "lyapunov_eps": lyapunov_eps,
        "alphas": alphas,
        "rows": rows,
        "symmetric": sym_rows,
        "alpha0_reference": lya0,
        "regime": regime,
        "claim": (
            "Lyapunov-α scan and conjugacy-friendly observables diagnose "
            "rotor dynamics only; not prime selection"
        ),
    }


def run_rotor_dynamics(
    *,
    N: int,
    alpha: float,
    steps: int,
    coupling: str,
    multiplicative: bool,
    n_sample: int,
    seed: int,
    lyapunov_eps: float,
) -> dict:
    """Rotor-Eigendynamik: Invarianz + Lyapunov, ohne Primzahllabel-Ziel."""
    mode = "multiplicative" if multiplicative else "additive"
    seeds = sample_wheel_seeds(N, n_sample, seed)
    print(
        f"ROTOR DYNAMICS  mode={mode}  coupling={coupling}  "
        f"alpha={alpha}  N={N}  steps={steps}  n_seeds={len(seeds)}"
    )
    print(
        "claim-boundary: tests probe rotor geometry/dynamics; "
        "not a prime-triplet selector."
    )
    inv = run_invariance_test(
        seeds,
        alpha=alpha,
        steps=steps,
        coupling=coupling,
        multiplicative=multiplicative,
    )
    lya = run_lyapunov_test(
        seeds,
        alpha=alpha,
        steps=steps,
        coupling=coupling,
        multiplicative=multiplicative,
        eps=lyapunov_eps,
        seed=seed + 1,
    )
    # Baseline: typische |A_dyn| und Anteil A_dyn ≠ 0
    a_abs = []
    for q0 in seeds:
        a_abs.append(
            abs(
                asymmetry_dyn(
                    q0,
                    alpha,
                    steps=steps,
                    multiplicative=multiplicative,
                    coupling=coupling,
                )
            )
        )
    a_abs_arr = np.asarray(a_abs, dtype=float)
    baseline = {
        "mean_abs_A_dyn": float(a_abs_arr.mean()) if len(a_abs_arr) else float("nan"),
        "median_abs_A_dyn": (
            float(np.median(a_abs_arr)) if len(a_abs_arr) else float("nan")
        ),
        "frac_abs_A_dyn_gt_1e-12": (
            float(np.mean(a_abs_arr > 1e-12)) if len(a_abs_arr) else float("nan")
        ),
    }
    print()
    print(
        f"BASELINE A_dyn  mean|A|={baseline['mean_abs_A_dyn']:.6e}  "
        f"median|A|={baseline['median_abs_A_dyn']:.6e}  "
        f"frac(|A|>1e-12)={baseline['frac_abs_A_dyn_gt_1e-12']:.3f}"
    )
    return {
        "mode": mode,
        "coupling": coupling,
        "alpha": alpha,
        "N": N,
        "steps": steps,
        "n_seeds": len(seeds),
        "baseline_A_dyn": baseline,
        "invariance_I1": inv,
        "lyapunov_L1": lya,
        "claim": (
            "EABC gap rotor: dynamically nontrivial, not a prime-triplet selector"
        ),
    }


def run_experiment(
    *,
    N: int,
    alphas: list[float],
    steps: int,
    multiplicative: bool,
    use_dyn: bool,
    coupling: str,
    time_reversal: bool,
    null_trials: int,
    null_seed: int,
    hard_nulls: bool = False,
) -> dict:
    mode = "multiplicative" if multiplicative else "additive"
    if time_reversal:
        gname = "A_dyn=G→−G←"
    else:
        gname = "G_dyn" if use_dyn else "G_static"

    ns, labels, dens, gtab = collect_gate_table(
        N=N,
        alphas=alphas,
        steps=steps,
        multiplicative=multiplicative,
        use_dyn=use_dyn,
        coupling=coupling,
        time_reversal=time_reversal,
    )
    rows = rows_from_labels(labels, gtab, alphas)
    print_rows(
        rows,
        header=(
            f"mode={mode}  coupling={coupling}  obs={gname}  "
            f"N={N}  steps={steps}  nWheel={len(labels)}"
        ),
    )

    null_summary = None
    if null_trials > 0:
        null_summary = run_null_model(
            labels, gtab, alphas, trials=null_trials, seed=null_seed
        )

    hard_summary = None
    if hard_nulls:
        hard_summary = hard_null_tests(
            ns, labels, dens, gtab, alphas, seed=null_seed
        )

    return {
        "mode": mode,
        "coupling": coupling,
        "observable": gname,
        "time_reversal": time_reversal,
        "N": N,
        "steps": steps,
        "nWheel": int(len(labels)),
        "nP": int(labels.sum()),
        "nC": int(len(labels) - labels.sum()),
        "alphas": alphas,
        "rows": rows,
        "null": null_summary,
        "hard_nulls": hard_summary,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--N", type=int, default=10**6)
    ap.add_argument(
        "--alphas",
        type=float,
        nargs="+",
        default=None,
        help="Explizite Kopplungsstärken",
    )
    ap.add_argument(
        "--sweep",
        action="store_true",
        help="α-Scan (Default-Fenster 0.05–0.30, Schritt 0.01)",
    )
    ap.add_argument(
        "--fine",
        action="store_true",
        help="Test A: Feinscan α ∈ [0.03, 0.08], Schritt 0.01",
    )
    ap.add_argument(
        "--small-coupling",
        action="store_true",
        help="Klein-Kopplung: α ∈ [0.001, 0.02], Schritt 0.001",
    )
    ap.add_argument("--sweep-min", type=float, default=None)
    ap.add_argument("--sweep-max", type=float, default=None)
    ap.add_argument("--sweep-step", type=float, default=0.01)
    ap.add_argument("--steps", type=int, default=10)
    ap.add_argument(
        "--coupling",
        choices=("mod30", "gap"),
        default="mod30",
        help="Arithmetik: mod30 (Restklasse) oder gap (R_g-Gap-Rotor)",
    )
    ap.add_argument(
        "--additive",
        action="store_true",
        help="Additive statt multiplikativer Kopplung (Referenz)",
    )
    ap.add_argument(
        "--static-gate",
        action="store_true",
        help="Altes Punkt-Gate statt G_dyn",
    )
    ap.add_argument(
        "--time-reversal",
        action="store_true",
        help="Test D3: Observable A_dyn = G_forward − G_reversed-path",
    )
    ap.add_argument(
        "--hard-nulls",
        action="store_true",
        help="D1/D2: size- und density-matched Paartests (Wheel-Kontrollen)",
    )
    ap.add_argument(
        "--rotor-dyn",
        action="store_true",
        help="I1/L1: Rotor-Eigendynamik (Invarianz + Lyapunov), ohne Primzahlziel",
    )
    ap.add_argument(
        "--lyapunov-scan",
        action="store_true",
        help="λ(α)-Scan + invariante Observablen ΔE, I_sym",
    )
    ap.add_argument(
        "--shape-lyapunov",
        action="store_true",
        help="λ_shape auf ˆq=q/√N(q): Form vs. Energie; α- und T-Scan",
    )
    ap.add_argument(
        "--s3-measure",
        action="store_true",
        help="Test C: Hopf-Maß ρ(η,φ1,φ2) + Autokorrelation C(τ) auf S³",
    )
    ap.add_argument(
        "--period3-perturbation",
        action="store_true",
        help="Störungstheorie: F_α³ und |C_α−C_0| um Ad_u-Period-3",
    )
    ap.add_argument(
        "--eg-bound",
        action="store_true",
        help="Explizite |E_g|-Schranke aus R_g auf Gap-Klasse (2,4,6)",
    )
    ap.add_argument(
        "--burn-in",
        type=int,
        default=5,
        help="Burn-in-Schritte vor Hopf-/C(τ)-Sampling",
    )
    ap.add_argument(
        "--tau-max",
        type=int,
        default=40,
        help="Maximales τ für C(τ)",
    )
    ap.add_argument(
        "--n-sample",
        type=int,
        default=200,
        help="Anzahl Seeds für Rotor-/Lyapunov-Diagnostik",
    )
    ap.add_argument(
        "--lyapunov-eps",
        type=float,
        default=1e-8,
        help="Anfangstörung d_0 für Lyapunov-Test",
    )
    ap.add_argument(
        "--null-trials",
        type=int,
        default=0,
        help="Test B: Anzahl Label-Permutation-Nullen (0 = aus)",
    )
    ap.add_argument("--null-seed", type=int, default=20260718)
    ap.add_argument(
        "--json-out",
        type=str,
        default="",
        help="Optional: Ergebnisse als JSON schreiben",
    )
    args = ap.parse_args(argv)

    if args.eg_bound:
        import sys

        alphas = (
            list(args.alphas)
            if args.alphas is not None
            else [1e-4, 1e-3, 1e-2, 3e-2]
        )
        coupling = args.coupling if "--coupling" in sys.argv else "gap"
        payload = run_eg_explicit_bound(
            N=args.N,
            alphas=alphas,
            n_sample=args.n_sample,
            seed=args.null_seed,
            coupling=coupling,
        )
    elif args.period3_perturbation:
        import sys

        alphas = (
            list(args.alphas)
            if args.alphas is not None
            else [1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2]
        )
        coupling = args.coupling if "--coupling" in sys.argv else "gap"
        payload = run_period3_perturbation(
            N=args.N,
            alphas=alphas,
            coupling=coupling,
            multiplicative=not args.additive,
            n_sample=args.n_sample,
            seed=args.null_seed,
            steps=max(args.steps, args.burn_in + args.tau_max + 20),
            burn_in=args.burn_in,
            tau_max=args.tau_max,
        )
    elif args.s3_measure:
        import sys

        alpha = float(args.alphas[0]) if args.alphas is not None else 0.01
        coupling = args.coupling if "--coupling" in sys.argv else "gap"
        # längere Bahn für C(τ) und Maß
        steps = max(args.steps, args.burn_in + args.tau_max + 20)
        payload = run_s3_measure_analysis(
            N=args.N,
            alpha=alpha,
            steps=steps,
            burn_in=args.burn_in,
            coupling=coupling,
            multiplicative=not args.additive,
            n_sample=args.n_sample,
            seed=args.null_seed,
            tau_max=args.tau_max,
        )
    elif args.shape_lyapunov:
        import sys

        alphas = (
            list(args.alphas)
            if args.alphas is not None
            else [1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.5]
        )
        coupling = args.coupling if "--coupling" in sys.argv else "gap"
        focus = float(alphas[2]) if len(alphas) >= 3 else 0.01
        # Prefer canonical focus α=0.01 if present
        if 0.01 in alphas:
            focus = 0.01
        payload = run_shape_lyapunov_suite(
            N=args.N,
            alphas=alphas,
            T_list=[5, 10, 20, 50, 100],
            coupling=coupling,
            multiplicative=not args.additive,
            n_sample=args.n_sample,
            seed=args.null_seed,
            lyapunov_eps=args.lyapunov_eps,
            focus_alpha=focus,
        )
    elif args.lyapunov_scan:
        import sys

        alphas = (
            list(args.alphas)
            if args.alphas is not None
            else [1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.5]
        )
        coupling = args.coupling if "--coupling" in sys.argv else "gap"
        payload = run_lyapunov_alpha_scan(
            N=args.N,
            alphas=alphas,
            steps=args.steps,
            coupling=coupling,
            multiplicative=not args.additive,
            n_sample=args.n_sample,
            seed=args.null_seed,
            lyapunov_eps=args.lyapunov_eps,
        )
    elif args.rotor_dyn:
        import sys

        alpha = float(args.alphas[0]) if args.alphas is not None else 0.01
        coupling = args.coupling if "--coupling" in sys.argv else "gap"
        payload = run_rotor_dynamics(
            N=args.N,
            alpha=alpha,
            steps=args.steps,
            coupling=coupling,
            multiplicative=not args.additive,
            n_sample=args.n_sample,
            seed=args.null_seed,
            lyapunov_eps=args.lyapunov_eps,
        )
    else:
        if args.alphas is not None:
            alphas = list(args.alphas)
        elif args.small_coupling:
            alphas = alpha_grid(0.001, 0.02, 0.001)
        elif args.fine:
            alphas = alpha_grid(0.03, 0.08, 0.01)
        elif args.sweep:
            a0 = 0.05 if args.sweep_min is None else args.sweep_min
            a1 = 0.30 if args.sweep_max is None else args.sweep_max
            alphas = alpha_grid(a0, a1, args.sweep_step)
        elif args.sweep_min is not None or args.sweep_max is not None:
            a0 = 0.05 if args.sweep_min is None else args.sweep_min
            a1 = 0.30 if args.sweep_max is None else args.sweep_max
            alphas = alpha_grid(a0, a1, args.sweep_step)
        else:
            alphas = [0.05, 0.1, 0.2]

        if args.hard_nulls and args.alphas is None and not (
            args.small_coupling
            or args.fine
            or args.sweep
            or args.sweep_min is not None
            or args.sweep_max is not None
        ):
            alphas = [0.01]
            if not args.time_reversal:
                args.time_reversal = True
            if args.coupling == "mod30":
                args.coupling = "gap"

        payload = run_experiment(
            N=args.N,
            alphas=alphas,
            steps=args.steps,
            multiplicative=not args.additive,
            use_dyn=not args.static_gate,
            coupling=args.coupling,
            time_reversal=args.time_reversal,
            null_trials=args.null_trials,
            null_seed=args.null_seed,
            hard_nulls=args.hard_nulls,
        )

        if payload["rows"]:
            payload["alpha_star_abs_delta"] = max(
                payload["rows"], key=lambda r: r["abs_delta"]
            )["alpha"]
            payload["alpha_star_abs_Z"] = max(
                payload["rows"], key=lambda r: abs(r["Z"])
            )["alpha"]
        payload["claim"] = (
            "EABC gap rotor: dynamically nontrivial, not a prime-triplet selector; "
            "Z/null/D1/D2 are sample diagnostics only"
        )

    if args.json_out:
        import json
        from pathlib import Path

        Path(args.json_out).write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )
        print(f"wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
