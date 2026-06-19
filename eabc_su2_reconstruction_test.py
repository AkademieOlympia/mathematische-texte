#!/usr/bin/env python3
"""
EABC SU(2)-Rekonstruktionstest
==============================

Prüft die Forschungsfrage aus collatz_eabc_zirkulationshypothese.md §4.8.12:
Kann endliche EABC-/Ikosaeder-Orientierungsgeometrie eine SU(2)-artige
Drehimpulsstruktur — insbesondere (l, m) mit Degenerierung 2l+1 —
approximieren oder rekonstruieren?

Operationalisierung (§4.8.12):
    1. Baseline: Ikosaeder-Ecken V_12 (Default) oder Tetraeder V_4
    2. Hebung: C_E: N -> Z_4 und Zuordnung Z_4 -> V (strukturiert vs. Zufall)
    3. SU(2)-Vergleich: Partitionen von V in „l-Niveaus" mit Kardinalitaet 2l+1
    4. Resonanz: R_ij = u_i . u_j, Casimir-ähnliche Observable, Wigner-D-Analog
    5. Klassifikation: Fall (1) SU(2)-artig, (2) nur spezielle Zuordnung, (3) keine

Wahl V (Default icosahedron):
    12 Ikosaeder-Ecken auf S^2 — nicht 30 Kantenmittelpunkte — weil
    eabc_icosahedron_test.py bereits die geometrische EABC-Zuordnung 3+3+3+3
    auf Ecken definiert (mod-12-Restklassen, A_5-Gitter). Kantenmittelpunkte
    hätten keine etablierte EABC-Hebung.

Hebung Z_4 -> V (Default structured):
    assign_eabc_geometric liefert je 3 Ecken pro Buchstabe E,A,B,C.
    Für Prim p mit Klasse c und Stromindex k:
        u(p) = pool[c][k mod 3]
    Nullmodell: uniform zufällige Ecke pro Schritt (unabhängig von c).

Epistemik: diagnostisch/explorativ — kein Claim, dass SU(2) reproduziert wird.

Querverweise:
    - collatz_eabc_zirkulationshypothese.md §4.8.8–§4.8.12
    - eabc_icosahedron_test.py
    - Bahndrehimpuls.tex (Drehimpuls auf C_4, Auswahlregeln)
"""

from __future__ import annotations

import argparse
import json
from math import isqrt, lgamma
from pathlib import Path

import numpy as np

from eabc_icosahedron_test import (
    assign_eabc_geometric,
    icosahedron_edges,
    icosahedron_vertices,
    normalize,
    triangular_faces,
)

# ---------------------------------------------------------------------------
# Konstanten
# ---------------------------------------------------------------------------

LABEL_TO_CODE: dict[str, int] = {"E": 0, "A": 1, "B": 2, "C": 3}
CODE_TO_LABEL: dict[int, str] = {v: k for k, v in LABEL_TO_CODE.items()}
RESIDUE_MAP: dict[int, str] = {1: "E", 5: "A", 7: "B", 11: "C"}
DEFAULT_LIMIT = 10_000
DEFAULT_SEED = 42
AXIS_Z = np.array([0.0, 0.0, 1.0])

CASE_SU2_LIKE = "case_1_su2_like_structure"
CASE_SPECIAL = "case_2_special_assignment_only"
CASE_NO_REPRO = "case_3_no_reproduction"


# ---------------------------------------------------------------------------
# Orientierungsgitter V
# ---------------------------------------------------------------------------


def tetrahedron_vertices() -> np.ndarray:
    """4 reguläre Tetraeder-Ecken auf S^2 (Alternative zu Ikosaeder)."""
    raw = np.array(
        [
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ],
        dtype=float,
    )
    return np.array([normalize(v) for v in raw])


def build_orientation_set(lift: str) -> tuple[np.ndarray, dict]:
    """
    Liefert Einheitsrichtungen V und Metadaten.

    lift='icosahedron': 12 Ecken (Default, EABC-geometrisch 3+3+3+3)
    lift='tetrahedron': 4 Ecken (Minimaltest, nur l=0,1 teilweise)
    """
    if lift == "icosahedron":
        verts = icosahedron_vertices()
        edges = icosahedron_edges(verts)
        faces = triangular_faces(verts, edges)
        assignment, residues = assign_eabc_geometric(verts, edges, faces)
        meta = {
            "lift": lift,
            "n_vertices": len(verts),
            "symmetry_note": "A_5 (60), endliche Rotationsgruppe",
            "eabc_balance": {letter: list(assignment.values()).count(letter) for letter in "EABC"},
            "residues_mod12": {int(i): int(r) for i, r in residues.items()},
        }
        return verts, meta

    if lift == "tetrahedron":
        verts = tetrahedron_vertices()
        # Tetraeder: zyklische EABC-Zuordnung (keine mod-12-Geometrie auf S^2)
        assignment = {i: "EABC"[i % 4] for i in range(len(verts))}
        meta = {
            "lift": lift,
            "n_vertices": len(verts),
            "symmetry_note": "A_4 (12), Tetraedergruppe",
            "eabc_balance": {letter: list(assignment.values()).count(letter) for letter in "EABC"},
            "assignment_rule": "zyklisch i mod 4 -> E,A,B,C",
        }
        return verts, meta

    raise ValueError(f"Unbekannter lift-Modus: {lift!r} (erwartet icosahedron|tetrahedron)")


def z4_vertex_pools(vertices: np.ndarray, lift: str) -> dict[str, list[int]]:
    """
    Strukturierte Hebung Z_4 -> Teilmengen von V.

    icosahedron: inverse der geometrischen assign_eabc_geometric-Zuordnung.
    tetrahedron: je eine Ecke pro Buchstabe (E,A,B,C).
    """
    pools: dict[str, list[int]] = {letter: [] for letter in "EABC"}

    if lift == "icosahedron":
        edges = icosahedron_edges(vertices)
        faces = triangular_faces(vertices, edges)
        assignment, _ = assign_eabc_geometric(vertices, edges, faces)
        for idx, letter in assignment.items():
            pools[letter].append(idx)
        for letter in "EABC":
            pools[letter].sort()
        return pools

    # tetrahedron: je 1 Ecke pro Klasse
    for i in range(len(vertices)):
        pools["EABC"[i % 4]].append(i)
    return pools


# ---------------------------------------------------------------------------
# Primstrom / EABC-Wort
# ---------------------------------------------------------------------------


def sieve_primes(limit: int) -> list[int]:
    if limit < 2:
        return []
    flags = bytearray([1]) * (limit + 1)
    flags[0] = flags[1] = 0
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p * p :: p] = bytearray(len(flags[p * p :: p]))
    return [i for i, ok in enumerate(flags) if ok]


def class_of(p: int) -> str | None:
    return RESIDUE_MAP.get(p % 12)


def build_prime_eabc_stream(n_primes: int) -> tuple[np.ndarray, np.ndarray]:
    """Erste n_primes EABC-Labels und Primzahlen > 3."""
    limit = max(50_000, n_primes * 20)
    while True:
        primes_all = sieve_primes(limit)
        codes: list[int] = []
        prime_list: list[int] = []
        for p in primes_all:
            if p <= 3:
                continue
            c = class_of(p)
            if c is not None:
                codes.append(LABEL_TO_CODE[c])
                prime_list.append(p)
            if len(codes) >= n_primes:
                return (
                    np.array(codes[:n_primes], dtype=np.int8),
                    np.array(prime_list[:n_primes], dtype=np.int64),
                )
        limit = int(limit * 1.5)


def lift_orientations(
    eabc_codes: np.ndarray,
    pools: dict[str, list[int]],
    *,
    structured: bool = True,
    seed: int = DEFAULT_SEED,
) -> np.ndarray:
    """
    Hebung auf Vertex-Indizes in V.

    structured=True: pool[letter][k mod |pool|] mit Stromindex k.
    structured=False: uniform Zufall auf alle |V|-Ecken (Nullmodell).
    """
    n = len(eabc_codes)
    n_verts = sum(len(p) for p in pools.values())
    indices = np.zeros(n, dtype=np.int32)

    if structured:
        for k, code in enumerate(eabc_codes):
            letter = CODE_TO_LABEL[int(code)]
            pool = pools[letter]
            indices[k] = pool[k % len(pool)]
        return indices

    rng = np.random.default_rng(seed)
    all_indices = np.concatenate([np.array(p, dtype=np.int32) for p in pools.values()])
    for k in range(n):
        indices[k] = int(rng.choice(all_indices))
    return indices


# ---------------------------------------------------------------------------
# SU(2)-Degenerierung 2l+1
# ---------------------------------------------------------------------------


def su2_level_dimensions(l_max: int) -> list[int]:
    return [2 * ell + 1 for ell in range(l_max + 1)]


def partition_vertices_su2_shells(
    vertices: np.ndarray,
    axis: np.ndarray | None = None,
) -> list[dict]:
    """
    Partitioniert V in l-Schalen durch Sortierung nach Polwinkel (u·axis).

    Greedy: nacheinander 1, 3, 5, 7, ... Ecken (SU(2)-Degenerierung 2l+1)
    bis alle Ecken verteilt sind. Letzte Schale kann unvollständig sein.
    """
    axis = normalize(axis if axis is not None else AXIS_Z)
    dots = vertices @ axis
    order = np.argsort(-dots)  # Nordpol → Südpol

    shells: list[dict] = []
    pos = 0
    ell = 0
    n = len(vertices)
    while pos < n:
        target = 2 * ell + 1
        take = min(target, n - pos)
        idx = order[pos : pos + take].tolist()
        shells.append(
            {
                "l": ell,
                "target_dim": target,
                "actual_dim": take,
                "vertex_indices": idx,
                "matches_su2": take == target,
            }
        )
        pos += take
        ell += 1
    return shells


def static_degeneracy_match(shells: list[dict]) -> dict:
    """Vergleicht statische Schalengrößen mit SU(2)-Vorhersage 2l+1."""
    actual = [s["actual_dim"] for s in shells]
    target = [s["target_dim"] for s in shells]
    full_match = all(s["matches_su2"] for s in shells if s["actual_dim"] > 0)
    partial_levels = sum(1 for s in shells if s["matches_su2"])
    return {
        "shell_dims_actual": actual,
        "shell_dims_su2_target": target,
        "full_su2_partition": full_match,
        "n_matching_levels": partial_levels,
        "n_shells": len(shells),
    }


def empirical_shell_occupancy(
    vertex_indices: np.ndarray,
    shells: list[dict],
) -> dict:
    """Zählt Besetzung der l-Schalen im gehobenen Primstrom."""
    index_to_shell = {}
    for sh in shells:
        for vi in sh["vertex_indices"]:
            index_to_shell[vi] = sh["l"]

    counts = {sh["l"]: 0 for sh in shells}
    for vi in vertex_indices:
        counts[index_to_shell[int(vi)]] += 1

    total = len(vertex_indices)
    fracs = {str(l): counts[l] / total for l in counts}
    su2_weights = {str(sh["l"]): sh["actual_dim"] for sh in shells}
    w_sum = sum(su2_weights.values())
    expected_uniform = {k: v / w_sum for k, v in su2_weights.items()}

    chi2 = 0.0
    for sh in shells:
        l = sh["l"]
        obs = counts[l]
        exp = expected_uniform[str(l)] * total
        if exp > 0:
            chi2 += (obs - exp) ** 2 / exp

    return {
        "counts": {str(k): v for k, v in counts.items()},
        "fractions": fracs,
        "expected_su2_uniform": expected_uniform,
        "chi2_vs_su2_uniform": float(chi2),
        "total": total,
    }


# ---------------------------------------------------------------------------
# Korrelationsmatrix, Winkelhistogramm, Observablen
# ---------------------------------------------------------------------------


def representation_correlation_matrix(vertices: np.ndarray) -> np.ndarray:
    """R_ij = u_i · u_j (Resonanzmatrix aus §4.8.10)."""
    return vertices @ vertices.T


def transition_correlation_matrix(
    vertices: np.ndarray,
    stream_indices: np.ndarray,
) -> np.ndarray:
    """Empirische Übergangsgewichte T_ij = Häufigkeit i -> j."""
    n = len(vertices)
    T = np.zeros((n, n), dtype=np.float64)
    for a, b in zip(stream_indices[:-1], stream_indices[1:]):
        T[int(a), int(b)] += 1.0
    if T.sum() > 0:
        T /= T.sum()
    return T


def angle_histogram(vertices: np.ndarray, stream_indices: np.ndarray, n_bins: int = 18) -> dict:
    """Winkelverteilung aufeinanderfolgender Orientierungen arccos(u_k · u_{k+1})."""
    u = vertices[stream_indices]
    dots = np.sum(u[:-1] * u[1:], axis=1)
    dots = np.clip(dots, -1.0, 1.0)
    angles = np.arccos(dots)
    hist, edges = np.histogram(angles, bins=n_bins, range=(0.0, np.pi))
    return {
        "bin_edges": edges.tolist(),
        "counts": hist.astype(int).tolist(),
        "mean_angle_rad": float(np.mean(angles)),
        "std_angle_rad": float(np.std(angles)),
    }


def casimir_like_observable(R: np.ndarray) -> float:
    """Casimir-Analog: mittleres R_ij^2 über alle Paare (diagonal = 1)."""
    return float(np.mean(R ** 2))


def wigner_d_projection_diagnostic(
    vertices: np.ndarray,
    l_max: int = 3,
    axis: np.ndarray | None = None,
) -> dict:
    """
    Wigner-D-Matrix-Analog: projiziere Y_l^m (reell, Condon-Shortley) auf V.

    Für jedes l vergleichen wir die Korrelationsstruktur
        S_l(i,j) = sum_m Y_l^m(u_i) Y_l^m(u_j)
    mit R_ij = u_i·u_j (Frobenius-Korrelation pro l-Schale).
    """
    axis = normalize(axis if axis is not None else AXIS_Z)
    R = representation_correlation_matrix(vertices)
    n = len(vertices)

    # sphärische Koordinaten relativ zu axis
    def sph_coords(v: np.ndarray) -> tuple[float, float]:
        v = normalize(v)
        cos_theta = np.clip(float(np.dot(v, axis)), -1.0, 1.0)
        theta = np.arccos(cos_theta)
        # Azimuth: Projektion in Ebene orthogonal zu axis
        proj = v - cos_theta * axis
        norm = np.linalg.norm(proj)
        if norm < 1e-12:
            phi = 0.0
        else:
            ref = normalize(np.cross(axis, [1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.cross(axis, [0.0, 1.0, 0.0]))
            tangent = normalize(np.cross(axis, ref))
            d = proj / norm
            phi = np.arctan2(float(np.dot(d, tangent)), float(np.dot(d, ref)))
        return theta, phi

    thetas = np.zeros(n)
    phis = np.zeros(n)
    for i, v in enumerate(vertices):
        thetas[i], phis[i] = sph_coords(v)

    try:
        from scipy.special import lpmv
    except ImportError:
        lpmv = None

    def assoc_legendre(l: int, m: int, x: np.ndarray) -> np.ndarray:
        x = np.clip(np.asarray(x, dtype=float), -1.0, 1.0)
        if lpmv is not None:
            return lpmv(m, l, x)
        m_abs = abs(m)
        if m_abs > l:
            return np.zeros_like(x)
        if m_abs == 0:
            p_mm = np.ones_like(x)
        else:
            df = float(np.prod([2 * k - 1 for k in range(1, m_abs + 1)]))
            p_mm = df * (1.0 - x ** 2) ** (m_abs / 2.0)
        if l == m_abs:
            return p_mm
        p_lm1 = x * (2 * m_abs + 1) * p_mm
        if l == m_abs + 1:
            return p_lm1
        p_lm2 = p_mm
        for ell in range(m_abs + 2, l + 1):
            p_l = ((2 * ell - 1) * x * p_lm1 - (ell + m_abs - 1) * p_lm2) / (ell - m_abs)
            p_lm2, p_lm1 = p_lm1, p_l
        return p_lm1

    def real_sph_harm(l: int, m: int, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
        cos_t = np.cos(theta)
        m_abs = abs(m)
        norm = np.sqrt(
            (2 * l + 1) / (4 * np.pi) * np.exp(lgamma(l - m_abs + 1) - lgamma(l + m_abs + 1))
        )
        plm = assoc_legendre(l, m_abs, cos_t)
        if m > 0:
            return norm * plm * np.cos(m * phi)
        if m < 0:
            return norm * plm * np.sin(m_abs * phi)
        return norm * plm

    levels = []
    for ell in range(l_max + 1):
        Y = np.zeros((n, 2 * ell + 1))
        for mp, m in enumerate(range(-ell, ell + 1)):
            Y[:, mp] = real_sph_harm(ell, m, thetas, phis)
        S_l = Y @ Y.T
        denom = np.linalg.norm(S_l) * np.linalg.norm(R)
        corr = float(np.sum(S_l * R) / denom) if denom > 0 else 0.0
        levels.append({"l": ell, "dim": 2 * ell + 1, "frobenius_corr_with_R": corr})

    return {"levels": levels, "l_max": l_max}


# ---------------------------------------------------------------------------
# Klassifikation und Haupttest
# ---------------------------------------------------------------------------


def classify_reconstruction(
    static: dict,
    structured_occ: dict,
    random_occ: dict,
    wigner: dict,
) -> str:
    """
    Fall (1): statische 2l+1-Partition vollständig + strukturierte Hebung näher an SU(2)
    Fall (2): statische Partition teilweise, nur strukturierte Hebung verbessert Chi²
    Fall (3): keine SU(2)-Degenerationsstruktur reproduziert
    """
    chi2_struct = structured_occ["chi2_vs_su2_uniform"]
    chi2_rand = random_occ["chi2_vs_su2_uniform"]
    structured_better = chi2_struct < chi2_rand * 0.95

    wigner_corr = max((lv["frobenius_corr_with_R"] for lv in wigner["levels"]), default=0.0)

    if static["full_su2_partition"] and structured_better and wigner_corr > 0.3:
        return CASE_SU2_LIKE
    if static["n_matching_levels"] >= 2 and structured_better:
        return CASE_SPECIAL
    return CASE_NO_REPRO


def run_su2_reconstruction_test(
    n_primes: int = DEFAULT_LIMIT,
    lift: str = "icosahedron",
    seed: int = DEFAULT_SEED,
    l_max: int = 3,
) -> dict:
    """Führt den vollständigen diagnostischen SU(2)-Rekonstruktionstest durch."""
    vertices, v_meta = build_orientation_set(lift)
    pools = z4_vertex_pools(vertices, lift)
    eabc_codes, primes = build_prime_eabc_stream(n_primes)

    idx_struct = lift_orientations(eabc_codes, pools, structured=True)
    idx_random = lift_orientations(eabc_codes, pools, structured=False, seed=seed)

    shells = partition_vertices_su2_shells(vertices)
    static = static_degeneracy_match(shells)
    occ_struct = empirical_shell_occupancy(idx_struct, shells)
    occ_random = empirical_shell_occupancy(idx_random, shells)

    R = representation_correlation_matrix(vertices)
    T = transition_correlation_matrix(vertices, idx_struct)
    angles = angle_histogram(vertices, idx_struct)
    casimir = casimir_like_observable(R)
    wigner = wigner_d_projection_diagnostic(vertices, l_max=min(l_max, 5))

    classification = classify_reconstruction(static, occ_struct, occ_random, wigner)

    su2_dims = su2_level_dimensions(l_max)
    return {
        "meta": {
            "script": "eabc_su2_reconstruction_test.py",
            "section": "§4.8.12",
            "epistemology": "diagnostisch/explorativ — kein SU(2)-Reproduktions-Claim",
            "n_primes": n_primes,
            "lift": lift,
            "seed": seed,
        },
        "orientation_set": v_meta,
        "lift_rule": {
            "structured": "pool[EABC( p mod 12 )][stream_index mod |pool|]",
            "null_model": "uniform random vertex per step",
            "pools": {k: v for k, v in pools.items()},
        },
        "static_degeneracy": static,
        "shells": shells,
        "su2_reference_dims_2l_plus_1": su2_dims,
        "empirical_occupancy": {
            "structured": occ_struct,
            "random_null": occ_random,
        },
        "correlation": {
            "R_matrix_shape": list(R.shape),
            "R_casimir_like": casimir,
            "R_offdiag_mean": float(np.mean(R[np.triu_indices(len(R), k=1)])),
            "transition_matrix_nonzero": int(np.count_nonzero(T)),
        },
        "angle_histogram": angles,
        "wigner_d_analog": wigner,
        "classification": classification,
        "prime_sample_head": primes[:10].tolist(),
        "stream_indices_head": idx_struct[:10].tolist(),
    }


def print_summary(report: dict) -> None:
    """Kompakte stdout-Zusammenfassung."""
    meta = report["meta"]
    static = report["static_degeneracy"]
    occ = report["empirical_occupancy"]["structured"]
    wigner = report["wigner_d_analog"]

    print("=== EABC SU(2)-Rekonstruktionstest (§4.8.12) ===")
    print(f"Lift: {meta['lift']}, Primzahl-Schritte: {meta['n_primes']}")
    print(f"|V| = {report['orientation_set']['n_vertices']}")
    print()
    print("Statische l-Schalen (Greedy 2l+1 auf Polwinkel):")
    print(f"  Ist-Dimensionen:  {static['shell_dims_actual']}")
    print(f"  SU(2)-Ziel:       {static['shell_dims_su2_target']}")
    print(f"  Volle 2l+1-Partition: {static['full_su2_partition']}")
    print()
    print("Empirische Schalenbesetzung (strukturierte Hebung):")
    print(f"  counts: {occ['counts']}")
    print(f"  chi² vs. SU(2)-uniform: {occ['chi2_vs_su2_uniform']:.4f}")
    print(f"  Null chi² (Zufall):     {report['empirical_occupancy']['random_null']['chi2_vs_su2_uniform']:.4f}")
    print()
    print("Wigner-D-Analog (Frobenius-Korrelation R_ij mit Y_l-Projektion):")
    for lv in wigner["levels"]:
        print(f"  l={lv['l']}: corr={lv['frobenius_corr_with_R']:.4f}")
    print()
    print(f"Casimir-ähnlich ⟨R²⟩: {report['correlation']['R_casimir_like']:.4f}")
    print(f"Klassifikation: {report['classification']}")
    print()
    print("Hinweis: diagnostisch — keine Behauptung SU(2)-Reproduktion.")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="SU(2)-Rekonstruktionstest für EABC-Orientierungshebung (§4.8.12)"
    )
    parser.add_argument(
        "--limit",
        "--n-primes",
        dest="n_primes",
        type=int,
        default=DEFAULT_LIMIT,
        help=f"Anzahl Primzahlen > 3 im EABC-Strom (default: {DEFAULT_LIMIT})",
    )
    parser.add_argument(
        "--lift",
        choices=("icosahedron", "tetrahedron"),
        default="icosahedron",
        help="Orientierungsgitter V (default: icosahedron = 12 Ecken)",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--l-max",
        type=int,
        default=3,
        help="Maximales l für Wigner-D-Analog (default: 3)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="JSON-Ausgabedatei (optional)",
    )
    args = parser.parse_args()

    report = run_su2_reconstruction_test(
        n_primes=args.n_primes,
        lift=args.lift,
        seed=args.seed,
        l_max=args.l_max,
    )
    print_summary(report)

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
        print(f"JSON gespeichert: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
