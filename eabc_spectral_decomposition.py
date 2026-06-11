#!/usr/bin/env python3
"""
EABC spektrale Zerlegung des Kanal-Bias δ auf Dirichlet-Charaktere mod 420.

Kanäle: r = [11, 101, 191, 221, 311, 401]
Counts: [765, 831, 786, 809, 809, 767], N = 4767

Berechnet:
  - χ₃ mod 7 (kubisch) und quadratische χ₄, χ₅, χ₇ sowie Produkte
  - Orthogonale Charakter-Projektion von δ
  - Rest-Bias-Aufschlüsselung (40,29 % nach χ₃)
  - Spezialzerlegung 11 vs 221 innerhalb χ₃ = 1
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

# ── Daten ────────────────────────────────────────────────────────────────────
R_VALS = np.array([11, 101, 191, 221, 311, 401])
COUNTS = np.array([765, 831, 786, 809, 809, 767], dtype=float)
N = COUNTS.sum()
MEAN = N / 6
DELTA_COUNT = COUNTS - MEAN
DELTA_PROB = COUNTS / N - 1 / 6

OUT_DIR = Path(__file__).resolve().parent
REPORT_TXT = OUT_DIR / "eabc_spectral_decomposition_report.txt"
REPORT_JSON = OUT_DIR / "eabc_spectral_decomposition.json"


# ── Charaktere ───────────────────────────────────────────────────────────────
OMEGA = np.exp(2j * np.pi / 3)

CHI3_TABLE = {1: 1.0, 2: OMEGA**2, 3: OMEGA, 4: 1.0, 5: OMEGA, 6: OMEGA**2}


def chi4(a: int) -> int:
    """Quadratischer Charakter mod 4."""
    a = int(a) % 4
    if a == 0:
        return 0
    return 1 if a == 1 else -1


def chi5(a: int) -> int:
    """Legendre-Symbol mod 5."""
    a = int(a) % 5
    if a == 0:
        return 0
    return 1 if a in (1, 4) else -1


def chi7_quad(a: int) -> int:
    """Quadratischer (Legendre-)Charakter mod 7."""
    a = int(a) % 7
    if a == 0:
        return 0
    return 1 if pow(a, 3, 7) == 1 else -1  # a^{(7-1)/2} mod 7


def chi3_mod7(a: int) -> complex:
    """Kubischer primitiver Charakter mod 7."""
    return CHI3_TABLE[int(a) % 7]


def eval_char(r: int, factors: dict[str, bool]) -> complex:
    """Produkt-Charakter χ₄^{e4} χ₅^{e5} χ₇^{e7q} χ₃^{e3} auf r."""
    val = 1.0 + 0.0j
    if factors.get("chi4"):
        val *= chi4(r)
    if factors.get("chi5"):
        val *= chi5(r)
    if factors.get("chi7q"):
        val *= chi7_quad(r)
    if factors.get("chi3"):
        val *= chi3_mod7(r)
    return val


# Alle relevanten Charaktere (Name → Faktoren)
CHARACTER_DEFS: dict[str, dict[str, bool]] = {
    "1": {},
    "χ₄": {"chi4": True},
    "χ₅": {"chi5": True},
    "χ₇": {"chi7q": True},
    "χ₃": {"chi3": True},
    "χ̄₃": {"chi3": True},  # wird konjugiert
    "χ₄χ₅": {"chi4": True, "chi5": True},
    "χ₄χ₇": {"chi4": True, "chi7q": True},
    "χ₅χ₇": {"chi5": True, "chi7q": True},
    "χ₄χ₅χ₇": {"chi4": True, "chi5": True, "chi7q": True},
    "χ₄χ₃": {"chi4": True, "chi3": True},
    "χ₄χ̄₃": {"chi4": True, "chi3": True},
    "χ₇χ₃": {"chi7q": True, "chi3": True},
    "χ₇χ̄₃": {"chi7q": True, "chi3": True},
    "χ₄χ₇χ₃": {"chi4": True, "chi7q": True, "chi3": True},
    "χ₄χ₇χ̄₃": {"chi4": True, "chi7q": True, "chi3": True},
}


def character_vector(name: str) -> np.ndarray:
    factors = CHARACTER_DEFS[name]
    conj_chi3 = name == "χ̄₃" or name.endswith("χ̄₃")
    vals = []
    for r in R_VALS:
        v = eval_char(r, factors)
        if conj_chi3:
            v = np.conj(v) if factors.get("chi3") else v
        vals.append(v)
    return np.array(vals, dtype=complex)


# ── Hilfsfunktionen ──────────────────────────────────────────────────────────
def inner_product(a: np.ndarray, b: np.ndarray) -> complex:
    return np.sum(a * np.conj(b)) / 6


def norm_sq(v: np.ndarray) -> float:
    return float(np.real(inner_product(v, v)))


def gram_schmidt(vectors: list[tuple[str, np.ndarray]]) -> list[tuple[str, np.ndarray]]:
    """Orthonormalisiere; verwerfe linear abhängige / Null-Vektoren."""
    basis: list[tuple[str, np.ndarray]] = []
    for name, v in vectors:
        w = v.astype(complex).copy()
        for _, e in basis:
            w -= inner_product(w, e) * e
        ns = norm_sq(w)
        if ns < 1e-14:
            continue
        basis.append((name, w / math.sqrt(ns)))
    return basis


def dft_mode(v: np.ndarray, k: int) -> complex:
    n = np.arange(6)
    return np.sum(v * np.exp(-2j * np.pi * k * n / 6)) / 6


@dataclass
class ProjectionRow:
    name: str
    coeff: complex
    magnitude: float
    variance_share_pct: float
    chi_values: list[complex]
    is_conjugate_pair: bool = False


def project_onto_basis(
    delta: np.ndarray, basis: list[tuple[str, np.ndarray]]
) -> tuple[list[ProjectionRow], np.ndarray, float]:
    """Projektion δ = Σ c_k e_k; Varianzanteil = |c_k|² / Var(δ)."""
    var_total = norm_sq(delta.astype(complex))
    rows: list[ProjectionRow] = []
    recon = np.zeros(6, dtype=complex)
    for name, e in basis:
        c = inner_product(delta.astype(complex), e)
        recon += c * e
        rows.append(
            ProjectionRow(
                name=name,
                coeff=c,
                magnitude=abs(c),
                variance_share_pct=100 * abs(c) ** 2 / var_total if var_total > 0 else 0.0,
                chi_values=[],
            )
        )
    return rows, recon, var_total


def main() -> None:
    lines: list[str] = []

    def p(s: str = "") -> None:
        print(s)
        lines.append(s)

    p("=" * 78)
    p("EABC Spektrale Zerlegung — Dirichlet-Charaktere mod 420")
    p("=" * 78)
    p(f"Kanäle r = {list(R_VALS)}")
    p(f"Counts   = {list(COUNTS.astype(int))},  N = {int(N)}")
    p(f"δ (prob) = [{', '.join(f'{x:+.5f}' for x in DELTA_PROB)}]")
    p()

    # ── 1. Charakterwerte auf Kanälen ────────────────────────────────────────
    p("─" * 78)
    p("1. Charakterwerte auf den 6 Kanälen")
    p("─" * 78)
    core_names = ["χ₄", "χ₅", "χ₇", "χ₃"]
    header = f"{'r':>5} {'r%4':>4} {'r%5':>4} {'r%7':>4}"
    for cn in core_names:
        header += f" {cn:>8}"
    p(header)
    for i, r in enumerate(R_VALS):
        row = f"{r:>5} {r % 4:>4} {r % 5:>4} {r % 7:>4}"
        for cn in core_names:
            v = character_vector(cn)[i]
            if np.isreal(v) and abs(np.imag(v)) < 1e-12:
                row += f" {int(np.real(v)):>8}"
            else:
                row += f" {v:>8.4f}"
        p(row)
    p()

    # Konstant-Check χ₅, χ mod 3
    chi5_vec = character_vector("χ₅")
    p(f"χ₅-Konstanz auf Kanälen: {np.allclose(chi5_vec, chi5_vec[0])}  (Werte: {chi5_vec.astype(int)})")
    chi3_mod3 = np.array([1 if r % 3 == 1 else -1 for r in R_VALS])
    p(f"Alle r ≡ {int(R_VALS[0] % 3)} (mod 3) → quadratischer χ mod 3 konstant")
    p()

    # ── 2. χ₃-Projektion (59,71 %) ─────────────────────────────────────────
    chi3 = character_vector("χ₃")
    chi3bar = np.conj(chi3)
    c_chi3 = inner_product(DELTA_PROB.astype(complex), chi3)
    c_chi3bar = inner_product(DELTA_PROB.astype(complex), chi3bar)
    var_total = norm_sq(DELTA_PROB.astype(complex))

    e0 = np.ones(6, dtype=complex) / math.sqrt(6)
    chi3_o = chi3 - inner_product(chi3, e0) * e0
    e1 = chi3_o / math.sqrt(norm_sq(chi3_o))
    chi3b_o = chi3bar - inner_product(chi3bar, e0) * e0 - inner_product(chi3bar, e1) * e1
    e2 = chi3b_o / math.sqrt(norm_sq(chi3b_o))

    var_chi3_pair = abs(inner_product(DELTA_PROB.astype(complex), e1)) ** 2 + abs(
        inner_product(DELTA_PROB.astype(complex), e2)
    ) ** 2
    pct_chi3 = 100 * var_chi3_pair / var_total

    F2 = dft_mode(DELTA_PROB, 2)

    p("─" * 78)
    p("2. χ₃ mod 7 und k=2-DFT (Theorem 1)")
    p("─" * 78)
    p(f"  ⟨δ, χ₃⟩     = {c_chi3:.10f}")
    p(f"  |⟨δ, χ₃⟩|    = {abs(c_chi3):.10f}")
    p(f"  F₂ (DFT k=2) = {F2:.10f}")
    p(f"  |F₂|         = {abs(F2):.10f}")
    p(f"  |⟨δ,χ₃⟩|/|F₂| = {abs(c_chi3) / abs(F2):.6f}")
    p(f"  Varianzanteil χ₃+χ̄₃: {pct_chi3:.2f}%  (Rest: {100 - pct_chi3:.2f}%)")
    p()

    # ── 3. Vollständige orthogonale Zerlegung ────────────────────────────────
    p("─" * 78)
    p("3. Orthogonale Zerlegung δ = Σ c_χ · ê_χ  (Gram-Schmidt)")
    p("─" * 78)

    # Basis-Reihenfolge: trivial, χ₃-Richtungen, dann quadratische Reste
    raw_vectors = [
        ("1", character_vector("1")),
        ("χ₃", chi3),
        ("χ̄₃", chi3bar),
        ("χ₄", character_vector("χ₄")),
        ("χ₇", character_vector("χ₇")),
        ("χ₄χ₇", character_vector("χ₄χ₇")),
        ("χ₄χ₃", character_vector("χ₄χ₃")),
        ("χ₇χ₃", character_vector("χ₇χ₃")),
        ("χ₄χ₇χ₃", character_vector("χ₄χ₇χ₃")),
        ("χ₅", character_vector("χ₅")),
        ("χ₄χ₅", character_vector("χ₄χ₅")),
        ("χ₅χ₇", character_vector("χ₅χ₇")),
        ("χ₄χ₅χ₇", character_vector("χ₄χ₅χ₇")),
    ]
    basis = gram_schmidt(raw_vectors)
    proj_rows, recon, _ = project_onto_basis(DELTA_PROB, basis)
    residual = DELTA_PROB.astype(complex) - recon
    var_residual = norm_sq(residual)

    p(f"  Dimension orthogonale Basis: {len(basis)}")
    p(f"  Var(δ) = {var_total:.10f}")
    p(f"  Σ|c_k|² = {sum(r.magnitude**2 for r in proj_rows):.10f}")
    p(f"  Residual-Varianz = {var_residual:.2e}")
    p()
    p(f"  {'Charakter':<12} {'|c|':>12} {'arg(c)°':>10} {'Var %':>10}")
    p(f"  {'-'*12} {'-'*12} {'-'*10} {'-'*10}")
    for row in proj_rows:
        p(
            f"  {row.name:<12} {row.magnitude:>12.8f} "
            f"{np.degrees(np.angle(row.coeff)):>10.2f} {row.variance_share_pct:>10.4f}"
        )
    p()

    # ── 4. Rest-Bias (40,29 %) in χ₄, χ₅, χ₇ ───────────────────────────────
    p("─" * 78)
    p("4. Rest-Bias-Aufschlüsselung (orthogonal zu χ₃+χ̄₃)")
    p("─" * 78)

    # Rest-Komponente: δ minus χ₃-Projektion
    delta_rest = DELTA_PROB.astype(complex) - (
        inner_product(DELTA_PROB.astype(complex), e1) * e1
        + inner_product(DELTA_PROB.astype(complex), e2) * e2
    )
    var_rest = norm_sq(delta_rest)

    quad_chars = ["χ₄", "χ₅", "χ₇", "χ₄χ₅", "χ₄χ₇", "χ₅χ₇", "χ₄χ₅χ₇"]
    p(f"  Var(Rest) = {var_rest:.10f}  ({100 * var_rest / var_total:.2f}% von Var(δ))")
    p()
    p(f"  {'Charakter':<12} {'|⟨Rest,χ⟩|':>14} {'Var-Anteil %':>14} {'Bemerkung':<20}")
    p(f"  {'-'*12} {'-'*14} {'-'*14} {'-'*20}")

    rest_breakdown = []
    for cn in quad_chars:
        chi_v = character_vector(cn)
        ns = norm_sq(chi_v)
        if ns < 1e-14:
            note = "konstant/0"
            mag = 0.0
            vpct = 0.0
        else:
            # Orthogonale Projektion auf χ-Richtung (nicht orthonormalisiert)
            c_raw = inner_product(delta_rest, chi_v)
            # Anteil an Rest-Varianz via orthogonalisiertem χ
            chi_orth = chi_v.copy()
            for _, e in [("e0", e0), ("χ₃", e1), ("χ̄₃", e2)]:
                chi_orth -= inner_product(chi_orth, e) * e
            ns_o = norm_sq(chi_orth)
            if ns_o < 1e-14:
                note = "∈ span(χ₃)"
                mag = 0.0
                vpct = 0.0
            else:
                e_chi = chi_orth / math.sqrt(ns_o)
                c_orth = inner_product(delta_rest, e_chi)
                mag = abs(c_orth)
                vpct = 100 * mag**2 / var_rest
                note = ""
        rest_breakdown.append((cn, mag, vpct, note))
        p(f"  {cn:<12} {mag:>14.8f} {vpct:>14.4f} {note:<20}")

    p()

    # ── 5. Spezial: 11 vs 221 (χ₃ = 1-Klasse) ───────────────────────────────
    p("─" * 78)
    p("5. Zerlegung innerhalb χ₃ = 1  (Kanäle 11 und 221)")
    p("─" * 78)
    idx_11 = int(np.where(R_VALS == 11)[0][0])
    idx_221 = int(np.where(R_VALS == 221)[0][0])
    d11 = DELTA_PROB[idx_11]
    d221 = DELTA_PROB[idx_221]
    diff = d11 - d221

    p(f"  δ₁₁  = {d11:+.6f}  (count={int(COUNTS[idx_11])})")
    p(f"  δ₂₂₁ = {d221:+.6f}  (count={int(COUNTS[idx_221])})")
    p(f"  δ₁₁ − δ₂₂₁ = {diff:+.6f}")
    p()
    p(f"  {'Charakter':<12} {'χ(11)':>10} {'χ(221)':>10} {'χ(11)−χ(221)':>14}")
    p(f"  {'-'*12} {'-'*10} {'-'*10} {'-'*14}")
    for cn in ["χ₄", "χ₅", "χ₇", "χ₃"]:
        v11 = character_vector(cn)[idx_11]
        v221 = character_vector(cn)[idx_221]
        if np.isreal(v11):
            p(f"  {cn:<12} {np.real(v11):>10.0f} {np.real(v221):>10.0f} {np.real(v11 - v221):>14.0f}")
        else:
            p(f"  {cn:<12} {v11:>10.4f} {v221:>10.4f} {v11 - v221:>14.4f}")

    # 2-Punkt-Analyse: nur 11 und 221
    chi4_2 = np.array([chi4(11), chi4(221)], dtype=float)
    delta_2 = np.array([d11, d221])
    # Auf χ₄: c = (1/2) Σ δᵢ χᵢ, Var = 2*c² (für ±1 Muster)
    c_chi4_2 = np.dot(delta_2, chi4_2) / 2
    var_2 = np.dot(delta_2, delta_2) / 2  # Var auf 2-Punkt-Menge
    p()
    p(f"  2-Punkt-Projektion (11, 221):")
    p(f"    ⟨δ, χ₄⟩_{'{2pt}'} = {c_chi4_2:+.8f}")
    p(f"    Erklärte Varianz |c|² / Var₂ = {100 * c_chi4_2**2 / var_2:.2f}%")
    p(f"    χ₅, χ₇ unterscheiden 11 und 221 nicht (beide ≡ 1 mod 5, ≡ 4 mod 7)")
    p()

    # ── 6. Direkte (nicht-ortho) Projektionstabelle ──────────────────────────
    p("─" * 78)
    p("6. Direkte Charakter-Projektionen c_χ = (1/6)⟨δ, χ⟩")
    p("─" * 78)
    all_names = [
        "1", "χ₃", "χ̄₃", "χ₄", "χ₅", "χ₇",
        "χ₄χ₅", "χ₄χ₇", "χ₅χ₇", "χ₄χ₅χ₇",
        "χ₄χ₃", "χ₇χ₃", "χ₄χ₇χ₃",
    ]
    direct_rows = []
    p(f"  {'Charakter':<14} {'|c|':>12} {'Var %':>10}")
    p(f"  {'-'*14} {'-'*12} {'-'*10}")
    for cn in all_names:
        chi_v = character_vector(cn)
        if cn == "χ̄₃":
            chi_v = np.conj(character_vector("χ₃"))
        c = inner_product(DELTA_PROB.astype(complex), chi_v)
        # Varianzbeitrag nur sinnvoll für orthogonal; hier |c|² als Rohgröße
        vpct = 100 * abs(c) ** 2 / var_total
        direct_rows.append({"name": cn, "coeff": complex(c), "magnitude": abs(c), "var_pct_raw": vpct})
        p(f"  {cn:<14} {abs(c):>12.8f} {vpct:>10.4f}")

    p()
    p("─" * 78)
    p("ZUSAMMENFASSUNG")
    p("─" * 78)
    p(f"  χ₃+χ̄₃ erklärt {pct_chi3:.2f}% der Varianz von δ")
    p(f"  Rest-Bias ({100 - pct_chi3:.2f}%) liegt in quadratischen Charakteren")
    p(f"  χ₅ ist auf allen 6 Kanälen trivial (=1) → kein Beitrag")
    p(f"  11 vs 221: Differenz wird allein durch χ₄ (mod 4) getragen")
    p(f"  |⟨δ,χ₃⟩| = |F₂| exakt (k=2-DFT ≡ kubischer Charakter mod 7)")
    p()

    # JSON export
    export = {
        "channels": R_VALS.tolist(),
        "counts": COUNTS.astype(int).tolist(),
        "delta_prob": DELTA_PROB.tolist(),
        "var_total": var_total,
        "chi3_variance_pct": pct_chi3,
        "rest_variance_pct": 100 - pct_chi3,
        "F2": complex(F2),
        "c_chi3": complex(c_chi3),
        "orthogonal_projections": [
            {
                "name": r.name,
                "coeff_real": r.coeff.real,
                "coeff_imag": r.coeff.imag,
                "magnitude": r.magnitude,
                "variance_pct": r.variance_share_pct,
            }
            for r in proj_rows
        ],
        "rest_breakdown": [
            {"character": cn, "magnitude": mag, "variance_pct": vpct, "note": note}
            for cn, mag, vpct, note in rest_breakdown
        ],
        "direct_projections": direct_rows,
        "diff_11_221": {
            "delta_11": float(d11),
            "delta_221": float(d221),
            "diff": float(diff),
            "chi4_projection_2pt": float(c_chi4_2),
        },
    }
    REPORT_JSON.write_text(json.dumps(export, indent=2, default=str), encoding="utf-8")
    REPORT_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    p(f"Bericht gespeichert: {REPORT_TXT}")
    p(f"JSON gespeichert:    {REPORT_JSON}")


if __name__ == "__main__":
    main()
