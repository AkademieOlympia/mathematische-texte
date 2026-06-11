#!/usr/bin/env python3
"""
EABC Superalgebra-Fit: {Q, Q̄} ≈ Σ α_k Π_k in verschiedenen Π-Basen.

Struktureller Analogie-Test (keine physikalische SUSY).
Erweitert eabc_superalgebra_test.py um lineare Hüllen-Fits und Rank-Analyse.
"""

from __future__ import annotations

import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.linalg import eigvals, matrix_rank, norm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from Stiefel import build_transition_matrix, channel_mod420, dirac_operator, quadruplets  # noqa: E402
from eabc_superalgebra_test import (  # noqa: E402
    N_CH,
    anticommutator,
    build_candidates,
    cyclic_shift,
    frob,
    projector_mod12,
    test_fermionic,
)

OUT_DIR = Path(__file__).resolve().parent
REPORT = OUT_DIR / "eabc_superalgebra_fit_report.txt"

TOL_FIT = 1e-10
TOL_ALPHA = 1e-6


def vec(M: np.ndarray) -> np.ndarray:
    return M.ravel()


def build_pi_bases(T: np.ndarray, L: np.ndarray) -> dict[str, dict[str, np.ndarray]]:
    S = cyclic_shift()
    I = np.eye(N_CH, dtype=complex)

    p_shift = {
        "I": I,
        "S": S,
        "S^2": np.linalg.matrix_power(S, 2),
        "S^3": np.linalg.matrix_power(S, 3),
        "S^4": np.linalg.matrix_power(S, 4),
        "S^5": np.linalg.matrix_power(S, 5),
    }

    p_mod12 = {
        "Π_E": projector_mod12(1),
        "Π_A": projector_mod12(5),
        "Π_B": projector_mod12(7),
        "Π_C": projector_mod12(11),
    }

    p_T = {
        "I": I,
        "T": T.astype(complex),
        "T^T": T.T.astype(complex),
        "L": L.astype(complex),
        "L^T": L.T.astype(complex),
    }

    p_combined: dict[str, np.ndarray] = {}
    for d in (p_shift, p_mod12, p_T):
        for name, op in d.items():
            if name not in p_combined:
                p_combined[name] = op

    return {
        "P_shift": p_shift,
        "P_mod12": p_mod12,
        "P_T": p_T,
        "P_combined": p_combined,
    }


@dataclass
class FitResult:
    basis_name: str
    rel_err: float
    r2: float
    alpha: np.ndarray
    pi_names: list[str]
    dominant_pi: str
    dominant_contrib: float
    matrix_rank: int
    in_span: bool
    cond: float
    alpha_str: str = ""


@dataclass
class CandidateFitSummary:
    name: str
    fermionic: bool
    qq_norm: float
    qbarqbar_norm: float
    a_norm: float
    best_basis: str
    best_rel_err: float
    best_r2: float
    best_dominant_pi: str
    best_alpha_str: str
    fits: dict[str, FitResult] = field(default_factory=dict)


def fit_in_basis(A: np.ndarray, basis: dict[str, np.ndarray], basis_name: str) -> FitResult:
    pi_names = list(basis.keys())
    Pi_list = [basis[n] for n in pi_names]
    b = vec(A)
    M = np.column_stack([vec(P) for P in Pi_list])

    alpha, _, rank, s = np.linalg.lstsq(M, b, rcond=None)
    A_fit = sum(alpha[k] * Pi_list[k] for k in range(len(pi_names)))

    ss_res = float(norm(b - M @ alpha) ** 2)
    ss_tot = float(norm(b) ** 2)
    r2 = 1.0 - ss_res / max(ss_tot, 1e-30)
    rel_err = frob(A - A_fit) / max(frob(A), 1e-15)

    contribs = [abs(alpha[k]) * frob(Pi_list[k]) for k in range(len(pi_names))]
    dom_idx = int(np.argmax(contribs))
    dominant_pi = pi_names[dom_idx]
    dominant_contrib = contribs[dom_idx]

    mat_rank = int(matrix_rank(M, tol=1e-10))
    in_span = rel_err < TOL_FIT

    cond = float(s[0] / s[-1]) if len(s) > 0 and s[-1] > 1e-15 else float("inf")

    # kompakte α-Darstellung: Top-Koeffizienten nach Beitrag
    order = np.argsort(-np.array(contribs))
    parts: list[str] = []
    for idx in order:
        if contribs[idx] < TOL_ALPHA * max(contribs[dom_idx], 1e-15):
            continue
        a = alpha[idx]
        if abs(a.imag) < 1e-8:
            parts.append(f"{pi_names[idx]}:{a.real:+.4f}")
        else:
            parts.append(f"{pi_names[idx]}:{a.real:+.4f}{a.imag:+.4f}j")
        if len(parts) >= 4:
            break

    return FitResult(
        basis_name=basis_name,
        rel_err=rel_err,
        r2=r2,
        alpha=alpha,
        pi_names=pi_names,
        dominant_pi=dominant_pi,
        dominant_contrib=dominant_contrib,
        matrix_rank=mat_rank,
        in_span=in_span,
        cond=cond,
        alpha_str=", ".join(parts) if parts else "(alle ≈ 0)",
    )


def analyze_stability(summaries: list[CandidateFitSummary], basis_name: str) -> dict:
    """Prüft konsistente dominante Π und α-Richtungen über Kandidaten."""
    dom_counts: dict[str, int] = {}
    rel_errs: list[float] = []
    alpha_dirs: list[np.ndarray] = []

    for s in summaries:
        fit = s.fits.get(basis_name)
        if fit is None:
            continue
        dom_counts[fit.dominant_pi] = dom_counts.get(fit.dominant_pi, 0) + 1
        rel_errs.append(fit.rel_err)
        norm_a = norm(fit.alpha)
        if norm_a > 1e-15:
            alpha_dirs.append(fit.alpha / norm_a)

    if not dom_counts:
        return {"stable": "nein", "dominant": "-", "consistency": 0.0}

    top_dom = max(dom_counts, key=dom_counts.get)
    consistency = dom_counts[top_dom] / len(summaries)

    # Winkel zwischen α-Richtungen (fermionische Teilmenge separat behandeln wir extern)
    angle_spread = 0.0
    if len(alpha_dirs) >= 2:
        angles = []
        for i in range(len(alpha_dirs)):
            for j in range(i + 1, len(alpha_dirs)):
                inner = abs(np.vdot(alpha_dirs[i], alpha_dirs[j]))
                angles.append(inner)
        angle_spread = float(np.mean(angles)) if angles else 0.0

    med_err = float(np.median(rel_errs)) if rel_errs else float("inf")

    if consistency >= 0.75 and med_err < 0.2:
        stable = "ja"
    elif consistency >= 0.5 or med_err < 0.5:
        stable = "teilweise"
    else:
        stable = "nein"

    return {
        "stable": stable,
        "dominant": top_dom,
        "consistency": consistency,
        "angle_spread": angle_spread,
        "median_err": med_err,
    }


def run_fits(T: np.ndarray) -> tuple[list[CandidateFitSummary], dict[str, dict[str, np.ndarray]]]:
    L = np.eye(N_CH) - T
    D = dirac_operator(T)
    candidates = build_candidates(T, L, D)
    bases = build_pi_bases(T, L)

    summaries: list[CandidateFitSummary] = []

    for cand in candidates:
        nqq, nqbqb, ferm = test_fermionic(cand.Q, cand.Qbar)
        A = anticommutator(cand.Q, cand.Qbar)
        a_norm = frob(A)

        fits: dict[str, FitResult] = {}
        for bname, basis in bases.items():
            fits[bname] = fit_in_basis(A, basis, bname)

        best_basis = min(fits, key=lambda k: fits[k].rel_err)
        best = fits[best_basis]

        summaries.append(CandidateFitSummary(
            name=cand.name,
            fermionic=ferm,
            qq_norm=nqq,
            qbarqbar_norm=nqbqb,
            a_norm=a_norm,
            best_basis=best_basis,
            best_rel_err=best.rel_err,
            best_r2=best.r2,
            best_dominant_pi=best.dominant_pi,
            best_alpha_str=best.alpha_str,
            fits=fits,
        ))

    return summaries, bases


def write_report(
    T: np.ndarray,
    eigT: np.ndarray,
    gap: float,
    summaries: list[CandidateFitSummary],
    bases: dict[str, dict[str, np.ndarray]],
    elapsed: float,
) -> None:
    lines: list[str] = []
    lines.append("=" * 78)
    lines.append("EABC Superalgebra-Fit: {Q,Q̄} ≈ Σ α_k Π_k")
    lines.append("=" * 78)
    lines.append("")
    lines.append("EPISTEMOLOGISCHER HINWEIS")
    lines.append("-" * 40)
    lines.append(
        "Struktureller Analogie-Test auf H = C^6. Kein Anspruch auf physikalische "
        "Supersymmetrie oder Super-Poincaré-Einbettung."
    )
    lines.append("")

    lines.append("1. Π-BASEN")
    lines.append("-" * 40)
    for bname, basis in bases.items():
        lines.append(f"  {bname} ({len(basis)} Operatoren): {', '.join(basis.keys())}")
    lines.append("")

    lines.append("2. FIT-METHODE")
    lines.append("-" * 40)
    lines.append("  A = {Q, Q̄}")
    lines.append("  vec(A) = M @ α,  Spalten von M = vec(Π_k)")
    lines.append("  Least Squares (komplex), rel.Fehler = ||A - Σα_k Π_k||_F / ||A||_F")
    lines.append("  R² = 1 - ||Residuum||² / ||A||²_F")
    lines.append("  Dominantes Π_k: max_k |α_k| · ||Π_k||_F")
    lines.append("")

    lines.append("3. ERGEBNISSE — ALLE KANDIDATEN (bester Fit über alle Basen)")
    lines.append("-" * 40)
    hdr = f"{'Kandidat':<26} {'ferm?':>5} {'rel.Fehler':>12} {'R²':>8} {'Basis':>12} {'dom. Π':>8}"
    lines.append(hdr)
    sorted_all = sorted(summaries, key=lambda s: s.best_rel_err)
    for s in sorted_all:
        lines.append(
            f"{s.name:<26} {'JA' if s.fermionic else 'NEIN':>5} "
            f"{s.best_rel_err:12.4e} {s.best_r2:8.6f} {s.best_basis:>12} {s.best_dominant_pi:>8}"
        )
        lines.append(f"    α: {s.best_alpha_str}")
    lines.append("")

    fermionic = [s for s in summaries if s.fermionic]
    lines.append("4. FIT NUR FERMIONISCHE KANDIDATEN ({Q,Q}≈0, {Q̄,Q̄}≈0)")
    lines.append("-" * 40)
    if not fermionic:
        lines.append("  Keine fermionischen Kandidaten.")
    else:
        lines.append(hdr)
        for s in sorted(fermionic, key=lambda x: x.best_rel_err):
            lines.append(
                f"{s.name:<26} {'JA':>5} "
                f"{s.best_rel_err:12.4e} {s.best_r2:8.6f} {s.best_basis:>12} {s.best_dominant_pi:>8}"
            )
            lines.append(f"    α: {s.best_alpha_str}")
    lines.append("")

    lines.append("5. DETAIL — FIT PRO Π-BASIS (rel.Fehler)")
    lines.append("-" * 40)
    basis_names = list(bases.keys())
    header = f"{'Kandidat':<26}" + "".join(f"{b:>14}" for b in basis_names)
    lines.append(header)
    for s in sorted_all:
        row = f"{s.name:<26}"
        for b in basis_names:
            row += f"{s.fits[b].rel_err:14.4e}"
        lines.append(row)
    lines.append("")

    lines.append("6. RANK-ANALYSE: LIEGT A IN DER LINEAREN HÜLLE?")
    lines.append("-" * 40)
    lines.append(
        f"{'Kandidat':<26} {'Basis':>12} {'rank(M)':>8} {'cond':>12} {'in Hülle?':>10} {'rel.Fehler':>12}"
    )
    for s in sorted_all:
        for b in basis_names:
            f = s.fits[b]
            in_h = "JA" if f.in_span else "NEIN"
            cond_s = f"{f.cond:.2e}" if np.isfinite(f.cond) else "inf"
            lines.append(
                f"{s.name:<26} {b:>12} {f.matrix_rank:8d} {cond_s:>12} {in_h:>10} {f.rel_err:12.4e}"
            )
    lines.append("")
    lines.append("  Basis-Dimensionen und Ränge der Design-Matrizen M:")
    for bname, basis in bases.items():
        Pi_list = list(basis.values())
        M = np.column_stack([vec(P) for P in Pi_list])
        lines.append(f"    {bname}: {len(basis)} Spalten, rank(M) = {matrix_rank(M, tol=1e-10)}")
    lines.append("")

    lines.append("7. VERGLEICH DER BASEN")
    lines.append("-" * 40)
    for b in basis_names:
        errs = [s.fits[b].rel_err for s in summaries]
        lines.append(
            f"  {b}: min={min(errs):.4e}, median={np.median(errs):.4e}, "
            f"max={max(errs):.4e}, mean={np.mean(errs):.4e}"
        )
    best_basis_global = min(
        basis_names,
        key=lambda b: np.median([s.fits[b].rel_err for s in summaries]),
    )
    lines.append(f"  → Beste Basis (Median rel.Fehler): {best_basis_global}")
    lines.append("")

    lines.append("8. STABILE PARAMETRISIERUNG?")
    lines.append("-" * 40)
    for b in basis_names:
        stab = analyze_stability(summaries, b)
        lines.append(
            f"  {b}: {stab['stable'].upper()} "
            f"(dom. Π={stab['dominant']}, Konsistenz={stab['consistency']:.0%}, "
            f"median err={stab['median_err']:.4e})"
        )
    stab_combined = analyze_stability(summaries, "P_combined")
    stab_ferm = analyze_stability(fermionic, "P_combined") if fermionic else {"stable": "nein"}
    lines.append("")
    lines.append(f"  Gesamturteil (P_combined, alle Q): {stab_combined['stable']}")
    if fermionic:
        lines.append(f"  Gesamturteil (P_combined, nur fermionisch): {stab_ferm['stable']}")
    lines.append("")

    best = sorted_all[0]
    lines.append("9. BESTER KANDIDAT & ANALOGIE")
    lines.append("-" * 40)
    lines.append(f"  Bester Fit (global): {best.name}")
    lines.append(f"    rel.Fehler = {best.best_rel_err:.4e}, R² = {best.best_r2:.6f}")
    lines.append(f"    Basis = {best.best_basis}, dominantes Π = {best.best_dominant_pi}")
    lines.append(f"    α: {best.best_alpha_str}")
    lines.append("")

    if fermionic:
        best_f = min(fermionic, key=lambda s: s.best_rel_err)
        lines.append(f"  Bester fermionischer Fit: {best_f.name}")
        lines.append(f"    rel.Fehler = {best_f.best_rel_err:.4e}, Basis = {best_f.best_basis}")
        lines.append(f"    dominantes Π = {best_f.best_dominant_pi}, α: {best_f.best_alpha_str}")
    lines.append("")

    lines.append("10. SCHLUSSFOLGERUNG")
    lines.append("-" * 40)
    if best.best_rel_err < 0.05:
        analogie = (
            f"Starke lineare Analogie: {{Q,Q̄}} ≈ Σα_k Π_k für {best.name} "
            f"(rel.Fehler {best.best_rel_err:.4e}) in {best.best_basis}."
        )
    elif best.best_rel_err < 0.3:
        analogie = (
            f"Teilweise Analogie: {best.name} mit dominierendem {best.best_dominant_pi} "
            f"(rel.Fehler {best.best_rel_err:.4e}), aber keine universelle Schließung."
        )
    else:
        analogie = (
            "Keine der getesteten Konstruktionen liegt nahe genug in der Π-Hülle; "
            "Superalgebra-Schließung {Q,Q̄} ∝ 2ΣΠ_r bleibt numerisch unbestätigt."
        )
    lines.append(analogie)
    lines.append("")
    lines.append(f"λ₂ ≈ {eigT[1]}, Spektrallücke Δ ≈ {gap:.6f}")
    lines.append(f"Laufzeit: {elapsed:.2f} s")
    lines.append(f"Skript: {Path(__file__).name}")
    lines.append(f"Report: {REPORT.name}")

    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    t0 = time.perf_counter()
    print("EABC Superalgebra-Fit")
    print("Lade Vierlingsdaten (N=10^8) ...")

    starts = quadruplets(100_000_000)
    channels = np.array(
        [c for p in starts if (c := channel_mod420(int(p))) is not None],
        dtype=np.int64,
    )
    T = build_transition_matrix(channels)
    eigT = eigvals(T)
    eigT = eigT[np.argsort(-np.abs(eigT))]
    gap = float(1.0 - np.abs(eigT[1]))

    summaries, bases = run_fits(T)
    elapsed = time.perf_counter() - t0
    write_report(T, eigT, gap, summaries, bases, elapsed)

    best = min(summaries, key=lambda s: s.best_rel_err)
    ferm = [s.name for s in summaries if s.fermionic]
    print(f"K = {len(channels)}, bester Fit: {best.name} ({best.best_rel_err:.4e})")
    print(f"Fermionische Kandidaten: {ferm}")
    print(f"Report: {REPORT}")


if __name__ == "__main__":
    main()
