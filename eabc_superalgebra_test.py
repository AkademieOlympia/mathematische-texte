#!/usr/bin/env python3
"""
EABC Superalgebra-Analogie-Test (strukturell, keine physikalische SUSY).

Testet auf H = C^6 (Kanäle mod 420) Kandidaten für fermionische Generatoren
Q, Q̄ und arithmetische Translationsoperatoren Π_r gegen die N=1
Super-Poincaré-Antikommutator-Struktur

    {Q, Q̄} ≟ 2 Σ^r Π_r ,   {Q,Q} = {Q̄,Q̄} = 0 .

Epistemologisch: strukturelle Analogie-Test, kein Anspruch auf physikalische
Supersymmetrie oder Super-Poincaré-Einbettung.
"""

from __future__ import annotations

import math
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from numpy.linalg import eigvals, norm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from Stiefel import (  # noqa: E402
    HALF_CHANNELS,
    INDEX,
    build_transition_matrix,
    channel_mod420,
    dirac_operator,
    quadruplets,
)

OUT_DIR = Path(__file__).resolve().parent
REPORT = OUT_DIR / "eabc_superalgebra_report.txt"

R_VALS = np.array(HALF_CHANNELS, dtype=int)
N_CH = 6
OMEGA = np.exp(2j * np.pi / 3)
CHI3_TABLE = {1: 1.0, 2: OMEGA**2, 3: OMEGA, 4: 1.0, 5: OMEGA, 6: OMEGA**2}
COUNTS = np.array([765, 831, 786, 809, 809, 767], dtype=float)
DELTA_PROB = COUNTS / COUNTS.sum() - 1.0 / 6


# ── Hilfsfunktionen ──────────────────────────────────────────────────────────

def anticommutator(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return A @ B + B @ A


def commutator(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return A @ B - B @ A


def frob(M: np.ndarray) -> float:
    return float(norm(M, "fro"))


def rel_err(M: np.ndarray, target: np.ndarray) -> float:
    denom = max(frob(target), 1e-15)
    return frob(M - target) / denom


def best_scale(M: np.ndarray, target: np.ndarray) -> complex:
    """Skalar c mit minimaler ||M - c·target||_F (komplexes least squares)."""
    denom = np.vdot(target.ravel(), target.ravel())
    if abs(denom) < 1e-15:
        return 0.0 + 0.0j
    return np.vdot(target.ravel(), M.ravel()) / denom


def fmt_matrix(M: np.ndarray, prec: int = 4) -> str:
    lines = []
    for row in M:
        parts = []
        for z in row:
            if abs(z.imag) < 10 ** (-prec - 1):
                parts.append(f"{z.real:>{prec+6}.{prec}f}")
            else:
                parts.append(f"{z.real:>{prec+3}.{prec}f}{z.imag:+.{prec}f}j")
        lines.append("  [" + ", ".join(parts) + "]")
    return "\n".join(lines)


def dft_vector(k: int) -> np.ndarray:
    n = np.arange(N_CH)
    return np.exp(2j * np.pi * k * n / N_CH) / math.sqrt(N_CH)


def projector(k: int) -> np.ndarray:
    psi = dft_vector(k)
    return np.outer(psi, psi.conj())


def cyclic_shift() -> np.ndarray:
    S = np.zeros((N_CH, N_CH), dtype=complex)
    for i in range(N_CH):
        S[(i + 1) % N_CH, i] = 1.0
    return S


def chi3_on_channels() -> np.ndarray:
    return np.array([CHI3_TABLE[int(r % 7)] for r in R_VALS], dtype=complex)


def mult_op(values: np.ndarray) -> np.ndarray:
    return np.diag(values)


def mod12_class(r: int) -> str:
    m = r % 12
    return {1: "E", 5: "A", 7: "B", 11: "C"}.get(m, f"?{m}")


def projector_mod12(cls: int) -> np.ndarray:
    mask = np.array([1.0 if (r % 12) == cls else 0.0 for r in R_VALS])
    if mask.sum() == 0:
        return np.zeros((N_CH, N_CH), dtype=complex)
    v = mask / math.sqrt(mask.sum())
    return np.outer(v, v.conj())


# Pauli-Matrizen (für Σ^r in der reduzierten 2×2-Analogie)
SIGMA_0 = np.eye(2, dtype=complex)
SIGMA_1 = np.array([[0, 1], [1, 0]], dtype=complex)
SIGMA_2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
SIGMA_3 = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMAS = [SIGMA_0, SIGMA_1, SIGMA_2, SIGMA_3]


@dataclass
class Candidate:
    name: str
    Q: np.ndarray
    Qbar: np.ndarray
    construction: str


@dataclass
class TestResult:
    name: str
    QQ_norm: float
    QbarQbar_norm: float
    fermionic: bool
    best_pi_name: str
    best_pi_scale: complex
    best_pi_rel_err: float
    target_P_shift10_scale: complex
    target_P_shift10_rel_err: float
    anticomm: np.ndarray
    sigma_match: dict[str, float]


def build_candidates(T: np.ndarray, L: np.ndarray, D: np.ndarray) -> list[Candidate]:
    S = cyclic_shift()
    chi = chi3_on_channels()
    chibar = np.conj(chi)
    psi1 = dft_vector(1)
    psi5 = dft_vector(5)

    def ketbra(u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """|u><v| (not necessarily Hermitian). Assumes normalized u,v for nilpotency checks."""
        return np.outer(u, v.conj())

    cands: list[Candidate] = []

    # (a) FFT-Moden k=1 / k=5 als Projektoren
    cands.append(Candidate(
        "FFT_proj_k1_k5",
        projector(1),
        projector(5),
        "Q = P_{k=1}, Q̄ = P_{k=5} (DFT-Projektoren)",
    ))

    # (a0) nilpotente Außenprodukt-Konstruktion zwischen k=1 und k=5
    # Da DFT-Moden orthonormal sind: <ψ5|ψ1> = 0 ⇒ Q^2 = Q̄^2 = 0.
    cands.append(Candidate(
        "FFT_nilpotent_k1_to_k5",
        ketbra(psi1, psi5),
        ketbra(psi5, psi1),
        "Q = |ψ1><ψ5|, Q̄ = |ψ5><ψ1| (nilpotent wegen Orthogonalität)",
    ))

    # (a') FFT mit Shift-Kopplung (nilpotent-ähnlich)
    cands.append(Candidate(
        "FFT_shift_k1_k5",
        projector(1) @ S,
        projector(5) @ S.conj().T,
        "Q = P_1 S, Q̄ = P_5 S†",
    ))

    # (a'') Re/Im-Teile der k=1,k=5 Moden als „Majorana"-Paar
    M1 = np.outer(np.real(psi1), np.real(psi1)) + 1j * np.outer(np.imag(psi1), np.real(psi1))
    M5 = np.outer(np.real(psi5), np.real(psi5)) + 1j * np.outer(np.imag(psi5), np.real(psi5))
    cands.append(Candidate(
        "FFT_majorana_k1_k5",
        M1,
        M5,
        "Q, Q̄ aus Re/Im-Kopplung der k=1/k=5 DFT-Moden",
    ))

    # (b) Charakter-Multiplikationsoperatoren χ₃ / χ̄₃
    cands.append(Candidate(
        "chi3_mult",
        mult_op(chi),
        mult_op(chibar),
        "Q = M_{χ₃}, Q̄ = M_{χ̄₃} (Diagonal)",
    ))

    # (b0) nilpotente Außenprodukt-Konstruktion zwischen χ₃ und χ̄₃
    # (wenn χ̄₃ orthogonal zu χ₃ auf den 6 Punkten ist, dann Q^2=0 und Q̄^2=0)
    cands.append(Candidate(
        "chi3_nilpotent_kets",
        ketbra(chi, chibar),
        ketbra(chibar, chi),
        "Q = |χ₃><χ̄₃|, Q̄ = |χ̄₃><χ₃| (nilpotent-ähnlich)",
    ))

    # (b') χ₃ mit zyklischem Shift verflochten
    cands.append(Candidate(
        "chi3_shift",
        mult_op(chi) @ S,
        mult_op(chibar) @ S.conj().T,
        "Q = M_{χ₃} S, Q̄ = M_{χ̄₃} S†",
    ))

    # (b'') Off-diagonal χ₃-Kopplung benachbarter Kanäle
    B = np.zeros((N_CH, N_CH), dtype=complex)
    for i in range(N_CH):
        j = (i + 1) % N_CH
        B[i, j] = chi[i]
        B[j, i] = chibar[j]
    cands.append(Candidate(
        "chi3_neighbor",
        B,
        B.conj().T,
        "Q_{ij} = χ₃(i)δ_{j,i+1} + χ̄₃(j)δ_{i,j+1} (zyklisch)",
    ))

    # (c) Clifford-Struktur auf 6-Ring: antisymmetrisch + Shift
    J = (S - S.conj().T) / 2  # Imaginärteil / antisymmetrisch
    K = (S + S.conj().T) / 2  # symmetrischer Shift
    cands.append(Candidate(
        "clifford_JS",
        J,
        1j * K,
        "Q = (S-S†)/2, Q̄ = i(S+S†)/2",
    ))

    # (c') Gamma-Analogie: {γ,γ}=2 auf k=2-Unterraum
    P2 = projector(2)
    P4 = projector(4)
    gamma = P2 - P4  # Hermitische Differenz k=2 vs k=4
    cands.append(Candidate(
        "clifford_gamma_k2k4",
        gamma @ S,
        1j * gamma,
        "Q = (P₂-P₄)S, Q̄ = i(P₂-P₄)",
    ))

    # (d) Aus Dirac-Operator D = [[0,L],[L†,0]] in Stiefel.py
    L_T = L.T
    cands.append(Candidate(
        "dirac_L_and_transpose",
        L.astype(complex),
        L_T.astype(complex),
        "Q = L = I-T, Q̄ = Lᵀ (Direkt aus dem Dirac-Block D abgeleitet)",
    ))

    cands.append(Candidate(
        "dirac_coupling_normalized",
        L.astype(complex) / (frob(L) + 1e-15),
        L_T.astype(complex) / (frob(L) + 1e-15),
        "Q = L/‖L‖, Q̄ = Lᵀ/‖L‖ (Dirac-Block-Kopplung, normalisiert)",
    ))

    # (d') Nilpotenter Anteil aus Schur-Zerlegung von L
    U, s, Vh = np.linalg.svd(L, full_matrices=False)
    if s[0] > 1e-12:
        Qn = U[:, :1] @ Vh[:1, :]  # Rang-1
        cands.append(Candidate(
            "dirac_svd_rank1",
            Qn,
            Qn.conj().T,
            "Q = u₁v₁† (Rang-1-SVD von L=I-T; i.A. nicht fermionisch nilpotent)",
        ))

    # (e) δ-Vektor / Besetzungs-Bias als „Ordnungsparameter"
    delta = DELTA_PROB.astype(complex)
    cands.append(Candidate(
        "delta_outer",
        np.outer(delta, np.ones(N_CH)),
        np.outer(np.ones(N_CH), delta.conj()),
        "Q = δ·1ᵀ, Q̄ = 1·δ† (Rang-1 aus δ)",
    ))

    return cands


def build_pi_operators(T: np.ndarray, L: np.ndarray) -> dict[str, np.ndarray]:
    S = cyclic_shift()
    ops: dict[str, np.ndarray] = {}

    ops["S (zykl. Shift)"] = S
    ops["S + S†"] = S + S.conj().T
    ops["i(S - S†)"] = 1j * (S - S.conj().T)
    ops["Re(S)"] = (S + S.conj().T) / 2
    ops["Im(S)"] = (S - S.conj().T) / (2j)

    for r in range(1, 4):
        Sr = np.linalg.matrix_power(S, r)
        ops[f"S^{r} + S^{-r}"] = Sr + Sr.conj().T

    ops["L = I - T"] = L.astype(complex)
    ops["Lᵀ"] = L.T.astype(complex)
    ops["T - I"] = (T - np.eye(N_CH)).astype(complex)
    ops["T (Übergangsmatrix)"] = T.astype(complex)

    # Translation-Generatoren als Shift-Potenzen (Π_r = S^r)
    # Wir benutzen später typischerweise r=0..3 als Analogon zu σ^μ-Komponenten.
    for r in range(6):
        ops[f"Π_shift^{r} (S^{r})"] = np.linalg.matrix_power(S, r)

    # mod-12 Projektionen E,A,B,C
    for cls, label in [(1, "E"), (5, "A"), (7, "B"), (11, "C")]:
        ops[f"Π_{label} (mod 12)"] = projector_mod12(cls)

    # Kombination Σ^r Π_r mit Pauli-Struktur auf 2D-Unterraum {A,C}
    # Nur A(5) und C(11) besetzt auf R
    Pi_A = projector_mod12(5)
    Pi_C = projector_mod12(11)
    ops["Σ⁰(Π_A+Π_C)"] = Pi_A + Pi_C
    ops["Σ¹(Π_C-Π_A)"] = Pi_C - Pi_A
    ops["Σ³(Π_C+Π_A)"] = Pi_C + Pi_A  # = Σ⁰ hier

    return ops


def test_fermionic(Q: np.ndarray, Qbar: np.ndarray) -> tuple[float, float, bool]:
    nqq = frob(anticommutator(Q, Q))
    nqbqb = frob(anticommutator(Qbar, Qbar))
    tol = 1e-8
    return nqq, nqbqb, (nqq < tol and nqbqb < tol)


def test_against_pi(
    A: np.ndarray,
    pi_ops: dict[str, np.ndarray],
) -> tuple[str, complex, float]:
    best_name = ""
    best_c = 0.0 + 0.0j
    best_err = float("inf")
    for name, Pi in pi_ops.items():
        c = best_scale(A, Pi)
        err = rel_err(A, c * Pi)
        if err < best_err:
            best_err = err
            best_c = c
            best_name = name
    return best_name, best_c, best_err


def test_sigma_structure(
    A: np.ndarray,
    pi_ops: dict[str, np.ndarray],
) -> dict[str, float]:
    """Prüft 2 Σ^r Π_r mit r=0..3 auf dem {A,C}-Unterraum."""
    Pi_A = pi_ops["Π_A (mod 12)"]
    Pi_C = pi_ops["Π_C (mod 12)"]
    bases = [
        ("Σ⁰·I", Pi_A + Pi_C),
        ("Σ¹·(Π_C-Π_A)", Pi_C - Pi_A),
        ("Σ²·i(Π_C-Π_A)", 1j * (Pi_C - Pi_A)),
        ("Σ³·(Π_A+Π_C)", Pi_A + Pi_C),
    ]
    out: dict[str, float] = {}
    for label, target in bases:
        c = best_scale(A, 2 * target)
        out[label] = rel_err(A, c * 2 * target)
    return out


def run_tests(T: np.ndarray) -> tuple[list[TestResult], dict[str, np.ndarray], np.ndarray]:
    L = np.eye(N_CH) - T
    D = dirac_operator(T)
    pi_ops = build_pi_operators(T, L)
    candidates = build_candidates(T, L, D)

    results: list[TestResult] = []
    for cand in candidates:
        nqq, nqbqb, ferm = test_fermionic(cand.Q, cand.Qbar)
        A = anticommutator(cand.Q, cand.Qbar)
        best_name, best_c, best_err = test_against_pi(A, pi_ops)

        # zusätzlich: Zieloperator P_target = 2 Σ_{r=0..3} Π_shift^r
        S = cyclic_shift()
        P_target = 2 * sum(np.linalg.matrix_power(S, r) for r in range(4))
        c_target = best_scale(A, P_target)
        err_target = rel_err(A, c_target * P_target)

        sigma = test_sigma_structure(A, pi_ops)
        results.append(TestResult(
            name=cand.name,
            QQ_norm=nqq,
            QbarQbar_norm=nqbqb,
            fermionic=ferm,
            best_pi_name=best_name,
            best_pi_scale=best_c,
            best_pi_rel_err=best_err,
            target_P_shift10_scale=c_target,
            target_P_shift10_rel_err=err_target,
            anticomm=A,
            sigma_match=sigma,
        ))

    return results, pi_ops, D


def write_report(
    T: np.ndarray,
    eigT: np.ndarray,
    gap: float,
    results: list[TestResult],
    candidates: list[Candidate],
    pi_ops: dict[str, np.ndarray],
    D: np.ndarray,
    elapsed: float,
) -> None:
    lines: list[str] = []
    lines.append("=" * 78)
    lines.append("EABC Superalgebra-Analogie-Test")
    lines.append("=" * 78)
    lines.append("")
    lines.append("EPISTEMOLOGISCHER HINWEIS")
    lines.append("-" * 40)
    lines.append(
        "Dies ist ein struktureller Analogie-Test zwischen N=1 Super-Poincaré-"
    )
    lines.append(
        "Antikommutator-Relationen und EABC-Kanaloperatoren auf H = C^6."
    )
    lines.append(
        "Kein Anspruch auf physikalische Supersymmetrie oder echte Super-Poincaré-"
    )
    lines.append("Einbettung.")
    lines.append("")

    lines.append("1. VORHANDENER DIRAC-OPERATOR (Stiefel.py)")
    lines.append("-" * 40)
    lines.append("Stiefel.py definiert:")
    lines.append("  L = I - T   (6×6, T = zeilenstochastische Übergangsmatrix)")
    lines.append("  D = [[0, L], [Lᵀ, 0]]   (12×12, block-off-diagonal / „Dirac“)")
    lines.append("Es gibt keinen separaten Majorana-Operator; der 12×12-Block-D")
    lines.append("ist der einzige Dirac-ähnliche Operator. Kleinste Eigenwerte von D")
    lines.append("entsprechen √(eigenvalues of L Lᵀ)).")
    lines.append("")

    lines.append("2. EABC-DATEN")
    lines.append("-" * 40)
    lines.append(f"Kanäle R = {list(R_VALS)}")
    lines.append(f"Counts   = {list(COUNTS.astype(int))}")
    lines.append(f"δ (prob) = {np.round(DELTA_PROB, 5).tolist()}")
    lines.append("")
    lines.append("Übergangsmatrix T:")
    lines.append(fmt_matrix(T.astype(complex)))
    lines.append("")
    lines.append(f"Eigenwerte(T) (|λ| absteigend): {eigT.tolist()}")
    lines.append(f"λ₂ = {eigT[1]}")
    lines.append(f"Spektrallücke Δ = 1 - |λ₂| = {gap:.6f}")
    lines.append("")
    lines.append("χ₃|_R = k=2 DFT-Modus (Theorem): identisch numerisch.")
    chi = chi3_on_channels()
    psi2 = dft_vector(2) * math.sqrt(N_CH)
    lines.append(f"  max|χ₃ - ψ₂| = {np.max(np.abs(chi - psi2)):.2e}")
    lines.append("")

    lines.append("3. DEFINITIONEN DER GETESTETEN OPERATOREN")
    lines.append("-" * 40)
    for cand in candidates:
        lines.append(f"  [{cand.name}]")
        lines.append(f"    {cand.construction}")
    lines.append("")
    lines.append("  Π_r-Kandidaten (Translations-/Impuls-Analogien):")
    for name in pi_ops:
        lines.append(f"    - {name}")
    lines.append("")

    lines.append("4. FERMIONISCHE SCHLIESSUNG {Q,Q} = {Q̄,Q̄} = 0")
    lines.append("-" * 40)
    lines.append(f"{'Kandidat':<22} {'‖{Q,Q}‖':>12} {'‖{Q̄,Q̄}‖':>12} {'fermionisch?':>12}")
    for r in results:
        lines.append(
            f"{r.name:<22} {r.QQ_norm:12.4e} {r.QbarQbar_norm:12.4e} "
            f"{'JA' if r.fermionic else 'NEIN':>12}"
        )
    lines.append("")

    lines.append("5. ANTIKOMMUTATOR {Q,Q̄} vs Π-OPERATOREN")
    lines.append("-" * 40)
    lines.append(f"{'Kandidat':<22} {'beste Π':<24} {'rel.Fehler':>12} {'rel.Fehler vs 2ΣS^r':>18}")
    sorted_results = sorted(results, key=lambda r: r.best_pi_rel_err)
    for r in sorted_results:
        lines.append(
            f"{r.name:<22} {r.best_pi_name:<24} "
            f"{r.best_pi_scale.real:+.4f}{r.best_pi_scale.imag:+.4f}j "
            f"{r.best_pi_rel_err:12.4e} {r.target_P_shift10_rel_err:18.4e}"
        )
    lines.append("")

    lines.append("6. PAULI-ÄHNLICHE STRUKTUR 2 Σ^r Π_r")
    lines.append("-" * 40)
    for r in sorted_results[:5]:
        lines.append(f"  [{r.name}]")
        for label, err in sorted(r.sigma_match.items(), key=lambda x: x[1]):
            lines.append(f"    {label}: rel.Fehler = {err:.4e}")
    lines.append("")

    # wähle bester nicht-trivialer Kandidat:
    # (Q,Q̄ sollen nicht beide fast-Null sein)
    def is_trivial(cand: Candidate) -> bool:
        return (frob(cand.Q) < 1e-12) or (frob(cand.Qbar) < 1e-12)

    nontrivial = [r for r in sorted_results if not is_trivial(next(c for c in candidates if c.name == r.name))]
    if nontrivial:
        fermionic_nontrivial = [r for r in nontrivial if r.fermionic]
    else:
        fermionic_nontrivial = []

    # Für den eigentlichen Superalgebra-Schluss ist fermionische Schließung entscheidend.
    # Daher wählen wir den besten fermionischen Kandidaten bzgl. Abweichung gegen
    # den Zieloperator P_target = 2 Σ_{r=0..3} S^r.
    best = (
        min(fermionic_nontrivial, key=lambda r: r.target_P_shift10_rel_err)
        if fermionic_nontrivial
        else nontrivial[0] if nontrivial else sorted_results[0]
    )
    best_cand = next(c for c in candidates if c.name == best.name)
    lines.append("7. BESTER (NICHTTRIVIALER) KANDIDAT — NUMERISCHE {Q,Q̄}-MATRIX")
    lines.append("-" * 40)
    lines.append(f"Kandidat: {best.name}")
    lines.append(f"Konstruktion: {best_cand.construction}")
    lines.append(f"Beste Π-Analogie: {best.best_pi_name}, c = {best.best_pi_scale}")
    lines.append(f"Relativer Fehler: {best.best_pi_rel_err:.4e}")
    lines.append(f"Rel. Fehler gegen 2Σ_(r=0..3) S^r: {best.target_P_shift10_rel_err:.4e}")
    lines.append("")
    lines.append("{Q, Q̄} =")
    lines.append(fmt_matrix(best.anticomm))
    lines.append("")
    lines.append(f"c · {best.best_pi_name} =")
    Pi_best = pi_ops[best.best_pi_name]
    lines.append(fmt_matrix(best.best_pi_scale * Pi_best))
    lines.append("")

    lines.append("8. SCHLUSSFOLGERUNG (ANALOGIE-TEST)")
    lines.append("-" * 40)
    target_err = best.target_P_shift10_rel_err
    if best.fermionic and target_err < 0.05:
        closure = "JA (teilweise)"
        closure_detail = (
            "Ein fermionischer Kandidat erfüllt nilpotente Schließung "
            f"und {{Q,Q̄}} ≈ c·(2ΣᵣΠᵣ) bis rel. Fehler {target_err:.4f}."
        )
    elif target_err < 0.15:
        closure = "TEILWEISE"
        closure_detail = (
            f"{{Q,Q̄}} ist näherungsweise proportional zu (2ΣᵣΠᵣ) "
            f"(rel. Fehler {target_err:.4f}), aber fermionische Schließung ist "
            "nur für einzelne Konstruktionen zufriedenstellend."
        )
    else:
        closure = "NEIN"
        closure_detail = (
            "Keine getestete Konstruktion erfüllt simultan fermionische "
            "Schließung und {{Q,Q̄}} ≈ c·(2ΣᵣΠᵣ) mit kleinem Fehler."
        )

    lines.append(f"Schließt die Algebra? {closure}")
    lines.append(closure_detail)
    lines.append("")
    lines.append(
        f"Beste Analogie zu {{Q,Q̄}} ~ P: {best.name} → {best.best_pi_name}"
    )
    lines.append("")
    lines.append("Verbindung zu bekannten EABC-Ergebnissen:")
    lines.append(f"  λ₂ ≈ {eigT[1].real:.4f} ± {abs(eigT[1].imag):.4f}i")
    lines.append(f"  Δ ≈ {gap:.4f}  (kein Markov-Gedächtnis)")
    lines.append(f"  χ₃-Projektion: |⟨δ,χ₃⟩| = |F₂| = {abs(np.vdot(DELTA_PROB, chi)):.6f}")
    lines.append(f"  Varianzanteil χ₃+χ̄₃ ≈ 59.7 %")
    lines.append("")
    lines.append(f"Laufzeit: {elapsed:.2f} s")
    lines.append(f"Skript: {Path(__file__).name}")
    lines.append(f"Report: {REPORT.name}")

    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    t0 = time.perf_counter()
    print("EABC Superalgebra-Analogie-Test")
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

    L = np.eye(N_CH) - T
    D = dirac_operator(T)
    candidates = build_candidates(T, L, D)
    results, pi_ops, D = run_tests(T)

    elapsed = time.perf_counter() - t0
    write_report(T, eigT, gap, results, candidates, pi_ops, D, elapsed)

    print(f"K = {len(channels)}, λ₂ = {eigT[1]}, Δ = {gap:.6f}")
    best = min(results, key=lambda r: r.best_pi_rel_err)
    print(f"Bester Kandidat: {best.name} → {best.best_pi_name} (err={best.best_pi_rel_err:.4e})")
    ferm = [r.name for r in results if r.fermionic]
    print(f"Fermionische Kandidaten: {ferm if ferm else 'keiner'}")
    print(f"Report geschrieben: {REPORT}")


if __name__ == "__main__":
    main()
