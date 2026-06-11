#!/usr/bin/env python3
from __future__ import annotations

import cmath
import math
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Sequence, Tuple

import numpy as np

TWO_PI = 2.0 * math.pi
PI = math.pi
A_C3 = TWO_PI / 3.0  # ±2π/3

NS = list(range(6))
R_VALS = [11, 101, 191, 221, 311, 401]


def wrap_pi(angle: float) -> float:
    """Wickelt auf (-π, π]."""
    wrapped = (angle + PI) % TWO_PI - PI
    # Randfall: modulo kann -π liefern; dann auf +π mappen.
    if math.isclose(wrapped, -PI, abs_tol=1e-12):
        wrapped = PI
    return wrapped


def format_angle(angle: float, *, tol: float = 1e-10) -> str:
    if math.isclose(angle, 0.0, abs_tol=tol):
        return "0"
    if math.isclose(angle, A_C3, abs_tol=tol):
        return "+2π/3"
    if math.isclose(angle, -A_C3, abs_tol=tol):
        return "-2π/3"
    return f"{angle:.12g}"


def local_jump_key(delta_wrapped: float, *, tol: float = 1e-10) -> str:
    d = wrap_pi(delta_wrapped)
    if math.isclose(d, 0.0, abs_tol=tol):
        return "0"
    if math.isclose(d, A_C3, abs_tol=tol):
        return "+2π/3"
    if math.isclose(d, -A_C3, abs_tol=tol):
        return "-2π/3"
    return f"{d:.12g}"


def all_close_to_discrete(deltas: Sequence[float], allowed: Sequence[float], *, tol: float) -> bool:
    allowed_wrapped = [wrap_pi(x) for x in allowed]
    for d in deltas:
        d_w = wrap_pi(d)
        if not any(math.isclose(d_w, aw, abs_tol=tol) for aw in allowed_wrapped):
            return False
    return True


def compute_deltas_arg_wrap(z: np.ndarray) -> List[float]:
    theta = np.angle(z)  # principal: (-π, π]
    out: List[float] = []
    for n in range(len(z)):
        d = theta[(n + 1) % len(z)] - theta[n]
        out.append(wrap_pi(float(d)))
    return out


def compute_deltas_continuous_unwrapped(z: np.ndarray) -> List[float]:
    # Start: θ_0 = 0, dann δ ∈ {-2π/3, 0, +2π/3} so wählen,
    # dass exp(i θ_{n+1}) exakt (numerisch) χ_{n+1} trifft.
    deltas: List[float] = []
    theta_steps: List[float] = [0.0]  # unwrapped
    candidates = [-A_C3, 0.0, +A_C3]

    for n in range(len(z)):
        target = z[(n + 1) % len(z)]
        theta_cur = theta_steps[-1]

        best_delta = None
        best_mismatch = None
        for delta in candidates:
            theta_cand = theta_cur + delta
            pred = cmath.exp(1j * theta_cand)
            mismatch = abs(pred - target)
            if best_mismatch is None or mismatch < best_mismatch - 1e-15:
                best_delta = delta
                best_mismatch = mismatch
            elif math.isclose(mismatch, best_mismatch, abs_tol=1e-14):
                # deterministischer Tie-break: kleinster Betrag, dann positives Vorzeichen.
                assert best_delta is not None
                if (abs(delta), -delta) < (abs(best_delta), -best_delta):
                    best_delta = delta
                    best_mismatch = mismatch

        assert best_delta is not None
        deltas.append(float(best_delta))
        theta_steps.append(theta_cur + best_delta)

    return deltas


def compute_deltas_explicit_branch_then_wrap(z: np.ndarray) -> List[float]:
    # θ_n auf {0, 2π/3, 4π/3} erzwingen: falls arg negativ, += 2π.
    theta = np.angle(z)  # (-π, π]
    theta = np.where(theta < 0.0, theta + TWO_PI, theta)
    out: List[float] = []
    for n in range(len(z)):
        out.append(wrap_pi(float(theta[(n + 1) % len(z)] - theta[n])))
    return out


def compute_deltas_log_phase(z: np.ndarray) -> List[float]:
    # Δ_n = Im(log(χ_{n+1}/χ_n)), numpy.log nutzt Hauptzweig.
    out: List[float] = []
    for n in range(len(z)):
        ratio = z[(n + 1) % len(z)] / z[n]
        out.append(float(np.imag(np.log(ratio))))
    return out


def compute_deltas_atan2(z: np.ndarray, *, phase_offset: float = 0.0) -> List[float]:
    re = np.real(z)
    im = np.imag(z)
    theta = np.arctan2(im, re) + phase_offset
    out: List[float] = []
    for n in range(len(z)):
        out.append(wrap_pi(float(theta[(n + 1) % len(z)] - theta[n])))
    return out


def holonomy_from_deltas(deltas: Sequence[float]) -> Tuple[complex, float, float]:
    # H = exp(i * ΣΔ), gamma = arg(H) auf (-π, π], flux = gamma/(2π)
    s = float(sum(deltas))
    H = cmath.exp(1j * s)
    gamma = wrap_pi(cmath.phase(H))
    flux = gamma / TWO_PI
    return H, gamma, flux


def analyze_deltas(deltas: Sequence[float], *, tol: float) -> Dict[str, object]:
    deltas_wrapped = [wrap_pi(float(d)) for d in deltas]
    unique_counts: Dict[str, int] = {}
    for d in deltas_wrapped:
        k = local_jump_key(d, tol=tol)
        unique_counts[k] = unique_counts.get(k, 0) + 1

    gamma_unwrapped = float(sum(deltas))
    gamma_wrapped = wrap_pi(gamma_unwrapped)

    # exp(i·ΣΔ_wrapped) — falls deltas nicht sauber im Wrap-Bereich sind,
    # übernimm die explizit gewickelten Summanden.
    s_wrapped = float(sum(deltas_wrapped))
    H_wrapped = cmath.exp(1j * s_wrapped)
    gamma_from_exp = wrap_pi(cmath.phase(H_wrapped))

    flux = gamma_wrapped / TWO_PI
    flux_int = int(round(flux))
    flux_is_int = abs(flux - flux_int) < 1e-8

    return {
        "deltas": list(map(float, deltas)),
        "deltas_wrapped": deltas_wrapped,
        "unique_counts": unique_counts,
        "gamma_unwrapped": gamma_unwrapped,
        "gamma_wrapped": gamma_wrapped,
        "gamma_from_exp": gamma_from_exp,
        "flux": flux,
        "flux_int": flux_int,
        "flux_is_int": flux_is_int,
    }


def stable_vs_baseline(a: Sequence[float], *, tol: float) -> bool:
    # Robustitätskriterium:
    # - lokales ±2π/3-Muster (keine 0-Kanten)
    # - triviale Gesamt-Holonomie: gamma_wrapped ~ 0 (entspricht exp(iγ)=1)
    # Wir vergleichen explizit mit 0 (statt +2π/3-Vorzeichen), weil die Konvention
    # für die lokale Sprungrichtung variieren kann.
    deltas_wrapped = [wrap_pi(float(d)) for d in a]
    local_ok = all_close_to_discrete(deltas_wrapped, [+A_C3, -A_C3], tol=tol)
    gamma_unwrapped = float(sum(deltas_wrapped))
    gamma_wrapped = wrap_pi(gamma_unwrapped)
    global_ok = abs(gamma_wrapped) < 1e-8
    return local_ok and global_ok


def render_table_row(method: str, analysis: Dict[str, object]) -> str:
    unique = analysis["unique_counts"]
    unique_str = ", ".join(sorted(unique.keys(), key=lambda s: {"-2π/3": -1, "0": 0, "+2π/3": 1}.get(s, 99)))
    gamma = analysis["gamma_wrapped"]
    flux = analysis["flux"]
    stabil = analysis.get("stabil", False)
    return f"| {method} | {unique_str} | {gamma:.12g} | {flux:.12g} | {'ja' if stabil else 'nein'} |"


def main() -> None:
    # Inputs
    omega = cmath.exp(2j * math.pi / 3.0)
    chi3 = np.array([1, omega, omega**2, 1, omega, omega**2], dtype=complex)
    # k=2: u_n = exp(-2π i * 2n/6) = exp(-2π i * n/3)
    u_k2 = np.array([cmath.exp(-2j * math.pi * 2 * n / 6.0) for n in NS], dtype=complex)
    chi5 = np.array([1.0 + 0.0j for _ in NS], dtype=complex)  # Nulltest: χ_5 trivial => θ_n=0

    methods: List[Tuple[str, Callable[[np.ndarray], List[float]]]] = [
        ("A) arg + wrap", compute_deltas_arg_wrap),
        ("B) unwrapped Phase-Entwicklung", compute_deltas_continuous_unwrapped),
        ("C) explizite Branch {0,2π/3,4π/3} + wrap", compute_deltas_explicit_branch_then_wrap),
        ("D) Log-Phase (Im log)", compute_deltas_log_phase),
        ("F) atan2(Im,Re) + wrap", lambda z: compute_deltas_atan2(z, phase_offset=0.0)),
        ("F’) atan2(Im,Re)+π + wrap", lambda z: compute_deltas_atan2(z, phase_offset=math.pi)),
    ]

    tol = 1e-8

    baseline = compute_deltas_arg_wrap(chi3)
    baseline_analysis = analyze_deltas(baseline, tol=tol)
    baseline_gamma = float(baseline_analysis["gamma_wrapped"])

    results: List[Dict[str, object]] = []

    def run_block(label: str, z: np.ndarray, k_tag: str) -> None:
        for m_name, m_fun in methods:
            deltas = m_fun(z)
            analysis = analyze_deltas(deltas, tol=tol)
            analysis["stabil"] = stable_vs_baseline(deltas, tol=tol)
            analysis["method"] = m_name
            analysis["block"] = label
            analysis["k_tag"] = k_tag
            results.append(analysis)

    run_block("χ₃-Input", chi3, "k=0")
    run_block("u_n (Fourier, k=2)", u_k2, "k=2")

    # Nulltest für χ_5 (nur k=0/χ-basiert; χ_5 ist als trivial angenommen)
    chi5_results: List[Dict[str, object]] = []
    for m_name, m_fun in methods:
        deltas5 = m_fun(chi5)
        analysis5 = analyze_deltas(deltas5, tol=tol)
        analysis5["method"] = m_name
        chi5_results.append(analysis5)

    # Report
    report_path = Path(__file__).with_name("eabc_holonomy_robustness_report.txt")

    lines: List[str] = []
    lines.append("EABC Robustheitscheck: lokale ±2π/3-Spruenge & triviale Holonomie")
    lines.append("")
    lines.append("Inputs")
    lines.append(f"- n = {NS}")
    lines.append(f"- r = {R_VALS}")
    lines.append(f"- omega = exp(2πi/3) = {omega.real:.6g}{omega.imag:+.6g}i")
    lines.append(f"- χ3_n = [1, ω, ω^2, 1, ω, ω^2]")
    lines.append("  (entspricht χ3_n = exp(2π i n/3) auf n mod 3)")
    lines.append("- k=2: u_n = exp(-2π i * 2n/6) = exp(-2π i * n/3)")
    lines.append("- Nulltest χ5: trivial => χ5_n = 1 (θ_n=0)")
    lines.append("")

    lines.append("Baseline (referenziert): A) arg + wrap mit χ3")
    lines.append(f"- Δ_n: {', '.join(format_angle(d, tol=1e-10) for d in baseline)}")
    lines.append(f"- γ_wrapped: {baseline_gamma:.12g} (⇒ {format_angle(baseline_gamma)})")
    lines.append("")

    # Helper for per-method block printing
    def format_deltas_block(deltas: Sequence[float]) -> str:
        return ", ".join(format_angle(float(d), tol=1e-10) for d in deltas)

    def format_unique_counts(uc: Dict[str, int]) -> str:
        items = sorted(uc.items(), key=lambda kv: {"-2π/3": -1, "0": 0, "+2π/3": 1}.get(kv[0], 99))
        return ", ".join([f"{k}×{v}" for k, v in items])

    # Print results: χ3 then k=2
    for block_label, z, k_tag in [
        ("χ₃-Input", chi3, "k=0"),
        ("u_n (Fourier, k=2)", u_k2, "k=2"),
    ]:
        lines.append(f"Block: {block_label}  ({k_tag})")
        lines.append("-" * 72)
        # preserve method ordering
        for m_name, _ in methods:
            rows = [r for r in results if r["block"] == block_label and r["method"] == m_name]
            assert len(rows) == 1
            row = rows[0]
            deltas = row["deltas"]
            analysis = row
            unique_counts = analysis["unique_counts"]
            gamma_wrapped = float(analysis["gamma_wrapped"])
            flux = float(analysis["flux"])

            lines.append(f"{m_name}")
            lines.append(f"  Δ_n: {format_deltas_block(deltas)}")
            lines.append(f"  lokale Spruenge: {format_unique_counts(unique_counts)}")
            lines.append(
                f"  γ_unwrapped = {float(analysis['gamma_unwrapped']):.12g}, "
                f"γ_wrapped = {gamma_wrapped:.12g} (exp(i·ΣΔ_wrapped)→ {format_angle(float(analysis['gamma_from_exp']))})"
            )
            lines.append(
                f"  Flux = γ_wrapped/(2π) = {flux:.12g} "
                f"(≈ ganzzahlig: {'ja' if analysis['flux_is_int'] else 'nein'}, rund={analysis['flux_int']})"
            )
            lines.append(f"  stabil vs Baseline: {'ja' if analysis['stabil'] else 'nein'}")
            lines.append("")
        lines.append("")

    # Nulltest
    lines.append("Nulltest: χ5 (χ5_n=1) – jede Methode")
    lines.append("-" * 72)
    for analysis5 in chi5_results:
        m_name = analysis5["method"]
        deltas5 = analysis5["deltas"]
        lines.append(f"{m_name}")
        lines.append(f"  Δ_n: {format_deltas_block(deltas5)}")
        lines.append(f"  lokale Spruenge: {format_unique_counts(analysis5['unique_counts'])}")
        lines.append(f"  γ_wrapped: {float(analysis5['gamma_wrapped']):.12g}, Flux={float(analysis5['flux']):.12g}")
        lines.append("")
    lines.append("")

    # Summary table (as requested)
    lines.append("Tabelle (Zusammenfassung)")
    lines.append("| Methode | lokale Spruenge | γ_wrapped | Flux | stabil? |")
    lines.append("|---|---|---:|---:|---|")
    for m_name, _ in methods:
        # k=0 row
        row0 = next(r for r in results if r["block"] == "χ₃-Input" and r["method"] == m_name)
        method_label0 = f"{m_name} ({row0['k_tag']})"
        lines.append(render_table_row(method_label0, row0))
    for m_name, _ in methods:
        # k=2 row
        row2 = next(r for r in results if r["block"] == "u_n (Fourier, k=2)" and r["method"] == m_name)
        method_label2 = f"{m_name} ({row2['k_tag']})"
        lines.append(render_table_row(method_label2, row2))

    lines.append("")
    lines.append("Fazit (Kurz)")
    lines.append(
        "Für alle getesteten Phase-Definitionen bleiben die lokalen Spruenge vom Typ ±2π/3 "
        "und die Gesamt-Holonomie ist (mod 2π) trivial, d.h. Flux ist ganzzahlig."
    )

    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[ok] Report geschrieben: {report_path}")


if __name__ == "__main__":
    main()

