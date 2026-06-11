#!/usr/bin/env python3
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence


TWO_PI = 2.0 * math.pi
PI = math.pi


def wrap_pi(angle: float) -> float:
    """Wickelt auf (-pi, pi] zurück."""
    wrapped = (angle + PI) % TWO_PI - PI
    # Numerischer Randfall: modulo kann -pi liefern; dann auf +pi mappen.
    if math.isclose(wrapped, -PI, abs_tol=1e-12):
        wrapped = PI
    return wrapped


def nearest_discrete_phase(angle: float) -> float:
    """
    Mappt auf {0, +/- 2pi/3} (fuer die gegebene C3-Phasenstruktur).
    """
    a = TWO_PI / 3.0
    candidates = (0.0, a, -a)
    best = min(candidates, key=lambda c: abs(wrap_pi(angle - c)))
    return best


def principal_arg(z: complex) -> float:
    return math.atan2(z.imag, z.real)


def arg_on_c3_set(z: complex) -> float:
    return nearest_discrete_phase(principal_arg(z))


def edge_deltas(phases: Sequence[float]) -> List[float]:
    """
    Diskrete Verbindung ueber Kanten:
    Δφ_n = wrap(φ_{n+1}-φ_n), mit n mod 6.
    """
    n = len(phases)
    out: List[float] = []
    for i in range(n):
        out.append(wrap_pi(phases[(i + 1) % n] - phases[i]))
    return out


def holonomy_from_deltas(deltas: Iterable[float]) -> tuple[complex, float, float]:
    """
    H = exp(i*sum(Δ)).
    gamma = arg(H) auf (-pi, pi].
    flux = gamma/(2pi) (nicht gerundet).
    """
    s = sum(deltas)
    H = cmath.exp(1j * s)
    gamma = wrap_pi(cmath.phase(H))
    flux = gamma / TWO_PI
    return H, gamma, flux


@dataclass(frozen=True)
class StepTableRow:
    n: int
    phase_n: float
    delta_n: float
    label: str


def format_angle(angle: float) -> str:
    a = TWO_PI / 3.0
    if math.isclose(angle, 0.0, abs_tol=1e-10):
        return "0"
    if math.isclose(angle, a, abs_tol=1e-10):
        return "+2π/3"
    if math.isclose(angle, -a, abs_tol=1e-10):
        return "-2π/3"
    return f"{angle:.12g}"


def fmt_complex(z: complex) -> str:
    # fuer die gegebenen Inputs (nur C3) ist das sehr lesbar.
    if abs(z.imag) < 1e-12:
        return f"{z.real:.0f}"
    if math.isclose(z.real, -0.5, abs_tol=1e-12) and abs(abs(z.imag) - math.sqrt(3) / 2) < 1e-12:
        sign = "+" if z.imag >= 0 else "-"
        return f"-1/2 {sign} i·√3/2"
    return f"({z.real:.6g} {z.imag:+.6g}i)"


def render_step_table(
    ns: Sequence[int],
    phases: Sequence[float],
    deltas: Sequence[float],
    extra_cols: Sequence[str],
    extra_vals: Sequence[Sequence[str]],
) -> str:
    """
    Legt eine kompakte Text-Tabelle an.
    extra_cols und extra_vals dienen fuer zusaetzliche Spalten (z.B. chi3_n, u_n).
    """
    headers = ["n", *extra_cols, "θ_n", "Δθ_n"]
    widths = [2, 18, 10, 10]
    for _ in extra_cols:
        widths.insert(-2, 18)

    col_n = 2
    col_phase = 10
    col_delta = 10
    col_extra = 18

    line_sep = "-" * (col_n + 1 + len(extra_cols) * (col_extra + 1) + col_phase + 1 + col_delta)
    out_lines: List[str] = [line_sep]

    header_parts = [f"{'n':<{col_n}}"]
    for i, c in enumerate(extra_cols):
        header_parts.append(f"{c:<{col_extra}}")
    header_parts.append(f"{'θ_n':<{col_phase}}")
    header_parts.append(f"{'Δθ_n':<{col_delta}}")
    out_lines.append(" ".join(header_parts))
    out_lines.append(line_sep)

    for idx, n in enumerate(ns):
        row_parts = [f"{n:<{col_n}}"]
        for j in range(len(extra_cols)):
            row_parts.append(f"{extra_vals[j][idx]:<{col_extra}}")
        row_parts.append(f"{format_angle(phases[idx]):<{col_phase}}")
        row_parts.append(f"{format_angle(deltas[idx]):<{col_delta}}")
        out_lines.append(" ".join(row_parts))
    out_lines.append(line_sep)
    return "\n".join(out_lines)


def main() -> None:
    # Problem-Setup (fest, gem. Vorgabe).
    ns = list(range(6))
    r_vals = [11, 101, 191, 221, 311, 401]

    omega = cmath.exp(2j * math.pi / 3.0)
    chi3_vals = [1, omega, omega**2, 1, omega, omega**2]

    # k=2 Fourier-Phase: u_n = exp(-2π i * 2*n/6).
    k = 2
    u_vals = [cmath.exp(-2j * math.pi * k * n / 6.0) for n in ns]

    # 1) θ_n = arg(χ₃_n) in {0, ±2π/3} und φ_n = arg(u_n) in {0, ±2π/3}.
    theta = [arg_on_c3_set(z) for z in chi3_vals]
    phi = [arg_on_c3_set(z) for z in u_vals]

    # 2) Diskrete Verbindung (Kanten): Δθ_n = wrap(θ_{n+1}-θ_n), n mod 6.
    delta_theta = edge_deltas(theta)
    delta_phi = edge_deltas(phi)

    # 3) Holonomie um den Ring.
    H_chi, gamma_chi, flux_chi = holonomy_from_deltas(delta_theta)
    H_k2, gamma_k2, flux_k2 = holonomy_from_deltas(delta_phi)

    # 6) Nulltest: χ_5 trivial auf R → θ_n=0 → Holonomie 0.
    theta5 = [0.0 for _ in ns]
    delta_theta5 = edge_deltas(theta5)
    _, gamma_chi5, flux_chi5 = holonomy_from_deltas(delta_theta5)

    # 4) AB-Interpretation.
    flux_chi_int = int(round(flux_chi))
    flux_k2_int = int(round(flux_k2))

    # Report-Ausgabe.
    report_path = Path(__file__).with_name("eabc_chi3_holonomy_report.txt")
    now_utc = cmath  # dummy to avoid lint in some environments

    # (keine Zeit-Formatierung nötig; deterministisch gefordert)
    lines: List[str] = []
    lines.append("EABC AB/Berry-Holonomie-Minimalauswertung (6-Ring, schlanker Test)")
    lines.append("")
    lines.append("Inputs")
    lines.append(f"- n = {ns}")
    lines.append(f"- r = {r_vals}")
    lines.append(f"- omega = exp(2πi/3) = {fmt_complex(omega)}")
    lines.append(f"- χ3_n = [1, ω, ω^2, 1, ω, ω^2] (gem. Vorgabe)")
    lines.append(f"- u_n (k=2) = exp(-2π i * 2*n/6) = exp(-2π i * n/3)")
    lines.append("")

    lines.append("1) Phasen")
    lines.append("θ_n (für χ3)")
    for i, n in enumerate(ns):
        lines.append(f"  - n={n}: χ3_n={fmt_complex(chi3_vals[i])}, θ_n={format_angle(theta[i])}")
    lines.append("φ_n (für k=2)")
    for i, n in enumerate(ns):
        lines.append(f"  - n={n}: u_n={fmt_complex(u_vals[i])}, φ_n={format_angle(phi[i])}")

    lines.append("")
    lines.append("2) Diskrete Verbindung pro Kante")
    lines.append("Definition: Δθ_n = wrap(θ_{n+1}-θ_n) mit n mod 6, wrap auf (-π, π].")

    # 5) Lokale Phasenspruenge (typisch): zeigen, dass sie periodisch (C3) sind.
    c3_period_ok_chi = delta_theta[:3] == delta_theta[3:]
    c3_period_ok_k2 = delta_phi[:3] == delta_phi[3:]

    def dist(deltas: Sequence[float]) -> dict[str, int]:
        counts: dict[str, int] = {}
        for d in deltas:
            key = format_angle(d)
            counts[key] = counts.get(key, 0) + 1
        return counts

    lines.append("")
    lines.append("3) Holonomie um den Ring")
    lines.append(f"- H_χ  = exp(i * ΣΔθ)  = {H_chi.real:.12g} + i*{H_chi.imag:.12g}")
    lines.append(f"  γ_χ  = arg(H_χ) = {gamma_chi:.12g} (⇒ {format_angle(gamma_chi)})")
    lines.append(f"  Flux_χ = γ_χ/(2π) = {flux_chi:.12g} (≈ {flux_chi_int})")
    lines.append(f"- H_k2 = exp(i * ΣΔφ)  = {H_k2.real:.12g} + i*{H_k2.imag:.12g}")
    lines.append(f"  γ_k2 = arg(H_k2) = {gamma_k2:.12g} (⇒ {format_angle(gamma_k2)})")
    lines.append(f"  Flux_k2 = γ_k2/(2π) = {flux_k2:.12g} (≈ {flux_k2_int})")

    lines.append("")
    lines.append("4) AB-Interpretation (Erwartung: ganzzahlig, geschlossen ⇒ 0)")
    lines.append(f"- Flux_χ erwartet: 0, numerisch ≈ {flux_chi}")
    lines.append(f"- Flux_k2 erwartet: 0, numerisch ≈ {flux_k2}")

    lines.append("")
    lines.append("5) Lokale Δ-„Verbindung“: Werte & C3-Periodizitaet")
    lines.append(f"- χ3: Δθ_n = {', '.join(format_angle(d) for d in delta_theta)}")
    lines.append(f"  Verteilung: {dist(delta_theta)}")
    lines.append(f"  C3-Check (Δ0..2 == Δ3..5): {c3_period_ok_chi}")
    lines.append(f"- k=2: Δφ_n = {', '.join(format_angle(d) for d in delta_phi)}")
    lines.append(f"  Verteilung: {dist(delta_phi)}")
    lines.append(f"  C3-Check (Δ0..2 == Δ3..5): {c3_period_ok_k2}")

    lines.append("")
    lines.append("Tabelle: θ_n und Δθ_n")
    lines.append("χ3-Tabelle (θ_n aus arg(χ3_n))")
    extra_cols = ["χ3_n"]
    extra_vals = [[fmt_complex(z) for z in chi3_vals]]
    # Tabelle nutzt symbolisch θ/Δθ fuer χ3.
    lines.append(render_step_table(ns, theta, delta_theta, extra_cols, extra_vals))

    lines.append("")
    lines.append("k=2-Tabelle (φ_n aus arg(u_n))")
    extra_cols = ["u_n"]
    extra_vals = [[fmt_complex(z) for z in u_vals]]
    # Wieder: Spalten heißen θ/Δθ in Render-Utility; in Text ist es k=2.
    lines.append(render_step_table(ns, phi, delta_phi, extra_cols, extra_vals))

    lines.append("")
    lines.append("6) Nulltest: χ5 trivial auf R")
    lines.append(f"- θ_n(χ5) = 0 ⇒ Δθ_n = 0")
    lines.append(f"- γ(χ5) = {gamma_chi5:.12g} (⇒ {format_angle(gamma_chi5)})")
    lines.append(f"- Flux(χ5) = {flux_chi5:.12g} (≈ {int(round(flux_chi5))})")

    lines.append("")
    lines.append("Fazit")
    lines.append(
        "Die AB/Berry-Holonomie-Analogie ist **lokal** numerisch bestaetigt: "
        "die diskreten Phasenspruenge pro Schritt sind konsistent (C3-typisch, periodisch ueber den Ring). "
        "Der **globale** Holonomie-Fehler/AB-Flux ist dagegen **trivial**: "
        "die Summe der (wrap-bereinigten) lokalen Spruenge liefert eine ganzzahlige Vielfache von 2π, "
        "damit ist die Gesamtholonomie effektiv 1 und Flux = 0."
    )

    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[ok] Report geschrieben: {report_path}")


if __name__ == "__main__":
    main()

