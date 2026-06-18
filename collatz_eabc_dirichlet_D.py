#!/usr/bin/env python3
"""
Dirichlet-Erzeugerfunktion der EABC-Schalen-Anomalie D(n).

Kanonsiche Theorie: collatz_eabc_quaternion_mass_hypothese.md §13
(EABC-Spektralgeometrische Erzeugerhypothese).

D̂(s) = Σ_{n≥1} D(n) / n^s  mit skalarer Anomalie |D(n)| = ||I(μ_n) − I_ref(n)||.

I_ref (epistemisch dokumentiert):
  rolling      — bevorzugt (lokale geometrische Glättung, §12.4)
  mu_infinity  — bevorzugt (globales μ_∞-Proxy)

Ausführung:
    python3 collatz_eabc_dirichlet_D.py
    python3 collatz_eabc_dirichlet_D.py --max-n 500 --output collatz_eabc_dirichlet_D.json
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

from collatz_eabc_shell_defekt_test import (
    I_REF_VARIANTS,
    ShellInvariants,
    bernoulli_V_n,
    compute_shell_invariants,
    cumulative_baseline,
    defect_magnitude,
    mean_I_vector,
    omega_baseline,
    rolling_baseline,
    shell_defekt_report,
    tau_baseline,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_dirichlet_D.json"

# Bernoulli B_{2m} als exakte Brüche (m = 1..8)
BERNOULLI_EVEN: dict[int, Fraction] = {
    2: Fraction(1, 6),
    4: Fraction(-1, 30),
    6: Fraction(1, 42),
    8: Fraction(-1, 30),
    10: Fraction(5, 66),
    12: Fraction(-691, 2730),
    14: Fraction(7, 6),
    16: Fraction(-3617, 510),
}

DEFAULT_S_VALUES: tuple[float, ...] = (
    2.0,
    1.0,
    0.0,
    -1.0,
    -2.0,
    -4.0,
    -6.0,
    -8.0,
    -10.0,
)

PREFERRED_IREF: tuple[str, ...] = ("rolling", "mu_infinity")


@dataclass(frozen=True, slots=True)
class DTerm:
    n: int
    D: float
    is_prime: bool
    omega: int
    tau: int


def bernoulli_B(two_m: int) -> Fraction | None:
    """B_{2m} für gerades 2m in der Tabelle."""
    return BERNOULLI_EVEN.get(two_m)


def zeta_bernoulli_value(m: int) -> float | None:
    """
    ζ(1−2m) = −B_{2m}/(2m) für m ≥ 1.
    Beispiel: m=1 → ζ(−1) = −B_2/2 = −1/12.
    """
    if m < 1:
        return None
    two_m = 2 * m
    B = bernoulli_B(two_m)
    if B is None:
        return None
    return float(-B / two_m)


def collect_D_terms(
    max_n: int,
    iref: str,
    rolling_window: int = 5,
) -> list[DTerm]:
    """Berechnet skalare |D(n)| für n = 2..max_n mit gewähltem I_ref."""
    if iref not in I_REF_VARIANTS:
        raise ValueError(f"Unbekanntes I_ref: {iref!r}")

    invariants: list[ShellInvariants] = []
    for n in range(2, max_n + 1):
        inv = compute_shell_invariants(n)
        if inv is not None:
            invariants.append(inv)

    if not invariants:
        return []

    K_scale = max(statistics.mean(abs(r.K_trace) for r in invariants), 0.01)
    I_global = mean_I_vector(invariants)

    terms: list[DTerm] = []
    for i, inv in enumerate(invariants):
        if iref == "rolling":
            I_ref = rolling_baseline(invariants, i, window=rolling_window)
        elif iref == "cumulative":
            I_ref = cumulative_baseline(invariants, i)
        elif iref == "omega":
            I_ref = omega_baseline(invariants, inv)
        elif iref == "tau":
            I_ref = tau_baseline(invariants, inv)
        elif iref == "mu_infinity":
            I_ref = I_global
        else:
            raise ValueError(f"Unbekanntes I_ref: {iref!r}")

        d_mag = defect_magnitude(
            inv.H_n - I_ref[0],
            inv.chi_n - I_ref[1],
            inv.K_trace - I_ref[2],
            K_scale,
        )
        terms.append(
            DTerm(
                n=inv.n,
                D=d_mag,
                is_prime=inv.is_prime,
                omega=inv.omega,
                tau=inv.tau,
            )
        )
    return terms


def dirichlet_partial(terms: list[DTerm], s: float) -> float:
    """Partialsumme D̂_N(s) = Σ D(n) / n^s."""
    total = 0.0
    for t in terms:
        total += t.D / (t.n**s)
    return total


def zeta_partial(N: int, s: float, start_n: int = 2) -> float:
    """Σ_{n=start_n}^{N} 1/n^s — Referenz für Repackaging-Check."""
    return sum(1.0 / (n**s) for n in range(start_n, N + 1))


def arithmetic_partial(
    terms: list[DTerm],
    s: float,
    field: str,
) -> float:
    """Σ field(n) / n^s über die Schalen-indizes."""
    total = 0.0
    for t in terms:
        val = getattr(t, field)
        total += float(val) / (t.n**s)
    return total


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3 or len(xs) != len(ys):
        return None
    mx = statistics.mean(xs)
    my = statistics.mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mx) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - my) ** 2 for y in ys))
    if den_x < 1e-15 or den_y < 1e-15:
        return None
    return num / (den_x * den_y)


def growth_table(terms: list[DTerm], s_values: tuple[float, ...]) -> list[dict[str, Any]]:
    if not terms:
        return []
    N = terms[-1].n
    rows: list[dict[str, Any]] = []
    for s in s_values:
        D_hat = dirichlet_partial(terms, s)
        zeta_hat = zeta_partial(N, s, start_n=terms[0].n)
        omega_hat = arithmetic_partial(terms, s, "omega")
        tau_hat = arithmetic_partial(terms, s, "tau")
        ratio_zeta = D_hat / zeta_hat if abs(zeta_hat) > 1e-15 else None
        rows.append(
            {
                "s": s,
                "D_hat_N": round(D_hat, 6),
                "zeta_partial_N": round(zeta_hat, 6),
                "omega_dirichlet_partial": round(omega_hat, 6),
                "tau_dirichlet_partial": round(tau_hat, 6),
                "D_hat_over_zeta_partial": round(ratio_zeta, 6) if ratio_zeta is not None else None,
            }
        )
    return rows


def bernoulli_comparison(
    terms: list[DTerm],
    s_targets: tuple[float, ...] = (-2.0, -4.0, -1.0, -3.0),
) -> dict[str, Any]:
    """
    Explorativer Vergleich D̂_N(s) mit Bernoulli-Zahlen und ζ(1−2m)-Brücke.
    Keine Gleichheitsbehauptung — nur Skalen und Verhältnisse.
    """
    out: dict[str, Any] = {}
    for s in s_targets:
        D_hat = dirichlet_partial(terms, s)
        entry: dict[str, Any] = {
            "D_hat_N": round(D_hat, 6),
            "note": "explorativ; Partialsumme divergiert typischerweise bei s ≤ 0",
        }
        # s = −2m → Vergleich mit B_{2m}
        if s < 0 and s == int(s) and int(s) % 2 == 0:
            two_m = -int(s)
            m = two_m // 2
            B = bernoulli_B(two_m)
            zeta_odd = zeta_bernoulli_value(m) if two_m % 4 == 2 else None
            if B is not None:
                entry["B_2m"] = {"2m": two_m, "value": float(B), "fraction": str(B)}
                if abs(float(B)) > 1e-15:
                    entry["D_hat_over_B_2m"] = round(D_hat / float(B), 6)
            if zeta_odd is not None:
                entry["zeta_at_s_plus_1"] = {
                    "s_bernoulli": 1 - two_m,
                    "value": zeta_odd,
                    "note": f"ζ({1 - two_m}) = −B_{two_m}/{two_m}",
                }
                if abs(zeta_odd) > 1e-15:
                    entry["D_hat_over_zeta_bernoulli"] = round(D_hat / zeta_odd, 6)
            if two_m % 4 == 0:
                entry["zeta_trivial_zero"] = {
                    "s": s,
                    "value": 0.0,
                    "note": "ζ(s)=0 an geraden negativen s — Bernoulli-Brücke liegt bei ungeraden s",
                }
        key = str(int(s)) if s == int(s) else str(s)
        out[key] = entry
    return out


def first_terms_table(terms: list[DTerm], count: int = 12) -> list[dict[str, Any]]:
    return [
        {
            "n": t.n,
            "D": round(t.D, 6),
            "D_over_n2": round(t.D / (t.n**2), 8),
            "is_prime": t.is_prime,
            "omega": t.omega,
            "tau": t.tau,
        }
        for t in terms[:count]
    ]


def repackaging_assessment(
    terms: list[DTerm],
    growth: list[dict[str, Any]],
) -> dict[str, Any]:
    """Ehrliche Einordnung: ähnelt D̂ bekannter arithmetischer Reihen?"""
    if not terms:
        return {"verdict": "no_data"}

    D_vals = [t.D for t in terms]
    omega_vals = [float(t.omega) for t in terms]
    tau_vals = [float(t.tau) for t in terms]
    log_n = [math.log(t.n) for t in terms]

    rho_omega = pearson(D_vals, omega_vals)
    rho_tau = pearson(D_vals, tau_vals)
    rho_log = pearson(D_vals, log_n)

    s2_row = next((r for r in growth if r["s"] == 2.0), None)
    s0_row = next((r for r in growth if r["s"] == 0.0), None)

    ratios = [
        r["D_hat_over_zeta_partial"]
        for r in growth
        if r.get("D_hat_over_zeta_partial") is not None
    ]
    ratio_stable = (
        len(ratios) >= 2
        and max(ratios) > 0
        and (max(ratios) - min(ratios)) / max(abs(r) for r in ratios) < 0.05
    )

    looks_like_zeta = ratio_stable and s2_row is not None
    looks_like_arithmetic = (
        (rho_omega is not None and abs(rho_omega) > 0.5)
        or (rho_tau is not None and abs(rho_tau) > 0.5)
    )

    if looks_like_zeta and looks_like_arithmetic:
        verdict = (
            "D̂_N(s) skaliert annähernd wie ζ-Partialsummen; |D(n)| korreliert mit ω/τ — "
            "wahrscheinlich Repackaging bekannter arithmetischer Struktur, keine neue "
            "EABC-Erzeugerfunktion sichtbar."
        )
    elif looks_like_zeta:
        verdict = (
            "D̂_N(s) ≈ konstant · Σ 1/n^s über untersuchte s — eher ζ-Repackaging; "
            "EABC-spezifische Zusatzstruktur schwach."
        )
    elif looks_like_arithmetic:
        verdict = (
            "|D(n)| trägt starke ω-/τ-Korrelation — Anomalie reflektiert Teilerstruktur, "
            "nicht offensichtlich neue Dirichlet-Analytik."
        )
    else:
        verdict = (
            "Kein klares ζ- oder divisor-sum-Repackaging in den ersten N — weiterer "
            "Skalenvergleich und größeres N nötig; kein Bernoulli-Match behauptet."
        )

    return {
        "pearson_D_vs_omega": round(rho_omega, 6) if rho_omega is not None else None,
        "pearson_D_vs_tau": round(rho_tau, 6) if rho_tau is not None else None,
        "pearson_D_vs_log_n": round(rho_log, 6) if rho_log is not None else None,
        "D_hat_sum_D_at_s0": s0_row["D_hat_N"] if s0_row else None,
        "zeta_ratio_stable_across_s": ratio_stable,
        "verdict": verdict,
    }


def dirichlet_D_report(
    max_n: int = 200,
    rolling_window: int = 5,
    s_values: tuple[float, ...] = DEFAULT_S_VALUES,
    iref_keys: tuple[str, ...] = PREFERRED_IREF,
) -> dict[str, Any]:
    by_iref: dict[str, Any] = {}
    for iref in iref_keys:
        terms = collect_D_terms(max_n, iref, rolling_window=rolling_window)
        growth = growth_table(terms, s_values)
        by_iref[iref] = {
            "I_ref_note": I_REF_VARIANTS[iref],
            "term_count": len(terms),
            "first_terms": first_terms_table(terms),
            "growth_by_s": growth,
            "bernoulli_exploration": bernoulli_comparison(terms),
            "repackaging": repackaging_assessment(terms, growth),
        }

    # Kurzvergleich mit shell_defekt_report (Konsistenz)
    shell_snap = shell_defekt_report(max_n=min(max_n, 50), rolling_window=rolling_window)

    return {
        "meta": {
            "hypothesis_doc": "collatz_eabc_quaternion_mass_hypothese.md §13",
            "D_definition": "|D(n)| = ||I(μ_n) − I_ref(n)|| (euklidische Norm, K-skaliert)",
            "I_ref_preferred": list(iref_keys),
            "I_ref_all_variants": I_REF_VARIANTS,
            "max_n": max_n,
            "rolling_window": rolling_window,
            "s_values": list(s_values),
            "partial_sum_note": (
                "D̂_N(s) sind endliche Partialsummen; Konvergenz bei s≤1 nicht behauptet. "
                "Bernoulli-Vergleich bei s=−2,−4 ist explorativ (ζ trivial null dort)."
            ),
            "bernoulli_bridge": "ζ(1−2m) = −B_{2m}/(2m); vgl. collatz_eabc_bernoulli_uebersetzung.md",
            "shell_defekt_crosscheck_n50": {
                "rolling_ratio_prime_composite": shell_snap["prime_vs_composite"]["rolling"].get(
                    "ratio_mean_prime_over_composite"
                ),
            },
        },
        "by_I_ref": by_iref,
        "epistemic_status": (
            "Experiment — testet EABC-Spektralgeometrische Erzeugerhypothese (§13). "
            "Keine Gleichheit D̂(−2m) = B_{2m} behauptet."
        ),
    }


def run(
    max_n: int = 200,
    output: Path = DEFAULT_OUTPUT,
    rolling_window: int = 5,
) -> dict[str, Any]:
    report = dirichlet_D_report(max_n=max_n, rolling_window=rolling_window)
    report["output_path"] = str(output)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC Dirichlet-Erzeuger D̂(s) für Schalen-Defekt D(n)")
    parser.add_argument("--max-n", type=int, default=200, help="Obergrenze für n (Schalen)")
    parser.add_argument("--rolling-window", type=int, default=5)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(max_n=args.max_n, output=args.output, rolling_window=args.rolling_window)
    primary = report["by_I_ref"]["rolling"]
    b = primary["bernoulli_exploration"]
    print(
        f"collatz_eabc_dirichlet_D: N≤{args.max_n}, rolling I_ref, "
        f"D̂(-2)={b['-2']['D_hat_N']}, D̂(-4)={b['-4']['D_hat_N']} → {args.output}"
    )


if __name__ == "__main__":
    main()
