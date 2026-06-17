#!/usr/bin/env python3
"""
Σ→p-Schalen-Defekt-Experiment: I(μ_n), Referenz I_ref(n), Anomalie D(n).

Kanonsiche Theorie: collatz_eabc_quaternion_mass_hypothese.md §11–§12
(EABC-Spektralgeometrische Hauptvermutung).

Perspektivwechsel: nicht p ↦ Σ_p, sondern n ↦ Σ_n ↦ μ_n; Primzahlen als
mögliche Spektralanomalien, an denen |D(n)| = |I(μ_n) − I_ref(n)| außergewöhnlich ist.

I_ref-Varianten (epistemisch kritisch — gleiche Arithmetik, andere Sprache?):
  rolling     — Nachbar-Mittel über n±w
  cumulative  — Präfix-Mittel Ī(n) = mean_{m<n} I(μ_m)
  omega       — ω(n)-Stratum (verschiedene Primfaktoren)
  tau         — τ(n)-Stratum (Teileranzahl)
  mu_infinity — globales Mittel als μ_∞-Proxy

Ausführung:
    python3 collatz_eabc_shell_defekt_test.py
    python3 collatz_eabc_shell_defekt_test.py --max-n 200 --output collatz_eabc_shell_defekt.json
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from collections import Counter, defaultdict
from dataclasses import dataclass
from math import isqrt
from pathlib import Path
from typing import Any

from collatz_eabc_hurwitz_orbit_test import (
    K_p_covariance,
    chi_p_from_counts,
    classify_shell_point,
    entropy_from_counts,
    hurwitz_shell_elements,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_shell_defekt.json"


def _sieve_primes(limit: int) -> set[int]:
    if limit < 2:
        return set()
    flags = [True] * (limit + 1)
    flags[0] = flags[1] = False
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            step = p
            start = p * p
            flags[start : limit + 1 : step] = [False] * (((limit - start) // step) + 1)
    return {i for i, ok in enumerate(flags) if ok}


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = isqrt(n)
    for k in range(3, d + 1, 2):
        if n % k == 0:
            return False
    return True


def omega_distinct(n: int) -> int:
    """ω(n): Anzahl verschiedener Primfaktoren."""
    if n < 2:
        return 0
    count = 0
    x = n
    p = 2
    while p * p <= x:
        if x % p == 0:
            count += 1
            while x % p == 0:
                x //= p
        p += 1 if p == 2 else 2
    if x > 1:
        count += 1
    return count


def tau_divisor_count(n: int) -> int:
    """τ(n): Anzahl positiver Teiler."""
    if n < 1:
        return 0
    count = 0
    k = 1
    while k * k <= n:
        if n % k == 0:
            count += 1 if k * k == n else 2
        k += 1
    return count


def mean_I_vector(rows: list[ShellInvariants]) -> tuple[float, float, float]:
    if not rows:
        return 0.0, 0.0, 0.0
    return (
        statistics.mean(r.H_n for r in rows),
        statistics.mean(r.chi_n for r in rows),
        statistics.mean(r.K_trace for r in rows),
    )


@dataclass(frozen=True, slots=True)
class ShellInvariants:
    n: int
    shell_size: int
    H_n: float
    chi_n: float
    K_trace: float
    distinct_gamma: int
    is_prime: bool
    omega: int
    tau: int

    def I_vector(self) -> tuple[float, float, float]:
        return (self.H_n, self.chi_n, self.K_trace)

    def to_dict(self) -> dict[str, Any]:
        return {
            "n": self.n,
            "shell_size": self.shell_size,
            "H_n": round(self.H_n, 6),
            "chi_n": round(self.chi_n, 6),
            "K_trace": round(self.K_trace, 6),
            "distinct_gamma": self.distinct_gamma,
            "is_prime": self.is_prime,
            "omega": self.omega,
            "tau": self.tau,
        }


def compute_shell_invariants(n: int) -> ShellInvariants | None:
    raw = hurwitz_shell_elements(n)
    if not raw:
        return None
    points = [classify_shell_point(*coords) for coords in raw]
    gamma_counts: Counter = Counter(pt.gamma for pt in points)
    size = len(points)
    H_n = entropy_from_counts(gamma_counts, size)
    chi_n = chi_p_from_counts(gamma_counts, size)
    K_trace = K_p_covariance(points)["trace"]
    return ShellInvariants(
        n=n,
        shell_size=size,
        H_n=H_n,
        chi_n=chi_n,
        K_trace=K_trace,
        distinct_gamma=len(gamma_counts),
        is_prime=is_prime(n),
        omega=omega_distinct(n),
        tau=tau_divisor_count(n),
    )


def rolling_baseline(
    rows: list[ShellInvariants], idx: int, window: int = 5
) -> tuple[float, float, float]:
    """Mittel von I über Nachbarn n±window (ohne n selbst)."""
    n = rows[idx].n
    lo, hi = max(2, n - window), n + window
    peers = [r for r in rows if lo <= r.n <= hi and r.n != n]
    if not peers:
        return rows[idx].I_vector()
    H = statistics.mean(r.H_n for r in peers)
    chi = statistics.mean(r.chi_n for r in peers)
    K = statistics.mean(r.K_trace for r in peers)
    return H, chi, K


def omega_baseline(rows: list[ShellInvariants], target: ShellInvariants) -> tuple[float, float, float]:
    """Mittel von I über alle m mit gleichem ω(m), m ≠ n."""
    peers = [r for r in rows if r.omega == target.omega and r.n != target.n]
    if not peers:
        return mean_I_vector(rows)
    return mean_I_vector(peers)


def tau_baseline(rows: list[ShellInvariants], target: ShellInvariants) -> tuple[float, float, float]:
    """Mittel von I über alle m mit gleichem τ(m), m ≠ n."""
    peers = [r for r in rows if r.tau == target.tau and r.n != target.n]
    if not peers:
        return mean_I_vector(rows)
    return mean_I_vector(peers)


def cumulative_baseline(rows: list[ShellInvariants], idx: int) -> tuple[float, float, float]:
    """Präfix-Mittel Ī(n) = mean_{m<n} I(μ_m) über alle vorherigen Schalen."""
    peers = rows[:idx]
    if not peers:
        return rows[idx].I_vector()
    return mean_I_vector(peers)


def defect_magnitude(dH: float, dchi: float, dK: float, K_scale: float = 1.0) -> float:
    """Skalierte euklidische Norm des Defektvektors."""
    return math.sqrt(dH * dH + dchi * dchi + (dK / K_scale) ** 2)


def defect_components(
    inv: ShellInvariants,
    I_ref: tuple[float, float, float],
    K_scale: float,
) -> dict[str, Any]:
    dH = inv.H_n - I_ref[0]
    dchi = inv.chi_n - I_ref[1]
    dK = inv.K_trace - I_ref[2]
    return {
        "H": round(dH, 6),
        "chi": round(dchi, 6),
        "K_trace": round(dK, 6),
        "magnitude": round(defect_magnitude(dH, dchi, dK, K_scale), 6),
    }


def shell_defekt_row(
    inv: ShellInvariants,
    baselines: dict[str, tuple[float, float, float]],
    K_scale: float,
    V_n: float,
) -> dict[str, Any]:
    row: dict[str, Any] = {**inv.to_dict(), "V_n_bernoulli_proxy": round(V_n, 6)}
    for name, I_ref in baselines.items():
        row[f"I_ref_{name}"] = [round(x, 6) for x in I_ref]
        row[f"D_{name}"] = defect_components(inv, I_ref, K_scale)
    return row


I_REF_VARIANTS: dict[str, str] = {
    "rolling": "Nachbar-Mittel über n±w (lokale Glättung; Primzahlen oft zwischen dichten zusammengesetzten n)",
    "cumulative": "Präfix-Mittel Ī(n)=mean_{m<n}I(μ_m) (spektralgeometrische Referenz aus der Folge selbst)",
    "omega": "ω(n)-Stratum: Mittel über gleiche Anzahl verschiedener Primfaktoren",
    "tau": "τ(n)-Stratum: Mittel über gleiche Teileranzahl (divisor-count baseline)",
    "mu_infinity": "Globales Mittel als μ_∞-Proxy über alle untersuchten Schalen",
}


def bernoulli_V_n(inv: ShellInvariants, I_inf: tuple[float, float, float], K_scale: float) -> float:
    """Bernoulli-Proxy V_n = ||I(μ_n) − I_∞|| — Diskret-vs.-kontinuierlich-Defekt (Heuristik)."""
    return defect_magnitude(
        inv.H_n - I_inf[0],
        inv.chi_n - I_inf[1],
        inv.K_trace - I_inf[2],
        K_scale,
    )


def pearson_correlation(xs: list[float], ys: list[float]) -> float | None:
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


def bernoulli_D_correlation(rows: list[dict[str, Any]], baseline_keys: list[str]) -> dict[str, Any]:
    """Explorativ: Korrelation V_n (μ_n−μ_∞-Proxy) mit |D(n)| je I_ref."""
    V_vals = [r["V_n_bernoulli_proxy"] for r in rows]
    out: dict[str, Any] = {}
    for key in baseline_keys:
        dkey = f"D_{key}"
        d_vals = [r[dkey]["magnitude"] for r in rows]
        rho = pearson_correlation(V_vals, d_vals)
        out[key] = {
            "pearson_Vn_vs_abs_D": round(rho, 6) if rho is not None else None,
            "note": "Explorativ; collatz_eabc_bernoulli_sensor (Branch eabc-bernoulli-sensor) nicht lokal — V_n als I−I_∞-Proxy",
        }
    return out


def pick_best_iref(comparisons: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Wählt I_ref mit größtem ratio_mean_prime_over_composite (falls >1), sonst ehrlich 'keiner'."""
    ranked: list[tuple[str, float, bool]] = []
    for name, cmp in comparisons.items():
        if cmp.get("insufficient_data"):
            continue
        ratio = cmp.get("ratio_mean_prime_over_composite", 0.0)
        ranked.append((name, ratio, bool(cmp.get("primes_larger_on_mean"))))
    ranked.sort(key=lambda t: t[1], reverse=True)
    if not ranked:
        return {"candidate": None, "reason": "insufficient_data"}
    best_name, best_ratio, primes_larger = ranked[0]
    stable = primes_larger and best_ratio >= 1.05
    return {
        "candidate": best_name if stable else None,
        "best_ratio_name": best_name,
        "best_ratio": round(best_ratio, 6),
        "primes_larger_on_mean": primes_larger,
        "stable_prime_anomaly": stable,
        "ranking": [
            {"iref": n, "ratio_mean_prime_over_composite": round(r, 6), "primes_larger_on_mean": pl}
            for n, r, pl in ranked
        ],
        "epistemic_note": (
            "ω-/τ-Strata und μ_∞-Proxy können Prim-Überhang zeigen, weil ω(p)=1 und τ(p)=2 "
            "arithmetisch fixiert sind — nicht zwingend neue Geometrie. Rolling-Baseline zeigt bei kleinem n "
            "typischerweise keinen robusten Prim-Überhang."
        ),
    }


def prime_vs_composite_comparison(rows: list[dict[str, Any]], key: str = "D_rolling") -> dict[str, Any]:
    primes = [r[key]["magnitude"] for r in rows if r["is_prime"]]
    composites = [r[key]["magnitude"] for r in rows if not r["is_prime"]]
    if not primes or not composites:
        return {"key": key, "insufficient_data": True}

    mean_p = statistics.mean(primes)
    mean_c = statistics.mean(composites)
    med_p = statistics.median(primes)
    med_c = statistics.median(composites)
    ratio_mean = mean_p / mean_c if mean_c > 1e-12 else float("inf")

    ranked = sorted(rows, key=lambda r: r[key]["magnitude"], reverse=True)
    top10 = [
        {
            "n": r["n"],
            "is_prime": r["is_prime"],
            "magnitude": r[key]["magnitude"],
            "omega": r["omega"],
            "tau": r.get("tau"),
        }
        for r in ranked[:10]
    ]
    prime_in_top10 = sum(1 for t in top10 if t["is_prime"])

    return {
        "key": key,
        "prime_count": len(primes),
        "composite_count": len(composites),
        "mean_abs_D_prime": round(mean_p, 6),
        "mean_abs_D_composite": round(mean_c, 6),
        "median_abs_D_prime": round(med_p, 6),
        "median_abs_D_composite": round(med_c, 6),
        "ratio_mean_prime_over_composite": round(ratio_mean, 6),
        "primes_larger_on_mean": mean_p > mean_c,
        "top10_outliers": top10,
        "primes_in_top10": prime_in_top10,
        "epistemic_note": (
            "Exploratives Muster — kein Beweis, dass Primzahlen systematisch größere |D(n)| tragen. "
            "Kleine N, Schalengröße korreliert mit n; I_ref-Wahl beeinflusst Verhältnis."
        ),
    }


def bernoulli_discretization_hint(
    rows: list[ShellInvariants], I_inf: tuple[float, float, float]
) -> dict[str, Any]:
    """Heuristik: μ_n − μ_∞ als Diskretisierungskorrektur (Bernoulli-Analogie)."""
    if len(rows) < 4:
        return {"hint": "insufficient_data"}
    prime_rows = [r for r in rows if r.is_prime]
    comp_rows = [r for r in rows if not r.is_prime]
    return {
        "mu_infinity_proxy": {
            "H": round(I_inf[0], 6),
            "chi": round(I_inf[1], 6),
            "K_trace": round(I_inf[2], 6),
        },
        "prime_mean_delta_H": round(statistics.mean(r.H_n - I_inf[0] for r in prime_rows), 6)
        if prime_rows
        else 0.0,
        "composite_mean_delta_H": round(statistics.mean(r.H_n - I_inf[0] for r in comp_rows), 6)
        if comp_rows
        else 0.0,
        "heuristic": "μ_n − μ_∞ als endliche-Schalen-Korrektur (Bernoulli-Heuristik, nicht bewiesen)",
        "V_n_note": "V_n = ||I(μ_n)−I_∞|| pro Zeile; Korrelation mit |D(n)| in bernoulli_D_correlation",
    }


def shell_defekt_report(max_n: int = 50, rolling_window: int = 5) -> dict[str, Any]:
    invariants: list[ShellInvariants] = []
    for n in range(2, max_n + 1):
        inv = compute_shell_invariants(n)
        if inv is not None:
            invariants.append(inv)

    if not invariants:
        return {"max_n": max_n, "rows": [], "verdict": "no_shells"}

    K_scale = max(statistics.mean(abs(r.K_trace) for r in invariants), 0.01)
    I_global = mean_I_vector(invariants)
    baseline_keys = list(I_REF_VARIANTS.keys())

    rows: list[dict[str, Any]] = []
    for i, inv in enumerate(invariants):
        baselines = {
            "rolling": rolling_baseline(invariants, i, window=rolling_window),
            "cumulative": cumulative_baseline(invariants, i),
            "omega": omega_baseline(invariants, inv),
            "tau": tau_baseline(invariants, inv),
            "mu_infinity": I_global,
        }
        V_n = bernoulli_V_n(inv, I_global, K_scale)
        rows.append(shell_defekt_row(inv, baselines, K_scale, V_n))

    comparisons = {
        name: prime_vs_composite_comparison(rows, f"D_{name}") for name in baseline_keys
    }
    best = pick_best_iref(comparisons)
    bernoulli = bernoulli_discretization_hint(invariants, I_global)
    b_corr = bernoulli_D_correlation(rows, baseline_keys)

    cmp_roll = comparisons["rolling"]
    verdict_parts: list[str] = []
    if cmp_roll.get("primes_larger_on_mean"):
        verdict_parts.append(
            f"rolling: Prim-Mittel |D|={cmp_roll['mean_abs_D_prime']:.4f} > "
            f"zusammengesetzt {cmp_roll['mean_abs_D_composite']:.4f} "
            f"(Ratio {cmp_roll['ratio_mean_prime_over_composite']:.3f})"
        )
    else:
        verdict_parts.append(
            f"rolling: kein Prim-Überhang (Prim {cmp_roll['mean_abs_D_prime']:.4f} "
            f"vs zusammengesetzt {cmp_roll['mean_abs_D_composite']:.4f})"
        )
    if best.get("candidate"):
        verdict_parts.append(
            f"bester I_ref-Kandidat: {best['candidate']} (Ratio {best['best_ratio']:.3f})"
        )
    else:
        verdict_parts.append(
            f"kein stabiler Prim-I_ref-Kandidat (beste Ratio: {best.get('best_ratio_name')} "
            f"{best.get('best_ratio')})"
        )
    verdict_parts.append(f"Top-10-|D| rolling: {cmp_roll.get('primes_in_top10', 0)}/10 prim")

    return {
        "meta": {
            "hypothesis_doc": "collatz_eabc_quaternion_mass_hypothese.md §12",
            "perspective": "EABC-Spektralgeometrische Hauptvermutung: n ↦ Σ_n ↦ μ_n fundamental",
            "invariants": "I(μ_n) = (H_n, χ_n, K_trace)",
            "defect": "D(n) = I(μ_n) − I_ref(n)",
            "I_ref_variants": I_REF_VARIANTS,
            "max_n": max_n,
            "rolling_window": rolling_window,
            "shell_count": len(rows),
            "K_scale_for_magnitude": round(K_scale, 6),
        },
        "bernoulli_heuristic": bernoulli,
        "bernoulli_D_correlation": b_corr,
        "prime_vs_composite": comparisons,
        "best_I_ref": best,
        "rows": rows,
        "verdict": "; ".join(verdict_parts),
        "epistemic_status": (
            "Experiment — testet EABC-Spektralgeometrische Hauptvermutung (§12). "
            "Rolling-Baseline bei n≤50 typischerweise ohne Prim-Überhang; ω-/τ-Strata und μ_∞-Proxy "
            "können Ratio>1 zeigen, oft repackagierte Arithmetik (ω(p)=1). Keine Falsifikation der "
            "Emergenz-Idee, aber I_ref muss sorgfältig gewählt werden."
        ),
    }


def multi_scale_summary(max_ns: list[int], rolling_window: int = 5) -> dict[str, Any]:
    """Vergleicht Prim-Separation über mehrere n-Obergrenzen."""
    scales: dict[str, Any] = {}
    for n in max_ns:
        rep = shell_defekt_report(max_n=n, rolling_window=rolling_window)
        scales[str(n)] = {
            "best_I_ref": rep.get("best_I_ref"),
            "rolling_ratio": rep["prime_vs_composite"]["rolling"].get("ratio_mean_prime_over_composite"),
            "omega_ratio": rep["prime_vs_composite"]["omega"].get("ratio_mean_prime_over_composite"),
            "tau_ratio": rep["prime_vs_composite"]["tau"].get("ratio_mean_prime_over_composite"),
            "mu_infinity_ratio": rep["prime_vs_composite"]["mu_infinity"].get(
                "ratio_mean_prime_over_composite"
            ),
            "verdict": rep["verdict"],
        }
    return {"scales": scales, "rolling_window": rolling_window}


def format_table(report: dict[str, Any]) -> str:
    lines = [
        "Σ→p-Schalen-Anomalie D(n) = I(μ_n) − I_ref(n)  [§12]",
        "=" * 58,
        f"n ∈ [2,{report['meta']['max_n']}], Schalen: {report['meta']['shell_count']}",
        f"Verdict: {report['verdict']}",
        "",
    ]
    best = report.get("best_I_ref", {})
    if best.get("ranking"):
        lines.append("I_ref Prim/Composite-Ratios (mean |D|):")
        for item in best["ranking"]:
            flag = " *" if item["iref"] == best.get("best_ratio_name") else ""
            lines.append(
                f"  {item['iref']}: ratio={item['ratio_mean_prime_over_composite']:.3f}"
                f"  prim>{'comp' if item['primes_larger_on_mean'] else 'comp?'}{flag}"
            )
        lines.append("")
    pvc = report["prime_vs_composite"]["rolling"]
    lines.extend(
        [
            f"rolling: Prim {pvc['mean_abs_D_prime']:.4f}  comp {pvc['mean_abs_D_composite']:.4f}  "
            f"ratio={pvc['ratio_mean_prime_over_composite']:.3f}",
            "",
            "Top Ausreißer (rolling |D|):",
        ]
    )
    for t in pvc["top10_outliers"][:6]:
        tag = "prim" if t["is_prime"] else f"ω={t['omega']}"
        lines.append(f"  n={t['n']} ({tag}): |D|={t['magnitude']:.4f}")
    return "\n".join(lines)


def run(
    max_n: int = 100,
    rolling_window: int = 5,
    output: Path | None = None,
    include_multi_scale: bool = True,
) -> dict[str, Any]:
    report = shell_defekt_report(max_n=max_n, rolling_window=rolling_window)
    if include_multi_scale:
        scale_ns = sorted({50, 100, min(200, max_n)} | {max_n})
        report["multi_scale"] = multi_scale_summary(scale_ns, rolling_window=rolling_window)
    out = output or DEFAULT_OUTPUT
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report["output_path"] = str(out)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Σ→p-Schalen-Defekt / Spektralgeometrie-Experiment")
    parser.add_argument("--max-n", type=int, default=100, help="Oberes n für Σ_n")
    parser.add_argument("--rolling-window", type=int, default=5, help="Fenster für rollendes I_ref")
    parser.add_argument("--no-multi-scale", action="store_true", help="Kein multi_scale-Block im JSON")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(
        max_n=args.max_n,
        rolling_window=args.rolling_window,
        output=args.output,
        include_multi_scale=not args.no_multi_scale,
    )
    print(format_table(report))
    print(f"\nJSON: {report['output_path']}")


if __name__ == "__main__":
    main()
