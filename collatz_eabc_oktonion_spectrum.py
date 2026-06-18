#!/usr/bin/env python3
"""
EABC-Oktanion-Assoziator-Spektrum M_n(t) und M_n^E(t).

Kanonsiche Theorie: collatz_eabc_oktonion_assoziator_spektralhypothese.md §3–§8
                   collatz_eabc_oktonion_singularitaet.md §3.8
                   collatz_eabc_plattenuebergang.md §2.7

Größen:
  α(x,y,z)   = N([x,y,z])
  α_E(x,y,z) = N(Γ_E((xy)z) - Γ_E(x(yz)))
  M_n(t)     = #{(x,y,z) ∈ 𝔞_n : α = t}
  M_n^E(t)   = #{(x,y,z) ∈ 𝔞_n : α_E = t}

Sampling (n ≤ max_n): Stichproben aus 𝔞_n auf Z^8-Stub; kein volles Σ_n-Zählen.

Ausführung:
    python3 collatz_eabc_oktonion_spectrum.py
    python3 collatz_eabc_oktonion_spectrum.py --max-n 30 --samples 40
"""

from __future__ import annotations

import argparse
import json
import math
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from collatz_eabc_oktonion_associator import (
    associator,
    eabc_associator_vector,
    o_norm_sq,
    shell_z8,
)

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_oktonion_spectrum.json"

Oct = tuple[int, ...]


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def triple_norm_factorizations(n: int) -> list[tuple[int, int, int]]:
    """Alle geordneten Tripel (a,b,c) mit a*b*c=n und a,b,c>=1."""
    if n < 1:
        return []
    out: list[tuple[int, int, int]] = []
    for a in range(1, int(n ** (1 / 3)) + 2):
        if n % a != 0:
            continue
        rest = n // a
        for b in range(1, int(math.isqrt(rest)) + 1):
            if rest % b != 0:
                continue
            c = rest // b
            out.append((a, b, c))
    return out


def alpha_algebraic(x: Oct, y: Oct, z: Oct) -> int:
    """α(x,y,z) = N([x,y,z])."""
    return o_norm_sq(associator(x, y, z))


def alpha_eabc(x: Oct, y: Oct, z: Oct) -> int:
    """α_E(x,y,z) = N(Γ_E((xy)z) - Γ_E(x(yz)))."""
    v = eabc_associator_vector(x, y, z)
    return sum(c * c for c in v)


def _sample_shell(n: int, k: int, rng: random.Random) -> list[Oct]:
    pts = list(shell_z8(n))
    if not pts:
        return []
    if len(pts) <= k:
        return pts
    return rng.sample(pts, k)


def sample_triples_for_n(
    n: int,
    samples: int,
    rng: random.Random,
    max_shell_points: int = 300,
) -> list[tuple[Oct, Oct, Oct]]:
    """Stichproben-Tripel mit N(x)N(y)N(z)=n."""
    facs = triple_norm_factorizations(n)
    if not facs:
        return []
    triples: list[tuple[Oct, Oct, Oct]] = []
    shells: dict[tuple[int, int, int], tuple[list[Oct], list[Oct], list[Oct]]] = {}
    for a, b, c in facs:
        key = (a, b, c)
        if key not in shells:
            xs = _sample_shell(a, max_shell_points, rng)
            ys = _sample_shell(b, max_shell_points, rng)
            zs = _sample_shell(c, max_shell_points, rng)
            shells[key] = (xs, ys, zs)
        xs, ys, zs = shells[key]
        if not xs or not ys or not zs:
            continue
        per_fac = max(1, samples // len(facs))
        for _ in range(per_fac):
            triples.append((rng.choice(xs), rng.choice(ys), rng.choice(zs)))
    return triples


def build_histogram(values: list[int]) -> dict[str, int]:
    return {str(k): v for k, v in sorted(Counter(values).items())}


def spectrum_entropy(hist: dict[str, int]) -> float | None:
    total = sum(hist.values())
    if total == 0:
        return None
    h = 0.0
    for c in hist.values():
        p = c / total
        if p > 0:
            h -= p * math.log(p)
    return h


def spectrum_stats(hist: dict[str, int]) -> dict[str, Any]:
    if not hist:
        return {
            "support_size": 0,
            "total_count": 0,
            "entropy": None,
            "mean_level": None,
            "var_level": None,
        }
    levels = [int(k) for k in hist for _ in range(hist[k])]
    total = len(levels)
    return {
        "support_size": len(hist),
        "total_count": total,
        "entropy": spectrum_entropy(hist),
        "mean_level": statistics.mean(levels),
        "var_level": statistics.pvariance(levels) if len(levels) > 1 else 0.0,
    }


def I_M_n_E(hist: dict[str, int]) -> dict[str, Any]:
    """I(M_n^E) = (H, support, mean, var)."""
    st = spectrum_stats(hist)
    return {
        "H": st["entropy"],
        "support_size": st["support_size"],
        "mean_alpha_E": st["mean_level"],
        "var_alpha_E": st["var_level"],
    }


def normalize_hist(hist: dict[str, int]) -> dict[str, float]:
    total = sum(hist.values())
    if total == 0:
        return {}
    return {k: v / total for k, v in hist.items()}


def kl_divergence(p: dict[str, float], q: dict[str, float], eps: float = 1e-12) -> float | None:
    """KL(p || q) über Vereinigung der Stützpunkte."""
    if not p:
        return None
    keys = set(p) | set(q)
    if not keys:
        return None
    kl = 0.0
    for k in keys:
        pk = p.get(k, 0.0)
        qk = q.get(k, eps)
        if pk > 0:
            kl += pk * math.log(pk / qk)
    return kl


def pooled_composite_reference(
    rows: list[dict[str, Any]], spectrum_key: str = "M_n_E"
) -> dict[str, float]:
    pooled: Counter[str] = Counter()
    for r in rows:
        if r.get("is_prime"):
            continue
        hist = r.get(spectrum_key, {})
        for k, v in hist.items():
            pooled[k] += v
    return normalize_hist(dict(pooled))


def partial_zeta(hist: dict[str, int], s: float) -> float | None:
    """S_n(s) = Σ_t m_n(t) / t^s mit m_n normalisiert."""
    p = normalize_hist(hist)
    if not p:
        return None
    return sum(prob / (float(t) ** s) for t, prob in p.items() if float(t) > 0)


def spectrum_for_n(
    n: int,
    samples: int,
    rng: random.Random,
) -> dict[str, Any]:
    triples = sample_triples_for_n(n, samples, rng)
    alpha_vals: list[int] = []
    alpha_e_vals: list[int] = []
    for x, y, z in triples:
        alpha_vals.append(alpha_algebraic(x, y, z))
        alpha_e_vals.append(alpha_eabc(x, y, z))

    hist_alg = build_histogram(alpha_vals)
    hist_e = build_histogram(alpha_e_vals)

    return {
        "n": n,
        "is_prime": _is_prime(n),
        "triple_factorizations": len(triple_norm_factorizations(n)),
        "sample_triples": len(triples),
        "M_n": hist_alg,
        "M_n_E": hist_e,
        "I_M_n": spectrum_stats(hist_alg),
        "I_M_n_E": I_M_n_E(hist_e),
        "S_n_s1": partial_zeta(hist_alg, 1.0),
        "S_n_s2": partial_zeta(hist_alg, 2.0),
        "S_n_E_s1": partial_zeta(hist_e, 1.0),
        "S_n_E_s2": partial_zeta(hist_e, 2.0),
    }


def prime_vs_composite_kl(
    rows: list[dict[str, Any]],
) -> dict[str, Any]:
    ref = pooled_composite_reference(rows, "M_n_E")
    ref_alg = pooled_composite_reference(rows, "M_n")

    prime_kls: list[float] = []
    comp_kls: list[float] = []
    prime_details: list[dict[str, Any]] = []

    for r in rows:
        p_dist = normalize_hist(r.get("M_n_E", {}))
        if not p_dist:
            continue
        kl = kl_divergence(p_dist, ref)
        if kl is None:
            continue
        entry = {"n": r["n"], "kl_M_n_E_vs_composite_ref": kl}
        if r.get("is_prime"):
            prime_kls.append(kl)
            prime_details.append(entry)
        else:
            comp_kls.append(kl)

    def _mean(xs: list[float]) -> float | None:
        return statistics.mean(xs) if xs else None

    ratio = None
    mp, mc = _mean(prime_kls), _mean(comp_kls)
    if mp is not None and mc is not None and mc > 0:
        ratio = mp / mc

    distinguishable = (
        mp is not None
        and mc is not None
        and mp > mc * 1.1
        and len(prime_kls) >= 2
    )

    return {
        "composite_reference_support": len(ref),
        "composite_reference_M_n_support": len(ref_alg),
        "prime_count_with_spectrum": len(prime_kls),
        "composite_count_with_spectrum": len(comp_kls),
        "mean_kl_prime_vs_composite_ref": mp,
        "mean_kl_composite_vs_composite_ref": mc,
        "ratio_mean_kl_prime_over_composite": ratio,
        "prime_kl_details": prime_details,
        "distinguishable_prime_spectrum": distinguishable,
        "verdict": (
            "kein systematisches Prim-Spektrum (KL_prim ≈ KL_comp)"
            if mp is not None and mc is not None and abs(mp - mc) < 0.05 * max(mc, 1e-9)
            else (
                "mögliche Prim-Spektral-Separation (explorativ)"
                if distinguishable
                else "schwaches oder gemischtes Signal (explorativ)"
            )
        ),
    }


def falsification_verdict(
    rows: list[dict[str, Any]],
    kl_report: dict[str, Any],
) -> dict[str, Any]:
    mp = kl_report.get("mean_kl_prime_vs_composite_ref")
    mc = kl_report.get("mean_kl_composite_vs_composite_ref")
    trivial_symmetry = (
        mp is not None
        and mc is not None
        and abs(mp - mc) < 0.02
    )
    substance = kl_report.get("distinguishable_prime_spectrum", False)

    return {
        "trivial_symmetries_detected": trivial_symmetry,
        "substance_reproducible_prime_anomalies": substance,
        "overall": (
            "falsifiziert (triviale Symmetrie: Prim ≈ Composite KL-Profil)"
            if trivial_symmetry
            else (
                "substanz (explorativ): reproduzierbare Prim-Spektral-Anomalie"
                if substance
                else "offen: weder klar falsifiziert noch substanziell gestützt"
            )
        ),
        "epistemic_note": (
            "Nichtassoziativität allein ≠ Prim; Test prüft M_n^E mit Normniveau-Bindung. "
            "Sampling auf Z^8, nicht Hurwitz O_H."
        ),
    }


def example_compare(n_comp: int, n_prime: int, samples: int, seed: int) -> dict[str, Any]:
    rng = random.Random(seed)
    c = spectrum_for_n(n_comp, samples, rng)
    p = spectrum_for_n(n_prime, samples, rng)
    kl_cp = kl_divergence(
        normalize_hist(c.get("M_n_E", {})),
        normalize_hist(p.get("M_n_E", {})),
    )
    return {
        "composite_n": n_comp,
        "prime_n": n_prime,
        "M_n_E_composite": c.get("M_n_E"),
        "M_n_E_prime": p.get("M_n_E"),
        "I_composite": c.get("I_M_n_E"),
        "I_prime": p.get("I_M_n_E"),
        "kl_composite_vs_prime": kl_cp,
        "distinguishable": kl_cp is not None and kl_cp > 0.1,
    }


def spectrum_report(
    max_n: int = 30,
    samples: int = 40,
    seed: int = 54,
) -> dict[str, Any]:
    rng = random.Random(seed)
    rows = [spectrum_for_n(n, samples, rng) for n in range(2, max_n + 1)]
    kl = prime_vs_composite_kl(rows)
    fals = falsification_verdict(rows, kl)
    ex = example_compare(6, 7, samples, seed)

    return {
        "meta": {
            "hypothesis_doc": [
                "collatz_eabc_oktonion_assoziator_spektralhypothese.md",
                "collatz_eabc_oktonion_singularitaet.md §3.8",
                "collatz_eabc_plattenuebergang.md §2.7",
            ],
            "script": "collatz_eabc_oktonion_spectrum.py",
            "associator_module": "collatz_eabc_oktonion_associator.py",
            "lattice": "Z^8 stub (nicht Hurwitz O_H)",
            "max_n": max_n,
            "samples_per_n": samples,
            "seed": seed,
            "definitions": {
                "alpha": "N([x,y,z])",
                "alpha_E": "N(Γ_E((xy)z) - Γ_E(x(yz)))",
                "M_n(t)": "#{ (x,y,z) in 𝔞_n : alpha=t }",
                "M_n_E(t)": "#{ (x,y,z) in 𝔞_n : alpha_E=t }",
            },
            "dirichlet_link": "collatz_eabc_dirichlet_D.py (Quaternion D̂(s) Referenz)",
            "not_computed": [
                "full 𝔞_n enumeration",
                "Hurwitz O_H lattice",
                "D̂_E(s) partial sums (future)",
            ],
        },
        "example_n6_vs_n7": ex,
        "rows": rows,
        "prime_vs_composite_kl": kl,
        "falsification": fals,
    }


def format_summary(report: dict[str, Any]) -> str:
    ex = report["example_n6_vs_n7"]
    kl = report["prime_vs_composite_kl"]
    fals = report["falsification"]
    lines = [
        "EABC-Oktanion-Assoziator-Spektrum M_n / M_n^E",
        "=" * 44,
        f"n ≤ {report['meta']['max_n']}, samples={report['meta']['samples_per_n']}",
        "",
        f"Beispiel n=6 (zsg.) vs n=7 (prim):",
        f"  M_n^E(6) = {ex['M_n_E_composite']}",
        f"  M_n^E(7) = {ex['M_n_E_prime']}",
        f"  KL(6||7) = {ex['kl_composite_vs_prime']}",
        f"  unterscheidbar? {ex['distinguishable']}",
        "",
        f"KL Prim vs Composite-Ref: mean={kl['mean_kl_prime_vs_composite_ref']}",
        f"KL Composite vs Composite-Ref: mean={kl['mean_kl_composite_vs_composite_ref']}",
        f"Ratio Prim/Comp = {kl['ratio_mean_kl_prime_over_composite']}",
        f"Verdict KL: {kl['verdict']}",
        "",
        f"Falsifikation: {fals['overall']}",
    ]
    return "\n".join(lines)


def run(
    max_n: int = 30,
    samples: int = 40,
    seed: int = 54,
    output: Path | None = None,
) -> dict[str, Any]:
    report = spectrum_report(max_n=max_n, samples=samples, seed=seed)
    out = output or DEFAULT_OUTPUT
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report["output_path"] = str(out)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC-Oktanion-Assoziator-Spektrum")
    parser.add_argument("--max-n", type=int, default=30)
    parser.add_argument("--samples", type=int, default=40)
    parser.add_argument("--seed", type=int, default=54)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = run(
        max_n=args.max_n,
        samples=args.samples,
        seed=args.seed,
        output=args.output,
    )
    print(format_summary(report))
    print(f"\nJSON: {report['output_path']}")


if __name__ == "__main__":
    main()
