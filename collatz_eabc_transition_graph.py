#!/usr/bin/env python3
"""
Gerichteter EABC-Übergangsgraph und Zyklus-Holonomie.

Kanonsiche Theorie:
  collatz_eabc_zyklus_holonomie.md  — χ_E(N) auf Primfolge-Gleitfenstern, Hol_E
  collatz_eabc_transport.md         — G_E, Transport T_n
  collatz_eabc_holonomie.md         — χ_E^quad auf Prim-Vierlingen (Vergleich)

χ_E(N) = Σ Ω(Q_n) / #{Ω≠0}  auf Gleitfenstern Q_n=(X_n,…,X_{n+3})
χ_E^quad(N) — arithmetische Prim-Vierlinge (collatz_eabc_holonomie_test.py)

Ausführung:
    python3 collatz_eabc_transition_graph.py
    python3 collatz_eabc_transition_graph.py --max-p 1000000 --null-trials 2000
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
from math import isqrt
from pathlib import Path
from typing import Any, Literal

from collatz_eabc_holonomie_test import chi_E as chi_E_quadruplet
from collatz_eabc_holonomie_test import enumerate_quadruplets
from eabc_from_lean import EClass, class_of, t

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_transition_graph.json"

LABELS = ("E", "A", "B", "C")
IDX = {label: i for i, label in enumerate(LABELS)}
ABCE_WORD = "ABCE"
CEAB_WORD = "CEAB"
T_INV = {t(c): c for c in EClass}

Omega = Literal[-1, 0, 1]


def _sieve_primes(limit: int) -> list[int]:
    if limit < 2:
        return []
    flags = [True] * (limit + 1)
    flags[0] = flags[1] = False
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            step = p
            start = p * p
            flags[start : limit + 1 : step] = [False] * (((limit - start) // step) + 1)
    return [i for i, ok in enumerate(flags) if ok]


def prime_eabc_sequence(max_p: int) -> list[dict[str, Any]]:
    """Primzahlen p > 3 mit EABC-Klasse."""
    rows: list[dict[str, Any]] = []
    for p in _sieve_primes(max_p):
        if p <= 3:
            continue
        cls = class_of(p)
        if cls is None:
            continue
        rows.append({"p": p, "class": cls.value})
    return rows


def classes_from_sequence(seq: list[dict[str, Any]]) -> list[str]:
    return [row["class"] for row in seq]


def omega_window(word: str) -> Omega:
    """Ω(Q_n)=+1 (ABCE), -1 (CEAB), 0 sonst."""
    if word == ABCE_WORD:
        return 1
    if word == CEAB_WORD:
        return -1
    return 0


def sliding_windows(classes: list[str]) -> list[dict[str, Any]]:
    """Alle 4-Fenster Q_n auf der Klassenfolge."""
    rows: list[dict[str, Any]] = []
    for i in range(len(classes) - 3):
        word = "".join(classes[i : i + 4])
        rows.append({"index": i, "word": word, "omega": omega_window(word)})
    return rows


def chi_E_sliding(classes: list[str]) -> dict[str, Any]:
    """
    χ_E(N) = Σ Ω(Q_n) / #{Ω≠0} auf Primfolge-Gleitfenstern.

    Kanonisch: collatz_eabc_zyklus_holonomie.md §4.
    """
    windows = sliding_windows(classes)
    omega_sum = sum(w["omega"] for w in windows)
    abce = sum(1 for w in windows if w["omega"] == 1)
    ceab = sum(1 for w in windows if w["omega"] == -1)
    nonzero = abce + ceab
    chi = omega_sum / nonzero if nonzero else 0.0
    return {
        "window_count": len(windows),
        "abce_windows": abce,
        "ceab_windows": ceab,
        "omega_sum": omega_sum,
        "nonzero_windows": nonzero,
        "chi_E": chi,
        "formula": "sum(Omega(Q_n)) / #{Omega != 0}",
        "carrier": "prime_sequence_sliding_windows",
    }


def chi_E_quadruplet_report(max_p: int) -> dict[str, Any]:
    """χ_E^quad(N) auf arithmetischen Prim-Vierlingen (holonomie_test)."""
    quads = enumerate_quadruplets(max_p)
    chi = chi_E_quadruplet(max_p, quads)
    return {
        "quadruplet_count": chi["quadruplet_count"],
        "abce_count": chi["abce_count"],
        "ceab_count": chi["ceab_count"],
        "chi_E_quad": chi["chi_E"],
        "formula": chi["formula"],
        "carrier": "prime_quadruplets",
    }


def transition_counts(classes: list[str]) -> list[list[int]]:
    """T_ij = Anzahl Übergänge i → j."""
    mat = [[0 for _ in LABELS] for _ in LABELS]
    for i in range(len(classes) - 1):
        a = IDX[classes[i]]
        b = IDX[classes[i + 1]]
        mat[a][b] += 1
    return mat


def row_stochastic_matrix(counts: list[list[int]]) -> list[list[float]]:
    total = sum(sum(row) for row in counts)
    if total == 0:
        return [[0.0] * 4 for _ in range(4)]
    return [[c / total for c in row] for row in counts]


def transition_probabilities(counts: list[list[int]]) -> list[list[float]]:
    """Zeilenstochastische Matrix P_ij = T_ij / sum_j T_ij."""
    probs: list[list[float]] = []
    for row in counts:
        s = sum(row)
        if s == 0:
            probs.append([0.0] * 4)
        else:
            probs.append([c / s for c in row])
    return probs


def stationary_distribution(counts: list[list[int]], tol: float = 1e-12) -> dict[str, Any]:
    """π via Potenziteration auf zeilenstochastischer Matrix."""
    p_mat = transition_probabilities(counts)
    pi = [0.25, 0.25, 0.25, 0.25]
    for _ in range(10_000):
        new_pi = [0.0] * 4
        for j in range(4):
            for i in range(4):
                new_pi[j] += pi[i] * p_mat[i][j]
        if max(abs(new_pi[k] - pi[k]) for k in range(4)) < tol:
            pi = new_pi
            break
        pi = new_pi
    return {
        "pi": {LABELS[k]: pi[k] for k in range(4)},
        "pi_vector": pi,
    }


def markov_irreducible(counts: list[list[int]]) -> bool:
    """Irreduzibilität via BFS auf gerichtetem Graph mit positivem Gewicht."""
    n = 4
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if counts[i][j] > 0:
                adj[i].append(j)
    seen = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                stack.append(v)
    return len(seen) == n


def edge_frequencies(counts: list[list[int]]) -> dict[str, Any]:
    total = sum(sum(row) for row in counts)
    edges: dict[str, float] = {}
    for i, src in enumerate(LABELS):
        for j, dst in enumerate(LABELS):
            key = f"{src}->{dst}"
            edges[key] = counts[i][j] / total if total else 0.0
    t_aligned = sum(counts[i][(i + 1) % 4] for i in range(4))
    t_inv_aligned = sum(counts[i][(i - 1) % 4] for i in range(4))
    return {
        "total_transitions": total,
        "edge_frequencies": edges,
        "t_forward_fraction": t_aligned / total if total else 0.0,
        "t_inverse_fraction": t_inv_aligned / total if total else 0.0,
        "t_bias": (t_aligned - t_inv_aligned) / total if total else 0.0,
    }


def count_word_windows(classes: list[str], word: str) -> int:
    w = len(word)
    if len(classes) < w:
        return 0
    count = 0
    for i in range(len(classes) - w + 1):
        if "".join(classes[i : i + w]) == word:
            count += 1
    return count


def count_t_cycle_windows(classes: list[str], forward: bool = True) -> int:
    """Vier Schritte mit c_{k+1} = t(c_k) bzw. t^{-1}(c_k)."""
    if len(classes) < 4:
        return 0
    count = 0
    for i in range(len(classes) - 3):
        ok = True
        for k in range(3):
            c0 = EClass(classes[i + k])
            c1 = EClass(classes[i + k + 1])
            expected = t(c0) if forward else T_INV[c0]
            if c1 is not expected:
                ok = False
                break
        if ok:
            count += 1
    return count


def chi_from_bias(pos: int, neg: int) -> float:
    denom = pos + neg
    return (pos - neg) / denom if denom else 0.0


def chi_transport(classes: list[str]) -> dict[str, Any]:
    """Legacy-Alias: χ_trans = χ_E auf Gleitfenstern (identische Formel)."""
    sliding = chi_E_sliding(classes)
    t_fwd = count_t_cycle_windows(classes, forward=True)
    t_inv = count_t_cycle_windows(classes, forward=False)
    return {
        "abce_windows": sliding["abce_windows"],
        "ceab_windows": sliding["ceab_windows"],
        "chi_E": sliding["chi_E"],
        "chi_trans": sliding["chi_E"],
        "t_forward_windows": t_fwd,
        "t_inverse_windows": t_inv,
        "chi_t_cycle": chi_from_bias(t_fwd, t_inv),
        "formula": sliding["formula"],
        "formula_t_cycle": "(#t_forward_4cycles - #t_inverse_4cycles) / (#fwd + #inv)",
    }


def _z_score(observed: float, samples: list[float]) -> float:
    if len(samples) < 2:
        return 0.0
    mu = statistics.mean(samples)
    sigma = statistics.pstdev(samples)
    if sigma == 0:
        return 0.0
    return (observed - mu) / sigma


def _apply_label_permutation(classes: list[str], perm: dict[str, str]) -> list[str]:
    return [perm[c] for c in classes]


def _random_label_permutation(rng: random.Random) -> dict[str, str]:
    labels = list(LABELS)
    rng.shuffle(labels)
    return dict(zip(LABELS, labels))


def shuffle_null_chi_E(
    classes: list[str],
    trials: int = 1000,
    seed: int = 42,
) -> dict[str, Any]:
    """Marginal-Shuffle-Null: permutiere Klassenfolge, erhalte Häufigkeiten."""
    rng = random.Random(seed)
    base = chi_E_sliding(classes)
    samples: list[float] = []
    for _ in range(trials):
        perm = classes[:]
        rng.shuffle(perm)
        samples.append(chi_E_sliding(perm)["chi_E"])
    return {
        "null_type": "marginal_shuffle",
        "trials": trials,
        "seed": seed,
        "observed_chi_E": base["chi_E"],
        "null_mean": statistics.mean(samples),
        "null_std": statistics.pstdev(samples) if trials > 1 else 0.0,
        "z_score": _z_score(base["chi_E"], samples),
    }


def isotropy_null_chi_E(
    classes: list[str],
    trials: int = 1000,
    seed: int = 42,
) -> dict[str, Any]:
    """Isotropie-Null: zufällige Relabeling σ∈S_4 auf {E,A,B,C}."""
    rng = random.Random(seed + 1)
    base = chi_E_sliding(classes)
    samples: list[float] = []
    for _ in range(trials):
        relabeled = _apply_label_permutation(classes, _random_label_permutation(rng))
        samples.append(chi_E_sliding(relabeled)["chi_E"])
    return {
        "null_type": "isotropy_relabel",
        "trials": trials,
        "seed": seed + 1,
        "observed_chi_E": base["chi_E"],
        "null_mean": statistics.mean(samples),
        "null_std": statistics.pstdev(samples) if trials > 1 else 0.0,
        "z_score": _z_score(base["chi_E"], samples),
    }


def hol_E_supported(
    chi_series: list[dict[str, Any]],
    shuffle_null: dict[str, Any],
    isotropy_null: dict[str, Any],
    tol: float = 0.01,
) -> dict[str, Any]:
    """Bewertung: Hol_E≠0 vs. stabile Fluktuationen vs. Null."""
    if not chi_series:
        return {"supported": False, "reason": "empty series"}
    last_chi = chi_series[-1]["chi_E"]
    nonzero_limit = abs(last_chi) > tol
    stable = len(chi_series) >= 2 and all(
        abs(chi_series[i]["chi_E"] - chi_series[-1]["chi_E"]) < 0.05
        for i in range(len(chi_series))
    )
    shuffle_sig = abs(shuffle_null.get("z_score", 0.0)) > 2.0
    isotropy_sig = abs(isotropy_null.get("z_score", 0.0)) > 2.0
    strong = nonzero_limit and (shuffle_sig or isotropy_sig)
    weak = not strong and (shuffle_sig or isotropy_sig or stable)
    return {
        "hol_E_nonzero_at_limit": nonzero_limit,
        "stable_fluctuations": stable,
        "shuffle_significant": shuffle_sig,
        "isotropy_significant": isotropy_sig,
        "hol_E_ne_zero_supported": strong,
        "weak_hypothesis_supported": weak,
        "verdict": (
            f"Hol_E≠0 {'gestützt' if strong else 'nicht gestützt'} "
            f"(χ_E={last_chi:.4f}); "
            f"schwächere Lesart {'gestützt' if weak and not strong else 'nicht gestützt'}."
        ),
    }


def transition_matrix_report(max_p: int) -> dict[str, Any]:
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    counts = transition_counts(classes)
    freqs = edge_frequencies(counts)
    stat = stationary_distribution(counts)
    return {
        "max_p": max_p,
        "prime_count": len(seq),
        "transition_counts": {LABELS[i]: {LABELS[j]: counts[i][j] for j in range(4)} for i in range(4)},
        "transition_probabilities": {
            LABELS[i]: {LABELS[j]: round(transition_probabilities(counts)[i][j], 6) for j in range(4)}
            for i in range(4)
        },
        "joint_frequencies": row_stochastic_matrix(counts),
        "edge_summary": freqs,
        "stationary": stat,
        "irreducible": markov_irreducible(counts),
    }


def chi_sliding_vs_quadruplet(
    max_p: int,
    null_trials: int = 1000,
    seed: int = 42,
) -> dict[str, Any]:
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    sliding = chi_E_sliding(classes)
    quad = chi_E_quadruplet_report(max_p)
    trans = chi_transport(classes)
    shuffle = shuffle_null_chi_E(classes, trials=null_trials, seed=seed)
    isotropy = isotropy_null_chi_E(classes, trials=null_trials, seed=seed)
    same_sign = (sliding["chi_E"] >= 0) == (quad["chi_E_quad"] >= 0)
    return {
        "max_p": max_p,
        "chi_E": sliding["chi_E"],
        "chi_E_quad": quad["chi_E_quad"],
        "chi_t_cycle": trans["chi_t_cycle"],
        "sliding_detail": sliding,
        "quadruplet_detail": quad,
        "transport_detail": trans,
        "null_models": {
            "shuffle": shuffle,
            "isotropy": isotropy,
        },
        "comparison": {
            "same_sign": same_sign,
            "difference": sliding["chi_E"] - quad["chi_E_quad"],
            "verdict": (
                "χ_E (Gleitfenster) und χ_E^quad (Vierlinge) sind verwandte chirale Observablen "
                "in verschiedenen Trägern — nicht identisch. "
                f"Beobachtet: χ_E={sliding['chi_E']:.4f}, "
                f"χ_E^quad={quad['chi_E_quad']:.4f}, "
                f"Δ={sliding['chi_E'] - quad['chi_E_quad']:.4f}."
            ),
        },
    }


# Backward-compatible alias
chi_transport_vs_quadruplet = chi_sliding_vs_quadruplet


def hol_E_estimates(
    limits: list[int],
    null_trials: int = 500,
    seed: int = 42,
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    last_cmp: dict[str, Any] | None = None
    for lim in limits:
        cmp = chi_sliding_vs_quadruplet(lim, null_trials=null_trials, seed=seed)
        last_cmp = cmp
        rows.append(
            {
                "limit": lim,
                "chi_E": cmp["chi_E"],
                "chi_E_quad": cmp["chi_E_quad"],
                "chi_t_cycle": cmp["chi_t_cycle"],
                "z_shuffle": cmp["null_models"]["shuffle"]["z_score"],
                "z_isotropy": cmp["null_models"]["isotropy"]["z_score"],
            }
        )
    support = hol_E_supported(
        rows,
        last_cmp["null_models"]["shuffle"] if last_cmp else {},
        last_cmp["null_models"]["isotropy"] if last_cmp else {},
    )
    return {
        "limits": limits,
        "series": rows,
        "hol_E_estimate": rows[-1]["chi_E"] if rows else 0.0,
        "hol_E_quad_estimate": rows[-1]["chi_E_quad"] if rows else 0.0,
        "hol_E_support": support,
        "note": "Hol_E = lim χ_E(N); empirische Schätzung = letzter χ_E-Wert — kein Beweis des Grenzwerts.",
    }


def run(
    max_p: int = 100_000,
    null_trials: int = 1000,
    output: Path = DEFAULT_OUTPUT,
    seed: int = 42,
) -> dict[str, Any]:
    matrix = transition_matrix_report(max_p)
    comparison = chi_sliding_vs_quadruplet(max_p, null_trials=null_trials, seed=seed)
    limits = sorted({min(max_p, x) for x in (1_000, 5_000, 10_000, 50_000, 100_000, max_p)})
    hol = hol_E_estimates(limits, null_trials=min(null_trials, 300), seed=seed)

    report: dict[str, Any] = {
        "meta": {
            "module": "collatz_eabc_transition_graph.py",
            "theory": "collatz_eabc_zyklus_holonomie.md",
            "transport": "collatz_eabc_transport.md",
            "holonomy_hierarchy": "collatz_eabc_holonomie.md",
            "max_p": max_p,
            "null_trials": null_trials,
        },
        "transition_matrix": matrix,
        "cycle_holonomy": comparison,
        "hol_E_estimates": hol,
        "epistemic_labels": {
            "transition_matrix": "Definition + Experiment",
            "chi_E_sliding": "Definition",
            "chi_E_quad": "Definition (Vierlingsträger)",
            "shuffle_null": "Experiment",
            "isotropy_null": "Experiment",
            "hol_E_limit": "Definition (Grenzwert offen)",
            "hol_E_ne_zero": "Hypothese",
        },
    }
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC-Übergangsgraph und Zyklus-Holonomie")
    parser.add_argument("--max-p", type=int, default=100_000, help="Obergrenze für Primzahlen")
    parser.add_argument("--null-trials", type=int, default=1000, help="Null-Durchläufe")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    report = run(max_p=args.max_p, null_trials=args.null_trials, output=args.output, seed=args.seed)
    cmp = report["cycle_holonomy"]
    mat = report["transition_matrix"]
    hol = report["hol_E_estimates"]
    print("=== EABC-Zyklus-Holonomie ===")
    print(f"Primzahlen >3 bis {args.max_p}: {mat['prime_count']}")
    print(f"χ_E (Gleitfenster) = {cmp['chi_E']:.4f}  |  "
          f"χ_E^quad = {cmp['chi_E_quad']:.4f}  |  "
          f"χ_t_cycle = {cmp['chi_t_cycle']:.4f}")
    print(f"Null z(shuffle) = {cmp['null_models']['shuffle']['z_score']:.2f}  |  "
          f"z(isotropy) = {cmp['null_models']['isotropy']['z_score']:.2f}")
    print(cmp["comparison"]["verdict"])
    print(hol["hol_E_support"]["verdict"])
    print(f"Hol_E-Schätzung = {hol['hol_E_estimate']:.4f}")
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
