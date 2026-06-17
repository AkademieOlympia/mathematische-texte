#!/usr/bin/env python3
"""
Gerichteter EABC-Übergangsgraph und Zyklus-Holonomie-Kandidat.

Kanonsiche Theorie: collatz_eabc_transport.md

Transport T: κ(p_n) → κ(p_{n+1}) entlang der Primfolge (p > 3).
Fundamentales Objekt: gerichtete Kante τ(p_n) = (κ(p_n), κ(p_{n+1})).

Vergleich:
  χ_E(N)        — Prim-Vierlinge (collatz_eabc_holonomie.md)
  χ_trans(N)    — ABCE/CEAB-Fenster auf 4 aufeinanderfolgenden Prim-Klassen
  χ_t_cycle(N)  — vier t- bzw. t^{-1}-alignierte Schritte

Ausführung:
    python3 collatz_eabc_transition_graph.py
    python3 collatz_eabc_transition_graph.py --max-p 100000 --null-trials 2000
"""

from __future__ import annotations

import argparse
import json
import random
import statistics
from math import isqrt
from pathlib import Path
from typing import Any

from collatz_eabc_holonomie_test import chi_E, enumerate_quadruplets
from eabc_from_lean import EClass, class_of, t

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_transition_graph.json"

LABELS = ("E", "A", "B", "C")
IDX = {label: i for i, label in enumerate(LABELS)}
ABCE_WORD = "ABCE"
CEAB_WORD = "CEAB"
T_INV = {t(c): c for c in EClass}


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
    abce = count_word_windows(classes, ABCE_WORD)
    ceab = count_word_windows(classes, CEAB_WORD)
    t_fwd = count_t_cycle_windows(classes, forward=True)
    t_inv = count_t_cycle_windows(classes, forward=False)
    return {
        "abce_windows": abce,
        "ceab_windows": ceab,
        "chi_trans": chi_from_bias(abce, ceab),
        "t_forward_windows": t_fwd,
        "t_inverse_windows": t_inv,
        "chi_t_cycle": chi_from_bias(t_fwd, t_inv),
        "formula_trans": "(#ABCE_windows - #CEAB_windows) / (#ABCE + #CEAB)",
        "formula_t_cycle": "(#t_forward_4cycles - #t_inverse_4cycles) / (#fwd + #inv)",
    }


def chi_quadruplet_report(max_p: int) -> dict[str, Any]:
    quads = enumerate_quadruplets(max_p)
    chi = chi_E(max_p, quads)
    return {
        "quadruplet_count": chi["quadruplet_count"],
        "abce_count": chi["abce_count"],
        "ceab_count": chi["ceab_count"],
        "chi_E": chi["chi_E"],
        "formula": chi["formula"],
    }


def shuffle_null_chi(
    classes: list[str],
    trials: int = 1000,
    seed: int = 42,
) -> dict[str, Any]:
    """Permutiere Klassenfolge, erhalte Marginal-Häufigkeiten."""
    rng = random.Random(seed)
    base = chi_transport(classes)
    chi_trans_samples: list[float] = []
    chi_t_samples: list[float] = []
    for _ in range(trials):
        perm = classes[:]
        rng.shuffle(perm)
        null = chi_transport(perm)
        chi_trans_samples.append(null["chi_trans"])
        chi_t_samples.append(null["chi_t_cycle"])
    return {
        "trials": trials,
        "seed": seed,
        "observed_chi_trans": base["chi_trans"],
        "observed_chi_t_cycle": base["chi_t_cycle"],
        "null_chi_trans_mean": statistics.mean(chi_trans_samples),
        "null_chi_trans_std": statistics.pstdev(chi_trans_samples) if trials > 1 else 0.0,
        "null_chi_t_cycle_mean": statistics.mean(chi_t_samples),
        "null_chi_t_cycle_std": statistics.pstdev(chi_t_samples) if trials > 1 else 0.0,
        "z_chi_trans": _z_score(base["chi_trans"], chi_trans_samples),
        "z_chi_t_cycle": _z_score(base["chi_t_cycle"], chi_t_samples),
    }


def _z_score(observed: float, samples: list[float]) -> float:
    if len(samples) < 2:
        return 0.0
    mu = statistics.mean(samples)
    sigma = statistics.pstdev(samples)
    if sigma == 0:
        return 0.0
    return (observed - mu) / sigma


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


def chi_transport_vs_quadruplet(max_p: int, null_trials: int = 1000, seed: int = 42) -> dict[str, Any]:
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    trans = chi_transport(classes)
    quad = chi_quadruplet_report(max_p)
    null = shuffle_null_chi(classes, trials=null_trials, seed=seed)
    same_sign = (trans["chi_trans"] >= 0) == (quad["chi_E"] >= 0)
    close = abs(trans["chi_trans"] - quad["chi_E"]) < 0.05
    return {
        "max_p": max_p,
        "chi_trans": trans["chi_trans"],
        "chi_t_cycle": trans["chi_t_cycle"],
        "chi_E": quad["chi_E"],
        "transport_detail": trans,
        "quadruplet_detail": quad,
        "null_model": null,
        "comparison": {
            "same_sign": same_sign,
            "difference": trans["chi_trans"] - quad["chi_E"],
            "close_within_0_05": close,
            "verdict": (
                "χ_trans und χ_E sind verwandte chirale Observablen in verschiedenen Trägern "
                "(Primfolge-Fenster vs. Prim-Vierlinge) — nicht identisch. "
                f"Beobachtet: χ_trans={trans['chi_trans']:.4f}, χ_E={quad['chi_E']:.4f}, "
                f"Δ={trans['chi_trans'] - quad['chi_E']:.4f}."
            ),
        },
    }


def hol_E_estimates(limits: list[int], null_trials: int = 500) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for lim in limits:
        cmp = chi_transport_vs_quadruplet(lim, null_trials=null_trials)
        rows.append(
            {
                "limit": lim,
                "chi_E": cmp["chi_E"],
                "chi_trans": cmp["chi_trans"],
                "chi_t_cycle": cmp["chi_t_cycle"],
                "z_chi_trans": cmp["null_model"]["z_chi_trans"],
            }
        )
    return {
        "limits": limits,
        "series": rows,
        "hol_E_candidate": rows[-1]["chi_E"] if rows else 0.0,
        "hol_trans_candidate": rows[-1]["chi_trans"] if rows else 0.0,
        "note": "Hol_E = lim χ_E(N); empirische Stützung via letzte Stufe — kein Beweis des Grenzwerts.",
    }


def run(
    max_p: int = 100_000,
    null_trials: int = 1000,
    output: Path = DEFAULT_OUTPUT,
    seed: int = 42,
) -> dict[str, Any]:
    matrix = transition_matrix_report(max_p)
    comparison = chi_transport_vs_quadruplet(max_p, null_trials=null_trials, seed=seed)
    limits = sorted({min(max_p, x) for x in (1_000, 5_000, 10_000, 50_000, max_p)})
    hol = hol_E_estimates(limits, null_trials=min(null_trials, 300))

    report: dict[str, Any] = {
        "meta": {
            "module": "collatz_eabc_transition_graph.py",
            "theory": "collatz_eabc_transport.md",
            "holonomy_hierarchy": "collatz_eabc_holonomie.md §5–6",
            "max_p": max_p,
            "null_trials": null_trials,
        },
        "transition_matrix": matrix,
        "cycle_holonomy": comparison,
        "hol_E_estimates": hol,
        "epistemic_labels": {
            "transition_matrix": "Definition + Experiment",
            "chi_trans": "Definition + Experiment",
            "chi_E": "Definition + Experiment",
            "shuffle_null": "Experiment",
            "hol_E_limit": "Definition (Grenzwert offen)",
        },
    }
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    report["output_path"] = str(output)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC-Übergangsgraph und Zyklus-Holonomie")
    parser.add_argument("--max-p", type=int, default=100_000, help="Obergrenze für Primzahlen")
    parser.add_argument("--null-trials", type=int, default=1000, help="Shuffle-Null-Durchläufe")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    report = run(max_p=args.max_p, null_trials=args.null_trials, output=args.output, seed=args.seed)
    cmp = report["cycle_holonomy"]
    mat = report["transition_matrix"]
    hol = report["hol_E_estimates"]
    print("=== EABC-Übergangsgraph ===")
    print(f"Primzahlen >3 bis {args.max_p}: {mat['prime_count']}")
    print(f"t-Bias Kanten: {mat['edge_summary']['t_bias']:.4f}  "
          f"(fwd {mat['edge_summary']['t_forward_fraction']:.4f}, "
          f"inv {mat['edge_summary']['t_inverse_fraction']:.4f})")
    print(f"χ_trans = {cmp['chi_trans']:.4f}  |  χ_t_cycle = {cmp['chi_t_cycle']:.4f}  |  "
          f"χ_E = {cmp['chi_E']:.4f}")
    print(f"Null z(χ_trans) = {cmp['null_model']['z_chi_trans']:.2f}")
    print(cmp["comparison"]["verdict"])
    print(f"Hol_E-Kandidat (letzte Stufe) = {hol['hol_E_candidate']:.4f}")
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
