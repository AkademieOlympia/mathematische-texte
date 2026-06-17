#!/usr/bin/env python3
"""
Gerichteter EABC-Übergangsgraph, Pfadorientierung (4-Block) und Zyklus-Holonomie (5-Block).

Kanonsiche Theorie:
  collatz_eabc_zyklus_holonomie.md  — χ_path (4-Block), χ_hol^(5) (5-Block), Hol_E
  collatz_eabc_transport.md         — G_E, Transport T_n
  collatz_eabc_holonomie.md         — χ_E^quad auf Prim-Vierlingen (Vergleich)

χ_path(N)   = Σ χ_path(Q_n^(4)) / #{χ_path≠0}   — Pfadorientierung ABCE/CEAB
χ_hol^(5)(N) = Σ Ω^(5)(Q_n^(5)) / #{Ω^(5)≠0}    — geschlossener Zyklus ABCEA/CEABC
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
ABCEA_WORD = "ABCEA"
CEABC_WORD = "CEABC"
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


def omega_path(word: str) -> Omega:
    """4-Block Pfadorientierung: ABCE=+1, CEAB=-1, 0 sonst."""
    if word == ABCE_WORD:
        return 1
    if word == CEAB_WORD:
        return -1
    return 0


# Legacy alias
omega_window = omega_path


def omega_5(word: str) -> Omega:
    """5-Block Zyklus-Holonomie: ABCEA=+1, CEABC=-1, 0 sonst."""
    if word == ABCEA_WORD:
        return 1
    if word == CEABC_WORD:
        return -1
    return 0


def sliding_windows(classes: list[str], width: int = 4) -> list[dict[str, Any]]:
    """Gleitfenster fester Breite auf der Klassenfolge."""
    if width < 1:
        return []
    omega_fn = omega_5 if width == 5 else omega_path
    rows: list[dict[str, Any]] = []
    for i in range(len(classes) - width + 1):
        word = "".join(classes[i : i + width])
        rows.append({"index": i, "word": word, "omega": omega_fn(word)})
    return rows


def _chi_from_windows(
    windows: list[dict[str, Any]],
    *,
    carrier: str,
    formula: str,
    positive_key: str,
    negative_key: str,
) -> dict[str, Any]:
    omega_sum = sum(w["omega"] for w in windows)
    positive = sum(1 for w in windows if w["omega"] == 1)
    negative = sum(1 for w in windows if w["omega"] == -1)
    nonzero = positive + negative
    chi = omega_sum / nonzero if nonzero else 0.0
    return {
        "window_count": len(windows),
        positive_key: positive,
        negative_key: negative,
        "omega_sum": omega_sum,
        "nonzero_windows": nonzero,
        "chi": chi,
        "formula": formula,
        "carrier": carrier,
    }


def chi_path_sliding(classes: list[str]) -> dict[str, Any]:
    """
    χ_path(N) = Σ χ_path(Q_n^(4)) / #{χ_path≠0} auf Primfolge-Gleitfenstern.

    Kanonisch: collatz_eabc_zyklus_holonomie.md §4.
    """
    windows = sliding_windows(classes, width=4)
    report = _chi_from_windows(
        windows,
        carrier="prime_sequence_sliding_windows_4",
        formula="sum(chi_path(Q_n^(4))) / #{chi_path != 0}",
        positive_key="abce_windows",
        negative_key="ceab_windows",
    )
    report["chi_path"] = report["chi"]
    report["chi_E"] = report["chi"]  # legacy alias
    return report


def chi_hol_sliding(classes: list[str]) -> dict[str, Any]:
    """
    χ_hol^(5)(N) = Σ Ω^(5)(Q_n^(5)) / #{Ω^(5)≠0} auf Primfolge-Gleitfenstern.

    Kanonisch: collatz_eabc_zyklus_holonomie.md §5.
    """
    windows = sliding_windows(classes, width=5)
    report = _chi_from_windows(
        windows,
        carrier="prime_sequence_sliding_windows_5",
        formula="sum(Omega^(5)(Q_n^(5))) / #{Omega^(5) != 0}",
        positive_key="abcea_windows",
        negative_key="ceabc_windows",
    )
    report["chi_hol"] = report["chi"]
    report["chi_hol_5"] = report["chi"]
    return report


# Legacy alias
chi_E_sliding = chi_path_sliding


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
    """Legacy-Alias: χ_trans = χ_path auf Gleitfenstern (identische Formel)."""
    path = chi_path_sliding(classes)
    t_fwd = count_t_cycle_windows(classes, forward=True)
    t_inv = count_t_cycle_windows(classes, forward=False)
    return {
        "abce_windows": path["abce_windows"],
        "ceab_windows": path["ceab_windows"],
        "chi_path": path["chi_path"],
        "chi_E": path["chi_path"],
        "chi_trans": path["chi_path"],
        "t_forward_windows": t_fwd,
        "t_inverse_windows": t_inv,
        "chi_t_cycle": chi_from_bias(t_fwd, t_inv),
        "formula": path["formula"],
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


def _null_chi(
    classes: list[str],
    chi_fn,
    trials: int,
    seed: int,
    null_type: str,
    use_isotropy: bool,
) -> dict[str, Any]:
    rng = random.Random(seed)
    base = chi_fn(classes)
    chi_key = "chi_path" if chi_fn is chi_path_sliding else "chi_hol"
    observed = base[chi_key]
    samples: list[float] = []
    for _ in range(trials):
        if use_isotropy:
            perm_classes = _apply_label_permutation(classes, _random_label_permutation(rng))
        else:
            perm_classes = classes[:]
            rng.shuffle(perm_classes)
        samples.append(chi_fn(perm_classes)[chi_key])
    return {
        "null_type": null_type,
        "trials": trials,
        "seed": seed,
        f"observed_{chi_key}": observed,
        "null_mean": statistics.mean(samples),
        "null_std": statistics.pstdev(samples) if trials > 1 else 0.0,
        "z_score": _z_score(observed, samples),
    }


def shuffle_null_chi_path(
    classes: list[str],
    trials: int = 1000,
    seed: int = 42,
) -> dict[str, Any]:
    """Marginal-Shuffle-Null für χ_path (4-Block)."""
    return _null_chi(classes, chi_path_sliding, trials, seed, "marginal_shuffle", False)


def isotropy_null_chi_path(
    classes: list[str],
    trials: int = 1000,
    seed: int = 42,
) -> dict[str, Any]:
    """Isotropie-Null für χ_path (4-Block)."""
    return _null_chi(classes, chi_path_sliding, trials, seed + 1, "isotropy_relabel", True)


def shuffle_null_chi_hol(
    classes: list[str],
    trials: int = 1000,
    seed: int = 42,
) -> dict[str, Any]:
    """Marginal-Shuffle-Null für χ_hol^(5) (5-Block)."""
    return _null_chi(classes, chi_hol_sliding, trials, seed + 2, "marginal_shuffle", False)


def isotropy_null_chi_hol(
    classes: list[str],
    trials: int = 1000,
    seed: int = 42,
) -> dict[str, Any]:
    """Isotropie-Null für χ_hol^(5) (5-Block)."""
    return _null_chi(classes, chi_hol_sliding, trials, seed + 3, "isotropy_relabel", True)


# Legacy aliases
shuffle_null_chi_E = shuffle_null_chi_path
isotropy_null_chi_E = isotropy_null_chi_path


def hol_E_supported(
    chi_series: list[dict[str, Any]],
    shuffle_null: dict[str, Any],
    isotropy_null: dict[str, Any],
    tol: float = 0.01,
    chi_key: str = "chi_hol",
) -> dict[str, Any]:
    """Bewertung: Hol_E≠0 vs. stabile Fluktuationen vs. Null (5-Block)."""
    if not chi_series:
        return {"supported": False, "reason": "empty series"}
    last_chi = chi_series[-1][chi_key]
    nonzero_limit = abs(last_chi) > tol
    stable = len(chi_series) >= 2 and all(
        abs(chi_series[i][chi_key] - chi_series[-1][chi_key]) < 0.05
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
            f"(χ_hol^(5)={last_chi:.4f}); "
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


def chi_path_vs_hol(
    max_p: int,
    null_trials: int = 1000,
    seed: int = 42,
) -> dict[str, Any]:
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    path = chi_path_sliding(classes)
    hol = chi_hol_sliding(classes)
    quad = chi_E_quadruplet_report(max_p)
    trans = chi_transport(classes)
    null_path = {
        "shuffle": shuffle_null_chi_path(classes, trials=null_trials, seed=seed),
        "isotropy": isotropy_null_chi_path(classes, trials=null_trials, seed=seed),
    }
    null_hol = {
        "shuffle": shuffle_null_chi_hol(classes, trials=null_trials, seed=seed),
        "isotropy": isotropy_null_chi_hol(classes, trials=null_trials, seed=seed),
    }
    path_stronger = abs(path["chi_path"]) > abs(hol["chi_hol"])
    return {
        "max_p": max_p,
        "chi_path": path["chi_path"],
        "chi_hol": hol["chi_hol"],
        "chi_E_quad": quad["chi_E_quad"],
        "chi_t_cycle": trans["chi_t_cycle"],
        "path_detail": path,
        "hol_detail": hol,
        "quadruplet_detail": quad,
        "transport_detail": trans,
        "null_models": {
            "path": null_path,
            "hol": null_hol,
        },
        "comparison": {
            "same_sign_path_hol": (path["chi_path"] >= 0) == (hol["chi_hol"] >= 0),
            "difference_path_hol": path["chi_path"] - hol["chi_hol"],
            "abs_path_vs_abs_hol": {
                "abs_chi_path": abs(path["chi_path"]),
                "abs_chi_hol": abs(hol["chi_hol"]),
                "path_stronger_signal": path_stronger,
            },
            "verdict": (
                "χ_path (4-Block-Pfad) und χ_hol^(5) (5-Block-Zyklus) sind verwandte Observablen "
                "mit unterschiedlicher geometrischer Lesart — Pfad vs. geschlossene Holonomie. "
                f"Beobachtet: χ_path={path['chi_path']:.4f}, "
                f"χ_hol^(5)={hol['chi_hol']:.4f}, "
                f"|χ_path|={'>' if path_stronger else '≤'} |χ_hol|."
            ),
        },
    }


# Backward-compatible aliases
chi_sliding_vs_quadruplet = chi_path_vs_hol
chi_transport_vs_quadruplet = chi_path_vs_hol


def hol_E_estimates(
    limits: list[int],
    null_trials: int = 500,
    seed: int = 42,
) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    last_cmp: dict[str, Any] | None = None
    for lim in limits:
        cmp = chi_path_vs_hol(lim, null_trials=null_trials, seed=seed)
        last_cmp = cmp
        rows.append(
            {
                "limit": lim,
                "chi_path": cmp["chi_path"],
                "chi_hol": cmp["chi_hol"],
                "chi_E_quad": cmp["chi_E_quad"],
                "chi_t_cycle": cmp["chi_t_cycle"],
                "z_shuffle_path": cmp["null_models"]["path"]["shuffle"]["z_score"],
                "z_isotropy_path": cmp["null_models"]["path"]["isotropy"]["z_score"],
                "z_shuffle_hol": cmp["null_models"]["hol"]["shuffle"]["z_score"],
                "z_isotropy_hol": cmp["null_models"]["hol"]["isotropy"]["z_score"],
                # legacy keys
                "chi_E": cmp["chi_path"],
                "z_shuffle": cmp["null_models"]["hol"]["shuffle"]["z_score"],
                "z_isotropy": cmp["null_models"]["hol"]["isotropy"]["z_score"],
            }
        )
    support = hol_E_supported(
        rows,
        last_cmp["null_models"]["hol"]["shuffle"] if last_cmp else {},
        last_cmp["null_models"]["hol"]["isotropy"] if last_cmp else {},
    )
    return {
        "limits": limits,
        "series": rows,
        "hol_E_estimate": rows[-1]["chi_hol"] if rows else 0.0,
        "chi_path_estimate": rows[-1]["chi_path"] if rows else 0.0,
        "hol_E_quad_estimate": rows[-1]["chi_E_quad"] if rows else 0.0,
        "hol_E_support": support,
        "note": (
            "Hol_E = lim χ_hol^(5)(N); empirische Schätzung = letzter χ_hol-Wert — "
            "kein Beweis des Grenzwerts. χ_path ist Pfadorientierung (4-Block), nicht Holonomie."
        ),
    }


def run(
    max_p: int = 100_000,
    null_trials: int = 1000,
    output: Path = DEFAULT_OUTPUT,
    seed: int = 42,
) -> dict[str, Any]:
    matrix = transition_matrix_report(max_p)
    comparison = chi_path_vs_hol(max_p, null_trials=null_trials, seed=seed)
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
            "terminology": {
                "chi_path": "4-Block Pfadorientierung ABCE/CEAB",
                "chi_hol": "5-Block Zyklus-Holonomie ABCEA/CEABC",
            },
        },
        "transition_matrix": matrix,
        "path_vs_holonomy": comparison,
        "cycle_holonomy": comparison,  # legacy key
        "hol_E_estimates": hol,
        "epistemic_labels": {
            "transition_matrix": "Definition + Experiment",
            "chi_path_sliding": "Definition (4-Block Pfad)",
            "chi_hol_sliding": "Definition (5-Block Holonomie)",
            "chi_E_quad": "Definition (Vierlingsträger)",
            "shuffle_null": "Experiment",
            "isotropy_null": "Experiment",
            "hol_E_limit": "Definition (Grenzwert auf χ_hol^(5))",
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
    cmp = report["path_vs_holonomy"]
    mat = report["transition_matrix"]
    hol = report["hol_E_estimates"]
    print("=== EABC Pfad (4-Block) vs. Holonomie (5-Block) ===")
    print(f"Primzahlen >3 bis {args.max_p}: {mat['prime_count']}")
    print(f"χ_path (4-Block) = {cmp['chi_path']:.4f}  |  "
          f"χ_hol^(5) = {cmp['chi_hol']:.4f}  |  "
          f"χ_E^quad = {cmp['chi_E_quad']:.4f}")
    print(f"Null z_path(shuffle) = {cmp['null_models']['path']['shuffle']['z_score']:.2f}  |  "
          f"z_hol(shuffle) = {cmp['null_models']['hol']['shuffle']['z_score']:.2f}")
    print(cmp["comparison"]["verdict"])
    print(hol["hol_E_support"]["verdict"])
    print(f"Hol_E-Schätzung (χ_hol^(5)) = {hol['hol_E_estimate']:.4f}")
    print(f"JSON: {report['output_path']}")


if __name__ == "__main__":
    main()
