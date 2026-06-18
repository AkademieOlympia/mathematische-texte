#!/usr/bin/env python3
"""
EABC-Bell-Ungleichung: P_same auf Transportfenstern, CHSH-Analog, G_E-Vergleich.

Theorie: collatz_eabc_bell_holonomie.md

  σ(n) = 1[X_{n+1}=t(X_n)]                    Kantenlesart
  P_same^win(i,j) auf ABCE-Fenstern           gemeinsamer Träger (Theorem: Summe ≥ 1)
  P_same^marg(i,j)                            marginalisierte Kantenränder
  S_EABC = E(a,b)-E(a,b')+E(a',b)+E(a',b')   CHSH auf Zyklus E-A-B-C

Ausführung:
    python3 collatz_eabc_bell_inequality_test.py
    python3 collatz_eabc_bell_inequality_test.py --max-p 1000000
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

from collatz_eabc_transition_graph import (
    ABCE_WORD,
    classes_from_sequence,
    omega_hol,
    omega_pfad,
    prime_eabc_sequence,
    transition_counts,
    transition_probabilities,
)
from eabc_from_lean import EClass, t

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "collatz_eabc_bell_inequality.json"
THEORY = "collatz_eabc_bell_holonomie.md"
DEFAULT_MAX_P = 1_000_000

LABELS = ("E", "A", "B", "C")
IDX = {label: i for i, label in enumerate(LABELS)}
TRIPLE = ("E", "A", "C")
PAIR_KEYS = (("E", "A"), ("E", "C"), ("A", "C"))
QM_CHSH_MAX = 2 * math.sqrt(2)
CLASSICAL_CHSH_MAX = 2.0


def _valid_class(ch: str) -> bool:
    return ch in IDX


def t_alignment(classes: list[str]) -> list[int]:
    """σ(n) = 1 iff X_{n+1} = t(X_n)."""
    aligned: list[int] = []
    for i in range(len(classes) - 1):
        if not _valid_class(classes[i]) or not _valid_class(classes[i + 1]):
            aligned.append(0)
            continue
        c0 = EClass(classes[i])
        c1 = EClass(classes[i + 1])
        aligned.append(1 if t(c0) is c1 else 0)
    return aligned


def _observable_on_abce_window(
    classes: list[str], align: list[int], start: int, role: str
) -> int | None:
    """O_E, O_A, O_C on ABCE window starting at start (needs 5th class for O_C)."""
    if start + 4 >= len(classes):
        return None
    word = "".join(classes[start : start + 4])
    if word != ABCE_WORD:
        return None
    offsets = {"A": 0, "C": 2, "E": 3}
    off = offsets.get(role)
    if off is None:
        return None
    idx = start + off
    if idx >= len(align):
        return None
    return align[idx]


def p_same_win_report(classes: list[str]) -> dict[str, Any]:
    """
    P_same^win on ABCE windows for triple (E,A,C).
    Theorem: B_win >= 1 always (on windows with valid O_C).
    """
    align = t_alignment(classes)
    abce_starts: list[int] = []
    for i in range(len(classes) - 4):
        if "".join(classes[i : i + 4]) == ABCE_WORD:
            abce_starts.append(i)

    pair_match: dict[tuple[str, str], int] = {k: 0 for k in PAIR_KEYS}
    pair_total = 0

    for start in abce_starts:
        o_e = _observable_on_abce_window(classes, align, start, "E")
        o_a = _observable_on_abce_window(classes, align, start, "A")
        o_c = _observable_on_abce_window(classes, align, start, "C")
        if o_e is None or o_a is None or o_c is None:
            continue
        pair_total += 1
        if o_e == o_a:
            pair_match[("E", "A")] += 1
        if o_e == o_c:
            pair_match[("E", "C")] += 1
        if o_a == o_c:
            pair_match[("A", "C")] += 1

    p_same: dict[str, float | None] = {}
    for key in PAIR_KEYS:
        p_same[f"{key[0]}_{key[1]}"] = (
            pair_match[key] / pair_total if pair_total else None
        )

    b_win = sum(pair_match[k] / pair_total for k in PAIR_KEYS) if pair_total else None
    per_window_min = None
    if pair_total:
        mins = []
        for start in abce_starts:
            o_e = _observable_on_abce_window(classes, align, start, "E")
            o_a = _observable_on_abce_window(classes, align, start, "A")
            o_c = _observable_on_abce_window(classes, align, start, "C")
            if o_e is None or o_a is None or o_c is None:
                continue
            mins.append(int(o_e == o_a) + int(o_e == o_c) + int(o_a == o_c))
        per_window_min = min(mins) if mins else None

    return {
        "abce_window_count": pair_total,
        "abce_starts_requiring_5th_class": len(abce_starts),
        "pair_match_counts": {f"{a}_{b}": pair_match[(a, b)] for a, b in PAIR_KEYS},
        "P_same_win": p_same,
        "B_win": b_win,
        "per_window_min_pair_matches": per_window_min,
        "theorem_B_win_ge_1": (
            b_win is not None and b_win >= 1.0 - 1e-12 and per_window_min is not None
        ),
        "carrier": "ABCE_sliding_windows_with_5th_class",
    }


def _marginal_bit_distribution(
    classes: list[str], align: list[int], start_class: str
) -> dict[str, float]:
    """P(σ=1 | X_n = start_class) over all edges."""
    counts = {0: 0, 1: 0}
    for i in range(len(align)):
        if classes[i] == start_class:
            counts[align[i]] += 1
    total = counts[0] + counts[1]
    if total == 0:
        return {"p0": 0.0, "p1": 0.0, "n": 0}
    return {
        "p0": counts[0] / total,
        "p1": counts[1] / total,
        "n": total,
    }


def p_same_marginal(classes: list[str]) -> dict[str, Any]:
    """Marginal P_same without common window — can violate B >= 1."""
    align = t_alignment(classes)
    dists = {c: _marginal_bit_distribution(classes, align, c) for c in TRIPLE}

    def overlap(i: str, j: str) -> float:
        di, dj = dists[i], dists[j]
        return di["p0"] * dj["p0"] + di["p1"] * dj["p1"]

    p_same = {f"{a}_{b}": overlap(a, b) for a, b in PAIR_KEYS}
    b_marg = sum(p_same.values())

    return {
        "marginals": dists,
        "P_same_marg": p_same,
        "B_marg": b_marg,
        "marginal_B_ge_1": b_marg >= 1.0 - 1e-12,
        "carrier": "marginalized_edge_starts",
    }


def p_same_hol_path(classes: list[str]) -> dict[str, Any]:
    """Coincidence of Pfad vs Holonomie orientation on overlapping indices."""
    match = 0
    total = 0
    for n in range(len(classes) - 4):
        w4 = "".join(classes[n : n + 4])
        w5 = "".join(classes[n : n + 5])
        op = omega_pfad(w4)
        oh = omega_hol(w5)
        if op == 0 or oh == 0:
            continue
        o_pfad = (1 + op) // 2
        o_hol = (1 + oh) // 2
        total += 1
        if o_pfad == o_hol:
            match += 1
    p_hol = match / total if total else None
    return {
        "overlap_count": total,
        "match_count": match,
        "P_same_hol": p_hol,
        "carrier": "Pfad_vs_Holonomie_overlap",
    }


def _binary_from_pfad(classes: list[str], n: int) -> int | None:
    w = "".join(classes[n : n + 4])
    o = omega_pfad(w)
    if o == 0:
        return None
    return (1 + o) // 2


def _binary_from_hol(classes: list[str], n: int) -> int | None:
    w = "".join(classes[n : n + 5])
    o = omega_hol(w)
    if o == 0:
        return None
    return (1 + o) // 2


def _correlation(bits_a: list[int], bits_b: list[int]) -> float | None:
    if not bits_a:
        return None
    agree = sum(1 for a, b in zip(bits_a, bits_b) if a == b)
    return 2 * agree / len(bits_a) - 1


def _chsh_sum(e_ab: float, e_abp: float, e_apb: float, e_apbp: float) -> float:
    """S = E(a,b) - E(a,b') + E(a',b) + E(a',b')."""
    return e_ab - e_abp + e_apb + e_apbp


def chsh_eabc_cycle_report(classes: list[str]) -> dict[str, Any]:
    """
    CHSH on common ABCE carrier (collatz_eabc_bell_holonomie.md §7.2):
      a  = σ(n)           E-edge
      a' = Õ_Pfad(n)      path orientation
      b  = σ(n+3)         C-edge (needs X_{n+4})
      b' = Õ_Hol(n)       holonomy
    """
    align = t_alignment(classes)
    a_bits: list[int] = []
    a_prime: list[int] = []
    b_bits: list[int] = []
    b_prime: list[int] = []

    for n in range(len(classes) - 4):
        if "".join(classes[n : n + 4]) != ABCE_WORD:
            continue
        o_pfad = _binary_from_pfad(classes, n)
        o_hol = _binary_from_hol(classes, n)
        if o_pfad is None or o_hol is None:
            continue
        idx_e = n + 3
        idx_c = n + 2
        if idx_e >= len(align) or idx_c >= len(align):
            continue
        a_bits.append(align[idx_e])
        a_prime.append(o_pfad)
        b_bits.append(align[idx_c])
        b_prime.append(o_hol)

    empty = {
        "mapping": {
            "a": "sigma(n+3) at E-node of ABCE",
            "a_prime": "O_Pfad(n)",
            "b": "sigma(n+2) at C-node of ABCE",
            "b_prime": "O_Hol(n)",
        },
        "E_a_b": None,
        "E_a_bp": None,
        "E_ap_b": None,
        "E_ap_bp": None,
        "S_EABC": None,
        "abs_S_EABC": None,
        "sample_size": 0,
        "classical_bound_2": CLASSICAL_CHSH_MAX,
        "qm_reference_2sqrt2": QM_CHSH_MAX,
        "exceeds_classical_CHSH": None,
        "reaches_qm_reference": None,
        "carrier": "ABCE_common_window_cycle",
    }

    if not a_bits:
        return empty

    e_ab = _correlation(a_bits, b_bits)
    e_abp = _correlation(a_bits, b_prime)
    e_apb = _correlation(a_prime, b_bits)
    e_apbp = _correlation(a_prime, b_prime)
    assert None not in (e_ab, e_abp, e_apb, e_apbp)
    s_val = _chsh_sum(e_ab, e_abp, e_apb, e_apbp)

    return {
        **empty,
        "E_a_b": e_ab,
        "E_a_bp": e_abp,
        "E_ap_b": e_apb,
        "E_ap_bp": e_apbp,
        "S_EABC": s_val,
        "abs_S_EABC": abs(s_val),
        "sample_size": len(a_bits),
        "exceeds_classical_CHSH": abs(s_val) > CLASSICAL_CHSH_MAX + 1e-12,
        "reaches_qm_reference": abs(s_val) > QM_CHSH_MAX - 1e-6,
    }


def chsh_eabc_marginal_report(classes: list[str]) -> dict[str, Any]:
    """CHSH with mismatched carriers (non-factorizable contexts) — may exceed 2."""
    align = t_alignment(classes)

    a_bits, a_prime = [], []
    for n in range(len(classes) - 3):
        if classes[n] not in ("E", "A"):
            continue
        o_pfad = _binary_from_pfad(classes, n)
        if o_pfad is None:
            continue
        a_bits.append(align[n])
        a_prime.append(o_pfad)

    b_bits, b_prime = [], []
    for n in range(len(classes) - 4):
        if classes[n] not in ("B", "C"):
            continue
        o_hol = _binary_from_hol(classes, n)
        if o_hol is None:
            continue
        b_bits.append(align[n])
        b_prime.append(o_hol)

    m = min(len(a_bits), len(b_bits), len(a_prime), len(b_prime))
    if m == 0:
        return {
            "S_EABC_marginal": None,
            "abs_S_EABC_marginal": None,
            "sample_size": 0,
            "carrier": "marginal_mismatched_carriers",
        }

    a_bits, a_prime = a_bits[:m], a_prime[:m]
    b_bits, b_prime = b_bits[:m], b_prime[:m]
    e_ab = _correlation(a_bits, b_bits)
    e_abp = _correlation(a_bits, b_prime)
    e_apb = _correlation(a_prime, b_bits)
    e_apbp = _correlation(a_prime, b_prime)
    assert None not in (e_ab, e_abp, e_apb, e_apbp)
    s_val = _chsh_sum(e_ab, e_abp, e_apb, e_apbp)

    return {
        "E_a_b": e_ab,
        "E_a_bp": e_abp,
        "E_ap_b": e_apb,
        "E_ap_bp": e_apbp,
        "S_EABC_marginal": s_val,
        "abs_S_EABC_marginal": abs(s_val),
        "sample_size": m,
        "exceeds_classical_CHSH": abs(s_val) > CLASSICAL_CHSH_MAX + 1e-12,
        "carrier": "marginal_mismatched_carriers",
    }


# Legacy alias
chsh_eabc_report = chsh_eabc_cycle_report


def transition_graph_comparison(classes: list[str]) -> dict[str, Any]:
    """Compare t-alignment rates with row-stochastic G_E."""
    align = t_alignment(classes)
    counts = transition_counts(classes)
    probs = transition_probabilities(counts)

    t_align_rate: dict[str, float] = {}
    for i, lab in enumerate(LABELS):
        row_sum = sum(counts[i])
        if row_sum == 0:
            t_align_rate[lab] = 0.0
            continue
        t_dst = LABELS[(i + 1) % 4]
        t_align_rate[lab] = counts[i][IDX[t_dst]] / row_sum

    empirical_sigma: dict[str, float] = {}
    for lab in LABELS:
        d = _marginal_bit_distribution(classes, align, lab)
        empirical_sigma[lab] = d["p1"]

    return {
        "t_align_rate_from_G_E": t_align_rate,
        "empirical_sigma_p1": empirical_sigma,
        "transition_probabilities": {
            f"{LABELS[i]}->{LABELS[j]}": probs[i][j] for i in range(4) for j in range(4)
        },
        "t_forward_edge_probs": {
            lab: probs[IDX[lab]][(IDX[lab] + 1) % 4] for lab in LABELS
        },
    }


def run(max_p: int) -> dict[str, Any]:
    seq = prime_eabc_sequence(max_p)
    classes = classes_from_sequence(seq)
    win = p_same_win_report(classes)
    marg = p_same_marginal(classes)
    hol = p_same_hol_path(classes)
    chsh_cycle = chsh_eabc_cycle_report(classes)
    chsh_marg = chsh_eabc_marginal_report(classes)
    ge = transition_graph_comparison(classes)

    return {
        "meta": {
            "module": "collatz_eabc_bell_inequality_test.py",
            "theory": THEORY,
            "max_p": max_p,
            "prime_eabc_count": len(seq),
            "chsh_formula": "S = E(a,b) - E(a,b') + E(a',b) + E(a',b')",
            "labels": {
                "B_win": "Theorem: >= 1 on ABCE windows",
                "B_marg": "Marginal: may be < 1",
                "S_EABC_cycle": "CHSH on common ABCE window, LHV bound |S|<=2",
                "S_EABC_marginal": "CHSH mismatched carriers, not LHV-fair",
            },
        },
        "p_same_win": win,
        "p_same_marginal": marg,
        "p_same_hol_path": hol,
        "chsh_eabc_cycle": chsh_cycle,
        "chsh_eabc_marginal": chsh_marg,
        "chsh_eabc": chsh_cycle,
        "transition_graph": ge,
        "summary": {
            "B_win": win["B_win"],
            "B_win_satisfies_bound": win["theorem_B_win_ge_1"],
            "B_marg": marg["B_marg"],
            "B_marg_satisfies_bound": marg["marginal_B_ge_1"],
            "P_same_hol": hol["P_same_hol"],
            "S_EABC": chsh_cycle["S_EABC"],
            "abs_S_EABC": chsh_cycle["abs_S_EABC"],
            "S_exceeds_2": chsh_cycle.get("exceeds_classical_CHSH"),
            "S_reaches_2sqrt2": chsh_cycle.get("reaches_qm_reference"),
            "abs_S_EABC_marginal": chsh_marg.get("abs_S_EABC_marginal"),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="EABC Bell / CHSH experiment")
    parser.add_argument("--max-p", type=int, default=DEFAULT_MAX_P)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    report = run(args.max_p)
    args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")

    s = report["summary"]
    ch = report["chsh_eabc_cycle"]
    print(f"max_p={args.max_p}  primes>{3}: {report['meta']['prime_eabc_count']}")
    print(f"B_win={s['B_win']:.4f}  satisfies >=1: {s['B_win_satisfies_bound']}")
    print(f"B_marg={s['B_marg']:.4f}  satisfies >=1: {s['B_marg_satisfies_bound']}")
    print(f"P_same_hol={s['P_same_hol']}")
    print(
        f"CHSH cycle: S={s['S_EABC']:.4f}  |S|={s['abs_S_EABC']:.4f}  "
        f"n={ch['sample_size']}  exceeds 2: {s['S_exceeds_2']}"
    )
    print(f"  E(a,b)={ch['E_a_b']:.4f}  E(a,b')={ch['E_a_bp']:.4f}  "
          f"E(a',b)={ch['E_ap_b']:.4f}  E(a',b')={ch['E_ap_bp']:.4f}")
    print(f"  classical bound 2, QM ref 2√2≈{QM_CHSH_MAX:.4f}")
    print(f"written: {args.output}")


if __name__ == "__main__":
    main()
