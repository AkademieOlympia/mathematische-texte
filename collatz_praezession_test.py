#!/usr/bin/env python3
"""Explorativer Datentest: corr(I(Q), Π(Q)) auf Primvierlingsfenstern.

Modell: collatz_praezession_info.tex
- I(Q) = log(1 + g) mit g aus Spannenregel span = 10 + 12·g (Integrationsstrom)
- Π(Q) = 4/log(p₁) unter uniformer Drift ε(p)=1/log p (witness.py-Kandidat)

Kein Collatz-Anspruch, keine Behauptung über asymptotische Primzahlstruktur.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Literal

from wolfram import (
    build_mod12_integration_stream,
    count_eabc_family_gaps,
    run_automaton_fast_on_stream,
)

DatasetKind = Literal["integration_stream", "marked"]


@dataclass(frozen=True, slots=True)
class QuadrupletRecord:
    kind: DatasetKind
    signature: str
    primes: tuple[int, int, int, int]
    span: int
    g: int
    information: float
    information_log_span: float
    perihel_proxy: float
    projection_proxy: float | None


@dataclass(frozen=True, slots=True)
class CorrelationSummary:
    n: int
    pearson_i_pi: float | None
    spearman_i_pi: float | None
    pearson_i_projection: float | None
    g_min: int
    g_max: int
    g_unique: int
    interpretation: str


def g_from_span(span: int) -> int:
    """g ≥ 0 aus span = 10 + 12·g; wirft bei ungültiger Spanne."""
    if span < 10 or (span - 10) % 12 != 0:
        raise ValueError(f"span={span} passt nicht zu 10+12g")
    return (span - 10) // 12


def information_content(g: int) -> float:
    """I(Q) = log(1 + g) — diskrete Leiter-Information."""
    return math.log1p(g)


def information_log_span(span: int) -> float:
    """Alternative I(Q) = log(span)."""
    return math.log(span)


def perihel_proxy_uniform(p_start: int) -> float:
    """Π-Proxy: 4·ε mit ε(p)=1/log p ⇒ Π = 4/log p (witness.py-Kandidat)."""
    if p_start <= 1:
        return 0.0
    return 4.0 / math.log(p_start)


def projection_proxy_marked(p_start: int) -> float | None:
    """|S| aus Projektionszeuge für p ≡ 5,11 (mod 12); sonst None."""
    residue = p_start % 12
    if residue not in (5, 11):
        return None
    try:
        from susy_fourlinge.witness import quadruplet_witness
    except ImportError:
        return None
    witness = quadruplet_witness(p_start, t=0.0)
    return abs(witness.start_witness.sum_xz)


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    mx = statistics.mean(xs)
    my = statistics.mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den_x = math.sqrt(sum((x - mx) ** 2 for x in xs))
    den_y = math.sqrt(sum((y - my) ** 2 for y in ys))
    if den_x == 0.0 or den_y == 0.0:
        return None
    return num / (den_x * den_y)


def rank_data(values: list[float]) -> list[float]:
    """Durchschnittsrang bei Bindungen (Spearman)."""
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i
        while j + 1 < len(indexed) and indexed[j + 1][1] == indexed[i][1]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[indexed[k][0]] = avg_rank
        i = j + 1
    return ranks


def spearman(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    return pearson(rank_data(xs), rank_data(ys))


def enumerate_integration_quadruplets(limit: int) -> list[QuadrupletRecord]:
    stream = build_mod12_integration_stream(limit)
    result = run_automaton_fast_on_stream(
        stream, max_n=limit, quadruplet_limit=0
    )
    records: list[QuadrupletRecord] = []
    for quad in result.quadruplets:
        g, _ = count_eabc_family_gaps(quad.primes)
        span = quad.primes[3] - quad.primes[0]
        p_start = quad.primes[0]
        records.append(
            QuadrupletRecord(
                kind="integration_stream",
                signature=quad.signature,
                primes=quad.primes,
                span=span,
                g=g,
                information=information_content(g),
                information_log_span=information_log_span(span),
                perihel_proxy=perihel_proxy_uniform(p_start),
                projection_proxy=projection_proxy_marked(p_start),
            )
        )
    return records


def enumerate_marked_quadruplets(limit: int) -> list[QuadrupletRecord]:
    """Markierte Q(p)=(p,p+2,p+6,p+8); span=8, g_leiter=0 — Kontrollstichprobe."""
    try:
        from susy_fourlinge.witness import PrimeSieve
    except ImportError:
        return []

    sieve = PrimeSieve(limit)
    records: list[QuadrupletRecord] = []
    for p in sieve.iter_quadruplet_starts():
        primes = (p, p + 2, p + 6, p + 8)
        span = 8
        g = 0
        word = "ABCE" if p % 12 == 5 else "CEAB"
        records.append(
            QuadrupletRecord(
                kind="marked",
                signature=word.lower(),
                primes=primes,
                span=span,
                g=g,
                information=information_content(g),
                information_log_span=information_log_span(span),
                perihel_proxy=perihel_proxy_uniform(p),
                projection_proxy=projection_proxy_marked(p),
            )
        )
    return records


def summarize_correlation(records: list[QuadrupletRecord]) -> CorrelationSummary:
    xs = [r.information for r in records]
    ys = [r.perihel_proxy for r in records]
    proj_pairs = [
        (r.information, r.projection_proxy)
        for r in records
        if r.projection_proxy is not None
    ]
    gs = [r.g for r in records]
    g_unique = len(set(gs))
    r_pearson = pearson(xs, ys)
    r_spearman = spearman(xs, ys)
    r_proj = (
        pearson([a for a, _ in proj_pairs], [b for _, b in proj_pairs])
        if len(proj_pairs) >= 2
        else None
    )

    if g_unique <= 1:
        interpretation = (
            "I(Q) ist konstant (g hat keine Variation); "
            "Korrelation mit Π ist nicht aussagekräftig."
        )
    elif r_pearson is None:
        interpretation = "Korrelation nicht berechenbar (degenerierte Varianz)."
    elif abs(r_pearson) < 0.1:
        interpretation = (
            f"Keine signifikante lineare Kopplung: Pearson r≈{r_pearson:.4f}. "
            "Das Modell Π=βK wird durch diese Proxy-Wahl nicht gestützt."
        )
    elif abs(r_pearson) < 0.3:
        interpretation = (
            f"Schwache Korrelation (Pearson r≈{r_pearson:.4f}); "
            "explorativ, kein Beweis der Kopplung I→Π."
        )
    else:
        interpretation = (
            f"Mittlere bis starke Korrelation (Pearson r≈{r_pearson:.4f}); "
            "weiterer Test nötig — kausal nicht folgert."
        )

    return CorrelationSummary(
        n=len(records),
        pearson_i_pi=r_pearson,
        spearman_i_pi=r_spearman,
        pearson_i_projection=r_proj,
        g_min=min(gs) if gs else 0,
        g_max=max(gs) if gs else 0,
        g_unique=g_unique,
        interpretation=interpretation,
    )


def run_analysis(limit: int) -> dict:
    integration = enumerate_integration_quadruplets(limit)
    marked = enumerate_marked_quadruplets(limit)
    summary_int = summarize_correlation(integration)
    summary_marked = summarize_correlation(marked) if marked else None

    return {
        "limit": limit,
        "model": {
            "I": "log(1+g)",
            "g_from": "span=10+12g via count_eabc_family_gaps (wolfram.py)",
            "Pi_proxy": "4/log(p_start) with epsilon=1/log p",
            "reference": "collatz_praezession_info.tex",
        },
        "integration_stream": {
            "count": len(integration),
            "correlation": asdict(summary_int),
            "g_histogram": _histogram([r.g for r in integration]),
        },
        "marked_quadruplets": {
            "count": len(marked),
            "correlation": asdict(summary_marked) if summary_marked else None,
            "note": (
                "Markierte Vierlinge haben stets g=0 und span=8; "
                "I(Q)=0 konstant — nur Kontrollstichprobe."
            ),
        },
    }


def _histogram(values: list[int]) -> dict[str, int]:
    hist: dict[str, int] = {}
    for v in values:
        key = str(v)
        hist[key] = hist.get(key, 0) + 1
    return dict(sorted(hist.items(), key=lambda kv: int(kv[0])))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="corr(I, Π) auf Primvierlingsfenstern (explorativ)"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=10**6,
        help="Obergrenze für Integrationsstrom (Standard: 10^6)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("collatz_praezession_test.json"),
        help="JSON-Ausgabedatei",
    )
    args = parser.parse_args()

    report = run_analysis(args.limit)
    args.output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    corr = report["integration_stream"]["correlation"]
    print(f"Limit: {args.limit:,}")
    print(f"Integrationsstrom-Vierlinge: {corr['n']}")
    print(f"g: min={corr['g_min']}, max={corr['g_max']}, unique={corr['g_unique']}")
    print(f"Pearson corr(I, Π): {corr['pearson_i_pi']}")
    print(f"Spearman corr(I, Π): {corr['spearman_i_pi']}")
    print(f"Pearson corr(I, |S|): {corr['pearson_i_projection']}")
    print()
    print(corr["interpretation"])
    print(f"\nJSON: {args.output}")


if __name__ == "__main__":
    main()
