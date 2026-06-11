"""
Differenzmuster zwischen aufeinanderfolgenden Vierlingen (4-Tupel, z. B. Primzahlen)
sowie **Mod-12-Flavor** (E/A/B/C) und die **Kleinsche Vierergruppe** für
``delta_signature`` (gruppenwertiges „Delta“ zwischen zwei Vierling-Signaturen).

    delta_patterns = []
    for Q1, Q2 in zip(quadruplets, quadruplets[1:]):
        delta = group_diff(Q2, Q1)
        delta_patterns.append(delta)

    n_distinct = count_unique(delta_patterns)  # Anzahl verschiedener Muster
    # oder: histogram = dict(Counter(delta_patterns))

    Gruppen-Signaturen (Klein) zwischen benachbarten Vierlingen::

        deltas = consecutive_delta_signatures(quadruplets)

    Häufigkeiten mit ``collections.Counter``::

        from collections import Counter
        cnt = Counter(deltas)  # oder: counter_signature_deltas(quadruplets)

    Shannon-Entropie der Verteilung der Deltas::

        H = shannon_entropy_counts(cnt)  # H = -sum p log p (nats)

    Übergangsgraph ``d_i \\to d_{i+1}`` (braucht ``networkx``)::

        G = delta_transition_digraph(deltas)
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Hashable, Iterable, Mapping, Sequence
from typing import Literal

import numpy as np

Flavor = Literal["E", "A", "B", "C"]

# Kleinsche Vierergruppe auf den Flavor-Labels (E neutrales Element)
GROUP_MUL: dict[tuple[Flavor, Flavor], Flavor] = {
    ("E", "E"): "E",
    ("E", "A"): "A",
    ("E", "B"): "B",
    ("E", "C"): "C",
    ("A", "E"): "A",
    ("A", "A"): "E",
    ("A", "B"): "C",
    ("A", "C"): "B",
    ("B", "E"): "B",
    ("B", "A"): "C",
    ("B", "B"): "E",
    ("B", "C"): "A",
    ("C", "E"): "C",
    ("C", "A"): "B",
    ("C", "B"): "A",
    ("C", "C"): "E",
}

# Alias zum Lesen wie im Nutzer-Snippet
group_mul = GROUP_MUL


def mod12_class(p: int) -> Flavor | None:
    """Restklasse von ``p`` mod 12 in der EABC-Zuordnung (1→E, 5→A, 7→B, 11→C)."""
    r = int(p) % 12
    if r == 1:
        return "E"
    if r == 5:
        return "A"
    if r == 7:
        return "B"
    if r == 11:
        return "C"
    return None


def inv(x: Flavor) -> Flavor:
    """Inverses im Klein-Modell: jedes Element ist **selbstinvers**."""
    return x


def quadruplet_signature(p: int) -> tuple[Flavor | None, Flavor | None, Flavor | None, Flavor | None]:
    """
    Signatur des „dichten“ Musters (Start ``p``, Abstände 0, 2, 6, 8) in E/A/B/C.

    Für ein echtes Primzahl-Vierling-Tupel ``(p, p+2, p+6, p+8)`` sind alle
    Einträge typischerweise gesetzt; sonst kann ``None`` vorkommen.
    """
    return tuple(mod12_class(int(p) + d) for d in (0, 2, 6, 8))


def signature_from_primes(Q: Sequence[int]) -> tuple[Flavor | None, Flavor | None, Flavor | None, Flavor | None]:
    """Signatur aus vier Primzahlen (oder beliebigen Integers) ``(p0, p2, p6, p8)``."""
    Q = tuple(map(int, Q))
    if len(Q) != 4:
        raise ValueError("Q muss Länge 4 haben")
    return tuple(mod12_class(x) for x in Q)


FlavorSet = frozenset({"E", "A", "B", "C"})


def delta_signature(
    Q1: Sequence[Flavor],
    Q2: Sequence[Flavor],
) -> tuple[Flavor, ...]:
    """
    Komponentenweise Gruppenoperation: für ``x`` aus ``Q1``, ``y`` aus ``Q2`` wird
    ``GROUP_MUL[(y, inv(x))] = GROUP_MUL[(y, x)]`` gebildet (Selbstinversität).

    ``Q1``, ``Q2`` sind Tupel aus ``'E'|'A'|'B'|'C'`` der Länge 4 (Signaturen).
    """
    t1 = tuple(Q1)
    t2 = tuple(Q2)
    if len(t1) != len(t2):
        raise ValueError("Q1 und Q2 gleiche Länge nötig")
    out: list[Flavor] = []
    for x, y in zip(t1, t2):
        if x not in FlavorSet or y not in FlavorSet:
            raise ValueError(
                f"delta_signature: erwartete Flavors in {sorted(FlavorSet)}, "
                f"bekam {x!r}, {y!r}"
            )
        xf: Flavor = x  # durch FlavorSet geprüft
        yf: Flavor = y
        out.append(GROUP_MUL[(yf, inv(xf))])
    return tuple(out)


def _signature_from_quadruplet_entry(
    q: Sequence[int] | int,
) -> tuple[Flavor | None, Flavor | None, Flavor | None, Flavor | None]:
    """
    Ein Listeneintrag: **Startzahl** ``p`` → ``quadruplet_signature(p)``,
    oder **4-Tupel** → ``signature_from_primes``.
    """
    if isinstance(q, (list | tuple)):
        return signature_from_primes(q)
    if isinstance(q, int) and not isinstance(q, bool):
        return quadruplet_signature(q)
    raise TypeError(
        "Quadruplet-Eintrag muss int (Start) oder Sequenz der Länge 4 sein, "
        f"nicht {type(q).__name__!r}"
    )


def consecutive_delta_signatures(
    quadruplets: Sequence[Sequence[int] | int],
) -> list[tuple[Flavor, ...]]:
    """
    Entspricht (korrekt für 4-Tupel **oder** nur Startzahl)::

        deltas = []
        for i in range(len(quadruplets) - 1):
            Q1 = _signature_from_quadruplet_entry(quadruplets[i])
            Q2 = _signature_from_quadruplet_entry(quadruplets[i + 1])
            deltas.append(delta_signature(Q1, Q2))

    Hinweis: ``quadruplet_signature`` erwartet einen **einen** Start ``p``;
    ist ``quadruplets[i]`` ein Tupel ``(p, p+2, p+6, p+8)``, wird intern
    ``signature_from_primes`` verwendet.

    Zählen der Vorkommen: ``from collections import Counter`` und ``cnt = Counter(deltas)``
    — oder ``counter_signature_deltas(quadruplets)``.
    """
    n = len(quadruplets)
    if n < 2:
        return []
    deltas: list[tuple[Flavor, ...]] = []
    for i in range(n - 1):
        sig1 = _signature_from_quadruplet_entry(quadruplets[i])
        sig2 = _signature_from_quadruplet_entry(quadruplets[i + 1])
        if any(v is None for v in sig1) or any(v is None for v in sig2):
            raise ValueError(
                f"Ungültige mod12-Signatur bei Index {i} oder {i + 1}: {sig1!r}, {sig2!r}"
            )
        deltas.append(delta_signature(sig1, sig2))
    return deltas


def counter_signature_deltas(
    quadruplets: Sequence[Sequence[int] | int],
) -> Counter[tuple[Flavor, ...]]:
    """
    ``Counter(consecutive_delta_signatures(quadruplets))`` — wie::

        from collections import Counter
        deltas = consecutive_delta_signatures(quadruplets)
        cnt = Counter(deltas)

    Als ``dict``: ``dict(cnt)`` oder ``count_patterns(deltas)``.
    """
    return Counter(consecutive_delta_signatures(quadruplets))


def group_diff(Q2: Sequence[int], Q1: Sequence[int]) -> tuple[int, ...]:
    """
    Komponentenweise Differenz ``Q2 - Q1`` (gleiche Länge).

    Für Vierlinge ``Q_k = (p, p+2, p+6, p+8)`` beschreibt das z. B. den „Sprung“
    zum nächsten Vierling in jeder Komponente.
    """
    if len(Q2) != len(Q1):
        raise ValueError("Q1 und Q2 müssen gleiche Länge haben")
    return tuple(int(b) - int(a) for a, b in zip(Q1, Q2))


def iter_consecutive_pairs(quadruplets: Sequence[Sequence[int]]):
    """Paare ``(Q_i, Q_{i+1})`` aus einer Liste von Tupeln."""
    q = list(quadruplets)
    for i in range(len(q) - 1):
        yield tuple(map(int, q[i])), tuple(map(int, q[i + 1]))


def deltas_from_quadruplets(quadruplets: Sequence[Sequence[int]]) -> list[tuple[int, ...]]:
    """Alle ``group_diff(Q2, Q1)`` für aufeinanderfolgende Vierlinge."""
    return [group_diff(Q2, Q1) for Q1, Q2 in iter_consecutive_pairs(quadruplets)]


def count_unique(patterns: Iterable[Hashable]) -> int:
    """Anzahl **verschiedener** Muster (Kardinalität der Menge)."""
    return len(set(patterns))


def count_patterns(patterns: Iterable[Hashable]) -> dict[Hashable, int]:
    """Häufigkeit jedes Musters — ``dict(Counter(...))``, Analogon zu ``Counter(deltas)``."""
    return dict(Counter(patterns))


def shannon_entropy_counts(
    cnt: Mapping[Hashable, int] | Counter,
    *,
    base: float | None = None,
) -> float:
    """
    Shannon-Entropie ``H = -\\sum_i p_i \\log p_i`` aus **Zähl**daten (``Counter``/``dict``).

    Entspricht (mit ``p`` nur für positive Häufigkeiten)::

        p = np.array(list(cnt.values())) / sum(cnt.values())
        H = -np.sum(p * np.log(p))

    - ``base is None``: natürlicher Log (Einheit **nat**).
    - ``base=2``: Logarithmus zur Basis 2 (Einheit **bit**).

    Leeres ``cnt`` oder Summe 0 → ``0.0``.
    """
    values = np.asarray(list(cnt.values()), dtype=float)
    s = float(values.sum())
    if s <= 0:
        return 0.0
    p = values / s
    p = p[p > 0]
    if base is None:
        return float(-np.sum(p * np.log(p)))
    b = float(base)
    if b <= 0 or b == 1.0:
        raise ValueError("base muss > 0 und ungleich 1 sein")
    return float(-np.sum(p * np.log(p) / np.log(b)))


def delta_transition_digraph(deltas: Sequence[Hashable]):
    """
    Gerichteter **Übergangs**graph der Delta-Folge: für jedes aufeinanderfolgende
    Paar ``(deltas[i], deltas[i+1])`` eine Kante; Kantenattribut ``weight`` = Vielfachheit.

    Erfordert **networkx** (``pip install networkx``).

    **Nicht** so zählen: ``G[d] = G.get(d, 0) + 1`` — ein ``nx.DiGraph`` ist kein
    Häulichkeits-``dict``. Für Vorkommen einzelner Deltas: ``Counter(deltas)``;
    für den Übergangsgraphen diese Funktion.
    """
    try:
        import networkx as nx
    except ImportError as exc:
        raise ImportError(
            "delta_transition_digraph benötigt networkx (z. B. pip install networkx)."
        ) from exc

    seq = list(deltas)
    G = nx.DiGraph()
    if not seq:
        return G
    for d in seq:
        G.add_node(d)
    for (u, v), w in Counter(zip(seq, seq[1:])).items():
        G.add_edge(u, v, weight=int(w))
    return G
