#!/usr/bin/env python3
"""
Suche nach einem arithmetischen Peres-Mermin/Kochen-Specker-Quadrat.

Eingabeformat: eine CSV-Kantenliste mit mindestens

    context_id, orbit_id, context_sign

Jeder Kontext muss genau drei verschiedene Orbits enthalten. Das Vorzeichen
`context_sign` ist das algebraische Produkt des Kontextes (+1 oder -1), nicht
ein fixer Einzelwert des Orbits. Optional koennen Zeilen ueber Phase/Residuen
vorgefiltert werden, wenn die entsprechenden Spalten vorhanden sind.

Beispiel:

    python3 arithmetisches_ks_suchprotokoll.py daten.csv \
        --phase-column phase_k --phase-k 1 3 \
        --residue-column residue --modulus 210 --residue 5 11 101 191
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


POSITIVE_PATTERN = (+1, +1, +1)
COLUMN_PATTERN = (+1, +1, -1)


@dataclass(frozen=True)
class Context:
    context_id: str
    sign: int
    orbits: frozenset[str]


@dataclass(frozen=True)
class MagicSquare:
    rows: tuple[Context, Context, Context]
    columns: tuple[Context, Context, Context]
    grid: tuple[tuple[str, str, str], tuple[str, str, str], tuple[str, str, str]]


def parse_sign(value: str) -> int:
    text = str(value).strip()
    if text in {"+", "+1", "1", "plus", "pos", "positive"}:
        return +1
    if text in {"-", "-1", "minus", "neg", "negative"}:
        return -1
    raise ValueError(f"Ungueltiges Vorzeichen: {value!r}")


def parse_int_set(values: Iterable[str] | None) -> set[int] | None:
    if not values:
        return None
    return {int(v) for v in values}


def row_passes_filters(
    row: dict[str, str],
    *,
    phase_column: str | None,
    allowed_phases: set[int] | None,
    residue_column: str | None,
    modulus: int | None,
    allowed_residues: set[int] | None,
) -> bool:
    if phase_column and allowed_phases is not None:
        if phase_column not in row or int(row[phase_column]) % 4 not in allowed_phases:
            return False

    if residue_column and allowed_residues is not None:
        if residue_column not in row:
            return False
        residue = int(row[residue_column])
        if modulus:
            residue %= modulus
        if residue not in allowed_residues:
            return False

    return True


def load_contexts(
    path: Path,
    *,
    context_column: str,
    orbit_column: str,
    sign_column: str,
    phase_column: str | None,
    allowed_phases: set[int] | None,
    residue_column: str | None,
    modulus: int | None,
    allowed_residues: set[int] | None,
) -> list[Context]:
    grouped_orbits: dict[str, set[str]] = defaultdict(set)
    grouped_signs: dict[str, set[int]] = defaultdict(set)

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {context_column, orbit_column, sign_column}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Fehlende CSV-Spalten: {sorted(missing)}")

        for row in reader:
            if not row_passes_filters(
                row,
                phase_column=phase_column,
                allowed_phases=allowed_phases,
                residue_column=residue_column,
                modulus=modulus,
                allowed_residues=allowed_residues,
            ):
                continue

            context_id = row[context_column].strip()
            orbit_id = row[orbit_column].strip()
            grouped_orbits[context_id].add(orbit_id)
            grouped_signs[context_id].add(parse_sign(row[sign_column]))

    contexts: list[Context] = []
    for context_id, orbits in grouped_orbits.items():
        signs = grouped_signs[context_id]
        if len(orbits) != 3 or len(signs) != 1:
            continue
        contexts.append(
            Context(
                context_id=context_id,
                sign=next(iter(signs)),
                orbits=frozenset(orbits),
            )
        )

    return contexts


def intersection_grid(
    rows: tuple[Context, Context, Context],
    columns: tuple[Context, Context, Context],
) -> tuple[tuple[str, str, str], tuple[str, str, str], tuple[str, str, str]] | None:
    row_union: set[str] = set().union(*(row.orbits for row in rows))
    col_union: set[str] = set().union(*(col.orbits for col in columns))

    if len(row_union) != 9 or row_union != col_union:
        return None

    grid_rows: list[tuple[str, str, str]] = []
    seen: set[str] = set()
    for row in rows:
        cells: list[str] = []
        for column in columns:
            overlap = row.orbits.intersection(column.orbits)
            if len(overlap) != 1:
                return None
            orbit = next(iter(overlap))
            if orbit in seen:
                return None
            seen.add(orbit)
            cells.append(orbit)
        grid_rows.append(tuple(cells))

    if len(seen) != 9:
        return None

    return tuple(grid_rows)  # type: ignore[return-value]


def find_magic_squares(contexts: list[Context], max_solutions: int) -> Iterable[MagicSquare]:
    positive = [ctx for ctx in contexts if ctx.sign == +1]
    negative = [ctx for ctx in contexts if ctx.sign == -1]

    for rows in itertools.combinations(positive, 3):
        row_tuple = tuple(rows)  # type: ignore[assignment]
        if tuple(ctx.sign for ctx in row_tuple) != POSITIVE_PATTERN:
            continue

        for positive_columns in itertools.combinations(positive, 2):
            for negative_column in negative:
                columns = (*positive_columns, negative_column)
                for ordered_columns in itertools.permutations(columns, 3):
                    if tuple(ctx.sign for ctx in ordered_columns) != COLUMN_PATTERN:
                        continue

                    grid = intersection_grid(row_tuple, ordered_columns)  # type: ignore[arg-type]
                    if grid is None:
                        continue

                    yield MagicSquare(
                        rows=row_tuple,  # type: ignore[arg-type]
                        columns=ordered_columns,  # type: ignore[arg-type]
                        grid=grid,
                    )
                    max_solutions -= 1
                    if max_solutions <= 0:
                        return


def square_to_dict(square: MagicSquare) -> dict[str, object]:
    return {
        "row_contexts": [ctx.context_id for ctx in square.rows],
        "column_contexts": [ctx.context_id for ctx in square.columns],
        "row_signs": [ctx.sign for ctx in square.rows],
        "column_signs": [ctx.sign for ctx in square.columns],
        "grid": square.grid,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Findet 3x3-Kontextquadrate mit Zeilenprodukt +++ und Spaltenprodukt ++-."
    )
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--context-column", default="context_id")
    parser.add_argument("--orbit-column", default="orbit_id")
    parser.add_argument("--sign-column", default="context_sign")
    parser.add_argument("--phase-column", default=None)
    parser.add_argument("--phase-k", nargs="*", default=None, help="Erlaubte k in lambda=exp(2*pi*i*k/4)")
    parser.add_argument("--residue-column", default=None)
    parser.add_argument("--modulus", type=int, default=None)
    parser.add_argument("--residue", nargs="*", default=None, help="Erlaubte Residuen, optional modulo --modulus")
    parser.add_argument("--max-solutions", type=int, default=10)
    args = parser.parse_args()

    contexts = load_contexts(
        args.csv_path,
        context_column=args.context_column,
        orbit_column=args.orbit_column,
        sign_column=args.sign_column,
        phase_column=args.phase_column,
        allowed_phases=parse_int_set(args.phase_k),
        residue_column=args.residue_column,
        modulus=args.modulus,
        allowed_residues=parse_int_set(args.residue),
    )

    solutions = list(find_magic_squares(contexts, max_solutions=args.max_solutions))
    result = {
        "input": str(args.csv_path),
        "valid_contexts": len(contexts),
        "solutions": [square_to_dict(square) for square in solutions],
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
