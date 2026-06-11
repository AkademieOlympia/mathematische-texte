#!/usr/bin/env python3
"""
Leiter-DFA (ladder finite automaton) for E-ABC span prediction.

Formalises the mod-12 residue-class ladder between consecutive quadruplet
primes.  Each composite integrable slot (n ≡ 1,5,7,11 mod 12, n composite)
skipped along a transition contributes +12 to the quadruplet span:

    span_predicted = 10 + 12 · g
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from enum import Enum, auto

# ── Shared constants (aligned with wolfram.py) ───────────────────────────────

UNIT_GROUP = frozenset({1, 5, 7, 11})
FAMILY_LETTERS = ("a", "b", "c", "e")
RESIDUE_TO_LETTER = {1: "a", 5: "b", 7: "c", 11: "e"}
LETTER_TO_RESIDUE = {letter: residue for residue, letter in RESIDUE_TO_LETTER.items()}
EABC_LABELS = {"a": "A", "b": "B", "c": "C", "e": "E"}
EABC_RESIDUES = {"A": 1, "B": 5, "C": 7, "E": 11}

QUADRUPLET_SIGNATURES: dict[str, tuple[int, ...]] = {
    "abce": (1, 5, 7, 11),
    "ceab": (7, 11, 1, 5),
}
SIGNATURE_PREFIXES: dict[str, tuple[int, ...]] = {
    sig: pattern[:3] for sig, pattern in QUADRUPLET_SIGNATURES.items()
}
NEXT_RESIDUE_FOR_SIGNATURE = {
    sig: pattern[3] for sig, pattern in QUADRUPLET_SIGNATURES.items()
}

MINIMAL_SPAN = 10
SPAN_STEP = 12

# ── Primality ────────────────────────────────────────────────────────────────


def is_prime(n: int) -> bool:
    """Trial-division primality test (suitable for demo-scale numbers)."""
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def sieve_primes(limit: int) -> list[int]:
    """Return all primes ≤ limit via Eratosthenes."""
    if limit < 2:
        return []
    is_p = bytearray(b"\x01") * (limit + 1)
    is_p[0] = is_p[1] = 0
    root = int(limit**0.5)
    for i in range(2, root + 1):
        if is_p[i]:
            start = i * i
            is_p[start : limit + 1 : i] = b"\x00" * ((limit - start) // i + 1)
    return [i for i in range(2, limit + 1) if is_p[i]]


def build_mod12_stream(max_n: int) -> list[tuple[int, int]]:
    """Integrable primes (mod 12 ∈ units) up to max_n."""
    return [(p, p % 12) for p in sieve_primes(max_n) if p % 12 in UNIT_GROUP]


def first_residue_slot_after(prime: int, target_residue: int) -> int:
    """Smallest n > prime with n ≡ target_residue (mod 12)."""
    offset = (target_residue - prime % 12) % 12
    if offset == 0:
        offset = 12
    return prime + offset


def smallest_prime_after_with_residue(prime: int, residue: int) -> int:
    """Smallest prime strictly greater than prime with p ≡ residue (mod 12)."""
    candidate = first_residue_slot_after(prime, residue)
    while not is_prime(candidate):
        candidate += 12
    return candidate


# ── Leiter-DFA ───────────────────────────────────────────────────────────────


class LeiterState(Enum):
    """Finite-state machine states for the ladder walk."""

    IDLE = auto()
    ON_LADDER = auto()


@dataclass(frozen=True)
class LadderSlot:
    """One position visited on the mod-12 ladder."""

    position: int
    is_composite: bool
    accepted: bool  # composite slot that counts toward g


@dataclass(frozen=True)
class LeiterScanResult:
    """Output of a single ladder scan between two quadruplet primes."""

    family: str
    residue: int
    p_low: int
    p_high: int
    slots: tuple[LadderSlot, ...]
    composite_count: int
    trace: tuple[tuple[LeiterState, int], ...]


@dataclass(frozen=True)
class LeiterDFA:
    """
    Ladder finite automaton for one E-ABC family (A≡1, B≡5, C≡7, E≡11 mod 12).

    States
    ------
    IDLE      – waiting for scan(p_low, p_high)
    ON_LADDER – walking positions n ≡ residue (mod 12), step +12

    Transitions
    -----------
    IDLE --start--> ON_LADDER at first slot > p_low on this ladder
    ON_LADDER --+12--> ON_LADDER while n < p_high
    ON_LADDER --accept--> output composite slot (counts toward g)
    ON_LADDER --halt--> when n ≥ p_high
    """

    family: str
    residue: int

    def scan(self, p_low: int, p_high: int) -> LeiterScanResult:
        """
        Walk the residue-class ladder from p_low toward p_high.

        Counts composite slots strictly between p_low and p_high on the
        ladder n ≡ self.residue (mod 12).
        """
        slots: list[LadderSlot] = []
        trace: list[tuple[LeiterState, int]] = [(LeiterState.IDLE, p_low)]

        pos = first_residue_slot_after(p_low, self.residue)
        trace.append((LeiterState.ON_LADDER, pos))

        while pos < p_high:
            composite = not is_prime(pos)
            slots.append(LadderSlot(position=pos, is_composite=composite, accepted=composite))
            pos += 12
            if pos < p_high:
                trace.append((LeiterState.ON_LADDER, pos))

        return LeiterScanResult(
            family=self.family,
            residue=self.residue,
            p_low=p_low,
            p_high=p_high,
            slots=tuple(slots),
            composite_count=sum(1 for s in slots if s.accepted),
            trace=tuple(trace),
        )

    def state_diagram_ascii(self) -> str:
        """ASCII state diagram for this family's Leiter-DFA."""
        return "\n".join(
            [
                f"  Leiter-DFA Familie {self.family} (≡{self.residue} mod 12)",
                "",
                "       ┌─────────┐",
                "       │  IDLE   │",
                "       └────┬────┘",
                f"            │ start(p_low, p_high): erster Slot > p_low",
                f"            │ mit n ≡ {self.residue} (mod 12)",
                "            ▼",
                "    ┌───────────────┐",
                f"    │  ON_LADDER    │◄──┐",
                f"    │  pos ≡ {self.residue:<2}       │   │ +12 solange pos < p_high",
                "    └───────┬───────┘───┘",
                "            │ pos zusammengesetzt → AKZEPT (g+=1)",
                "            │ pos ≥ p_high        → HALT",
                "            ▼",
                "         [ENDE]",
            ]
        )


LEITER_DFAS: dict[str, LeiterDFA] = {
    fam: LeiterDFA(family=fam, residue=EABC_RESIDUES[fam]) for fam in ("A", "B", "C", "E")
}


# ── Quadruplet analysis ──────────────────────────────────────────────────────


@dataclass(frozen=True)
class TransitionReport:
    """One quadruplet transition analysed by the appropriate Leiter-DFA."""

    index: int
    p_from: int
    p_to: int
    target_family: str
    scan: LeiterScanResult
    g_step: int


@dataclass(frozen=True)
class QuadrupletReport:
    """Full span/g analysis for a prime quadruplet."""

    signature: str
    primes: tuple[int, int, int, int]
    residues: tuple[int, int, int, int]
    transitions: tuple[TransitionReport, ...]
    g_total: int
    span_predicted: int
    span_actual: int
    span_hit: bool


def detect_signature(residues: tuple[int, ...]) -> str | None:
    for signature, pattern in QUADRUPLET_SIGNATURES.items():
        if residues == pattern:
            return signature
    return None


def analyze_quadruplet(
    primes: tuple[int, int, int, int],
    signature: str | None = None,
) -> QuadrupletReport:
    """
    Analyse all three transitions of a quadruplet with the Leiter-DFA model.

    For each transition p_i → p_{i+1}, the target family is determined by
    p_{i+1} mod 12; the corresponding Leiter-DFA scans the ladder and counts
    composite slots.  span_predicted = 10 + 12·g.
    """
    residues = tuple(p % 12 for p in primes)
    if signature is None:
        signature = detect_signature(residues)
        if signature is None:
            raise ValueError(f"Keine gültige Signatur für Reste {residues}")

    transitions: list[TransitionReport] = []
    g_total = 0

    for i in range(3):
        p_from, p_to = primes[i], primes[i + 1]
        letter = RESIDUE_TO_LETTER[p_to % 12]
        family = EABC_LABELS[letter]
        dfa = LEITER_DFAS[family]
        scan = dfa.scan(p_from, p_to)
        g_total += scan.composite_count
        transitions.append(
            TransitionReport(
                index=i + 1,
                p_from=p_from,
                p_to=p_to,
                target_family=family,
                scan=scan,
                g_step=scan.composite_count,
            )
        )

    span_actual = primes[3] - primes[0]
    span_predicted = MINIMAL_SPAN + SPAN_STEP * g_total
    return QuadrupletReport(
        signature=signature,
        primes=primes,
        residues=residues,
        transitions=tuple(transitions),
        g_total=g_total,
        span_predicted=span_predicted,
        span_actual=span_actual,
        span_hit=span_predicted == span_actual,
    )


# ── Signature forecast from 3-prefix ─────────────────────────────────────────


@dataclass(frozen=True)
class SignatureForecast:
    """Prediction from three consecutive integrable-prime residues."""

    signature: str
    prefix_residues: tuple[int, int, int]
    prefix_primes: tuple[int, int, int]
    required_family: str
    required_residue: int
    candidate_p4: int
    span_predicted: int
    g_total: int


def predict_from_prefix(
    prefix_primes: tuple[int, int, int],
) -> SignatureForecast | None:
    """
    Given three consecutive integrable primes, predict signature, 4th family,
    candidate p₄, span, and g.
    """
    prefix_residues = tuple(p % 12 for p in prefix_primes)
    predicted_sig: str | None = None
    for signature, pref in SIGNATURE_PREFIXES.items():
        if prefix_residues == pref:
            predicted_sig = signature
            break
    if predicted_sig is None:
        return None

    required_residue = NEXT_RESIDUE_FOR_SIGNATURE[predicted_sig]
    required_family = EABC_LABELS[RESIDUE_TO_LETTER[required_residue]]
    candidate_p4 = smallest_prime_after_with_residue(prefix_primes[2], required_residue)
    quad = (prefix_primes[0], prefix_primes[1], prefix_primes[2], candidate_p4)
    report = analyze_quadruplet(quad, signature=predicted_sig)
    return SignatureForecast(
        signature=predicted_sig,
        prefix_residues=prefix_residues,
        prefix_primes=prefix_primes,
        required_family=required_family,
        required_residue=required_residue,
        candidate_p4=candidate_p4,
        span_predicted=report.span_predicted,
        g_total=report.g_total,
    )


# ── Quadruplet discovery (batch verification) ────────────────────────────────


@dataclass(frozen=True)
class PrimeQuadruplet:
    signature: str
    primes: tuple[int, int, int, int]
    residues: tuple[int, int, int, int]
    span: int
    stream_index: int


def scan_quadruplets(stream: list[tuple[int, int]], limit: int = 0) -> list[PrimeQuadruplet]:
    """Find consecutive abce/ceab quadruplets in integration-stream order."""
    hits: list[PrimeQuadruplet] = []
    seen: set[tuple[int, int, int, int]] = set()
    for i in range(len(stream) - 3):
        window = stream[i : i + 4]
        residues = tuple(r for _, r in window)
        for signature, pattern in QUADRUPLET_SIGNATURES.items():
            if residues != pattern:
                continue
            primes = tuple(p for p, _ in window)
            if primes in seen:
                break
            seen.add(primes)
            hits.append(
                PrimeQuadruplet(
                    signature=signature,
                    primes=primes,
                    residues=residues,
                    span=primes[3] - primes[0],
                    stream_index=i + 4,
                )
            )
            if limit > 0 and len(hits) >= limit:
                return hits
            break
    return hits


def load_quadruplets_from_wolfram(limit: int = 100) -> list[PrimeQuadruplet]:
    """Import first `limit` quadruplets via wolfram.run_automaton."""
    from wolfram import run_automaton  # noqa: WPS433 — optional integration path

    result = run_automaton(max_n=50_000, quadruplet_limit=limit, fast_mode=True)
    return [
        PrimeQuadruplet(
            signature=q.signature,
            primes=q.primes,
            residues=q.residues,
            span=q.span,
            stream_index=q.stream_index,
        )
        for q in result.quadruplets[:limit]
    ]


def discover_quadruplets(limit: int = 100, *, from_wolfram: bool = False) -> list[PrimeQuadruplet]:
    if from_wolfram:
        return load_quadruplets_from_wolfram(limit)
    stream = build_mod12_stream(50_000)
    return scan_quadruplets(stream, limit=limit)


# ── German output formatters ─────────────────────────────────────────────────


def format_leiter_diagrams() -> str:
    lines = ["── Leiter-DFA Zustandsdiagramme (pro Familie) ──", ""]
    for fam in ("A", "B", "C", "E"):
        lines.append(LEITER_DFAS[fam].state_diagram_ascii())
        lines.append("")
    return "\n".join(lines)


def format_transition_scan(tr: TransitionReport) -> str:
    lines = [
        f"  Übergang {tr.index}: p{tr.index}={tr.p_from} → p{tr.index + 1}={tr.p_to}  "
        f"(Ziel-Familie {tr.target_family}, g={tr.g_step})",
    ]
    if not tr.scan.slots:
        lines.append("    Leiter leer (keine Zwischen-Slots)")
        return "\n".join(lines)

    slot_parts: list[str] = []
    for slot in tr.scan.slots:
        mark = " ◆ zusammengesetzt" if slot.accepted else "   (Primzahl)"
        slot_parts.append(f"    n={slot.position}{mark}")
    lines.extend(slot_parts)
    return "\n".join(lines)


def format_quadruplet_report(report: QuadrupletReport, *, index: int | None = None) -> str:
    idx = f"#{index}  " if index is not None else ""
    p = report.primes
    lines = [
        f"── Vierling {idx}({p[0]}, {p[1]}, {p[2]}, {p[3]}) ──",
        f"  Signatur: {report.signature}  Reste: {report.residues}",
        f"  E-ABC:    {''.join(EABC_LABELS[RESIDUE_TO_LETTER[r]] for r in report.residues)}",
        "",
        "  Leiter-Scans:",
    ]
    for tr in report.transitions:
        lines.append(format_transition_scan(tr))
        lines.append("")

    hit = "✓" if report.span_hit else "✗"
    lines += [
        "  g-Aufschlüsselung:",
        f"    g₁={report.transitions[0].g_step}, "
        f"g₂={report.transitions[1].g_step}, "
        f"g₃={report.transitions[2].g_step}  →  g={report.g_total}",
        "",
        f"  Spanne: vorhergesagt = {MINIMAL_SPAN} + 12·{report.g_total} "
        f"= {report.span_predicted}",
        f"          tatsächlich  = {report.span_actual}  [{hit}]",
    ]
    return "\n".join(lines)


def format_forecast(fc: SignatureForecast) -> str:
    p = fc.prefix_primes
    return "\n".join(
        [
            "── Signatur-Prognose aus 3-Präfix ──",
            f"  Präfix-Reste: {fc.prefix_residues}  →  Signatur {fc.signature}",
            f"  Präfix-Primzahlen: ({p[0]}, {p[1]}, {p[2]})",
            f"  4. Familie: {fc.required_family} (Rest ≡ {fc.required_residue} mod 12)",
            f"  Kandidat p₄: {fc.candidate_p4}",
            f"  Spanne: {fc.span_predicted}  (g={fc.g_total})",
        ]
    )


def format_forward_forecasts_from_wolfram(*, limit: int = 100) -> str:
    """Batch-Vorwärtsmetrik via wolfram.evaluate_forward_forecasts (gleiche Ausgabe)."""
    from wolfram import (  # noqa: WPS433 — shared integration path
        evaluate_forward_forecasts,
        format_forecast_section,
        run_automaton,
    )

    result = run_automaton(max_n=50_000, quadruplet_limit=limit, fast_mode=True)
    report = evaluate_forward_forecasts(result)
    lines = [
        "── Vorwärtsprognosen (wolfram.evaluate_forward_forecasts) ──",
        format_forecast_section(report, verbose=True),
    ]
    return "\n".join(lines)


def format_verification(quadruplets: list[PrimeQuadruplet], spot_checks: tuple[int, ...]) -> str:
    hits = 0
    total = len(quadruplets)
    for q in quadruplets:
        report = analyze_quadruplet(q.primes, signature=q.signature)
        if report.span_hit:
            hits += 1

    lines = [
        "── Verifikation (100 Vierlinge) ──",
        f"  Spannen-Prognose: {hits}/{total} korrekt ({100.0 * hits / total:.1f} %)",
        "",
        "  Stichproben:",
    ]
    for idx in spot_checks:
        if idx < 1 or idx > len(quadruplets):
            continue
        q = quadruplets[idx - 1]
        report = analyze_quadruplet(q.primes, signature=q.signature)
        mark = "✓" if report.span_hit else "✗"
        lines.append(
            f"    #{idx:3d}  {q.signature}  {q.primes}  "
            f"span={q.span}  pred={report.span_predicted}  g={report.g_total}  {mark}"
        )
    return "\n".join(lines)


# ── Demo / CLI ───────────────────────────────────────────────────────────────

DEMO_QUADRUPLETS = {
    6: (109, 113, 127, 131),
    24: (2803, 2819, 2833, 2837),
    87: (25471, 25523, 25537, 25541),
    100: (31387, 31391, 31393, 31397),
}


def run_demo(*, from_wolfram: bool = False) -> None:
    print()
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  Leiter-DFA · E-ABC Spanne-Prognose  (span = 10 + 12·g)         ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    print()
    print("  Familien:  A ≡ 1,  B ≡ 5,  C ≡ 7,  E ≡ 11  (mod 12)")
    print("  Regel:     span = 10 + 12 · g")
    print("             g = Σ zusammengesetzter E-ABC-Slots auf den Leitern")
    print()
    print(format_leiter_diagrams())

    for idx in (6, 87, 100):
        primes = DEMO_QUADRUPLETS[idx]
        report = analyze_quadruplet(primes)
        print(format_quadruplet_report(report, index=idx))
        print()

    q6 = DEMO_QUADRUPLETS[6]
    fc = predict_from_prefix((q6[0], q6[1], q6[2]))
    if fc:
        print(format_forecast(fc))
        print()

    print(format_forward_forecasts_from_wolfram(limit=100))
    print()

    quads = discover_quadruplets(limit=100, from_wolfram=from_wolfram)
    print(format_verification(quads, spot_checks=(6, 24, 87, 100)))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Leiter-DFA für E-ABC Spanne-Prognose (span = 10 + 12·g)",
    )
    parser.add_argument(
        "--quadruplet",
        type=int,
        metavar="N",
        help="Analysiere Vierling #N aus dem Integrations-Stream (1-basiert)",
    )
    parser.add_argument(
        "--primes",
        type=int,
        nargs=4,
        metavar="P",
        help="Analysiere explizite vier Primzahlen p₁ p₂ p₃ p₄",
    )
    parser.add_argument(
        "--from-wolfram",
        action="store_true",
        help="Vierlinge via wolfram.run_automaton laden (sonst lokaler Sieb-Scan)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)

    if args.primes:
        primes = tuple(args.primes)
        report = analyze_quadruplet(primes)
        print(format_quadruplet_report(report))
        return

    if args.quadruplet:
        quads = discover_quadruplets(limit=args.quadruplet, from_wolfram=args.from_wolfram)
        if args.quadruplet > len(quads):
            print(f"Fehler: nur {len(quads)} Vierlinge gefunden.", file=sys.stderr)
            sys.exit(1)
        q = quads[args.quadruplet - 1]
        report = analyze_quadruplet(q.primes, signature=q.signature)
        print(format_quadruplet_report(report, index=args.quadruplet))
        return

    run_demo(from_wolfram=args.from_wolfram)


if __name__ == "__main__":
    main()
