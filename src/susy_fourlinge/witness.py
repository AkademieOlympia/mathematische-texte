"""Projektionszeuge für markierte Primzahlvierlinge Q(p) = (p, p+2, p+6, p+8)."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

Flavor = Literal["E", "A", "B", "C"]
ChiralityWord = Literal["ABCE", "CEAB"]
Edge = tuple[Flavor, Flavor]

FLAVORS: tuple[Flavor, ...] = ("E", "A", "B", "C")
MONOTONE_CYCLE: tuple[Edge, ...] = (("E", "A"), ("A", "B"), ("B", "C"), ("C", "E"))

FLAVOR_ANGLE: dict[Flavor, float] = {
    "E": math.pi / 6,
    "A": 5 * math.pi / 6,
    "B": 7 * math.pi / 6,
    "C": 11 * math.pi / 6,
}

START_EDGE: dict[ChiralityWord, Edge] = {
    "ABCE": ("A", "B"),
    "CEAB": ("C", "E"),
}

# Unnormierte Basisrichtungen; ||BASE_AB||^2 = ||BASE_CE||^2 = 7.
_BASE_AB = (-math.sqrt(3.0), 1.0, -math.sqrt(3.0))
_BASE_CE = (math.sqrt(3.0), -1.0, math.sqrt(3.0))
_EDGE_BASE: dict[Edge, tuple[float, float, float]] = {
    ("A", "B"): _BASE_AB,
    ("B", "C"): _BASE_AB,
    ("C", "E"): _BASE_CE,
    ("E", "A"): _BASE_CE,
}

_SQRT7 = math.sqrt(7.0)
_VECTOR_NORM = _SQRT7


@dataclass(frozen=True, slots=True)
class Vec3:
    x: float
    y: float
    z: float

    def norm(self) -> float:
        return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def sum_xz(self) -> float:
        return self.x + self.z


@dataclass(frozen=True, slots=True)
class EdgeWitness:
    edge: Edge
    weight: float
    projection: Vec3
    sum_xz: float


@dataclass(frozen=True, slots=True)
class QuadrupletWitness:
    start: int
    word: ChiralityWord
    start_edge: Edge
    start_witness: EdgeWitness
    delta_weight: float


@dataclass(frozen=True, slots=True)
class SieveStats:
    limit: int
    total: int
    abce: int
    ceab: int
    bias: float
    signed_bias: float


@dataclass(frozen=True, slots=True)
class CentroidStats:
    limit: int
    total: int
    abce: int
    ceab: int
    max_centroid_error: float


@dataclass(frozen=True, slots=True)
class QuadrupletEllipseParams:
    """Kanonische Bamberger-Primvierlingsellipse aus der Normalform M+(-4,-2,+2,+4)."""

    a: float = 4.0
    b: float = 2.0
    f: float = 2.0 * math.sqrt(3.0)
    e: float = math.sqrt(3.0) / 2.0
    step_pattern: tuple[int, int, int] = (2, 4, 2)
    offsets: tuple[int, int, int, int] = (-4, -2, 2, 4)


@dataclass(frozen=True, slots=True)
class WitnessEllipseBridge:
    """Brücke Ebene B (Projektionszeuge) ↔ Ebene B′ (Primvierlingsellipse)."""

    start: int
    center: int
    word: ChiralityWord
    ellipse: QuadrupletEllipseParams
    witness: QuadrupletWitness
    rho_pv: float
    sum_xz: float


QUADRUPLET_OFFSETS: tuple[int, int, int, int] = (-4, -2, 2, 4)
STEP_PATTERN: tuple[int, int, int] = (2, 4, 2)
RHO_PV: float = 1.5


def family_vector(flavor: Flavor) -> Vec3:
    """EABC-Tetraeder-Einbettung (Heeger): s = 1/sqrt(3)."""
    s = 1.0 / math.sqrt(3.0)
    if flavor == "E":
        return Vec3(s, s, s)
    if flavor == "A":
        return Vec3(s, -s, -s)
    if flavor == "B":
        return Vec3(-s, s, -s)
    if flavor == "C":
        return Vec3(-s, -s, s)
    raise ValueError(f"Unbekannte Flavor-Klasse: {flavor!r}")


def edge_midpoint_angle(u: Flavor, v: Flavor) -> float:
    """Winkelmitte der beiden Flavor-Positionen auf dem Einheitskreis."""
    return math.atan2(
        math.sin(FLAVOR_ANGLE[u]) + math.sin(FLAVOR_ANGLE[v]),
        math.cos(FLAVOR_ANGLE[u]) + math.cos(FLAVOR_ANGLE[v]),
    )


def delta_weight(t: float) -> float:
    """ΔW(t) = cos t - sin t."""
    return math.cos(t) - math.sin(t)


def _rotate_y(base: tuple[float, float, float], phase: float) -> Vec3:
    x0, y0, z0 = base
    c = math.cos(phase)
    s = math.sin(phase)
    return Vec3(c * x0 + s * z0, y0, -s * x0 + c * z0)


def _edge_weight_core(t: float, phase: float, sign: int) -> float:
    delta = math.cos(phase) - math.sin(phase)
    return 2.0 + sign * delta / 2.0


def edge_weight(
    u: Flavor,
    v: Flavor,
    t: float,
    *,
    sign: int = 1,
    start_edge: bool = False,
) -> float:
    """
    Kantengewicht W_UV(t).

    Monotone Folge E→A→B→C→E: sign=+1, phase = t + Winkelmitte(U,V).
    Startkante ABCE (A,B): sign=+1, phase = t.
    Startkante CEAB (C,E): sign=-1, phase = t.
    """
    phase = t if start_edge else t + edge_midpoint_angle(u, v)
    return _edge_weight_core(t, phase, sign)


def edge_projection(
    u: Flavor,
    v: Flavor,
    t: float,
    *,
    sign: int = 1,
    start_edge: bool = False,
) -> Vec3:
    """Projektionsvektor P_UV(t) mit |P| * W = sqrt(7)/2."""
    edge = (u, v)
    if edge not in _EDGE_BASE:
        raise ValueError(f"Unbekannte gerichtete Kante: {edge!r}")

    weight = edge_weight(u, v, t, sign=sign, start_edge=start_edge)
    phase = t if start_edge else t + edge_midpoint_angle(u, v)
    direction = _rotate_y(_EDGE_BASE[edge], phase)
    # ||direction|| = sqrt(7), daher |P| = sqrt(7)/(2W)
    scale = 1.0 / (2.0 * weight)
    return Vec3(direction.x * scale, direction.y * scale, direction.z * scale)


def chirality_word(p: int) -> ChiralityWord:
    """Klassifikation über den markierten Startpunkt p."""
    residue = int(p) % 12
    if residue == 5:
        return "ABCE"
    if residue == 11:
        return "CEAB"
    raise ValueError(f"Startpunkt {p} hat ungültige Restklasse {residue} mod 12")


def start_edge_sign(word: ChiralityWord) -> int:
    return 1 if word == "ABCE" else -1


def projection_sum_witness(projection: Vec3) -> float:
    """Einfaches Zeuge S = x + z."""
    return projection.sum_xz()


def classify_by_projection(projection: Vec3) -> ChiralityWord:
    """S < 0 ⇒ ABCE, S > 0 ⇒ CEAB."""
    witness = projection_sum_witness(projection)
    if witness < 0:
        return "ABCE"
    if witness > 0:
        return "CEAB"
    raise ValueError("Projektionszeuge S = x + z ist degeneriert (S = 0)")


def edge_witness(
    u: Flavor,
    v: Flavor,
    t: float,
    *,
    sign: int = 1,
    start_edge: bool = False,
) -> EdgeWitness:
    projection = edge_projection(u, v, t, sign=sign, start_edge=start_edge)
    return EdgeWitness(
        edge=(u, v),
        weight=edge_weight(u, v, t, sign=sign, start_edge=start_edge),
        projection=projection,
        sum_xz=projection.sum_xz(),
    )


def monotone_transition_weights(t: float) -> dict[str, float]:
    """Volle monotone Folge E→A→B→C→E."""
    labels = ("W_EA", "W_AB", "W_BC", "W_CE")
    return {
        label: edge_weight(u, v, t, sign=1, start_edge=False)
        for label, (u, v) in zip(labels, MONOTONE_CYCLE)
    }


def monotone_transition_witnesses(t: float) -> dict[str, EdgeWitness]:
    """Projektions-Gate-Zeuge für die volle monotone Folge E→A→B→C→E."""
    labels = ("EA", "AB", "BC", "CE")
    return {
        label: edge_witness(u, v, t, sign=1, start_edge=False)
        for label, (u, v) in zip(labels, MONOTONE_CYCLE)
    }


def quadruplet_witness(p: int, t: float = 0.0) -> QuadrupletWitness:
    """Zeuge für den markierten Vierling mit Startpunkt p bei Parameter t."""
    word = chirality_word(p)
    edge = START_EDGE[word]
    sign = start_edge_sign(word)
    start = edge_witness(*edge, t, sign=sign, start_edge=True)
    return QuadrupletWitness(
        start=p,
        word=word,
        start_edge=edge,
        start_witness=start,
        delta_weight=delta_weight(t),
    )


def _odd_index(n: int) -> int:
    return (n - 3) // 2


def _simple_sieve_odd(limit: int) -> bytearray:
    if limit < 3:
        return bytearray()
    size = _odd_index(limit) + 1
    is_prime_odd = bytearray(b"\x01") * size
    max_p = math.isqrt(limit)
    for i in range(_odd_index(3), _odd_index(max_p) + 1):
        if is_prime_odd[i]:
            p = 2 * i + 3
            start = p * p
            step = 2 * p
            for n in range(start, limit + 1, step):
                is_prime_odd[_odd_index(n)] = 0
    return is_prime_odd


class PrimeSieve:
    """Odd-only Primtabelle für schnelle Vierlingssuche."""

    def __init__(self, limit: int) -> None:
        self.limit = limit
        self._odd = _simple_sieve_odd(limit)

    def is_prime(self, n: int) -> bool:
        if n == 2:
            return True
        if n < 2 or n > self.limit or n % 2 == 0:
            return False
        return bool(self._odd[_odd_index(n)])

    def iter_quadruplet_starts(self) -> list[int]:
        """Alle Startpunkte p mit Primvierling (p, p+2, p+6, p+8)."""
        starts: list[int] = []
        for p in range(5, self.limit - 7, 2):
            if p % 12 not in (5, 11):
                continue
            if all(self.is_prime(q) for q in (p, p + 2, p + 6, p + 8)):
                starts.append(p)
        return starts


def quadruplet_center(p: int) -> int:
    """Symmetrieanker M = p + 4 des Primvierlings Q(p)."""
    return p + 4


def prime_quadruplet_centroid(p: int) -> float:
    """Arithmetischer Schwerpunkt M(Q) = (p + (p+2) + (p+6) + (p+8)) / 4."""
    q = quadruplet_members(p)
    return sum(q) / 4


def centered_prime_quadruplet(p: int) -> tuple[int, int, int, int]:
    """Zentrierte Normalform Q(p) - M(Q) = (-4, -2, 2, 4)."""
    m = quadruplet_center(p)
    return (p - m, p + 2 - m, p + 6 - m, p + 8 - m)


def prime_quadruplet_ellipse_parameters() -> dict[str, float]:
    """Kanonische Ellipsenparameter (a, b, e) und projektives Gewicht rho_PV."""
    a = 4.0
    b = 2.0
    e = math.sqrt(1.0 - (b * b) / (a * a))
    return {
        "a_pv": a,
        "b_pv": b,
        "e_pv": e,
        "rho_pv": RHO_PV,
    }


def canonical_ellipse_params() -> QuadrupletEllipseParams:
    """Ellipsenparameter (a,b,f,e) fest durch die arithmetische Normalform."""
    return QuadrupletEllipseParams()


def kepler_ellipse_point(p: int, theta: float) -> tuple[float, float]:
    """Kepler-Ellipse E_PV(θ) = M + a·cos θ + i·b·sin θ als (Re, Im).

    Chiralität (ABCE vs. CEAB) steuert die π-Verschiebung gemäß `chirality_word`.
    """
    m = float(quadruplet_center(p))
    a, b = 4.0, 2.0
    word = chirality_word(p)
    if word == "CEAB":
        theta += math.pi
    return (m + a * math.cos(theta), b * math.sin(theta))


def kepler_phase_tick(t: int) -> float:
    """Diskreter Phasentick θ(t) = (π/2)·t für t ∈ {0,1,2,3}."""
    if t not in (0, 1, 2, 3):
        raise ValueError(f"Phasentick t muss in {{0,1,2,3}} liegen, nicht {t}")
    return (math.pi / 2.0) * t


def quadruplet_normal_form(p: int) -> tuple[int, tuple[int, int, int, int]]:
    """Liefert (M, (-4,-2,+2,+4)) mit Q(p) = M + Offsets."""
    return quadruplet_center(p), QUADRUPLET_OFFSETS


def quadruplet_members(p: int) -> tuple[int, int, int, int]:
    """Q(p) = (p, p+2, p+6, p+8)."""
    return (p, p + 2, p + 6, p + 8)


def projective_start_weight_ce(t: float) -> float:
    """Projektives Kantengewicht ρ_PV = W_CE(t) an der Startkante (CEAB)."""
    return edge_weight("C", "E", t, sign=-1, start_edge=True)


def witness_ellipse_bridge(p: int, t: float = 0.0) -> WitnessEllipseBridge:
    """Verknüpft kanonische Ellipsenparameter mit dem Projektionszeugen."""
    witness = quadruplet_witness(p, t)
    return WitnessEllipseBridge(
        start=p,
        center=quadruplet_center(p),
        word=witness.word,
        ellipse=canonical_ellipse_params(),
        witness=witness,
        rho_pv=projective_start_weight_ce(t),
        sum_xz=witness.start_witness.sum_xz,
    )


def centroid_error(p: int) -> float:
    """|arithmetischer Schwerpunkt von Q(p) minus M|."""
    members = (p, p + 2, p + 6, p + 8)
    centroid = sum(members) / 4
    return abs(centroid - quadruplet_center(p))


def centroid_statistics(limit: int) -> CentroidStats:
    """Zählt Primvierlinge bis N und prüft den Schwerpunktfehler gegen M = p + 4."""
    sieve = PrimeSieve(limit)
    max_error = 0.0
    abce = 0
    ceab = 0
    for p in sieve.iter_quadruplet_starts():
        max_error = max(max_error, centroid_error(p))
        word = chirality_word(p)
        if word == "ABCE":
            abce += 1
        else:
            ceab += 1
    total = abce + ceab
    return CentroidStats(
        limit=limit,
        total=total,
        abce=abce,
        ceab=ceab,
        max_centroid_error=max_error,
    )


def ideal_phase(m: int) -> float:
    """Ideale EABC-Uhr: θ_m = (π/2)·m."""
    return (math.pi / 2) * m


def perihel_shift(theta_m: float, theta_m4: float) -> float:
    """Π = Θ_{m+4} − Θ_m − 2π."""
    return theta_m4 - theta_m - 2 * math.pi


def uniform_drift_phase(m: int, epsilon: float) -> float:
    """θ_m = m·(π/2 + ε)."""
    return m * (math.pi / 2 + epsilon)


def perihel_shift_uniform_epsilon(epsilon: float, m: int = 0) -> float:
    """Π unter uniformer Drift ε; idealerweise Π = 4ε."""
    return perihel_shift(
        uniform_drift_phase(m, epsilon),
        uniform_drift_phase(m + 4, epsilon),
    )


def drift_epsilon_inv_log(p: int) -> float:
    """Kandidat ε(p) ~ 1/log p (nicht bewiesen)."""
    if p <= 1:
        return 0.0
    return 1.0 / math.log(p)


@dataclass(frozen=True, slots=True)
class PerihelWitnessStats:
    """Numerischer Vergleich Π_real vs Π_random (Evidenz, kein Beweis)."""

    limit: int
    sample_count: int
    real_mean_pi: float
    random_mean_pi: float
    real_std_pi: float
    random_std_pi: float
    ideal_zero_max_abs: float


def test_perihel_real_vs_random(
    limit: int = 10**5,
    *,
    random_count: int | None = None,
) -> PerihelWitnessStats:
    """
    Vergleicht Π unter ε(p)=1/log p für echte Primvierlinge vs. Zufalls-ε.

    Nur numerische Evidenz — keine Behauptung über asymptotische Primzahlstruktur.
    """
    import random
    import statistics

    sieve = PrimeSieve(limit)
    starts = sieve.iter_quadruplet_starts()
    real_pi = [perihel_shift_uniform_epsilon(drift_epsilon_inv_log(p)) for p in starts]
    count = random_count if random_count is not None else len(real_pi)
    random_pi = [
        perihel_shift_uniform_epsilon(random.uniform(-0.05, 0.05) / random.randint(2, 100))
        for _ in range(count)
    ]

    ideal_checks = [
        abs(perihel_shift(ideal_phase(m), ideal_phase(m + 4))) for m in range(64)
    ]
    return PerihelWitnessStats(
        limit=limit,
        sample_count=len(real_pi),
        real_mean_pi=statistics.mean(real_pi) if real_pi else 0.0,
        random_mean_pi=statistics.mean(random_pi) if random_pi else 0.0,
        real_std_pi=statistics.pstdev(real_pi) if len(real_pi) > 1 else 0.0,
        random_std_pi=statistics.pstdev(random_pi) if len(random_pi) > 1 else 0.0,
        ideal_zero_max_abs=max(ideal_checks) if ideal_checks else 0.0,
    )


def sieve_statistics(limit: int) -> SieveStats:
    sieve = PrimeSieve(limit)
    abce = 0
    ceab = 0
    for p in sieve.iter_quadruplet_starts():
        word = chirality_word(p)
        if word == "ABCE":
            abce += 1
        else:
            ceab += 1
    total = abce + ceab
    signed_bias = (abce - ceab) / total if total else 0.0
    bias = abs(signed_bias)
    return SieveStats(
        limit=limit,
        total=total,
        abce=abce,
        ceab=ceab,
        bias=bias,
        signed_bias=signed_bias,
    )
