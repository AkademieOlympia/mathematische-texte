from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math
import random

SIGMAS = ("E", "A", "B", "C")

# V4-Multiplikation fuer den EABC-Sektor.
MUL_TABLE = {
    ("E", "E"): "E", ("E", "A"): "A", ("E", "B"): "B", ("E", "C"): "C",
    ("A", "E"): "A", ("A", "A"): "E", ("A", "B"): "C", ("A", "C"): "B",
    ("B", "E"): "B", ("B", "A"): "C", ("B", "B"): "E", ("B", "C"): "A",
    ("C", "E"): "C", ("C", "A"): "B", ("C", "B"): "A", ("C", "C"): "E",
}

# Quaternionisch inspirierte Achsenprojektion:
# A -> i, B -> j, C -> k, E -> 0.
AXES = {
    "E": (0.0, 0.0, 0.0),
    "A": (1.0, 0.0, 0.0),
    "B": (0.0, 1.0, 0.0),
    "C": (0.0, 0.0, 1.0),
}

PAIR_INDICES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
NONTRIVIAL = ("A", "B", "C")


State = tuple[str, str, str, str]


@dataclass(frozen=True)
class ShellRecord:
    shell: int
    state: State
    energy: float
    axis_norm: float
    signature: str


def mul(x: str, y: str) -> str:
    return MUL_TABLE[(x, y)]


def total_signature(state: State) -> str:
    sig = "E"
    for x in state:
        sig = mul(sig, x)
    return sig


def shift(state: State) -> State:
    return (state[1], state[2], state[3], state[0])


def shift_inv(state: State) -> State:
    return (state[3], state[0], state[1], state[2])


def pair_flip(state: State, i: int, j: int, g: str) -> State:
    data = list(state)
    data[i] = mul(g, data[i])
    data[j] = mul(g, data[j])
    return tuple(data)  # type: ignore[return-value]


def axis_sum(state: State) -> tuple[float, float, float]:
    sx = sy = sz = 0.0
    for sigma in state:
        ax, ay, az = AXES[sigma]
        sx += ax
        sy += ay
        sz += az
    return sx, sy, sz


def axis_norm(state: State) -> float:
    sx, sy, sz = axis_sum(state)
    return math.sqrt(sx * sx + sy * sy + sz * sz)


def count_nontrivial(state: State) -> int:
    return sum(1 for x in state if x != "E")


def hamming(a: State, b: State) -> int:
    return sum(x != y for x, y in zip(a, b))


def eichler_variant_score(state: State) -> float:
    """
    Niedrige Werte markieren "balancierte" Vierlinge:
    - ABCE wird bevorzugt
    - AABB-artige Muster bleiben ebenfalls guenstig
    """
    cnt = Counter(state)
    spread = sum((cnt[sigma] - 1) ** 2 for sigma in SIGMAS)
    abce_bonus = 2.0 if all(cnt[sigma] == 1 for sigma in SIGMAS) else 0.0
    pair_balance = 0.6 if tuple(sorted(cnt.values(), reverse=True)) == (2, 2) else 0.0
    return spread - abce_bonus - pair_balance


def local_shell_energy(state: State, shell: int, previous: State | None = None) -> float:
    """
    Uebertraegt die Eichler-Idee auf einen einfachen Schalenbrenner:
    balancierte EABC-Zustaende werden in spaeteren Schalen leicht bevorzugt.
    """
    shell_target = shell % 4
    occupancy_penalty = abs(count_nontrivial(state) - shell_target)
    anisotropy_penalty = 0.35 * axis_norm(state)
    eichler_term = 0.9 * eichler_variant_score(state)
    transport_penalty = 0.25 * hamming(previous, state) if previous is not None else 0.0
    return occupancy_penalty + anisotropy_penalty + eichler_term + transport_penalty


def local_shell_energy_complex(
    state: State,
    shell: int,
    complexity: int = 4,
    previous: State | None = None,
) -> float:
    if complexity <= 0:
        raise ValueError("complexity muss positiv sein.")

    # Das Schalen-Target skaliert mit der gewaehlten Komplexitaet.
    shell_target = shell % complexity
    occupancy_penalty = abs(count_nontrivial(state) - shell_target)

    # Der Eichler-Schutz bleibt konstant, damit Vergleiche moeglich bleiben.
    eichler_term = 0.9 * eichler_variant_score(state)

    anisotropy_penalty = 0.35 * axis_norm(state)
    transport_penalty = 0.25 * hamming(previous, state) if previous is not None else 0.0

    return occupancy_penalty + anisotropy_penalty + eichler_term + transport_penalty


def neighbor_states(state: State) -> list[State]:
    neighbors = {state, shift(state), shift_inv(state)}
    for i, j in PAIR_INDICES:
        for g in NONTRIVIAL:
            candidate = pair_flip(state, i, j, g)
            if total_signature(candidate) == "E":
                neighbors.add(candidate)
    return sorted(neighbors)


def soft_choice(candidates: list[tuple[float, State]], temperature: float) -> State:
    if temperature <= 1e-9:
        return min(candidates, key=lambda item: item[0])[1]

    weights = [math.exp(-(energy - candidates[0][0]) / temperature) for energy, _ in candidates]
    return random.choices([state for _, state in candidates], weights=weights, k=1)[0]


def schalenbrenner_run(
    start: State = ("E", "A", "B", "C"),
    shells: int = 12,
    temperature: float = 0.35,
    seed: int = 42,
) -> list[ShellRecord]:
    random.seed(seed)

    if total_signature(start) != "E":
        raise ValueError("Startzustand muss im neutralen EABC-Sektor liegen.")

    history: list[ShellRecord] = []
    state = start
    previous: State | None = None

    for shell in range(shells):
        candidates = []
        for candidate in neighbor_states(state):
            energy = local_shell_energy(candidate, shell=shell, previous=previous)
            candidates.append((energy, candidate))

        candidates.sort(key=lambda item: item[0])
        state = soft_choice(candidates, temperature=temperature)
        energy = local_shell_energy(state, shell=shell, previous=previous)
        history.append(
            ShellRecord(
                shell=shell,
                state=state,
                energy=energy,
                axis_norm=axis_norm(state),
                signature=total_signature(state),
            )
        )
        previous = state

    return history


def schalenbrenner_run_ext(
    start: State = ("E", "A", "B", "C"),
    shells: int = 12,
    temperature: float = 0.35,
    seed: int = 42,
    complexity: int = 4,
) -> list[ShellRecord]:
    random.seed(seed)

    if total_signature(start) != "E":
        raise ValueError("Startzustand muss im neutralen EABC-Sektor liegen.")

    history: list[ShellRecord] = []
    state = start
    previous: State | None = None

    for shell in range(shells):
        candidates = []
        for candidate in neighbor_states(state):
            energy = local_shell_energy_complex(
                candidate,
                shell=shell,
                complexity=complexity,
                previous=previous,
            )
            candidates.append((energy, candidate))

        candidates.sort(key=lambda item: item[0])
        state = soft_choice(candidates, temperature=temperature)
        energy = local_shell_energy_complex(
            state,
            shell=shell,
            complexity=complexity,
            previous=previous,
        )
        history.append(
            ShellRecord(
                shell=shell,
                state=state,
                energy=energy,
                axis_norm=axis_norm(state),
                signature=total_signature(state),
            )
        )
        previous = state

    return history


def print_run(history: list[ShellRecord]) -> None:
    print("Eichler-/Quaternionen-Prototyp fuer den Schalenbrenner")
    print("=" * 72)
    for record in history:
        print(
            f"Schale {record.shell:2d} | Zustand {record.state} | "
            f"E={record.energy:6.3f} | |axis|={record.axis_norm:5.3f} | "
            f"sig={record.signature}"
        )


if __name__ == "__main__":
    run = schalenbrenner_run()
    print_run(run)