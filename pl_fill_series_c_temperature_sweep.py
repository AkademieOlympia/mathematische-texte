
#!/usr/bin/env python3
"""
PL-FILL v1.1 — PREREGISTERED SERIES C: THERMAL T-SWEEP

Goal:
Test whether thermally activated uphill acceptance opens Pachner sectors
that are inaccessible in the low-T / descent-dominated regime.

Frozen controls:
    alpha = 1.0
    gamma = 1.0
    lambda_assoc = 0.0
    T in {0.05, 0.10, 0.20, 0.50}
    identical initial state
    identical deterministic Hurwitz edge rule
    identical move generator
    associator backend SEALED throughout blind evolution

Important:
This is a stochastic relaxation experiment. It is not claimed to sample an
exact Gibbs distribution unless proposal symmetry or a Hastings correction
has separately been established.
"""

import copy
import json
import math
import random
import hashlib
import numpy as np
from dataclasses import dataclass, asdict

TEMPERATURES = (0.05, 0.10, 0.20, 0.50)
ALPHA = 1.0
GAMMA = 1.0
N_STEPS = 200
N_REPLICATES = 32
BASE_SEED = 20260909
REFERENCE_BARRIER = 1.56


@dataclass(frozen=True)
class SeriesCPreregistration:
    temperatures: tuple = TEMPERATURES
    alpha: float = ALPHA
    gamma: float = GAMMA
    lambda_assoc: float = 0.0
    n_steps: int = N_STEPS
    n_replicates: int = N_REPLICATES
    base_seed: int = BASE_SEED
    reference_barrier: float = REFERENCE_BARRIER
    acceptance_rule: str = "downhill accept; uphill accept with exp(-DeltaS/T)"
    assoc_curv_rule: str = "SEALED during run; evaluate only after state freeze/hash"
    epistemic_status: str = "stochastic relaxation heuristic, not exact Gibbs sampling"


def barrier_probability(delta_s, T):
    if delta_s <= 0:
        return 1.0
    if T <= 0:
        return 0.0
    return math.exp(-delta_s / T)


def canonical_freeze_hash(engine):
    """
    Freeze visible state only. Engine must already have deterministic,
    side-effect-free edge assignments.
    """
    payload = {
        "coords": {
            str(k): [float(x) for x in np.asarray(v)]
            for k, v in sorted(engine.coords.items(), key=lambda kv: str(kv[0]))
        },
        "tets": {
            str(k): [str(v) for v in engine.tets[k].vertices]
            for k in sorted(engine.tets, key=str)
        },
        "edges": {}
    }

    for key in sorted(engine.edge_connections, key=lambda e: (str(e[0]), str(e[1]))):
        U = np.asarray(engine.edge_connections[key], dtype=complex)
        payload["edges"][f"{key[0]}|{key[1]}"] = [
            [[float(z.real), float(z.imag)] for z in row] for row in U
        ]

    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def run_one(engine_factory, T, replicate, prereg):
    seed = prereg.base_seed + int(round(T * 1000))*10000 + replicate
    random.seed(seed)
    np.random.seed(seed % (2**32 - 1))

    engine = engine_factory(seed=seed, temperature=T)

    if getattr(engine, "face_assoc", None) is not None:
        raise RuntimeError("Blindness violation before thermal sweep.")

    S0 = float(engine.evaluate_blind_action())

    accepted = 0
    uphill_proposed = 0
    uphill_accepted = 0
    n14_prop = n14_acc = 0
    n23_prop = n23_acc = 0
    max_n3 = len(engine.tets)
    trace = []

    for step in range(1, prereg.n_steps + 1):
        if getattr(engine, "face_assoc", None) is not None:
            raise RuntimeError("Blindness violation during thermal sweep.")

        res = engine.step_stochastic_move()
        move = str(res["move"])
        ds = float(res["delta_S"])
        acc = bool(res["accepted"])

        if ds > 0:
            uphill_proposed += 1
            if acc:
                uphill_accepted += 1

        if move.startswith("1_to_4"):
            n14_prop += 1
            n14_acc += int(acc)
        elif move.startswith("2_to_3"):
            n23_prop += 1
            n23_acc += int(acc)

        accepted += int(acc)
        max_n3 = max(max_n3, len(engine.tets))

        trace.append({
            "step": step,
            "move": move,
            "delta_S": ds,
            "accepted": acc,
            "N3_after": len(engine.tets),
        })

    Sf = float(engine.evaluate_blind_action())
    digest = canonical_freeze_hash(engine)

    return {
        "T": T,
        "replicate": replicate,
        "seed": seed,
        "S_initial": S0,
        "S_final": Sf,
        "delta_S_total": Sf - S0,
        "acceptance_rate": accepted / prereg.n_steps,
        "uphill_proposed": uphill_proposed,
        "uphill_accepted": uphill_accepted,
        "uphill_acceptance_rate":
            uphill_accepted / uphill_proposed if uphill_proposed else None,
        "one_to_four_proposed": n14_prop,
        "one_to_four_accepted": n14_acc,
        "two_to_three_proposed": n23_prop,
        "two_to_three_accepted": n23_acc,
        "N3_final": len(engine.tets),
        "N3_max": max_n3,
        "freeze_sha256": digest,
        "assoc_curv_status": "SEALED_NOT_EVALUATED",
        "trace": trace,
    }


def run_series_c(engine_factory, output="pl_fill_series_c_temperature_sweep.json"):
    prereg = SeriesCPreregistration()
    results = []

    for T in prereg.temperatures:
        for rep in range(prereg.n_replicates):
            results.append(run_one(engine_factory, T, rep, prereg))

    artifact = {
        "preregistration": asdict(prereg),
        "reference_acceptance_probabilities": {
            str(T): barrier_probability(prereg.reference_barrier, T)
            for T in prereg.temperatures
        },
        "results": results,
        "status": "BLIND_SERIES_C_COMPLETE_ASSOC_CURV_SEALED"
    }

    with open(output, "w", encoding="utf-8") as f:
        json.dump(artifact, f, indent=2, sort_keys=True)

    return artifact


if __name__ == "__main__":
    p = SeriesCPreregistration()
    print("=" * 72)
    print("PL-FILL v1.1 — PREREGISTERED SERIES C: THERMAL T-SWEEP")
    print("=" * 72)
    print("T values:", p.temperatures)
    print("alpha =", p.alpha, "gamma =", p.gamma, "lambda =", p.lambda_assoc)
    print("replicates per T =", p.n_replicates, "steps per replicate =", p.n_steps)
    print("-" * 72)
    print("Reference uphill acceptance for DeltaS = 1.56:")
    for T in p.temperatures:
        print(f"T={T:0.2f} -> P_acc={barrier_probability(1.56,T):.10f}")
    print("-" * 72)
    print("ASSOC-CURV remains SEALED until each final state is frozen and hashed.")
    print("Numerical sweep requires the final patched engine_factory.")
