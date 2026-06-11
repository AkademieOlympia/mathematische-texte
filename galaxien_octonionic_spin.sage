# -*- mode: sage -*-
# Ausführung:  sage galaxien_octonionic_spin.sage
# Erfordert SageMath (Oktaven: Octonions(QQ)).

from sage.all import *
import random


def octonionic_spin_simulation(n_galaxies=200_000, drift_scale=None, seed=None):
    """
    Oktave O über QQ; e1,e2,e3 ~ quaternionische Scheibe; e7: globaler Drall.

    Wendet p(cw) = 0,5 + | Drift | an. Für |e7|=1 (orthonormale Fano-Basis in Sage) ist
    |drift_scale * e7| = |drift_scale|.

    (local_v beinhaltet die pro-Galaxie quaternionalen Anteile – aktuell nur motivierend
    in der Rechnung, die Entscheidung nutzt den festen p aus dem globalen Drall.)
    """
    if drift_scale is None:
        drift_scale = QQ(1) / 200
    if seed is not None:
        set_random_seed(int(seed))
        random.seed(int(seed))
    if n_galaxies <= 0:
        return float("nan")

    O = Octonions(QQ)
    basis = list(O.basis())

    # e0 = 1; e1, e2, e3, …, e7
    e1, e2, e3 = basis[1], basis[2], basis[3]
    global_drift = basis[7] * drift_scale
    nd = abs(RR(drift_scale))
    p_cw = min(RR(0.999), max(RR(0.001), RR(0.5) + nd))

    count_cw = 0
    count_ccw = 0

    for _ in range(n_galaxies):
        local_v = (
            random.uniform(-1, 1) * e1
            + random.uniform(-1, 1) * e2
            + random.uniform(-1, 1) * e3
        )
        # Wechselwirkung (bisher nicht in p einbezogen) – hält e1,e2,e3 + Drift in der Kette
        _ = local_v + global_drift  # Oktave; nur auswerten, kein Abbruch
        if random.random() < p_cw:
            count_cw += 1
        else:
            count_ccw += 1

    return float((count_cw - count_ccw) / n_galaxies)


if __name__ == "__main__":
    n = 200_000
    # 0,005 = „4,4σ“-Kugel-Drall im Modellbild
    result_bias = octonionic_spin_simulation(
        n,
        drift_scale=QQ(5) / 1000,
        seed=1,
    )
    print(f"Simulation mit n={n}")
    print(f"Okt. induzierter Bias: {result_bias:.5f}")
