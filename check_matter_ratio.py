from pathlib import Path

import numpy as np


def _prime_sieve(limit: int) -> np.ndarray:
    """Return a boolean sieve where True marks prime numbers up to limit."""
    if limit < 2:
        return np.zeros(limit + 1, dtype=bool)

    is_prime = np.ones(limit + 1, dtype=bool)
    is_prime[:2] = False

    for p in range(2, int(np.sqrt(limit)) + 1):
        if is_prime[p]:
            is_prime[p * p : limit + 1 : p] = False

    return is_prime


def _prime_quartet_shells(limit: int) -> set[int]:
    is_prime = _prime_sieve(limit)
    shells: set[int] = set()

    for p in range(11, limit - 7):
        if is_prime[p] and is_prime[p + 2] and is_prime[p + 6] and is_prime[p + 8]:
            shells.update((p, p + 2, p + 6, p + 8))

    return shells


def check_matter_ratio_at_timestamp(filepath: str | Path = "zeros6.npy") -> float:
    print("--- #Energiedoku: Analyse des Materieverhaeltnisses ---")

    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Die Datei wurde nicht gefunden: {path}")

    zeros = np.load(path).astype(float).ravel()
    if zeros.size == 0:
        raise ValueError(f"Die Datei enthaelt keine Werte: {path}")

    finite_mask = np.isfinite(zeros)
    if not np.all(finite_mask):
        removed = int(zeros.size - np.count_nonzero(finite_mask))
        print(f"Warnung: {removed} nicht-endliche Werte wurden ignoriert.")
        zeros = zeros[finite_mask]

    if zeros.size == 0:
        raise ValueError("Nach dem Entfernen nicht-endlicher Werte bleiben keine Daten uebrig.")

    norms = np.rint(zeros).astype(np.int64)

    # Primzahl-Sieb fuer die Identifikation der baryonischen Periodizitaet.
    sieve_limit = max(20, int(np.max(norms)) + 10)
    quartet_shells = _prime_quartet_shells(sieve_limit)

    # Zaehlung der besetzten Zustaende im Gitter.
    omega_b_nodes = sum(1 for n in norms if int(n) in quartet_shells)
    omega_c_nodes = int(norms.size) - omega_b_nodes

    print(f"Zustaende im Selberg-Sieb (Dunkle Materie Omega_c): {omega_c_nodes}")
    print(f"Zustaende in periodischen Vierlingen (Baryonisch Omega_b): {omega_b_nodes}")

    if omega_b_nodes == 0:
        print("--> Numerisch ermitteltes Verhaeltnis: nicht definiert (Omega_b = 0)")
        return float("inf")

    ratio = omega_c_nodes / omega_b_nodes
    print(f"--> Numerisch ermitteltes Verhaeltnis: {ratio:.4f}")

    # Abweichungen vom Planck-Wert koennen als Grad unvollstaendiger
    # Kristallisation im Szenario C interpretiert werden.
    return ratio


if __name__ == "__main__":
    check_matter_ratio_at_timestamp()
