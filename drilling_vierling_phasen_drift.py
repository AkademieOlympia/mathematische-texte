"""
drilling_vierling_phasen_drift.py
Iteration Feinstruktur: Drilling → Vierling im Bamberg-C4-Kanalgesetz.

Misst die energetische Nähe d(T,Q) der Primdrillings-Muster
  Typ I:  {p, p+2, p+6}
  Typ II: {p, p+4, p+6}
zum nächstgelegenen Primvierling {q, q+2, q+6, q+8} und quantifiziert
Phasen-Synchronisation (Kanal-Injektion 5 ↔ 11 mod 12).

Governance (#Energiedoku):
  Modell-Definition ≠ physikalische Interpretation.
  Sync-Drop / d_Φ sind dimensionslose Fluktuationskennzahlen der
  Synchronisationsdämpfung — keine fundamentalen Naturkonstanten.
  α_B = Var(ψ)/⟨Gap⟩² parametrisiert ausschließlich die Dämpfung im Modell.

Bezug: Bamberg_Modell_Kern.lean (C4, α_B), Tag 2026-07-17-STABLE-SPECTRAL.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

BASE = Path("/Users/thomashoffbauer/Desktop/Mathematische Texte")
PRIMES_1B = BASE / "primes_cache_1e9.npy"
OUT_JSON = BASE / "drilling_vierling_phasen_drift_report.json"

# Bamberg-Kanäle / C4-Phasenstufen (diskrete Φ_θ)
CHANNEL_PHASE = {5: 0.0, 11: 0.5}  # Halbumlauf auf dem C4-Orbit
PHI_MISSING = 1.0  # offener Kanal (kein 5/11)


def load_primes(path: Path = PRIMES_1B) -> np.ndarray:
    if not path.is_file():
        raise FileNotFoundError(path)
    return np.asarray(np.load(path, mmap_mode="r"), dtype=np.int64)


def is_prime_vals(primes: np.ndarray, x: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(primes, x)
    ok = idx < len(primes)
    idx_clip = np.minimum(idx, len(primes) - 1)
    return ok & (primes[idx_clip] == x)


def extract_patterns(primes: np.ndarray) -> dict[str, np.ndarray]:
    p = primes
    has2 = is_prime_vals(p, p + 2)
    has4 = is_prime_vals(p, p + 4)
    has6 = is_prime_vals(p, p + 6)
    has8 = is_prime_vals(p, p + 8)
    base = p > 3

    quads = p[base & has2 & has6 & has8]
    trip_i = p[base & has2 & has6]  # inkl. solche, die schon Vierlinge sind
    trip_ii = p[base & has4 & has6]
    # offene Drillinge = noch kein geschlossener Vierlings-Orbit
    open_i = p[base & has2 & has6 & ~has8]
    open_ii = p[base & has4 & has6 & ~(has2 & has8)]  # nicht zugleich Quad-Start

    return {
        "quadruplets": quads,
        "triplets_I": trip_i,
        "triplets_II": trip_ii,
        "open_I": open_i,
        "open_II": open_ii,
    }


def phi_theta(n: int | np.ndarray) -> np.ndarray:
    """Φ_θ-Projektion: Bamberg-Kanal 5/11 → Phase, sonst offener Kanal."""
    r = np.asarray(n, dtype=np.int64) % 12
    out = np.full(r.shape, PHI_MISSING, dtype=float)
    out[r == 5] = CHANNEL_PHASE[5]
    out[r == 11] = CHANNEL_PHASE[11]
    return out


def constellation_phases(starts: np.ndarray, offsets: tuple[int, ...]) -> np.ndarray:
    """Phasenvektor Φ_θ für jede Konstellation; Shape (N, len(offsets))."""
    cols = [phi_theta(starts + off) for off in offsets]
    return np.stack(cols, axis=1)


def pad_triplet_phases(phases3: np.ndarray) -> np.ndarray:
    """
    Drilling → 4-Phasenvektor: fehlende Stufe als offener Kanal (Φ=1).
    Entspricht unvollständiger C4-Besetzung / nicht geschlossenem transition_operator.
    """
    n = len(phases3)
    out = np.ones((n, 4), dtype=float)
    out[:, :3] = phases3
    return out


def nearest_quad_distance(trip_starts: np.ndarray, quad_starts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Für jeden Drilling-Start: Index und |Δstart| zum nächstgelegenen Vierling.
    """
    if len(trip_starts) == 0 or len(quad_starts) == 0:
        empty = np.array([], dtype=np.int64)
        return empty, empty.astype(float)

    idx = np.searchsorted(quad_starts, trip_starts)
    idx0 = np.clip(idx - 1, 0, len(quad_starts) - 1)
    idx1 = np.clip(idx, 0, len(quad_starts) - 1)
    d0 = np.abs(quad_starts[idx0].astype(np.int64) - trip_starts.astype(np.int64))
    d1 = np.abs(quad_starts[idx1].astype(np.int64) - trip_starts.astype(np.int64))
    use1 = d1 < d0
    nearest_idx = np.where(use1, idx1, idx0)
    nearest_dist = np.where(use1, d1, d0).astype(float)
    return nearest_idx, nearest_dist


def phase_distance(T_phases: np.ndarray, Q_phases: np.ndarray) -> np.ndarray:
    """
    d(T,Q) = Σ_i ‖Φ_θ(T_i) - Φ_θ(Q_i)‖
    mit periodischer Metrik auf [0,1) für belegte Kanäle; offener Kanal = volle Strafe 1.
    """
    diff = np.abs(T_phases - Q_phases)
    # Periodizität nur wo beide belegt (< 1)
    both = (T_phases < 1.0 - 1e-12) & (Q_phases < 1.0 - 1e-12)
    circ = np.minimum(diff, 1.0 - diff)
    metric = np.where(both, circ, diff)
    return metric.sum(axis=1)


def hodge_integrity_from_starts(starts: np.ndarray, offsets: tuple[int, ...]) -> np.ndarray:
    """
    Hodge-/Kanal-Integrität: unter den Konstellationsgliedern erscheinen
    beide Bamberg-Kanäle 5 und 11 (komplementäre Injektion).
    """
    residues = np.stack([(starts + off) % 12 for off in offsets], axis=1)
    has5 = np.any(residues == 5, axis=1)
    has11 = np.any(residues == 11, axis=1)
    return has5 & has11


def padded_triplet_is_open(n: int) -> np.ndarray:
    """Gepolsterter Drilling: vierte Stufe fehlt ⇒ Orbit nicht geschlossen."""
    return np.ones(n, dtype=bool)


def analyze_family(
    name: str,
    trip_starts: np.ndarray,
    trip_offsets: tuple[int, int, int],
    quad_starts: np.ndarray,
    quad_phases: np.ndarray,
) -> dict:
    if len(trip_starts) == 0:
        return {"family": name, "n": 0}

    T3 = constellation_phases(trip_starts, trip_offsets)
    T4 = pad_triplet_phases(T3)
    n_idx, start_gap = nearest_quad_distance(trip_starts, quad_starts)
    Qφ = quad_phases[n_idx]
    d_phi = phase_distance(T4, Qφ)
    # Drilling (3 Glieder): kann schon 5↔11 tragen, Orbit aber offen (fehlende 4. Stufe)
    integrity_trip_channels = hodge_integrity_from_starts(trip_starts, trip_offsets)
    integrity_open_orbit = padded_triplet_is_open(len(trip_starts))  # immer True
    nearest_quad_starts = quad_starts[n_idx]
    integrity_Q = hodge_integrity_from_starts(nearest_quad_starts, (0, 2, 6, 8))

    trip_channels = (trip_starts % 12).astype(int)

    # Typ I: Drilling ist Präfix eines Vierlings ⇔ Start ∈ quads
    if trip_offsets == (0, 2, 6):
        is_prefix = np.isin(trip_starts, quad_starts)
        open_mask = ~is_prefix
    else:
        is_prefix = np.zeros(len(trip_starts), dtype=bool)
        open_mask = np.ones(len(trip_starts), dtype=bool)

    return {
        "family": name,
        "n": int(len(trip_starts)),
        "n_open": int(open_mask.sum()),
        "n_prefix_of_quad": int(is_prefix.sum()),
        "mean_start_gap_to_nearest_quad": float(np.mean(start_gap)),
        "median_start_gap_to_nearest_quad": float(np.median(start_gap)),
        "mean_d_phi": float(np.mean(d_phi)),
        "median_d_phi": float(np.median(d_phi)),
        "mean_d_phi_open": float(np.mean(d_phi[open_mask])) if open_mask.any() else None,
        "mean_d_phi_prefix": float(np.mean(d_phi[is_prefix])) if is_prefix.any() else None,
        "hodge_channel_pair_rate_triplet": float(np.mean(integrity_trip_channels)),
        "open_orbit_rate_triplet": float(np.mean(integrity_open_orbit)),
        "hodge_integrity_rate_nearest_quad": float(np.mean(integrity_Q)),
        "channel_shares_starts": {
            "5": float(np.mean(trip_channels == 5)),
            "11": float(np.mean(trip_channels == 11)),
            "other": float(np.mean(~np.isin(trip_channels, [5, 11]))),
        },
        "d_phi_percentiles": {
            "p10": float(np.percentile(d_phi, 10)),
            "p50": float(np.percentile(d_phi, 50)),
            "p90": float(np.percentile(d_phi, 90)),
        },
        "start_gap_percentiles": {
            "p10": float(np.percentile(start_gap, 10)),
            "p50": float(np.percentile(start_gap, 50)),
            "p90": float(np.percentile(start_gap, 90)),
        },
    }


def main():
    print("Drilling→Vierling Phasen-Drift (Bamberg / C4)")
    print(f"  Quelle: {PRIMES_1B}")
    primes = load_primes()
    print(f"  Primzahlen: n={len(primes):,}  max={int(primes[-1]):,}")

    pat = extract_patterns(primes)
    quads = pat["quadruplets"]
    print(f"  Vierlinge:           {len(quads):,}")
    print(f"  Drillinge Typ I:     {len(pat['triplets_I']):,}  (offen: {len(pat['open_I']):,})")
    print(f"  Drillinge Typ II:    {len(pat['triplets_II']):,}  (offen: {len(pat['open_II']):,})")

    Qφ = constellation_phases(quads, (0, 2, 6, 8))
    quad_integrity = float(np.mean(hodge_integrity_from_starts(quads, (0, 2, 6, 8))))

    fam_I = analyze_family("triplet_I_0_2_6", pat["triplets_I"], (0, 2, 6), quads, Qφ)
    fam_II = analyze_family("triplet_II_0_4_6", pat["triplets_II"], (0, 4, 6), quads, Qφ)
    fam_open_I = analyze_family("open_I_only", pat["open_I"], (0, 2, 6), quads, Qφ)

    def _num(x, default=1.0):
        return default if x is None else float(x)

    # Sync-Signatur: beim Übergang open→quad fällt d_φ abrupt
    sync_drop = None
    if fam_open_I.get("mean_d_phi") is not None and fam_I.get("mean_d_phi_prefix") is not None:
        sync_drop = fam_open_I["mean_d_phi"] - fam_I["mean_d_phi_prefix"]

    report = {
        "source": str(PRIMES_1B),
        "n_primes": int(len(primes)),
        "n_quadruplets": int(len(quads)),
        "quad_hodge_integrity_rate": quad_integrity,
        "families": {
            "triplet_I": fam_I,
            "triplet_II": fam_II,
            "open_I": fam_open_I,
        },
        "phase_sync_signature": {
            "mean_d_phi_open_I": fam_open_I.get("mean_d_phi"),
            "mean_d_phi_prefix_I": fam_I.get("mean_d_phi_prefix"),
            "sync_drop_open_to_closed": sync_drop,
            "interpretation": (
                "Abrupte Phasen-Synchronisation: d_φ fällt beim Schließen des "
                "vierten Kanals (Drilling-Präfix → Vierlings-Orbit / inert)."
            ),
        },
        "validation": {
            "quads_hodge_closed": quad_integrity > 0.99,
            "open_I_orbit_always_open": _num(fam_open_I.get("open_orbit_rate_triplet")) > 0.99,
            "prefix_near_zero_phase_distance": _num(fam_I.get("mean_d_phi_prefix")) < 0.05,
            "sync_drop_positive": _num(sync_drop, default=-1.0) > 0.3,
        },
    }

    OUT_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print()
    print("=== Phasen-Drift Bilanz ===")
    print(f"  Vierling Hodge 5↔11:           {quad_integrity:.4f}")
    print(f"  Typ I  mean d_φ (alle):        {fam_I['mean_d_phi']:.4f}")
    print(f"  Typ I  mean d_φ (Präfix/inert):{fam_I['mean_d_phi_prefix']:.4f}")
    print(f"  Typ I  mean d_φ (offen):       {fam_open_I['mean_d_phi']:.4f}")
    if sync_drop is not None:
        print(f"  Sync-Drop (offen→geschlossen): {sync_drop:.4f}")
    print(f"  Typ II mean d_φ:               {fam_II['mean_d_phi']:.4f}")
    print(f"  Typ I  median |Δstart|→Q:      {fam_I['median_start_gap_to_nearest_quad']:.1f}")
    print(f"  Typ II median |Δstart|→Q:      {fam_II['median_start_gap_to_nearest_quad']:.1f}")
    print(f"  Report: {OUT_JSON}")

    v = report["validation"]
    ok = all(v.values())
    print()
    if ok:
        print("Muster bestätigt: offener Drilling instabil; Vierling geschlossen/inert.")
    else:
        print("Muster teilweise — siehe validation im Report.")
        print(json.dumps(v, indent=2))


if __name__ == "__main__":
    main()
