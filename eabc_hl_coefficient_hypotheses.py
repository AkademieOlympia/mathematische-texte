#!/usr/bin/env python3
"""
Hardy–Littlewood-Koeffizienten für Primzahlvierlinge (p, p+2, p+6, p+8)
auf den sechs mod-420-Halbkanälen R = {11, 101, 191, 221, 311, 401}.

Berechnet:
  - Singulärreihe S(Q) und lokale Faktoren mod 2, 3, 4, 5, 7, 420
  - Kanalspezifische HL-Gewichte (falls unterschiedlich)
  - Vergleich mit beobachteten Anteilen p_i = N_i/N
  - χ²-Test und KL-Divergenz HL vs. Beobachtung
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

# ── Beobachtete Daten (N_max = 10^8) ─────────────────────────────────────────
R_VALS = [11, 101, 191, 221, 311, 401]
COUNTS = np.array([765, 831, 786, 809, 809, 767], dtype=float)
N_TOTAL = float(COUNTS.sum())
OBSERVED = COUNTS / N_TOTAL
DELTA = OBSERVED - 1 / 6

OFFSETS = (0, 2, 6, 8)
SMALL_PRIMES = (2, 3, 5, 7)
MODULI = (2, 3, 4, 5, 7, 420)

OUT_DIR = Path(__file__).resolve().parent
REPORT_PY = OUT_DIR / "eabc_hl_coefficient_hypotheses_report.txt"
REPORT_JSON = OUT_DIR / "eabc_hl_coefficient_hypotheses.json"


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0:
        return False
    r = int(math.isqrt(n))
    for d in range(3, r + 1, 2):
        if n % d == 0:
            return False
    return True


def primes_up_to(limit: int) -> list[int]:
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    r = int(math.isqrt(limit))
    for p in range(2, r + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * ((limit - start) // p + 1)
    return [i for i in range(2, limit + 1) if sieve[i]]


def admissible_mod(q: int, residue: int) -> bool:
    """Alle vier Positionen p, p+2, p+6, p+8 sind ≠ 0 mod q."""
    r = int(residue) % q
    return all((r + h) % q != 0 for h in OFFSETS)


def omega(q: int) -> int:
    """Anzahl zulässiger Startrestklassen mod q."""
    if is_prime(q) and q > max(OFFSETS):
        # Für Primzahl p > max(OFFSETS): verboten sind n ≡ -h (mod p), h ∈ OFFSETS
        forb = {(q - h) % q for h in OFFSETS}
        return q - len(forb)
    return sum(1 for r in range(q) if admissible_mod(q, r))


def forbidden_residues_mod(q: int) -> set[int]:
    return {r for r in range(q) if not admissible_mod(q, r)}


def singular_factor_bateman_horn(p: int) -> float:
    """
    Lokaler HL-Faktor für Konstellation (0,2,6,8) bei Primzahl p:
      w_p = 1 - ν_p / (p-1)^4,
    wobei ν_p = Anzahl verbotener Restklassen mod p.
    """
    nu = p - omega(p)
    return 1.0 - nu / (p - 1) ** 4


def singular_series_global(prime_limit: int = 200_000) -> dict:
    """Globale Singulärreihe S(Q) = ∏_{p>2} w_p (p=2 separat)."""
    prod = 1.0
    factors: dict[int, float] = {}
    for p in primes_up_to(prime_limit):
        if p == 2:
            continue
        w = singular_factor_bateman_horn(p)
        prod *= w
        factors[p] = w
    return {
        "S_partial": prod,
        "prime_limit": prime_limit,
        "small_prime_factors": {p: factors[p] for p in SMALL_PRIMES if p in factors},
    }


def local_weight_mod_q(q: int, residue: int) -> float:
    """
    Relatives HL-Gewicht für Start p ≡ residue mod q:
    1/ω(q) wenn zulässig, sonst 0.
    (Gleichverteilung auf zulässigen Klassen mod q.)
    """
    if not admissible_mod(q, residue):
        return 0.0
    w = omega(q)
    return 1.0 / w if w > 0 else 0.0


def hl_weight_channel(r: int, moduli: Iterable[int] = MODULI) -> float:
    """Produkt der lokalen Gleichverteilungsgewichte über die Moduln."""
    w = 1.0
    for q in moduli:
        w *= local_weight_mod_q(q, r)
    return w


def hl_weight_channel_with_small_singular(r: int) -> float:
    """
    Gewicht inkl. expliziter Singulärreihen-Faktoren für p ∈ {2,3,5,7}:
    w(r) = (∏_{q|420} 1/ω(q)) · (∏_{p∈{3,5,7}} (1 - ν_p/(p-1)^4))
    (p=2: nur Zulässigkeit, Faktor identisch für alle ungeraden r)
    """
    base = hl_weight_channel(r, moduli=(2, 3, 4, 5, 7))
    for p in (3, 5, 7):
        if admissible_mod(p, r):
            base *= singular_factor_bateman_horn(p)
        else:
            return 0.0
    return base


def channel_residue_table(r: int) -> dict:
    row: dict = {"r": r}
    for h in OFFSETS:
        for q in MODULI:
            key = f"p+{h} mod {q}" if h else f"p mod {q}"
            row[key] = (r + h) % q
    for q in MODULI:
        row[f"adm_{q}"] = admissible_mod(q, r)
        row[f"omega_{q}"] = omega(q)
        row[f"loc_w_{q}"] = local_weight_mod_q(q, r)
    row["r_mod4"] = r % 4
    row["r_mod5"] = r % 5
    row["r_mod7"] = r % 7
    return row


@dataclass
class ChannelHL:
    r: int
    count: float
    observed: float
    delta: float
    hl_uniform: float
    hl_weight_raw: float
    hl_weight_norm: float
    admissible_all: bool
    r_mod7: int
    r_mod4: int
    singular_small: float


def normalize_weights(weights: np.ndarray) -> np.ndarray:
    s = weights.sum()
    return weights / s if s > 0 else weights


def chi2_stat(observed: np.ndarray, expected: np.ndarray) -> float:
    mask = expected > 0
    return float(np.sum((observed[mask] - expected[mask]) ** 2 / expected[mask]))


def kl_divergence(p: np.ndarray, q: np.ndarray) -> float:
    mask = (p > 0) & (q > 0)
    return float(np.sum(p[mask] * np.log(p[mask] / q[mask])))


def enumerate_admissible_mod420() -> list[int]:
    adm: list[int] = []
    for r in range(420):
        if all(admissible_mod(q, r) for q in (2, 3, 4, 5, 7)):
            adm.append(r)
    return sorted(adm)


def main() -> None:
    lines: list[str] = []
    p = lines.append

    p("=" * 72)
    p("Hardy–Littlewood-Koeffizienten: Primzahlvierlinge mod 420")
    p("=" * 72)
    p(f"R = {R_VALS}")
    p(f"Counts = {COUNTS.tolist()}, N = {int(N_TOTAL)}")
    p(f"Offsets = {OFFSETS}")
    p("")

    # ── Globale Singulärreihe ────────────────────────────────────────────────
    sg = singular_series_global()
    p("1. Globale Singulärreihe S(Q) für (0,2,6,8)")
    p(f"   S(Q) ≈ {sg['S_partial']:.12f}  (partielle Produkt bis p≤{sg['prime_limit']})")
    p("   Kleine Primfaktoren w_p = 1 - ν_p/(p-1)^4:")
    for pr, wf in sg["small_prime_factors"].items():
        nu = pr - omega(pr)
        p(f"     p={pr}: ν_p={nu}, ω_p={omega(pr)}, w_p={wf:.10f}")
    p("")

    # ── Zulässige Klassen mod 420 ────────────────────────────────────────────
    adm420 = enumerate_admissible_mod420()
    p("2. Zulässige Restklassen mod 420")
    p(f"   ω_420 = {len(adm420)} Klassen: {adm420}")
    p(f"   R ⊆ admissible: {set(R_VALS) <= set(adm420)}")
    p(f"   R = admissible Klassen: {sorted(R_VALS) == adm420}")
    p("")

    # ── Lokale Struktur pro Kanal ────────────────────────────────────────────
    p("3. Lokale Kongruenzbedingungen pro Kanal")
    p(f"{'r':>5} {'p%4':>4} {'p%5':>4} {'p%7':>4} "
      f"{'ω₂':>4} {'ω₃':>4} {'ω₄':>4} {'ω₅':>4} {'ω₇':>4} "
      f"{'w_loc':>10} {'S_small':>10}")
    channels: list[ChannelHL] = []
    raw_weights = []
    small_sing = []
    for r, c in zip(R_VALS, COUNTS):
        tbl = channel_residue_table(r)
        w_raw = hl_weight_channel(r)
        w_small = hl_weight_channel_with_small_singular(r)
        raw_weights.append(w_raw)
        small_sing.append(w_small)
        ch = ChannelHL(
            r=r,
            count=c,
            observed=c / N_TOTAL,
            delta=c / N_TOTAL - 1 / 6,
            hl_uniform=1 / 6,
            hl_weight_raw=w_raw,
            hl_weight_norm=0.0,
            admissible_all=all(tbl[f"adm_{q}"] for q in MODULI),
            r_mod7=r % 7,
            r_mod4=r % 4,
            singular_small=w_small,
        )
        channels.append(ch)
        p(
            f"{r:5d} {r%4:4d} {r%5:4d} {r%7:4d} "
            f"{omega(2):4d} {omega(3):4d} {omega(4):4d} {omega(5):4d} {omega(7):4d} "
            f"{w_raw:10.6f} {w_small:10.6f}"
        )
    p("")

    raw_w = np.array(raw_weights)
    small_w = np.array(small_sing)
    norm_raw = normalize_weights(raw_w)
    norm_small = normalize_weights(small_w)

    p("4. HL-Gewichte (normalisiert auf R)")
    p(f"   Roh-Produkt 1/ω(q): alle gleich = {np.allclose(raw_w, raw_w[0])}")
    p(f"   Mit Singulärfaktoren p∈{{3,5,7}}: alle gleich = {np.allclose(small_w, small_w[0])}")
    p(f"   Theoretische Anteile (uniform): {[1/6]*6}")
    p(f"   Theoretische Anteile (raw HL):  {norm_raw.round(6).tolist()}")
    p(f"   Beobachtete Anteile:            {OBSERVED.round(6).tolist()}")
    p("")

    # ── Statistischer Vergleich ──────────────────────────────────────────────
    exp_uniform = np.full(6, 1 / 6)
    chi2_uni = chi2_stat(COUNTS, exp_uniform * N_TOTAL)
    chi2_hl = chi2_stat(COUNTS, norm_raw * N_TOTAL)
    kl_uni = kl_divergence(OBSERVED, exp_uniform)
    kl_hl = kl_divergence(OBSERVED, norm_raw)

    p("5. Statistischer Vergleich (N = 4767)")
    p(f"   χ² (HL uniform 1/6):     {chi2_uni:.4f}  (df=5, p≈{1-chi2_cdf(chi2_uni,5):.4f})")
    p(f"   χ² (kanalspez. HL):      {chi2_hl:.4f}  (identisch wenn gleiche Gewichte)")
    p(f"   KL(Obs || Uniform):      {kl_uni:.8f}")
    p(f"   KL(Obs || HL-Vorhersage):{kl_hl:.8f}")
    p(f"   max |δ_i|:               {np.max(np.abs(DELTA)):.6f}")
    p(f"   Kanal 101 Überschuss:    {DELTA[1]*100:.3f}%")
    p(f"   Kanal 11 Unterschuss:    {DELTA[0]*100:.3f}%")
    p("")

    # ── Mod-7 / χ₃ Struktur ──────────────────────────────────────────────────
    p("6. Mod-7-Struktur und χ₃-Verbindung")
    mod7_groups: dict[int, list[int]] = {}
    for ch in channels:
        mod7_groups.setdefault(ch.r_mod7, []).append(ch.r)
    for g, rs in sorted(mod7_groups.items()):
        sub_counts = [COUNTS[R_VALS.index(r)] for r in rs]
        sub_frac = sum(sub_counts) / N_TOTAL
        p(f"   p≡{g} mod 7: Kanäle {rs}, Anteil {sub_frac:.4f} (Erwartung {len(rs)/6:.4f})")
    p("   χ₃-Zyklus auf R: r mod 7 = (4,3,2,4,3,2)")
    p("")

    # ── Verbotene kleine Primteiler ────────────────────────────────────────────
    p("7. Ausgeschlossene/erzwungene Primteiler (mod kleine q)")
    for q in (2, 3, 4, 5, 7):
        forb = forbidden_residues_mod(q)
        p(f"   mod {q}: verbotene p-Klassen = {sorted(forb)}")
        p(f"          erzwungen für R: p ≡ {[r % q for r in R_VALS]}")
    p("")

    # ── Hypothesen-Zusammenfassung ───────────────────────────────────────────
    p("8. Hypothesen-Bewertung (numerisch)")
    hypotheses = [
        ("H1", "Identische HL-Koeffizienten → Gleichverteilung",
         np.allclose(norm_raw, exp_uniform), "HL sagt 1/6 für jeden Kanal"),
        ("H2", "Mod-7-Struktur → unterschiedliche HL-Dichten",
         not np.allclose(norm_raw, exp_uniform), "ω_7=3, je 2 Kanäle/Klasse; HL-Gewichte gleich"),
        ("H3", "Mod-4/mod-5 → χ₄-Komponente via HL",
         not np.allclose(norm_raw, exp_uniform), "ω_4=2, ω_5=1; alle r∈R gleiche lokale Faktoren"),
        ("H4", "Korrelationen zwischen Vierlingen",
         None, "HL-Modell i.i.d.; Übergangsmatrix aus Stiefel.py testbar"),
        ("H5", "Finite-size bis 10^8",
         chi2_uni < 11.07, f"χ²={chi2_uni:.2f} konsistent mit Fluktuation (p≈{1-chi2_cdf(chi2_uni,5):.3f})"),
    ]
    for hid, title, verdict, note in hypotheses:
        vstr = "STÜTZT" if verdict is True else ("WIDERLEGT" if verdict is False else "OFFEN")
        p(f"   {hid}: {title}")
        p(f"        → {vstr}: {note}")
    p("")

    # ── Speichern ────────────────────────────────────────────────────────────
    report_txt = "\n".join(lines)
    REPORT_PY.write_text(report_txt, encoding="utf-8")

    payload = {
        "R": R_VALS,
        "counts": COUNTS.tolist(),
        "N": int(N_TOTAL),
        "offsets": list(OFFSETS),
        "singular_series": sg,
        "admissible_mod420": adm420,
        "channels": [asdict(ch) for ch in channels],
        "hl_weights_normalized": norm_raw.tolist(),
        "observed": OBSERVED.tolist(),
        "chi2_uniform": chi2_uni,
        "chi2_hl": chi2_hl,
        "kl_uniform": kl_uni,
        "kl_hl": kl_hl,
        "all_hl_weights_equal": bool(np.allclose(raw_w, raw_w[0])),
        "residue_tables": [channel_residue_table(r) for r in R_VALS],
    }
    REPORT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(report_txt)
    print(f"\nGespeichert: {REPORT_PY}")
    print(f"Gespeichert: {REPORT_JSON}")


def chi2_cdf(x: float, k: int) -> float:
    """χ²-CDF via reguläre untere unvollständige Gamma-Funktion (Reihe)."""
    if x <= 0:
        return 0.0
    a = k / 2.0
    z = x / 2.0
    term = 1.0 / a
    total = term
    n = 1
    while term > 1e-15:
        term *= z / (a + n)
        total += term
        n += 1
    return float(total * math.exp(-z + a * math.log(z) - math.lgamma(a)))


if __name__ == "__main__":
    main()
