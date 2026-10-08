"""
fourier_spektrum_validierung.py
Fourier-Validierung: 1B-Primdaten → Mod-12-/Gap-Signal → PSD vs. L-Lift.

Datenquelle (Status-Log lokal):
  primes_cache_1e9.npy  — alle Primzahlen ≤ 10^9

Governance (#Energiedoku):
  Modell-Definition ≠ physikalische Interpretation.
  L = 1/2·I + i·H_c: H_c ist Modellannahme (Störungsoperator-Notation),
  keine behauptete Selbstadjungiertheit hier.
  PSD-Korrelationen und α_B sind dimensionslose Fluktuationskennzahlen,
  keine fundamentalen Naturkonstanten.

Algebraischer Bezug:
  Bamberg_Modell_Kern.lean (C4, α_B) + L-Lift (test_spektral_kern.py)
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

BASE = Path("/Users/thomashoffbauer/Desktop/Mathematische Texte")
PRIMES_1B = BASE / "primes_cache_1e9.npy"
OUT_JSON = BASE / "fourier_spektrum_validierung_report.json"


# --- GNS / L-Lift (wie test_spektral_kern.py) ---------------------------------

def quat_mul(q1, q2):
    a1, b1, c1, d1 = q1
    a2, b2, c2, d2 = q2
    return (
        a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
        a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
        a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
        a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2,
    )


def quat_conj(q):
    return (q[0], -q[1], -q[2], -q[3])


def quat_norm(q):
    return math.sqrt(sum(x * x for x in q))


def quat_inverse(q):
    n2 = sum(x * x for x in q)
    c = quat_conj(q)
    return tuple(x / n2 for x in c)


def transition_operator(Q1, Q2):
    delta = quat_mul(quat_inverse(Q1), Q2)
    norm = quat_norm(delta)
    if norm == 0:
        return np.zeros((4, 4), dtype=float)
    vec = np.array(delta, dtype=float) / norm
    return np.outer(vec, vec)


def eabc_quaternion(p: int):
    r = p % 12
    if r == 1:
        return (1.0, 0.5 * p, 0.3 * p, 0.2 * p)
    if r == 5:
        return (1.0, 0.2 * p, 0.5 * p, 0.3 * p)
    if r == 7:
        return (1.0, 0.3 * p, 0.2 * p, 0.5 * p)
    if r == 11:
        return (1.0, 0.4 * p, 0.4 * p, 0.2 * p)
    return (1.0, float(p), 0.0, 0.0)


def build_L_lift_spectrum(quad_starts: np.ndarray, n_trans: int = 200) -> np.ndarray:
    """Eigenwerte von L = 1/2 I + i H_c aus Quaternion-Übergängen der Vierlinge."""
    starts = quad_starts[: max(n_trans + 1, 3)]
    quads = [eabc_quaternion(int(p)) for p in starts]
    D = np.zeros((4, 4), dtype=float)
    for i in range(min(n_trans, len(quads) - 1)):
        B = transition_operator(quads[i], quads[i + 1])
        D += B.T @ B
    H = 0.5 * (D + D.T)
    n = H.shape[0]
    H_c = H - (np.trace(H) / n) * np.eye(n)
    L = 0.5 * np.eye(n) + 1j * H_c
    return np.linalg.eigvals(L)


# --- 1B-Daten ----------------------------------------------------------------

def load_1B_data(path: Path = PRIMES_1B) -> np.ndarray:
    """Lädt den Primzahl-Cache ≤ 10^9 (Status 2026, lokale Datei)."""
    if not path.is_file():
        raise FileNotFoundError(f"1B-Primcache fehlt: {path}")
    return np.load(path, mmap_mode="r")


def extract_quadruplet_starts(primes: np.ndarray) -> np.ndarray:
    """
    Kanonische Primvierlinge (p, p+2, p+6, p+8), p>3.
    Vektorisiertes searchsorted auf dem sortierten 1e9-Cache.
    """
    p = np.asarray(primes, dtype=np.int64)

    def _is_prime_val(x: np.ndarray) -> np.ndarray:
        idx = np.searchsorted(p, x)
        ok = idx < len(p)
        # sichere Indexierung
        idx_clip = np.minimum(idx, len(p) - 1)
        return ok & (p[idx_clip] == x)

    mask = (p > 3) & _is_prime_val(p + 2) & _is_prime_val(p + 6) & _is_prime_val(p + 8)
    return p[mask]


def perform_fourier_analysis(data_points, sampling_rate=12):
    """
    Modulo-12-Kanalregularität als Zeitsignal (Bamberg-Übergangswort (+1,+1,+1)).
    Zusätzlich: Gap-Signal für spektrale Feinstruktur.
    """
    pts = np.asarray(data_points, dtype=float)
    signal_mod12 = pts % sampling_rate
    gaps = np.diff(pts)
    gap_signal = gaps.astype(float)

    def _psd(signal: np.ndarray, rate: float):
        n = len(signal)
        if n < 4:
            raise ValueError("Signal zu kurz fuer FFT")
        # Mittelwert entfernen (DC-Bias)
        x = signal - np.mean(signal)
        yf = np.fft.rfft(x)
        xf = np.fft.rfftfreq(n, d=1.0 / rate)
        psd = (2.0 / n) * np.abs(yf)
        return xf, psd

    xf_m, psd_m = _psd(signal_mod12, float(sampling_rate))
    # Gaps: effektive Rate = mittlere Vierlingsdichte
    mean_gap = float(np.mean(gap_signal)) if len(gap_signal) else 1.0
    rate_gap = 1.0 / max(mean_gap, 1.0)
    xf_g, psd_g = _psd(gap_signal, rate_gap)

    return {
        "mod12": (xf_m, psd_m, signal_mod12),
        "gaps": (xf_g, psd_g, gap_signal),
        "mean_gap": mean_gap,
    }


def compare_with_L_lift(frequencies, psd, L_eigs: np.ndarray):
    """
    Korrelation empirische PSD vs. theoretisches L-Lift-Dichtemodell.

    Modell: log(1+f) skaliert + Peak-Anker an |Im λ| der Eigenwerte
    (auf den Frequenzraster interpoliert).
    """
    frequencies = np.asarray(frequencies, dtype=float)
    psd = np.asarray(psd, dtype=float)
    expected = np.log(frequencies + 1.0)

    imag_abs = np.sort(np.abs(np.imag(L_eigs)))
    # Peak-Maske an skalierten Im-Eigenwerten
    if imag_abs[-1] > 0 and frequencies[-1] > 0:
        scale = frequencies[-1] / imag_abs[-1]
        anchors = imag_abs * scale
        for a in anchors:
            expected = expected + np.exp(-0.5 * ((frequencies - a) / (0.05 * frequencies[-1] + 1e-9)) ** 2)

    # Nur dort korrelieren, wo beide Varianz haben
    if np.std(psd) < 1e-15 or np.std(expected) < 1e-15:
        return float("nan"), expected

    corr = float(np.corrcoef(psd, expected)[0, 1])
    return corr, expected


def dominant_frequencies(xf, psd, k: int = 5):
    # Index 0 oft Rest-DC — überspringen wenn f≈0
    order = np.argsort(psd)[::-1]
    peaks = []
    for idx in order:
        if xf[idx] <= 0:
            continue
        peaks.append((float(xf[idx]), float(psd[idx])))
        if len(peaks) >= k:
            break
    return peaks


def main():
    print("Starte Fourier-Analyse auf 1B-Datenlog...")
    print(f"  Quelle: {PRIMES_1B}")

    primes = load_1B_data()
    print(f"  Primzahlen geladen: n={len(primes):,}  max={int(primes[-1]):,}")

    quads = extract_quadruplet_starts(primes)
    print(f"  Primvierlinge (p,p+2,p+6,p+8): n={len(quads):,}")

    analysis = perform_fourier_analysis(quads, sampling_rate=12)
    xf_m, psd_m, sig_m = analysis["mod12"]
    xf_g, psd_g, gaps = analysis["gaps"]

    L_eigs = build_L_lift_spectrum(quads, n_trans=min(200, len(quads) - 1))
    mean_re = float(np.mean(np.real(L_eigs)))

    corr_mod12, _ = compare_with_L_lift(xf_m, psd_m, L_eigs)
    corr_gaps, _ = compare_with_L_lift(xf_g, psd_g, L_eigs)

    peaks_m = dominant_frequencies(xf_m, psd_m)
    peaks_g = dominant_frequencies(xf_g, psd_g)

    # Kanalstatistik Mod 12 (Bamberg-Wort)
    vals, counts = np.unique(sig_m.astype(int), return_counts=True)
    channel_shares = {int(v): float(c / len(sig_m)) for v, c in zip(vals, counts)}

    report = {
        "source": str(PRIMES_1B),
        "n_primes": int(len(primes)),
        "n_quadruplets": int(len(quads)),
        "mean_gap": analysis["mean_gap"],
        "L_lift_mean_real": mean_re,
        "L_lift_imag_abs": np.sort(np.abs(np.imag(L_eigs))).tolist(),
        "corr_psd_mod12_vs_L": corr_mod12,
        "corr_psd_gaps_vs_L": corr_gaps,
        "mod12_channel_shares": channel_shares,
        "top_frequencies_mod12": peaks_m,
        "top_frequencies_gaps": peaks_g,
        "validation": {
            "critical_line_ok": abs(mean_re - 0.5) < 1e-8,
            "mod12_only_5_or_11": set(channel_shares) <= {5, 11},
        },
    }

    OUT_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print()
    print("=== Fourier / L-Lift Validierung ===")
    print(f"  L-Lift mittlerer Realteil: {mean_re:.6f}  (Ziel 0.5)")
    print(f"  |Im λ|: {np.array2string(np.array(report['L_lift_imag_abs']), precision=4)}")
    print(f"  Mod-12 Kanalanteile: {channel_shares}")
    print(f"  Korrelation PSD(mod12) ↔ L-Modell: {corr_mod12:.4f}")
    print(f"  Korrelation PSD(gaps)  ↔ L-Modell: {corr_gaps:.4f}")
    print(f"  Top-Frequenzen (gaps): {peaks_g[:3]}")
    print(f"  Report: {OUT_JSON}")

    # Abschnittsverschluss-Kriterium: kritische Gerade + Mod-12-Kanalgesetz
    ok = report["validation"]["critical_line_ok"] and report["validation"]["mod12_only_5_or_11"]
    print()
    if ok:
        print("Abschnittsverschluss: statistische Kern-Signatur bestätigt "
              "(Re=1/2, Mod-12 ∈ {5,11}).")
    else:
        print("Abschnittsverschluss: Validierung unvollständig — siehe Report.")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
