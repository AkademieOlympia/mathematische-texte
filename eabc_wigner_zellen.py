"""
eabc_wigner_zellen.py
=====================
Arithmetische Wigner-Zellen bei Primvierlingen (p, p+2, p+6, p+8).

Klassifikation nach r = p mod 210 ∈ {11, 101, 191} und Muster ABCE/CEAB:

  Die sechs erlaubten Startklassen mod 420 für Primvierlinge sind:
      11, 101, 191, 221, 311, 401

  Mod-12-Signaturen (EABC-Flavor):
      p mod 420 = 11  → (C,E,A,B) = CEAB    [11%12=11=C, 13%12=1=E, 17%12=5=A, 19%12=7=B]
      p mod 420 = 221 → (A,B,C,E) = ABCE    [221%12=5=A, 223%12=7=B, 227%12=11=C, 229%12=1=E]
      p mod 420 = 101 → (A,B,C,E) = ABCE    [101%12=5=A, 103%12=7=B, 107%12=11=C, 109%12=1=E]
      p mod 420 = 311 → (C,E,A,B) = CEAB    [311%12=11=C, 313%12=1=E, 317%12=5=A, 319%12=7=B]
      p mod 420 = 191 → (C,E,A,B) = CEAB    [191%12=11=C, 193%12=1=E, 197%12=5=A, 199%12=7=B]
      p mod 420 = 401 → (A,B,C,E) = ABCE    [401%12=5=A, 403%12=7=B, 407%12=11=C, 409%12=1=E]

  Drei Wigner-Zellen (r = p mod 210):
      Zelle r=11:  ABCE ↔ p mod 420 = 221,  CEAB ↔ p mod 420 = 11
      Zelle r=101: ABCE ↔ p mod 420 = 101,  CEAB ↔ p mod 420 = 311
      Zelle r=191: ABCE ↔ p mod 420 = 401,  CEAB ↔ p mod 420 = 191

Ausgaben: Tabellen N=10^7, N=10^8, Fenster 72e6 < p ≤ 77e6
          Plot: eabc_wigner_zellen_plot.png
"""

import math
import time
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ---------------------------------------------------------------------------
# Sieb des Eratosthenes (kompakt, ~0.9 GB für N=10^8 mit uint8)
# ---------------------------------------------------------------------------

def sieb_primzahlen(N: int) -> np.ndarray:
    """Gibt alle Primzahlen ≤ N als numpy-Array zurück (Sieb)."""
    is_prime = np.ones(N + 1, dtype=np.bool_)
    is_prime[0] = is_prime[1] = False
    is_prime[4::2] = False
    for i in range(3, int(math.isqrt(N)) + 1, 2):
        if is_prime[i]:
            is_prime[i * i :: 2 * i] = False
    return np.where(is_prime)[0].astype(np.int64)


# ---------------------------------------------------------------------------
# Primvierlinge finden (p, p+2, p+6, p+8 alle prim)
# ---------------------------------------------------------------------------

def find_quadruplets(N: int) -> np.ndarray:
    """Gibt Startprimzahlen p aller Primvierlinge (p,p+2,p+6,p+8) bis N zurück."""
    t0 = time.time()
    print(f"  Siebe bis N={N:,} …", flush=True)
    is_prime = np.ones(N + 1, dtype=np.bool_)
    is_prime[0] = is_prime[1] = False
    is_prime[4::2] = False
    for i in range(3, int(math.isqrt(N)) + 1, 2):
        if is_prime[i]:
            is_prime[i * i :: 2 * i] = False

    # Bedingung: is_prime[p] & is_prime[p+2] & is_prime[p+6] & is_prime[p+8]
    lim = N - 8
    arr = np.arange(lim + 1, dtype=np.int64)
    mask = (
        is_prime[:lim + 1]
        & is_prime[2:lim + 3]
        & is_prime[6:lim + 7]
        & is_prime[8:lim + 9]
    )
    quads = arr[mask]
    print(f"  Fertig in {time.time()-t0:.2f}s — {len(quads):,} Vierlinge gefunden.", flush=True)
    return quads


# ---------------------------------------------------------------------------
# Mod-420-Klassifikation und Wigner-Zellen-Auswertung
# ---------------------------------------------------------------------------

# Erlaubte Startklassen mod 420
ERLAUBT_MOD420 = {11, 101, 191, 221, 311, 401}

# Für jede mod-210-Zelle r: (ABCE_mod420, CEAB_mod420)
ZELLEN = {
    11:  (221, 11),
    101: (101, 311),
    191: (401, 191),
}


def klassifiziere(quads: np.ndarray, p_min: int = 0, p_max: int = None):
    """
    Klassifiziert Vierlinge im Bereich (p_min, p_max] nach Wigner-Zellen.

    Gibt zurück:
        counts: dict  r -> {'ABCE': int, 'CEAB': int}
        total_abce, total_ceab, gesamt (exkl. p=5)
    """
    if p_max is None:
        sub = quads[quads > p_min]
    else:
        sub = quads[(quads > p_min) & (quads <= p_max)]

    # p=5 separat heraushalten (5 mod 420 = 5, kein reguläres Muster)
    sub = sub[sub != 5]

    counts = {r: {'ABCE': 0, 'CEAB': 0} for r in (11, 101, 191)}
    mod420 = sub % 420

    for r, (abce_val, ceab_val) in ZELLEN.items():
        counts[r]['ABCE'] = int(np.sum(mod420 == abce_val))
        counts[r]['CEAB'] = int(np.sum(mod420 == ceab_val))

    return counts


def drucke_tabelle(label: str, counts: dict):
    """Gibt die Wigner-Zellen-Tabelle auf der Konsole aus."""
    print(f"\n{'='*62}")
    print(f"  {label}")
    print(f"{'='*62}")
    print(f"  {'Zelle r':>8}  {'ABCE':>6}  {'CEAB':>6}  {'D_r':>6}  {'Q_r':>6}  {'W_r':>8}")
    print(f"  {'-'*54}")

    total_abce = total_ceab = 0
    for r in (11, 101, 191):
        a = counts[r]['ABCE']
        b = counts[r]['CEAB']
        d = a - b
        q = a + b
        w = d / q if q > 0 else float('nan')
        total_abce += a
        total_ceab += b
        print(f"  r={r:3d}       {a:>6}  {b:>6}  {d:>+7}  {q:>6}  {w:>+8.4f}")

    D_11  = counts[11]['ABCE']  - counts[11]['CEAB']
    D_101 = counts[101]['ABCE'] - counts[101]['CEAB']
    D_191 = counts[191]['ABCE'] - counts[191]['CEAB']
    S = D_11 + D_101 - D_191
    Q_total = total_abce + total_ceab
    sqrt_Q = math.sqrt(Q_total) if Q_total > 0 else float('nan')

    print(f"  {'-'*54}")
    print(f"  Gesamt:         {total_abce:>6}  {total_ceab:>6}  {total_abce-total_ceab:>+7}  {Q_total:>6}")
    print(f"  S = D_11 + D_101 - D_191 = {D_11:+} + {D_101:+} - ({D_191:+}) = {S:+}")
    print(f"  √Q = {sqrt_Q:.3f},  S/√Q = {S/sqrt_Q:+.4f}" if sqrt_Q > 0 else "")
    print(f"{'='*62}")
    return counts, S, sqrt_Q


# ---------------------------------------------------------------------------
# Visualisierung
# ---------------------------------------------------------------------------

def erstelle_plot(alle_counts: list, labels: list, outpath: str):
    """Balkendiagramm der W_r-Werte mit Fehlerbalken (1/√Q_r)."""
    fig, ax = plt.subplots(figsize=(9, 5))

    r_vals = [11, 101, 191]
    x = np.arange(len(r_vals))
    bar_width = 0.25
    colors = ['#2196F3', '#4CAF50', '#FF5722']

    for idx, (counts, label, color) in enumerate(zip(alle_counts, labels, colors)):
        w_vals = []
        err_vals = []
        for r in r_vals:
            a = counts[r]['ABCE']
            b = counts[r]['CEAB']
            q = a + b
            w = (a - b) / q if q > 0 else 0.0
            err = 1.0 / math.sqrt(q) if q > 0 else 0.0
            w_vals.append(w)
            err_vals.append(err)

        offset = (idx - 1) * bar_width
        bars = ax.bar(x + offset, w_vals, bar_width,
                      label=label, color=color, alpha=0.85,
                      yerr=err_vals, capsize=5, error_kw={'linewidth': 1.5})

    ax.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_xticks(x)
    ax.set_xticklabels([f'r = {r}\n(mod 210)' for r in r_vals], fontsize=12)
    ax.set_xlabel('Wigner-Zelle r', fontsize=13)
    ax.set_ylabel('W_r = D_r / Q_r', fontsize=13)
    ax.set_title('Arithmetische Wigner-Zellen bei Primvierlingen\n'
                 r'$W_r = (N_{\mathrm{ABCE}} - N_{\mathrm{CEAB}}) / (N_{\mathrm{ABCE}} + N_{\mathrm{CEAB}})$',
                 fontsize=13)
    ax.legend(fontsize=11)
    ax.set_ylim(-0.12, 0.12)
    ax.grid(axis='y', alpha=0.4)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()
    print(f"\n  Plot gespeichert: {outpath}")


# ---------------------------------------------------------------------------
# Hauptprogramm
# ---------------------------------------------------------------------------

def main():
    import os
    workspace = os.path.dirname(os.path.abspath(__file__))

    N_VOLL = 10_000_000   # Für schnelle Verifikation (wie vom User angegeben)
    N_GROSS = 100_000_000  # Optional (falls Zeit vorhanden)

    # ---------- Sieb bis 10^8 (brauchen wir für Fenster und N=10^8) ----------
    print("\n" + "="*62)
    print("  EABC Wigner-Zellen – Primvierlinge")
    print("="*62)

    COMPUTE_LARGE = True  # Auf False setzen für reine Schnellverifikation

    N_sieb = N_GROSS if COMPUTE_LARGE else N_VOLL
    quads_all = find_quadruplets(N_sieb)

    # ---------- N = 10^7 ----------
    counts_7 = klassifiziere(quads_all, p_max=N_VOLL)
    drucke_tabelle(f"N = 10^7 = {N_VOLL:,}", counts_7)

    # ---------- Verifikation gegen Nutzerzahlen ----------
    print("\n  VERIFIKATION gegen Nutzerzahlen (N=10^7):")
    nutzer = {11: (161, 157), 101: (148, 140), 191: (140, 152)}
    ok = True
    for r in (11, 101, 191):
        a_ref, b_ref = nutzer[r]
        a_got = counts_7[r]['ABCE']
        b_got = counts_7[r]['CEAB']
        match = "✓" if (a_got == a_ref and b_got == b_ref) else "✗ ABWEICHUNG"
        print(f"  r={r:3d}: erwartet ABCE={a_ref}, CEAB={b_ref} | "
              f"berechnet ABCE={a_got}, CEAB={b_got}  {match}")
        if a_got != a_ref or b_got != b_ref:
            ok = False
    print(f"  Gesamt: {'Alle Werte bestätigt ✓' if ok else 'Abweichungen vorhanden ✗'}")

    alle_counts = [counts_7]
    labels = ['N = 10⁷']

    if COMPUTE_LARGE:
        # ---------- N = 10^8 ----------
        counts_8 = klassifiziere(quads_all, p_max=N_GROSS)
        drucke_tabelle(f"N = 10^8 = {N_GROSS:,}", counts_8)
        alle_counts.append(counts_8)
        labels.append('N = 10⁸')

        # ---------- Fenster 72e6 < p ≤ 77e6 ----------
        P_LOW, P_HIGH = 72_000_000, 77_000_000
        counts_f = klassifiziere(quads_all, p_min=P_LOW, p_max=P_HIGH)
        drucke_tabelle(f"Fenster 72·10^6 < p ≤ 77·10^6", counts_f)
        alle_counts.append(counts_f)
        labels.append('Fenster 72–77 Mio')

        # ---------- S-Wert für N=10^8 ----------
        D_11  = counts_8[11]['ABCE']  - counts_8[11]['CEAB']
        D_101 = counts_8[101]['ABCE'] - counts_8[101]['CEAB']
        D_191 = counts_8[191]['ABCE'] - counts_8[191]['CEAB']
        S = D_11 + D_101 - D_191
        Q = sum(counts_8[r]['ABCE'] + counts_8[r]['CEAB'] for r in (11, 101, 191))
        print(f"\n  S (N=10^8) = D_11 + D_101 - D_191 = {D_11:+} + {D_101:+} - ({D_191:+}) = {S:+}")
        print(f"  Q = {Q:,},  √Q = {math.sqrt(Q):.3f},  S/√Q = {S/math.sqrt(Q):+.4f}")

    # ---------- Plot ----------
    outpath = os.path.join(workspace, "eabc_wigner_zellen_plot.png")
    erstelle_plot(alle_counts, labels, outpath)

    print("\n  Fertig.")


if __name__ == "__main__":
    main()
