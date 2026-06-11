# SageMath 10.9 — Korrelation zwischen Interferenzmodell und Riemann-Nullstellen (gamma)
#
# Terminal:
#   sage "Resonanz.py"
#   sage -python "Resonanz.py"
#
# Laufende Sage-Sitzung (Notebook oder sage:) — gleicher Kernel wie „im Hintergrund“:
#   load("/Pfad/zu/Resonanz.py")     einmal ausführen
#   attach("/Pfad/zu/Resonanz.py")   bei jeder Speicherung der Datei neu laden (Entwicklung)
#
# load: einmalig. attach: auto-Neuladen nach Änderung am File (beachte: kann beim Arbeiten neu auslösen).
#
# Reines „python Resonanz.py“ scheitert am fehlenden sage-Modul (Absicht).

import sys
from pathlib import Path

try:
    from sage.all import cos, var  # type: ignore[import-untyped]
except ImportError:
    print(
        "SageMath 10.9 ist nicht geladen (Modul sage fehlt).\n"
        "Aufruf z.B.:\n"
        '  sage "Resonanz.py"\n'
        '  sage -python "Resonanz.py"\n',
        file=sys.stderr,
    )
    sys.exit(2)


def korrelations_test_riemann(nullstellen_liste, schwellenwert=0.7):
    """
    Vergleicht ein einfaches Resonanzmodell |cos(t/90) + cos(t/180)| mit den
    gegebenen Ordinaten gamma der Zeta-Nullstellen auf der kritischen Geraden.

    nullstellen_liste: Liste oder Folge reeller gamma-Werte
    schwellenwert: Schwelle für „Hochenergie“-Treffer
    """
    if not nullstellen_liste:
        print("Keine Nullstellen übergeben (leere Liste).")
        return float("nan")

    # Symbolisches Modell in t (reellwertiger „Energiedichte“-Proxy)
    t = var("t")
    f_resonanz = abs(cos(t / 90) + cos(t / 180))

    treffer = 0
    n = len(nullstellen_liste)

    print(f"Starte Test für {n} Nullstellen (Schwelle {float(schwellenwert)})...")

    for gamma in nullstellen_liste:
        wert = f_resonanz.subs({t: gamma})
        try:
            lokale_energie = float(wert)
        except (TypeError, ValueError):
            lokale_energie = float(wert.numerical_approx())

        if lokale_energie > schwellenwert:
            treffer += 1

    quote = float(treffer) / n * 100.0
    print(f"Ergebnis: {quote:.2f}% der Nullstellen liegen oberhalb der Schwelle.")
    return quote


if __name__ == "__main__":
    # Optional: gamma aus zeros6.npy (Projektroot), höchstens 10000 Werte; sonst Demo-Liste
    npy = Path(__file__).resolve().parent / "zeros6.npy"
    if npy.is_file():
        import numpy as np

        gammas = np.asarray(np.load(npy), dtype=float).ravel()[:10000].tolist()
        print(f"Lade gamma aus {npy.name}, verwendet: {len(gammas)} Werte.")
    else:
        gammas = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178]
        print("Keine zeros6.npy gefunden — Demo mit wenigen Literatur-gamma.")

    _ = korrelations_test_riemann(gammas)
