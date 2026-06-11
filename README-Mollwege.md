# Mollwege.py – Mollweide-Visualisierung (Primzahlen mod 12)

Primzahlen werden nach Restklassen mod 12 (E/A/B/C/P) auf die Sphäre abgebildet (Hopf), in die Mollweide-Ebene projiziert und als animierte Frames gerendert.

## Abhängigkeiten

```bash
pip install -r requirements-mollwege.txt
```

## Gamma-Datei (Zeta-Nullstellen)

Optional: **zeros_gamma.txt** – eine Zeile pro Imaginärteil γ einer Riemann-Zeta-Nullstelle ρ = 1/2 + iγ (in aufsteigender Reihenfolge).  
Fehlt die Datei, werden synthetische Gammas verwendet (nur zur Demo, keine echten Nullstellen).

Zum Erzeugen von `zeros_gamma.txt` z.B. in **Zeta_Funktion_33_Nullpunkt.ipynb** die berechneten γ-Werte zeilenweise in eine Textdatei schreiben.

## Nutzung

```bash
# Standard (1000 Primzahlen, Frames in ./frames/, zeros_gamma.txt)
python Mollwege.py

# Weniger Frames, anderer Ordner, eigene Gamma-Datei
python Mollwege.py -n 500 -o meine_frames -g meine_gammas.txt

# Alle Optionen anzeigen
python Mollwege.py --help
```

## Video erzeugen

```bash
ffmpeg -framerate 30 -i frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p eabc_overlay.mp4
```
