# RH + Planck „beides“: Frequenz- und Wellenlängen-Darstellung

Dieses Script macht genau das, was du mit „beides“ meinst:

- **Planckkurven über Frequenz \(\nu\)** (log-y)
- **Planckkurven über Wellenlänge \(\lambda\)** (log-y, und zusätzlich normalisiert)
- optional: **dein RH-Template** als Score-Heatmap **Score(T, N)**,
  einmal für **\(\nu\)**-Spektren und einmal für **\(\lambda\)**-Spektren.

## Dateien

- `rh_planck_both_axes.py`
- Input:
  - `zeros6.gz` (Riemann-Nullstellen \(\gamma_n\))
  - optional `spectrum.csv`:
    - entweder: `nu_Hz,S`
    - oder: `lambda_m,S`
    - oder: 2 Spalten ohne Header (Heuristik entscheidet)

## Quickstart

### 1) Nur klassische Planckkurven (ohne Spektrum)

```bash
python rh_planck_both_axes.py --zeros zeros6.gz --out out_rh_both
```

Erzeugt im Output:

- `planck_over_nu_log.png`
- `planck_over_lambda_log.png`  (bis 300 µm, damit der Abfall sichtbar ist)
- `planck_over_lambda_normalized.png` („klassische“ Lehrbuchform)

### 2) Mit Spektrum + Overlays

```bash
python rh_planck_both_axes.py --zeros zeros6.gz --spectrum spectrum.csv --out out_rh_both
```

Erzeugt zusätzlich:

- `spectrum_overlay_nu.png` **oder** `spectrum_overlay_lambda.png`

### 3) Mit Score(T,N)-Heatmap

```bash
python rh_planck_both_axes.py \
  --zeros zeros6.gz --spectrum spectrum.csv --do_heatmap \
  --Tmin 100 --Tmax 2000 --Tpoints 260 \
  --Nlist 200,500,1000,2000,5000,10000 \
  --out out_rh_both
```

Erzeugt:

- `score_heatmap_nu.png` oder `score_heatmap_lambda.png`

## Warum dein „hinterer Abfall“ vorher fehlte
Bei 273–1000 K liegt der Peak im IR. Im Bereich bis 30 µm siehst du im Fern-IR nur den Beginn des Rayleigh-Jeans-Flügels (\(\sim \lambda^{-4}\)).
Deshalb erweitert das Tool den Bereich im log-Plot standardmäßig bis 300 µm.

