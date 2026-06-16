# Reproduzierbare Pipeline: `zeros6.npy`

## Was ist `zeros6.npy`?

Eindimensionales `float64`-Array der Länge **2 001 052** mit den Imaginärteilen

\[
\gamma_1 < \gamma_2 < \cdots < \gamma_N
\]

der nicht-trivialen Nullstellen von \(\zeta(s)\) auf der kritischen Gerade \(\Re(s)=\tfrac12\) (also \(\rho_n = \tfrac12 + i\gamma_n\)).

Referenzwerte (Workspace-Stand):

| Größe | Wert |
|--------|------|
| \(N\) | 2 001 052 |
| \(\gamma_1\) | ≈ 14.134725142 |
| \(\gamma_N\) | ≈ 1 132 490.658714 |

Verwendung im Repo u. a. in `Resonanz.py`, `Resonanz.ipynb`, `BA heiliger Gral.py`, `EABC Kontinuum.py` (Projektroot).

**Kein Beweis der Riemannschen Vermutung:** Die Datei ist nur ein **numerischer Datensatz** für Korrelations- und Spektraltests.

## Trennung vom Hoffbauer-Sieb (0,9992)

Die empirische Restmetrik **≈ 0,9992** aus `Hoffbauer Sieb.py` / `hoffbauer_geometrie_test.csv` ist eine **geometrische Sieb-Statistik auf Primzahlen**, nicht \(\zeta(s)\) und nicht \(\gamma_n\).

Diese Pipeline liefert ausschließlich die \(\gamma\)-Folge. Ein direkter Nachweis „\(\gamma \to 0{,}9992\)“ ist damit **nicht** gemeint und **nicht** Teil dieses Ordners.

## Reproduktion

1. **Quelle beschaffen** (nicht im Git, typisch mehrere MB):

   - Projektkonvention laut `zeros6_meta.json`: Textdatei **`zeros6.gz`** (eine reelle Zahl pro Zeile), nach `riemann_zeros_pipeline/zeros6.gz` legen.
   - Öffentliche Referenztabellen (Montgomery–Odlyzko-Nullstellen):  
     [Odlyzko zeta tables — zeros1](https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros1)  
     Bei Binärformat (rohe `float64`): `generate_zeros6.py --binary --input /pfad/zur/datei`.
   - Plausibilität einzelner \(\gamma_n\): [LMFDB Riemann zeros](https://www.lmfdb.org/riemann/).

2. **Erzeugen** (schreibt standardmäßig ins **Repo-Root**, damit bestehende Skripte `zeros6.npy` finden):

   ```bash
   cd riemann_zeros_pipeline
   python generate_zeros6.py
   # oder: python generate_zeros6.py --input /pfad/zu/zeros6.gz
   ```

3. **Prüfen** (SHA256 + Inhalt):

   ```bash
   python verify_zeros6.py
   ```

## Erwartete SHA256-Prüfsummen

| Datei | SHA256 |
|--------|--------|
| `zeros6.npy` (Repo-Root) | `9062f83edd01af3803e848c3ff894a328e595f46bd2085b6356e3d30acc367be` |
| `zeros6.npz` (Repo-Root) | `09f4675cfb0cac2211f1bff9844cd863a59f854ec99bd11308225db54ce16054` |

`zeros6.npy` ist ~16 MB und wird **nicht** versioniert (siehe `.gitignore`). Im Git liegt u. a. `zeros6.npz` mit gleicher \(\gamma\)-Folge; nach `np.load(... )['zeros']` identisch.

## Claim → Datensatz → Skript → Reproduktion

| Schritt | Artefakt |
|---------|-----------|
| Claim | „Wir nutzen die ersten 2 001 052 bekannten \(\gamma_n\) wie in `zeros6.npy`.“ |
| Datensatz | `zeros6.npy` / `zeros6.npz` |
| Skript | `generate_zeros6.py` |
| Reproduktion | `python generate_zeros6.py` + `python verify_zeros6.py` |
