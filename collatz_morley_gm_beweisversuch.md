# G_M: experimenteller Kandidat für einen signierten Morley-Krümmungssensor

**Stand:** Juni 2026  
**Branch:** `collatz/morley-phi-sensor` (nach Merge PR #48 auf `main`)  
**Artefakte:** `collatz_morley_tm_numerik.py m3`, `collatz_morley_m3_beweisversuch.json`, `m2-correlations` → `collatz_morley_m2_correlations.json`

---

## Epistemische Einordnung

> **G_M ist derzeit ein experimenteller Kandidat** für einen signierten Morley-Krümmungssensor.

> **M2 testet, ob sign(G_M) mit sign(K_G) korreliert.**

> G_M wird als signierter Morley-Fehler eingeführt. Er ist kein bewiesener Krümmungssensor, sondern eine experimentelle Testgröße, mit der geprüft wird, ob der geodätische Morley-Operator zwischen positiver und negativer Krümmung unterscheiden kann.

Dieses Dokument folgt dem Tao-Stil: Definition → M1 → M2 → M3 → ehrliche Grenzen. Status: **Experiment / Conjecture** — offen, ob $G_M$ definitionsunabhängig tragfähig ist.

Es werden **keine** Collatz-Aussagen abgeleitet. Morley ist ein eigenständiges Testmodul für $T_M^{(g)}$ auf riemannschen Flächen.

---

## Recap: M1 → M2 → M3

| Stufe | Frage | Befund (numerisch) |
|-------|-------|-------------------|
| **M1** | Stimmen vier Realisierungen von $T_M^{(g)}$ bis 2. Ordnung überein? | Ja — Abstände $d_{ij}=O(\varepsilon^3)$ |
| **M2a** | Korreliert $F_M$ mit Krümmungsstärke? | $F_M\approx 0$ bei $K{=}0$; $F_M>0$ bei $K{\neq}0$; $F_M(S^2)\approx F_M(H^2)$ |
| **M2b** | Korreliert $\operatorname{sign}(G_M)$ mit $\operatorname{sign}(K_G)$? | $\operatorname{sign}(G_M)=\operatorname{sign}(K_G)$ auf $S^2/H^2$; $G_M\approx 0$ auf $\mathbb{R}^2$ |
| **M3** | Skalieren $F_M$, $G_M$ wie $c\,|K|^\alpha A^\beta$ bzw. $c\,K^\alpha A^\beta$? | Siehe Tabellen unten |

**Definitionen** (kanonische Variante `local_chart`):

$$F_M(\Delta) = \sum_i (\theta_i^M - \pi/3)^2, \qquad
G_M(\Delta) = \sum_i (\theta_i^M - \pi/3).$$
$\varepsilon$-Familie: $\{0{,}05,\;0{,}08,\;0{,}12,\;0{,}18,\;0{,}25,\;0{,}35\}$.

---

## M3 — Numerische Ergebnisse

### Ebene $K_G=0$ (Kontrolle, aus Fit ausgeschlossen)

| Größe | Max über $\varepsilon$ | Lesart |
|-------|------------------------|--------|
| $F_M$ | $2{,}2\times 10^{-30}$ | $\approx 0$ (Morley-Satz) |
| $|G_M|$ | $6{,}7\times 10^{-16}$ | $\approx 0$ |

### Dualer log-log-Fit (gepoolt $S^2 + H^2$, $n=12$)

| Testgröße | Modell | $c$ | $\alpha$ | $\beta$ | $R^2$ |
|-----------|--------|-----|----------|---------|-------|
| $F_M$ | $c_F\,|K_G|^\alpha A^\beta$ | $3{,}88\times 10^{-4}$ | $\approx 0$ | $\approx 2{,}00$ | $0{,}9998$ |
| $G_M$ | $c_G\,K_G^\alpha A^\beta$ | $3{,}41\times 10^{-2}$ | $\approx 0$ | $\approx 1{,}00$ | $0{,}9998$ |

**Lesart:** $F_M \propto A^2$ (quadratischer Morley-Fehler); $G_M \propto A$ (linearer, **signierter** Fehler). Die Exponenten $\alpha$ sind auf den Kontrollflächen mit $|K_G|\equiv 1$ **nicht identifizierbar** — nur $\beta$ ist belastbar.

### $G_M$-Vorzeichen-Cross-Check

| Test | Ergebnis |
|------|----------|
| Stabil über alle $\varepsilon$ | ja |
| $G_M>0$ auf $S^2$ (jedes $\varepsilon$) | ja |
| $G_M<0$ auf $H^2$ (jedes $\varepsilon$) | ja |
| $G_M\approx 0$ auf $\mathbb{R}^2$ | ja |
| $\operatorname{sign}(G_M)=\operatorname{sign}(K_G)$ für alle vier M1-Varianten | ja |

---

## Walter-Sensor $W_M$ — unabhängige Flächeninformation (M2-Erweiterung)

Morley ($F_M$, $G_M$) misst **Winkelstruktur** des Morley-Dreiecks; Marion Walter misst **Flächenstruktur** des zentralen Hexagons $H_W(\Delta)$:

$$W_M(\Delta) = \frac{\mathrm{Area}(H_W(\Delta))}{\mathrm{Area}(\Delta)} - \frac{1}{10}.$$

| Aspekt | $G_M$ (Winkel) | $W_M$ (Fläche) |
|--------|----------------|----------------|
| Euklidische Referenz | Morley-Satz → $\approx 0$ | Marion-Walter-Satz → $\approx 0$ |
| Gekrümmt | chart-basierter Morley-Operator | chart-nähe Trisektion + Cevianen |
| Epistemischer Status | experimenteller Vorzeichen-Kandidat | experimenteller Flächen-Kandidat |
| Theorem auf $S^2/H^2$? | nein | **nein** |

**Numerik (Median, `local_chart`, `collatz_morley_m2_sensor.json`):**

| Raum | $W_M$ |
|------|-------|
| $\mathbb{R}^2$ | $\approx 0$ |
| $S^2$ | $\approx -6{,}8\times 10^{-4}$ |
| $H^2$ | $\approx +6{,}8\times 10^{-4}$ |

**Vorzeichen-Test:** $\operatorname{sign}(G_M)=\operatorname{sign}(K_G)$ (ja), aber $\operatorname{sign}(W_M)=-\operatorname{sign}(G_M)$ auf $S^2/H^2$ — $W_M$ ist **nicht** parallel zum Morley-Vorzeichensensor, sondern liefert orthogonale Flächeninformation. Das spricht **für** die Unabhängigkeit der Sensoren, nicht für einen einfachen gemeinsamen Krümmungsindikator.

**Implementierung:** `walter_hexagon_area`, `walter_form_wm` in `collatz_morley_tm_numerik.py`; Tests in `tests/test_morley_walter_wm.py`.

---

## Sensorvektor $\Phi_M = (G_M, W_M)$ (Reviewer-Synthese, PR #48/49)

> **Boxed (Kernidee):**  
> Nicht $G_M \approx W_M$, sondern **$\Phi_M = (G_M, W_M)$** als zweidimensionaler lokaler Geometriezustand.

| Komponente | Defekt-Typ | Euklidische Nullreferenz |
|------------|------------|--------------------------|
| $G_M$ | Winkel-Formfehler $\Delta\theta$ | Morley-Satz |
| $W_M$ | Flächenfehler $\Delta A$ | Marion-Walter-Satz ($1/10$) |
| $F_M$ | Winkel-Formfehler (Betrag) | Morley-Satz — Stärke, kein Vorzeichen |

**Terminologie (Korrektur):** $F_M$ ist ein **Winkel-Formfehler** (quadratische Summe), **kein** Flächensignal. $W_M$ ist der **Flächensensor**.

**Anti-Parallelität auf $S^2/H^2$:** $\operatorname{sign}(W_M)=-\operatorname{sign}(G_M)=-\operatorname{sign}(K_G)$ — das ist ein **Merkmal**, kein Bug. $W_M$ misst Flächenprojektionsdefekt, $G_M$ Winkeldefekt; analog zu $K$ vs. $H$ in der DG.

**Kein Triple-Konsens:** $(F_M, G_M, W_M)$ tragen nicht dieselbe Information — das spricht für einen echten Sensorvektor.

| Objekt | Status |
|--------|--------|
| $F_M$ | experimenteller Morley-Sensor (Betrag) |
| $G_M$ | sign-korrelierter Krümmungs**kandidat** |
| $W_M$ | unabhängiger Walter-Flächensensor |
| Zusammenhang mit $K_G$ | empirisch, nicht bewiesen |

**Test A (Pearson-$\rho$):** `collatz_morley_tm_numerik.py m2-correlations` — $\rho(G_M,K_G)$, $\rho(W_M,K_G)$, $\rho(F_M,K_G)$ und Kreuzkorrelationen; siehe `collatz_morley_m2_correlations.json`.

---

## Babylon-Kalibrierung (3-4-5)

**Referenz:** [`Tri - Okto.tex`](Tri%20-%20Okto.tex), Definition *Orthogonal calibration*:

$$3^2 + 4^2 = 5^2$$

liefert die elementare orthogonale Referenz. Der erweiterte Oktogon-Defekt $D_8^{\mathrm{orth}}$ nutzt $|a^2+b^2-c^2|$ als Orthogonalitätsdefekt — analog prüft `babylon_orthogonalize` $| \|\mathrm{leg}_3\|^2 + \|\mathrm{leg}_4\|^2 - \|\mathrm{hyp}\|^2 |$.

**Zuordnung:**

| Skala | Sensor | Lesart |
|-------|--------|--------|
| 3 | $G_M$ | Morley-Winkeldefekt |
| 4 | $W_M$ | Walter-Flächendefekt |
| 5 | $\|3\alpha\hat u_3 + 4\beta\hat u_4\|$ | kombinierte Babylon-Sensorlänge |

$\Phi_M = (G_M, W_M)$ wird per G_M-first Gram-Schmidt in orthogonale Babylon-Achsen $(\hat u_3, \hat u_4)$ zerlegt. Da $\rho(G_M,W_M)\approx -1$ (Test A), wäre eine Identitätsannahme $G_M \propto W_M$ **falsch** — Babylon erzwingt orthogonale Beine statt Kollinearität.

**CLI:** `collatz_morley_tm_numerik.py m2-babylon` → `collatz_morley_m2_babylon.json`

**Epistemisch:** heuristisch / experimentell — **kein Theorem**, Link zu Tri-Okto, **nicht** Teil des Collatz-Beweisprogramms.

**Test C (Stub):** Orientiertes $W_M^{\mathrm{or}} = \sigma(\Delta)\,( \mathrm{Area}(H_W)/\mathrm{Area}(\Delta) - 1/10)$ mit $\sigma(\Delta)=\operatorname{sign}(\mathrm{Area}(\Delta))$ — `walter_form_wm_oriented` (nächste PR).

---

## Argumente **für** den experimentellen Kandidaten

1. **Vorzeichenkonsistenz:** $\operatorname{sign}(G_M)=\operatorname{sign}(K_G)$ auf $S^2$ und $H^2$ bei jedem $\varepsilon$ der Familie — nicht nur im Median.
2. **Variantenrobustheit:** Das Vorzeichenmuster gilt für alle vier M1-Realisierungen (`geodesic_angles`, `local_chart`, `exp_euclidean`, `parallel_transport`).
3. **Skalierung:** $G_M \propto A$ mit $R^2\approx 0{,}9998$ — der signierte Fehler skaliert linear mit der Dreiecksfläche, orthogonal zu $F_M \propto A^2$.
4. **Kombinierte Testgrößen:** $(F_M, G_M)$ trennt Stärke ($|K|$-Proxy via $F_M$) und Vorzeichen ($G_M$) — stärker als $F_M$ allein.

---

## Argumente **gegen** / Grenzen

1. **Keine Eindeutigkeit:** $G_M=\sum(\theta_i^M-\pi/3)$ ist eine **willkürliche** erste Wahl; orientierte Varianten (Flächenvorzeichen, Umlauf, geodätischer Exzess) sind noch nicht getestet.
2. **Nur numerisch:** Keine analytische Herleitung, dass $G_M$ eine Krümmungsinvariante approximiert.
3. **Eine kanonische Variante im Fit:** Der Exponentenfit nutzt `local_chart`; Varianten-Invarianz von $c_G$, $\beta$ ist nicht gefittet.
4. **Konstantes $|K_G|$:** Auf $S^2/H^2$ mit $K_G=\pm 1$ ist $\alpha$ degeneriert — echte $|K|$-Skalierung braucht variable Krümmung oder mehr Modellräume.
5. **Hyperbolisches $H^2$:** Numerisches Patch-Modell, keine vollständige Levi-Civita-Implementierung.
6. **Orthogonale Sensoren:** $W_M$ (Fläche) und $G_M$ (Winkel) können entgegengesetzte Vorzeichen auf $S^2/H^2$ zeigen — kombinierte $(F_M, G_M, W_M)$-Analyse nötig, kein einzelner Krümmungsindikator.
7. **Keine globale Aussage:** M1-$O(\varepsilon^3)$-Evidenz rechtfertigt keine Aussage über globale Krümmungserkennung.

---

## Vergleich: κ-Zweig vs. Morley-Asymmetrie

Der κ-Zweig (arithmetische Kodierung, $\mathcal{L}_{\mathrm{arith}}^*$) testet **Invarianz unter BE-Verbot** und scheitert dort strukturell. Der Morley-Zweig testet **geometrische Definitionsrobustheit** von $T_M^{(g)}$ und zeigt numerische Stabilität über Varianten. Die Asymmetrie ist beabsichtigt: arithmetische Symbolik und Riemann-Geometrie sind **parallele Spuren**, keine Vereinheitlichung. $G_M$ als experimenteller Vorzeichen-Kandidat adressiert eine Frage, die κ nicht stellt — und umgekehrt.

---

## Fazit (boxed)

> **Starke Evidenz wäre:** analytischer Zusammenhang $G_M = c\,K_G\,A + O(A^2)$ aus einer eindeutigen $T_M^{(g)}$-Definition; Invarianz von $G_M$ unter allen M1-Varianten **mit gleichem** $c_G$; Bestätigung auf variablem $K_G$ (z. B. Ikosaeder-20-Flächen mit unterschiedlichen lokalen Krümmungen); orientierte $G_M$-Definition, die unter Drehung/Spiegelung konsistent bleibt.
>
> **Offen bleibt:** Eindeutigkeit der Testgröße; analytische Herleitung; globale Aussage jenseits kleiner $\varepsilon$-Dreiecke; Kopplung an Collatz oder κ — **bewusst nicht beansprucht**.

---

## Nächste Schritte

1. ~~**Test A:** Pearson-Korrelationen $(F_M, G_M, W_M)$ vs. $K_G$~~ — implementiert (`m2-correlations`).
2. ~~**Babylon-Kalibrierung (3-4-5)**~~ — implementiert (`m2-babylon`, `babylon_orthogonalize`).
3. **Test B:** Vollständige PCA auf $(F_M, G_M, W_M)$ — 2 vs. 3 Freiheitsgrade?
4. **Test C:** Orientiertes $W_M^{\mathrm{or}}$ (`walter_form_wm_oriented`).
4. **Orientiertes $G_M$:** Testgröße mit explizitem Flächenvorzeichen / Umlaufsinn.
5. **Ikosaeder (20 Flächen):** Systematischer M3-Lauf auf allen Dreiecksflächen.
6. **Varianten-Invarianz von $G_M$:** $|c_G|$ und $\beta$ über alle vier Realisierungen.
7. **Variable Krümmung:** Modellräume mit unterschiedlichem $|K_G|$ zur Identifikation von $\alpha$.

---

*Siehe auch:* `collatz_morley_stufen_m.md` (Stufen-Roadmap), `collatz_morley_m2_sensor.json` (M2-Rohdaten).
