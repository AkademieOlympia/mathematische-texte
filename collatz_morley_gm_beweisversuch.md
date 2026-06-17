# Beweisversuch: G_M als signierter Morley-Krümmungssensor (experimentell)

**Stand:** Juni 2026  
**Branch:** `collatz/morley-m3-gm-beweisversuch`  
**Artefakte:** `collatz_morley_tm_numerik.py m3`, `collatz_morley_m3_beweisversuch.json`

---

## Epistemische Einordnung

> **Kein Beweis.** Dieses Dokument ist ein **experimenteller Beweisversuch** im Tao-Stil: Definition → M1 → M2 → M3 → ehrliche Grenzen. Status: **Experiment / Conjecture** — offen, ob $G_M$ als geometrischer Krümmungssensor definitionsunabhängig tragfähig ist.

Es werden **keine** Collatz-Aussagen abgeleitet. Morley ist ein eigenständiges Testmodul für $T_M^{(g)}$ auf riemannschen Flächen.

---

## Recap: M1 → M2 → M3

| Stufe | Frage | Befund (numerisch) |
|-------|-------|-------------------|
| **M1** | Stimmen vier Realisierungen von $T_M^{(g)}$ bis 2. Ordnung überein? | Ja — Abstände $d_{ij}=O(\varepsilon^3)$ |
| **M2a** | Misst $F_M$ Krümmungsstärke? | $F_M\approx 0$ bei $K{=}0$; $F_M>0$ bei $K{\neq}0$; $F_M(S^2)\approx F_M(H^2)$ |
| **M2b** | Misst $G_M$ das Vorzeichen von $K$? | $\operatorname{sign}(G_M)=\operatorname{sign}(K_G)$ auf $S^2/H^2$; $G_M\approx 0$ auf $\mathbb{R}^2$ |
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

| Sensor | Modell | $c$ | $\alpha$ | $\beta$ | $R^2$ |
|--------|--------|-----|----------|---------|-------|
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

## Argument **für** den signierten Sensor

1. **Vorzeichenkonsistenz:** $\operatorname{sign}(G_M)=\operatorname{sign}(K_G)$ auf $S^2$ und $H^2$ bei jedem $\varepsilon$ der Familie — nicht nur im Median.
2. **Variantenrobustheit:** Das Vorzeichenmuster gilt für alle vier M1-Realisierungen (`geodesic_angles`, `local_chart`, `exp_euclidean`, `parallel_transport`).
3. **Skalierung:** $G_M \propto A$ mit $R^2\approx 0{,}9998$ — der signierte Fehler skaliert linear mit der Dreiecksfläche, orthogonal zu $F_M \propto A^2$.
4. **Kombinierter Sensor:** $(F_M, G_M)$ trennt Stärke ($|K|$-Proxy via $F_M$) und Vorzeichen ($G_M$) — stärker als $F_M$ allein.

---

## Argument **gegen** / Grenzen

1. **Keine Eindeutigkeit:** $G_M=\sum(\theta_i^M-\pi/3)$ ist eine **willkürliche** erste Wahl; orientierte Varianten (Flächenvorzeichen, Umlauf, geodätischer Exzess) sind noch nicht getestet.
2. **Nur numerisch:** Kein analytischer Beweis, dass $G_M$ eine Krümmungsinvariante approximiert.
3. **Eine kanonische Variante im Fit:** Der Exponentenfit nutzt `local_chart`; Varianten-Invarianz von $c_G$, $\beta$ ist nicht gefittet.
4. **Konstantes $|K_G|$:** Auf $S^2/H^2$ mit $K_G=\pm 1$ ist $\alpha$ degeneriert — echte $|K|$-Skalierung braucht variable Krümmung oder mehr Modellräume.
5. **Hyperbolisches $H^2$:** Numerisches Patch-Modell, keine vollständige Levi-Civita-Implementierung.
6. **Kein Theorem:** M1-$O(\varepsilon^3)$-Evidenz rechtfertigt keine Aussage über globale Krümmungsmessung.

---

## Vergleich: κ-Zweig vs. Morley-Asymmetrie

Der κ-Zweig (arithmetische Kodierung, $\mathcal{L}_{\mathrm{arith}}^*$) testet **Invarianz unter BE-Verbot** und scheitert dort strukturell. Der Morley-Zweig testet **geometrische Definitionsrobustheit** von $T_M^{(g)}$ und zeigt numerische Stabilität über Varianten. Die Asymmetrie ist beabsichtigt: arithmetische Symbolik und Riemann-Geometrie sind **parallele Spuren**, keine Vereinheitlichung. $G_M$ als Vorzeichensensor adressiert eine Frage, die κ nicht stellt — und umgekehrt.

---

## Fazit (boxed)

> **Starke Evidenz wäre:** analytischer Zusammenhang $G_M = c\,K_G\,A + O(A^2)$ aus einer eindeutigen $T_M^{(g)}$-Definition; Invarianz von $G_M$ unter allen M1-Varianten **mit gleichem** $c_G$; Bestätigung auf variablem $K_G$ (z. B. Ikosaeder-20-Flächen mit unterschiedlichen lokalen Krümmungen); orientierte $G_M$-Definition, die unter Drehung/Spiegelung konsistent bleibt.
>
> **Offen bleibt:** Eindeutigkeit der Sensordefinition; analytischer Beweis; globale Aussage jenseits kleiner $\varepsilon$-Dreiecke; Kopplung an Collatz oder κ — **bewusst nicht beansprucht**.

---

## Nächste Schritte

1. **Orientiertes $G_M$:** Sensor mit explizitem Flächenvorzeichen / Umlaufsinn definieren und gegen die einfache Summe testen.
2. **Ikosaeder (20 Flächen):** Systematischer M3-Lauf auf allen Dreiecksflächen des regulären Ikosaeders (`eabc_icosahedron_test.py`).
3. **Varianten-Invarianz von $G_M$:** Nicht nur Vorzeichen, sondern $|c_G|$ und $\beta$ über alle vier Realisierungen vergleichen.
4. **Variable Krümmung:** Modellräume mit unterschiedlichem $|K_G|$ zur Identifikation von $\alpha$.

---

*Siehe auch:* `collatz_morley_stufen_m.md` (Stufen-Roadmap), `collatz_morley_m2_sensor.json` (M2-Rohdaten).
