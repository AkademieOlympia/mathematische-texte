# **Morley-Erweiterung: diskreter Krümmungs- und Konformitätssensor**

**Branch:** `collatz/kappa-invarianz-stufe3` · **PR:** #41 (Draft)  
**Stand:** Juni 2026  
**Kein Collatz-Beweis.** Geometrische Brücke — epistemisch abgegrenzt von Stufe-3-Kern (κ-Invarianz).

**Methodik:** Tao-Stil — Definition / Zeuge / Experiment / Theorem / Conjecture strikt trennen (`collatz_formalisierung_tao_stil.md`).

**Vorgänger:** `MorleyWalter.tex`, `collatz_dc_morley_walter.pdf`, `collatz_schlussartikel_arxiv.tex` (§ Morley–Walter-Geometrie, illustrativ).

---

## Kernthese (revidiert)

> **Boxed (Leitsatz):**  
> *Morley nicht als Metrik, sondern als diskreten Krümmungs- und Konformitätssensor.*

Morley liefert **keine** neue Riemann-Metrik $g_{ij}$ auf einer glatten Fläche. Der Morley-Kern lebt auf **Dreiecken** (Winkel, Flächenverhältnis), nicht im Tangentialraum — näher an **Regge-Kalkül**, **diskreter Gauß-Krümmung**, **Voronoi/Delaunay-Triangulierungen**, **FEM-Netzen** und **konformen Triangulierungen** als an einer infinitesimalen Metrik.

---

## Kontext: Morley vs. Riemann

| Ebene | Objekt | Lokalität | Typische Literatur |
|-------|--------|-----------|-------------------|
| **Riemann** | Metrik $g_{ij}$, Christoffel, Krümmungstensor | infinitesimal (Tangentenraum) | glatte Differentialgeometrie |
| **Morley** | Winkeldrittel, gleichseitiges Morley-Dreieck $\mathrm{Mor}(\Delta)$ | **winkel-lokal** auf einem Dreieck | klassische euklidische Geometrie |
| **Diskret** | Defizitwinkel, Regge-Winkel, dualer Volumenfluss | **triangulations-lokal** | Regge (1961), discrete exterior calculus |

**Konsequenz:** Eine naive „Morley-Metrik“ $K_M$ als einzelne skalare Größe war zu grob. Stattdessen trennen wir zwei **orthogonal lesbare** Sensoren:

1. **Morley-Form** $F_M$ — Gleichseitigkeitsabweichung (rein winkelbasiert),
2. **Morley-Skala** $S_M$ — Größenänderung (Flächenverhältnis).

Die kombinierte **Morley-Krümmung** $K_M$ gewichtet beide mit Parametern $\alpha,\beta$ und einem euklidischen Referenzwert $S_0$.

---

## Definitionen

Sei $\Delta = ABC$ ein nicht-degeneriertes Dreieck in der euklidischen Ebene (oder auf einer konformen Karte einer Riemann-Fläche). Die **Morley-Winkel** $\theta_i^M$ ($i=1,2,3$) seien die Innenwinkel des Morley-Dreiecks $\mathrm{Mor}(\Delta)$ — das gleichseitige Dreieck aus den Schnittpunkten benachbarter Winkeldrittelnden (Morley-Satz).

### Morley-Form (Gleichseitigkeitsabweichung)

$$F_M(\Delta) \;:=\; \sum_{i=1}^{3} \left(\theta_i^M - \frac{\pi}{3}\right)^2.$$

- $F_M = 0$ genau dann, wenn $\mathrm{Mor}(\Delta)$ perfekt gleichseitig ist (immer im euklidischen Fall — der Sensor misst **Abweichung von der Referenzform** unter Störung/Krümmung).
- Auf einer **krummen** Fläche (oder unter numerischem Rauschen) dient $F_M$ als **Formdefekt** des Morley-Kerns.

**Epistemisches Label:** **Definition** (definitorisch).

### Morley-Skala (Größenänderung)

$$S_M(\Delta) \;:=\; \frac{\operatorname{Area}\bigl(\mathrm{Mor}(\Delta)\bigr)}{\operatorname{Area}(\Delta)}.$$

**Ebenenformel** (klassisch, für Dreieck mit Winkeln $A,B,C$):

$$\operatorname{Area}\bigl(\mathrm{Mor}(\Delta)\bigr)
= \operatorname{Area}(\Delta)\cdot
\frac{\sin(A/3)\,\sin(B/3)\,\sin(C/3)}{\sin A\,\sin B\,\sin C}.$$

Damit

$$S_M(\Delta) = \frac{\sin(A/3)\,\sin(B/3)\,\sin(C/3)}{\sin A\,\sin B\,\sin C}.$$

- $S_M$ ist **dimensionslos** und invariant unter Ähnlichkeitstransformationen.
- Für gleichseitige Ausgangsdreiecke: $S_M = 1/3$ (bekannte Konstante).
- $S_M$ kodiert, wie stark der Morley-Kern **schrumpft** oder **wächst** relativ zur Hülle.

**Epistemisches Label:** **Definition** (definitorisch); Ebenenformel **Theorem** (klassisch).

### Kombinierte Morley-Krümmung

$$K_M(\Delta) \;:=\; \alpha\, F_M(\Delta) + \beta\,\bigl(S_M(\Delta) - S_0\bigr)^2,$$

wobei $\alpha,\beta \geq 0$ Gewichte und $S_0$ ein **euklidischer Referenzwert** ist (z. B. $S_0 = 1/3$ für gleichseitige Normalisierung, oder lokal gemittelt über eine Triangulierung).

**Epistemisches Label:** **Definition** (definitorisch). Die Wahl von $(\alpha,\beta,S_0)$ ist **Modellparameter**, keine geometrische Notwendigkeit.

---

## Verbindung zur Riemann-Fläche

Auf einer glatten orientierten Riemann-Fläche $(M,g)$ gilt für kleine geodätische Dreiecke mit Winkelsumme $\alpha+\beta+\gamma$ und Fläche $A$:

$$\alpha + \beta + \gamma = \pi + K_G\, A + O(A^2),$$

wobei $K_G$ die **Gauß-Krümmung** im Schwerpunkt ist (Gauß–Bonnet lokale Form).

### Offene asymptotische Frage

> **Conjecture (Morley–Gauß-Asymptotik):**  
> Existieren Konstanten $c_1, c_2, \ldots$ (abhängig von der Normalisierung von $K_M$), sodass für kleine Dreiecke
> $$K_M(\Delta) = c_1\, K_G + c_2\, K_G^2 + O(A)$$
> oder mindestens $K_M(\Delta) \sim c_1 K_G$ im Grenzfall $A\to 0$?

**Status:** **offen** — keine Beweise in diesem Dokument. Plausibel ist die Lesart:

- $F_M$ reagiert auf **Winkeldefizit** (erste Ordnung in $K_G A$),
- $S_M - S_0$ reagiert auf **Flächenverzerrung** (ebenfalls erste Ordnung, gekoppelt über konforme Faktorisierung).

**Epistemisches Label:** **Conjecture**.

---

## Komplexe Ebene und Morley-Konformitätsfeld

Auf der konformen Ebene mit $ds^2 = \lambda(z,\bar z)\,|dz|^2$ sei $\Delta = (z_1,z_2,z_3)$ ein Dreieck und $(m_1,m_2,m_3)$ die entsprechenden Morley-Eckpunkte.

### Diskrete konforme Ableitung

$$\mu_M \;:=\; \frac{m_2 - m_1}{z_2 - z_1}$$

(und analog auf den anderen Kanten) wirkt als **diskrete konforme Ableitung** entlang der Kanten — ein **Morley-Konformitätsfeld** auf der Triangulierung.

| Objekt | Kontinuierlich | Diskret (Morley) |
|--------|----------------|------------------|
| Konforme Struktur | $\lambda(z,\bar z)$ | Kantenquotienten $\mu_M$ |
| Quasikonform | Beltrami $\mu = \partial_{\bar z} f / \partial_z f$ | $\mu_M$ als Triangulationsfeld |
| Ableitung | $\partial_z$ | Differenzenquotient Morley/Hülle |

**Verbindungen (Literaturrahmen, nicht bewiesen hier):**

- Beltrami-Differentiale und quasikonforme Abbildungen,
- diskrete komplexe Analysis (Duffin, Mercat, circle packing),
- conformal weldings / discrete Riemann surfaces.

**Epistemisches Label:** **heuristisch** — $\mu_M$ als Konformitätssensor ist motiviert, aber nicht als Theorem etabliert.

---

## EABC-Perspektive (revidiert)

### Alt (zu naiv)

Morley-Kern als **zusätzlicher Punkt** $E$ in der EABC-Konfiguration — geometrisch suggestiv, aber strukturell unpräzise.

### Neu: Morley als Operator

Der Morley-Kern ist ein **Operator** auf Dreieckskonfigurationen:

$$(A,B,C) \;\longmapsto\; (E_A, E_B, E_C),$$

wobei $(E_A,E_B,E_C)$ die drei Eckpunkte von $\mathrm{Mor}(\Delta)$ sind (oder deren EABC-Kodierung unter DC-Interpretation).

**Transformation:**

$$T_M : \Delta \;\longmapsto\; \mathrm{Mor}(\Delta), \qquad \Delta_{k+1} = T_M(\Delta_k),$$

iterierbar: $\Delta_0 \to \Delta_1 \to \Delta_2 \to \cdots$

| Iteration | Lesart |
|-----------|--------|
| $T_M^1$ | ein Morley-Schritt (Winkeldrittel-Phasenextraktion) |
| $T_M^k$ | $k$-fache Morley-Renormierung |
| Grenzverhalten | offen — Attraktor? periodische Ordnung? |

**Brücke zu Collatz / κ (spekulativ, ehrlich):**

- DC-Blöcke als Dreieckskonfigurationen (`collatz_schlussartikel_arxiv.tex`, `collatz_dc_morley_walter.pdf`),
- EABC als mod-12-Klassifikation von Kanten/Spitzen,
- $T_M$ als **geometrischer Renormierungsoperator** parallel zu κ-Familien und Sprachverdünnung — **keine** bewiesene Implikation für Collatz-Trajektorien.

**Epistemisches Label:** **heuristisch / spekulativ** — operatorische Iteration als Forschungsbrücke, nicht als Lemma.

---

## Strukturbrücke $(M, g, T, \mathcal{M}_T)$

Die folgende Kette ordnet die Ebenen (kein Theorem, **organisatorisches Schema**):

$$\boxed{
\text{Riemann-Fläche } (M,g)
\;\to\;
\text{Triangulierung}
\;\to\;
\text{Morley-Operator } T_M
\;\to\;
\text{diskrete komplexe Dynamik}
}$$

| Komponente | Rolle |
|------------|-------|
| $(M,g)$ | glatte Hintergrundgeometrie (falls vorhanden) |
| Triangulierung | Diskretisierung; Regge-/Delaunay-Niveau |
| $T_M$ | Morley-Operator auf Dreiecken |
| $\mathcal{M}_T$ | Trajektorien unter Iteration $T_M^k$; Sensoren $F_M, S_M, K_M, \mu_M$ |

**Nähe zu:** diskreter Differentialgeometrie, Regge-Kalkül, diskreter konformer Geometrie — **nicht** zu einer neuen „Morley-Metrik“ im Riemannschen Sinn.

---

## Zwei Sensoren vs. naive $K_M$

| Größe | Misst | Sensitivität |
|-------|-------|--------------|
| $F_M$ | Winkeldefekt des Morley-Kerns | rein **angular** |
| $S_M$ | Flächenverhältnis Hülle/Kern | **metrisch** (Skala) |
| $K_M$ | gewichtete Kombination | tunable via $\alpha,\beta,S_0$ |

Die Trennung verhindert die Verwechslung „Morley = neue Metrik“. In der Praxis:

- **Formdominiert:** $F_M$ groß, $S_M \approx S_0$ → stark nichtgleichseitige Morley-Realisierung bei erhaltener Skala,
- **skalendominiert:** $S_M$ weit von $S_0$, $F_M$ klein → ähnliche Winkelstruktur, aber andere Größe,
- **kombiniert:** $K_M$ als einheitlicher **Krümmungssensor** auf der Triangulierung.

---

## Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| Definitionen $F_M$, $S_M$, $K_M$ | **Definition** |
| Ebenenformel für $\operatorname{Area}(\mathrm{Mor}(\Delta))$ | **Theorem** (klassisch) |
| $K_M \sim c_1 K_G + c_2 K_G^2 + \cdots$ asymptotisch | **Conjecture** (offen) |
| $\mu_M$ als diskretes Konformitätsfeld | **heuristisch** |
| $T_M$-Iteration, EABC-Operator $(A,B,C)\mapsto(E_A,E_B,E_C)$ | **heuristisch / spekulativ** |
| Brücke $T_M \to$ Collatz / κ-Invarianz | **spekulativ** (keine Beweisansprüche) |
| Morley als Krümmungs-/Konformitätssensor (nicht Metrik) | **Leitsatz** (organisatorisch) |

---

## Verhältnis zu Stufe 3 (κ-Invarianz)

Stufe 3 (`collatz_stufe3_kappa_invarianz.md`) fokussiert **Kodierungsinvarianz** von $R(k)$, $h_F$, $\mathcal{L}_{\mathrm{arith}}^*$ — nicht Geometrie als Hauptangriff.

Dieses Dokument hält die **geometrische Brücke** offen, aber **epistemisch abgegrenzt**:

| Stufe 3 (aktiv) | Morley-Erweiterung (Brücke) |
|-----------------|----------------------------|
| κ-Familien, Äquivalenzklassen | $T_M$-Iteration auf DC-Dreiecken |
| $R(k)$, $h_F$ | $F_M$, $S_M$, $K_M$ auf Triangulierungen |
| $\mathcal{L}_{\mathrm{arith}}^*$ | $\mu_M$-Felder, diskrete Konformität |

**Keine Prioritätskonkurrenz:** Morley-Sensorik ist **Rang 6+** in der Forschungshierarchie (nach `collatz_generalangriff_2026.md`).

---

## Living spreadsheet (Morley-Erweiterung)

| Zeile | Objekt | Ebene | Stand Juni 2026 |
|-------|--------|-------|-----------------|
| 1 | $F_M$ (Morley-Form) | Definition | **definiert** |
| 2 | $S_M$ (Morley-Skala) | Definition | **definiert** |
| 3 | $K_M = \alpha F_M + \beta(S_M-S_0)^2$ | Definition | **definiert** |
| 4 | $K_M \sim K_G$ asymptotisch | Conjecture | **offen** |
| 5 | $\mu_M$ Konformitätsfeld | heuristisch | **skizziert** |
| 6 | $T_M$ als EABC-Operator | heuristisch | **skizziert** |
| 7 | Sensor statt Metrik (Leitsatz) | organisatorisch | **revidiert** |
| 8 | Numerik / Lean-Formalismus | Experiment | **fehlt** |

---

## Konkrete nächste Schritte

1. **Numerik:** $F_M$, $S_M$, $K_M$ für DC-Dreieckskonfigurationen aus `collatz_dc_morley_walter.pdf` berechnen.
2. **Asymptotik:** kleine Dreiecke auf Flächen konstanter $K_G$ (Sphäre, Pseudosphäre) — Test der Conjecture $K_M \sim c_1 K_G$.
3. **$\mu_M$-Feld:** Triangulierung der mod-12-DC-Geometrie; Vergleich mit Beltrami-Norm.
4. **$T_M$-Iteration:** Attraktorverhalten für typische DC-Dreiecke (Experiment, kein Theorem).
5. **Lean:** optional `MorleyForm`, `MorleyScale` als reine Definitionen — **kein** Collatz-`sorry`-Abbau.

**Bewusst nicht:** Behauptung einer neuen Riemann-Metrik oder Collatz-Beweis via Morley.

---

## Artefakte und Verweise

| Datei | Rolle |
|-------|-------|
| `collatz_morley_metrik_erweiterung.md` | **kanonisches** Morley-Erweiterungsdokument (dieses) |
| `MorleyWalter.tex` / `collatz_dc_morley_walter.pdf` | Bamberg-Modell, Morley–Walter-Deutung |
| `collatz_schlussartikel_arxiv.tex` | DC-Dreiecksinterpretation (illustrativ) |
| `collatz_stufe3_kappa_invarianz.md` | Stufe 3 — κ-Invarianz (Hauptangriff) |
| `collatz_generalangriff_2026.md` | Forschungsreport, Geometrie abgegrenzt |
| `collatz_offene_punkte.md` | Synthese, Prioritäten |
