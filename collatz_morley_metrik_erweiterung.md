# **Morley-Erweiterung: Operator, Sensor und Geometrie zweiter Ordnung**

**Branch:** `collatz/morley-sensorik-docs` · **PR:** #43 (stacked auf #42)  
**Stand:** Juni 2026  
**Kein Collatz-Beweis.** Geometrische Brücke — epistemisch abgegrenzt von Stufe-3-Kern (κ-Invarianz).

**Methodik:** Tao-Stil — Definition / Zeuge / Experiment / Theorem / Conjecture strikt trennen (`collatz_formalisierung_tao_stil.md`).

**Vorgänger:** `MorleyWalter.tex`, `collatz_dc_morley_walter.pdf`, `collatz_schlussartikel_arxiv.tex` (§ Morley–Walter-Geometrie, illustrativ).

---

## Kernthese (revidiert)

> **Boxed (Leitsatz):**  
> *Der eigentliche neue Akteur ist nicht eine Morley-Metrik, sondern der Morley-Operator $T_M$.*

Morley liefert **keine** neue Riemann-Metrik $g_{ij}$ auf einer glatten Fläche. Der Morley-Kern lebt auf **Dreiecken** (Winkel, Flächenverhältnis), nicht im Tangentialraum — näher an **Regge-Kalkül**, **diskreter Gauß-Krümmung**, **Voronoi/Delaunay-Triangulierungen**, **FEM-Netzen** und **konformen Triangulierungen** als an einer infinitesimalen Metrik.

Sobald Morley als **Operator** $T_M$ gelesen wird (nicht als statisches Objekt), öffnen sich natürliche Fragen nach **Dynamik**, **Invarianten**, **Flüssen** und der **EABC-Brücke**. Die Sensoren $F_M$, $S_M$, $K_M$ und $\mu_M$ messen dann **Zustände** entlang von $T_M$-Trajektorien — sie ersetzen den Operator nicht.

> **Boxed (Nebenleitsatz):**  
> *Morley nicht als Metrik, sondern als diskreter Krümmungs- und Konformitätssensor auf Dreieckszuständen.*

---

## Kontext: Morley vs. Riemann

| Ebene | Objekt | Lokalität | Typische Literatur |
|-------|--------|-----------|-------------------|
| **Riemann** | Metrik $g_{ij}$, Christoffel, Krümmungstensor | infinitesimal (Tangentenraum) | glatte Differentialgeometrie |
| **Morley** | Winkeldrittel, gleichseitiges Morley-Dreieck $\mathrm{Mor}(\Delta)$ | **winkel-lokal** auf einem Dreieck | klassische euklidische Geometrie |
| **Diskret** | Defizitwinkel, Regge-Winkel, dualer Volumenfluss | **triangulations-lokal** | Regge (1961), discrete exterior calculus |

**Fundamentaler Unterschied (erste vs. zweite Ordnung):**

| Paradigma | Minimales Objekt | Lokalität | Analogie (statistische Physik) |
|-----------|------------------|-----------|--------------------------------|
| **Riemann** | $(p, v)$ — Punkt + Tangente | infinitesimale Umgebung | **Teilchen** |
| **Morley** | $(A,B,C)$ — Dreieck mit Fläche | winkel- und flächenlokal | **Zelle** |

Riemann-Geometrie baut auf **erster Ordnung** (Tangentenvektor, infinitesimale Nachbarschaft). Morley-Geometrie baut auf **zweiter Ordnung** (Dreieck als minimales Objekt mit nichtverschwindender Fläche). Das ist kein bloßer Diskretisierungsschritt, sondern ein **Wechsel der Grundobjekte**.

**Konsequenz:** Eine naive „Morley-Metrik“ $K_M$ als einzelne skalare Größe war zu grob. Stattdessen trennen wir zwei **orthogonal lesbare** Sensoren:

1. **Morley-Form** $F_M$ — Gleichseitigkeitsabweichung (rein winkelbasiert),
2. **Morley-Skala** $S_M$ — Größenänderung (Flächenverhältnis).

Die kombinierte **Morley-Krümmung** $K_M$ gewichtet beide mit Parametern $\alpha,\beta$ und einem euklidischen Referenzwert $S_0$.

---

## Was bewiesen ist — und was nicht übertragen wird

> **Boxed (Warnung):**  
> *Den euklidischen Morley-Satz ohne Definition auf gekrümmte Flächen zu übertragen, ist epistemisch unzulässig.*

### Bewiesen (nur euklidisch)

In der **Ebene** gilt (Morley-Satz, **klassisch / bewiesen**):

$$T_M(\Delta) = \mathrm{Mor}(\Delta), \qquad \mathrm{Mor}(\Delta)\ \text{ist gleichseitig.}$$

$T_M$ wirkt dort als **Symmetrisierungsoperator**: aus beliebigem Dreieck wird in einem Schritt ein gleichseitiger Morley-Kern. Dieses Resultat ist **strikt euklidisch** — es folgt aus geraden Linien, Winkeldrittelung und dem klassischen Morley-Beweis.

**Epistemisches Label:** **klassisch / bewiesen** (Morley).

### Nicht ohne Definition übertragbar

Auf einer Riemann-Fläche $(M,g)$ ist „Winkeldrittelung“ und „Schnittpunkt der Drittelnden“ **nicht** eindeutig, solange man keine Realisierungsvorschrift wählt. Der Schritt „Morley auf $M$“ ist daher **kein** Theorem, sondern ein **Definitionsproblem**.

**Epistemisches Label:** **Definition / offen** (Variante muss gewählt werden).

---

## Definitionen

Sei $\Delta = ABC$ ein nicht-degeneriertes Dreieck in der euklidischen Ebene (oder auf einer konformen Karte einer Riemann-Fläche). Die **Morley-Winkel** $\theta_i^M$ ($i=1,2,3$) seien die Innenwinkel des Morley-Dreiecks $\mathrm{Mor}(\Delta)$ — das gleichseitige Dreieck aus den Schnittpunkten benachbarter Winkeldrittelnden (Morley-Satz).

### Morley-Form (Gleichseitigkeitsabweichung)

$$F_M(\Delta) \;:=\; \sum_{i=1}^{3} \left(\theta_i^M - \frac{\pi}{3}\right)^2.$$

- $F_M = 0$ genau dann, wenn $\mathrm{Mor}(\Delta)$ perfekt gleichseitig ist (in der Ebene **immer** — Morley-Satz).
- Auf einer **krummen** Fläche (nach Wahl von $T_M^{(g)}$, s. unten) misst $F_M$ den **Formdefekt** des geodätischen Morley-Kerns relativ zur Referenz „gleichseitig“.

**Epistemisches Label:** **Definition** (definitorisch); $F_M=0$ in $\mathbb{R}^2$ **klassisch / bewiesen**.

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

### Erste konkrete Vermutung (publikationsreif formuliert)

**Definition (Arbeitswahl):** $T_M^{(g)}$ sei der **geodätische Morley-Operator** auf einer triangulierten Riemann-Fläche — d. h. eine explizit gewählte Realisierung von Winkeldrittelung und Morley-Konstruktion auf geodätischen Dreiecken (s. § Geometrie zweiter Ordnung).

> **Conjecture (Morley–Gauß-Form, erste Ordnung):**  
> Für hinreichend kleine geodätische Dreiecke $\Delta$ mit Schwerpunkt $p$ und Fläche $A = \operatorname{Area}(\Delta)$ existiert eine Konstante $c$ (abhängig von der Normalisierung von $T_M^{(g)}$, aber **unabhängig** von der Dreiecksorientierung im Grenzfall), sodass
> $$F_M(\Delta) = c\, K_G(p)\, A + O(A^2).$$

**Lesart:** Wenn wahr, ist der Morley-Operator ein **diskreter Krümmungssensor** — die Gleichseitigkeitsabweichung $F_M$ kodiert in erster Ordnung die Gauß-Krümmung, nicht eine neue Riemann-Metrik.

**Folgefragen (offen):** Analoge Entwicklung für $S_M - S_0$; kombiniertes $K_M = \alpha F_M + \beta(S_M-S_0)^2$ als zweite Lesart.

**Epistemisches Label:** **Conjecture**.

### Ältere asymptotische Frage (Kombinationssensor)

> **Conjecture (Morley–Gauß-Asymptotik, kombiniert):**  
> Existieren Konstanten $c_1, c_2, \ldots$ (abhängig von der Normalisierung von $K_M$), sodass für kleine Dreiecke
> $$K_M(\Delta) = c_1\, K_G + c_2\, K_G^2 + O(A)$$
> oder mindestens $K_M(\Delta) \sim c_1 K_G$ im Grenzfall $A\to 0$?

**Status:** **offen** — die **primäre** Vermutung für eine Arbeit ist $F_M \sim c\,K_G A$ (oben); $K_M$ ist ein abgeleiteter Kombinationssensor.

**Epistemisches Label:** **Conjecture** (sekundär).

---

## Geometrie zweiter Ordnung: der Morley-Operator $T_M$

Dieser Abschnitt hebt den **Operator** über die **Sensoren** — und damit über jede statische Metrik-Lesart. Der **erste Schritt** auf gekrümmten Flächen ist nicht „Morley-Fluss“, sondern die **Definition** des geodätischen Morley-Operators.

### Formraum $\mathcal{S}$ — primäre Größe ist die Dreiecksform

Ein Dreieck hat **zwei Freiheitsgrade für die Form** (modulo Skala und Rotation). Der **Formraum** sei

$$\mathcal{S} \;:=\; \{\text{Dreiecksformen}\} / \sim_{\text{Ähnlichkeit}}.$$

Der euklidische Morley-Operator wirkt zuerst auf **Form**, nicht auf Kantenlängen:

$$T_M : \mathcal{S} \longrightarrow \mathcal{S}, \qquad [\Delta] \longmapsto [\mathrm{Mor}(\Delta)].$$

Der **gleichseitige Punkt** in $\mathcal{S}$ ist ausgezeichnet (Morley-Attraktor in der Ebene). Analogien: Renormierungsoperatoren, Newton-Iteration, diskrete Ricci-Flüsse — allesamt **Abbildungen auf einem Zustandsraum**, nicht punktuelle Metriken.

**Epistemisches Label:** **Definition** (Formraum); Dynamik auf $\mathcal{S}$ **offen**.

### Dreieckraum und geodätischer Operator $T_M^{(g)}$

Sei $(M,g)$ eine orientierte Riemann-Fläche (oder triangulierte Fläche mit geodätischen bzw. Regge-Dreiecken). Der **Dreieckraum** sei

$$\mathcal{T}(M) \;:=\; \{\text{nicht-degenerierte geodätische Dreiecke auf } M\}$$

(modulo Kongruenz, falls explizit quotientiert). Der **geodätische Morley-Operator** ist die zu wählende Abbildung

$$T_M^{(g)} : \mathcal{T}(M) \longrightarrow \mathcal{T}(M), \qquad \Delta \longmapsto \mathrm{Mor}^{(g)}(\Delta).$$

In der Ebene schreiben wir $T_M^{(g)} = T_M$ und $T_{\mathbb{R}^2}(\Delta) = \mathrm{Mor}(\Delta)$ (**klassisch / bewiesen**).

### Vier Realisierungsvarianten (Definitionsproblem)

Auf $M \neq \mathbb{R}^2$ gibt es **mehrere** sinnvolle Weisen, „Winkeldrittelung“ und Morley-Schnittpunkte zu definieren. Diese müssen **nicht** dasselbe Objekt liefern:

| Variante | Idee | Typische Werkzeuge |
|----------|------|-------------------|
| **Geodätische Winkel** | Innenwinkel aus geodätischen Kanten; Drittelung im Tangentialraum, dann Rücktransport | Winkelsumme, Gauß–Bonnet lokal |
| **Lokale Karten** | Morley-Konstruktion in konformer/isothermer Karte, Rückabbildung auf $M$ | $z$-Koordinate, $\lambda\,|dz|^2$ |
| **Paralleltransport** | Winkeltrisektionsrichtungen parallel zu den Kanten transportieren | Levi-Civita-Verbindung |
| **Exponentialabbildung** | Konstruktion im tangentialeuklidischen Modell $\exp_p(T_pM)$, dann Einbettung | Normal coordinates |

> **Boxed (Definitionspriorität):**  
> *Bevor Morley-Dynamik, Fixpunkte oder $F_M \sim K_G A$ untersucht werden, muss $T_M^{(g)}$ explizit gewählt und die Variantenabweichung dokumentiert werden.*

**Epistemisches Label:** **Definition / offen** (Variante muss gewählt werden).

### Morley-Dynamik (nach Definition von $T_M^{(g)}$)

Die **diskrete Morley-Dynamik** ist

$$\Delta_{n+1} \;=\; T_M^{(g)}(\Delta_n), \qquad n = 0,1,2,\ldots$$

Auf dem Formraum $\mathcal{S}$ (Quotient nach Skala/Rotation) induziert $T_M^{(g)}$ ebenfalls eine Dynamik — die **primäre** Lesart für Invarianten und Attraktoren.

Der **Morley-Fluss** (heuristische Benennung) ist die Trajektorie $(\Delta_n)_{n\geq 0}$ — ein **diskretes Analogon** zu Ricci-Fluss, Newton-Iteration oder Renormierung, aber auf $\mathcal{T}(M)$ bzw. $\mathcal{S}$ statt auf dem Raum der Metriken:

| Kontinuierlich / klassisch | Morley-Dynamik (diskret) |
|----------------------------|--------------------------|
| $\partial_t g = -2\,\mathrm{Ric}(g)$ (Ricci-Fluss) | $\Delta_{n+1} = T_M^{(g)}(\Delta_n)$ |
| Newton-Iteration | $T_M^k$ als Iteration auf $\mathcal{S}$ |
| Renormierungsoperator | $T_M^{(g)}$ als Form-Symmetrisierung |
| evolviert Metrik | evolviert **Dreiecksform** |

Sobald $T_M^{(g)}$ definiert und iteriert wird, entstehen **automatisch** Fragen der dynamischen Systemtheorie:

- **Fixpunkte** $\Delta^*$ mit $T_M^{(g)}(\Delta^*) = \Delta^*$
- **Periodische Orbits** $T_M^{(g)\,k}(\Delta) = \Delta$
- **Attraktoren** und Basins auf $\mathcal{S}$
- **Invarianten** unter $T_M^{(g)}$ bzw. $T_M^{(g)\,k}$
- **Entropie** und Komplexität der Orbitstruktur
- **Kontraktivität** in einer geeigneten Metrik auf $\mathcal{S}$

**Epistemisches Label:** **offen** — keine Antworten in diesem Dokument.

### Morley-Satz als Symmetrisierungsoperator (nur $\mathbb{R}^2$)

Im euklidischen Fall (Morley-Satz, **klassisch / bewiesen**): Für jedes $\Delta \subset \mathbb{R}^2$ ist $\mathrm{Mor}(\Delta)$ **gleichseitig**. Damit ist $T_{\mathbb{R}^2}$ auf der Äquivalenzklasse gleichseitiger Dreiecke **konstant** — Morley wirkt in der Ebene als **Symmetrisierungsoperator**: aus beliebiger Form wird (in einem Schritt) maximale Winkelsymmetrie im Kern.

**Lesart:** In der Ebene ist der Attraktor in $\mathcal{S}$ **trivial und universal** — der gleichseitige Punkt. Auf gekrümmten Flächen (nach Wahl von $T_M^{(g)}$) bricht diese Universalität; genau dort wird $F_M$ zum **Krümmungssensor** (Conjecture oben).

### Krümmungssonde: $T_M^{(g)}(\Delta)$ vs. euklidische Referenz

Für ein kleines Dreieck $\Delta$ auf $(M,g)$ sei $T_{\mathrm{euclid}}(\Delta)$ die Morley-Realisierung, wenn $\Delta$ **als euklidisch** (gleiche Kantenlängen/Winkel in der Tangentialebene) konstruiert wird. Die **Abweichung**

$$\delta_M(\Delta) \;:=\; T_M^{(g)}(\Delta) - T_{\mathrm{euclid}}(\Delta)$$

(mit geeigneter Norm auf Winkel- und Skalenabweichung, z. B. via $F_M$, $S_M$) misst **lokale Geometrie**: der Morley-Kern als **Krümmungssonde**. Die Sensoren $F_M$, $S_M$, $K_M$ aus § Definitionen sind **Koordinaten** dieser Abweichung, nicht der Operator selbst.

**Epistemisches Label:** **heuristisch** — präzise Norm und Grenzfall $A\to 0$ **offen**; Verbindung zur Conjecture $F_M = c\,K_G A + O(A^2)$ (oben).

### Komplexe Geometrie: Kreuzverhältnis und konforme Invarianten

Auf $\mathbb{C}$ identifiziere $\Delta = (z_1,z_2,z_3)$ mit $z_i \in \mathbb{C}$. Ein nicht-degeneriertes Dreieck trägt einen **Kreuzverhältnis-Typ** (Möbius-Äquivalenzklasse). Die Morley-Abbildung

$$\mathcal{M} : (z_1,z_2,z_3) \longmapsto (m_1,m_2,m_3)$$

ist eine **punktweise diskrete Transformation** auf Dreieckskonfigurationen — kein globales $f:\mathbb{C}\to\mathbb{C}$, aber eine **abbildungsdefinierte** Evolution von Kreuzverhältnissen. Natürliche Fragen:

- Welche **konformen Invarianten** (Kreuzverhältnis, Schwarzian) sind unter $\mathcal{M}$ approximativ erhalten?
- Liegt $\mathcal{M}$ näher an **klassischer Riemannscher Flächentheorie** (konforme Struktur) als an einer neuen Metrik?

Verbindungen (Literaturrahmen, **nicht bewiesen**):

- Schwarzsche Ableitung (dritte invariante Ableitung unter Möbius),
- Beltrami-Differential $\mu = \partial_{\bar z}f / \partial_z f$,
- Teichmüller-Theorie und quasikonforme Abbildungen,
- diskrete Riemann-Flächen (circle packing, discrete conformal weldings).

Das Kantenquotientenfeld $\mu_M$ (§ Komplexe Ebene unten) ist die **lokale Ablese** von $\mathcal{M}$ entlang Kanten.

**Epistemisches Label:** **heuristisch**.

### EABC revidiert: $E = T_M$, nicht ein Punkt

**Alt (zu naiv):** Morley-Kern als zusätzlicher Punkt $E$ in der statischen EABC-Viererkonfiguration $(A,B,C,E)$.

**Neu (operatorisch):**

$$(A,B,C) \xrightarrow{\,T_M\,} (E_A, E_B, E_C),$$

wobei $(E_A,E_B,E_C)$ die Eckpunkte von $\mathrm{Mor}(\Delta)$ sind. **$E$ ist kein Punkt, sondern die Wirkung des Operators** — eine **Ordnungs-/Ausgleichsstruktur** (Winkeldrittelung, Symmetrisierung). In EABC-Terminologie: das vierte Element ist **$T_M$ selbst** (bzw. seine Ausgabe auf der aktuellen Dreieckskonfiguration), nicht eine feste Spitze im mod-12-Diagramm.

| Lesart | Objekt | Status |
|--------|--------|--------|
| statisch | Punkt $E$ | **verworfen** (zu naiv) |
| operatorisch | $T_M : (A,B,C)\mapsto(E_A,E_B,E_C)$ | **Leitrahmen** |
| iteriert | $\Delta_{k+1}=T_M(\Delta_k)$ | **Morley-Fluss** |

**Brücke zu Collatz / κ (spekulativ):** DC-Blöcke als Dreieckskonfigurationen; $T_M$ als **geometrischer Renormierungsoperator** parallel zu κ-Familien — **keine** bewiesene Implikation für Collatz-Trajektorien.

**Epistemisches Label:** **heuristisch / spekulativ**.

### Invarianten des Morley-Operators (gesichert vs. offen)

> **Boxed (Invarianten des Morley-Operators):**  
> **Gesichert** (ein Schritt $T_M$, euklidisch):
>
> 1. $F_M(\Delta)=0$ — Morley-Kern gleichseitig (**Theorem**, Morley-Satz).
> 2. $S_M(\Delta)$ invariant unter Ähnlichkeit (**Definition** / klassische Rechnung).
> 3. Winkeldrittel-Struktur — **3-fache** Phasensymmetrie pro Anwendung (**Definition**).
>
> **Offen** (Iteration $T_M^k$, gekrümmte Flächen, EABC):
>
> 4. Allgemeine Fixpunkte / Attraktoren auf $\mathcal{T}(M)$ (**offen**).
> 5. Funktorialität von $T_M$ auf mod-12-EABC-Labels (**spekulativ**).
> 6. Dualitätsverträglichkeit mit Ikosaeder–Dodekaeder-Renormierung (**heuristisch**, s. § Ikosaeder–Dodekaeder-Dualität).

### Hauptforschungsfrage und Teilfragen

> **Boxed (Hauptfrage):**  
> *Welche Invarianten besitzt der Morley-Operator auf triangulierten Riemann-Flächen?*

| Nr. | Teilfrage | Status |
|-----|-----------|--------|
| 1 | **Fixpunkt-Dreiecke:** $\exists\,\Delta^*$ mit $T_M(\Delta^*)=\Delta^*$? | **offen** |
| 2 | **Periodische Morley-Zyklen:** $T_M^k(\Delta)=\Delta$ für $k>1$? | **offen** |
| 3 | **Kontraktivität:** Ist $T_M$ kontraktiv in natürlicher Metrik auf $\mathcal{T}(M)$? | **offen** |
| 4 | **Krümmungsrekonstruktion:** Lässt sich $K_G$ aus $\delta_M(\Delta)$ bzw. $F_M,S_M$ rekonstruieren? | **Conjecture** (s. asymptotische Frage) |
| 5 | **Diskretes konformes Feld:** Entsteht unter $T_M$-Iteration ein stabiles $\mu_M$-Feld? | **heuristisch / offen** |
| 6 | **Konstante Krümmung:** Verhalten auf Sphäre ($K_G>0$), Ebene ($K_G=0$), hyperbolisch ($K_G<0$)? | **offen** (Ebene: Morley-Satz = Symmetrisierung) |

**Epistemisches Label der Gesamtfrage:** **offen** — Forschungsprogramm in Richtung diskrete Differentialgeometrie, **nicht** Collatz-Beweis.

---

## Ikosaeder–Dodekaeder-Dualität und Triangulierungsdualität

Verbindung des $T_M$-Rahmens mit der **platonischen Dualität** $I \leftrightarrow D$ und den EABC-Polyeder-Experimenten im Repo.

### Platonic dual pair (klassisch)

| Körper | Ecken | Kanten | Flächen | Eckfigur |
|--------|-------|--------|---------|----------|
| **Ikosaeder** $I$ | 12 | 30 | 20 (Dreiecke) | 5 Dreiecke/Ecke |
| **Dodekaeder** $D$ | 20 | 30 | 12 (Pentagone) | 3 Pentagone/Ecke |

**Dualitätsregel:** Ecken $\leftrightarrow$ Flächen tauschen die Rollen; Kanten bleiben dual verknüpft. Auf $S^2$ ist die ikosaedrische Triangulierung (20 sphärische Dreiecke) **dual** zur dodekaedrischen Voronoi-Zellstruktur (12 Pentagon-Zellen).

**Epistemisches Label:** **Theorem** (klassische Polyedergeometrie).

### Repo-Artefakte

| Artefakt | Inhalt |
|----------|--------|
| `eabc_icosahedron_test.py` | `build_dodecahedron_dual()`, `push_assignment_to_dual()`, `run_duality_test()`, `run_core_lemma_tests()` (Test 3: $\Delta(M)$-Stabilität) |
| `collatz_ikosaeder_spannung.pdf` | Ikosaeder-Spannung im EABC-Kontext |
| `collatz_kepler_gedankenexperiment.tex` | Zeugen-Tabelle: Anisotropie $=0$ (geometrisch) vs. $>0$ (Zufall) |
| `MorleyWalter.tex` | Kepler-Hülle: dodekaedrische Schließung **nach** Morley-Kern und Walter-Schicht |
| `kapitel-eabc-phaenomenologie.tex` | H4/120-Zeller, fünffache Ikosaeder-Symmetrie |
| `Zwilling VI.tex` | Mechanische Frustration Ikosaeder (12 Ecken, 20 Flächen) vs. Dual (20 Ecken, 12 Flächen) |

### Numerischer Zeuge: $\Delta(M)$ unter Dualität

`eabc_icosahedron_test.py` konstruiert das Dodekaeder-Dual via Flächenzentren $\to$ Ecken und Ikosaeder-Ecken $\to$ Pentagon-Flächen. **Test 3** in `run_core_lemma_tests`:

- $\Delta(M)_{\text{Ikosaeder}} = \Delta(M)_{\text{Dodekaeder-Dual}}$ am geometrischen Fixpunkt (**Experiment**),
- $D$- und $H$-Normen sind **nicht** dual-invariant — nur der Isotropie-Defekt $\Delta(M)$ bleibt stabil.

**Lesart:** Isotropie ist dualitätsverträglich; Holonomie-/Richtungssensoren ($D$, $H$) brechen die Vertex↔Face-Symmetrie.

**Epistemisches Label:** **Experiment** (reproduzierbar, kein Collatz-Befund).

### Triangulierungsdualität und $T_M$

$$\text{Delaunay-Triangulierung} \;\longleftrightarrow\; \text{Voronoi-Dualität}.$$

| Struktur | Rolle im Morley-Rahmen |
|----------|------------------------|
| Ikosaeder-Triangulation | Träger von $\Delta \in \mathcal{T}(S^2)$; $T_M$ wirkt **pro Dreieck** |
| Dodekaeder-Dual | Voronoi-Zellen; Morley **nicht** direkt auf Pentagone definiert |
| Fächer-Triangulation | `triangulate_faces()` — macht Dual Morley-zugänglich ($P \to 3$ Dreiecke) |

**Symmetrie-Ordnungen (heuristisch):**

| Operator / Struktur | Ordnung | Lesart |
|---------------------|---------|--------|
| $T_M$ (Winkeldrittel) | **3** | gleichseitiger Morley-Kern |
| Ikosaeder-Eckfigur | **5** | $A_5$-Symmetrie auf $S^2$ |
| Dodekaeder-Flächen | **5** | Kepler-Hülle (`MorleyWalter.tex`) |

Morley (3-fach) und Ikosaeder (5-fach) sind **verschiedene** Symmetrieoperatoren — keine Identität, aber komplementäre Hierarchie: Zelle $\to$ Polyeder-Hülle.

**Epistemisches Label:** Triangulierungsdualität **Theorem**; Morley↔Ikosaeder-Kopplung **heuristisch**.

### EABC und mod 12 — ehrliche Einordnung

In `eabc_icosahedron_test.py`:

- $E\equiv 1$, $A\equiv 5$, $B\equiv 7$, $C\equiv 11 \pmod{12}$,
- Ikosaeder-Koordinaten-Orbits $\{2,4,8,10\} \bmod 12$ (`ICO_RAW_RESIDUES`) $\to$ EABC via `nearest_eabc_residue`,
- geometrische Zuweisung: **3+3+3+3**-Balance (12 Ecken).

| Zahl | Ikosaeder | EABC / mod 12 | Bewertung |
|------|-----------|---------------|-----------|
| 12 | Ecken | $\mathbb{Z}/12\mathbb{Z}$ | **Zahlähnlichkeit**, keine Isomorphie |
| 20 | Dreiecksflächen | — | keine direkte Viererstruktur |
| 30 | Kanten | geteilt mit Dual | Dualitäts-Invariante |
| 5 | Eckfigur | $A\equiv 5$ | **Notation**, kein Theorem |

**Gestützte Aussage (schwächer):** geometrische EABC-Zuweisung minimiert Anisotropie; Zufallszuweisung nicht (`collatz_offene_punkte.md`).

**Epistemisches Label:** Anisotropie-Zeuge **Experiment**; 12↔mod-12 **spekulativ**.

### Operatorische Parallelität: $T_M$ und $R_{\mathrm{EABC}}$

| Operator | Domäne | Bild |
|----------|--------|------|
| $T_M$ | $\mathcal{T}(M)$ | $\mathrm{Mor}(\Delta)$ — **3-er** Symmetrisierung |
| $R_{\mathrm{EABC}}$ | Labelings auf Polyedergraph | geometrisch zulässige 3+3+3+3-Zuweisung (`apply_r_eabc`) |

> **Conjecture (dualitätsverträgliche Renormierung, offen):**  
> Existiert eine Abbildung $\Phi$, die $T_M$-Trajektorien auf ikosaedrischen Teil-Triangulierungen mit $R_{\mathrm{EABC}}$ zu $\Delta(M)=0$-Zuständen koppelt — analog zur $\Delta(M)$-Stabilität unter $I \leftrightarrow D$?

**Status:** **offen**.

### Sphärische Morley-Zellen (Forschungsfrage)

Der klassische Morley-Satz gilt auf $S^2$ **nicht** unverändert. Dennoch:

- **Experiment (offen):** $F_M$, $S_M$ auf ikosaedrischen Sphärendreiecken,
- **heuristisch:** Hierarchie Morley $\to$ Walter $\to$ Kepler/Dodekaeder (`MorleyWalter.tex`).

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

## EABC-Perspektive (Kurzverweis)

Die operatorische EABC-Lesart — $(A,B,C)\xrightarrow{T_M}(E_A,E_B,E_C)$, $E$ als Wirkung nicht als Punkt, Morley-Fluss $\Delta_{n+1}=T_M(\Delta_n)$ — ist in § **Geometrie zweiter Ordnung** ausgeführt; Ikosaeder–Dodekaeder-Dualität und $R_{\mathrm{EABC}}$ in § **Ikosaeder–Dodekaeder-Dualität**. Brücke zu Collatz/κ bleibt **spekulativ**; vgl. `collatz_dc_morley_walter.pdf`, `eabc_icosahedron_test.py`.

---

## Strukturbrücke $(M, g, T_M, \mathcal{T}(M))$

Die folgende Kette ordnet die Ebenen (kein Theorem, **organisatorisches Schema**):

$$\boxed{
\text{Riemann-Fläche } (M,g)
\;\to\;
\mathcal{T}(M) \;\text{(ikosaedrisch / Delaunay)}
\;\to\;
T_M : \mathcal{T}(M) \to \mathcal{T}(M)
\;\to\;
\text{Polyeder-Dualität } I \leftrightarrow D
\;\to\;
\text{Morley-Fluss } (\Delta_n)
\;\to\;
\text{Sensoren } F_M, S_M, K_M, \mu_M
}$$

| Komponente | Rolle |
|------------|-------|
| $(M,g)$ | glatte Hintergrundgeometrie (falls vorhanden) |
| $\mathcal{T}(M)$ | Dreieckraum — **Grundobjekt zweiter Ordnung** |
| $T_M$ | Morley-Operator; Symmetrisierung (Ebene) / Krümmungssonde (krumm) |
| $(\Delta_n)$ | Morley-Fluss; Fixpunkte, Zyklen, Attraktoren **offen** |
| $F_M, S_M, K_M, \mu_M$ | Sensoren auf Zuständen — **nicht** der Operator |

**Nähe zu:** diskreter Differentialgeometrie, Regge-Kalkül, diskreter konformer Geometrie, Ricci-Fluss (als kontinuierliches Analogon) — **nicht** zu einer neuen „Morley-Metrik“ im Riemannschen Sinn.

---

## Zwei Sensoren vs. naive $K_M$

| Größe | Misst | Sensitivität |
|-------|-------|--------------|
| $F_M$ | Winkeldefekt des Morley-Kerns | rein **angular** |
| $S_M$ | Flächenverhältnis Hülle/Kern | **metrisch** (Skala) |
| $K_M$ | gewichtete Kombination | tunable via $\alpha,\beta,S_0$ |

Die Trennung verhindert die Verwechslung „Morley = neue Metrik“. Der **Operator** $T_M$ ist der primäre Akteur; $K_M$ ist ein **abgeleiteter Sensor**. In der Praxis:

- **Formdominiert:** $F_M$ groß, $S_M \approx S_0$ → stark nichtgleichseitige Morley-Realisierung bei erhaltener Skala,
- **skalendominiert:** $S_M$ weit von $S_0$, $F_M$ klein → ähnliche Winkelstruktur, aber andere Größe,
- **kombiniert:** $K_M$ als einheitlicher **Krümmungssensor** auf der Triangulierung.

---

## Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| Definitionen $F_M$, $S_M$, $K_M$ | **Definition** |
| Ebenenformel für $\operatorname{Area}(\mathrm{Mor}(\Delta))$ | **Theorem** (klassisch) |
| Morley-Operator $T_M$ auf $\mathcal{T}(M)$ | **Definition** |
| Morley-Fluss $\Delta_{n+1}=T_M(\Delta_n)$ | **Definition** (Dynamik) |
| Fixpunkte, Zyklen, Kontraktivität, universelle Konvergenz | **offen** |
| Morley-Satz: $T_{\mathbb{R}^2}$ symmetrisiert zu gleichseitig | **Theorem** (klassisch) |
| $\delta_M = T_M(\Delta)-T_{\mathrm{euclid}}(\Delta)$ als Krümmungssonde | **heuristisch** |
| $\mathcal{M}:(z_1,z_2,z_3)\mapsto(m_1,m_2,m_3)$; Schwarz/Beltrami/Teichmüller | **heuristisch** |
| $E = T_M$ (EABC operatorisch, nicht Punkt $E$) | **heuristisch / spekulativ** |
| Hauptfrage: Invarianten von $T_M$ auf triangulierten Flächen | **offen** (Forschungsprogramm) |
| $K_M \sim c_1 K_G + c_2 K_G^2 + \cdots$ asymptotisch | **Conjecture** (offen) |
| $\mu_M$ als diskretes Konformitätsfeld | **heuristisch** |
| Brücke $T_M \to$ Collatz / κ-Invarianz | **spekulativ** (keine Beweisansprüche) |
| Ikosaeder–Dodekaeder-Dualität; $\Delta(M)$ dual-stabil | **Theorem** (Dualität) / **Experiment** ($\Delta(M)$) |
| 12 Ikosaeder-Ecken $\leftrightarrow$ mod 12 | **spekulativ** |
| $T_M$ vs. $R_{\mathrm{EABC}}$ Renormierung | **heuristisch / offen** |
| Operator statt Metrik (Leitsatz) | **organisatorisch** |
| Sensoren $F_M,S_M,K_M$ (Nebenleitsatz) | **organisatorisch** |

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
| 6 | $T_M$, Morley-Fluss, Invariantenfragen | offen / heuristisch | **§ Geometrie zweiter Ordnung** |
| 7 | $E=T_M$ (EABC operatorisch) | heuristisch | **revidiert** |
| 8 | Operator statt Metrik (Leitsatz) | organisatorisch | **revidiert** |
| 9 | Sensoren $F_M,S_M,K_M$ | Definition | **definiert** |
| 10 | Ikosaeder–Dodekaeder + $T_M$ | heuristisch / Experiment | **§ Ikosaeder–Dodekaeder-Dualität** |
| 11 | Numerik / Lean-Formalismus | Experiment | **fehlt** |

---

## Konkrete nächste Schritte

1. **$T_M$-Dynamik:** Fixpunkte, kurze Zyklen, Attraktorverhalten auf Ebene vs. Sphäre vs. Pseudosphäre (Experiment).
2. **Numerik:** $F_M$, $S_M$, $K_M$, $\delta_M$ für DC-Dreieckskonfigurationen aus `collatz_dc_morley_walter.pdf`.
3. **Asymptotik:** kleine Dreiecke auf Flächen konstanter $K_G$ — Test der Conjecture $K_M \sim c_1 K_G$.
4. **$\mu_M$-Feld:** Triangulierung der mod-12-DC-Geometrie; Vergleich mit Beltrami-Norm; Stabilität unter $T_M$.
5. **Ikosaeder-Dualität:** $F_M$, $S_M$ auf ikosaedrischen $S^2$-Dreiecken; Abgleich mit `eabc_icosahedron_test.py` ($\Delta(M)$, $R_{\mathrm{EABC}}$).
6. **Lean:** optional `MorleyOperator`, `MorleyForm`, `MorleyScale` als Definitionen — **kein** Collatz-`sorry`-Abbau.

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
| `eabc_icosahedron_test.py` | Ikosaeder–Dodekaeder-Dualität, $\Delta(M)$-Test, $R_{\mathrm{EABC}}$ |
| `collatz_ikosaeder_spannung.pdf` | Ikosaeder-Spannung (EABC) |
