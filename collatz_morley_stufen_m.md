# Morley-Stufen M1 → M2 → M3

**Branch:** `main` (nach Merge PR #44/#45, Commit `b2ec630`)  
**Stand:** Juni 2026  
**Kein Collatz-Beweis.** Operatives Forschungsprotokoll für den geodätischen Morley-Operator $T_M^{(g)}$.

**Methodik:** Tao-Stil — Definition → Variantenvergleich → Modellräume → numerische Evidenz → Conjectures (`collatz_formalisierung_tao_stil.md`).

**Kanonische Theorie:** `collatz_morley_metrik_erweiterung.md` (abgeschlossen nach Merge #42/#43). Dieses Dokument ist die **Stufen-Roadmap** für den Morley-Zweig; κ und Morley bleiben **parallele Spuren** (Asymmetrietabelle unten), keine künstliche Vereinheitlichung.

---

## Zwei Forschungslogiken (parallel, nicht verschmolzen)

> **Boxed (κ-Zweig):**  
> *Arithmetische Symbolik mit Invarianzproblem* — Kodierungsstabilität, $\mathcal{L}_{\mathrm{arith}}^*$, Stufe 3 (`collatz_stufe3_kappa_invarianz.md`).

> **Boxed (Morley-Zweig):**  
> *Geometrischer Operator mit überraschender lokaler Stabilität* — Definitionsrobustheit von $T_M^{(g)}$, experimentelle Testgröße $F_M$.

> **Boxed (Roadmap):**  
> *PR #44 mergen (M1) → M2 ($F_M$, $G_M$ auf Ebene / $S^2$ / $H^2$) → M3 (dualer Exponentenfit)*  
> — M1 **abgeschlossen** (PR #44); M2 **abgeschlossen** (PR #46); M3 **abgeschlossen** (PR #47, experimentell).

---

## Epistemischer Rahmen (Reviewer-Vorgabe)

**Kein Theorem.** Beobachtete Skalierung $O(\varepsilon^3)$ in M1 ist **numerische Evidenz**, kein bewiesener Satz. Die Formel $T_{M,1}=T_{M,2}+O(\varepsilon^3)$ wird **nicht** als Theorem behauptet.

> **Boxed (korrekte Lesart nach M1):**  
> *M1 liefert starke numerische Evidenz für eine universelle lokale Struktur von $T_M^{(g)}$ — kein Theorem.*

---

## κ vs. Morley — strukturelle Asymmetrie (Hauptfokus)

Die beiden Forschungszweige testen **verschiedene Robustheitsbegriffe** auf scheinbar verwandten Objekten. Das Ergebnis ist **asymmetrisch**, nicht symmetrisch:

| Eigenschaft | κ-Zweig (Stufe 3) | Morley-Zweig (M1–M2) |
|-------------|-------------------|----------------------|
| Primärer Stresstest | BE-Verbot vs. $T_M^{(g)}$ | BE-Verbot vs. vier Operatorvarianten |
| Stresstest-Objekt | $\kappa_2$ ($\nu_2$-Rotation) | geodätische Winkel, lokale Karte, exp-euklidisch, Paralleltransport |
| Beobachtetes Ergebnis | Struktur **zerfällt** (BE nicht invariant) | Struktur **bleibt** (2. Ordnung numerisch gemeinsam) |
| Invarianz-Stärke | schwach (listenabhängig) | stark (numerisch, vier Realisierungen) |
| Hauptrisiko | Kodierungsartefakt | numerische Überschätzung der Konsistenz |
| Nächste offene Frage | $\mathcal{L}_{\mathrm{arith}}^*$ definieren | Korreliert $T_M^{(g)}$ mit lokaler Krümmung? |

**Lesart:** Der κ-Zweig fragt nach **arithmetischer** Kodierungsrobustheit; der Morley-Zweig nach **geometrischer** Definitionsrobustheit. Ein Scheitern in κ schließt keinen Erfolg in Morley aus — und umgekehrt.

> **Boxed (Gesamturteil):**  
> *Der Morley-Zweig besitzt erstmals eine eigene mathematische Fragestellung, die unabhängig von Collatz interessant ist.*

Querverweis (minimal): `collatz_stufe3_kappa_invarianz.md`, Abschnitt „Geometrischer Nebenzweig“.

---

## Theoretischer Hintergrund (Conjecture, nicht Theorem)

Für kleine geodätische Dreiecke $\Delta_\varepsilon$ nahe $p$ mit Gauß-Krümmung $K_G(p)$ und Fläche $A(\Delta_\varepsilon)$:

$$T_M^{(g)}(\Delta_\varepsilon) = T_{\mathrm{euclid}}(\Delta_\varepsilon) + Q(K_G(p))\,\varepsilon^2 + O(\varepsilon^3).$$

- $T_{\mathrm{euclid}}$: euklidischer Morley-Referenzoperator (Morley-Satz).
- $Q(K_G)$: noch **undefiniert** — möglicherweise von $|K_G|$, möglicherweise vorzeichenabhängig.
- $O(\varepsilon^3)$: **numerisch beobachtet** in M1 (Variantenabstände); **nicht** als Theorem behauptet.

M1 testet die **Definitionskonsistenz** der linken Seite; M2 testet die **Krümmungskopplung** von $F_M$ und $G_M$ als experimentelle Testgrößen.

---

## Warum die Reihenfolge?

Der frühe numerische Lauf (`collatz_morley_tm_numerik.py`, Commit auf `main`) testete direkt $F_M \sim c\,K_G A$. Das war **vorzeitig**: Auf gekrümmten Flächen ist $T_M^{(g)}$ kein Theorem, sondern ein **Definitionsproblem** mit mehreren Realisierungen. Solange die Varianten nicht verglichen sind, ist jede $F_M$-Regression **interpretationslos**.

> **Boxed (Leitregel):**  
> *M1 (Operator) → M2 (Modellräume) → M3 (Sensor). Nicht umgekehrt.*

---

## M1 — Variantenvergleich (abgeschlossen, PR #44)

**Frage:** Liefern vier sinnvolle Realisierungen von $T_M^{(g)}$ auf kleinen geodätischen Dreiecken $\Delta_\varepsilon$ **dieselbe** lokale Geometrie bis zur 2. Ordnung?

| Nr. | Variante | Implementierung (numerisch) |
|-----|----------|----------------------------|
| 1 | Geodätische Winkel | Trisektion im Tangentialraum, Schnitt über Großkreise |
| 2 | Lokale Karte | $\log_p$ am Schwerpunkt, euklidischer Morley, $\exp_p$ zurück |
| 3 | Exponential + euklidisch | wie (2), explizit als tangentialeuklidisches Modell |
| 4 | Paralleltransport | Trisektionsrichtungen entlang Kanten transportiert, dann Schnitt |

**Experiment:** Familie $\Delta_\varepsilon$ auf $S^2$ (Seitenwinkel $\varepsilon$), paarweise Abstände

$$d_{ij}(\varepsilon) = \max_k \operatorname{dist}\bigl(T_{M,i}^{(g)}(\Delta_\varepsilon)_k,\; T_{M,j}^{(g)}(\Delta_\varepsilon)_k\bigr).$$

**Skalierung:** log-log-Regression $\log d_{ij} \sim \alpha \log \varepsilon$.

| $\alpha$ | Lesart |
|----------|--------|
| $\approx 3$ | **beobachtet** — $O(\varepsilon^3)$-Abweichung → 2.-Ordnung-Geometrie gemeinsam |
| $\approx 2$ | Übereinstimmung auf Krümmungsebene (Szenario C) |
| $\approx 1$ | systematische Erstordnungs-Abweichung → Definition offen |
| $\approx 0$ | $O(1)$ — Projekt tot (Szenario A) |

M1-Ergebnis (PR #44): Szenario **D** — log-log-Steigung $\approx 3$ für alle nicht-trivialen Variantenpaare.

> M1 zeigt numerisch, dass natürliche Varianten des geodätischen Morley-Operators bis zur zweiten Ordnung dieselbe lokale Struktur zu besitzen scheinen. Damit ist der Operator hinreichend stabil, um in M2 als experimentelle Testgröße geprüft zu werden.

> **Boxed (M1-Bedeutung):**  
> *Die zweite Ordnung scheint universell zu sein* — Abweichungen zwischen Realisierungen liegen **jenseits** der Krümmungsterme ($O(\varepsilon^3)$), während die 2.-Ordnung-Struktur gemeinsam bleibt.

> **Boxed (M1 → M2):**  
> *Der Morley-Zweig ist jetzt eigenständig genug für M2.*

**Artefakte:** `collatz_morley_tm_numerik.py m1`, `collatz_morley_m1_konsistenz.json`, `tests/test_morley_m1_konsistenz.py`.

**Epistemisches Label:** **Experiment** — numerische Evidenz, kein Theorem.

---

## M2 — Krümmungs-Signatur (abgeschlossen, PR #46)

> **G_M ist derzeit ein experimenteller Kandidat** für einen signierten Morley-Krümmungssensor.

> **M2 testet, ob sign(G_M) mit sign(K_G) korreliert.**

M2 prüft $F_M$ und $G_M$ auf Kontroll- und Testflächen mit bekannter $K_G \in \{0, +1, -1\}$, **dieselbe $\varepsilon$-Familie** $\{0{,}05, 0{,}08, 0{,}12, 0{,}18, 0{,}25, 0{,}35\}$.

**Definitionen:**
$$F_M(\Delta) = \sum_i (\theta_i^M - \pi/3)^2, \qquad
G_M(\Delta) = \sum_i (\theta_i^M - \pi/3).$$

> G_M wird als signierter Morley-Fehler eingeführt. Er ist kein bewiesener Krümmungssensor, sondern eine experimentelle Testgröße, mit der geprüft wird, ob der geodätische Morley-Operator zwischen positiver und negativer Krümmung unterscheiden kann.

**Kernpunkt:** $F_M$ ist eine **nichtnegative Testgröße** (Krümmungsstärke-Proxy), **kein** natürlicher Vorzeichen-Indikator. Da $F_M \geq 0$, kann $F_M(K{=}{+}1) \approx F_M(K{=}{-}1)$ bei gleicher Fläche gelten.

> **Boxed (M2-Leitfragen):**  
> 1. Ist Morley bei $K \neq 0$ verzerrt?  
> 2. Korreliert $\operatorname{sign}(G_M)$ mit $\operatorname{sign}(K_G)$?

> **Boxed (M2-Roadmap):**  
> - **M2a:** $F_M$ als Krümmungsstärke-Testgröße  
> - **M2b:** $G_M$ als experimenteller Kandidat für Vorzeichen-Kopplung

**Testmatrix** (Median über $\varepsilon$-Familie, `local_chart`, Quelle: `collatz_morley_m2_sensor.json`):

| Raum | $K$ | $F_M$ | $G_M$ |
|------|-----|-------|-------|
| Ebene | $0$ | $\approx 0$ ($6{,}2\times 10^{-31}$) | $\approx 0$ ($2{,}2\times 10^{-16}$) |
| Sphäre | $+1$ | $> 0$ ($4{,}1\times 10^{-7}$) | $> 0$ ($+1{,}04\times 10^{-3}$) |
| Hyperbolisch | $-1$ | $> 0$ ($4{,}1\times 10^{-7}$) | $< 0$ ($-1{,}04\times 10^{-3}$) |

**Durchbruchskriterium:** $\operatorname{sign}(G_M) = \operatorname{sign}(K_G)$ — **erreicht** (Juni 2026): $G_M > 0$ auf $S^2$, $G_M < 0$ auf $H^2$, $G_M \approx 0$ auf $\mathbb{R}^2$, bei jedem $\varepsilon$ der Familie.

**Vorbehalt:** $G_M$ muss definitionsunabhängig genug sein; die einfache Summe ist der erste Test. Orientierte Varianten (Fläche, Umlaufzahl, geodätischer Exzess) folgen später.

> **Boxed (M2-Befund):**  
> *M2a: $F_M$ korreliert mit Krümmungsstärke ($F_M \approx 0$ bei $K{=}0$, $F_M > 0$ bei $K \neq 0$; $S^2/H^2 \approx 1{,}001$). M2b: $\operatorname{sign}(G_M)$ korreliert mit $\operatorname{sign}(K_G)$ — kombinierte Testgrößen $(F_M, G_M)$ empfohlen. Nur numerische Evidenz, kein bewiesener Krümmungssensor.*

### M2a — Stärke $F_M$

| Geometrie | $F_M$ | $F_M/A$ | $F_M/A^2$ |
|-----------|-------|---------|-----------|
| Ebene $K=0$ | $6{,}2\times 10^{-31}$ | $2{,}8\times 10^{-29}$ | $1{,}0\times 10^{-27}$ |
| Sphäre $K=+1$ | $4{,}1\times 10^{-7}$ | $1{,}2\times 10^{-5}$ | $3{,}8\times 10^{-4}$ |
| Hyperbolisch $K=-1$ | $4{,}1\times 10^{-7}$ | $1{,}2\times 10^{-5}$ | $3{,}9\times 10^{-4}$ |

**Lesart:** $F_M \approx 0$ auf der Ebene; $F_M > 0$ auf $S^2$ und $H^2$; $F_M(S^2) \approx F_M(H^2)$ → Abhängigkeit von $|K_G|$, nicht Vorzeichen.

### M2b — Vorzeichen $G_M$

| Geometrie | $G_M$ (Median) | Vorzeichen |
|-----------|----------------|------------|
| Ebene $K=0$ | $2{,}2\times 10^{-16}$ | $\approx 0$ |
| Sphäre $K=+1$ | $+1{,}04\times 10^{-3}$ | positiv |
| Hyperbolisch $K=-1$ | $-1{,}04\times 10^{-3}$ | negativ |

**Lesart:** $G_M \approx 0$ auf $\mathbb{R}^2$ (Morley-Satz); auf gekrümmten Flächen **stabile Vorzeichentrennung** $S^2 > 0$, $H^2 < 0$ bei jedem $\varepsilon$ der Familie. $|G_M(S^2)| \approx |G_M(H^2)|$ — Betrag analog zu $F_M$, Vorzeichen orthogonal.

### M2c — Exponentenfit $F_M$ (nach ausgefüllter Tabelle)

$$F_M \stackrel{?}{=} c\,|K_G|^\alpha A^\beta \qquad (\alpha,\beta \text{ aus Daten}).$$

Erste Anpassung (nur $K_G \neq 0$): $c \approx 3{,}9\times 10^{-4}$, $\alpha \approx 0$, $\beta \approx 2$, $R^2 \approx 0{,}9998$ — $F_M \propto A^2$, schwache $|K_G|$-Kopplung in dieser $\varepsilon$-Bandbreite.

**Artefakte:** `collatz_morley_tm_numerik.py m2`, `collatz_morley_m2_sensor.json`, `tests/test_morley_m2_sensor.py`.

**Epistemisches Label:** **Experiment** — numerische Evidenz, kein Theorem.

> **Boxed (M2-Gate):**  
> *M2a/b/c-Tabelle, $F_M$- und $G_M$-Vorzeichenstruktur dokumentiert (Juni 2026). M3 (dualer Fit) abgeschlossen — siehe `collatz_morley_gm_beweisversuch.md`.*

---

## Ikosaeder — natürliche Testumgebung (4. Stufe, geplant)

Nach Ebene / $S^2$ / $H^2$ ist der **reguläre Ikosaeder** die natürliche **4. Teststufe**: diskretes Polyeder mit 20 dreieckigen Flächen, bekannter Kombinatorik und Ikosaeder–Dodekaeder-Dualität. `eabc_icosahedron_test.py` konstruiert Ikosaeder-Ecken, Kanten und Dreiecksflächen auf der Einheitskugel und koppelt EABC-Labels an mod-12-Restklassen. Das liefert eine **Brücke** von diskreter Triangulation zu $T_M^{(g)}$ auf $S^2$: jedes Ikosaeder-Dreieck ist ein Kandidat für $F_M(\Delta)$, $G_M(\Delta)$ und Morley-Fluss $\Delta_{n+1} = T_M(\Delta_n)$. Noch **nicht** durchgeführt: systematischer M2-Lauf auf allen 20 Flächen.

---

## M3 — Dualer Exponentenfit (abgeschlossen, PR #47, experimentell)

M2 ist abgeschlossen (Vorzeichenstruktur geklärt). M3 testet **getrennte** Skalierungsgesetze:

$$F_M \stackrel{?}{=} c_F\,|K_G|^\alpha A^\beta, \qquad
G_M \stackrel{?}{=} c_G\,K_G^\alpha A^\beta.$$

Beachte: bei $G_M$ steht $K_G$ **ohne** Betrag — Vorzeichenkopplung ist Teil der Conjecture. Ebene ($K_G=0$) ist Kontrolle und aus den Fits ausgeschlossen.

**Ergebnisse** (`collatz_morley_tm_numerik.py m3`, Quelle: `collatz_morley_m3_beweisversuch.json`):

| Sensor | $c$ | $\alpha$ | $\beta$ | $R^2$ | $n$ |
|--------|-----|----------|---------|-------|-----|
| $F_M$ | $3{,}88\times 10^{-4}$ | $\approx 0$ | $\approx 2$ | $0{,}9998$ | 12 |
| $G_M$ | $3{,}41\times 10^{-2}$ | $\approx 0$ | $\approx 1$ | $0{,}9998$ | 12 |

**Lesart:** $F_M \propto A^2$ (quadratischer Fehler); $G_M \propto A$ (linear, signiert). $\alpha$ ist auf $S^2/H^2$ mit $|K_G|\equiv 1$ nicht identifizierbar.

**$G_M$-Cross-Check:** $\operatorname{sign}(G_M)=\operatorname{sign}(K_G)$ stabil über alle $\varepsilon$ und alle vier M1-Varianten — numerische Korrelation, kein bewiesener Krümmungssensor.

**Dokumentation:** [`collatz_morley_gm_beweisversuch.md`](collatz_morley_gm_beweisversuch.md) — Argumente für/gegen, Grenzen, nächste Schritte.

**CLI:** `collatz_morley_tm_numerik.py m3` → `collatz_morley_m3_beweisversuch.json`.

**Artefakte:** `tests/test_morley_m3_beweisversuch.py`.

**Epistemisches Label:** **Conjecture / Experiment** — nur numerische Evidenz.

> **Boxed (M3-Befund):**  
> *Dualer Fit bestätigt $F_M \propto A^2$ und $G_M \propto A$ mit hohem $R^2$. Vorzeichen von $G_M$ ist über $\varepsilon$ und Varianten stabil — $G_M$ bleibt experimenteller Kandidat, definitionsabhängig und nur numerisch.*

---

## Abgrenzung

| Thema | Status |
|-------|--------|
| Collatz / κ-Invarianz | **nicht** Ziel dieser Stufen |
| Morley-Metrik auf glatter Fläche | **verworfen** (Operator statt Metrik) |
| Eigenes Projekt „Morley-Operatoren auf riemannschen Flächen“ | möglich, wenn M1–M3 tragfähig |

---

## Nächste Schritte (kurz)

1. ~~M1~~ — abgeschlossen (PR #44).
2. ~~M2: Krümmungs-Signatur ($F_M$ + $G_M$)~~ — erweitert (PR #46).
3. ~~**M3:** Dualer Exponentenfit $F_M \sim c_F |K|^\alpha A^\beta$, $G_M \sim c_G K^\alpha A^\beta$~~ — experimentell (PR #47).
4. Optional: Ikosaeder-20-Flächen-Lauf (4. Stufe nach $R^2/S^2/H^2$).
5. Optional: orientiertes $G_M$, Varianten-Invarianz von $c_G$ und $\beta$.
6. Optional: Lean-Definitionen (`MorleyOperator`) — ohne Collatz-`sorry`.
