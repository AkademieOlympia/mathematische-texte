# Morley-Stufen M1 → M2 → M3

**Branch:** `collatz/morley-m2-sensor`  
**Stand:** Juni 2026  
**Kein Collatz-Beweis.** Operatives Forschungsprotokoll für den geodätischen Morley-Operator $T_M^{(g)}$.

**Methodik:** Tao-Stil — Definition → Variantenvergleich → Modellräume → numerische Evidenz → Conjectures (`collatz_formalisierung_tao_stil.md`).

**Kanonische Theorie:** `collatz_morley_metrik_erweiterung.md` (abgeschlossen nach Merge #42/#43). Dieses Dokument ist nur die **Stufen-Roadmap**.

---

## Epistemischer Rahmen (Reviewer-Vorgabe)

**Kein Theorem.** Beobachtete Skalierung $O(\varepsilon^3)$ in M1 ist **numerische Evidenz**, kein bewiesener Satz. Die Formel $T_{M,1}=T_{M,2}+O(\varepsilon^3)$ wird **nicht** als Theorem behauptet.

> **Boxed (korrekte Lesart nach M1):**  
> *M1 liefert starke numerische Evidenz für eine universelle lokale Struktur von $T_M^{(g)}$ — kein Theorem.*

### Vier Szenarien vor M1

| Szenario | Skalierung $d_{ij}(\varepsilon)$ | Folge |
|----------|-----------------------------------|-------|
| **A** | $O(1)$ | Operatoren grundsätzlich verschieden → Projekt tot |
| **B** | $O(\varepsilon)$ | Kein kanonischer Operator |
| **C** | $O(\varepsilon^2)$ | Krümmungsebene — definitionsabhängig |
| **D** | $O(\varepsilon^3)$ | **beobachtet** — Abweichungen jenseits der Krümmungsterme |

M1-Ergebnis (PR #44): Szenario **D** — log-log-Steigung $\approx 3$ für alle nicht-trivialen Variantenpaare.

> **Boxed (revidierte Conjecture):**  
> *Alle vernünftigen Morley-Konstruktionen besitzen dieselbe universelle 2.-Ordnung-Geometrie.*

### PR #44 — Reviewer-Paragraph

> Die Resultate von M1 rechtfertigen die weitere Untersuchung des geodätischen Morley-Operators. Die beobachtete Skalierung $O(\varepsilon^3)$ deutet darauf hin, dass verschiedene natürliche Konstruktionen von $T_M^{(g)}$ dieselbe lokale Geometrie bis zur zweiten Ordnung teilen. Dies ist kein Beweis eines Konsistenzsatzes, aber eine starke numerische Evidenz dafür, dass der Operator nicht vollständig konstruktionsabhängig ist.

### Abgrenzung κ vs. Morley

| Zweig | Ergebnis |
|-------|----------|
| **κ** (Stufe 3) | Invarianz zerbricht unter $\kappa_2$ |
| **Morley** (M1) | Invarianz bis Ordnung 2 (numerisch) |

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

**Artefakte:** `collatz_morley_tm_numerik.py m1`, `collatz_morley_m1_konsistenz.json`, `tests/test_morley_m1_konsistenz.py`.

**Epistemisches Label:** **Experiment** — numerische Evidenz, kein Theorem.

---

## M2 — Modellräume (aktuell)

Kontroll- und Testflächen mit bekannter $K_G$, **dieselbe $\varepsilon$-Familie**:

| Raum | $K_G$ | Rolle |
|------|-------|-------|
| $\mathbb{R}^2$ | $0$ | Kontrolle: Morley-Satz, $F_M \approx 0$ |
| $S^2$ (Einheitskugel) | $+1$ | positive Krümmung |
| $H^2$ (Hyperboloid) | $-1$ | negative Krümmung |

**Kanoniche Variante:** `local_chart` (nach M1-Evidenz).

### M2a — Vorzeichenstruktur (Priorität)

**Vor** jedem Exponentenfit $F_M = c\,|K_G|^\alpha A^\beta$:

| Test | Erwartung |
|------|-----------|
| $K_G=0$ (Ebene) | $F_M \approx 0$ |
| $K_G>0$ (Kugel) | $F_M > 0$ |
| $K_G<0$ (hyperbolisch) | $F_M > 0$ |
| **stärker:** $F_M(K{>}0) \neq F_M(K{<}0)$ | Sensor erkennt Krümmungsvorzeichen |

Erste Daten: $F_M>0$ auf $S^2$ und $H^2$, $F_M\approx 0$ auf $\mathbb{R}^2$. Bei gleichem $\varepsilon$ gilt $F_M(S^2)\approx F_M(H^2)$ → Abhängigkeit von $|K_G|$, nicht vom Vorzeichen. **Noch nicht gezeigt:** $F_M(K{>}0) \neq F_M(K{<}0)$ (stärkerer Vorzeichentest).

### M2b — Exponentenfit (nach M2a)

$$F_M \stackrel{?}{=} c\,|K_G|^\alpha A^\beta \qquad (\alpha,\beta \text{ aus Daten, nicht vorausgesetzt}).$$

Erste Anpassung: $\alpha \approx 0$, $\beta \approx 2$ — $F_M \propto A^2$, schwache $|K_G|$-Kopplung in dieser $\varepsilon$-Bandbreite.

**Artefakte:** `collatz_morley_tm_numerik.py m2`, `collatz_morley_m2_sensor.json`, `tests/test_morley_m2_sensor.py`.

**Epistemisches Label:** **Experiment**.

---

## M3 — Morley-Sensor (nach M1 + M2)

Erst wenn M2-Vorzeichenstruktur geklärt:

$$F_M(\Delta) = \sum_i (\theta_i^M - \pi/3)^2, \qquad
\text{Test: } F_M(\Delta) \stackrel{?}{=} c\,K_G(p)\,A + O(A^2).$$

Sekundär: $F_M \propto K_G^2 A^2$ als Alternativ-Conjecture.

**CLI:** `collatz_morley_tm_numerik.py m3` (hinter Flag, nicht Default).

**Epistemisches Label:** **Conjecture** — nur Evidenz, kein Beweis.

---

## Abgrenzung

| Thema | Status |
|-------|--------|
| Collatz / κ-Invarianz | **nicht** Ziel dieser Stufen |
| Morley-Metrik auf glatter Fläche | **verworfen** (Operator statt Metrik) |
| Eigenes Projekt „Morley-Operatoren auf riemannschen Flächen“ | möglich, wenn M1–M3 tragfähig |

---

## Nächste Schritte (kurz)

1. ~~M1 auswerten~~ — $O(\varepsilon^3)$ beobachtet (PR #44, merged).
2. M2a: Vorzeichenstruktur dokumentieren ($F_M$ vs. $K_G$-Vorzeichen).
3. M2b: $\alpha,\beta$ in breiterer $\varepsilon$-Bandbreite verfeinern.
4. M3: $F_M$-Skalierung nur mit geklärter M2-Struktur.
5. Optional: Lean-Definitionen (`MorleyOperator`) — ohne Collatz-`sorry`.
