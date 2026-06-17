# Morley-Stufen M1 → M2 → M3

**Branch:** `main` (M1: PR #44, Commit `c2fa37c`; M2: `collatz/morley-m2-sensor`)  
**Stand:** Juni 2026  
**Kein Collatz-Beweis.** Operatives Forschungsprotokoll für den geodätischen Morley-Operator $T_M^{(g)}$.

**Methodik:** Tao-Stil — Definition → Variantenvergleich → Modellräume → numerische Evidenz → Conjectures (`collatz_formalisierung_tao_stil.md`).

**Kanonische Theorie:** `collatz_morley_metrik_erweiterung.md` (abgeschlossen nach Merge #42/#43). Dieses Dokument ist nur die **Stufen-Roadmap**.

---

## Warum die Reihenfolge?

Der frühe numerische Lauf (`collatz_morley_tm_numerik.py`, Commit auf `main`) testete direkt $F_M \sim c\,K_G A$. Das war **vorzeitig**: Auf gekrümmten Flächen ist $T_M^{(g)}$ kein Theorem, sondern ein **Definitionsproblem** mit mehreren Realisierungen. Solange die Varianten nicht konsistent sind, ist jede $F_M$-Regression **interpretationslos**.

> **Boxed (Leitregel):**  
> *M1 (Operator) → M2 (Modellräume) → M3 (Sensor). Nicht umgekehrt.*

---

## M1 — Konsistenzsatz erster Ordnung (abgeschlossen)

**Frage:** Liefern vier sinnvolle Realisierungen von $T_M^{(g)}$ auf kleinen geodätischen Dreiecken $\Delta_\varepsilon$ **dieselbe** Morley-Konfiguration bis auf höhere Ordnung?

| Nr. | Variante | Implementierung (numerisch) |
|-----|----------|----------------------------|
| 1 | Geodätische Winkel | Trisektion im Tangentialraum, Schnitt über Großkreise |
| 2 | Lokale Karte | $\log_p$ am Schwerpunkt, euklidischer Morley, $\exp_p$ zurück |
| 3 | Exponential + euklidisch | wie (2), explizit als tangentialeuklidisches Modell |
| 4 | Paralleltransport | Trisektionsrichtungen entlang Kanten transportiert, dann Schnitt |

**Experiment:** Familie $\Delta_\varepsilon$ auf $S^2$ (Seitenwinkel $\varepsilon$), paarweise Abstände

$$d_{ij}(\varepsilon) = \max_k \operatorname{dist}\bigl(T_{M,i}^{(g)}(\Delta_\varepsilon)_k,\; T_{M,j}^{(g)}(\Delta_\varepsilon)_k\bigr).$$

**Skalierung:** log-log-Regression $\log d_{ij} \sim \alpha \log \varepsilon$.

### Boxed — Hauptergebnis M1 (numerisch)

$$\boxed{d\bigl(T_{M,i}^{(g)}, T_{M,j}^{(g)}\bigr) = O(\varepsilon^3)\quad\text{— nicht }O(\varepsilon)\text{ oder }O(\varepsilon^2).}$$

**Interpretation (heuristisch, numerisch gestützt):** $T_M^{(g)}$ ist **wohldefiniert bis zur zweiten Ordnung** in $\varepsilon$. Die vier Realisierungen unterscheiden sich wie Chart-, Exponential- und Normal-Koordinaten — erst in dritter Ordnung. Variante 2 und 3 sind numerisch identisch; alle übrigen Paare zeigen log-log-Steigung $\alpha \approx 3{,}01$ bei $R^2 > 0{,}999$.

**Epistemisches Label:** **Experiment** — kein Theorem; „Konsistenzsatz erster Ordnung“ im **heuristischen** Sinn (Operator konsistent bis $O(\varepsilon^2)$ in den Daten).

**Konsequenz für M3:** $F_M \sim c\,K_G A$ ist **plausibler** als vor M1: alle Varianten teilen denselben Hauptterm; verbleibende $O(\varepsilon^3)$-Differenzen stören die Sensorik nicht systematisch in erster Ordnung.

**Artefakte:** `collatz_morley_tm_numerik.py m1`, `collatz_morley_m1_konsistenz.json`, `tests/test_morley_m1_konsistenz.py`.

**PR #44:** *Morley M1: Operator-Konsistenz der vier $T_M^{(g)}$-Varianten* — auf `main` gemergt (`c2fa37c`).

### Tao living spreadsheet (Zeile M1)

| ID | Ebene | Aussage | Status | Artefakt |
|----|-------|---------|--------|----------|
| **M1** | Experiment | $d(T_{M,i}^{(g)}, T_{M,j}^{(g)}) = O(\varepsilon^3)$ auf $S^2$; $T_M^{(g)}$ wohldefiniert bis 2. Ordnung (heuristisch) | **numerisch gestützt** | `collatz_morley_m1_konsistenz.json` |

---

## κ vs. Morley — methodischer Kontrast

| Aspekt | κ-Zweig (Stufe 2–3) | Morley-Zweig (M1–M3) |
|--------|---------------------|----------------------|
| Kernproblem | Kodierungsinvarianz von $R(k)$, $h_F$ | Definition von $T_M^{(g)}$ auf gekrümmten Flächen |
| Strukturabhängigkeit | **κ-abhängig** — verschiedene κ-Familien | **4 Definitionen**, Differenzen nur $O(\varepsilon^3)$ |
| Gegenbeispiel / Zeuge | $w=\mathrm{BE} \notin \mathcal{L}_{\mathrm{arith}}$ | Paarweise Variantenabstände skalieren mit $\varepsilon^3$ |
| BE-Invarianz | **nicht** invariant unter Block-Encoding | nicht relevant (geometrisch) |
| Stabilität (Juni 2026) | κ-Invarianz offen; BE-Bruch belegt | **Morley-Zweig stabiler** — M1 bestanden |
| Epistemik | Zeuge ≠ Theorem | Experiment → Conjecture (M3) |

---

## M2 — Modellräume (aktuell)

Kontroll- und Testflächen mit bekannter $K_G$, **dieselbe** Dreiecksfamilie $\Delta_\varepsilon$:

| Raum | $K_G$ | Rolle |
|------|-------|-------|
| $\mathbb{R}^2$ | $0$ | Kontrolle: Morley-Satz, $F_M=0$ exakt |
| $S^2$ (Einheitskugel) | $+1$ | positive Krümmung |
| $H^2$ (Hyperboloid) | $-1$ | negative Krümmung |

### M2a — $F_M/A$ vs. $K_G$

Für jedes $\varepsilon$: Flächeninhalt $A(\varepsilon)$, Morley-Form $F_M(\varepsilon)$ mit kanonischer Variante `local_chart`. Vergleich der Mediane $F_M/A$ über $\varepsilon$ zwischen $S^2$ und $H^2$.

**Vorläufiger Befund (Juni 2026):** $F_M/A$ auf $S^2$ und $H^2$ **gleiche Größenordnung** ($\sim 10^{-5}$ bei $\varepsilon \in [0{,}05, 0{,}35]$); euklidische Kontrolle $F_M \approx 0$.

### M2b — Exponentenfit

$$\boxed{F_M \stackrel{?}{=} c\,|K_G|^{\alpha}\,A^{\beta}}\qquad\text{(}\alpha,\beta\text{ aus Daten — nicht vorausgesetzt).}$$

Log-lineare Regression auf $K_G \neq 0$-Punkte: $\log F_M \sim \log c + \alpha\log|K_G| + \beta\log A$.

**Vorläufiger Befund:** $\beta \approx 2{,}0$ ($R^2 \approx 0{,}9998$); $\alpha \approx 0$ bei festem $|K_G|=1$ (**nicht identifizierbar** ohne Krümmungs-Skalenvariation). Naives $F_M \sim c\,K_G A$ trennt $S^2$/$H^2$ bei gleichem $|F_M|$ schlecht ($R^2 \ll 1$) — Vorzeichen von $K_G$ allein reicht nicht.

**Artefakte:** `collatz_morley_tm_numerik.py m2`, `collatz_morley_m2_sensor.json`, `tests/test_morley_m2_sensor.py`.

**Epistemisches Label:** **Experiment**.

**Gate M3:** `m1_gate_passed` in M2-JSON; M3-Subcommand bleibt `m3`.

---

## M3 — Morley-Sensor (nach M1 + M2)

Erst wenn $T_M^{(g)}$ feststeht und M2 die Modellräume abdeckt:

$$F_M(\Delta) = \sum_i (\theta_i^M - \pi/3)^2, \qquad
\text{Test: } F_M(\Delta) \stackrel{?}{=} c\,K_G(p)\,A + O(A^2).$$

Sekundär: $F_M \propto K_G^2 A^2$ bzw. $F_M \propto A^\beta$ (M2b-Hinweis $\beta \approx 2$) als Alternativ-Conjectures.

**CLI:** `collatz_morley_tm_numerik.py m3` (hinter Gate, nicht Default).

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

1. ~~M1 auswerten~~ — **erledigt** ($O(\varepsilon^3)$).
2. ~~Kanoniche Variante~~ — `local_chart` (bzw. äquivalent `exp_euclidean`).
3. ~~M2: $H^2$-Patch~~ — **implementiert**; M2b mit variablen Radien $R$ für $\alpha$-Identifikation.
4. M3: $F_M$-Skalierung mit gewählter Variante und M2-Kontext ($\beta$).
5. Optional: Lean-Definitionen (`MorleyOperator`) — ohne Collatz-`sorry`.
