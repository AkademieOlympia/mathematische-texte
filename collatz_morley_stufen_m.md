# Morley-Stufen M1 → M2 → M3

**Branch:** `collatz/morley-m1-konsistenz`  
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

## M1 — Operator-Konsistenz (jetzt)

**Frage:** Liefern vier sinnvolle Realisierungen von $T_M^{(g)}$ auf kleinen geodätischen Dreiecken $\Delta_\varepsilon$ **dieselbe** Morley-Konfiguration bis auf $O(\varepsilon^2)$?

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
| $\approx 2$ | Varianten konsistent im Grenzfall → $T_M^{(g)}$ wohl definiert |
| $\approx 1$ | systematische Erstordnungs-Abweichung → Definition offen |
| sonst | Zwischenordnung / numerisches Rauschen |

**Artefakte:** `collatz_morley_tm_numerik.py m1`, `collatz_morley_m1_konsistenz.json`, `tests/test_morley_m1_konsistenz.py`.

**Epistemisches Label:** **Experiment** (Definitionsvergleich).

---

## M2 — Modellräume (nach M1)

Kontroll- und Testflächen mit bekannter $K_G$:

| Raum | $K_G$ | Status |
|------|-------|--------|
| $\mathbb{R}^2$ | $0$ | Kontrolle: Morley-Satz, $F_M=0$ exakt |
| $S^2$ (Radius $R$) | $+1/R^2$ | M1-Daten auf Einheitskugel |
| $H^2$ (Radius $R$) | $-1/R^2$ | **noch nicht implementiert** |

**Zweck:** Prüfen, ob die in M1 gewählte Variante auf allen Konstant-Krümmungs-Modellen stabil ist.

**Epistemisches Label:** **Experiment**.

---

## M3 — Morley-Sensor (nach M1 + M2)

Erst wenn $T_M^{(g)}$ feststeht:

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

1. M1 auswerten: $\alpha \approx 2$ für alle Paare?
2. Bei Konsistenz: eine Variante als **kanonische** $T_M^{(g)}$-Wahl dokumentieren.
3. M2: $H^2$-Patch implementieren.
4. M3: $F_M$-Skalierung nur mit gewählter Variante.
5. Optional: Lean-Definitionen (`MorleyOperator`) — ohne Collatz-`sorry`.
