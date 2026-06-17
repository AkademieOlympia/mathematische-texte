# Morley-Stufen M1 → M2 → M3

**Branch:** `main` (nach Merge PR #45, Commit `7b2c92c`)  
**Stand:** Juni 2026  
**Kein Collatz-Beweis.** Operatives Forschungsprotokoll für den geodätischen Morley-Operator $T_M^{(g)}$.

**Methodik:** Tao-Stil — Definition → Variantenvergleich → Modellräume → numerische Evidenz → Conjectures (`collatz_formalisierung_tao_stil.md`).

**Kanonische Theorie:** `collatz_morley_metrik_erweiterung.md` (abgeschlossen nach Merge #42/#43). Dieses Dokument ist die **Stufen-Roadmap** und Synthese der κ–Morley-Asymmetrie.

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
| Nächste offene Frage | $\mathcal{L}_{\mathrm{arith}}^*$ definieren | Misst $T_M^{(g)}$ tatsächlich Krümmung? |

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

M1 testet die **Definitionskonsistenz** der linken Seite; M2 testet die **Krümmungskopplung** von $F_M$ als Sensorgröße.

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

> **Boxed (M1-Bedeutung):**  
> *Die zweite Ordnung scheint universell zu sein* — nicht nur „vier Definitionen sehen ähnlich aus“, sondern: Abweichungen zwischen Realisierungen liegen **jenseits** der Krümmungsterme ($O(\varepsilon^3)$), während die 2.-Ordnung-Struktur gemeinsam bleibt.

**Artefakte:** `collatz_morley_tm_numerik.py m1`, `collatz_morley_m1_konsistenz.json`, `tests/test_morley_m1_konsistenz.py`.

**Epistemisches Label:** **Experiment** — numerische Evidenz, kein Theorem.

---

## M2 — Modellräume (aktuell, PR #45)

Kontroll- und Testflächen mit bekannter $K_G$, **dieselbe $\varepsilon$-Familie** $\{0{,}05, 0{,}08, 0{,}12, 0{,}18, 0{,}25, 0{,}35\}$:

| Raum | $K_G$ | Rolle |
|------|-------|-------|
| $\mathbb{R}^2$ | $0$ | Kontrolle: Morley-Satz, $F_M \approx 0$ |
| $S^2$ (Einheitskugel) | $+1$ | positive Krümmung |
| $H^2$ (Hyperboloid) | $-1$ | negative Krümmung |

**Kanoniche Variante:** `local_chart` (nach M1-Evidenz).

### M2-Protokoll: zwei Welten

**Welt 1:** $F_M \propto |K_G|$ — Krümmungsbetrag, Vorzeichen irrelevant.  
**Welt 2:** $F_M(K{>}0) \neq F_M(K{<}0)$ — Vorzeichen-Sensor.

Empirisch: $F_M(\varepsilon, K_G)$ messen; fit $F_M = c\,|K_G|^\alpha A^\beta + o(A^\beta)$ — **Exponenten aus Daten**, nicht vorausgesetzt.

### M2a/b/c — Geometrietabelle (Pflicht vor M3)

Median über die $\varepsilon$-Familie, Variante `local_chart`, Stand Juni 2026:

| Geometrie | $F_M$ | $F_M/A$ | $F_M/A^2$ |
|-----------|-------|---------|-----------|
| Ebene $K=0$ | $6{,}2\times 10^{-31}$ | $2{,}8\times 10^{-29}$ | $1{,}0\times 10^{-27}$ |
| Sphäre $K=+1$ | $4{,}1\times 10^{-7}$ | $1{,}2\times 10^{-5}$ | $3{,}8\times 10^{-4}$ |
| Hyperbolisch $K=-1$ | $4{,}1\times 10^{-7}$ | $1{,}2\times 10^{-5}$ | $3{,}9\times 10^{-4}$ |

**Lesart M2a:** $F_M \approx 0$ auf der Ebene (Kontrolle); $F_M > 0$ auf $S^2$ und $H^2$.  
**Lesart M2b:** $F_M(S^2) \approx F_M(H^2)$ bei gleichem $\varepsilon$ (Median-Ratio $S^2/H^2 \approx 1{,}001$) → **Welt 1** (Abhängigkeit von $|K_G|$, nicht Vorzeichen); **Welt 2 nicht gestützt**.  
**Lesart M2c:** $F_M/A^2$ auf gekrümmten Flächen $\sim 10^{-4}$, auf der Ebene $\sim 10^{-27}$ — konsistent mit $F_M \propto A^2$ in dieser Bandbreite.

### M2b — Exponentenfit (nach ausgefüllter Tabelle)

$$F_M \stackrel{?}{=} c\,|K_G|^\alpha A^\beta \qquad (\alpha,\beta \text{ aus Daten}).$$

Erste Anpassung (nur $K_G \neq 0$): $c \approx 3{,}9\times 10^{-4}$, $\alpha \approx 0$, $\beta \approx 2$, $R^2 \approx 0{,}9998$ — $F_M \propto A^2$, schwache $|K_G|$-Kopplung in dieser $\varepsilon$-Bandbreite.

**Artefakte:** `collatz_morley_tm_numerik.py m2`, `collatz_morley_m2_sensor.json`, `tests/test_morley_m2_sensor.py`.

**Epistemisches Label:** **Experiment**.

> **Boxed (M2-Gate):**  
> *M3 erst starten, wenn die M2a/b/c-Tabelle vollständig ist und die Vorzeichenstruktur dokumentiert — erfüllt (Juni 2026). M3 bleibt Conjecture-Ebene.*

---

## Ikosaeder — natürliche Testumgebung

Ein regulärer Ikosaeder hat **20 dreieckige Flächen** — eine diskrete Triangulation mit bekannter Kombinatorik. `eabc_icosahedron_test.py` konstruiert Ikosaeder-Ecken, Kanten und Dreiecksflächen auf der Einheitskugel und koppelt EABC-Labels an mod-12-Restklassen. Das liefert eine **Brücke** von diskreter Triangulation zu $T_M^{(g)}$ auf $S^2$: jedes Ikosaeder-Dreieck ist ein Kandidat für $F_M(\Delta)$ und Morley-Fluss $\Delta_{n+1} = T_M(\Delta_n)$. Noch **nicht** durchgeführt: systematischer M2-Lauf auf allen 20 Flächen (nächster Schritt nach M2b-Verfeinerung).

---

## M3 — Morley-Sensor (nach M1 + M2, noch nicht aktiv)

Erst wenn M2-Vorzeichenstruktur und Geometrietabelle geklärt (s.o.):

$$F_M(\Delta) = \sum_i (\theta_i^M - \pi/3)^2, \qquad
\text{Test: } F_M(\Delta) \stackrel{?}{=} c\,K_G(p)\,A + O(A^2).$$

Sekundär: $F_M \propto K_G^2 A^2$ als Alternativ-Conjecture — **M2 deutet eher auf** $F_M \propto A^2$ **ohne starke $K_G$-Kopplung**.

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
2. ~~M2a/b/c: Geometrietabelle~~ — ausgefüllt (PR #45 + Synthese).
3. M2b verfeinern: $\alpha,\beta$ in breiterer $\varepsilon$-Bandbreite; Ikosaeder-20-Flächen-Lauf.
4. M3: $F_M$-Skalierung als Conjecture-Test — **erst nach** M2b-Verfeinerung.
5. Optional: Lean-Definitionen (`MorleyOperator`) — ohne Collatz-`sorry`.
