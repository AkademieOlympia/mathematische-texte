# **Stufe 2B: κ-Robustheit widerlegt die starke Interpretation von BE.**

**Branch:** `collatz/kappa-robustheit-stufe2b` · **PR:** #40  
**Vorgänger:** PR #39 (Stufe 2 — $\mathcal{L}_{\mathrm{arith}} \subsetneq \mathcal{L}$ experimentell)  
**Kein Collatz-Beweis.**

---

## Narrativ (Tao-Stil)

**Hypothese getestet:** Die Stufe-2-Befunde (BE verboten, kleine $R(k)$, wachsendes $|F_n|$) sind κ-unabhängig und spiegeln echte arithmetische Struktur wider.

**Abhängigkeit gefunden:** Die konkrete Verbotsstruktur — welches Wort minimal verboten ist, wie groß $|F_n|$ ist — hängt von der gewählten Kodierung $\kappa$ ab.

**Behauptung verfeinert:** Nicht „BE ist arithmetisch verboten“, sondern „die realisierbare Sprache $\mathcal{L}_{\mathrm{arith}}^\kappa$ ist stark ausgedünnt gegenüber der Grammatik $\mathcal{L}$“ — mit κ-sensitiven Details.

| PR | Aussage |
|----|---------|
| **#39** | $\mathcal{L}_{\mathrm{arith}} \subsetneq \mathcal{L}$ sieht **stark** aus (BE, $R(10)\approx 0{,}87\,\%$, $|F_n|$-Katalog). |
| **#40** | Aber die **exakte** Verbotsstruktur hängt von $\kappa$ ab. |

---

## Kernbefund

**BE ist kein kodierungsunabhängiger Zeuge.** Unter $\kappa_2$ wird BE realisiert ($\mathrm{BE} \in \mathcal{L}_{\mathrm{arith}}^{\kappa_2}$, aber $\mathrm{BE} \notin \mathcal{L}_{\mathrm{arith}}^{\kappa_1}$). BE ist ein **κ-abhängiges Symptom**, kein arithmetischer Invariant.

| Kodierung | BE verboten? | Min. Gegenbeispiel | $R(10)$ |
|-----------|--------------|-------------------|---------|
| $\kappa_1$ (naiv mod-12) | **ja** | BE (L=2) | 0,0087 |
| $\kappa_2$ ($\nu_2$-Rotation) | **nein** | EAEAA (L=5) | **0,0041** |
| $\kappa_3$ (Successor-Shift) | **ja** | BE (L=2) | 0,0088 |

**$\kappa_1 \sim \kappa_3$:** Nahezu identisch (Successor-Shift = gleicher Informationskanal). Das beweist **keine** Robustheit — es zeigt nur, dass zwei ähnliche Kodierungen ähnliche Artefakte erzeugen.

**$\kappa_2$** ist der echte Stresstest: **zerstört BE** als universellen Zeugen. Gleichzeitig bleibt $R(10)\approx 0{,}0041$ — sogar **kleiner** als bei $\kappa_1$.

> **Boxed (Hauptthese):**  
> *Die Verbotsstruktur ist κ-sensitiv, die Sprachverdünnung aber möglicherweise robust.*

> **Boxed (Schluss):**  
> *PR #40 schwächt BE, stärkt aber die Entropiefrage.*

Stufe 2B verschiebt den Befund von einer **verbotenen-Wort-Hypothese** zu einer **Entropie-/Dünnheits-Hypothese**. BE ist kein stabiler arithmetischer Invariant, sondern ein Artefakt oder Symptom der gewählten Kodierung. Robust erscheint eher die **starke Ausdünnung** der realisierbaren Sprache gegenüber der Grammatik $\mathcal{L}$.

**Neue saubere Conjecture (heuristisch):** Nicht die verbotenen Wörter selbst sind invariant, sondern das **starke Ausdünnen** von $\mathcal{L}_{\mathrm{arith}}$ gegenüber $\mathcal{L}$.

---

## Tao-Einordnung

| Ebene | Inhalt Stufe 2B | Status |
|-------|-----------------|--------|
| **Definition** | Drei dokumentierte $\kappa$-Varianten; $F_n := \mathcal{L}(n)\setminus\mathcal{L}_{\mathrm{arith}}^\kappa(n)$; $R(k)$; $h_F$ | formalisiert |
| **Experiment** | `collatz_kappa_robustheit.py` → `collatz_kappa_robustheit.json` | reproduzierbar |
| **Conjecture** | $R(k)\to 0$; $h_F$ endlich — κ-robust? | **offen** |
| **Theorem** | — | keine neuen Beweise |

---

## Definition — die drei Kodierungen

### $\kappa_1$ — naive mod-12 (Referenz, PR #38/#39)

Buchstabe$_i = \mathrm{classOf}(\mathrm{iterateU}(n,i))$; Abbruch bei $\mathrm{classOf}=\bot$ ($n\equiv 3,9 \pmod{12}$).

Identisch mit `collatz_l_arith_test.realized_words_by_length`. Dynamiktreu (Lean `naiveKappa_shift`), nicht injektiv.

### $\kappa_2$ — $\nu_2$-Rotation (Stresstest)

Buchstabe$_i = t^{\nu_2(3n_i+1)\bmod 4}(\mathrm{classOf}(n_i))$ mit $n_i=\mathrm{iterateU}(n,i)$ und Rotation $t$: $E\!\to\! A\!\to\! B\!\to\! C\!\to\! E$.

Nutzt die 2-adische Tiefe des geraden Blocks vor dem nächsten ungeraden Wert — gleiche mod-12-Basis, **andere Feinstruktur**. Experimentell; testet, ob Verbotsmuster von $\nu_2$-Information stammen.

### $\kappa_3$ — Successor-/Block-Kodierung

Buchstabe$_i = \mathrm{classOf}(\mathrm{iterateU}(n,i{+}1))$ (Nachfolger-Klasse statt aktueller).

Shift um einen odd-to-odd-Schritt; $\kappa_1$-äquivalente Informationslage → erwartete Alignment mit $\kappa_1$.

---

## Experiment — Test 1: Drei Kodierungen ($n\leq 10^6$)

### BE und minimale Gegenbeispiele

| $\kappa$ | BE verboten? | Min. Gegenbeispiel | $|F_2|$ | $|F_5|$ | $|F_8|$ |
|----------|--------------|-------------------|---------|---------|---------|
| $\kappa_1$ | ja | BE | 1 | 183 | 13 975 |
| $\kappa_2$ | **nein** | EAEAA (L=5) | 0 | 17 | 13 353 |
| $\kappa_3$ | ja | BE | 1 | 183 | 13 979 |

Unter $\kappa_2$ existieren **keine** verbotenen Wörter der Länge 2–4 ($|F_2|=|F_3|=|F_4|=0$); erst ab Länge 5 trennt sich $\mathcal{L}_{\mathrm{arith}}^{\kappa_2}$ von $\mathcal{L}$.

### Ratios $R(k) = |L_{\mathrm{arith}}^\kappa \cap L(k)| / |L(k)|$

| $k$ | $\kappa_1$ | $\kappa_2$ | $\kappa_3$ |
|-----|-----------|-----------|-----------|
| 4 | 0,539 | **0,966** | 0,539 |
| 8 | 0,056 | 0,065 | 0,056 |
| 10 | 0,0087 | **0,0041** | 0,0088 |

**Lesart:** $\kappa_2$ realisiert bei $k=4$ fast die gesamte Grammatik ($R(4)\approx 97\,\%$) — BE ist dort **nicht** ausgeschlossen. Bei $k=10$ fällt $R(10)$ trotzdem auf $\approx 0{,}4\,\%$: die **Dünnheit** bei großem $k$ überlebt den κ-Wechsel besser als der BE-Zeuge.

---

## Experiment — Test 2: Bootstrap $\kappa_1$

| Limit $n$ | BE verboten? | min. Wort | $R(10)$ | neue min. Wörter |
|-----------|--------------|-----------|---------|------------------|
| $10^6$ | ja | BE | 0,0087 | BE |
| $10^7$ | ja | BE | 0,0119 | — |

BE bleibt unter $\kappa_1$ bis $10^7$ **stabil** (experimentell gestützt, kein Theorem). $R(10)$ steigt leicht ($0{,}0087 \to 0{,}0119$) — mehr Starts, mehr Realisierungen; kein neues minimales Gegenbeispiel.

*Optional $10^8$:* `--include-1e8` (nicht im Standardlauf).

---

## Experiment — Test 3: Entropie $h_F$

Referenz $|F_n|$ bei $\kappa_1$, $n\leq 10^6$:

| $n$ | $|F_n|$ | $(1/n)\log|F_n|$ |
|-----|---------|------------------|
| 2 | 1 | 0,000 |
| 3 | 6 | 0,597 |
| 4 | 38 | 0,909 |
| 5 | 183 | 1,042 |
| 6 | 807 | 1,116 |
| 7 | 3 403 | 1,162 |
| 8 | 13 975 | **1,193** |

**Schätzer:** $h_F \approx \max_n (1/n)\log|F_n| \approx 1{,}19$ (nur $n\leq 8$ — **kein** asymptotischer Beweis).

Bei $\kappa_2$ (andere $|F_n|$-Folge): $h_F \approx 1{,}19$ bei $n=8$ — ähnliche Größenordnung, **andere** konkrete Wortlisten.

---

## Conjecture — verfeinerte Forschungsfragen

1. **Dünnheits-Conjecture (heuristisch):** $\displaystyle\lim_{k\to\infty} R(k) = 0$ für jede „vernünftige“ $\kappa$ — **offen**.
2. **Entropie-Conjecture (heuristisch):** $h_F$ endlich und κ-robust in der Größenordnung — **offen**.
3. **Verbotene-Wörter-Conjecture:** **widerlegt** in der starken Form („BE ist kodierungsunabhängig verboten“).

---

## Forschungsstatus nach Stufe 2B

| Aussage | Status |
|---------|--------|
| BE verboten bei $\kappa_1$ | **experimentell** stabil bis $10^7$ |
| BE kodierungsunabhängig verboten | **widerlegt** ($\kappa_2$ realisiert BE) |
| $R(k)$ klein für $\kappa_1,\kappa_2,\kappa_3$ | **experimentell gestützt** ($R(10)<1{,}2\,\%$) |
| konkrete $F_n$-Listen | **κ-abhängig** |
| $R(k)\to 0$ | **offen** |
| kodierungsfreie Definition von $\mathcal{L}_{\mathrm{arith}}$ | **zentral offen** |

---

## Ehrliche Grenzen

1. Drei $\kappa$ sind **dokumentierte Experimente**, keine Lösung von `faithfulKappaExists`.
2. $\kappa_2$ zeigt: BE kann κ-realisiert werden — der Hero-Zeuge aus PR #39 ist **κ-abhängig**.
3. $\kappa_1 \approx \kappa_3$ beweist keine Robustheit, nur Kanaläquivalenz.
4. Bootstrap und Vergleich auf ungeraden $n\leq$ jeweilige Grenze; $n\equiv 3,9\pmod{12}$ trägt nicht bei.
5. $F_n$-Vollliste nur $n\leq 8$; $|L(n)|$ explodiert.
6. $h_F$ aus endlichen $n$ — keine asymptotische Aussage.
7. **Kein** Theorem $\mathrm{BE}\notin\mathcal{L}_{\mathrm{arith}}$ für alle $n$ und alle $\kappa$.
8. Arithmetische **Erklärung** der Verbotsmuster — **offen**.

---

## Nächster Schritt: κ-Klassifikation (nicht $10^8$)

Vor größeren Bootstrap-Läufen: Kodierungen **klassifizieren** nach:

| Kriterium | Frage |
|-----------|-------|
| dynamiktreu? | Stimmen Präfixe mit Bahn überein? |
| injektiv? | Verliert $\kappa$ Information? |
| shiftäquivalent? | Gleicher Informationskanal wie $\kappa_1$? |
| rotationssensitiv? | Reagiert auf $\nu_2$-Rotation wie $\kappa_2$? |
| stabilisiert $R(k)$? | Bleibt Dünnheit bei wachsendem $k$? |
| erhält minimale Verbotslängen? | Wann trennt sich $\mathcal{L}_{\mathrm{arith}}^\kappa$ von $\mathcal{L}$? |

Ziel: eine **FaithfulKappa**-Näherung finden, die Dünnheit erklärt, ohne BE als Schein-Invariant zu missbrauchen.

---

## Artefakte

| Datei | Rolle |
|-------|-------|
| `collatz_kappa_robustheit.py` | Pipeline Tests 1–3 |
| `collatz_kappa_robustheit.json` | reproduzierbare Zahlen ($n\leq 10^6$; Bootstrap $10^6$/$10^7$) |
| `tests/test_kappa_robustheit.py` | pytest (inkl. κ₂-BE-Realisierung) |
| `collatz_generalangriff_2026.md` | Verweis Stufe 2B |
