# Stufe 2B — κ-Robustheit und Entropie von $\mathcal{L}_{\mathrm{arith}}$

**Branch:** `collatz/kappa-robustheit-stufe2b` · **PR:** #40  
**Vorgänger:** PR #39 (Stufe 2, geschlossen — $ \mathcal{L}_{\mathrm{arith}} \subsetneq \mathcal{L}$ experimentell)  
**Kein Collatz-Beweis.**

---

## Motivation (kritische Schwachstelle)

Alles in Stufe 2 hängt an der Operationalisierung

$$\text{„Realisierbar über naive $\kappa$-Präfixe“}.$$

| Fall | Mechanismus | Konsequenz |
|------|-------------|------------|
| **A** | $\kappa$ trifft arithmetische Struktur | $\mathcal{L}_{\mathrm{arith}} \subsetneq \mathcal{L}$ ist echte Grammatik-Lücke |
| **B** | $\kappa$ erzeugt künstliche Ausschlüsse | $R(k)$ misst nur $\kappa$-Artefakte |

> **Boxed (nächster Angriff):** κ-Robustheit — bleiben BE, $F_n$, $R(k)$ unter verschiedenen Kodierungen stabil?

> **Boxed (offene Frage):** Welche **arithmetischen Regeln** erzeugen die verbotenen Wörter? BE ist **Symptom**, nicht Erklärung.

---

## Tao-Einordnung

| Ebene | Inhalt Stufe 2B |
|-------|-----------------|
| **Definition** | Drei dokumentierte $\kappa$-Varianten; $F_n$, $R(k)$, $h_F$ |
| **Experiment** | `collatz_kappa_robustheit.py` → `collatz_kappa_robustheit.json` |
| **Conjecture** | $h_F = \limsup_{n\to\infty} \frac{1}{n}\log|F_n|$ endlich und $R(k)\to 0$ — **offen** |
| **Theorem** | — (keine neuen Beweise) |

---

## Die drei Kodierungen

### $\kappa_1$ — naive mod-12 (Referenz)

Buchstabe$_i = \mathrm{classOf}(\mathrm{iterateU}(n,i))$; Abbruch bei $\mathrm{classOf}=\bot$ ($n\equiv 3,9 \pmod{12}$).

Identisch mit `collatz_l_arith_test.realized_words_by_length`. Dynamiktreu (Lean), nicht injektiv (PR #38).

### $\kappa_2$ — $\nu_2$-Rotation (Blocktyp)

Buchstabe$_i = t^{\nu_2(3n_i+1)\bmod 4}(\mathrm{classOf}(n_i))$ mit $n_i=\mathrm{iterateU}(n,i)$ und $t$ die EABC-Rotation $E\!\to\! A\!\to\! B\!\to\! C\!\to\! E$.

Nutzt die 2-adische Tiefe des geraden Blocks vor dem nächsten ungeraden Wert — gleiche mod-12-Basis, andere Feinstruktur.

### $\kappa_3$ — Successor-/Block-Kodierung

Buchstabe$_i = \mathrm{classOf}(\mathrm{iterateU}(n,i{+}1))$ für $i<k{-}1$; letzter Buchstabe $= \mathrm{classOf}(\mathrm{iterateU}(n,k{-}1))$.

Kodiert die Klasse des **nächsten** odd-to-odd-Zustands (Nachfolger), nicht den aktuellen.

---

## Test 1 — Drei Kodierungen ($n\leq 10^6$)

| Größe | $\kappa_1$ | $\kappa_2$ | $\kappa_3$ |
|-------|-----------|-----------|-----------|
| BE verboten? | **ja** | **nein** | **ja** |
| Min. Gegenbeispiel | BE (L=2) | EAEAA (L=5) | BE (L=2) |
| $R(4)$ | 0,539 | **0,966** | 0,539 |
| $R(8)$ | 0,056 | 0,065 | 0,056 |
| $R(10)$ | 0,0087 | 0,0041 | 0,0088 |

**Interpretation:** $\kappa_1$ und $\kappa_3$ sind **aligniert** (gleiches minimales Gegenbeispiel, nahezu gleiche Ratios). $\kappa_2$ ist **nicht robust**: BE wird realisiert, $R(4)$ springt auf $\approx 97\,\%$ — starkes Indiz für **Fall B** bei $\nu_2$-Rotation. Die Stufe-2-Lücke ist **nicht** κ-unabhängig; treue $\kappa$ (FaithfulKappa) bleibt zwingend.

---

## Test 2 — Bootstrap ($\kappa_1$)

| Limit $n$ | BE verboten? | min. Wort | neue min. Wörter |
|-----------|--------------|-----------|------------------|
| $10^6$ | ja | BE | — |
| $10^7$ | ja | BE | — |

$F_n$-Größen und $R(4),R(8),R(10)$ stabil zwischen $10^6$ und $10^7$ (vgl. JSON). Kein neues minimales Gegenbeispiel unter $\kappa_1$.

*Optional $10^8$:* `--include-1e8` (Laufzeit/Speicher intensiv; nicht im Standardlauf).

---

## Test 3 — Entropie $h_F$

Referenz $|F_n|$ bei $\kappa_1$, $n\leq 10^6$ (Vollliste $n\leq 8$):

| $n$ | $|F_n|$ | $(1/n)\log|F_n|$ |
|-----|---------|------------------|
| 2 | 1 | 0,000 |
| 3 | 6 | 0,598 |
| 4 | 38 | 0,974 |
| 5 | 183 | 1,021 |
| 6 | 807 | 1,100 |
| 7 | 3402 | 1,154 |
| 8 | 13924 | **1,193** |

**Schätzer:** $\displaystyle h_F \approx \max_n \frac{1}{n}\log|F_n| \approx 1{,}19$ (nur $n\leq 8$ — **kein** asymptotischer Beweis).

Die Folge steigt monoton in $(1/n)\log|F_n|$ — Hinweis auf positives Wachstum der verbotenen Sprache, Conjecture-Ebene.

---

## Artefakte

| Datei | Rolle |
|-------|-------|
| `collatz_kappa_robustheit.py` | Pipeline Tests 1–3 |
| `collatz_kappa_robustheit.json` | reproduzierbare Zahlen |
| `tests/test_kappa_robustheit.py` | pytest |

---

## Ehrliche Grenzen

1. Drei $\kappa$ sind **Experimente**, keine Lösung von `faithfulKappaExists`.
2. $\kappa_2$ zeigt: **BE kann κ-realisiert werden** — Hero-Zeuge BE ist κ-abhängig.
3. Bootstrap $10^8$ optional; Standardlauf $10^6$/$10^7$.
4. $h_F$ aus $n\leq 8$ — extrapolieren mit Vorsicht.
5. **Kein** Theorem $\mathrm{BE}\notin\mathcal{L}_{\mathrm{arith}}$ für alle $n$.
6. Arithmetische **Erklärung** der verbotenen Wörter — **offen** (boxed oben).

---

## Nächste Schritte (nicht in diesem PR)

- Arithmetische Ausschlussregeln für BE / $F_n$ suchen (Symptom → Mechanismus).
- Weitere $\kappa$-Kandidaten näher an `FaithfulKappa` testen.
- $F_n$ für $n>8$ bei größerem Rechnerbudget (Vollliste $|L(n)|$ explodiert).
