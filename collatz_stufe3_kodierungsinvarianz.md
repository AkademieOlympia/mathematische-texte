# **Stufe 3: Kodierungsinvarianz**

**Stand:** Juni 2026 · **Vorgänger:** PR #39 (Stufe 2), PR #40 (Stufe 2B) — beide auf `main`  
**Kein Collatz-Beweis.** Dieses Dokument definiert die nächste Forschungsstufe; es greift Collatz nicht direkt an.

**Methodik:** Tao-Stil (IEANTN/PNT+-Parallele) — Definition / Zeuge / Experiment / Theorem / Conjecture strikt trennen (`collatz_formalisierung_tao_stil.md`).

---

## Ausgangslage: die Kombination PR #39 + PR #40

Die beiden PRs sind **zusammen** der eigentliche Befund — nicht isoliert.

| PR | Kern | Epistemisches Label |
|----|------|---------------------|
| **#39** | $\mathcal{L}_{\mathrm{arith}} \subsetneq \mathcal{L}$ — massive Sprachverdünnung; BE als erster Zeuge; $R(10)\approx 0{,}87\,\%$ bei $\kappa_1$ | **Experiment** (stark, reproduzierbar) |
| **#40** | BE ist **nicht** kodierungsunabhängig: $\kappa_2$ realisiert BE; $\kappa_1/\kappa_3$ verbieten es; konkrete $F_n$-Listen sind κ-abhängig | **Negativtest** (starke BE-Interpretation widerlegt) |

> **Boxed (Kombinationsbefund):**  
> *Konkrete Verbotslisten sind κ-abhängig; die Ausdünnung $R(k)$ erscheint dagegen möglicherweise robust.*

Das verschiebt die Forschung von „welches Wort ist verboten?“ zu „welche Eigenschaften der Realisierbarkeitslücke überleben Kodierungswechsel?“

### Referenztabelle $R(10)$ ($n\leq 10^6$)

| Kodierung | $R(10)$ | BE verboten? | min. Gegenbeispiel |
|-----------|---------|--------------|-------------------|
| $\kappa_1$ (naiv mod-12) | **0,0087** | ja | BE |
| $\kappa_2$ ($\nu_2$-Rotation) | **0,0041** | nein | EAEAA |
| $\kappa_3$ (Successor-Shift) | **0,0088** | ja | BE |

$\kappa_1 \approx \kappa_3$ zeigt **Kanaläquivalenz**, keine universelle Robustheit. $\kappa_2$ ist der Stresstest: zerstört BE, hält $R(10)$ sogar **unter** $\kappa_1$.

---

## Drei Leitfragen (Stufe 3)

### Frage A — Was ist κ-invariant unter „vernünftigen“ Kodierungen?

**Nicht invariant (widerlegt oder κ-sensitiv):**

- Konkrete verbotene Wörter ($\mathrm{BE}$, minimale Gegenbeispiele)
- Exakte $F_n$-Listen und $|F_n|$ bei festem $n$
- Hero-Zeugen wie BE als arithmetischer Invariant

**Möglicherweise invariant (heuristisch, offen):**

- **Wachstum von $R(k)$:** $R(k)$ klein und fallend mit $k$ (Conjecture: $R(k)\to 0$)
- **Existenz vieler verbotener Wörter:** $|F_n|$ wächst schnell mit $n$ (bei jeder getesteten $\kappa$)
- **Entropie-Ordnung:** $h_F$ in ähnlicher Größenordnung ($\approx 1{,}19$ bei $n\leq 8$, Schätzer)

| Eigenschaft | $\kappa_1$ | $\kappa_2$ | $\kappa_3$ | Invariant? |
|-------------|-----------|-----------|-----------|------------|
| BE verboten | ja | **nein** | ja | **nein** |
| $R(10)$ klein ($<1{,}2\,\%$) | ja | ja | ja | **experimentell ja** |
| $|F_n|$ groß für $n\leq 8$ | ja | ja (andere Liste) | ja | **qualitativ ja** |
| $h_F \approx 1{,}19$ | ja | ja | — | **heuristisch** |

**Operative Aufgabe:** Für jede neue κ-Variante systematisch prüfen, welche Spalte der Tabelle sich ändert — nicht nur BE.

---

### Frage B — Wann sind zwei Kodierungen äquivalent?

**Definition (vorgeschlagen, Conjecture-Ebene):**

$$\kappa_i \sim \kappa_j \quad:\Longleftrightarrow\quad R_i(k) \text{ und } R_j(k) \text{ stimmen asymptotisch überein}$$

(d.h. $\displaystyle\lim_{k\to\infty} R_i(k)/R_j(k) = 1$ oder beide $\to 0$ mit vergleichbarer Rate — präzise Form noch offen).

**Hypothese (experimentell gestützt, kein Theorem):**

- $\kappa_1 \sim \kappa_3$ — gleicher Informationskanal (Successor-Shift)
- $\kappa_2 \not\sim \kappa_1$ — andere Feinstruktur ($\nu_2$-Rotation), andere $F_n$-Listen, aber ähnliche Dünnheit bei großem $k$

**Klassifikationskriterien** (aus Stufe 2B, weiterzuführen):

| Kriterium | Frage |
|-----------|-------|
| dynamiktreu? | Präfixe = Bahn? |
| injektiv? | Informationsverlust? |
| shiftäquivalent? | Gleicher Kanal wie $\kappa_1$? |
| rotationssensitiv? | Reagiert auf $\nu_2$ wie $\kappa_2$? |
| stabilisiert $R(k)$? | Dünnheit bei wachsendem $k$? |

Ziel: Äquivalenzklassen von Kodierungen statt Einzeltests.

---

### Frage C — Existiert universelles $\mathcal{L}_{\mathrm{arith}}^*$?

**Definition (natürliche Tao-Frage):**

$$\mathcal{L}_{\mathrm{arith}}^* \;:=\; \bigcap_{\kappa \in \mathcal{K}_{\mathrm{reason}}} \mathcal{L}_{\mathrm{arith}}^\kappa$$

wobei $\mathcal{K}_{\mathrm{reason}}$ die Klasse „vernünftiger“ dynamiktreuer Kodierungen ist (noch zu präzisieren).

**Interpretation:**

- Wörter in $\mathcal{L}_{\mathrm{arith}}^*$ wären **kodierungsunabhängig realisierbar** — echte arithmetische Invarianten
- Wörter in $\mathcal{L} \setminus \mathcal{L}_{\mathrm{arith}}^*$ sind κ-abhängig verboten oder nur unter Spezialkodierungen realisierbar
- BE liegt **nicht** in $\mathcal{L}_{\mathrm{arith}}^*$ (widerlegt durch $\kappa_2$)

**Status:** **zentral offen.** Kein Algorithmus, keine Charakterisierung. Verbindung zu `FaithfulKappa` (Lean) und kodierungsfreier Definition von $\mathcal{L}_{\mathrm{arith}}$.

---

## Neue Hierarchie offener Fragen

Ersetzt die frühere Priorität (Attraktor, Geometrie, Präzession, κ als Hauptangriffe).

| Rang | Frage | Status | Artefakt |
|------|-------|--------|----------|
| **1** | **Kodierungsfreie Definition von $\mathcal{L}_{\mathrm{arith}}$** | zentral offen | Frage C; `FaithfulKappa` |
| **2** | **Verhalten von $R(k)$** | Conjecture: $R(k)\to 0$ | `collatz_kappa_robustheit.py`, `collatz_forbidden_words.py` |
| **3** | **Entropie $h_F$** | Schätzer $\approx 1{,}19$; κ-robust? | Stufe 2B Test 3 |
| **4** | **Dynamische Konsequenzen** | Pipeline-Ende: $F_n \to$ Trajektorien | `collatz_equivalenz_e_infty.tex` |
| **5** | **Lemma E** (Präperiodizität) | TeX-Skizze, nicht Lean | späterer Endpunkt, nicht Stufe-3-Kern |

**Bewusst zurückgestellt:** Attraktor $E$, Geometrie ($\Phi_{\mathrm{pref}}$, Präzession), treue $\kappa$ als **Hauptangriff** — epistemisch abgegrenzt oder in Stufe 1/2B bearbeitet.

---

## Living spreadsheet (Stufe 3)

| Zeile | Objekt | Ebene | Stand Juni 2026 |
|-------|--------|-------|-----------------|
| 1 | $\mathcal{L}_{\mathrm{arith}} \subsetneq \mathcal{L}$ | Experiment | **gesichert** (PR #39) |
| 2 | BE kodierungsunabhängig | Negativtest | **widerlegt** (PR #40, $\kappa_2$) |
| 3 | $R(10)$ klein für $\kappa_1,\kappa_2,\kappa_3$ | Experiment | **gestützt** |
| 4 | $R(k)\to 0$ | Conjecture | **offen** |
| 5 | $\kappa_1 \sim \kappa_3$ | Experiment | **Kanaläquivalenz** |
| 6 | $\kappa_2 \not\sim \kappa_1$ | Experiment | **gestützt** |
| 7 | $\mathcal{L}_{\mathrm{arith}}^*$ | Definition | **offen** |
| 8 | $h_F$ κ-robust | Conjecture | **heuristisch** |
| 9 | Dynamik: $F_n \to$ Collatz | Brücke | **fehlt** |
| 10 | Lemma E | Conjecture/Skizze | **offen** |

---

## Konkrete nächste Schritte (kein Collatz-Angriff)

1. **κ-Klassifikation:** Weitere Varianten nach den sechs Kriterien (Stufe 2B); Äquivalenzklassen kartieren.
2. **$R(k)$-Kurve:** $k=4,6,8,10,12,\ldots$ für alle drei $\kappa$; Rate des Abfalls schätzen.
3. **$\mathcal{L}_{\mathrm{arith}}^*$ approximieren:** Schnitt über dokumentierte $\kappa$; Kandidaten für kodierungsunabhängige Wörter.
4. **$h_F$ verfeinern:** Mehr $n$, Vergleich $\kappa_1$ vs. $\kappa_2$; asymptotische Schätzung — kein Theorem erwarten.
5. **Dynamische Brücke skizzieren:** Welche κ-invarianten Eigenschaften könnten Trajektorien einschränken? — nur Prop-Ebene, kein Beweis.

---

## Ehrliche Grenzen

1. Drei $\kappa$ sind **Experimente**, keine Lösung von `faithfulKappaExists`.
2. „Vernünftige“ Kodierung ist **noch nicht** formal definiert.
3. Asymptotische Äquivalenz $\kappa_i \sim \kappa_j$ ist eine **Arbeitsdefinition**, kein Theorem.
4. $\mathcal{L}_{\mathrm{arith}}^*$ kann leer oder trivial sein — muss getestet werden.
5. Kleine $R(k)$ bei endlichem $k$ impliziert **nicht** $R(k)\to 0$.
6. Dynamische Konsequenzen und Lemma E sind **nicht** Ziel von Stufe 3 — nur vorbereitende Hierarchie.
7. **Kein** Collatz-Beweis in dieser Stufe.

---

## Artefakte und Verweise

| Datei | Rolle |
|-------|-------|
| `collatz_stufe3_kodierungsinvarianz.md` | dieses Dokument |
| `collatz_stufe2b_kappa_robustheit.md` | Stufe 2B-Bericht (PR #40) |
| `collatz_kappa_robustheit.py` / `.json` | Drei-κ-Experimente |
| `collatz_forbidden_words.py` / `.json` | $F_n$-Katalog ($\kappa_1$) |
| `collatz_generalangriff_2026.md` | Forschungsreport, Stufe-3-Verweis |
| `collatz_offene_punkte.md` | Synthese, revidierte Prioritäten |
| `CollatzEabc.Kappa.lean` | `FaithfulKappa`, `kappaConjecture` |
