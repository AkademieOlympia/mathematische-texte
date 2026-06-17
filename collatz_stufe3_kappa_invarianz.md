# **Stufe 3: Kodierungsinvarianz**

**Branch:** `collatz/kappa-invarianz-stufe3` · **Stand:** Juni 2026  
**Vorgänger:** PR #39 (Stufe 2), PR #40 (Stufe 2B) — beide auf `main` (Commit `50dea21`)  
**Kein Collatz-Beweis.** Dieses Dokument definiert die nächste Forschungsstufe; es greift Collatz nicht direkt an.

**Methodik:** Tao-Stil (IEANTN/PNT+-Parallele) — Definition / Zeuge / Experiment / Theorem / Conjecture strikt trennen (`collatz_formalisierung_tao_stil.md`).

---

## Kernthese

Stufe 3 untersucht **nicht mehr einzelne verbotene Wörter**, sondern **kodierungsstabile Größen**: Sprachverdünnung, Entropie und Äquivalenzklassen von κ-Kodierungen.

> **Boxed (Leitsatz):**  
> *Nicht BE ist robust, sondern möglicherweise die Verdünnung.*

---

## Leitkette

Die Forschungslogik von Stufe 3 folgt einer festen Kette — von der Kodierungsfamilie bis zur dynamischen Brücke:

$$\boxed{\kappa\text{-Familie} \;\to\; \text{Äquivalenzklassen} \;\to\; R_\kappa(k) \;\to\; h_\kappa \;\to\; \mathcal{L}_{\mathrm{arith}}^* \;\to\; \text{Dynamik}}$$

| Stufe | Objekt | Ebene | Aufgabe |
|-------|--------|-------|---------|
| 1 | κ-Familie | Definition / Experiment | Zulässige Kodierungen kartieren |
| 2 | Äquivalenzklassen | Conjecture | $\kappa_i \sim \kappa_j$ nach asymptotischem $R(k)$ |
| 3 | $R_\kappa(k)$ | Experiment / Conjecture | Dünnheitskurven vergleichen |
| 4 | $h_\kappa$ | Experiment / Conjecture | Entropie der verbotenen Sprache |
| 5 | $\mathcal{L}_{\mathrm{arith}}^*$ | Definition | Schnitt über vernünftige κ — **noch nicht berechnen** |
| 6 | Dynamik | Brücke (offen) | Welche Invarianten schränken Trajektorien ein? |

**Wichtig:** $\mathcal{L}_{\mathrm{arith}}^*$ wird hier **nur definiert**, nicht approximiert oder algorithmisch berechnet. Erst die fünf definitorischen Fragen (unten) müssen geklärt sein.

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

## Fünf definitorische Fragen vor $\mathcal{L}_{\mathrm{arith}}^*$

Bevor der Schnitt $\mathcal{L}_{\mathrm{arith}}^* = \bigcap_{\kappa \in \mathcal{K}} \mathcal{L}_{\mathrm{arith}}^\kappa$ berechnet oder approximiert wird, müssen diese Fragen sauber beantwortet werden:

### Frage 1 — Was ist eine „zulässige" κ-Kodierung?

**Definition (Arbeitsentwurf):** Eine Kodierung $\kappa:\mathbb{Z}_2 \to \{E,A,B,C\}^{<\omega}$ (bzw. $\kappa_K$ für festes $K$) ist *zulässig*, wenn sie mindestens erfüllt:

| Kriterium | Bedeutung | Lean-Bezug |
|-----------|-----------|------------|
| **dynamiktreu** | Präfixe von $\kappa(U^k(n))$ = EABC-Blockschritte der Bahn | `kappaPrefix_get_shift` |
| **vollständig auf Testmenge** | Für ungerade Starts mit mod-$12 \in \{1,5,7,11\}$ liefert $\kappa_K$ ein vollständiges Wort | Experiment |
| **kein trivialer Informationsverlust** | Nicht konstant; nicht nur Länge-1-Wörter | heuristisch |

**Offen:** Ist Injektivität erforderlich? Reicht Dynamiktreue? Verbindung zu `FaithfulKappa` in `CollatzEabc.Kappa.lean`.

**Epistemisches Label:** **Definition** (unvollständig) — noch keine formale Klasse $\mathcal{K}_{\mathrm{reason}}$.

---

### Frage 2 — Wann sind zwei κ äquivalent?

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

**Epistemisches Label:** **Conjecture** — Äquivalenzklassen statt Einzeltests.

---

### Frage 3 — Welche Größen sind invariant?

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

**Epistemisches Label:** **Experiment** (für $R(10)$) / **Conjecture** (für asymptotisches Verhalten).

---

### Frage 4 — Welche Größen sind Artefakte?

**Artefakte** (κ-abhängig, nicht als Invariant interpretieren):

| Artefakt | Beispiel | Befund |
|----------|----------|--------|
| Konkretes verbotenes Wort | BE bei $\kappa_1$, EAEAA bei $\kappa_2$ | PR #40 |
| Exakte $F_n$-Liste | unterschiedliche Listen für $\kappa_1$ vs. $\kappa_2$ | PR #40 |
| Hero-Zeuge BE | „arithmetisch verboten" | **widerlegt** |
| $|F_n|$ bei festem $n$ | numerisch verschieden | κ-sensitiv |

**Nicht-Artefakt (Kandidat):** die **qualitative** Verdünnung — $R(k)$ klein und möglicherweise fallend.

**Epistemisches Label:** **Negativtest** (für BE) / **Experiment** (für $R$-Robustheit).

---

### Frage 5 — Welche Rolle spielt Dynamiktreue?

Dynamiktreue ist die **Mindestbedingung** für jede zulässige κ: ohne sie ist $\mathcal{L}_{\mathrm{arith}}^\kappa$ eine willkürliche Teilmenge der Grammatik, nicht eine arithmetische Realisierbarkeitssprache.

| Eigenschaft | Ohne Dynamiktreue | Mit Dynamiktreue |
|-------------|-------------------|-------------------|
| $\mathcal{L}_{\mathrm{arith}}^\kappa$ | beliebig manipulierbar | an Collatz-Bahn gebunden |
| $R(k)$ | ohne Bedeutung | misst echte Lücke |
| $\mathcal{L}_{\mathrm{arith}}^*$ | leer oder trivial | sinnvoller Schnitt |

**Lean:** Naive mod-12-$\kappa$ ist dynamiktreu (`kappaPrefix_get_shift`), aber nicht injektiv — reicht Dynamiktreue allein?

**Epistemisches Label:** **Definition** (notwendige Bedingung) — Verhältnis zu Injektivität **offen**.

---

## Drei Leitfragen (Stufe 3)

### Frage A — Was ist κ-invariant unter „vernünftigen" Kodierungen?

Siehe Frage 3 oben. **Operative Aufgabe:** Für jede neue κ-Variante systematisch prüfen, welche Spalte der Invarianten-Tabelle sich ändert — nicht nur BE.

**Priorität in der Hierarchie:** Rang 2 ($R(k)$) und Rang 3 ($h_F$).

---

### Frage B — Wann sind zwei Kodierungen äquivalent?

Siehe Frage 2 oben. **Ziel:** Äquivalenzklassen von Kodierungen statt Einzeltests.

**Experimenteller Befund:** $\kappa_1 \approx \kappa_3$ (Kanaläquivalenz); $\kappa_2$ in anderer Klasse.

---

### Frage C — Existiert universelles $\mathcal{L}_{\mathrm{arith}}^*$?

**Definition (natürliche Tao-Frage — noch nicht berechnen):**

$$\mathcal{L}_{\mathrm{arith}}^* \;:=\; \bigcap_{\kappa \in \mathcal{K}_{\mathrm{reason}}} \mathcal{L}_{\mathrm{arith}}^\kappa$$

wobei $\mathcal{K}_{\mathrm{reason}}$ die Klasse „vernünftiger" dynamiktreuer Kodierungen ist (Frage 1).

**Interpretation:**

- Wörter in $\mathcal{L}_{\mathrm{arith}}^*$ wären **kodierungsunabhängig realisierbar** — echte arithmetische Invarianten
- Wörter in $\mathcal{L} \setminus \mathcal{L}_{\mathrm{arith}}^*$ sind κ-abhängig verboten oder nur unter Spezialkodierungen realisierbar
- BE liegt **nicht** in $\mathcal{L}_{\mathrm{arith}}^*$ (widerlegt durch $\kappa_2$)

**Status:** **zentral offen.** Kein Algorithmus, keine Charakterisierung, **keine Berechnung in Stufe 3**. Verbindung zu `FaithfulKappa` (Lean) und kodierungsfreier Definition von $\mathcal{L}_{\mathrm{arith}}$.

**Epistemisches Label:** **Definition** — Existenz und Nicht-Trivialität sind offene Fragen.

---

## Prioritätshierarchie (Stufe 3)

Ersetzt die frühere Priorität (Attraktor, Geometrie, Präzession, κ als Hauptangriffe).

| Rang | Frage | Status | Artefakt |
|------|-------|--------|----------|
| **1** | **Kodierungsfreie Definition von $\mathcal{L}_{\mathrm{arith}}$** | zentral offen | Frage C; `FaithfulKappa` |
| **2** | **Verhalten von $R(k)$** | Conjecture: $R(k)\to 0$ | `collatz_kappa_robustheit.py`, `collatz_forbidden_words.py` |
| **3** | **Entropie $h_F$** | Schätzer $\approx 1{,}19$; κ-robust? | Stufe 2B Test 3 |
| **4** | **Dynamische Konsequenzen** | Pipeline-Ende: $F_n \to$ Trajektorien | `collatz_equivalenz_e_infty.tex` |
| **5** | **Lemma E** (Präperiodizität) | TeX-Skizze, nicht Lean | späterer Endpunkt, nicht Stufe-3-Kern |

**Bewusst zurückgestellt:** Attraktor $E$, Geometrie ($\Phi_{\mathrm{pref}}$, Präzession), treue $\kappa$ als **Hauptangriff** — epistemisch abgegrenzt oder in Stufe 1/2B bearbeitet.

### Geometrischer Nebenzweig: Morley-Operator (nicht Metrik)

Parallel zur arithmetischen Stufe-3-Pipeline existiert ein **definitorischer** geometrischer Rahmen in `collatz_morley_metrik_erweiterung.md` (§ **Geometrie zweiter Ordnung**): Grundobjekt ist das Dreieck $(A,B,C)$ statt $(p,v)$; zentraler Akteur $T_M : \mathcal{T}(M)\to\mathcal{T}(M)$ mit Morley-Fluss $\Delta_{n+1}=T_M(\Delta_n)$. Sensoren $F_M$, $S_M$, $K_M$, $\mu_M$ messen Zustände; EABC: $E=T_M$, nicht Punkt $E$. Offene Hauptfrage: Invarianten von $T_M$ auf triangulierten Flächen. Ikosaeder–Dodekaeder-Dualität (`eabc_icosahedron_test.py`) koppelt Triangulierung und Renormierung. Verzahnung mit $\kappa$-Robustheit **spekulativ** — Rang 6+, kein Stufe-3-Kern.

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
| 7 | $\mathcal{L}_{\mathrm{arith}}^*$ definiert | Definition | **offen** (nicht berechnet) |
| 8 | $h_F$ κ-robust | Conjecture | **heuristisch** |
| 9 | Dynamik: $F_n \to$ Collatz | Brücke | **fehlt** |
| 10 | Lemma E | Conjecture/Skizze | **offen** |
| 11 | Verdünnung κ-robust (nicht BE) | Conjecture | **Leitsatz** |

---

## Konkrete nächste Schritte (kein Collatz-Angriff)

1. **κ-Klassifikation:** Weitere Varianten nach den sechs Kriterien (Stufe 2B); Äquivalenzklassen kartieren.
2. **$R(k)$-Kurve:** $k=4,6,8,10,12,\ldots$ für alle drei $\kappa$; Rate des Abfalls schätzen.
3. **$\mathcal{K}_{\mathrm{reason}}$ präzisieren:** Formale Kriterien für zulässige κ — Voraussetzung für $\mathcal{L}_{\mathrm{arith}}^*$.
4. **$h_F$ verfeinern:** Mehr $n$, Vergleich $\kappa_1$ vs. $\kappa_2$; asymptotische Schätzung — kein Theorem erwarten.
5. **Dynamische Brücke skizzieren:** Welche κ-invarianten Eigenschaften könnten Trajektorien einschränken? — nur Prop-Ebene, kein Beweis.

**Bewusst nicht in Stufe 3:** Berechnung oder Approximation von $\mathcal{L}_{\mathrm{arith}}^*$.

---

## Ehrliche Grenzen

1. Drei $\kappa$ sind **Experimente**, keine Lösung von `faithfulKappaExists`.
2. „Vernünftige" Kodierung ist **noch nicht** formal definiert.
3. Asymptotische Äquivalenz $\kappa_i \sim \kappa_j$ ist eine **Arbeitsdefinition**, kein Theorem.
4. $\mathcal{L}_{\mathrm{arith}}^*$ kann leer oder trivial sein — muss getestet werden, aber **erst nach** Klärung der fünf definitorischen Fragen.
5. Kleine $R(k)$ bei endlichem $k$ impliziert **nicht** $R(k)\to 0$.
6. Dynamische Konsequenzen und Lemma E sind **nicht** Ziel von Stufe 3 — nur vorbereitende Hierarchie.
7. **Kein** Collatz-Beweis in dieser Stufe.
8. Der Leitsatz „Verdünnung, nicht BE" ist eine **Conjecture**, kein Theorem.

---

## Artefakte und Verweise

| Datei | Rolle |
|-------|-------|
| `collatz_stufe3_kappa_invarianz.md` | **kanonisches** Stufe-3-Dokument (dieses) |
| `collatz_stufe3_kodierungsinvarianz.md` | Redirect → dieses Dokument |
| `collatz_stufe2b_kappa_robustheit.md` | Stufe 2B-Bericht (PR #40) |
| `collatz_kappa_robustheit.py` / `.json` | Drei-κ-Experimente |
| `collatz_forbidden_words.py` / `.json` | $F_n$-Katalog ($\kappa_1$) |
| `collatz_generalangriff_2026.md` | Forschungsreport, Stufe-3-Verweis |
| `collatz_morley_metrik_erweiterung.md` | Morley-Sensorik ($F_M$, $S_M$, $K_M$, $\mu_M$, $T_M$); Brücke, nicht Hauptangriff |
| `collatz_offene_punkte.md` | Synthese, revidierte Prioritäten |
| `CollatzEabc.Kappa.lean` | `FaithfulKappa`, `kappaConjecture` |
