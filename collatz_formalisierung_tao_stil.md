# Methodik: Lean als Wahrheitsfilter

**Stand:** Juni 2026 · **Kontext:** EABC/Collatz-Generalangriff (`collatz_generalangriff_2026.md`,
PR #38–#39). **Epistemische Warnung:** Kein Collatz-Beweisanspruch — methodische Einordnung
der Formalisierungsarbeit.

---

## Tao-Stil: Lean nicht als Orakel

Terence Tao hat in der Formalisierung der Primzahlzahlensatz-Kette (PNT+) Lean nicht als
**Entdeckungsmaschine** oder Beweisorakel eingesetzt, sondern als **Buchhaltungsmaschine**
für fragile Mathematik: explizite Konstanten, Fehlerterme, Kettenabhängigkeiten und
jede implizite Annahme werden sichtbar gemacht, bevor sie in den nächsten Schritt übergeht.

Diese Haltung überträgt sich direkt auf EABC/Collatz — mit anderen Objekten, aber derselben
Disziplin: Definitionen, Zeugen, Heuristiken und Theoreme dürfen nicht unbemerkt ineinander
übergehen.

| Tao / PNT+ | EABC / Collatz |
|------------|----------------|
| explizite Fehlerterme | explizite Zeugenmengen ($E_\infty$, $E_{\mathrm{diag}}$, $E_{N,N}$) |
| abhängige Konstanten | abhängige Definitionen ($\kappa$, $\mathcal{L}$, $\mathcal{L}_{\mathrm{arith}}$) |
| PNT-Ketten | $\kappa \to L_{\mathrm{arith}} \to$ Lemma E |
| Formalisierung gegen Konstantenfehler | Formalisierung gegen Bedeutungsverschiebung |
| CI-Kultur (kontinuierliche Integrität) | PR-/Zeugen-/Negativtest-Kultur |

> **Boxed:** EABC braucht keine KI als Entdeckungsmaschine, sondern Lean als **Wahrheitsfilter**.
>
> **Wert der Formalisierung:** Sie verhindert, dass Definitionen, Zeugen, Heuristiken und
> Theoreme unbemerkt ineinander übergehen — und schneidet falsche Pfade ab, bevor sie
> strategische Ressourcen binden.

---

## Heutiger Befund im Tao-Stil (Juni 2026)

Die jüngsten Negativ- und Grenzresultate sind keine Niederlagen, sondern **Buchhaltung**:

| Pfad | Befund | Tao-Stil-Lesart |
|------|--------|-----------------|
| Präzession $I(Q)\to\Pi(Q)$ | **negativ** ($r\approx -0{,}18$) | Heuristik als Hauptangriff verworfen |
| $\Phi_{\mathrm{pref}}$ Länge 4 | **trivial** (nur ABCE/CEAB) | Diskriminante ohne dynamische Brücke unzureichend |
| Naive $\kappa_K$ | dynamiktreu, **nicht injektiv** | Brückenskizze, keine treue Kodierung (`FaithfulKappa` offen) |
| $\mathcal{L}_{\mathrm{arith}}\subsetneq\mathcal{L}$ | $w=\mathrm{BE}$ Gegenbeispiel | Grammatik ⊄ Realisierbarkeit operationalisiert (PR #39) |

Jeder Eintrag trennt sauber: *Was ist bewiesen? Was ist nur numerisch? Was ist die falsche
Zielsetzung?* — analog zu Taos expliziten Fehlertermen in der PNT-Kette.

---

## Operative Konsequenz für den Generalangriff

1. **Stufe 2 ($L_{\mathrm{arith}}$):** Realisierbarkeit als getrennte Sprache neben Grammatik
   (`CollatzEabc.ArithLanguage`, `collatz_l_arith_test.py`) — **stärkster operativer Zweig**
   (PR #39; Ratio $k=10$: $\approx 0{,}0087$).
2. **Stufe 1 ($\kappa$):** Treue Kodierung als explizites Lean-Objekt (`FaithfulKappa`, `kappaConjecture`);
   naive $\kappa$ dokumentiert und verworfen — soll verborgene Realisierbarkeitsregeln sichtbar machen.
3. **Stufe 3 (Lemma E):** Präperiodizität nur nach geklärter Definition von $E_{N,N}$ und
   Realisierbarkeit — keine Bedeutungsverschiebung zwischen $E_{\mathrm{diag}}$ und $E_\infty$.

**PR-/Test-Kultur:** Jeder Zweig liefert Zeugen (positive Evidenz), Negativtests (verworfene
Pfade) und Lean-Schnittstellen (`Prop` mit klarer Beweisstrategie) — vergleichbar mit
kontinuierlicher Integrität in großen Formalisierungsprojekten.

---

## Methodik 2 — Lean als Wahrheitsfilter (fünf Ebenen)

Die Stufe-2-Arbeit ($\mathcal{L}_{\mathrm{arith}}\subsetneq\mathcal{L}$, PR #39) ist kein
philosophischer Anhang, sondern **methodischer Rahmen** für die Interpretation numerischer
und formaler Befunde. Analog zu Taos PNT+-Kette werden Definition, Zeuge, Experiment, Theorem
und Vermutung strikt getrennt — bevor daraus strategische Prioritäten abgeleitet werden.

### Die fünf Ebenen

| Ebene | Inhalt (Stufe 2) | Status (Juni 2026) |
|-------|------------------|-------------------|
| **1 — Definition** | Alphabet $\{E,A,B,C\}$; Grammatik $\mathcal{L}(k)$; Verbot $BB$; endliche $C$-Ketten mit Kapazität $\mathrm{cap}\geq 1$; Zwang $C^*_{\max}\to EA$; arithmetische Realisierbarkeit über $\kappa$-Präfix | **formalisiert** (`CollatzEabc.ArithLanguage`, `collatz_l_arith_test.py`) |
| **2 — Zeuge** | Konkrete nicht-realisierbare Wörter: $w=\mathrm{BE}$ (minimal, Länge 2); weitere Gegenbeispiele aus Vollliste | **Beobachtung** — kein universelles Theorem |
| **3 — Experiment** | Ratios $|L_{\mathrm{arith}}|/|L|$; Häufigkeiten; Volllisten bzw. Stichproben | **durchgeführt** ($k=10$: $24\,818$ Realisierungen über naive $\kappa$, $n\leq 10^6$) |
| **4 — Theorem** | Bewiesene Ausschlüsse innerhalb $\mathcal{L}$: $BB\notin\mathcal{L}$; endliche $C$-Kapazität; $C^*_{\max}\to EA$ | **bewiesen** (lokale Grammatik, Lean/TeX) |
| **5 — Conjecture** | $\displaystyle\lim_{k\to\infty}\frac{|L_{\mathrm{arith}}(k)|}{|L(k)|}$ — asymptotisches Verhältnis Grammatik zu Realisierbarkeit | **offen** |

> **Epistemische Regel:** Eine Ebene darf die nächsthöhere nicht ersetzen. Ein Zeuge ($\mathrm{BE}$)
> ist kein Theorem; ein Experiment ($k=10$) ist keine asymptotische Aussage; ein Theorem über
> $\mathcal{L}$ impliziert nichts über $\mathcal{L}_{\mathrm{arith}}$.

### Kernbefund Stufe 2 ($k=10$)

$$\boxed{|L(10)| = 2\,860\,558;\quad |L_{\mathrm{arith}}(10)| = 24\,818;\quad
\frac{|L_{\mathrm{arith}}|}{|L|}\approx 0{,}0087\;\;(99{,}13\,\%\text{ ausgeschlossen})}$$

**Hauptbefund:** $L_{\mathrm{arith}}\ll L$ — strukturell, nicht nur punktuell. Die arithmetische
Sprache ist auf endlicher Suchtiefe ein **kleiner Anteil** der formalen Grammatik. Das Wort
$w=\mathrm{BE}$ ist der **erste sichtbare Zeuge** für $\mathcal{L}_{\mathrm{arith}}\subsetneq\mathcal{L}$,
nicht der Beweis einer allgemeinen Aussage.

**Tao-Interpretation:** Auf Experimentebene zeigt die Ratio, dass naive Grammatik-Zulässigkeit
massiv überzeichnet, was Collatz-Bahnen tatsächlich realisieren. Auf Zeugenebene liefert $\mathrm{BE}$
ein konkretes, überprüfbares Gegenbeispiel. Auf Theoremebene bleibt offen, ob $\mathrm{BE}\notin
\mathcal{L}_{\mathrm{arith}}$ für **alle** $n\in\mathbb{N}_{\mathrm{odd}}$ gilt — das wäre ein
Schritt von Beobachtung zu Aussage. Die Conjecture-Ebene fragt, ob das Verhältnis mit $k\to\infty$
gegen $0$ strebt oder eine positive asymptotische Dichte besitzt.

### Status-Tabelle (Stufe 2, ehrlich)

| Aspekt | Befund | Ebene |
|--------|--------|-------|
| Grammatik $\mathcal{L}(k)$ formal definiert | ja | Definition |
| $|L_{\mathrm{arith}}(10)|=24\,818$ über naive $\kappa$, $n\leq 10^6$ | ja | Experiment |
| $w=\mathrm{BE}$ tritt in der Suche nicht auf | ja | Beobachtung / Zeuge |
| $\mathrm{BE}\notin\mathcal{L}_{\mathrm{arith}}$ für alle ungeraden $n$ | **nein** | Theorem (offen) |
| Collatz / $E_\infty=\emptyset$ | **nein** | — |

### Strategische Einordnung

Nach dem Stufe-2-Befund ist $\mathcal{L}_{\mathrm{arith}}$ der **stärkste operative Zweig**
des Generalangriffs: Er operationalisiert die Brücke *Grammatik ⊄ Realisierbarkeit* mit
harten Zahlen, während treue $\kappa$ (Stufe 1) und Lemma E (Stufe 3) darauf aufbauen bzw.
später kommen. Die frühere Priorität $\kappa\to L_{\mathrm{arith}}\to$ Lemma E ist damit
**umgedreht**: zuerst die Realisierbarkeitslücke quantifizieren, dann Kodierung und
Präperiodizität anbinden.

---

## Referenzen

- Tao (2019) — logarithmische Dichte für „almost bounded“ Collatz-Bahnen (nicht punktweise Konvergenz)
- `collatz_generalangriff_2026.md` — Forschungsplan Stufe 1–3
- `collatz_offene_punkte.md` — Synthese offener Punkte und Negativresultate
- `CollatzEabc.Kappa`, `CollatzEabc.ArithLanguage` — Lean-Schnittstellen
