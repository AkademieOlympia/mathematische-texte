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

1. **Stufe 1 ($\kappa$):** Treue Kodierung als explizites Lean-Objekt (`FaithfulKappa`, `kappaConjecture`);
   naive $\kappa$ dokumentiert und verworfen — nicht stillschweigend ersetzt.
2. **Stufe 2 ($L_{\mathrm{arith}}$):** Realisierbarkeit als getrennte Sprache neben Grammatik
   (`CollatzEabc.ArithLanguage`, `collatz_l_arith_test.py`).
3. **Stufe 3 (Lemma E):** Präperiodizität nur nach geklärter Definition von $E_{N,N}$ und
   Realisierbarkeit — keine Bedeutungsverschiebung zwischen $E_{\mathrm{diag}}$ und $E_\infty$.

**PR-/Test-Kultur:** Jeder Zweig liefert Zeugen (positive Evidenz), Negativtests (verworfene
Pfade) und Lean-Schnittstellen (`Prop` mit klarer Beweisstrategie) — vergleichbar mit
kontinuierlicher Integrität in großen Formalisierungsprojekten.

---

## Referenzen

- Tao (2019) — logarithmische Dichte für „almost bounded“ Collatz-Bahnen (nicht punktweise Konvergenz)
- `collatz_generalangriff_2026.md` — Forschungsplan Stufe 1–3
- `collatz_offene_punkte.md` — Synthese offener Punkte und Negativresultate
- `CollatzEabc.Kappa`, `CollatzEabc.ArithLanguage` — Lean-Schnittstellen
