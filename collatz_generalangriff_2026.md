# Generalangriff Collatz (Juni 2026) — Forschungsreport

**Ziel:** Die kleinste noch fehlende Brücke **L** identifizieren, sodass
*(bewiesene EABC-/Lean-Struktur)* **+ L** ⟹ Collatz (odd-to-odd:
$\forall n$ ungerade $\exists K: U^K(n)=1$).

**Epistemische Warnung:** Dieses Dokument ist ein **Forschungsplan**, kein Collatz-Beweis.
Die Collatz-Vermutung bleibt offen.

---

## Ziel und Erfolgskriterium

> **Kleinste fehlende Brücke L:**  
> *(Bewiesene EABC-/Lean-Struktur)* **+ L** ⟹ Collatz.

**Erfolgskriterium:** L ist

1. **logisch schwächer** als volle Collatz-Vermutung (oder eine äquivalente, aber operationalisierbare Zerlegung),
2. **nicht** bereits durch lokale Grammatik, Dichte oder Markov-Heuristiken impliziert,
3. **Lean-formalisierbar** (mindestens als `Prop` mit klarer Beweisstrategie),
4. schließt die dokumentierte Lücke *endliche EABC-Grammatik ⇒ Ausschluss unendlich schlechter natürlicher Realisierungen*.

**Kernfeststellung:** Die bewiesene Struktur kontrolliert **Präfixe, Wörter und typische Bahnen** —
nicht den punktweisen Ausschluss von $E_\infty$. Jede minimal plausible L muss daher den Sprung von
*lokal zulässig* zu *global für jedes $n\in\mathbb{N}_{\mathrm{odd}}$ realisiert* leisten.

> **Boxed (Generalangriff):** Finde die **minimale Brücke** zwischen $E_\infty$ und $\mathbb{N}$.

**PR-Bedeutung (Generalangriff 2026):** EABC $\not\Rightarrow$ Collatz — aber
*(bewiesene EABC-/Lean-Struktur)* **+ kleine Brücke L** $\Rightarrow$ Collatz.

---

## Bekannter Stand

### Ebene A — bewiesen (arithmetischer Kern)

| Resultat | Beleg |
|----------|-------|
| mod-12 EABC-Klassifikation, ABCE/CEAB-Chiralität | `EABC.lean`, `Projektionszeuge.tex`, `CollatzEabc.Primvierling.*` |
| C-Ketten endlich, $B\to B$ verboten, $C^*_{\max}\to EA$ | `collatz_schlussartikel_arxiv.tex`, `CollatzEabc.Density` |
| LTE-Reset (gerades $N$), LTE-Worst-Familie | `CollatzEabc.Uniformity`, `collatz_uniformity.lean` |
| C-Ketten-Dichte $N/2^k$, Tail-Reihe $(1-p)^n\to 0$ | `collatz_density_appendix.lean` |

### Collatz-Struktur (Lean, sorry-frei bis Stufe D)

- `ExceptionSetInfinity` = $E_\infty$; Collatz $\Leftrightarrow$ $E_\infty=\emptyset$ (`CollatzEabc.Open`)
- `ExceptionSet` / `ExceptionSetDiag` = $E_{\mathrm{diag}}$ — **nicht** Collatz-äquivalent ($n=27$)
- Stufen A–D: Metrik, `ExceptionSetApprox`, Monotonie, $U$-Invarianz (`CollatzEabc.Z2Attraktor`)
- Stufe E: `collatzUniformityConjecture` — offen (`collatz_z2_attraktor.lean` Zeilen 431–467 mit `sorry`)

### Negativresultate / No-Go (Forschungsauftrag E)

| Idee | Status |
|------|--------|
| $\mathrm{dist}_2(n,1)\geq c\cdot 2^{-\log n}$ | **widerlegt** (LTE $4\cdot 3^r-1$, `dist_to_one_not_uniform_bound`) |
| Drift $\Lambda\approx -0{,}830$, mod-12-Mischung | nur **stationär/fast-sicher**, nicht punktweise |
| Bernoulli-Normschale als Lyapunov | **No-Go** |
| Präzession $I(Q)\to\Pi(Q)$ | **negativ** ($r\approx -0{,}18$; `collatz_praezession_info.tex`) |
| $\Phi_{\mathrm{pref}}$ Länge 4 | **trivial** (nur ABCE/CEAB) |
| $E_{\mathrm{diag}}\cap\mathbb{N}=\emptyset$ | **falsche Zielsetzung** |

### Φ_pref-Diskriminantentest (Auftrag D, explorativ)

`collatz_phi_pref_test.json`: bei Integrationsstrom-Wörtern (nicht Länge 4) **zwei chirale Röhren**
$T_{\mathrm{ABCE}}/T_{\mathrm{CEAB}}$ vs. diffuse Zufallsdaten — **geometrische Evidenz**, kein Collatz-Beweis.
Dynamische Brücke $\Phi=\Phi_{\mathrm{pref}}\circ\kappa$ und Kodierung $\kappa:\mathbb{Z}_2\to\{E,A,B,C\}^{<\omega}$
bleiben **offen** (`collatz_kepler_gedankenexperiment.tex`).

---

## Forschungsaufträge A–E als Arbeitsplan

| Auftrag | Inhalt | Konkrete nächste Schritte |
|---------|--------|--------------------------|
| **A** | Minimales Zwischenlemma | Lemma E formalisieren; Verbindung $E_{N,N}$ ↔ Trajektorientail |
| **B** | Symbolische Dynamik $\mathcal{L}\subset\{E,A,B,C\}^{\mathbb{N}}$ | $\mathcal{L}_{\mathrm{fin}}$ aus bewiesenen Regeln axiomatisieren; Realisierungsabbildung $w(n)$ definieren |
| **C** | $E_\infty$ vs. Approximationen | Dreier-Trennung $E_{\mathrm{diag}}, E_{\mathrm{tail}}, E_\infty$ (`collatz_equivalenz_e_infty.tex`); Lean-Brücke `ExceptionSetInfinity` ≠ `ExceptionSetDiag` |
| **D** | $\Phi_{\mathrm{pref}}$ mit echtem Diskriminantentest | Wörter aus Bahn-Präfixen $\kappa_K(n)$, nicht nur Primvierling-Länge 4; Kriterium aus `collatz_kepler_gedankenexperiment.tex` |
| **E** | Falsifikation | Drift/Uniformität/Präzession/Bernoulli nicht weiter als Hauptangriff; negative Tests dokumentieren |

---

## Strategische Verfeinerung (Juni 2026)

### Lücken-Diagnose

Die verbleibende Lücke liegt **nicht** bei mod-12, Chiralität oder Geometrie (diese Ebenen sind
gesichert bzw. epistemisch abgegrenzt), sondern beim **Übergang strukturell → punktweise**:
von lokal zulässigen EABC-Wörtern und typischen Bahnen zu einer Aussage über **jedes**
$n\in\mathbb{N}$ (odd-to-odd).

### Top-3 neu eingordnet

| Kandidat | Vermutete Stärke | Stufe |
|----------|------------------|-------|
| **Treue Kodierung $\kappa$** (L₄) | vermutlich **fundamentale Brücke** | **1** |
| **Realisierbarkeit** $\mathcal{L}_{\mathrm{arith}}\subsetneq\mathcal{L}$ (L₂) | möglicherweise **schwächer** | **2** |
| **Präperiodizität** Lemma E (L₁) | vermutlich **nahe an Collatz** | **3** |

### Warum $\kappa$ zuerst?

Drei Abbildungen aus demselben $n$:

- $n \mapsto$ Parität $\mapsto$ EABC-Blockschritt,
- $n \mapsto \Phi_{\mathrm{pref}}(\kappa(\cdot))$,
- $n \mapsto \mathbb{Z}_2$ (2-adische Metrik).

**Frage:** Geht dabei Information verloren? Eine **treue** Kodierung $\kappa$ würde die Collatz-Frage
in eine **symbolische** Frage über EABC-Wörter übersetzen — und damit Auftrag B (Grammatik) mit
Auftrag D ($\Phi_{\mathrm{pref}}$) verbinden.

Konkrete Stufe-1-Prüfpunkte:

1. Ist $\kappa$ **injektiv** (auf relevanten Präfixen)?
2. Ist $\kappa$ **dynamiktreu** (Präfixe von $\kappa(U^k(n))$ = EABC-Blockschritte der Bahn)?
3. Wo genau tritt **Informationsverlust** auf ($n \to$ Parität $\to$ EABC vs. $n \to \mathbb{Z}_2$)?

### $\mathcal{L}_{\mathrm{arith}}\subsetneq\mathcal{L}$

Nicht jedes grammatisch zulässige Wort $w\in\mathcal{L}$ ist **arithmetisch realisierbar**
als Bahnwort $w(n)$ eines $n\in\mathbb{N}_{\mathrm{odd}}$. Schlechte **unendliche** Wörter können
formal in $\mathcal{L}^{\mathbb{N}}$ liegen, ohne dass ein natürliches $n$ sie realisiert —
das ist die operationale Form der Brücke (vgl. L₂, `collatz_equivalenz_e_infty.tex`).

### Generalangriff-Prioritäten (Stufe 1–3)

| Rang | Stufe | PR / Zweig | Kernfrage |
|------|-------|------------|-----------|
| **1** | Treue Kodierung $\kappa$ | PR #38 `collatz/kappa-encoding` | Injektivität, Dynamiktreue, Informationsverlust |
| **2** | $\mathcal{L}_{\mathrm{arith}}$-Realisierbarkeit | **PR #39** `collatz/l-arith-stufe2` | Welche $w\in\mathcal{L}$ sind Bahnwörter? |
| **3** | Präperiodizität (Lemma E) | offen | Endliche Beobachtungsschlechtigkeit $\Rightarrow$ Präperiodizität |

**Nutzer-Priorität (Juni 2026):** PR **#39** und treue $\kappa$ (Stufe 1) sind die beiden
wichtigsten Zweige des Generalangriffs; Stufe 2 operationalisiert die offene Brücke aus
`collatz_equivalenz_e_infty.tex` Frage (ii).

### Stufe 1 — Implementierung (Juni 2026)

- **Lean:** `CollatzEabc.Kappa` — `kappaPrefix`, `FaithfulKappa`, `kappaConjecture`; Theorem `naiveKappa_shift`.
- **Python:** `collatz_kappa_test.py` — bei $K=8$, $n\leq 5000$: viele Starts mit `none` (mod $12\notin\{1,5,7,11\}$),
  Kollisionen unter definierten Starts; **Dynamik-Shift** für naive $\kappa$ verifiziert.
- **Fazit:** Treue $\kappa$ ist **stärker** als naive Präfix-Kodierung; `faithfulKappaExists K` bleibt offen.

---

## Kandidaten-Lemmas

### L₁ — Lemma E (Präperiodizität) ⭐ strategischer Kandidat

**Exakte Aussage** (`collatz_equivalenz_e_infty.tex`, Zeilen 204–211):

Sei $(n_k)$ eine odd-to-odd-Bahn mit $n_{k+1}=U(n_k)$. Falls
$$\forall N\;\exists k\geq N:\quad n_k\in E_{N,N},$$
dann ist $(n_k)$ **präperiodisch**.

| | |
|---|---|
| **Stärke** | **Schwächer** als Collatz (allein); zusammen mit Zyklus-Ausschluss und Divergenz-Ausschluss potenziell äquivalent |
| **Bewiesen** | Hypothese, Skizze in TeX; **nicht** in Lean |
| **Fehlt** | Vollständiger Beweis; Klarstellung, ob divergente Bahnen die Prämisse überhaupt erfüllen (vermutlich **nein** → Lemma E deckt eher Zyklus/Präperiodizität ab, nicht Divergenz) |
| **Lean** | **mittel–hoch** — $E_{N,N}$=`ExceptionSetApprox N N`, `iterateU` vorhanden; Präperiodizität als `∃ p, q, n_{p+i}=n_{q+i}` |

**Strategischer Wert:** Reduziert „dauerhaft schlecht auf endlicher Skala“ auf Analyse **präperiodischer EABC-Wörter**;
nichttriviale Zyklen sind durch Produktformeln stark eingeschränkt (Literatur + TeX-Epilog).

---

### L₂ — Arithmetisches Realisierbarkeits-Lemma (symbolische Collatz-Form)

**Exakte Aussage** (`collatz_equivalenz_e_infty.tex`):

Sei $\mathcal{L}\subset\{E,A,B,C\}^*$ die von lokalen Regeln erzeugte formale Sprache ($BB\notin\mathcal{L}$,
endliche $C$-Ketten, $C^*_{\max}\to EA$, …). Für jedes ungerade $n$ sei $w(n)\in\{E,A,B,C\}^{\mathbb{N}}$
das Bahnwort. Dann:
$$\nexists n\in\mathbb{N}_{\mathrm{odd}}:\; w(n)\in\mathcal{L}^{\mathbb{N}}\text{ und }U^K(n)\neq 1\;\forall K.$$

| | |
|---|---|
| **Stärke** | **Äquivalent** zu Collatz / $E_\infty=\emptyset$ |
| **Bewiesen** | Lokale Grammatik-Regeln (Teilmenge von $\mathcal{L}$) |
| **Fehlt** | Schritt 2 der TeX-Trennung: *Welche grammatisch zulässigen unendlichen Wörter sind arithmetisch realisierbar?* |
| **Lean** | **mittel** — `EabcLetter`, `classOf` in `Primvierling.Chirality`; Bahnkodierung $w(n)$ noch nicht formalisiert |

**Nutzen:** Macht die offene Brücke explizit: Grammatik ⊄ Realisierbarkeit.

---

### L₃ — Divergenz-Ausschluss auf odd-to-odd-Bahnen

**Exakte Aussage:**
$$\forall n\in\mathbb{N}_{\mathrm{odd}},\;\{U^k(n):k\in\mathbb{N}\}\text{ ist endlich (beschränkt).}$$

| | |
|---|---|
| **Stärke** | **Strikt schwächer** als Collatz (Collatz ⇒ Beschränktheit; Umkehrung falsch ohne Zyklus-Ausschluss) |
| **Bewiesen** | Tao (2019): logarithmische Dichte 1 für „almost bounded“ — **nicht** punktweise |
| **Fehlt** | Punktweiser Beweis; keine Brücke aus mod-12-Drift |
| **Lean** | **niedrig** — Mathlib hat keinen punktweisen Birkhoff/Kingman für Collatz |

---

### L₄ — Existenz treuer Kodierung $\kappa$

**Exakte Aussage** (`collatz_kepler_gedankenexperiment.tex`):

Es existiert $\kappa:\mathbb{Z}_2\to\{E,A,B,C\}^{<\omega}$ (bzw. $\kappa_K$ für jedes $K$), sodass für ungerade $n$
die Präfixe von $\kappa(U^k(n))$ mit den EABC-Blockschritten der Collatz-Bahn übereinstimmen, und
$$n\in E_\infty \iff \Phi_{\mathrm{pref}}(\kappa(U^k(n)))\in\mathcal{E}_\infty\subset\mathcal{M}\;\forall k.$$

| | |
|---|---|
| **Stärke** | **Schwächer** als Collatz, wenn $\mathcal{E}_\infty$ geometrisch charakterisiert werden kann |
| **Bewiesen** | $\Phi_{\mathrm{pref}}$ auf Wörtern (`PrefProjection.lean`); Diskriminantentest explorativ positiv |
| **Fehlt** | Konstruktion von $\kappa$; dynamische, nicht-triviale Brücke |
| **Lean** | **niedrig–mittel** — $\kappa$ ist das fehlende Objekt |

---

### L₅ — Zyklus-Ausschluss (klassisch, zerlegbar)

**Exakte Aussage:** Der einzige odd-to-odd-Zyklus ist $\{1\}$.

| | |
|---|---|
| **Stärke** | **Schwächer** als Collatz (Divergenz bleibt offen) |
| **Bewiesen** | Ausschlüsse modulo $2^k$ für endliches $k$ |
| **Fehlt** | „Für alle $k$“ / alle $n$ |
| **Lean** | **mittel** — projektspezifisch, nicht in Mathlib |

---

### L₆ — Uniformität zu $E_{\mathrm{diag}}$ (Stufe E, verworfen als Hauptweg)

**Exakte Aussage:** $\mathrm{dist}_2(U^k(n),E_{\mathrm{diag}})\to 0$ für alle $n$.

| | |
|---|---|
| **Stärke** | Unklar vs. Collatz; Zielobjekt $E_{\mathrm{diag}}$ ist **falsches** Collatz-Äquivalent |
| **Status** | Nutzer explizit **nicht** als Hauptangriff; LTE widerlegt naive $\mathrm{dist}_2(\cdot,1)$-Variante |

---

## Empfehlung: welches L ist minimal/plausibel?

**Epistemisch ehrlich:** Die **logisch schwächste** Aussage, die zusammen mit der bewiesenen Struktur
Collatz impliziert, ist weiterhin eine punktweise Aussage über **natürliche** Starts — im Kern
**$E_\infty=\emptyset$** (äquivalent zu Collatz). Die lokale Grammatik liefert dafür **keine** Implikation;
jede scheinbar schwächere Formulierung ist meist nur **Umformulierung** derselben Lücke.

**Strategisch (Juni 2026, nach Verfeinerung):** Drei-Stufen-Angriff — siehe Abschnitt
*Strategische Verfeinerung*.

**Stufe 1 — $\kappa$ (L₄):** Fundamentale Brücke zwischen $\mathbb{Z}_2$, EABC-Wörtern und
$\Phi_{\mathrm{pref}}$. Ohne treue Kodierung bleibt die dynamische Frage symbolisch unzugänglich;
erste positive Evidenz bei langen Integrationsstrom-Wörtern (`collatz_phi_pref_test.json`).

**Stufe 2 — $\mathcal{L}_{\mathrm{arith}}\subsetneq\mathcal{L}$ (L₂):** Klarste operationale Form
der offenen Brücke für Auftrag B; boxed Ziel in `collatz_equivalenz_e_infty.tex`.

**Stufe 3 — Lemma E (L₁):** Zwischenlemma in `collatz_equivalenz_e_infty.tex` und
`collatz_schlussartikel_arxiv.tex` (Epilog) — trennt endliche Beobachtungsschlechtigkeit von
unendlicher Nichtkonvergenz **für den präperiodischen Anteil**; vermutlich **näher an Collatz**
als Stufe 1–2, daher zuletzt im Angriffsplan.

Zerlegung (unverändert gültig):
Collatz $\Leftrightarrow$ (keine Divergenz) $\land$ (kein nichttrivialer Zyklus) $\land$
(Lemma E für präperiodischen Rest).

---

## Top-3 Kandidaten-Lemmas (Kurzliste)

| Rang | Stufe | Lemma | Warum |
|------|-------|-------|-------|
| **1** | **1** | **L₄ Treue Kodierung $\kappa$** | Fundamentale Brücke $n\leftrightarrow$ EABC-Wort $\leftrightarrow\mathbb{Z}_2$; Informationsverlust-Frage; verbindet Auftrag D ($\Phi_{\mathrm{pref}}$) mit dynamischer Collatz-Frage |
| **2** | **2** | **L₂ Arithmetische Realisierbarkeit** | Präzise Trennung $\mathcal{L}_{\mathrm{arith}}\subsetneq\mathcal{L}$; schlechte unendliche Wörter formal, aber nicht in $\mathbb{N}$ realisierbar |
| **3** | **3** | **L₁ Lemma E** (Präperiodizität) | Im Repo als Zwischenlemma skizziert; reduziert globale Frage auf präperiodische EABC-Analyse; vermutlich nahe an Collatz; Lean-Grundgerüst vorhanden |

---

## Referenzen

- `collatz_equivalenz_e_infty.tex` — Lemma E, $E_\infty$ vs. $E_{\mathrm{diag}}$, Realisierbarkeits-Lemma
- `CollatzEabc.Open` / `ExceptionSetInfinity` — Lean-Äquivalenz Collatz $\Leftrightarrow$ $E_\infty=\emptyset$
- `collatz_offene_punkte.md` — Synthese offener Punkte, Negativresultate
- `collatz_kepler_gedankenexperiment.tex` — $\kappa$, $\Phi_{\mathrm{pref}}$, Diskriminantentest
- `collatz_schlussartikel_arxiv.tex` — Epilog, Uniformität, EABC-Struktur

---

## Stufe 1 — Implementierung (κ, Juni 2026)

**Lean:** `collatz_eabc_core/CollatzEabc/Kappa.lean`

| Objekt | Inhalt |
|--------|--------|
| `EabcWord` | `List EabcLetter` |
| `classOfLetter` | mod-12 $\to$ `Option EabcLetter` |
| `kappaPrefix n K` | erste $K$ Schritte als `List (Option EabcLetter)` |
| `FaithfulKappa K` | Schnittstelle: volle Wörter, Klassenübereinstimmung, Shift+Append, Injektivität |
| `kappaConjecture` | $\forall K>0$, existiert treue $\kappa$ — **offen** |
| `kappaPrefix_get_shift` | Dynamiktreue der **naiven** $\kappa$ (sorry-frei) |

**Python:** `collatz_kappa_test.py` → `collatz_kappa_test.json`

**TeX:** `collatz_kappa_encoding.tex`

### Ehrliche Testergebnisse ($N=10^5$, $K=4$)

- **Dynamiktreue:** ja (0 Fehler; deckt sich mit Lean `kappaPrefix_get_shift`)
- **Injektivität:** nein ($\sim 1{,}47\times 10^8$ Kollisionspaare)
- **Vollständigkeit:** $33333/50000$ Starts ohne $\bot$-Einträge; $\approx 8{,}3\%$ undefinierte Schritte ($n\equiv 3,9\pmod{12}$)

Die naive mod-12-$\kappa$ ist eine **Brückenskizze**, keine treue Kodierung im Sinne von `FaithfulKappa`.

---

## Stufe 2 — $\mathcal{L}_{\mathrm{arith}}$ (L₂, Juni 2026)

> **Boxed (Stufe 2):** Welche EABC-Wörter sind **tatsächlich arithmetisch realisierbar**?
> $\mathcal{L}_{\mathrm{arith}}\subsetneq\mathcal{L}$ könnte die erste echte Brücke zwischen
> Collatz-Arithmetik, symbolischer Dynamik und Lemma E sein.

**Priorität:** PR **#39** (dieser Zweig) ist neben treuer $\kappa$ (Stufe 1 / PR #38) der
wichtigste Generalangriff-Zweig — beide hängen zusammen: treue $\kappa$ soll verborgene
Verbotsregeln von $\mathcal{L}_{\mathrm{arith}}$ sichtbar machen.

### Verbindung zu PR #38 ($\kappa$-Negativresultat)

| Befund (Stufe 1) | Konsequenz für Stufe 2 |
|------------------|------------------------|
| Naive $\kappa_K$ **dynamiktreu** | Realisierbarkeitstest über `kappaPrefixWord` ist konsistent mit Bahn |
| Naive $\kappa_K$ **nicht injektiv** | $\kappa_{\mathrm{naiv}}$ reicht nicht als finale Brücke |
| $\approx 8{,}3\,\%$ undefinierte Schritte ($n\equiv 3,9\pmod{12}$) | Diese Starts tragen nicht zu $\mathcal{L}_{\mathrm{arith}}$ bei |

Die Lücke $\mathcal{L}_{\mathrm{arith}}\subsetneq\mathcal{L}$ bleibt damit **operational**:
grammatisch zulässige Wörter können ohne passendes $n$ existieren.

### Implementierung

| Artefakt | Inhalt |
|----------|--------|
| **Python** | `collatz_l_arith_test.py` — Grammatik $L(k)$ (BB-Verbot, endliche $C$-Ketten, EA nach $C^*_{\max}$), Vollliste für $k\leq 8$, Stichprobe darüber |
| **JSON** | `collatz_l_arith_test.json` — Ratios, minimale Gegenbeispiele, ehrliche Grenzen |
| **Lean** | `CollatzEabc/ArithLanguage.lean` — `isGrammarValid` (BB), `RealizableWord` |
| **pytest** | `tests/test_l_arith.py` |

**Grammatikregeln** (vgl. `collatz_equivalenz_e_infty.tex`, `collatz_schlussartikel_arxiv.tex` §C-Ketten):
1. $BB\notin\mathcal{L}$
2. Endliche $C$-Ketten mit Kapazität $\mathrm{cap}\geq 1$
3. Zwang $C^*_{\max}\to EA$ nach maximaler $C$-Kette

### Ehrliche Testergebnisse ($n\leq 10^6$, ungerade Starts mit vollständigem $\kappa$-Präfix)

| $k$ | $|L(k)|$ | $|L_{\mathrm{arith}}\cap L(k)|$ | Ratio | Methode |
|-----|----------|----------------------------------|-------|---------|
| 4 | 89 | 48 | 0,539 | Vollliste |
| 5 | 410 | 144 | 0,351 | Vollliste |
| 6 | 2091 | 432 | 0,207 | Vollliste |
| 7 | 11589 | 1295 | 0,112 | Vollliste |
| 8 | 68753 | 3836 | 0,056 | Vollliste |
| 10+ | $\gg 10^6$ | — | Stichprobe | $|L(10)|\approx 2{,}9\times 10^6$ — keine Vollliste |

**Minimales Gegenbeispiel:** $w=\mathrm{BE}$ (Länge 2) — grammatisch zulässig, aber kein
ungerades $n\leq 10^6$ mit $\kappa$-Präfix $\mathrm{BE}$.

**Interpretation (kein Beweisanspruch):** Die Ratio $|L_{\mathrm{arith}}|/|L|$ fällt mit $k$
— viele grammatische Wörter sind (auf endlicher Suchtiefe) nicht realisierbar. Das stützt
$\mathcal{L}_{\mathrm{arith}}\subsetneq\mathcal{L}$, beweist aber weder Collatz noch
$L_{\mathrm{arith}}=\mathcal{L}$.

**Offen:** Vollständige Charakterisierung von $\mathcal{L}_{\mathrm{arith}}$; Zusammenhang mit
treuer $\kappa$; unendliche Wörter vs. endliche Präfixe.
