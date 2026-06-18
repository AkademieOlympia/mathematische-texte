# EABC-Evolutionspfad und offene analytische Fragen

**Status:** Forschungsnotiz — Evolutionstabelle, Wachstumsszenarien, Dirichlet-Zerlegung  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Vermutung | Hypothese | Experiment | Analogie

**Querverweise:**
- `collatz_eabc_zirkulationshypothese.md` — **kanonische Hypothese:** $N_\pm$, $C_E$, $D_E$, $S_E$
- `collatz_eabc_zirkulation_spektral.md` — diskrete 1-Form $\alpha$, $\mathrm{Spec}(L_E)$
- `collatz_eabc_fehlerterm_hypothese.md` — Teilhypothese Fehlerterm $D_E$, $\widetilde{D}_E$
- `collatz_eabc_holonomie_beweisversuch.md` — mod-$12$-Symmetrie, $\mathrm{Hol}_E=0$
- `collatz_eabc_bell_holonomie.md` — **sekundäre** Analogie Bell/CHSH
- `collatz_eabc_sagnac.md` — **Intuition only:** Sagnac-Bild
- `collatz_eabc_holonomie_fehlerterm.py` — $N_\pm$, $D_E$, $\widetilde{D}_E$
- `collatz_eabc_D_growth.py` — Wachstumsdiagnostik $D_E(X)$, Charakter-Stub
- `collatz_eabc_graph_laplacian.py` — $\mathrm{Spec}(L_E)$
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54 / PR #59

---

## 1. Evolution (vier Phasen)

Der EABC-Forschungsstrang entwickelt sich **formal** von Paarstatistik zu Spektralgeometrie:

$$\boxed{\;\text{Bell} \;\to\; \text{Sagnac} \;\to\; C_E \;\to\; \mathrm{Spec}(L_E).\;}$$

| Phase | Metapher (didaktisch) | Mathematischer Kern | Objekt |
|:-----:|----------------------|---------------------|--------|
| 1 | Bell / CHSH | Korrelationen auf Kanten | $E(a,b)$, Paarstatistik auf $G_E$ |
| 2 | Sagnac | Orientierung auf $C_4$ | $\gamma^\pm$, ABCEA vs. CEABC |
| 3 | Zirkulation | diskrete 1-Formen | $\alpha$, $\omega(e)$, $C_E(X)=\oint_\gamma\alpha$ |
| 4 | Spektral | Laplace-Spektrum, Fehlerterme | $L_E$, $\mathrm{Spec}(L_E)$, $D_E(X)$ |

**Epistemische Lesart:** Phasen 1–2 sind **Analogien** (Bell, Sagnac); Phasen 3–4 sind der **mathematische Kern** (Zirkulation, Spektralgeometrie). Details: `collatz_eabc_zirkulation_spektral.md` §1, §9.

**Shift der Fragestellung:**
$$E(a,b)\;\text{(Kanten-Korrelation)}\;\longrightarrow\;\oint_\gamma \alpha\;\text{(Zyklusstatistik)}\;\longrightarrow\;D_E(X)\;\longrightarrow\;\mathrm{Spec}(L_E).$$

**Label:** Evolutionstabelle = **Heuristik** (Begriffshierarchie); Bell/Sagnac = **Analogie**; $C_E$, $\mathrm{Spec}(L_E)$ = **Definition**.

---

## 2. $H_1(C_4)$-Interpretation

**Gerichteter 4-Zyklus** auf $V=\{E,A,B,C\}$:
$$E \to A \to B \to C \to E.$$

**Fundamentalgruppe / erste Homologie** (abgekürzt $C_4$-Gerüst):
$$H_1(C_4,\mathbb{Z}) \cong \mathbb{Z}.$$

**Zwei Orientierungen** desselben Zyklus:

| Orientierung | Wort | Klasse in $H_1(C_4)$ |
|--------------|------|----------------------|
| $\gamma^+$ | ABCEA | $+1$ |
| $\gamma^-$ | CEABC | $-1$ |

Beide Wörter tragen dasselbe Lückenmuster $(2,4,2,4)$ mod $12$; sie unterscheiden sich durch **zyklische Verschiebung** (Start $A\equiv 5$ vs. $C\equiv 11$).

**Observable:**
$$D_E(X) := N_+(X) - N_-(X) = \sum_{\gamma \le X} \mathrm{sgn}(\gamma),$$
wobei $\mathrm{sgn}(\gamma^+)=+1$, $\mathrm{sgn}(\gamma^-)=-1$.

$$\boxed{\;D_E(X)\;\text{misst die Überrepräsentation einer Orientierung in } H_1(C_4,\mathbb{Z})\text{ entlang der Primfolge.}\;}$$

Unter Symmetrie-Hypothese mod $12$: $\mathbb{E}[D_E(X)]=0$ im Hauptterm; der **interessante** Anteil ist der strukturierte Fehlerterm (Prime Race). Kanonische Definition: `collatz_eabc_zirkulationshypothese.md` §2–3.

**Label:** $H_1(C_4)$-Lesart = **Definition** / **Interpretation**; $\mathbb{E}[D_E]=0$ = **Vermutung** (Hauptterm).

---

## 3. Offene Frage: Wachstum von $D_E(X)$ für $X\to\infty$

**Zentrale Wachstumsfrage:** Wie schnell wächst $|D_E(X)|$ relativ zu $X$ und zur Fensterzahl $N_+(X)+N_-(X)$?

| Szenario | Asymptotik | Bedeutung |
|:--------:|------------|-----------|
| **A** | $D_E(X) = O(1)$ | Effekt verschwindet asymptotisch absolut |
| **B** | $D_E(X) = O(\log X)$ | langsame logarithmische Drift |
| **C** | $D_E(X) = O(\sqrt{X})$ | Chebyshev-Race-artig (klassische Größenordnung) |
| **D** | $D_E(X) = c\,X^\alpha$, $\alpha>0$ | starke strukturelle Asymmetrie |

**Zusatz-Observable:** $\widetilde{D}_E(X) := D_E(X)/\sqrt{N_+(X)+N_-(X)}$. Wächst $\widetilde{D}_E$ oder stabilisiert sie sich?

**Hauptterm vs. Fehlerterm:** $S_E(X)=D_E/(N_++N_-)\to 0$ (Vermutung) ist **verträglich** mit wachsendem $|D_E(X)|$, solange der Nenner schneller wächst.

**Experiment:** `collatz_eabc_D_growth.py` — Gitter $X\in[10^3,10^6]$, heuristische Modellwahl (Konstante, $\log X$, $\sqrt{X}$, Potenzgesetz).

**Label:** Wachstumsszenarien = **offene Frage**; numerische Klassifikation = **Experiment**.

---

## 4. Dirichlet-Zerlegung (nächster analytischer Schritt)

**Ziel:** $D_E(X)$ als **mod-$12$-Chebyshev-Race** über Dirichlet-Charaktere formulieren.

**Heuristische Zerlegung:**
$$D_E(X) = \sum_{\chi \bmod 12} a_\chi \sum_{p \le X} \chi(p),$$
wobei $\chi$ über die Charaktere auf $(\mathbb{Z}/12\mathbb{Z})^\times=\{1,5,7,11\}$ läuft (isomorph zu $C_2\times C_2$).

**Verbindungen:**
- **Dirichlet-$L$-Funktionen** $L(s,\chi)$ mod $12$ und ihre Nullstellen
- **Chebyshev-Races** ($\pi(x;4,3)-\pi(x;4,1)$ als Referenz)
- **Explizite Formel** — Oszillationen in $D_E$ aus Nullstellen-Beiträgen

$$\boxed{\;\text{Stärkster nächster Schritt: } D_E \text{ als mod-12-Chebyshev-Race via Dirichlet-Charaktere.}\;}$$

**Stub-Implementierung:** `collatz_eabc_D_growth.py::dirichlet_decomposition_stub` — Koeffizienten $a_\chi$ per lineare Projektion (experimentell, nicht Theorem).

**Label:** Dirichlet-Zerlegung = **Hypothese** / **Forschungsprogramm**; $a_\chi$-Numerik = **Experiment**.

---

## 5. Spektral: $C_4$-Laplace und Prim-gewichtete Eigenmoden

**Kombinatorischer $C_4$-Laplace** (unabhängig von Primdaten, reines Gerüst):

Auf dem gerichteten 4-Zyklus mit symmetrisierter Laplace $L_{C_4}^{\mathrm{sym}}$ gilt bekanntlich:
$$\mathrm{Spec}(L_{C_4}^{\mathrm{sym}}) = \{0,\,2,\,2,\,4\}$$
(Eigenwert $0$ konstant, doppelte Schwingungsmode bei $2$, maximale Mode bei $4$).

**Prim-induzierter Transportgraph** $G_E(X)$ aus $\kappa(p_n)\to\kappa(p_{n+1})$:
$$L_E(X) = D_{\mathrm{out}} - A_E(X),\qquad \mathrm{Spec}(L_E(X))$$
variabel mit Obergrenze $X$. Experiment: `collatz_eabc_graph_laplacian.py`.

**Offene spektrale Fragen:**
1. Tragen **prim-gewichtete Eigenmoden** von $L_E(X)$ Information über $D_E(X)$?
2. Korreliert die **Spektrallücke** $\lambda_1-\lambda_0$ mit $|D_E(X)|$ oder $S_E(X)$?
3. Lässt sich $D_E(X)$ als **Linearkombination** von Eigenmoden-Projektionen auf die Orientierungsklassen $\gamma^\pm$ schreiben?

$$\boxed{\;\text{Primfolge} \to G_E \to L_E \to \mathrm{Spec}(L_E)\;\leftrightarrow\;D_E(X)\;\text{(Fehlerterm-Kopplung offen).}\;}$$

**Label:** $\mathrm{Spec}(L_{C_4})=\{0,2,2,4\}$ = **Definition** (kombinatorisch); Kopplung $D_E\leftrightarrow$ Eigenmoden = **offene Frage**.

---

## 6. Schlussfolgerung und Forschungsagenda

| Priorität | Schritt | Werkzeug |
|:---------:|---------|----------|
| 1 | $D_E$ als mod-$12$-Chebyshev-Race formalisieren | Dirichlet-Charaktere, $L(s,\chi)$ |
| 2 | Wachstumsszenario A–D empirisch eingrenzen | `collatz_eabc_D_growth.py` |
| 3 | $a_\chi$ schätzen und mit $L$-Nullstellen vergleichen | Charakter-Stub, explizite Formel |
| 4 | Spektrallücke / Eigenmoden vs. $D_E$ | `collatz_eabc_graph_laplacian.py` |

$$\boxed{\;\text{Kern: } D_E \text{ ist Prime Race in } H_1(C_4) \text{ — analytisch via Dirichlet mod } 12 \text{, spektral via } \mathrm{Spec}(L_E).\;}$$

---

## 7. Chebyshev-Bias: gemeinsame Ursache?

### 7.1 Was ist die etablierte „Ursache“ des klassischen Chebyshev-Bias?

**Label: Theorem (bedingt auf GRH+LI)** — Rubinstein–Sarnak (1994) und Nachfolger.

Der Bias $\pi(x;4,3)>\pi(x;4,1)$ für „die meisten“ $x$ (logarithmische Dichte $\approx 0{,}996$ unter GRH+LI) ist **kein** verstecktes Primzahlgesetz und **kein** algebraischer Defekt der EABC-Klassifikation. Er entsteht aus einem **Prime Race** zwischen zwei Restklassen desselben Modulus:

$$\pi(x;4,3)-\pi(x;4,1) = -\frac{1}{2}\sum_{\substack{p\le x\\ p\equiv 1(4)}}\!\!\Bigl(1-\frac{p}{x}\Bigr) + \cdots$$

Die Vorzeichenstruktur wird durch die **Nullstellen** der zugehörigen Dirichlet-$L$-Funktionen $L(s,\chi_4)$ gesteuert; unter der **Linear Independence (LI)**-Hypothese über Nullstellen bestimmen die niedrigsten kritischen Nullstellen, welche Restklasse in endlichen Fenstern typischerweise vorne liegt. Mod $4$ gewinnt Klasse $3$, weil die erste Nullstelle von $L(s,\chi_4)$ negativ imaginär ist.

**Epistemische Kurzform:** Die „wahre Ursache“ ist **analytisch** (Nullstellen von $L$-Funktionen + LI), nicht **kombinatorisch** (keine versteckte Primregel).

### 7.2 Ist $D_E$ derselbe Race?

**Label: Definition** vs. **Heuristik** — $D_E$ ist ein **anderes** Observable.

| Aspekt | Klassischer Chebyshev-Bias | EABC-Zirkulation $D_E$ |
|--------|---------------------------|------------------------|
| Objekt | $\pi(x;q,a)-\pi(x;q,b)$ | $N_{\mathrm{ABCEA}}(X)-N_{\mathrm{CEABC}}(X)$ |
| Struktur | marginales Zählen mod $q$ | **5-Fenster** auf aufeinanderfolgenden Prim-EABC-Labels |
| Modulus | typisch $q=4$ | indirekt mod $12$ (Pfad $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$) |
| Symmetrie | zwei Restklassen | zwei **Orientierungen desselben** $C_4$-Zyklus, gleiches Lückenmuster $(2,4,2,4)$ |
| Unterschied | Klassen $1$ vs. $3$ | Startklasse $5$ vs. $11$ (Antipoden auf dem 4-Zyklus) |

$D_E$ misst keinen marginalen mod-$12$-Race $\pi(x;12,5)-\pi(x;12,11)$, sondern die Überrepräsentation geschlossener **5-Wörter** $\mathrm{ABCEA}$ gegenüber $\mathrm{CEABC}$ entlang der Primfolge. Beide Wörter erzwingen dieselbe lokale Gap-Struktur; der einzige kombinatorische Unterschied ist die zyklische Phase.

**Label:** „Prime Race in $H_1(C_4)$“ = **Definition** / **Interpretation**; Gleichsetzung mit $\pi(x;4,3)-\pi(x;4,1)$ = **Heuristik** (nicht bewiesen).

### 7.3 Kann eine Charakterzerlegung von $D_E$ die gemeinsame Mechanik offenlegen?

**Label: Hypothese** — Forschungsprogramm, kein Theorem.

Die naive Zerlegung aus §4,
$$D_E(X) \stackrel{?}{=} \sum_{\chi\bmod 12} a_\chi \sum_{p\le X}\chi(p),$$
ist **unzureichend**: $D_E$ hängt von der **gemeinsamen Verteilung** von $(\kappa(p_n),\ldots,\kappa(p_{n+4}))$ ab, nicht von marginalen Einzelsummen $\sum_p\chi(p)$. Eine korrekte analytische Formulierung müsste $5$-Korrelationsindikatoren oder Äquidistribution von $5$-Fenstern mod $12$ involvieren (offene Lücke in `collatz_eabc_holonomie_beweisversuch.md` Schritt 5).

**Gemeinsame Mechanik — wann plausibel?**

- **Plausibel (bedingt):** Wenn sich $D_E$ als Linearkombination von **Fehlertermen** klassischer mod-$12$-Races schreiben lässt, wären dieselben $L(s,\chi)$-Nullstellen die treibende Größe — analog Rubinstein–Sarnak, nur mit anderen Koeffizienten $a_\chi$ aus der Fenster-Geometrie.
- **Nicht plausibel:** Dass $D_E$ eine **neue** „wahre Ursache“ jenseits der $L$-Funktionen enthüllt; die moderne Antwort auf Chebyshev ist bereits vollständig in der Nullstellen-/LI-Sprache.

**Label:** geteilte $L$-Nullstellen-Mechanik = **Hypothese**; neue Ursache jenseits von $L$-Funktionen = **widerlegt** (für klassischen Bias).

### 7.4 Numerischer Schnellvergleich ($X\le 10^6$)

**Label: Experiment** — `collatz_eabc_holonomie_fehlerterm.py::chebyshev_bias_comparison`, ergänzende Auswertung.

| $X$ | $D_E$ | $\widetilde{D}_E$ | $\pi(X;4,3)-\pi(X;4,1)$ | $\pi(X;12,5)-\pi(X;12,11)$ |
|----:|------:|------------------:|-------------------------:|---------------------------:|
| $10^4$ | $+2$ | $0{,}54$ | $+9$ | $+2$ |
| $10^5$ | $+9$ | $1{,}08$ | $+24$ | $+12$ |
| $10^6$ | $+59$ | $2{,}66$ | $+146$ | $-43$ |

Beobachtungen:

1. $D_E>0$ (ABCEA-Führung) im gesamten getesteten Bereich — **kein** Vorzeichenwechsel (im Gegensatz zum oszillierenden klassischen Race).
2. $\widetilde{D}_E$ wächst weiter ($\approx 2{,}7$ bei $10^6$), nicht $O(1)$-stabilisiert — Wachstumsszenario C/D aus §3 noch offen.
3. Korrelation $D_E$ vs. mod-$4$-Differenz auf grobem Gitter: $r\approx 0{,}94$ (beide monoton positiv); vs. mod-$12$-Race $(5)-(11)$ bei $10^6$: **entgegengesetztes Vorzeichen** ($+59$ vs. $-43$).
4. Verhältnis $D_E/(\pi_4(3)-\pi_4(1))$ ist nicht stabil ($\approx 0{,}22$–$0{,}68$) — kein einfacher Proportionalitätsfaktor.

**Label:** Numerik = **Experiment**; daraus folgt **kein** Theorem über gemeinsame Ursache.

### 7.5 Direkte Antwort

$$\boxed{\;\text{Kann } D_E \text{ die „wahre Ursache“ des Chebyshev-Bias finden?}\;\Rightarrow\;\textbf{nein}\;\text{(sie ist bereits bekannt).}\;}$$

$$\boxed{\;\text{Kann } D_E \text{ dieselbe Mechanik (}L\text{-Nullstellen mod }12\text{) wie mod }4 \text{ tragen?}\;\Rightarrow\;\textbf{bedingt ja}\;\text{— wenn die Fensterzerlegung gelingt.}\;}$$

| Frage | Antwort | Was fehlt |
|-------|---------|-----------|
| Neue „wahre Ursache“ des klassischen Bias? | **Nein** | Rubinstein–Sarnak: $L$-Nullstellen + LI |
| Ist $D_E$ derselbe Race? | **Nein** | 5-Fenster-Orientierungsrace $\neq$ $\pi(x;q,a)-\pi(x,q,b)$ |
| Gemeinsame Mechanik über Charaktere mod $12$? | **Bedingt** | Beweis: $D_E=\sum a_\chi E_\chi(X)$ mit expliziten $a_\chi$; Äquidistribution von $5$-Fenstern |
| Numerische Evidenz für Identität mit mod-$4$-Bias? | **Nein** | Kein stabiles Verhältnis; mod-$12$-$(5,11)$-Race bei $10^6$ contra |

**Hinweis:** `collatz_eabc_D_growth.py` liefert Wachstumsdiagnostik und einen **experimentellen** Charakter-Stub ($a_\chi$ per lineare Projektion); die Zerlegung ist **kein Theorem** und ersetzt nicht die $5$-Fenster-Korrelationsstruktur.

---

## 8. Python-Symbolzuordnung

| LaTeX | Python | Modul |
|-------|--------|-------|
| $D_E(X)$ | `D_E` | `collatz_eabc_D_growth`, `collatz_eabc_holonomie_fehlerterm` |
| $\widetilde{D}_E(X)$ | `D_tilde_E` | `collatz_eabc_D_growth` |
| Wachstumsszenario | `growth_scenario` | `collatz_eabc_D_growth` |
| $a_\chi$ (Stub) | `dirichlet_coefficients` | `collatz_eabc_D_growth` |
| $\mathrm{Spec}(L_E)$ | `eigenvalues_symmetrized` | `collatz_eabc_graph_laplacian` |

---

*Evolutionspfad und offene Fragen ergänzen die kanonische Zirkulationshypothese. Numerische Wachstumsdiagnostik: `collatz_eabc_D_growth.py`.*
