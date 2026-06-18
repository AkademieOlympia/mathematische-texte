# EABC-Holonomie: analytischer Beweisversuch (konservativ)

**Status:** Lemma-Skizze + Vermutung + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Tao-Labels:** Definition | Lemma-Skizze | Vermutung | Experiment

**Querverweise:**
- `collatz_eabc_fehlerterm_hypothese.md` — **kanonische Endform:** $N_\pm$, Hauptvermutung, Fehlerterm-Hypothese, $\widetilde{D}_E$
- `collatz_eabc_zyklus_holonomie.md` — kanonische Hierarchie, $\chi_{\mathrm{Hol}}$, $\mathrm{Hol}_E$
- `collatz_eabc_transport.md` — Übergangsgraph $G_E$, Transport $T_n$
- `collatz_eabc_transition_graph.py` — $\chi_{\mathrm{Pfad}}$, $\chi_{\mathrm{Hol}}$, Nullmodelle
- `collatz_eabc_holonomie_fehlerterm.py` / `.json` — $D_E$, $\widetilde{D}_E$, Lückenmuster, Chebyshev-Vergleich
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54

---

## 0. Endform (Verweis)

Die **kanonische Kurzform** (Setup, $N_\pm$, Hauptvermutung vs. Fehlerterm-Hypothese, $\widetilde{D}_E$, zentrale Frage) steht in `collatz_eabc_fehlerterm_hypothese.md`. Dieses Dokument liefert die **analytische Ausarbeitung** (Lemma-Skizzen, HL-Argument, Numerik).

---

## Epistemische Einordnung

> **Dieses Dokument ist ein analytischer Beweisversuch**, nicht ein Theorem.

> Die **starke Hypothese** $\mathrm{Hol}_E\neq 0$ (`collatz_eabc_zyklus_holonomie.md` §7) wird hier durch eine **konservative Symmetrie-Lesart** in Frage gestellt: $\mathrm{Hol}_E=0$ als Hauptterm, interessante Struktur nur im **Fehlerterm** $D_E$.

---

## 1. EABC-Klassen modulo $12$

**Definition (Restklassen).** Für $p>3$ Prim:
$$\kappa(p)\in\{E,A,B,C\},\qquad
E\equiv 1,\; A\equiv 5,\; B\equiv 7,\; C\equiv 11\pmod{12}.$$

**Label:** $\kappa$, Restklassen = **Definition** (`eabc_from_lean.py::class_of`, `EABC.lean`).

Die vier zulässigen Prim-Restklassen $>3$ bilden genau die Menge $\{1,5,7,11\}\subset\mathbb{Z}/12\mathbb{Z}$.

---

## 2. ABCEA und CEABC als zyklische Verschiebungen

**Lemma-Skizze (4-Zyklus).** Die Rotation $t\colon E\mapsto A\mapsto B\mapsto C\mapsto E$ induziert auf $V_4=\{E,A,B,C\}$ einen **4-Zyklus**
$$5 \xrightarrow{+2} 7 \xrightarrow{+4} 11 \xrightarrow{+2} 1 \xrightarrow{+4} 5 \pmod{12},$$
d. h. in Buchstaben $A\to B\to C\to E\to A$.

Die geschlossenen 5-Wörter der Zyklus-Holonomie sind:
$$\mathrm{ABCEA} = (A,B,C,E,A),\qquad \mathrm{CEABC} = (C,E,A,B,C).$$

Beide sind **zyklische Verschiebungen desselben geschlossenen 4-Zyklus** $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$; sie unterscheiden sich nur durch den **Startpunkt** auf dem Zyklus.

**Label:** 4-Zyklus, Wortpaar = **Lemma-Skizze**.

---

## 3. Gemeinsames Lückenmuster $(2,4,2,4)$ modulo $12$

**Lemma-Skizze (Lücken).** Schreibe die Restklassen der Buchstaben entlang des Zyklus:
$$\mathrm{ABCEA}\colon\; 5,\;7,\;11,\;1,\;5 \quad\Rightarrow\quad \text{Lücken }(2,4,2,4),$$
$$\mathrm{CEABC}\colon\; 11,\;1,\;5,\;7,\;11 \quad\Rightarrow\quad \text{Lücken }(2,4,2,4).$$

Formal: für aufeinanderfolgende Restklassen $r_0,r_1,r_2,r_3,r_4$ mit $r_4\equiv r_0$ gilt
$$\Delta_i := (r_{i+1}-r_i)\bmod 12 \in \{2,4\},\qquad (\Delta_0,\Delta_1,\Delta_2,\Delta_3)=(2,4,2,4).$$

**Label:** Lückenmuster = **Lemma-Skizze** (reine mod-$12$-Kombinatorik).

**Experiment:** `collatz_eabc_holonomie_fehlerterm.py::verify_gap_patterns`.

---

## 4. Einziger Unterschied: Startklasse $5$ vs. $11$

**Lemma-Skizze.** $\mathrm{ABCEA}$ startet in Klasse $A\equiv 5$; $\mathrm{CEABC}$ startet in Klasse $C\equiv 11$.

Unter der $t$-Rotation sind $A$ und $C$ **Antipoden** auf dem 4-Zyklus (Abstand $2$ Schritte). Die beiden Holonomie-Orientierungen sind daher **dieselbe geometrische Schleife**, nur mit umgekehrter zyklischer Phase — analog zu zwei gegenläufigen Parametrisierungen derselben geschlossenen Kurve.

**Label:** Startklassen-Dualität = **Lemma-Skizze**.

---

## 5. Hardy–Littlewood-Äquidistribution und $\mathrm{Hol}_E=0$

**Definition (Zählgrößen).** Für $X>0$ (Endform: `collatz_eabc_fehlerterm_hypothese.md`):
$$N_+(X) := \#\{n:\,p_{n+4}\le X,\; C_n=\mathrm{ABCEA}\},$$
$$N_-(X) := \#\{n:\,p_{n+4}\le X,\; C_n=\mathrm{CEABC}\}.$$
Legacy: $N_{\mathrm{ABCEA}}:=N_+$, $N_{\mathrm{CEABC}}:=N_-$.

**Definition (normierte Holonomie-Observable).**
$$\chi_{\mathrm{Hol}}(X)=\frac{N_+(X)-N_-(X)}{N_+(X)+N_-(X)}
=\frac{D_E(X)}{N_+(X)+N_-(X)},$$
wobei $D_E(X):=N_+(X)-N_-(X)$.

**Definition (Grenzwert).**
$$\mathrm{Hol}_E := \lim_{X\to\infty}\chi_{\mathrm{Hol}}(X),$$
sofern der Grenzwert existiert (`collatz_eabc_zyklus_holonomie.md` §6).

**Hauptvermutung (konservativ).** Unter der **Hardy–Littlewood-Heuristik** für die Verteilung aufeinanderfolgender Prim-Restklassen modulo $12$ (bzw. äquivalent: Gleichverteilung der $4$-Schritt-Fenster mit festem Lückenmuster $(2,4,2,4)$ auf dem $4$-Zyklus) gilt asymptotisch
$$N_+(X)\sim N_-(X),$$
und damit
$$\mathrm{Hol}_E = 0.$$

**Begründung (Skizze).** Beide Wörter realisieren **dieselbe** lokale Gap-Struktur auf demselben 4-Zyklus; der einzige Unterschied ist der Startpunkt $5$ vs. $11$. Wenn aufeinanderfolgende Prim-Restklassen (in einem geeigneten Sinne) equidistributed sind, darf kein systematischer Vorzeichen-Überschuss einer Phase gegenüber der anderen verbleiben.

**Label:** HL-Äquidistribution $\Rightarrow$ $\mathrm{Hol}_E=0$ = **Vermutung** (konservativ); Zähldefinitionen = **Definition**.

---

## 6. Bias nur im Fehlerterm — Chebyshev-Analogie

**Lemma-Skizze (Fehlerterm).** Selbst wenn $\mathrm{Hol}_E=0$, bleibt
$$D_E(X)=N_+(X)-N_-(X)$$
nicht-trivial. Für endliche $X$ ist $\chi_{\mathrm{Hol}}(X)=D_E(X)/(N_++N_-)$ typischerweise $\neq 0$.

**Fehlerterm-Hypothese (stärker).** $D_E$ trägt nichttrivialen Chebyshev-artigen Bias, gesteuert durch Nullstellen der Dirichlet-$L$-Funktionen modulo $12$ (`collatz_eabc_fehlerterm_hypothese.md` §3).

**Heuristik (Chebyshev-Bias-Analogie).** Bei Chebyshev-Bias überwiegen oft Primzahlen $p\equiv 3\pmod 4$ gegenüber $p\equiv 1\pmod 4$ in endlichen Fenstern, obwohl die asymptotische Dichte gleich ist. Analog: $D_E(X)$ kann **oszillierende, vorzeichenbehaftete endliche Abweichungen** tragen, ohne dass $\mathrm{Hol}_E\neq 0$ folgt.

**Label:** Fehlerterm-Bias = **Lemma-Skizze** + **Heuristik**; kein Theorem.

**Experiment:** `collatz_eabc_holonomie_fehlerterm.py::chebyshev_bias_comparison`.

---

## 7. Satzskizze (5 Schritte)

**Satzskizze (konservative Holonomie-Vanishing).** *Unter HL-Äquidistribution der aufeinanderfolgenden Prim-Restklassen modulo $12$ gilt $\mathrm{Hol}_E=0$.*

| Schritt | Inhalt |
|--------:|--------|
| 1 | $\kappa(p)\in\{1,5,7,11\}\pmod{12}$ — vier zulässige Prim-Restklassen. |
| 2 | $\mathrm{ABCEA}$ und $\mathrm{CEABC}$ sind zyklische Realisierungen desselben 4-Zyklus $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$. |
| 3 | Beide tragen dasselbe Lückenmuster $(2,4,2,4)$ in $\mathbb{Z}/12\mathbb{Z}$. |
| 4 | Der einzige Unterschied ist die Startklasse $A\equiv 5$ vs. $C\equiv 11$ (Antipoden auf dem 4-Zyklus). |
| 5 | HL-Äquidistribution / symmetrische Zählung der Phasen $\Rightarrow$ $N_+(X)\sim N_-(X)$ $\Rightarrow$ $\chi_{\mathrm{Hol}}(X)\to 0$. |

**Offene Lücke:** Schritt 5 benötigt einen **beweisbaren** Äquidistributionssatz für $5$-Fenster auf der Primfolge — derzeit **Vermutung**, nicht Theorem.

**Label:** Satzskizze = **Lemma-Skizze**; HL-Schritt = **Vermutung**.

---

## 8. Interessante Größen: $D_E$ und $\widetilde{D}_E$

**Definition (Differenz und normalisierter Fehlerterm).**
$$D_E(X) := N_+(X) - N_-(X),$$
$$\widetilde{D}_E(X) := \frac{D_E(X)}{\sqrt{N_+(X)+N_-(X)}}.$$

**Zentrale Frage.** Verhält sich $\widetilde{D}_E$ wie reines Rauschen, oder zeigt stabile Vorzeichenasymmetrie / Oszillationsstruktur?

**Interpretation.**
- $D_E(X)$: absoluter Vorzeichen-Überschuss ABCEA gegenüber CEABC (Chebyshev-artig).
- $\widetilde{D}_E(X)$: $\sqrt{N}$-normalisierte Abweichung; bei $\mathrm{Hol}_E=0$ und „zufälligen" Phasenfluktuationen erwartet man $|\widetilde{D}_E(X)|=O(1)$, nicht lineares Wachstum in $X$.
- EABC-Holonomie als **sekundärer Bias** im Prime Race gegenläufiger Zyklusorientierungen — **nicht** permanenter Hauptterm.

**Label:** $D_E$, $\widetilde{D}_E$ = **Definition**; Grenzverhalten = **Vermutung** / **Experiment**.

**Experiment:** `collatz_eabc_holonomie_fehlerterm.py` — Tabelle für $X\in\{10^3,10^4,10^5,10^6\}$.

---

## 9. Boxed Schlussfolgerungen

$$\boxed{\;\mathrm{Hol}_E = 0\quad\text{als Hauptterm unter mod-12-Symmetrie und HL-Äquidistribution.}\;}$$

$$\boxed{\;\text{Interessante Forschungsfrage:}\;D_E(X)\;\text{— Bias und Oszillation im Fehlerterm (Chebyshev-Analogie), nicht }\mathrm{Hol}_E\neq 0.\;}$$

**Abgrenzung zur starken Hypothese** (`collatz_eabc_zyklus_holonomie.md` §7):

| Lesart | Aussage | Label |
|--------|---------|-------|
| **Stark** | $\mathrm{Hol}_E\neq 0$ | **Hypothese** (empirisch: $\chi_{\mathrm{Hol}}(10^6)\approx 0{,}12$) |
| **Konservativ (hier)** | $\mathrm{Hol}_E=0$, Struktur in $D_E$ | **Vermutung** + **Experiment** |

Die empirische Beobachtung $N_+>N_-$ bei endlichen $X$ widerspricht **nicht** der konservativen Lesart: sie betrifft den **Fehlerterm** $D_E$, nicht den Grenzwert $\mathrm{Hol}_E$. ABCE/CEAB-Asymmetrie = Fehlerterm-/Bias-Phänomen.

---

## Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $E\equiv 1$, $A\equiv 5$, $B\equiv 7$, $C\equiv 11\pmod{12}$ | **Definition** |
| ABCEA/CEABC als zyklische Verschiebungen desselben 4-Zyklus | **Lemma-Skizze** |
| Gemeinsames Lückenmuster $(2,4,2,4)$ | **Lemma-Skizze** |
| Startklassen-Dualität $5$ vs. $11$ | **Lemma-Skizze** |
| $N_+\sim N_-$ $\Rightarrow$ $\mathrm{Hol}_E=0$ (Hauptvermutung) | **Vermutung** |
| $D_E$-Bias, $L$-Funktionen mod $12$ (Fehlerterm-Hypothese) | **Hypothese** |
| $D_E$-Bias als Chebyshev-Analogie | **Heuristik** |
| Satzskizze (5 Schritte) | **Lemma-Skizze** |
| $D_E$, $\widetilde{D}_E$ | **Definition** |
| Numerik $X\le 10^6$ | **Experiment** |
| $\mathrm{Hol}_E\neq 0$ (stark) | **Hypothese** (`collatz_eabc_zyklus_holonomie.md`) |

---

*Kanonsiche Notiz: Dieser Beweisversuch **ersetzt nicht** die Zyklus-Holonomie-Definitionen in `collatz_eabc_zyklus_holonomie.md`. Er liefert eine **alternative analytische Lesart**: Vanishing des Hauptterms, strukturierte Fluktuation im Fehlerterm. Die numerische Entscheidung zwischen starker und konservativer Hypothese erfordert Grenzwert- und Fehlerterm-Analyse bei $X\to\infty$.*
