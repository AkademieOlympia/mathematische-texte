# EABC-Resonanzhypothese der Zetafunktion

**Stand:** Juni 2026 · **Epistemische Warnung:** Forschungsnotiz im Tao-Stil
(Definition / Zeuge / Experiment / Theorem / Conjecture / Heuristik / Negativtest).
**Kein Collatz-Beweis. Kein Beweis der Riemann-Hypothese.**

**Terminologie:** *ERPC* bezeichnet in diesem Repo **dieselbe EABC-Struktur**
(Restklassen $E,A,B,C$ mod $12$, chirale Sensoren $\sigma$, $\chi$, $\iota_{\mathrm{chir}}$).
Es gibt kein separates „ERPC-Projekt neben EABC“ — die Resonanz-/Bernoulli-/Zeta-Theorie
ist die **EABC-Deutung** trivialer Zetawerte und von-Staudt-Clausen-Signaturen.

**Querverweise:** `Staudt.tex`, `collatz_kepler_gedankenexperiment.tex` (Bernoulli-Uhr),
`eabc_from_lean.py`, `collatz_bernoulli_schalen.pdf`, `energiedoku_eabc_c4_kohaerenz.tex`,
`collatz_eabc_bernoulli_sensor.py`, `CollatzEabc.BernoulliClock.lean`.

---

## 1. Zentrale Fragestellung

Die **EABC-Resonanzhypothese der Zetafunktion** fragt, ob die Folge der
EABC-Zustandsvektoren $V_1,V_2,\ldots$ an der äquidistanten Achse trivialer Nullstellen
eine **diskrete Resonanzstruktur der Primzahlen** trägt — und ob nichttriviale Nullstellen
als **spektrale Projektion** dieser Resonanz lesbar sind.

| Ebene | Label in diesem Text |
|-------|----------------------|
| EABC mod $12$, $V_n=(E_n,A_n,B_n,C_n)$ | **Definition** |
| von Staudt--Clausen, $\mathrm{PrimeSig}(B_{2n})$ | **Definition** (klassisches Theorem) |
| Sensor $\Phi$, chirale Observablen | **Definition / Zeuge** |
| Test 1 (EABC-Asymmetrien in $V_n$) | **Experiment** |
| Resonanzstruktur auf der trivialen Nullstellenachse | **Conjecture** |
| Abbildung $V_n \to \rho_k$ bzw. $\Delta t_k$ | **Conjecture** (starke Form) |
| Bernoulli-Lyapunov für Collatz | **Negativtest** (No-Go) |

---

## 2. Äquidistante Basis: triviale Nullstellen

**Definition.** Die trivialen Nullstellen der Riemannschen Zetafunktion liegen bei
\[
s = -2,\,-4,\,-6,\,\ldots,\qquad s_n = -2n\quad (n\ge 1),
\]
einer **flachen, äquidistanten Basisschicht** auf der negativen reellen Achse.

Jede Stufe $n$ kodiert einen rationalen Spezialwert:
\[
\zeta(1-2n) = -\frac{B_{2n}}{2n},
\]
wobei $B_{2n}$ die geraden Bernoulli-Zahlen sind. Bernoulli-Zahlen sind damit die
**Übersetzer** von Zeta-Spezialwerten in rationale Arithmetik — nicht identisch mit
einzelnen Primzahlen, sondern Filtrationsobjekte.

> **Boxed (Präzision):** Bernoulli $\neq$ Prim. Die Korrespondenz verläuft zwischen
> **Filtrationen** $\{B_{2n}\}_{n\ge 1} \longleftrightarrow \{p\ \text{prim}: p-1\mid 2n\}_{n\ge 1}$.

---

## 3. von Staudt--Clausen und Primsignaturen

**Definition (von Staudt--Clausen).** Für gerade Bernoulli-Zahlen $B_{2n}$ ($n\ge 1$):
\[
B_{2n} + \sum_{\substack{p\ \text{prim}\\ p-1\mid 2n}} \frac{1}{p} \;\in\; \mathbb{Z},
\qquad
\mathrm{den}(B_{2n}) = \prod_{\substack{p\ \text{prim}\\ p-1\mid 2n}} p.
\]
Die **Primzahlsignatur** (PrimeSig) ist
\[
P_n := \mathrm{PrimeSig}(B_{2n}) := \{p\ \text{prim} : p-1\mid 2n\}.
\]
Referenz: `Staudt.tex`, Mathlib `NumberTheory.Bernoulli`.

**Negativer Befund (bereits dokumentiert):** Bernoulli-Normschalen als Lyapunov-Funktion
für Collatz sind **No-Go** (`collatz_generalangriff_2026.md`). Dieser Zweig ist **Sensorik**,
nicht Hauptangriff.

---

## 4. Native EABC mod $12$ und Zustandsvektor $V_n$

Die vier EABC-Familien sind **Definition** (nicht Brücke zu einem externen Modul):
\[
E\equiv 1,\quad A\equiv 5,\quad B\equiv 7,\quad C\equiv 11 \pmod{12}.
\]
Implementierung: `eabc_from_lean.py` (`class_of`, `residue`).

Für $n\ge 1$ zerlegen wir $P_n$ in mod-$12$-Klassen und zählen:
\[
E_n := \#\{p\in P_n : p\equiv 1\},\quad
A_n := \#\{p\in P_n : p\equiv 5\},
\]
\[
B_n := \#\{p\in P_n : p\equiv 7\},\quad
C_n := \#\{p\in P_n : p\equiv 11\}.
\]
Der **EABC-Zustandsvektor** pro trivialer Nullstelle $s_n=-2n$ ist
\[
V_n := (E_n, A_n, B_n, C_n).
\]
Primzahlen $2,3$ liegen typischerweise in $P_n$, fallen aber außerhalb der vier Klassen; dann gilt
\[
E_n+A_n+B_n+C_n + \#\{p\in P_n : p\in\{2,3\}\} = |P_n|.
\]

**Sensor:** $\Phi(n) := V_n$. Implementierung: `collatz_eabc_bernoulli_sensor.py`
$\to$ `collatz_eabc_bernoulli_sensor.json` (JSON-Schlüssel `E_n`…`C_n` bzw. `V.E`…`V.C`;
Python-Felder `e,a,b,c` mappen auf $E_n,\ldots,C_n$).

---

## 5. EABC-Resonanzhypothese (Conjecture)

> **Conjecture (Resonanzhypothese).** Die Folge $V_1,V_2,V_3,\ldots$ ist **nicht zufällig**:
> sie trägt eine **diskrete Resonanzstruktur der Primzahlen** entlang der äquidistanten
> Achse trivialer Nullstellen. Nichttriviale Nullstellen
> \[
> \rho_k = \tfrac12 + \mathrm{i}\,t_k
> \]
> sind die **spektrale Projektion** dieser Resonanz auf die kritische Linie.

Intuition: Die flache Gitterachse $s=-2n$ liefert über Bernoulli und $P_n$ eine
arithmetische „Grundschwingung“; Abweichungen und Korrelationen in $V_n$ kodieren
Information, die im Spektrum der nichttrivialen Nullstellen wieder auftaucht.

**Label: Conjecture** — empirisch falsifizierbar (Abschnitt 8), nicht als Theorem behauptet.

---

## 6. Starke Form: Abbildung auf das Spektrum (Conjecture)

> **Conjecture (starke Form).** Es existiert eine (deterministische oder statistische)
> Abbildung
> \[
> \Phi : (E_n,A_n,B_n,C_n) \;\longrightarrow\; \rho_k
> \qquad\text{oder}\qquad
> \Phi_{\mathrm{stat}} : V_n \;\longrightarrow\; \Delta t_k := t_{k+1}-t_k,
> \]
> die EABC-Zustände mit dem Spektrum der nichttrivialen Nullstellen koppelt.

Die starke Form behauptet **keinen** geschlossenen Beweis der Riemann-Hypothese, sondern
eine testbare Brücke zwischen arithmetischer Resonanz ($V_n$) und spektralen Abständen
($\Delta t_k$). Fehlt jede systematische Korrelation bei ausreichend großem $n$, ist die
Hypothese **widerlegt**.

---

## 7. Geometrische Interpretation (Heuristik)

Fünfgliedrige Lesart als Projektionskette:

1. **Äquidistante Basis:** Triviale Nullstellen $s=-2n$ bilden ein flaches, gleichabständiges Gitter auf $\mathrm{Re}(s)<0$.
2. **Krümmung durch Bernoulli:** $\zeta(1-2n)=-B_{2n}/(2n)$ koppelt jedes Gitterpunkt an rationale Bernoulli-Daten.
3. **Diskrete Prim-Einbettung:** $P_n$ legt Primzahlen als lokale „Massen“ auf die Stufe $n$.
4. **EABC-Faserung:** $V_n=(E_n,A_n,B_n,C_n)$ ordnet diese Massen der nativen mod-$12$-Symmetrie zu (chirale Bilanz $\sigma$, $\chi$, $\iota_{\mathrm{chir}}$).
5. **Spektrale Projektion:** Nichttriviale Nullstellen $\rho_k=\tfrac12+\mathrm{i}t_k$ erscheinen als gekrümmte Projektion der Resonanz auf $\mathrm{Re}(s)=\tfrac12$; RH entspricht der Fixierung dieser Projektion.

**Label: Heuristik / Conjecture** — geometrische Metapher, kein Beweisanspruch.

---

## 8. Falsifizierbarkeit und Experimente

**Falsifikationskriterium.** Berechne $V_n$ für viele $n$, vergleiche mit bekannten
Nullstellenabständen $\Delta t_k$. **Keine** systematische Korrelation $\Rightarrow$
Resonanzhypothese **falsch** (für die getestete Kopplung). **Systematische** Struktur
$\Rightarrow$ es existiert eine arithmetisch–spektrale Brücke — weiter zu präzisieren.

| # | Test | Status | Label |
|---|------|--------|-------|
| **1** | EABC-Asymmetrien in $V_n$ ($\sigma$, $\chi$, $\iota_{\mathrm{chir}}$) | **implementiert** | Experiment |
| **2** | $V_n$ vs. $\Delta t_k$ (Korrelation, Kreuzspektrum) | **zukünftig** | Experiment |
| **3** | Krümmung $K_B(n)$ vs. $\pi(x)-\mathrm{Li}(x)$ | **zukünftig** | Experiment |

Test 1 prüft reproduzierbar, ob EABC-Asymmetrien in der PrimeSig-Filtration
systematisch auftreten oder im Zufallsbereich liegen. Test 2 ist der **entscheidende**
Falsifikationstest der Resonanzhypothese.

**Chirale Observablen** (Definition):
\[
\sigma(n) := (E_n+A_n) - (B_n+C_n) \quad\text{(EA/BC-Bilanz)},
\]
\[
\chi(n) := (E_n+B_n) - (A_n+C_n) \quad\text{(Diagonalen-Bilanz)},
\]
\[
\iota_{\mathrm{chir}}(n) := \mathrm{sgn}\bigl(\sigma(n)\cdot\chi(n)\bigr)\in\{-1,0,+1\}.
\]

---

## 9. Boxed Kurzform

$$\boxed{
\text{Triviale Nullstellen}
\;\to\;
\text{Bernoulli}
\;\to\;
\text{Primsignaturen}
\;\to\;
\text{EABC-Zustände } V_n
\;\to\;
\text{Nichttriviale Nullstellen}
}$$

---

## 10. Parallele Tracks (nicht verschmolzen)

Der EABC-Bernoulli-/Resonanz-Zweig läuft **parallel** zu anderen Sensoren:

| Track | Sensor | Artefakt |
|-------|--------|----------|
| **EABC-Resonanz** | $\Phi$, $V_n$ | `collatz_eabc_bernoulli_sensor.json` |
| **Morley** | $F_M$, $S_M$, $K_M$ | `collatz_morley_tm_numerik.py` |
| **$\kappa$ / Grammatik** | $\mathcal{L}_{\mathrm{arith}}$, $F_n$ | `collatz_forbidden_words.py` |

Morley und $\kappa$ bleiben eigenständige Angriffslinien. Der Resonanz-Zweig nutzt
**dieselbe** EABC-Sprache, aber andere Eingangsgrößen (von Staudt statt Collatz-Bahnen).

---

## 11. Tao-Labels und nächste Schritte

| Aussage | Tao-Label |
|---------|-----------|
| EABC mod $12$, $V_n=(E_n,A_n,B_n,C_n)$ | Definition |
| PrimeSig aus von Staudt--Clausen | Definition (Theorem) |
| $\Phi$, $\sigma$, $\chi$, $\iota_{\mathrm{chir}}$ | Definition + Experiment |
| Resonanzstruktur $V_1,V_2,\ldots$ | Conjecture |
| $\Phi:V_n\to\rho_k$ bzw. $\Delta t_k$ | Conjecture |
| RH als Fixraum-Projektion | Conjecture / Heuristik |
| Bernoulli-Lyapunov für Collatz | Negativtest |

**Nächste Schritte:**
1. Test 1 auf größeres $n$ und mit expliziten Zufalls-Nullmodellen vergleichen.
2. Test 2: $\Delta t_k$-Daten einbinden, Korrelationen und Falsifikation dokumentieren.
3. Lean-Schnittstelle: PrimeSig als `Finset` über `Nat.Prime` (analog `BernoulliClock.lean`).

---

*Epistemische Einordnung gemischt je Abschnitt. Definitionen und Experimente sind streng;
Conjectures sind explizit falsifizierbar. Keine Ebene ersetzt die nächsthöhere.*
