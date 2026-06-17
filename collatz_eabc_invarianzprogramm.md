# EABC-Invarianzprogramm (mathematisch strikt)

**Kanonsiche Formulierung** des EABC-Invarianzprogramms — epistemisch scharf getrennt von
§19–§22 in `collatz_eabc_bernoulli_uebersetzung.md` (Holographie, Fünfstufenprogramm,
Fossilien-/Messreihen-Metapher). Dieses Dokument ist die **Definition / Forschungsfrage**-
Ebene im Tao-Sinn; die Übersetzungsdatei verweist hierher aus §22.

**Label-Schema (Tao):** Definition | Theorem | Conjecture | Experiment | Heuristik | Forschungsvision.

---

## 1. EABC-Klassifikation und Zählvektoren

### Definition 1 (κ, EABC-Klassen)

Sei $P_{>3}=\{p\in\mathbb{P}: p>3\}$. Die Abbildung
\[
\kappa: P_{>3} \longrightarrow \{E,A,B,C\}
\]
ordne jeder Primzahl $p>3$ ihre **native mod-$12$-Restklasse** zu:
\[
E\equiv 1,\quad A\equiv 5,\quad B\equiv 7,\quad C\equiv 11 \pmod{12}.
\]
**Label:** **Definition** (identisch mit `eabc_from_lean.class_of` / Lean `CollatzEabc.Kappa`).

### Definition 2 (Zählvektor $V(x)$)

Für $x\in\mathbb{R}$, $x\ge 5$, definiere den **EABC-Zählvektor**
\[
V(x)=\bigl(E(x),A(x),B(x),C(x)\bigr)\in\mathbb{N}_0^4,
\]
wobei $E(x)$ (bzw.\ $A(x)$, $B(x)$, $C(x)$) die Anzahl der Primzahlen $p\le x$ mit
$p>3$ und $\kappa(p)=E$ (bzw.\ $A$, $B$, $C$) ist.

Primzahlen $2$ und $3$ liegen **außerhalb** des Definitionsbereichs von $\kappa$ und
entfallen in $V(x)$.

**Label:** **Definition**.

### Definition 3 (Simplex-Anteil $S(x)$)

Sei $\pi_{>3}(x):=E(x)+A(x)+B(x)+C(x)=|P_{>3}\cap[2,x]|$. Der **Simplex-Anteil** ist
\[
S(x)=\frac{1}{\pi_{>3}(x)}\,V(x)\in\Delta_3,
\]
wobei $\Delta_3=\{(e,a,b,c)\in[0,1]^4 : e+a+b+c=1\}$ das **Standard-$3$-Simplex** in
$\mathbb{R}^4$ bezeichnet.

In der Literatur oft $\pi(x)$ für die klassische Primzahlzählfunktion; hier ist der
Nenner **explizit** $\pi_{>3}(x)$, damit $S(x)\in\Delta_3$ (Summe $1$).

**Label:** **Definition**.

---

## 2. EABC-Invarianten

### Definition 4 (EABC-Invariante)

Eine Abbildung $\mathcal{I}:\Delta_3\to\mathbb{R}$ heißt **EABC-Invariante**, wenn
**eine** der folgenden Bedingungen erfüllt ist:

1. $\displaystyle\lim_{x\to\infty}\mathcal{I}\bigl(S(x)\bigr)$ existiert (im Sinne endlicher
   Grenzwerte), oder
2. die Folge $\bigl(\mathcal{I}(S(x))\bigr)_{x\ge 5}$ hat **asymptotisch beschränkte
   Schwankungen** (d.\,h.\ $\sup_{x\ge x_0}|\mathcal{I}(S(x))|<\infty$ für ein $x_0$ und
   keine unbeschränkte Drift).

**Label:** **Definition** (Forschungsprogramm, nicht Behauptung eines konkreten Grenzwerts).

### Beispiel (χ-Observable)

\[
\chi(x)=\frac{(E(x)+C(x))-(A(x)+B(x))}{\pi_{>3}(x)}
=\mathcal{I}_\chi\bigl(S(x)\bigr),\quad
\mathcal{I}_\chi(e,a,b,c)=(e+c)-(a+b).
\]
Die Observable $\chi(x)$ ist **EABC-Invariante** im Sinne von Definition 4, **falls**
der Grenzwert $\lim_{x\to\infty}\chi(x)$ existiert **oder** $\chi(x)$ asymptotisch stabil
oszilliert (beschränkte Fluktuationen).

**Label:** **Definition** (Observable) + **Conjecture** (Grenzwert/Stabilität — offen).

**Experiment:** `collatz_eabc_invarianzprogramm.py` $\to$ `collatz_eabc_invarianzprogramm.json`.

---

## 3. Prim-Vierlinge und Signaturen

### Definition 5 (Quadrupel und Signatur σ)

Ein **Prim-Vierling** (Quadruplet) ist
\[
Q=(p,p+2,p+6,p+8),
\]
wobei alle vier Einträge prim sind. Die **EABC-Signatur** ist
\[
\sigma(Q)=\bigl(\kappa(p),\kappa(p+2),\kappa(p+6),\kappa(p+8)\bigr)\in\Sigma_4,
\]
wobei $\Sigma_4=\{E,A,B,C\}^4$ die Menge aller $4$-Tupel-Signaturen ist ($|\Sigma_4|=256$).

**Label:** **Definition** (Struktur); Vorkommen einzelner Signaturen = **Experiment**.

---

## 4. Forschungsfragen

### Forschungsfrage 1 (limitierende Verteilung μ)

Existiert eine Wahrscheinlichkeitsverteilung $\mu$ auf $\Sigma_4$, sodass die empirischen
Häufigkeiten der Signaturen $\sigma(Q)$ über Prim-Vierlinge $Q$ mit $p\le x$ gegen $\mu$
konvergieren, wenn $x\to\infty$?

**Label:** **Forschungsfrage** (offen).

### Forschungsfrage 2 (lineare Invariante mit Konstante)

Existiert eine nicht-triviale Abbildung $\mathcal{I}(E,A,B,C)$ (polynomial oder linear in
den Komponenten von $S(x)$), sodass
\[
\mathcal{I}\bigl(S(x)\bigr)=c+o(1)
\]
für eine Konstante $c\in\mathbb{R}$ und $x\to\infty$?

**Label:** **Forschungsfrage** (offen).

---

## 5. Arbeitshypothese (Conjecture)

> **Arbeitshypothese (EABC-Invarianzprogramm).**
> Es existieren **nicht-triviale EABC-Invarianten** $\mathcal{I}(E,A,B,C)$, die **nicht**
> allein aus der klassischen Primzahlzählfunktion $\pi(x)$ rekonstruierbar sind.

**Explizite Nicht-Behauptungen:**

- Es wird **nicht** behauptet, dass eine $8$D-Projektion oder EABC die Primzahlen **erzeugt**.
- Es wird **nicht** behauptet, dass $\chi(x)$ oder eine andere Observable bereits bewiesenermaßen
  konvergiert.
- Die Arbeitshypothese ist **falsifizierbar**: Instabilität oder Rekonstruierbarkeit aus $\pi(x)$
  allein widerlegt sie.

**Label:** **Conjecture**.

---

## 6. Epistemischer Kontrast zu §19–§22

| Abschnitt (Übersetzung) | Charakter | Dieses Dokument |
|-------------------------|-----------|-----------------|
| §19 Holographie | Conjecture / Forschungsvision | — |
| §20 Fünfstufenprogramm | Conjecture / Heuristik | — |
| §21 Fossilien-Metapher | Conjecture / Heuristik | — |
| §22 Messreihe / Invarianz-Hypothese | Conjecture (philosophisch) | **Strikte Definitionen 1–5, FF 1–2** |

§22 in `collatz_eabc_bernoulli_uebersetzung.md` liefert die **motivierende Lesart**;
dieses Dokument liefert die **mathematisch prüfbare Formulierung**.

---

---

## 8. EABC-Fluktuationsfeld

**Motivation:** Das Forschungsproblem betrifft die Struktur von Überschüssen und Defiziten zwischen
EABC-Restklassen — nicht die Existenz von $4n+1$-Primzahlen. Dirichlet liefert
$E,A,B,C\sim\pi(x)/4$ asymptotisch (führende Ordnung bekannt); die Information steckt in den
Abweichungen.

### Definition 6 (Fluktuationsvektor $\delta(x)$)

\[
\delta_E(x)=E(x)-\frac{\pi_{>3}(x)}{4},\quad
\delta_A(x)=A(x)-\frac{\pi_{>3}(x)}{4},\ \ldots
\]
\[
\delta(x)=\bigl(\delta_E,\delta_A,\delta_B,\delta_C\bigr)\in\mathbb{R}^4,\qquad
\sum_i \delta_i(x)=0.
\]
Damit liegt $\delta(x)$ in der Hyperebene
\[
\mathcal{E}=\Bigl\{v\in\mathbb{R}^4:\sum_i v_i=0\Bigr\}.
\]

**Label:** **Definition**.

### Definition 7 (Chirale Fluktuationsasymmetrie $\chi_{\mathrm{fluct}}$)

\[
\chi_{\mathrm{fluct}}(x)=(\delta_E+\delta_C)-(\delta_A+\delta_B)=(E+C)-(A+B).
\]
Die bestehende Observable aus §2 erfüllt
\[
\chi(x)=\frac{\chi_{\mathrm{fluct}}(x)}{\pi_{>3}(x)}.
\]

**Label:** **Definition** (vereinheitlichte Benennung mit $\chi$ aus Def.\,4 / Beispiel).

### Definition 8 (EABC-Energie $H(x)$)

\[
H(x)=\|\delta(x)\|_2^2=\sum_{i\in\{E,A,B,C\}}\delta_i(x)^2.
\]
Es gilt $H(x)=0$ genau dann, wenn die vier Klassen perfekt gleich verteilt sind.

**Label:** **Definition**.

### Definition 9 (Kovarianz $K(x)$)

Für eine diskrete Implementierung: **laufende bzw.\ Stichproben-Kovarianz** der Folge
$\bigl(\delta(x)\bigr)_{x\ge 5}$ an Primzähl-Punkten (Gitterschritte $x=5,6,\ldots$).
Die Matrix $K\in\mathbb{R}^{4\times 4}$ schätzt die gemeinsame Variabilität der
Klassenfluktuationen.

**Label:** **Definition** (Experiment: `fluctuation_covariance_at_grid` in der Implementierung).

### Forschungsproblem A (Skalierung von $H$ und $\chi$)

Existieren Grenzwerte, $\limsup$- oder $\liminf$-Größen für
\[
\frac{H(x)}{\pi_{>3}(x)}\qquad\text{oder}\qquad
\frac{\chi(x)}{\sqrt{\pi_{>3}(x)}}?
\]

**Label:** **Forschungsfrage**.

#### Experiment: Skalierung (FF A)

**Label:** **Experiment** (`collatz_eabc_fluktuation_skala_test.py` $\to$ `collatz_eabc_fluktuation_skala.json`).

Exploratives Gitter $x\in\{100,500,10^3,5{\cdot}10^3,10^4,5{\cdot}10^4,10^5,5{\cdot}10^5,10^6\}$;
ein Sieb bis $\max x$, inkrementelle Zählung an den Gitterpunkten ($\mathcal{O}(\pi(\max x))$).

| $x$ | $\pi_{>3}(x)$ | $H/\pi$ | $H/\sqrt\pi$ | $\chi/\sqrt\pi$ | $c_1,c_2,c_3$ |
|-----|---------------|---------|---------------|------------------|---------------|
| $10^2$ | 23 | 0,033 | 0,16 | $-9{,}1\cdot10^{-3}$ | $-0{,}25,-0{,}25,-0{,}25$ |
| $10^3$ | 166 | 0,259 | 3,34 | $-4{,}7\cdot10^{-3}$ | $-1{,}5,-1{,}5,-2{,}5$ |
| $10^4$ | 1227 | 0,056 | 1,96 | $-3{,}0\cdot10^{-4}$ | $-1{,}3,-2{,}3,-3{,}3$ |
| $10^5$ | 9590 | 0,088 | 8,59 | $-5{,}1\cdot10^{-5}$ | $-5{,}5,-6{,}0,-12{,}0$ |
| $10^6$ | 78496 | 0,083 | 23,3 | $-2{,}7\cdot10^{-6}$ | $-8{,}5,-36{,}5,-15{,}0$ |

**Beobachtungen (ehrlich, kein Beweis):**

1. **$H(x)$:** Log-log-Fit $\log H$ vs.\ $\log\pi$ liefert Steigung $\approx 1{,}02$ ($R^2\approx 0{,}96$) —
   führend $H\sim\pi_{>3}(x)$, also $H/\pi$ *asymptotisch* in der Größenordnung konstant ($\approx 0{,}05$–$0{,}10$ bei großem $x$).
   Gleichzeitig schwankt $H/\pi$ auf dem Gitter stark (Minimum $0{,}033$ bei $x=100$, lokales Maximum $0{,}26$ bei $x=10^3$); ein fester Grenzwert ist **nicht** etabliert.

2. **$\chi(x)/\sqrt{\pi}$:** Betrag fällt monoton von $\approx 9\cdot10^{-3}$ auf $\approx 3\cdot10^{-6}$ —
   tendenziell gegen $0$, nicht gegen eine von Null verschiedene Konstante. $\chi_{\mathrm{fluct}}$ wächst langsamer als $\sqrt\pi$ (Steigung $\approx 0{,}53$ in $\log|\chi_{\mathrm{fluct}}|$ vs.\ $\log\pi$).

3. **Moden $c_i$:** $|\!c_i\!|\sim\pi_{>3}^{\alpha}$ mit $\alpha\approx 0{,}45$–$0{,}53$; $c_i/\sqrt\pi$ bleibt klein ($\lesssim 0{,}2$) und oszilliert — keine stabile Normalisierung sichtbar.

4. **Beste einfache Hypothese im Skript-Score** (kleinster Residual auf Log-Steigungen): $H/\pi\approx\mathrm{const}$ —
   muss wegen der genannten Oszillationen als **heuristisch** gelesen werden; kein Ersatz für analytische Dirichlet-/Chebyshev-Information.

**Fazit FF A:** Weder ein sauberer Grenzwert von $H/\pi$ noch von $\chi/\sqrt\pi$ ist numerisch gesichert; $\chi/\sqrt\pi\to 0$ ist die klarste Tendenz. Weitere Punkte $x>10^6$ oder glatte Mittel über $x$-Intervalle wären nötig, um Grenzwert vs.\ Oszillation zu trennen.

### Forschungsproblem B (Spektrum von $K$)

Welche Eigenwerte und Eigenvektoren hat $K$? Lassen sich Eigenrichtungen mit den
EABC-Koordinaten $(E,A,B,C)$ oder mit Primzahl-Statistik interpretieren?

**Label:** **Forschungsfrage**.

### Forschungsproblem C (Fouriermoden auf Klein $V_4$)

Orthogonale Moden (Reihenfolge $E,A,B,C$):
\[
\Phi_0=(1,1,1,1),\quad
\Phi_1=(1,-1,1,-1),\quad
\Phi_2=(1,1,-1,-1),\quad
\Phi_3=(1,-1,-1,1).
\]
Wegen $\sum_i\delta_i=0$ ist $c_0=0$ und
\[
\delta(x)=c_1(x)\Phi_1+c_2(x)\Phi_2+c_3(x)\Phi_3,
\quad
c_i(x)=\frac{\delta(x)\cdot\Phi_i}{\|\Phi_i\|^2}=\frac{\delta(x)\cdot\Phi_i}{4}.
\]

**Label:** **Definition** + **Forschungsfrage** (Statistik der $c_i$).

### Conjecture (EABC-Spektralhypothese)

Die asymptotische Statistik der Modenkoeffizienten $c_i(x)$ enthält ehrliche Information über
$L$-Funktionen modulo $12$ bzw.\ über die Verteilung der $\zeta$-Nullstellen — ohne Behauptung
bereits bewiesener Grenzwerte.

**Label:** **Conjecture**.

**Experiment:** `collatz_eabc_invarianzprogramm.py` $\to$ JSON-Feld `fluctuation_field`.

---

## 9. Implementierung und Querverweise (aktualisiert)

| Artefakt | Rolle |
|----------|-------|
| `collatz_eabc_invarianzprogramm.py` | Numerik: $V$, $S$, $\chi$, $\delta$, $H$, $K$, $c_i$, $\sigma(Q)$ |
| `collatz_eabc_invarianzprogramm.json` | JSON-Output inkl.\ `fluctuation_field` |
| `collatz_eabc_fluktuation_skala_test.py` | FF A: Skalierungstest $H/\pi$, $\chi/\sqrt\pi$, $c_i$ |
| `collatz_eabc_fluktuation_skala.json` | JSON-Output Skalierungsexperiment |
| `tests/test_eabc_fluktuation_skala.py` | Unit-Tests Skalierungsskript |
| `tests/test_eabc_invarianzprogramm.py` | Unit-Tests ($\kappa$, Simplex, Fluktuationsfeld) |
| `eabc_from_lean.py` | Referenz-Implementierung von $\kappa$ |
| `collatz_eabc_holonomie.md` | Begriffshierarchie Stufe 1–6: $\kappa$, $\sigma(Q)$ (hier Def. 1, 5), $\omega$, $\chi_E$; Holonomie = Zielobjekt |
| `collatz_generalangriff_2026.md` | Strategischer Pointer |
| `collatz_eabc_quaternion_mass_hypothese.md` §12 | **EABC-Spektralgeometrische Hauptvermutung** ($D(n)=I(\mu_n)-I_{\mathrm{ref}}(n)$, Emergenz) |
| `collatz_eabc_quaternion_mass_hypothese.md` §13 | **EABC-Spektralgeometrische Erzeugerhypothese** ($\hat D(s)=\sum D(n)/n^s$, Bernoulli-Brücke) |
| `collatz_eabc_shell_defekt_test.py` | Experiment: $I_{\mathrm{ref}}$-Vergleich (rolling, cumulative, $\omega$-, $\tau$-Stratum, $\mu_\infty$) |
| `collatz_eabc_dirichlet_D.py` | Experiment: Partialsummen $\hat D_N(s)$, Vergleich zu $\zeta$ und Bernoulli |
| `collatz_eabc_bernoulli_uebersetzung.md` §17–§22 | Bernoulli-Sensor, $V_n$-Defekt, philosophischer Querverweis (Branch `collatz/eabc-bernoulli-sensor`) |

---

*Epistemische Einordnung: Definitionen 1–9 sind nicht verhandelbar; Forschungsfragen,
Forschungsprobleme A–C und Arbeitshypothese/Spektralhypothese sind explizit offen und falsifizierbar.*
