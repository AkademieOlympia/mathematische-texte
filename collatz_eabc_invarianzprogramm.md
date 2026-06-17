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

## 7. Implementierung und Querverweise

| Artefakt | Rolle |
|----------|-------|
| `collatz_eabc_invarianzprogramm.py` | Numerik: $V(x)$, $S(x)$, $\chi(x)$, $\sigma(Q)$, $\mu$-Schätzung |
| `collatz_eabc_invarianzprogramm.json` | JSON-Output (Experiment) |
| `tests/test_eabc_invarianzprogramm.py` | Unit-Tests ($\kappa$, Simplex, Referenzwerte) |
| `eabc_from_lean.py` | Referenz-Implementierung von $\kappa$ |
| `collatz_generalangriff_2026.md` | Strategischer Pointer |
| `collatz_eabc_bernoulli_uebersetzung.md` §22 | Philosophischer Querverweis hierher |

---

*Epistemische Einordnung: Definitionen 1–5 sind nicht verhandelbar; Forschungsfragen und
Arbeitshypothese sind explizit offen und falsifizierbar.*
