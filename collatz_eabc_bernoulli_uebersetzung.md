# EABC-Resonanzhypothese der Zetafunktion

**Stand:** Juni 2026 · **Epistemische Warnung:** Forschungsnotiz im Tao-Stil
(Definition / Zeuge / Experiment / Theorem / Conjecture / Heuristik / Negativtest).
**Kein Collatz-Beweis. Kein Beweis der Riemann-Hypothese.**

*Hinweis:* Der veraltete Begriff *ERPC* wurde durch **EABC-Zerlegungsprinzip** ersetzt.

**Querverweise:** `Staudt.tex`, `collatz_kepler_gedankenexperiment.tex` (Bernoulli-Uhr),
`eabc_from_lean.py`, `collatz_bernoulli_schalen.pdf`, `energiedoku_eabc_c4_kohaerenz.tex`,
`collatz_eabc_bernoulli_sensor.py`, `CollatzEabc.BernoulliClock.lean`,
`collatz_qed_arithmetik_resonanz.md` (Casimir als $\Delta$-Spektrum, nicht QED-Metapher),
`PAPER_HURWITZ_RESONANZ.md` (Hurwitz-Gitter, Quaternionen $H$; **Theorem**-Hintergrund für §17),
`collatz_kepler_gedankenexperiment.tex` (Kepler-/Projektions-Gedankenexperiment; **Heuristik**).

---

## 1. EABC-Zerlegungsprinzip

**Definition (Zerlegungsprinzip).** Jede natürliche Zahl $N$ besitzt zwei komplementäre Komponenten:
\[
N = (N_{\mathrm{glatt}},\, N_{\mathrm{EABC}}).
\]

| Komponente | Rolle | Typische Größen |
|------------|-------|-----------------|
| $N_{\mathrm{glatt}}$ | **metrisch / skalierend** — kontinuierliche Arithmetikschicht | $\log N$, Primdichte, $\Omega(N)$, Wachstum, analytische Näherungen ($\pi(x)$, $\mathrm{Li}(x)$, $\zeta(s)$, $\Gamma(s)$) |
| $N_{\mathrm{EABC}}$ | **Orientierung im EABC-Raum** — chirale Symmetrieschicht | Prime-Faktorzerlegung nach $E\equiv 1$, $A\equiv 5$, $B\equiv 7$, $C\equiv 11 \pmod{12}$; Zustandsvektor $V(N)=(e,a,b,c)$ mit gewichteten Exponenten pro Klasse |

$N_{\mathrm{glatt}}$ trägt das **analytische Skelett** (Skala, Dichte, glatte Approximationen);
$N_{\mathrm{EABC}}$ trägt die **Orientierung** (Restklassen-Bilanz, chirale Sensoren $\sigma$, $\chi$, $\iota_{\mathrm{chir}}$).
EABC ist damit **Symmetrie-/Orientierungsraum**, nicht bloße Primklassifikation.

---

## 2. EABC-Zustandsraum und Casimir-Spektrum (nicht EM)

Dieser Abschnitt präzisiert das EABC-Zerlegungsprinzip (§1) durch ein **Zustandsraum-Bild**:
Observablen entstehen aus **Asymmetrien erlaubter Konfigurationen**, nicht aus einzelnen Primobjekten.
Die Analogie bezieht sich auf **Casimir-Spektrum und Konfigurationsdichte**, **nicht** auf
virtuelle Teilchen oder die QED-Metapher in `collatz_qed_arithmetik_resonanz.md`.

### 2.1 Casimir-Kern: Observable als Differenz (Heuristik)

**Heuristik (Casimir, nicht Elektrodynamik).** Im klassischen Casimir-Effekt ist die messbare Größe
nicht die absolute Energie im Volumen, sondern die **Differenz** zwischen zwei Randbedingungen:
\[
\Delta E := E_{\mathrm{innen}} - E_{\mathrm{aussen}},
\]
wobei $E_{\mathrm{innen}}$ und $E_{\mathrm{aussen}}$ die (formalen) Energiesummen über erlaubte
Moden unter **verschiedenen** Zustandsraum-Geometrien bezeichnen. Die Kraft folgt aus
$\partial\Delta E/\partial d$ — die **Observable ist die Asymmetrie des Zustandsraums**, nicht ein
einzelner Modus.

**Label: Heuristik** — physikalisches Leitbild, kein arithmetischer Beweis.

### 2.2 Allgemeines Zustandsraumprinzip (Conjecture / Heuristik)

> **Conjecture (Zustandsraumprinzip, allgemein).** Viele arithmetische und physikalische Observablen
> sind **Differenzen in der Dichte erlaubter Konfigurationen** eines zugrunde liegenden Zustandsraums,
> nicht absolute Größen einzelner Elemente.

| Domäne | Zustandsraum (skizzenhaft) | Observable als Differenz |
|--------|---------------------------|--------------------------|
| Thermodynamik | Mikrozustände $\Omega$ | $S = k\log\Omega$ (logarithmische Dichte) |
| Pfadintegrale | Pfade mit Randbedingung | Relativgewicht gegen Referenzmaß |
| Casimir | Modenspektrum im Kasten vs. außen | $\Delta E = E_{\mathrm{innen}}-E_{\mathrm{aussen}}$ |
| Primzahlen | Globale Teilbarkeits-/Restklassenstruktur | Abweichungen von glatten Dichten ($\pi(x)-\mathrm{Li}(x)$, Bias) |

**Label: Conjecture / Heuristik** — Forschungsprogramm, keine abgeschlossene Theorie.

### 2.3 EABC-Version: $Q_4(N)$ und $\Delta Q_4(N)$ (Definition)

**Definition ($Q_4(N)$).** Für $N\ge 2$ sei
\[
Q_4(N) := \bigl(E(N), A(N), B(N), C(N)\bigr),
\]
wobei $E(N),A(N),B(N),C(N)$ die Anzahlen der Primzahlen $p\le N$ in den jeweiligen
EABC-Klassen $p\equiv 1,5,7,11\pmod{12}$ sind (Primzahlen $2,3$ und $p>12$ außerhalb der vier
Klassen werden separat gezählt, analog §7). $Q_4(N)$ ist damit der **aggregierte EABC-Zustandsvektor**
über alle erlaubten Primkonfigurationen bis zur Skala $N$ — nicht die Einzelprimzahl $p$.

**Definition ($\Delta Q_4(N)$).** Die **beobachtbare EABC-Asymmetrie** ist die Differenzstruktur
innerhalb $Q_4(N)$, insbesondere
\[
\Delta Q_4(N) := \bigl(\sigma(N), \chi(N)\bigr),
\qquad
\sigma(N):=(E+A)-(B+C),\quad \chi(N):=(E+B)-(A+C),
\]
sowie die spezifische Kongruenz-Asymmetrie $A(N)-C(N)$ (Diagonalen $A$ vs. $C$ mod $12$).
**Fundamental ist $\Delta Q_4(N)$, nicht $Q_4(N)$ allein** — analog zu $\Delta E$ im Casimir-Fall.

Implementierung: `collatz_eabc_bernoulli_sensor.py` (`q4_vector`, `delta_q4`).

**Label: Definition** — Zähldefinition; physikalische Casimir-Analogie bleibt **Heuristik**.

### 2.4 Casimir-Analogie über Kongruenz-Asymmetrie (Heuristik)

Die EABC-Lesart **vermeidet** das Bild virtueller Teilchen. Stattdessen:

| Casimir (physikalisch) | EABC (arithmetisch) | Label |
|------------------------|---------------------|-------|
| Modenspektrum $\{\omega_n\}$ | Primspektrum / PrimeSig-Filtration $\{P_n\}$ | Definition / Experiment |
| Energiesumme $\sum f(\omega_n)$ | Konfigurationsdichte / Zählvektor $Q_4(N)$, $V_n$ | Definition |
| Kraft aus $\partial\Delta E/\partial d$ | **Bias** aus $\Delta Q_4(N)$, $\sigma$, $\chi$ | Conjecture / Experiment |
| Randbedingung „innen“ vs. „außen“ | Kongruenz-Asymmetrie $A$ vs. $C$ mod $12$ | Definition |

Kette: **Spektrum $\to$ Energie $\to$ Kraft** (Casimir) entspricht heuristisch
**Primspektrum $\to$ Konfigurationsdichte $\to$ Bias** (EABC).

**Label: Heuristik** — strukturelle Analogie, kein Theorem über Primzahlen.

### 2.5 Zeta-Verbindung (Conjecture)

> **Conjecture (Zeta-Spektrum).** Nichttriviale Nullstellen $\rho_n=\tfrac12+\mathrm{i}t_n$ der
> Riemannschen Zetafunktion wirken — analog zu erlaubten Casimir-Eigenmoden — als **globales
> Spektrum**, das die Dichte arithmetischer Observablen mitprägt. Die explizite Formel
> \[
> \pi(x) \;\longleftrightarrow\; \{\rho_n\}
> \]
> zeigt: nicht einzelne Nullstellen, sondern die **gesamte spektrale Struktur** zählt.

Diese Conjecture verknüpft §2 mit der **EABC-Resonanzhypothese** (§8): Bernoulli-Stufe $n$ liefert
$V_n$ entlang trivialer Nullstellen; nichttriviale $\rho_k$ sind die postulierte spektrale Projektion.

**Label: Conjecture** — klassische explizite Formel ist **Theorem**; EABC-Deutung ist **Conjecture**.

### 2.6 Boxed: EABC-Zustandsraum-Hypothese

$$\boxed{
\text{Arithmetische Bias-Strukturen entstehen aus Asymmetrien des EABC-Zustandsraums.}
}$$

**Lesart:** Primzahlen sind **beobachtbare Ereignisse** (Realisierungen erlaubter Konfigurationen);
der **Zustandsraum** (EABC-Klassen, $Q_4(N)$, PrimeSig-Filtrationen) ist **fundamental**.
Nicht „Prim erzeugt EABC“, sondern: **EABC bestimmt, welche Primkonfigurationen möglich sind.**

**Label: Conjecture** — Forschungshypothese, epistemisch ehrlich nicht als etablierte Mathematik.

### 2.7 Boxed: Bernoulli-Kette (Übersetzer)

$$\boxed{
\text{triviale Nullstellen}
\;\to\;
\text{Bernoulli}
\;\to\;
\text{EABC-Konfigurationsraum}
\;\to\;
\text{nichttriviales Spektrum}
\;\to\;
\text{Primzahl-Bias}
}$$

Bernoulli-Zahlen sind **Übersetzer** von der äquidistanten Achse $s=-2,-4,-6,\ldots$ ($N_{\mathrm{glatt}}$)
in den **asymmetrischen** EABC-Zustandsraum ($N_{\mathrm{EABC}}$); vgl. §5 (Bernoulli-Brücke) und §7 ($V_n$).

**Label: Definition / Heuristik** — Bernoulli-Staudt ist **Theorem**; die Kette als Ganzes **Conjecture**.

### 2.8 Philosophischer Kern (Conjecture)

> **Conjecture (Forschungsprogramm).** Die gemeinsame Struktur zwischen Casimir-Physik und
> EABC-Arithmetik ist nicht „alles ist ähnlich“, sondern: **Geometrie des Konfigurationsraums +
> Modenselektion** erzeugt Observable als **Differenzen**. Primzahlen sind nicht die erzeugende
> Ursache der EABC-Struktur — EABC legt fest, **welche** Primkonfigurationen zulässig sind und
> welche Bias-Observablen $\Delta Q_4(N)$ messbar werden.

**Label: Conjecture** — philosophische Präzisierung, strenger als bloße Analogie.

---

## 3. Grundhypothese: glatte und EABC-Dynamik

> **Conjecture (Grundhypothese).** Klassische analytische Größen ($\pi(x)$, $\zeta(s)$, $\Gamma(s)$, $B_{2n}$)
> stammen überwiegend aus $N_{\mathrm{glatt}}$; Fluktuationen, Resonanzen und Bias aus $N_{\mathrm{EABC}}$:
> \[
> \text{Arithmetik} = \text{glatte Dynamik} + \text{EABC-Dynamik}.
> \]

**Label: Conjecture** — heuristische Leitlinie, nicht als Theorem behauptet.
Der Bernoulli-/Zeta-Zweig in diesem Dokument ist die **EABC-Deutung** trivialer Zetawerte
und von-Staudt-Clausen-Signaturen entlang dieser Zerlegung.

---

## 4. Zentrale Fragestellung

Die **EABC-Resonanzhypothese der Zetafunktion** fragt, ob die Folge der
EABC-Zustandsvektoren $V_1,V_2,\ldots$ an der äquidistanten Achse trivialer Nullstellen
eine **diskrete Resonanzstruktur der Primzahlen** trägt — und ob nichttriviale Nullstellen
als **spektrale Projektion** dieser Resonanz lesbar sind.

| Ebene | Label in diesem Text |
|-------|----------------------|
| $N=(N_{\mathrm{glatt}},N_{\mathrm{EABC}})$, Zerlegungsprinzip | **Definition** |
| Casimir-$\Delta E$, Zustandsraumprinzip | **Heuristik / Conjecture** |
| $Q_4(N)$, $\Delta Q_4(N)$, EABC-Zustandsraum-Hypothese | **Definition / Conjecture** |
| Arithmetik = glatte + EABC-Dynamik | **Conjecture** |
| EABC mod $12$, $V_n=(E_n,A_n,B_n,C_n)$ | **Definition** |
| von Staudt--Clausen, $\mathrm{PrimeSig}(B_{2n})$ | **Definition** (klassisches Theorem) |
| Sensor $\Phi$, chirale Observablen | **Definition / Zeuge** |
| Test 1 (EABC-Asymmetrien in $V_n$) | **Experiment** |
| Resonanzstruktur auf der trivialen Nullstellenachse | **Conjecture** |
| Abbildung $V_n \to \rho_k$ bzw. $\Delta t_k$ | **Conjecture** (starke Form) |
| Bernoulli-Lyapunov für Collatz | **Negativtest** (No-Go) |
| Hurwitz: normierte Divisionsalgebren $\mathbb{R},\mathbb{C},\mathbb{H},\mathbb{O}$ | **Theorem** (§17) |
| Peano-Arithmetik, $S(n)=n+1$ | **Definition / Theorem** (§17) |
| Peano als Projektion tieferer Defektdynamik | **Conjecture / Heuristik** (§17) |
| EABC-Tetraeder, Kepler-Füllung, Prim-Defekte | **Forschungsvision** (§17) |
| $\mathcal{K}(N)$, $D(N)$, $\Pi$, $\Phi_{\mathrm{def}}$ | **Conjecture** (§17, offen) |

---

## 5. Äquidistante Basis und Bernoulli-Brücke

**Definition.** Die trivialen Nullstellen der Riemannschen Zetafunktion liegen bei
\[
s = -2,\,-4,\,-6,\,\ldots,\qquad s_n = -2n\quad (n\ge 1),
\]
einer **flachen, äquidistanten Basisschicht** auf der negativen reellen Achse.

Jede Stufe $n$ kodiert einen rationalen Spezialwert:
\[
\zeta(1-2n) = -\frac{B_{2n}}{2n},
\]
wobei $B_{2n}$ die geraden Bernoulli-Zahlen sind. Bernoulli-Zahlen sind **Übersetzungsobjekte**:
sie überführen die äquidistante Struktur trivialer Nullstellen ($N_{\mathrm{glatt}}$)
in die EABC-Schicht ($N_{\mathrm{EABC}}$) — nicht Endobjekte, sondern Operatoren
von glatten Gitterpunkten $s=-2n$ zu rationalen Arithmetikdaten und Primsignaturen.

> **Boxed (Präzision):** Bernoulli $\neq$ Prim. Die Korrespondenz verläuft zwischen
> **Filtrationen** $\{B_{2n}\}_{n\ge 1} \longleftrightarrow \{p\ \text{prim}: p-1\mid 2n\}_{n\ge 1}$.

---

## 6. von Staudt--Clausen und Primsignaturen

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

## 7. Native EABC mod $12$ und Zustandsvektor $V_n$

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

## 8. EABC-Resonanzhypothese (Conjecture)

> **Conjecture (Resonanzhypothese, EABC-Version).**
> - **Triviale Nullstellen** $s=-2n$: äquidistante Ausgangsstruktur ($N_{\mathrm{glatt}}$).
> - **Bernoulli** $B_{2n}$: liest Primresonanzen über $\mathrm{PrimeSig}(B_{2n})=P_n$.
> - **EABC-Zerlegung** $V_n=(E_n,A_n,B_n,C_n)$: offenbart Orientierung der Resonanzen ($N_{\mathrm{EABC}}$).
> - **Nichttriviale Nullstellen** $\rho_k=\tfrac12+\mathrm{i}\,t_k$: **spektrale Projektion**
>   der EABC-durchdrungenen Struktur auf die kritische Linie.

Die Folge $V_1,V_2,V_3,\ldots$ ist **nicht zufällig**: sie trägt diskrete Primresonanz
entlang der trivialen Achse; Abweichungen in $V_n$ kodieren EABC-Dynamik, die im
Spektrum der nichttrivialen Nullstellen wieder auftaucht.

**Label: Conjecture** — empirisch falsifizierbar (Abschnitt 11), nicht als Theorem behauptet.

Die Resonanzhypothese setzt das **Zustandsraum-Bild** (§2) fort: $V_n$ ist die PrimeSig-Filtration
entlang trivialer Nullstellen; $\Delta Q_4(N)$ und $\sigma,\chi$ sind die aggregierte Asymmetrie
des EABC-Konfigurationsraums.

---

## 9. Starke Form: Abbildung auf das Spektrum (Conjecture)

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

## 10. Geometrische Interpretation (Heuristik)

Fünfgliedrige Lesart als Projektionskette (EABC-Zerlegungsprinzip, vgl. §2 Zustandsraum):

0. **Konfigurationsraum:** $Q_4(N)$, Observable $\Delta Q_4(N)$ — Asymmetrie erlaubter EABC-Konfigurationen.
1. **Äquidistante Basis ($N_{\mathrm{glatt}}$):** Triviale Nullstellen $s=-2n$ — gleichabständiges Gitter auf $\mathrm{Re}(s)<0$.
2. **Bernoulli-Übersetzung:** $\zeta(1-2n)=-B_{2n}/(2n)$ — Operator von glatten Gitterpunkten zu rationalen Daten.
3. **Primresonanz:** $P_n=\mathrm{PrimeSig}(B_{2n})$ — diskrete Einbettung der Primstruktur.
4. **EABC-Orientierung ($N_{\mathrm{EABC}}$):** $V_n=(E_n,A_n,B_n,C_n)$ — Symmetrie-/Orientierungsraum mod $12$ ($\sigma$, $\chi$, $\iota_{\mathrm{chir}}$).
5. **Spektrale Projektion:** Nichttriviale Nullstellen $\rho_k=\tfrac12+\mathrm{i}t_k$ — Spektrum der EABC-durchdrungenen Struktur auf $\mathrm{Re}(s)=\tfrac12$; RH als Fixierung dieser Projektion.

**Label: Heuristik / Conjecture** — geometrische Metapher, kein Beweisanspruch.

---

## 11. Falsifizierbarkeit und Experimente

**Falsifikationskriterium.** Berechne $V_n$ für viele $n$, vergleiche mit bekannten
Nullstellenabständen $\Delta t_k$. **Keine** systematische Korrelation $\Rightarrow$
Resonanzhypothese **falsch** (für die getestete Kopplung). **Systematische** Struktur
$\Rightarrow$ es existiert eine arithmetisch–spektrale Brücke — weiter zu präzisieren.

| # | Test | Status | Label |
|---|------|--------|-------|
| **1** | EABC-Asymmetrien in $V_n$ ($\sigma$, $\chi$, $\iota_{\mathrm{chir}}$) | **implementiert** | Experiment |
| **1b** | Lean-Kopplung: $V_n$ vs. `class_of` (formale EABC-Schicht) | **implementiert** | Experiment |
| **1c** | $\Delta Q_4(N)$ vs. $\pi(N)-\mathrm{Li}(N)$ (Zustandsraum-Bias) | **zukünftig** | Experiment |
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

## 12. Nomenklatur

| Begriff | Bedeutung in diesem Zweig |
|---------|---------------------------|
| **EABC** | Symmetrie-/Orientierungsraum (mod $12$) |
| **Zustandsraum** | Menge erlaubter EABC-Konfigurationen; Observable = $\Delta Q_4$ |
| **$Q_4(N)$** | Aggregierter EABC-Zählvektor über Primzahlen $p\le N$ |
| **$\Delta Q_4(N)$** | Chirale Asymmetrie $(\sigma,\chi)$, $A-C$-Bias |
| **glatt** ($N_{\mathrm{glatt}}$) | Metrische / analytische Schicht |
| **Bernoulli** | Übersetzungsobjekte (nicht Endobjekte) |
| **triviale Nullstellen** | Äquidistante Ausgangsstruktur |
| **nichttriviale Nullstellen** | Spektrum der EABC-durchdrungenen Struktur |

---

## 13. Boxed Kurzform

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

## 14. Parallele Tracks (nicht verschmolzen)

Der EABC-Bernoulli-/Resonanz-Zweig läuft **parallel** zu anderen Sensoren:

| Track | Sensor | Artefakt |
|-------|--------|----------|
| **EABC-Resonanz** | $\Phi$, $V_n$ | `collatz_eabc_bernoulli_sensor.json` |
| **Morley** | $F_M$, $S_M$, $K_M$ | `collatz_morley_tm_numerik.py` |
| **$\kappa$ / Grammatik** | $\mathcal{L}_{\mathrm{arith}}$, $F_n$ | `collatz_forbidden_words.py` |

Morley und $\kappa$ bleiben eigenständige Angriffslinien. Der Resonanz-Zweig nutzt
**dieselbe** EABC-Sprache, aber andere Eingangsgrößen (von Staudt statt Collatz-Bahnen).

---

## 15. Tao-Labels und nächste Schritte

| Aussage | Tao-Label |
|---------|-----------|
| $N=(N_{\mathrm{glatt}},N_{\mathrm{EABC}})$, Zerlegungsprinzip | Definition |
| Casimir-$\Delta E$, allg. Zustandsraumprinzip | Heuristik / Conjecture |
| $Q_4(N)$, $\Delta Q_4(N)$ | Definition |
| EABC-Zustandsraum-Hypothese, Bernoulli-Kette (§2) | Conjecture / Heuristik |
| Arithmetik = glatte + EABC-Dynamik | Conjecture |
| EABC mod $12$, $V_n=(E_n,A_n,B_n,C_n)$ | Definition |
| PrimeSig aus von Staudt--Clausen | Definition (Theorem) |
| Bernoulli als Übersetzungsobjekt | Definition |
| $\Phi$, $\sigma$, $\chi$, $\iota_{\mathrm{chir}}$ | Definition + Experiment |
| Resonanzstruktur $V_1,V_2,\ldots$ | Conjecture |
| $\Phi:V_n\to\rho_k$ bzw. $\Delta t_k$ | Conjecture |
| RH als Fixraum-Projektion | Conjecture / Heuristik |
| Bernoulli-Lyapunov für Collatz | Negativtest |
| Hurwitz-Theorem (Divisionsalgebren) | Theorem |
| Peano-Arithmetik | Definition / Theorem |
| Peano $\to$ Defektprojektion, Tetraeder, Kepler-Füllung | Conjecture / Heuristik / Forschungsvision |
| $\mathcal{K}(N)$, $D(N)$, $\Pi$, $\Phi_{\mathrm{def}}$ | Conjecture (offen) |

**Nächste Schritte:**
1. Test 1 auf größeres $n$ und mit expliziten Zufalls-Nullmodellen vergleichen.
2. Test 2: $\Delta t_k$-Daten einbinden, Korrelationen und Falsifikation dokumentieren.
3. Lean-Schnittstelle: PrimeSig als `Finset` über `Nat.Prime` (analog `BernoulliClock.lean`).

---

## 16. Experiment: Lean-Kopplung (Tao Experiment)

> **Hinweis:** „LEA-M“ im Gespräch war eine Fehlhörung von **Lean** — gemeint ist die
> formale EABC-Schicht (`EABC.lean` → `eabc_from_lean.py`), kein separates Modul.

**Label: Experiment** — kein Theorem, sondern **Konsistenz- und Strukturcheck** zwischen dem
analytischen Bernoulli-Sensor und der formalen EABC-Schicht in Lean.

### Motivation

Der Sensor `collatz_eabc_bernoulli_sensor.py` zählt $V_n=(E_n,A_n,B_n,C_n)$ über
$\mathrm{PrimeSig}(B_{2n})$ und nutzt dabei bereits `class_of` aus `eabc_from_lean.py`.
Dieses Experiment **verdoppelt** den Zählweg unabhängig: dieselben Primzahlen $p\in P_n$
werden erneut nur über die Lean-Spiegelung klassifiziert; Quadrupel- und Chiralitätszeuge
aus `EABC.lean` werden auf $P_n$ angewendet.

**Epistemik:** Übereinstimmung beweist weder die Resonanzhypothese noch RH — sie zeigt nur,
dass Python-Sensorik und Lean-Definitionen dieselbe mod-$12$-Semantik tragen.

### Lean-Module (formale Quelle)

| Modul | Inhalt |
|-------|--------|
| `EABC.lean` | `EClass`, `residue`, `classOf`, `Q`, `IsPrimeQuadruplet`, `Chirality`, `T`, `T4` |
| `CollatzEabc.Mod12Matrix` | `EabcIndex`, Restklassen $1,5,7,11$ |
| `CollatzEabc.PrefProjection` | `EabcLetter`, $\Phi_{\mathrm{pref}}$, Radial $\rho$ |
| `CollatzEabc.BernoulliClock` | Bernoulli-Tripel, Dreierphase (Gedankenmodell) |
| `eabc_from_lean.py` | Python-Spiegelung von `EABC.lean` |

### Prüfpunkte

Für $n=1,\ldots,N$:

1. **Zählvektor:** $V_n^{\mathrm{sensor}} = V_n^{\mathrm{lean}}$ (unabhängige `class_of`-Rechnung).
2. **Residuum:** für jedes $p\in P_n$ mit Klasse $X$: $p \equiv \mathrm{residue}(X) \pmod{12}$.
3. **Rotation:** $T^4(X)=X$ auf allen vorkommenden Klassen (Lean-Theorem `T4_has_order_4`).
4. **Vierlingszeuge:** $Q(p)=[p,p+2,p+6,p+8]$ in $P_n$; Chiralität ABCE bzw. CEAB falls zutreffend.
5. **Nullmodell:** Zufällige EABC-Verteilung bei fester $|P_n\cap\mathrm{EABC}|$; $z$-Abweichung von $\sigma$, $\chi$.

### Artefakte

| Datei | Rolle |
|-------|-------|
| `collatz_eabc_bernoulli_lean_test.py` | Experiment-Runner, JSON-Report |
| `tests/test_eabc_bernoulli_lean.py` | pytest (Konsistenz bis $n\le 500$) |
| `collatz_eabc_bernoulli_lean.json` | Standard-Ausgabe (`--output`) |

Ausführung:
```bash
python3 collatz_eabc_bernoulli_lean_test.py --max-n 200
pytest tests/test_eabc_bernoulli_lean.py -q
```

**Erwartung:** `summary.classification_match_all`, `residue_roundtrip_all`, `T4_identity_all`
sind `true` für alle getesteten $n$. Abweichung wäre ein **Implementierungsfehler**, nicht ein
Gegenbeispiel zur Resonanzhypothese.

---

## 17. EABC-Forschungsvision: Peano, Tetraeder, Kepler und Prim-Defekte

Dieser Abschnitt formuliert die **übergeordnete Forschungsvision** hinter §1–§2 und der
Bernoulli-/Zeta-Kette — epistemisch strikt getrennt von etablierter Mathematik.
Er ersetzt weder §1 (Zerlegungsprinzip) noch §2 (Zustandsraum) noch die Bernoulli-Abschnitte §5–§9;
er **vertieft** deren philosophische Lesart.

**Querverweise:** `PAPER_HURWITZ_RESONANZ.md` (Quaternionen $H$, Hurwitz-Maximalordnung),
`collatz_kepler_gedankenexperiment.tex`, `CollatzEabc.Mod12Matrix` / `EABC.lean` (formale
$E,A,B,C$-Semantik mod $12$).

### 17.1 Etablierte Mathematik (Theorem / Definition)

**Theorem (Hurwitz, normierte Divisionsalgebren).** Über $\mathbb{R}$ existieren genau vier
endlich-dimensionale assoziative normierte Divisionsalgebren:
\[
\mathbb{R}\;(1\text{D}),\quad
\mathbb{C}\;(2\text{D}),\quad
\mathbb{H}\;(4\text{D}),\quad
\mathbb{O}\;(8\text{D}),
\]
mit Dimensionen $1,2,4,8$. Es gibt **keine** weiteren endlich-dimensionalen normierten
Divisionsalgebren über $\mathbb{R}$. Die Kette endet bei den Oktanionen $\mathbb{O}$;
Assoziativität geht beim Übergang $\mathbb{H}\to\mathbb{O}$ verloren.

**Label: Theorem** — klassische Hurwitz-Theorie; keine EABC-Spezialisierung.

**Definition / Theorem (Peano-Arithmetik).** Die natürlichen Zahlen $\mathbb{N}=\{0,1,2,\ldots\}$
werden durch die Peano-Axiome charakterisiert; die Nachfolgerabbildung
\[
S:\mathbb{N}\to\mathbb{N},\qquad S(n)=n+1
\]
erzeugt die **eindimensionale lineare Dynamik** auf der Zahlengeraden: jeder Schritt verschiebt
um genau eine Einheit.

**Label: Definition / Theorem** — Standardarithmetik; $n\mapsto n+1$ ist die kanonische
1D-Dynamik.

### 17.2 Peano als Projektion (Conjecture / Heuristik)

> **Conjecture / Heuristik (Peano-Projektion).** Die Peano-Dynamik $S(n)=n+1$ ist nicht die
> fundamentale Dynamik der Arithmetik, sondern die **Projektion** einer tieferen Defektdynamik
> auf eine eindimensionale Achse.

In dieser Lesart erzeugt $n\mapsto n+1$ die sichtbare **Zahlengerade** $\mathbb{N}$.
**Primzahlen** erscheinen dort, wo die Projektion **innere Struktur verliert**: Stellen, an denen
die volle mehrkanalige Defektkonfiguration auf einen einzelnen Skalenwert kollabiert, ohne dass
die zugrunde liegende Symmetrie vollständig aufgelöst werden kann.

**Label: Conjecture / Heuristik** — Forschungsbild, kein Theorem über Primzahlen.

### 17.3 Vom 1D-Achse zum EABC-Tetraeder (Conjecture)

> **Conjecture (Tetraeder-Dynamik).** Die fundamentale Dynamik liegt auf einem **Tetraeder**
> mit Eckpunkten $(E,A,B,C)$ — optimaler Defektverteilung in vier Kanälen — statt auf der
> eindimensionalen Peano-Achse.

| Ebene | Struktur | Label |
|-------|----------|-------|
| Peano | 1 Kanal, $S(n)=n+1$ | Theorem (Projektionsziel) |
| EABC | 4 Kanäle $(E,A,B,C)$ mod $12$ | Definition (§1, §7) |
| Tetraeder | geometrische Trägerform der 4-Kanal-Dynamik | Conjecture |

Der Übergang **1-dimensional $\to$ 4-Kanal** ist der zentrale Schritt der Forschungsvision:
$\mathbb{N}$ ist die **Projektion**, $(E,A,B,C)$ die **volle Orientierung** (vgl. §1 Zerlegungsprinzip).

**Label: Conjecture** — geometrisches Postulat, nicht etablierte Zahlentheorie.

### 17.4 Oktanionen: Hurwitz-Kette und Spekulation (Theorem + Heuristik)

**Theorem (Hurwitz-Kette).** Die normierte Divisionsalgebren-Kette
\[
\mathbb{R}\to\mathbb{C}\to\mathbb{H}\to\mathbb{O}
\]
endet bei Dimension $8$; $\mathbb{O}$ ist die letzte endlich-dimensionale normierte Divisionsalgebra
über $\mathbb{R}$.

**Etablierte Strukturen (Definition / Theorem):**
- $\mathbb{H}$: assoziative Quaternionen, Fano-Ebene als Multiplikationstabelle der imaginären Einheiten;
- $\mathbb{O}$: nichtassoziativ; Automorphismengruppe $G_2=\mathrm{Aut}(\mathbb{O})$;
- $E_8$-Gitter, Kugelpackungen in Dimension $8$ — klassische Objekte der Geometrie und Lie-Theorie.

**Heuristik (Spekulation).** Die Hurwitz-Kette und ihre $8$-dimensionalen Nachbarschaften
($G_2$, Fano-Ebene, $E_8$, Gitterpackungen) könnten — **spekulativ** — mit einer tieferen
Primstruktur verknüpft sein, die in der 1D-Projektion nur als „Primdefekte" sichtbar wird.
Verbindung zu `PAPER_HURWITZ_RESONANZ.md` (Quaternionen-Basis $E\leftrightarrow 1$, $A\leftrightarrow i$,
$B\leftrightarrow j$, $C\leftrightarrow k$).

**Label: Theorem** für Hurwitz; **Heuristik** für jede Primzahl-Verknüpfung.

### 17.5 Kepler-Füllmechanismus (Forschungsvision)

> **Forschungsvision (Kepler-Füllung, nicht Astronomie).** „Kepler" bezeichnet hier **kein**
> astronomisches Modell, sondern **geometrisch optimale lokale Packung** — analog zu kristalliner
> Ordnung mit Defekten.

| Kristall-Defekt | Arithmetische Lesart (Vision) |
|-----------------|-------------------------------|
| Versetzung (dislocation) | lokale Prim-Störung |
| Leerstelle (vacancy) | fehlende Füllung auf der Projektionsachse |
| Korngrenze (grain boundary) | Übergang zwischen EABC-Konfigurationen |

**Primzahlen** sind in dieser Vision **Defekte** in einer optimalen Füllung: Die ideale
**Tetraeder-Füllung** würde vollständige lokale Schließung ergeben; Primzahlen markieren Stellen,
an denen die Schließung **unvollständig** bleibt — sichtbar erst nach Projektion $\Pi$ auf $\mathbb{N}$.

**Label: Forschungsvision** — physikalisch-geometrische Metapher, kein Beweis.

### 17.6 Mathematische Form: wachsender Zustandsraum und Defektvektor (Conjecture)

> **Conjecture (Kepler-Zustandsraum, offen).** Es existiert ein mit $N$ wachsender Zustandsraum
> $\mathcal{K}(N)$ („Kepler-Füllraum") und ein Defektvektor
> \[
> D(N) = (E(N), A(N), B(N), C(N)),
> \]
> sowie eine Projektion
> \[
> \Pi : \mathcal{K}(N) \longrightarrow \mathbb{N},
> \]
> die die natürlichen Zahlen als 1D-Bild erzeugt. Für $p\in\mathbb{N}$:
> \[
> p\ \text{prim}
> \quad\Longleftrightarrow\quad
> D(p) \in \mathcal{D}_{\mathrm{krit}},
> \]
> wobei $\mathcal{D}_{\mathrm{krit}}$ eine **noch unbekannte kritische Defektmenge** ist.

**Verbindung zu §2:** $D(N)$ stimmt heuristisch mit $Q_4(N)$ und $V_n$ überein, sobald
$\mathcal{K}(N)$ durch EABC-Konfigurationen und PrimeSig-Filtrationen operationalisiert wird.
Die formale Definition von $\mathcal{K}(N)$ und $\mathcal{D}_{\mathrm{krit}}$ ist **offen**.

**Label: Conjecture** — Forschungsprogramm, keine etablierte Charakterisierung der Primzahlen.

### 17.7 Offene Abbildung $\Phi_{\mathrm{def}}$ (Research program)

> **Research program (Defekt $\to$ Primzahl auf der Projektion).** Gesucht ist eine Abbildung
> \[
> \Phi_{\mathrm{def}} : D(N) \longrightarrow \pi(N)
> \qquad\text{oder}\qquad
> \Phi_{\mathrm{def}} : D(N) \longrightarrow \mathbf{1}_{\mathrm{prim}}(N),
> \]
> die den Defektzustand $D(N)$ auf Primzahl-Ereignisse entlang der Projektionsachse $\mathbb{N}$
> abbildet.

Dies ist **stärker** als der Bernoulli-Sensor $\Phi(n)=V_n$ (§7): $\Phi_{\mathrm{def}}$ soll
direkt die **Primheit** aus der Defektkonfiguration lesen, nicht nur PrimeSig-Aggregate entlang
trivialer Nullstellen. Keine konkrete Form von $\Phi_{\mathrm{def}}$ ist bekannt.

**Label: Research program / Conjecture** — offene Funktion, nicht implementiert.

### 17.8 Zeta-Verbindung und Bernoulli-Kette (Conjecture)

Verknüpfung mit §2.5 und §2.7 (Bernoulli-Kette):

| Zeta-Ebene | Defekt-Lesart (Vision) | Label |
|------------|------------------------|-------|
| Triviale Nullstellen $s=-2n$ | Idealer, äquidistanter Füllprozess ($N_{\mathrm{glatt}}$) | Definition |
| Bernoulli $B_{2n}$ | Übersetzer von glatter Füllung zu arithmetischen Defekten | Theorem / Heuristik |
| Nichttriviale Nullstellen $\rho_k$ | Spektrallinien des **gesamten** Defektsystems | Conjecture |

$$\boxed{
\text{Peano zählt Schichten;}
\quad
\text{EABC-Tetraeder organisiert Defekte;}
\quad
\text{Primzahlen markieren sichtbare Defektstellen;}
\quad
\text{Zeta misst das globale Spektrum.}
}$$

**Label: Conjecture** für die Defekt-Deutung nichttrivialer Nullstellen; Bernoulli-Staudt und
explizite Formel bleiben **Theorem** in ihrer klassischen Form.

### 17.9 Schluss: fundamentale Objekte (Forschungsvision)

> **Forschungsvision (epistemische Zusammenfassung).** Die fundamentalen Objekte sind
> **Defektkonfigurationen** in einem wachsenden EABC-Tetraeder-Raum $\mathcal{K}(N)$;
> die eindimensionale Projektion $\Pi$ erscheint als natürliche Zahlen $\mathbb{N}$.
> Primzahlen sind die **sichtbaren Defektstellen** dieser Projektion.

Dies ist **keine etablierte Theorie** — weder ein Beweis der Primzahlverteilung noch der
Riemann-Hypothese. Es ist eine **klarere Forschungsvision**, die §1 (Zerlegungsprinzip),
§2 (Zustandsraum, Casimir-$\Delta$), die Bernoulli-Brücke (§5–§9) und die Hurwitz-Kette
(§17.1, §17.4) in einem gemeinsamen Bild verbindet.

**Lean-Bezug (Definition, nicht Vision):** `CollatzEabc.Mod12Matrix` (`EabcIndex`, Restklassen
$1,5,7,11$), `EABC.lean` (`EClass`, `classOf`, `Q`, `Chirality`) formalisieren die **vier Kanäle**
$(E,A,B,C)$ — nicht $\mathcal{K}(N)$ oder $\mathcal{D}_{\mathrm{krit}}$.

**Label: Forschungsvision** — explizit von Theorem-Ebene getrennt.

---

*Epistemische Einordnung gemischt je Abschnitt. Definitionen und Experimente sind streng;
Conjectures sind explizit falsifizierbar. Keine Ebene ersetzt die nächsthöhere.*
