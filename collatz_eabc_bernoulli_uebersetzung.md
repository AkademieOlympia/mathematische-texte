# EABC-Theorie: Äquidistanz, Zeta-Krümmung und Bernoulli-Übersetzung

**Stand:** Juni 2026 · **Epistemische Warnung:** Forschungsnotiz im Tao-Stil
(Definition / Zeuge / Experiment / Theorem / Conjecture / Heuristik / Negativtest).
**Kein Collatz-Beweis. Kein Beweis der Riemann-Hypothese.**

**Terminologie:** *ERPC* bezeichnet in diesem Repo **dieselbe EABC-Struktur**
(Restklassen $E,A,B,C$ mod $12$, chirale Sensoren $\sigma$, $\chi$, $\iota_{\mathrm{chir}}$).
Es gibt kein separates „ERPC-Projekt neben EABC“ — die Bernoulli-/Zeta-/Äquidistanz-Theorie
ist die **EABC-Deutung** trivialer Zetawerte und von-Staudt-Clausen-Signaturen.

**Querverweise:** `Staudt.tex`, `collatz_kepler_gedankenexperiment.tex` (Bernoulli-Uhr),
`eabc_from_lean.py`, `collatz_bernoulli_schalen.pdf`, `energiedoku_eabc_c4_kohaerenz.tex`,
`collatz_eabc_bernoulli_sensor.py`, `CollatzEabc.BernoulliClock.lean`.

---

## 1. Einleitung und epistemischer Rahmen

Die EABC-Theorie ordnet natürliche Zahlen nach vier Restklassen modulo $12$:
$E\equiv 1$, $A\equiv 5$, $B\equiv 7$, $C\equiv 11$.
Diese Klassifikation ist **nativ EABC** (`EABC.lean`, `eabc_from_lean.py`) — kein externes
Anhängsel, sondern der arithmetische Kern des gesamten Collatz-Generalangriffs.

Dieses Dokument beschreibt einen **zweiten Sensorzweig innerhalb EABC**: die Übersetzung
äquidistanter Zeta-Geometrie (triviale Nullstellen, negative gerade Spezialwerte) in
Primzahlsignaturen über Bernoulli-Zahlen. Methodik: Tao-Stil (`collatz_formalisierung_tao_stil.md`).

| Ebene | Label in diesem Text |
|-------|----------------------|
| von Staudt--Clausen, mod-$12$-EABC | **Definition** |
| PrimeSig, $V(B_{2n})$, Sensor $\Phi_B$ | **Definition / Zeuge** |
| Test 1 (Bernoulli-Signatur, EABC-Asymmetrien) | **Experiment** |
| Kernhypothese (Bernoulli übersetzt Zeta $\to$ Primarithmetik) | **Conjecture** |
| RH als Projektion auf $\mathrm{Re}(s)=\tfrac12$ | **Conjecture / Heuristik** |
| Tests 2–3 (Zeta-Kopplung, Krümmung vs. $\pi(x)-\mathrm{Li}(x)$) | **zukünftig** |

---

## 2. Kritische Präzision: keine naive 1:1-Zuordnung

**Falsch (naiv):** $B_{2n}$ „ist“ eine einzelne Primzahl oder eine Bijektion Bernoulli $\leftrightarrow$ Prim.

**Richtig (EABC-Deutung):** Die Korrespondenz verläuft zwischen **Filtrationen**:
\[
\{B_{2n}\}_{n\ge 1} \;\longleftrightarrow\; \bigl\{p\ \text{Primzahl} : p-1\mid 2n\bigr\}_{n\ge 1}.
\]
Einzelne Bernoulli-Zahlen sind **Übersetzungsobjekte**, in denen Primzahlen über
Teilbarkeit, Nenner und Kongruenz **sichtbar** werden — nicht als identische Objekte.

> **Boxed (Präzision):** Bernoulli $\neq$ Prim. Filtration $\{B_{2n}\}$ $\leftrightarrow$
> Filtration $\{p : p-1\mid 2n\}$.

---

## 3. von Staudt--Clausen (Definition)

Für gerade Bernoulli-Zahlen $B_{2n}$ ($n\ge 1$) gilt der Satz von von Staudt--Clausen:
\[
B_{2n} + \sum_{\substack{p\ \text{prim}\\ p-1\mid 2n}} \frac{1}{p} \;\in\; \mathbb{Z},
\qquad
\mathrm{den}(B_{2n}) = \prod_{\substack{p\ \text{prim}\\ p-1\mid 2n}} p.
\]
Die **Primzahlsignatur** (PrimeSig) ist
\[
\mathrm{PrimeSig}(B_{2n}) := \{p\ \text{prim} : p-1\mid 2n\}.
\]
Referenz: `Staudt.tex`, Mathlib `NumberTheory.Bernoulli`.

**Negativer Befund (bereits dokumentiert):** Bernoulli-Normschalen als Lyapunov-Funktion
für Collatz sind **No-Go** (`collatz_generalangriff_2026.md`). Dieser Zweig ist **Sensorik**,
nicht Hauptangriff.

---

## 4. Vier Schichten innerhalb EABC

$$\boxed{
\text{Äquidistanz}
\;\to\;
\text{Bernoulli}
\;\to\;
\text{Primzahlsignatur}
\;\to\;
\text{Zeta-Spektrum}
}$$

| Schicht | Objekt | EABC-Rolle |
|---------|--------|------------|
| **Äquidistanz** | Triviale Nullstellen $s=-2,-4,\ldots$; äquidistante Stützstellen | Kontinuierliche „Gitter“-Seite |
| **Bernoulli** | $B_{2n}$, $\zeta(1-2n)=-B_{2n}/(2n)$ | Rationale Übersetzer |
| **Primzahlsignatur** | $\mathrm{PrimeSig}(B_{2n})$ | Diskrete Prim-Arithmetik |
| **Zeta-Spektrum** | Triviale + (heuristisch) nichttriviale Nullstellen | Spektrale Projektion |

Alle vier Schichten werden **im selben EABC-Rahmen** gelesen; die mod-$12$-Klassen
sind die native Kodierungsschicht für $\mathrm{PrimeSig}$.

---

## 5. EABC-Kernhypothese (Bernoulli-Übersetzung)

> **Boxed (Conjecture):** Bernoulli-Zahlen **übersetzen** äquidistante Zeta-Geometrie
> in Prim-Arithmetik — innerhalb EABC, nicht als externes Modul.

Präziser: Die Filtration $\{V(B_{2n})\}$ (Abschnitt 7) trägt EABC-Symmetrieinformation
über $\mathrm{PrimeSig}(B_{2n})$, analog zu chiralem Fluss bei Primvierlingen
(`energiedoku_chiraler_fluss.py`). **Label: Conjecture** — empirisch zu prüfen (Test 1),
nicht als Theorem behauptet.

---

## 6. RH-Interpretation (Conjecture / Heuristik)

Nichttriviale Nullstellen $\rho_n=\tfrac12+\mathrm{i}t_n$ werden **heuristisch** als
gekrümmte spektrale Projektionen gelesen. Die Riemann-Hypothese entspricht in dieser
Lesart der Fixierung der Projektion auf den Raum $\mathrm{Re}(s)=\tfrac12$.

**Label: Conjecture / Heuristik** — keine Behauptung eines RH-Beweises.
Kopplung $V(B_{2n})$ vs. Abstände $\Delta t_n$ ist **Test 2 (zukünftig)** und benötigt
Nullstellendaten.

---

## 7. Native EABC mod $12$ und $V(B_{2n})$

Die vier EABC-Familien (Definition, nicht Brücke):
\[
E\equiv 1,\quad A\equiv 5,\quad B\equiv 7,\quad C\equiv 11 \pmod{12}.
\]
Implementierung: `eabc_from_lean.py` (`class_of`, `residue`).

Für $n\ge 1$ definieren wir den **EABC-Zählvektor**
\[
V(B_{2n}) = (e_n,a_n,b_n,c_n),
\]
wobei $e_n$ (bzw. $a_n,b_n,c_n$) die Anzahl der Primzahlen $p\in\mathrm{PrimeSig}(B_{2n})$
in der jeweiligen Restklasse zählt. Primzahlen $2,3$ liegen typischerweise in
$\mathrm{PrimeSig}$, fallen aber außerhalb der vier Klassen; dann gilt
\[
e_n+a_n+b_n+c_n + \#\{p\in\mathrm{PrimeSig}: p\in\{2,3\}\} = |\mathrm{PrimeSig}(B_{2n})|.
\]

---

## 8. Sensor $\Phi_B$ und chirale Observablen

**Sensor:**
\[
\Phi_B(n) := V(B_{2n}).
\]
Implementierung: `collatz_eabc_bernoulli_sensor.py` $\to$ `collatz_eabc_bernoulli_sensor.json`.

**Chirale EABC-Observablen** (analog zu Vierlings-$\sigma$-Tests):
\[
\sigma_B(n) := (e_n+a_n) - (b_n+c_n) \quad\text{(EA/BC-Bilanz)},
\]
\[
\chi_B(n) := (e_n+b_n) - (a_n+c_n) \quad\text{(Diagonalen-Bilanz)},
\]
\[
\iota_{\mathrm{chir}}(n) := \mathrm{sgn}\bigl(\sigma_B(n)\cdot\chi_B(n)\bigr)\in\{-1,0,+1\}.
\]

**Label:** Definition (Observablen) + **Experiment** (Sequenzstatistik in JSON).

---

## 9. Drei Tests

| # | Test | Status | Label |
|---|------|--------|-------|
| **1** | Bernoulli-Signatur + EABC-Asymmetrien ($\sigma_B$, $\chi_B$, $\iota_{\mathrm{chir}}$) | **implementiert** (`collatz_eabc_bernoulli_sensor.py`) | Experiment |
| **2** | Zeta-Kopplung: $V(B_{2n})$ vs. $\Delta t_n$ | **zukünftig** (Nullstellendaten) | — |
| **3** | Krümmung $K_B(n)$ vs. $\pi(x)-\mathrm{Li}(x)$ | **zukünftig** | — |

Test 1 prüft reproduzierbar, ob EABC-Asymmetrien in der PrimeSig-Filtration
systematisch auftreten oder im Zufallsbereich liegen. **Kein Theorem** über deren
Bedeutung für Collatz oder RH.

---

## 10. Parallele Tracks (nicht verschmolzen)

Der EABC-Bernoulli-Zweig läuft **parallel** zu anderen Sensoren — ohne sie zu ersetzen:

| Track | Sensor | Artefakt |
|-------|--------|----------|
| **EABC-Bernoulli** | $\Phi_B$, $V(B_{2n})$ | `collatz_eabc_bernoulli_sensor.json` |
| **Morley** | $F_M$, $S_M$, $K_M$ | `collatz_morley_tm_numerik.py` |
| **$\kappa$ / Grammatik** | $\mathcal{L}_{\mathrm{arith}}$, $F_n$ | `collatz_forbidden_words.py` |

Morley und $\kappa$ bleiben eigenständige Angriffslinien (`collatz_morley_metrik_erweiterung.md`,
`collatz_stufe3_kappa_invarianz.md`). Der Bernoulli-Zweig nutzt **dieselbe** EABC-Sprache,
aber andere Eingangsgrößen (von Staudt statt Collatz-Bahnen).

---

## 11. Tao-Labels und nächste Schritte

| Aussage | Tao-Label |
|---------|-----------|
| EABC mod $12$, $V(B_{2n})$ | Definition |
| PrimeSig aus von Staudt--Clausen | Definition (Theorem der klassischen Zahlentheorie) |
| $\Phi_B$, $\sigma_B$, $\chi_B$, $\iota_{\mathrm{chir}}$ | Definition + Experiment |
| Bernoulli übersetzt Zeta $\to$ Primarithmetik | Conjecture |
| RH als Fixraum-Projektion | Conjecture / Heuristik |
| Bernoulli-Lyapunov für Collatz | Negativtest (No-Go) |

**Nächste saubere Schritte:**
1. Test 1 auf größeres $n$ und mit expliziten Zufalls-Nullmodellen vergleichen.
2. Lean-Schnittstelle: PrimeSig als `Finset` über `Nat.Prime` (analog `BernoulliClock.lean`).
3. Test 2 vorbereiten, sobald $\Delta t_n$-Daten eingebunden sind.

---

*Epistemische Einordnung gemischt je Abschnitt. Keine Ebene ersetzt die nächsthöhere.*
