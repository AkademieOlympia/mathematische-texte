# Hurwitz–EABC-Orbit-Spaltung (glatt-EABC auf $O_p$)

**Status:** Forschungshypothese + Experiment  
**Kanonsiche Erweiterung von:** `collatz_eabc_normabstieg_hypothese.md`, `collatz_eabc_euklidische_hebung.md` (PR #54)  
**Gauß-Vergleich:** `collatz_eabc_gauss_spaltung_hypothese.md`  
**Eisenstein-Vergleich:** `collatz_eabc_eisenstein_spaltung.md`  
**Invarianzprogramm:** `collatz_eabc_invarianzprogramm.md` ($V_4$, Fluktuationsfeld $\delta$)  
**Experiment:** `collatz_eabc_hurwitz_orbit_test.py` → `collatz_eabc_hurwitz_orbit.json`

---

## 1. Einzigkeit in $\mathbb{Z}[i]$ vs. Vielfalt in $\mathbb{H}_{\mathrm H}$

In den **Gaußschen Ganzzahlen** $\mathbb{Z}[i]$ zerlegt eine split-Primzahl $p\equiv 1\pmod 4$
eindeutig (bis auf Einheiten) als $p=N(a+bi)=a^2+b^2$ mit kanonischer Wahl $0<a\le b$.
Die **2D**-Normschale trägt genau **zwei** glatt-EABC-Beine $\Gamma(p)=(\kappa(a'),\kappa(b'))$.

In der **Hurwitz-Maximalordnung** $\mathbb{H}_{\mathrm H}$ gilt für $p>2$:
$$O_p := \{q\in\mathbb{H}_{\mathrm H} : N(q)=p\}$$
ist eine **große Menge** auf der 3-Sphäre $S^3$ (Projektivraum der Quaternionenrichtungen).
Die **24 Hurwitz-Einheiten** wirken durch Links-/Rechtsmultiplikation; es gibt **keine**
kanonische Einzelwahl von $q\in O_p$ analog zu $(a,b)$ in $\mathbb{Z}[i]$.

**Label:** Lagrange-Vier-Quadrate = **Theorem**; Verlust kanonischer Wahl = **ehrliche Limitation**.

---

## 2. $24$ Einheiten und Orbit-Struktur

Die Einheitengruppe $U(\mathbb{H}_{\mathrm H})$ hat Ordnung $24$ ($\pm 1,\pm i,\pm j,\pm k$
und $(\pm1\pm i\pm j\pm k)/2$). Für festes $p$ zerlegt $O_p$ in Bahnen unter
$q\mapsto u\,q\,v$ (Doppelbahnen); die Gesamtgröße $|O_p|$ wächst mit $p$
(siehe `hurwitz_orbit_summary_up_to_199.csv`, `Hurwitz 24.py`).

**Label:** **Theorem** (Hurwitz-Einheiten); Bahnzählung = **Experiment**.

---

## 3. Epistemische Korrektur: $\Gamma(q)$ ist **nicht** orbit-invariant

**Fehlannahme (verworfen):** Ein einzelnes $\Gamma(q)$ für $q\in O_p$ sei ein sinnvolles
Invariantenobjekt unter Einheiten-Bahnen.

**Korrektur:** $\Gamma$ ist eine **Punktfunktion** auf $O_p$. Unter $q\mapsto uqv$
($u,v\in U(\mathbb{H}_{\mathrm H})$) permutiert die Einheitengruppe Koordinatenachsen;
$\Gamma(q)$ ändert sich typischerweise. Das Gauß-Programm funktioniert, weil dort
**ein** kanonischer Repräsentant $(a,b)$ existiert — in $\mathbb{H}$ gibt es diesen nicht.

**Fundamentales Objekt:** die **Orbit-Verteilung**
$$\mu_p(\gamma) := \frac{\#\{q\in O_p : \Gamma(q)=\gamma\}}{|O_p|}, \qquad
\gamma\in\{E,A,B,C,0\}^4.$$
$\mu_p$ ist eine Wahrscheinlichkeit auf EABC-Signaturen, definiert auf der **gesamten**
Normschale $O_p$, nicht auf einem Einzelpunkt.

**Forschungskette (dokumentiert):**
$$p \;\longmapsto\; O_p \;\longmapsto\; \mu_p \;\longmapsto\; I(\mu_p).$$

---

## 3b. Verlust der kanonischen Wahl vs. Gauß-Einzigkeit

| Ring | Dimension | Kanonisches Objekt bei $N(\cdot)=p$ | EABC-Beine | Statistik |
|------|-----------|-------------------------------------|------------|-----------|
| $\mathbb{Z}[i]$ | $2$ | ein Paar $(a,b)$, $0<a\le b$ | $2$ | Punkt-$\Gamma(p)$ |
| $\mathbb{Z}[\omega]$ | $2$ | ein Paar $(a,b)$, $0<a\le b$ | $2$ | Punkt-$\Gamma(p)$ |
| $\mathbb{H}_{\mathrm H}$ | $4$ | **gesamtes** $O_p$ | $4$ pro $q$ | **$\mu_p$** |

Das Gauß-Experiment (`collatz_eabc_gauss_spaltung_hypothese.md`) misst **eindeutige**
$\Gamma(p)$ pro split-Primzahl — dort ist Einzigkeit **Theorem+Definition**.
Im Hurwitz-Fall **scheitert** dieselbe Punkt-Fragestellung; die ehrliche Antwort ist
Orbit-Statistik. Das ist keine Schwäche, sondern die korrekte $4$D-Fragestellung.

---

## 4. Definition ($q=(a,b,c,e)$, glatt-EABC, $\Gamma$)

Schreibe $q=a+bi+cj+ek\in\mathbb{H}_{\mathrm H}$ mit $N(q)=a^2+b^2+c^2+e^2=p$.
Die **vierte** Koordinate heißt $e$ (User-Nomenklatur; entspricht $a_3$ in $\mathbb{H}$).

**Glatter Kern** (wie Gauß/Eisenstein):
$$n = 2^{\alpha}\,3^{\beta}\, n',\qquad \gcd(n',6)=1,\qquad
\kappa(n')\in\{E,A,B,C\}\ \text{via `eabc_from_lean.class_of`}.$$

Auf jeder **nicht-Null**-Koordinate $|x|$: `strip_smooth` → $\kappa(x')$.
Null-Bein: Label `0` (kein EABC-Kern).

$$\Gamma(q) := \bigl(\kappa(a'),\,\kappa(b'),\,\kappa(c'),\,\kappa(e')\bigr)
\in (\{E,A,B,C,0\})^4.$$

**Implementierungswahl (Experiment):** $O_p$ wird als Vereinigung
- **ganzzahliger** Vier-Quadrat-Darstellungen $a^2+b^2+c^2+e^2=p$, und
- **halbganzzahliger** Hurwitz-Darstellungen mit ungeraden Zählern $m_i$,
  $\sum m_i^2 = 4p$, Koordinaten $m_i/2$.

Damit werden alle Norm-$p$-Elemente der Hurwitz-Maximalordnung erfasst
(vgl. `collatz_eabc_hurwitz_orbit_test.py::hurwitz_orbit_elements`).

---

## 5. $\mu_p$ als fundamentales Objekt

Für jede Primzahl $p$ und uniformes Gewicht auf $O_p$:
$$\mu_p(\gamma) := \frac{\#\{q\in O_p : \Gamma(q)=\gamma\}}{|O_p|}.$$

**Implementierung:** $O_p$ = Vereinigung ganzzahliger und halbganzzahliger Vier-Quadrat-
Darstellungen mit $N(q)=p$ (`collatz_eabc_hurwitz_orbit_test.py::hurwitz_orbit_elements`).

Aggregiert über $p\le X$ (gewichtet nach $|O_p|$):
$$\mu_X(\gamma) := \frac{\sum_{p\le X}\#\{q\in O_p:\Gamma(q)=\gamma\}}
{\sum_{p\le X}|O_p|}.$$

**Forschungsfrage:** Zeigt $\mu_p$ (asymptotisch in $p$) Struktur jenseits des
Produkts der Bein-Marginalen? Oder ist $\mu_p$ nahezu isotrop?

**Label:** **Experiment**; keine Beweisbehauptung.

---

## 5b. Invarianten-Kandidaten $I(\mu_p)$

Einzelne $\Gamma(q)$ sind **keine** Invarianten. Kandidaten auf der Verteilung $\mu_p$:

| Sensor | Definition | Interpretation |
|--------|------------|----------------|
| **$H(p)$** | Shannon-Entropie $-\sum_\gamma \mu_p(\gamma)\log\mu_p(\gamma)$ | niedrig = orientiert; hoch = isotrop |
| **$\bar\chi_p$** | $\mathbb{E}_{q\sim\mathrm{Unif}(O_p)}[\chi(q)]$ | Chiralitätsmittel auf dem Orbit |
| **Kanal-Spektrum** | Kovarianz / $\chi^2$ der Kanalpaare $(a',b')$, $(c',e')$ | Kopplung jenseits Marginalen |
| **Unabhängigkeit** | $\chi^2(\mu_p \,\|\, \prod_i \mathrm{marg}_i)$ | Abweichung vom Produktmodell |

**Morley-Parallel (getrennte Spur):** Im Morley-Modul (`collatz_morley_tm_numerik.py`) messen
$G_M=\sum(\theta_i-\pi/3)$ und $W_M=\mathrm{Area}(H_W)/\mathrm{Area}(\Delta)-1/10$
**Konfigurationssensoren** auf Dreiecken — analog als **orbitale Konfigurationssensoren**
auf $\mu_p$, nicht als Beweisbaustein. Morley und Hurwitz-EABC werden **nicht** in einem
Beweis verschmolzen; die Parallelität ist heuristisch (Konfiguration → Sensor).

**Label:** $H(p)$, $\bar\chi_p$ = **Experiment**; Morley-Analogie = **Heuristik**.

---

## 6. Korrelationen und Chiralität

**Kanalpaare** (Vierer-Symmetrie analog $V_4$):
- **AB-Kanal:** $(\kappa(a'),\kappa(b'))$
- **CE-Kanal:** $(\kappa(c'),\kappa(e'))$

**Chiralität** auf einem Punkt $q$:
$$\chi(q) := \#\{\text{Legs}\in\{E,C\}\} - \#\{\text{Legs}\in\{A,B\}\}$$
(nur über nicht-Null-Legs; $\chi\in\{-4,-2,0,2,4\}$).

**Orbit-mittlere Chiralität:** $\bar\chi_p := \frac{1}{|O_p|}\sum_{q\in O_p}\chi(q)$.

**Unabhängigkeitstest:** Vergleich $\mu_X(\Gamma)$ mit dem Produkt der vier Bein-Marginalen
(Pearson-$\chi^2$ gegen Unabhängigkeitshypothese).

---

## 7. Vier Komponenten $\leftrightarrow$ vier EABC-Klassen, $V_4$, Projektion

In $\mathbb{H}$ existiert die **Definition** (nicht Theorem der Zahlentheorie):
$$E\leftrightarrow 1,\quad A\leftrightarrow i,\quad B\leftrightarrow j,\quad C\leftrightarrow k.$$
Die Klein-Vierergruppe $V_4=\{E,A,B,C\}$ wirkt als Paarflip-Struktur auf Vierer-Tupeln
(`collatz_eabc_invarianzprogramm.md`, `PAPER_MASSE_KONSTITUTION.md`).

**Projektion** $\mathbb{H}_{\mathrm H}\to V_4^4$ via $\Gamma$ ist in $\mathbb{C}$ **nicht** sinnvoll
(nur zwei Beine); in $\mathbb{H}$ sind **vier** algebraische Achsen gleichberechtigt.

**Label:** $V_4$-Kodierung = **Definition**; $\Gamma$-Bias = **Conjecture/Experiment**.

---

## 8. Studienobjekt: $O_p\to\mu_p\to I(\mu_p)$, nicht $\Gamma(q)$

Im Gegensatz zu Gauß/Eisenstein wird **nicht** ein kanonisches Faktorpaar gewählt, sondern
die **volle Normschale** $O_p$ und die darauf induzierte Verteilung $\mu_p$.
Primzahlen erscheinen als irreduzible Hurwitz-Elemente mit Primnorm
(`collatz_eabc_euklidische_hebung.md` §3).

**Experimentelle Hypothesen:**
- $H(p)$ korreliert mit $p\bmod 12$ (EABC-Restklassen $1,5,7,11$)?
- $\bar\chi_p\approx 0$ für alle $p$ (kein globaler Chiralitäts-Drift)?
- $\mu_p$ unabhängig von Produkt-Marginalen (wie Gauß bei großem $X$)?

**Ehrliche Einordnung:** Das Gauß-$\Gamma(p)$-Programm liefert Punkt-Einzigkeit; das
Hurwitz-Programm testet, ob **Orbit-Statistik** noch EABC-Struktur trägt. Negative
Ergebnisse (hohes $H(p)$, $\mu_p\approx$ Produkt) wären **informativ**, nicht das
Scheitern des Gesamtprogramms.

---

## 9. Vergleich Gauß / Eisenstein / Hurwitz

| Aspekt | Gauß $\mathbb{Z}[i]$ | Eisenstein $\mathbb{Z}[\omega]$ | Hurwitz $\mathbb{H}_{\mathrm H}$ |
|--------|----------------------|--------------------------------|----------------------------------|
| Normbeine | $2$ | $2$ | $4$ |
| Kanonische Wahl | ja | ja | **nein** ($|O_p|\gg 1$) |
| $\Gamma$-Raum | $\{E,A,B,C\}^2$ | $\{E,A,B,C\}^2$ | $\{E,A,B,C,0\}^4$ |
| Split-Filter | $p\equiv 1\pmod 4$ | $p\equiv 1\pmod 3$ | alle $p$ (Lagrange) |
| Epistemik | Punkt-$\Gamma$ | Punkt-$\Gamma$ | **Orbit-$\Gamma$** |

---

## 10. Epistemische Einordnung

| Aussage | Label |
|---------|-------|
| Jedes $p$ ist Summe von vier Quadraten | **Theorem** (Lagrange) |
| $|O_p|$ und Einheiten-Bahnen | **Theorem** + **Experiment** |
| $\mu_p$ nahezu isotrop (hohes $H(p)$) | **offene Conjecture** |
| $H(p)$-Bias nach $p\bmod 12$ | **Experiment** |
| Orbit-$\chi$-Bias ($\bar\chi_p$) | **Experiment** |
| $\Gamma(q)$ als Invariante | **verworfen** |
| Verlust kanonischer Repräsentant | **ehrliche Limitation** |

**Falsifikation:** $\bar\chi_p\approx 0$, $H(p)$ hoch und ohne Restklassen-Struktur,
$\mu_p\approx$ Produkt der Marginalen, keine AB–CE-Kopplung → Hurwitz-Orbit trägt
**keine** zusätzliche EABC-Orientierung jenseits der mod-$12$-Marginalen.

---

## Querverweise

- `collatz_eabc_normabstieg_hypothese.md` — Normabstieg, Gauß-Brücke §8–§9
- `collatz_eabc_gauss_spaltung_hypothese.md` / `_test.py` — 2D Punkt-$\Gamma$ Referenz
- `collatz_eabc_eisenstein_spaltung.md` / `_test.py` — hexagonale Referenz
- `collatz_eabc_euklidische_hebung.md` — Hurwitz-Euklidizität, PR #54
- `collatz_hurwitz_polytop_eabc.tex` — Polytop-Geometrie, $p\bmod 12$
- `Hurwitz 24.py` — Einheiten und Bahnen
- `collatz_morley_tm_numerik.py` — $G_M$, $W_M$ als Konfigurationssensoren (Morley-Parallel)
- `collatz_morley_gm_beweisversuch.md` — Morley-Sensorik, getrennte Spur

---

*Epistemische Einordnung: $\Gamma(q)$ orbit-invariant = **falsch**; $\mu_p$, $H(p)$ = Experiment;
Gauß-Punkt-Einzigkeit ≠ Hurwitz-Orbit-Statistik.*
