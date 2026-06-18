# Hurwitz–EABC-Orbit-Spaltung (glatt-EABC auf $O_p$)

**Status:** Forschungshypothese + Experiment  
**Kanonsiche Erweiterung von:** `collatz_eabc_normabstieg_hypothese.md`, `collatz_eabc_euklidische_hebung.md` (PR #54)  
**Gauß-Vergleich:** `collatz_eabc_gauss_spaltung_hypothese.md`  
**Eisenstein-Vergleich:** `collatz_eabc_eisenstein_spaltung.md`  
**Invarianzprogramm:** `collatz_eabc_invarianzprogramm.md` ($V_4$, Fluktuationsfeld $\delta$)  
**Erweiterte Maßhypothese ($\Sigma_p$, $\mu_p$, $H_p$):** `collatz_eabc_quaternion_mass_hypothese.md`  
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

## 3. Verlust der kanonischen Wahl (epistemisch)

| Ring | Dimension | Kanonisches Objekt bei $N(\cdot)=p$ | EABC-Beine |
|------|-----------|-------------------------------------|------------|
| $\mathbb{Z}[i]$ | $2$ | ein Paar $(a,b)$, $0<a\le b$ | $2$ |
| $\mathbb{Z}[\omega]$ | $2$ | ein Paar $(a,b)$, $0<a\le b$ | $2$ |
| $\mathbb{H}_{\mathrm H}$ | $4$ | **gesamtes** $O_p$, nicht ein Punkt | $4$ pro $q$ |

Die Hurwitz-Variante testet daher **Orbit-Statistik**, nicht Punkt-Statistik.
Das ist keine Schwäche des Programms, sondern die korrekte Fragestellung in $4$D.

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
(vgl. `collatz_eabc_hurwitz_orbit_test.py::hurwitz_shell_elements`).

---

## 5. Forschungsfrage (Verteilung von $\Gamma$ auf $O_p$)

Für jede Primzahl $p$ und $q\in O_p$ uniform gewichtet:
$$\mu_p(\gamma) := \frac{\#\{q\in O_p : \Gamma(q)=\gamma\}}{|O_p|}.$$

**Forschungsfrage:** Ist $\mu_p$ (asymptotisch in $p$) uniform auf beobachteten Signaturen?
Aggregiert über $p\le X$:
$$\mu_X(\gamma) := \frac{\sum_{p\le X}\#\{q\in O_p:\Gamma(q)=\gamma\}}
{\sum_{p\le X}|O_p|}.$$

**Label:** **Experiment**; keine Beweisbehauptung.

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

## 8. Studienobjekt: Orbit $O_p$, nicht Einzelpunkt

Im Gegensatz zu Gauß/Eisenstein wird **nicht** ein kanonisches Faktorpaar gewählt, sondern
die **volle Normschale** $O_p$. Primzahlen erscheinen als irreduzible Hurwitz-Elemente
mit Primnorm (`collatz_eabc_euklidische_hebung.md` §3).

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
| $\Gamma$-Verteilung auf $O_p$ uniform | **offene Conjecture** |
| Orbit-$\chi$-Bias | **Experiment** |
| Verlust kanonischer Repräsentant | **ehrliche Limitation** |

**Falsifikation:** $\bar\chi_p\approx 0$ für alle $p$, $\mu_X\approx$ Produkt der Marginalen,
keine AB–CE-Kopplung → Hurwitz-Orbit trägt **keine** zusätzliche EABC-Orientierung jenseits
der mod-$12$-Marginalen.

---

## Querverweise

- `collatz_eabc_normabstieg_hypothese.md` — Normabstieg, Gauß-Brücke §8–§9
- `collatz_eabc_gauss_spaltung_hypothese.md` / `_test.py` — 2D Referenz
- `collatz_eabc_eisenstein_spaltung.md` / `_test.py` — hexagonale Referenz
- `collatz_eabc_euklidische_hebung.md` — Hurwitz-Euklidizität, PR #54
- `collatz_eabc_quaternion_mass_hypothese.md` — erweiterte $\Sigma_p$, $H_p$, $K_p$
- `collatz_hurwitz_polytop_eabc.tex` — Polytop-Geometrie
- `Hurwitz 24.py` — Einheiten und Bahnen

---

*Epistemische Einordnung: Hurwitz-Orbit = Theorem+Experiment; $\Gamma$-Bias auf $O_p$ = offene Conjecture.*
