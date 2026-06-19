# Zeta-Exponential-Gedankenexperiment

**Status:** Gedankenexperiment / Analogie (Schicht **C**) — **kein** Theorem über EABC, Primzahlen oder RH  
**Branch:** `collatz/eabc-h12e-zirkulationshypothese` (PR #79)  
**Tao-Labels:** Definition | Theorem (klassisch) | Analogie | Ikone

$$\boxed{\;\textbf{Schicht C:}\;\text{Dieses Dokument ist Interpretationsschicht — keine mathematische Konsequenz für }D_E,\,\Phi_E\text{ oder EABC-Hypothesen.}\;}$$

**Querverweise:**
- `collatz_eabc_zirkulationshypothese.md` — §4.9 (Kurzverweis); Level-2-Fluktuationsgeometrie §4.8 (Nullmodell-Hierarchie §4.8.2: $\zeta_F$ als **reguläres** Vergleichsensemble?); links-rechts-Paarung §2
- `collatz_eabc_epistemik_schichten.md` — Schichten A/B/C/R; Abgrenzung Ikone vs. Struktur; Hurwitz-/Basel-Ikone §2
- `collatz_eabc_diskrete_geometrie.md` — kanonische EABC-Definitionen ($G_E$, $\gamma^\pm$)
- `collatz_offene_punkte.md` — Stirling–Bernoulli–Zeta-Anbindung; Bernoulli-Uhr (`BernoulliClock.lean`)
- `collatz_bernoulli_schalen.pdf` — Bernoulli-Normschalen und $\zeta$-Anbindung (PDF)
- `collatz_oktonionen_beweis.pdf` — oktanionischer Beweis-/Programmkontext (PDF)
- `divisionsalgebren.tex` — Satz von Hurwitz, $\mathbb{R},\mathbb{C},\mathbb{H},\mathbb{O}$
- `collatz_mathlib_eabc_kandidaten.md` — $\zeta(2n)$ via Bernoulli (`riemannZeta_two_mul_nat`); Hurwitz-Namensfalle

---

## 1. Epistemischer Rahmen

Dieses Dokument sammelt ein **Gedankenexperiment**, das die klassische $\zeta$-Struktur mit **Exponentiallstreifen** $(e^x\pm e^{-x})$ in Beziehung setzt und daraus eine **heuristische Lesart** für das EABC-Programm ableitet.

| Ebene | Inhalt | Label |
|-------|--------|-------|
| **Klassisch** | Die drei Boxen in §2 sind **Standard** der analytischen Zahlentheorie | **Theorem** / **Definition** |
| **Heuristik** | Brücke zu logarithmischen Koordinaten, links-rechts-Paarung, Level-2-Geometrie | **Analogie** (Schicht **C**) |
| **Explizit ausgeschlossen** | RH-Beweis, Primzahltheorem, EABC-Holonomie, $\Phi_E\neq 0$ | **nicht behauptet** |

$$\boxed{\;\text{Zeta ist die Exponentialfunktion auf der logarithmischen Primzahl-Geometrie.}\;}$$

**Lesart des Leitsatzes:** Auf der **logarithmischen Achse** $n\mapsto\log n$ (bzw. $p\mapsto\log p$) erscheint die Dirichlet-Reihe $\zeta(s)=\sum n^{-s}$ als **gewichtete Exponentialsumme**; die Euler-Produktstruktur faktorisiert lokal in **Differenz-×-Summen-Paare** $(1-p^{-s})(1+p^{-s})$, und die Funktionalgleichung der abgeschlossenen $\xi$-Funktion ist **Spiegelung** $z\leftrightarrow -z$ um die kritische Mitte. Das ist **klassische Mathematik** — die EABC-Brücke in §4 überträgt dieselbe **Bildsprache** auf Transport, Chiralität und Level-2-Kovarianz, **ohne** Beweisanspruch.

**Label:** gesamter Abschnitt = **Analogie** / **Ikone** (Schicht **C**), außer explizit als klassisch markierte Formeln.

---

## 2. Drei klassische Boxen (analytische Zahlentheorie)

Die folgenden drei Relationen sind **exakt** (im jeweiligen Konvergenzbereich bzw. nach analytischer Fortsetzung). Sie werden hier **nicht** bewiesen, sondern als **Referenzgerüst** für das Gedankenexperiment fixiert.

### Box 1 — Zeta als Exponentialsumme auf der Log-Skala

$$\boxed{\;\zeta(s)=\sum_{n\ge 1} n^{-s}=\sum_{n\ge 1} e^{-s\log n}\qquad(\Re(s)>1).\;}$$

**Bedeutung:** Die Riemannsche Zetafunktion ist auf der **multiplikativen** Skala $\mathbb{N}_{>0}$ definiert; in **logarithmischen Koordinaten** $x_n:=\log n$ wird jeder Term eine **reine Exponentialfunktion** $e^{-sx_n}$. Die Summe aggregiert alle natürlichen „Log-Punkte“ mit Gewicht $e^{-sx}$.

**Label:** **Definition** / **Theorem** (klassisch).

### Box 2 — Lokaler Euler-Faktor und Streifen-Produkt

$$\boxed{\;(1-p^{-s})(1+p^{-s})=1-p^{-2s}\qquad(p\text{ prim}).\;}$$

**Algebra:** Differenzfaktor $(1-p^{-s})$ und Summenfaktor $(1+p^{-s})$ multiplizieren zu **Quadrat-Differenz** $1-p^{-2s}$.

**Streifen-Analogie (nach Normierung):** Setze $x:=\tfrac{s}{2}\log p$. Dann $p^{-s}=e^{-2x}$ und

$$(1-p^{-s})(1+p^{-s})=1-e^{-2s\log p}=1-p^{-2s}.$$

Auf der **reellen Streifen-Seite** (unabhängig von $p,s$):

$$\boxed{\;(e^x-e^{-x})(e^x+e^{-x})=e^{2x}-e^{-2x}=2\sinh(2x).\;}$$

**Strukturelle Parallele (heuristisch, nicht identisch):**

| Zeta / Euler (lokal) | Exponentialstreifen |
|----------------------|---------------------|
| $(1-p^{-s})$ — „Differenz“ gegen $1$ | $(e^x-e^{-x})\propto\sinh(x)$ — antisymmetrisch |
| $(1+p^{-s})$ — „Summe“ mit $1$ | $(e^x+e^{-x})\propto\cosh(x)$ — symmetrisch |
| Produkt $=1-p^{-2s}$ | Produkt $=2\sinh(2x)$ (nach Skalierung) |

Die Zuordnung $x\leftrightarrow\tfrac{s}{2}\log p$ ist **Bildsprache**: sie verknüpft **Prim-Euler-Faktoren** mit **hyperbolischen Streifen** — **kein** Theorem über Primzahlen oder EABC.

**Label:** Euler-Identität = **Theorem** (klassisch); Streifen-Tabelle = **Analogie** (Schicht **C**).

### Box 3 — Funktionalgleichung als $e^z\leftrightarrow e^{-z}$-Symmetrie

Sei $\xi(s)$ die **abgeschlossene** Zetafunktion (Riemann-Siegel-Completion),

$$\xi(s)=\tfrac12 s(s-1)\,\pi^{-s/2}\,\Gamma\!\left(\tfrac{s}{2}\right)\zeta(s),$$

die für alle $s\in\mathbb{C}$ meromorph ist und die Funktionalgleichung

$$\boxed{\;\xi(s)=\xi(1-s)\;}$$

erfüllt. Mit $s=\tfrac12+z$ folgt

$$\boxed{\;\xi\!\left(\tfrac12+z\right)=\xi\!\left(\tfrac12-z\right).\;}$$

**Lesart:** Um die **kritische Mitte** $\Re(s)=\tfrac12$ ist die abgeschlossene $\xi$-Funktion **gerade** in $z$ — formal dieselbe **Involution** $z\mapsto -z$, die in $e^z\leftrightarrow e^{-z}$-Paaren auftritt. Das ist **keine** Behauptung über die Riemann-Hypothese; RH wäre die **stärkere** Aussage, dass alle nichttrivialen Nullstellen auf $\Re(s)=\tfrac12$ liegen.

**Label:** Funktionalgleichung = **Theorem** (klassisch); $e^z$-Symmetrie-Lesart = **Analogie** (Schicht **C**).

---

## 3. Exponentialstreifen-Geometrie (Bildsprache)

### 3.1 Logarithmische Primzahl-Geometrie

$$\text{Primzahl }p \;\longmapsto\; \text{Log-Punkt } x_p:=\log p \;\longmapsto\; \text{lokaler Faktor } (1-p^{-s})(1+p^{-s}).$$

Die **globale** Zeta-Funktion entsteht durch **Summation** über alle $n$ (Box 1) bzw. **Produktion** über alle $p$ (Euler-Produkt). Geometrisch: eine **eindimensionale Log-Achse** mit diskreten Prim-Marken, auf der **Exponentiallgewichte** $e^{-s\log n}$ addiert werden.

### 3.2 Differenz–Summen-Paarung

Die Faktoren $(1-p^{-s})$ und $(1+p^{-s})$ bilden ein **links-rechts-artiges Paar**:

- der **Differenz**-Anteil misst Abweichung von der Einheit (antisymmetrische Komponente nach Log-Transformation);
- der **Summen**-Anteil misst kumulative Kopplung mit der Einheit (symmetrische Komponente).

Nach Multiplikation kollabiert das Paar zur **doppelten Frequenz** $p^{-2s}$ — analog zum Übergang $\sinh(x)\cosh(x)\propto\sinh(2x)$.

### 3.3 Spiegelung an der kritischen Mitte

Die Identität $\xi(\tfrac12+z)=\xi(\tfrac12-z)$ ist eine **zentrale Symmetrie** der abgeschlossenen Zeta-Funktion. Im Streifen-Bild: **Austausch** $e^z\leftrightarrow e^{-z}$ bei Fixierung der **mittleren Phase** ($\Re(s)=\tfrac12$).

**Label:** §3 = **Analogie** / **Ikone** — didaktische Verdichtung der Boxen 1–3.

---

## 4. Heuristische Brücke zu EABC (Interpretationsschicht)

$$\boxed{\;\text{Die folgende Brücke ist \textbf{methodisch}, nicht \textbf{beweisend}. Sie verknüpft Bildsprache — keine Theoreme.}\;}$$

### 4.1 Logarithmische Koordinaten

| Zeta-Bild | EABC-Programm | Label |
|-----------|---------------|-------|
| $n\mapsto\log n$, Gewichte $e^{-s\log n}$ | Primstrom als diskrete Folge; Skalierung über $Q(X)$, $\alpha_{\mathrm{loc}}$ | **Analogie** |
| Euler-Produkt über $p$ | HL-Kanäle, mod-$12$- / mod-$60060$-Schalen (§4.4–4.5 in `collatz_eabc_zirkulationshypothese.md`) | **Analogie** |
| Lokaler Faktor pro Prim | Transportkante $\tau_n=(\kappa(p_n),\kappa(p_{n+1}))$ auf $G_E$ | **Analogie** |

EABC arbeitet **nicht** mit $\zeta(s)$ als Observable; die Brücke betrifft nur die **Koordinatenwahl** (logarithmische Skalen, Prim-Marken).

### 4.2 Links–rechts-Paarung

| Zeta / Streifen | EABC | Label |
|-----------------|------|-------|
| $(1-p^{-s})$ vs. $(1+p^{-s})$ | $\gamma^+$ (ABCEA) vs. $\gamma^-$ (CEABC) | **Analogie** |
| $e^x-e^{-x}$ vs. $e^x+e^{-x}$ | antisymmetrische vs. symmetrische Streifenkomponente | **Analogie** |
| $D_E=N_+-N_-$, $Q=N_++N_-$ | Differenz–Summen-Paarung auf demselben $C_4$-Zyklus | **Definition** (EABC) + **Analogie** (Zeta-Bild) |

**Wichtig:** $D_E$ und $W_E=D_E/Q$ sind **definierte EABC-Observablen** (`collatz_eabc_zirkulationshypothese.md` §3). Die Zuordnung zu Euler-Faktoren **begründet keinen** analytischen Zusammenhang und **beweist nichts** über das Wachstum von $|D_E|$.

### 4.3 Level-2-Fluktuationsgeometrie

Level-2 verschiebt den Fokus von Mittelwerten ($D_E$, $W_E$) zur **Kovarianzgeometrie** $\Sigma_A$ im Raum $\Lambda^2(\mathbb{R}^4)$ (`collatz_eabc_zirkulationshypothese.md` §4.8):

| Zeta / Streifen (Bild) | Level-2 (EABC) | Label |
|------------------------|----------------|-------|
| Paarprodukt $(e^x-e^{-x})(e^x+e^{-x})$ erzeugt **zweite Harmonische** $\sinh(2x)$ | antisymmetrische Pfadsignatur $a\in\Lambda^2(\mathbb{R}^4)$; $\Sigma_A=\mathbb{E}[(a-\mu_A)(a-\mu_A)^T]$ | **Analogie** |
| Funktionalgleichung: Spiegelung $z\mapsto -z$ | $\Delta_F$: Abweichung $\Sigma_A^{\mathrm{prime}}$ vs. Nullmodell | **Experiment** (Level-2) |
| Information in **Quadrat** ($p^{-2s}$), nicht in linearem Term | $\mu_A\approx 0$, Struktur in $\Sigma_A$ (§4.8) | **Analogie** (methodisch) |

**Epistemische Grenze:** Ein Level-2-Befund ($\Delta_F>0$) ist **numerisches Experiment** — er wird durch dieses Zeta-Gedankenexperiment **weder erklärt noch bewiesen**.

### 4.4 Was ausdrücklich nicht folgt

| Nicht behauptet | Grund |
|-----------------|-------|
| RH oder Aussagen über $\zeta$-Nullstellen | Funktionalgleichung $\neq$ RH |
| $D_E$ als $\zeta$- oder $L$-Fehlerterm | heuristische Ähnlichkeit zu Prime Races, kein Identitätsbeweis |
| $\Phi_E\neq 0$ aus Zeta-Symmetrie | Schicht **C** exportiert nicht nach Schicht **B** |
| „Zeta beweist EABC-Struktur“ | **Gedankenexperiment**, kein Theorem |

---

## 5. Kurzform

$$\boxed{\;\zeta(s)=\sum_{n\ge 1} e^{-s\log n}.\;}$$

$$\boxed{\;(1-p^{-s})(1+p^{-s})=1-p^{-2s}\;\Leftrightarrow\;(e^x-e^{-x})(e^x+e^{-x})=2\sinh(2x)\;\text{(Streifen-Bild, normiert).}\;}$$

$$\boxed{\;\xi\!\left(\tfrac12+z\right)=\xi\!\left(\tfrac12-z\right)\;\Leftrightarrow\; e^z\leftrightarrow e^{-z}\text{-Symmetrie an der kritischen Mitte.}\;}$$

$$\boxed{\;\textbf{Zeta ist die Exponentialfunktion auf der logarithmischen Primzahl-Geometrie.}\;}$$

**EABC-Brücke (ein Satz):** Log-Koordinaten, Differenz–Summen-Paarung ($\gamma^\pm$, $D_E/Q$) und Level-2-Kovarianz ($\Sigma_A$) lesen sich **parallel** zur klassischen $\zeta$-Streifen-Geometrie — als **Interpretationsschicht** (Schicht **C**), nicht als Beweis.

| Aussage | Schicht | Tao-Label |
|---------|---------|-----------|
| Boxen 1–3 (klassisch) | — | **Theorem** / **Definition** |
| Streifen- und Spiegelungsbild | **C** | **Analogie** |
| EABC-Brücke §4 | **C** | **Analogie** / **Ikone** |
| $D_E$, $\Sigma_A$, $\Delta_F$ | **B** / Experiment | **Definition** / **Experiment** |
| §6–§8 (Ganzzahl/Halbzahl, Bernoulli) | **C** (+ klassisch in Boxen) | **Theorem** / **Analogie** |
| §9 (Fibonacci-Zeta) | **C** (+ klassisch in Boxen) | **Definition** / **Hypothese** / **Analogie** |

---

## 6. Ganzzahlige und halbzahlige Exponenten

Die Streifen-Normierung in Box 2 setzt bereits $x=\tfrac{s}{2}\log p$ — der **Faktor $\tfrac12$** in der Exponentenachse ist kein Zufall der Notation, sondern markiert eine **halbzahlige Gitterstruktur** auf der Log-Skala: $p^{-s}=e^{-2x}$ mit $x=\tfrac{s}{2}\log p$.

### 6.1 Exponentengitter $\mathbb{Z}$ und $\mathbb{Z}+\tfrac12$

Wir **beschränken** das Gedankenexperiment (methodisch, nicht als Konvergenzbehauptung für die volle Dirichlet-Reihe) auf Exponenten

$$s\in\mathcal{S}:=\mathbb{Z}\;\cup\;\Bigl(\tfrac12+\mathbb{Z}\Bigr)=\Bigl\{\ldots,-1,0,1,\ldots\Bigr\}\;\cup\;\Bigl\{\ldots,-\tfrac12,\tfrac12,\tfrac32,\ldots\Bigr\}.$$

| Klasse | Beispiel | Rolle im Streifen-Bild |
|--------|---------|------------------------|
| $s\in\mathbb{Z}$ | $s=2$: $\zeta(2)=\pi^2/6$ | ganzzahlige „Frequenz“ auf $\log n$; Euler-Faktor $1-p^{-2s}$ bei $s\in\mathbb{Z}$ |
| $s=\tfrac12+it$ | kritische Linie | reeller Anteil **halbzahlig**; $\xi(\tfrac12+z)=\xi(\tfrac12-z)$ (Box 3) |
| $s=\tfrac{n}{2}$, $n\in\mathbb{Z}$ | $s=1$: Rand des Konvergenzstreifens | $p^{-s}=p^{-n/2}$; Verdopplung $s\mapsto 2s$ in Box 2 |

**Lesart:** Die **Verdopplung** $s\mapsto 2s$ im Streifen-Produkt $(1-p^{-s})(1+p^{-s})=1-p^{-2s}$ ist auf $\mathcal{S}$ **geschlossen**: aus $s\in\mathcal{S}$ folgt $2s\in\mathbb{Z}$. Halbzahlige Exponenten liefern **ganzzahlige** Quadrat-Frequenzen in der Euler-Faktorisierung — parallel zur Übergangsregel $\sinh(x)\cosh(x)\propto\sinh(2x)$.

**Label:** algebraische Identitäten = **Theorem** (klassisch); Gitter-Einschränkung = **Definition** (Gedankenexperiment).

### 6.2 Klassische Spezialwerte (ganze gerade/ungerade Exponenten)

$$\boxed{\;\zeta(-2n)=0\qquad(n\in\mathbb{N}_{>0}).\;}$$

$$\boxed{\;\zeta(2n)=(-1)^{n+1}\,\frac{B_{2n}\,(2\pi)^{2n}}{2\,(2n)!}\qquad(n\in\mathbb{N}_{>0}),\;}$$

wobei $B_{2n}$ die **Bernoulli-Zahlen** ($B_0=1$, $B_2=1/6$, $B_4=-1/30$, …) bezeichnen.

Für **negative ungerade** ganze Exponente (nach analytischer Fortsetzung):

$$\boxed{\;\zeta(1-2n)=-\frac{B_{2n}}{2n}\qquad(n\in\mathbb{N}_{>0}).\;}$$

Die **Funktionalgleichung** $\xi(s)=\xi(1-s)$ verknüpft Werte bei $s$ und $1-s$: ganzzahlige Spezialwerte an der einen Seite bestimmen (über $\Gamma$-Faktoren und Polstellen) halbzahlige bzw. komplementäre Information — z. B. liegt die kritische Mitte $\Re(s)=\tfrac12$ genau **zwischen** $0$ und $1$.

**Label:** alle Boxen in §6.2 = **Theorem** (klassisch, Standardreferenz).

### 6.3 Halbzahlige Achse und kritische Linie

Mit $s=\tfrac12+it$ ist $\Re(s)$ **fix halbzahlig**; die Funktionalgleichung wird zur **Geradigkeit** $\xi(\tfrac12+z)=\xi(\tfrac12-z)$ (Box 3). Auf der Streifen-Seite: $p^{-s}=p^{-1/2}\cdot p^{-it}$ faktorisiert in einen **reellen Halbzahl-Anteil** $p^{-1/2}$ (Wurzel-Gewichtung auf der Log-Achse) und eine **reine Phase** $p^{-it}=e^{-it\log p}$.

**Wichtig:** Die Einschränkung auf $\mathcal{S}$ **ersetzt nicht** die analytische Theorie auf $\mathbb{C}$; sie ist eine **Lesart-Schicht** für das Gedankenexperiment — welche Exponenten im Streifen-Bild „resonieren“ und welche durch Bernoulli/Funktionalgleichung algebraisch geschlossen sind.

**Label:** §6.3 = **Analogie** (Schicht **C**), außer explizit geboxte Identitäten.

---

## 7. Oktanionische Einordnung (Schicht C)

$$\boxed{\;\text{Oktanionische Brücke = \textbf{Analogie} (Schicht C) — keine Behauptung über EABC-Beweise oder Collatz.}\;}$$

### 7.1 Hurwitz-Kette und Dimensionsleiter

Nach dem **Satz von Hurwitz** (normierte Divisionsalgebren) gibt es genau vier Stufen

$$\mathbb{R}\;(1)\;\subset\;\mathbb{C}\;(2)\;\subset\;\mathbb{H}\;(4)\;\subset\;\mathbb{O}\;(8),$$

jeweils mit Dimensionsverdopplung $1\to2\to4\to8$. Siehe `divisionsalgebren.tex`; oktanionischer Programmkontext: `collatz_oktonionen_beweis.pdf`.

| Algebra | Dimension | EABC-Programm (Lesart) | Label |
|---------|-----------|--------------------------|-------|
| $\mathbb{R}$ | $1$ | reelle Skalare, $Q(X)$, Fenstergröße | **Definition** (B) |
| $\mathbb{C}$ | $2$ | Phase, $\Phi_{\mathrm{pref}}$, komplexe Projektion | **Analogie** (C) |
| $\mathbb{H}$ | $4$ | **Quaternionischer Schatten:** $C_4\cong S^1$, $V=\{E,A,B,C\}$ mod $12$, $\gamma^\pm$ | **Definition** (B) + **Analogie** (C) |
| $\mathbb{O}$ | $8$ | **Schicht-C-Programm:** maximale nichtassoziative „Phase“; $G_2=\mathrm{Aut}(\mathbb{O})$ | **Analogie** (C) |

**EABC auf $C_4$** ist kanonisch **quaternionisch** (vier Restklassen, vierdimensional als $H_1(C_4,\mathbb{Z})$-Gerüst) — nicht oktanionisch im Beweissinne (`collatz_eabc_diskrete_geometrie.md`, `collatz_eabc_epistemik_schichten.md` §2). Die Oktanionen treten als **Nebenprogramm / Ikone** auf: mögliche **Hebung** von $4$ auf $8$ Freiheitsgrade in der Interpretationsschicht.

### 7.2 Halbzahlige Exponenten und Verdopplungsmuster

Die Hurwitz-Dimensionen verdoppeln sich: $2^k$ für $k=0,1,2,3$. **Vorsichtige** Parallele zum Exponentengitter §6:

| Hurwitz | Exponentengitter (Bild) | Gemeinsamkeit |
|---------|-------------------------|---------------|
| $1,2,4,8$ | Schrittweite $\tfrac12$ auf $s$-Achse (via $x=\tfrac{s}{2}\log p$) | **Verdopplung** als strukturelles Motiv |
| $\mathbb{H}\to\mathbb{O}$: Assoziativität bricht | $s\mapsto 2s$ in Box 2 | „Quadrat“-Information ($p^{-2s}$) statt linearer Term |

**Spinor / doppelte Überlagerung (nur Bildsprache):** In der Physik tritt $\mathrm{Spin}(n)$ als **zweifache Überlagerung** von $\mathrm{SO}(n)$ auf — eine **$2\pi$- vs. $4\pi$-Periodizität**. Die halbzahlige kritische Linie $\Re(s)=\tfrac12$ und die Normierung $x=\tfrac{s}{2}\log p$ erinnern **formal** an „halbe Schritte“ auf einer überlagerten Achse. Das ist **kein** Theorem über Spinoren, EABC oder $\zeta$-Nullstellen; es motiviert nur, warum $\mathcal{S}=\mathbb{Z}\cup(\tfrac12+\mathbb{Z})$ im Streifen-Gedankenexperiment **natürlicher** wirkt als beliebige $s\in\mathbb{C}$.

**Label:** gesamter §7 = **Analogie** / **Ikone** (Schicht **C**).

### 7.3 Mod $12$, $C_4$ und oktanionisches „Schicht-C-Deck“

| Objekt | Schicht | Rolle |
|--------|---------|-------|
| $G_E$, $C_4$, $\gamma^\pm$, $D_E$ | **B** | beweisbare/kanonische diskrete Struktur |
| Hurwitz-24-Einheiten, $24I_3$ | **A** / Ikone | `collatz_hurwitz_polytop_eabc.tex`, `collatz_eabc_epistemik_schichten.md` |
| $\mathbb{O}$, $G_2$, Symmetriebruch $\mathbb{O}\to\mathbb{H}$ | **C** | `divisionsalgebren.tex`, `AntiInflation.tex` (spekulativ) |

**Explizit ausgeschlossen:** EABC-Holonomie oder $D_E$-Wachstum **folgt nicht** aus oktanionischer Algebra; PR-Kontext „spektral-oktonion“ (falls vorhanden) ist **nicht** Teil des kanonischen B-Kerns.

---

## 8. Bernoulli-Subtraktion und regulierte analytische Größen

Dieser Abschnitt präzisiert, was „**Bernoulli subtrahieren**“ im Zeta-Streifen-Gedankenexperiment **konkret** bedeuten kann — und trennt **klassische** Regularisierung von **EABC-Analogie**.

### 8.1 Bernoulli in $\zeta(2n)$ und $\zeta(1-2n)$

Die Identitäten in §6.2 zeigen: **gerade positive** Spezialwerte $\zeta(2n)$ sind **reine Bernoulli-Zahlen** (mit $\pi$-Gewicht); **negative ungerade** Werte $\zeta(1-2n)$ sind **rationale Vielfache** $B_{2n}/(2n)$. Die Funktionalgleichung transportiert diese Information zwischen $s$ und $1-s$.

$$\boxed{\;\text{Spezialwerte an ganzzahligen Stellen }s\in\mathbb{Z}\text{ sind algebraisch durch }B_{2n}\text{ geschlossen (bis auf bekannte Transzendenz, z.\,B. }\pi\text{).}\;}$$

**Label:** **Theorem** (klassisch). Verknüpfung zu `collatz_bernoulli_schalen.pdf`, Bernoulli-Uhr in `collatz_offene_punkte.md` ($B_{2m}=-2m\,\zeta(1-2m)$), Lean `BernoulliClock.lean`.

### 8.2 Euler–Maclaurin: Subtraktion divergenter Anteile

Für glatte $f$ mit $\sum_{k\ge1}f(k)$ liefert **Euler–Maclaurin** eine asymptotische Entwicklung

$$\sum_{k=1}^{N} f(k)=\int_1^{N} f(x)\,dx+\frac{f(1)+f(N)}{2}+\sum_{j=1}^{m}\frac{B_{2j}}{(2j)!}\bigl(f^{(2j-1)}(N)-f^{(2j-1)}(1)\bigr)+R_m.$$

Die **Bernoulli-Gewichte** $B_{2j}/(2j)!$ sind die universellen Koeffizienten, mit denen **Rand- und Derivativterme** subtrahiert werden, um die Summe an die Integral-Hauptform anzugleichen. Für $f(k)=k^{-s}$ ($\Re(s)>1$) entsteht so die analytische Fortsetzung von $\zeta(s)$; die **Pole** und **triviale Nullstellen** bei negativen geraden $s$ sind das Ergebnis dieser **Subtraktions-/Kompensationsstruktur**.

**Lesart „Bernoulli subtrahieren“:** Man entfernt (formal oder asymptotisch) den **Polynom-/Potenz-Exponentialanteil** und behält den **endlichen Rest** bzw. die **regulierte Summe**.

**Label:** Euler–Maclaurin-Formel = **Theorem** (klassisch); Stirling–Bernoulli–Zeta-Kette in `collatz_offene_punkte.md` = **spekulativ** (separater Kontext).

### 8.3 Was „analytisches Resultat“ hier bedeutet

| Begriff | Bedeutung | Beispiel |
|---------|-----------|----------|
| **Spezialwert** | exakter Wert an $s\in\mathbb{Z}$ | $\zeta(2)=\pi^2/6$, $\zeta(-1)=-\tfrac1{12}$ |
| **Endlicher Teil** (Hadamard) | Konstante nach Abzug divergenter Hauptbeiträge in einer Regularisierung | $\zeta(-1)$ als endlicher Teil einer $\sum n$-Regularisierung |
| **Regulierte Summe** | Grenzwert nach Subtraktion der ersten $m$ Euler–Maclaurin-Terme | $\zeta(s)$ für $\Re(s)$ klein |
| **Algebraische Schließung** | Ausdruck nur durch $B_{2n}$, $\pi$, rationale Zahlen | $\zeta(2n)$ in §6.2 |

Im Streifen-Gedankenexperiment (§2–3) ist das **analytische Resultat** nach Bernoulli-Subtraktion: die **geschlossene** Spezialwertformel — nicht eine neue EABC-Observable.

### 8.4 Anbindung an Streifen-Produkt und $\zeta(s)/\zeta(2s)$

Aus Box 2: lokales Produkt $(1-p^{-s})(1+p^{-s})=1-p^{-2s}$. Global (formal, $\Re(s)>1$):

$$\frac{\zeta(s)}{\zeta(2s)}=\prod_{p}\frac{1-p^{-2s}}{1-p^{-s}}=\prod_{p}(1+p^{-s}).$$

An **ganzzahligen** $s=2n$ verbindet dies **gerade** Zeta-Spezialwerte mit einem **nur-Summen-Euler-Produkt** (Differenzfaktor weggekürzt). Die Bernoulli-Darstellung von $\zeta(2n)$ kontrolliert dann das **symmetrische** Streifen-Produkt $\prod_p(1+p^{-2n})$ indirekt über $\zeta(2n)$ und $\zeta(4n)$.

**Heuristik (Schicht C):** Die **Subtraktion** des antisymmetrischen Anteils $(1-p^{-s})$ im Quotienten $\zeta(s)/\zeta(2s)$ spiegelt die **Reduktion** von Differenz×Summe auf **reine Summe** — analog zum Übergang von $(e^x-e^{-x})(e^x+e^{-x})$ zu einer **einzigen** $\sinh(2x)$-Frequenz in Box 2. **Kein** Identitätsbeweis mit $D_E$ oder HL-Kanälen.

### 8.5 Bernoulli-Uhr und EABC (Grenze)

In der **Bernoulli-Uhr** (`collatz_offene_punkte.md`, `BernoulliClock.lean`) werden Tripel $(B_{2m-2},B_{2m},B_{2m+2})$ auf chiralen Zellen $\mathcal{E}_\pm$ (ABCE/CEAB) gelegt. Das ist **definitorische Geometrie** — **kein** numerischer Collatz-Befund, **keine** Lyapunov-Nutzung (Bernoulli-Normschale: No-Go in `collatz_generalangriff_2026.md`).

| Objekt | Schicht | Status |
|--------|---------|--------|
| $\zeta(2n)$, $B_{2n}$, Euler–Maclaurin | klassisch | **Theorem** |
| Bernoulli-Uhr-Tripel auf $\mathcal{E}_\pm$ | **C** | **Definition** (Gedankenmodell) |
| $D_E$, $\Phi_E$ aus Bernoulli-Subtraktion | — | **nicht behauptet** |

**Numerik (optional):** `collatz_eabc_zeta_bernoulli_check.py` verifiziert $\zeta(2n)$ gegen die Bernoulli-Formel in §6.2 für kleine $n$.

---

## 9. Fibonacci-Zeta: kontrolliertes Gegenmodell

$$\boxed{\;\text{Fibonacci-Zeta = geordneter Resonator; Riemann-Zeta = arithmetisch gestörter Resonator — \textbf{kontrolliertes Gegenmodell / Testfeld}, nicht „näher an Riemann“.}\;}$$

**Schlüsselunterscheidung (Schicht C):**

| Objekt | Spektrum | Lesart |
|--------|----------|--------|
| $\zeta(s)=\sum_{n\ge1} e^{-s\log n}$ | **arithmetisch irregulär** ($\log n$, Prim-Marken $\log p$) | gestörter Resonator |
| $\zeta_F(s)\approx 5^{s/2}\sum_{n\ge1} e^{-sn\log\varphi}$ | **goldener regulärer Frequenzkamm** ($n\log\varphi$ äquidistant) | kontrolliertes Gegenmodell |

Dieser Abschnitt fügt eine **dritte Skala** zum Dreieck aus §2–3 hinzu: neben $e^x$ (kontinuierlich) und $\zeta(s)$ (arithmetisch) die **goldene diskrete** Expansion/Kontraktion der Fibonacci-Zahlen $F_n$ — als **Testfeld**, in dem Resonanzstruktur explizit kontrollierbar ist.

### 9.1 Dreieck: kontinuierlich — arithmetisch — golden-diskret

| Skala | Objekt | Lesart |
|-------|--------|--------|
| **Kontinuierlich** | $e^x$ | reine Expansion/Kontraktion auf $\mathbb{R}$ |
| **Arithmetisch** | $\zeta(s)=\sum_{n\ge1} n^{-s}=\sum_{n\ge1} e^{-s\log n}$ | Exponentialsumme auf **irregulärer** Log-Achse ($\log n$, Prim-Marken $\log p$) |
| **Golden-diskret** | $F_n$ (Fibonacci) | **reguläre** diskrete Expansion/Kontraktion mit Basis $\varphi=(1+\sqrt5)/2$ |

**Label:** Tabelle = **Analogie** (Schicht **C**); $\zeta$-Definition = **Definition** (klassisch).

### 9.2 Binet-Formel und Matrix-Eigenwerte

Sei $\varphi=(1+\sqrt5)/2$ der goldene Schnitt und $\psi:=-\varphi^{-1}$. Dann gilt für $n\ge 0$ (**Binet**, klassisch):

$$\boxed{\;F_n=\frac{\varphi^n-\psi^n}{\sqrt5},\qquad \psi=-\varphi^{-1}.\;}$$

**Expansion vs. Kontraktion:** Der Term $\varphi^n$ wächst; der Term $\psi^n=(-\varphi^{-1})^n=(-1)^n\varphi^{-n}$ kontrahiert mit **Parität** $(-1)^n$. Die Fibonacci-Rekursion

$$F_{n+2}=F_{n+1}+F_n$$

entspricht der Potenz von

$$M=\begin{pmatrix}1&1\\1&0\end{pmatrix},\qquad
\begin{pmatrix}F_{n+1}\\F_n\end{pmatrix}=M^n\begin{pmatrix}1\\0\end{pmatrix},$$

mit Eigenwerten $\lambda_\pm=\varphi,\psi$ — dieselbe **Differenz zweier Exponentialmoden** wie in $\sinh/\cosh$, nun auf $\mathbb{Z}_{\ge0}$.

**Label:** Binet-Formel, Matrix-Eigenwerte = **Theorem** (klassisch).

**Numerik (optional):** `collatz_eabc_zeta_fibonacci_check.py` verifiziert $F_n$ (Binet vs. $M^n$) für $n\le 20$.

### 9.3 Fibonacci-Zeta $\zeta_F(s)$ — fast geometrisch, parity-gestört

**Definition** (formal, Konvergenz für $\Re(s)$ hinreichend groß):

$$\boxed{\;\zeta_F(s):=\sum_{n\ge1} F_n^{-s}.\;}$$

Mit Binet und $F_n\sim \varphi^n/\sqrt5$ folgt asymptotisch

$$\boxed{\;\zeta_F(s)\approx 5^{s/2}\sum_{n\ge1} e^{-sn\log\varphi}\;=\;5^{s/2}\,\frac{\varphi^{-s}}{1-\varphi^{-s}}\;}$$

— eine **geometrische Zeta** auf dem **regulären** Gitter $n\mapsto n\log\varphi$. Resonanzen (Polstellen der Hauptreihe) liegen auf der **imaginären Achse** bei

$$\boxed{\;s_k=\frac{2\pi i k}{\log\varphi},\qquad k\in\mathbb{Z}\setminus\{0\}\;}$$

(vertikales Resonanzgitter, goldener Fourier-Kamm $\omega_\varphi=\log\varphi$).

**Binet-Korrektur** — parity-gestörte geometrische Zeta:

$$\boxed{\;F_n^{-s}=5^{s/2}\,\varphi^{-ns}\,\bigl(1-(-1)^n\varphi^{-2n}\bigr)^{-s}.\;}$$

Der Faktor spiegelt den **Kontraktionsterm** $\psi^n=(-\varphi^{-1})^n=(-1)^n\varphi^{-n}$: Parität $(-1)^n$, Spiegel $\varphi^{-2n}$, komplexe Gewichtung $s$.

**Geboxte Struktur** (goldene Version von $e^x\pm e^{-x}$):

$$\boxed{\;\underbrace{\varphi^{-ns}}_{\text{Hauptterm (Expansion)}}\;+\;\underbrace{(-1)^n\varphi^{-2n}}_{\text{Spiegelterm (Kontraktion mit Parität)}}\;}$$

— formal parallel zu $(e^x+e^{-x})$ und $(e^x-e^{-x})$ in Box 2; $\zeta_F$ ist damit ein **diskretes golden-paritäts-geladenes Modell** für Zeta-artige Exponentialsummen.

**Label:** $\zeta_F$-Definition = **Definition**; asymptotische Entwicklung, Resonanzgitter = **Theorem** (klassisch); Gegenmodell-Lesart = **Hypothese/Modell** (Schicht **C**).

### 9.4 Vergleich auf der kritischen Linie — Gegenmodell, nicht Annäherung

$$\boxed{\;\text{Fibonacci-Zeta = geordneter Resonator}\;\|\;\text{Riemann-Zeta = arithmetisch gestörter Resonator.}\;}$$

| Aspekt | Riemann $\zeta(s)$ | Fibonacci-$\zeta_F(s)$ |
|--------|-------------------|------------------------|
| Gitter | $n\mapsto\log n$ (**irregulär**, Prim-Lücken) | $n\mapsto n\log\varphi$ (**regulär**, äquidistant) |
| „Frequenzkamm“ | $\omega_p=\log p$ (Prim-Fourier-Kamm) | $\omega_\varphi=\log\varphi$ (goldener Kamm) |
| Lokale Faktoren | Euler $(1-p^{-s})(1+p^{-s})$ | Binet-Haupt- + Spiegelterm (Paritätsstörung) |
| Epistemische Rolle | **Zielobjekt** (arithmetisch) | **kontrolliertes Gegenmodell / Testfeld** |

**Kritische Linie** ($s=\tfrac12+it$): Vergleich $\zeta(\tfrac12+it)$ vs. $\zeta_F(\tfrac12+it)$ — gedämpfter goldener Fourier-Kamm ($\varphi^{-n/2}$, Phase $e^{-itn\log\varphi}$) gegen **arithmetisch gestörten** Prim-Kamm ($p^{-1/2}$, Phase $e^{-it\log p}$). Die Frage ist **nicht** „liegt $\zeta_F$ näher an Riemann?“, sondern: **welche Observable unterscheidet regulären von gestörtem Kamm?**

**Nullmodell-Hierarchie** (`collatz_eabc_zirkulationshypothese.md` §4.8.2): $\zeta_F$ ist **Stufe 0** — das **regulärste** Referenzgitter ($\Sigma_A^{\mathrm{comb}}$ bzw. goldener Log-Kamm), **vor** Stufe 1 (Permutation). Prim-Log-Geometrie ist die **Störung** darüber; HL-Konsistenz (Stufe 3) testet arithmetische Realität. **Kein Theorem** über Identität mit Prim-Korrelationen.

**Label:** Gegenmodell / Stufe-0-Baseline = **Hypothese** / **Modell** (Schicht **C**).

### 9.5 Kritische Linie und Resonanzgitter

Auf dem Fibonacci-Gitter ($s=\tfrac12+it$):

$$\varphi^{-n(\tfrac12+it)}=\varphi^{-n/2}\,e^{-it\,n\log\varphi}.$$

Der **reelle Halbzahl-Anteil** $\varphi^{-n/2}$ und die **Phase** $e^{-it n\log\varphi}$ faktorisieren wie $p^{-1/2}\cdot p^{-it}$ auf der Prim-Achse — **Bildsprache** für eine „goldene kritische Linie“, **kein** RH-Anspruch.

**Resonanzgitter** (imaginäre Achse): $t_k=2\pi k/\log\varphi$ markiert die **vertikalen Resonanzlinien** der Hauptgeometrie $\varphi^{-s}/(1-\varphi^{-s})$. Auf $\Re(s)=\tfrac12$ liegen diese Punkte **diagnostisch** — sie strukturieren den Vergleich $\zeta$ vs. $\zeta_F$, beweisen aber **nichts** über Nullstellen von $\zeta$.

**Label:** **Analogie** (Schicht **C**); Resonanzformel = **Theorem** (klassisch).

### 9.6 EABC/Wigner: chirale Paarung (Schicht C)

| Fibonacci / golden | EABC (Lesart) | Label |
|------------------|---------------|-------|
| $\varphi^n$ (Expansion) vs. $(-\varphi^{-1})^n$ (Kontraktion) | chirale Paarung $\gamma^+$ (ABCEA) vs. $\gamma^-$ (CEABC) | **Analogie** (C) |
| Parität $(-1)^n$ | Vorzeichenwechsel auf $C_4$-Zyklus | **Analogie** (C) |
| Haupt- + Spiegelterm (Binet) | ABCE $\leftrightarrow$ CEAB — diskrete Wigner-Struktur | **Analogie** (C) |
| Fibonacci-Fenster $F_k\le p<F_{k+1}$ | natürliches Renormierungsgitter für chirale Wigner-Zeugen | **Hypothese/Modell** (C) |

**Explizit ausgeschlossen:** Fibonacci beweist **weder** die Riemann-Hypothese **noch** EABC-Holonomie oder $D_E$-Wachstum. Fibonacci ist ein **Renormierungsgitter** für chirale Wigner-Zeugen — **Interpretationsschicht**, nicht Beweiskette.

**Label:** gesamter §9.6 = **Analogie** / **Modell** (Schicht **C**).

### 9.7 Drei numerische Zeugen (Experiment, Schicht C)

Implementierung: `eabc_zeta_fibonacci_witnesses.py` → `eabc_zeta_fibonacci_witnesses.json`

| Zeuge | Observable | Lesart |
|-------|------------|--------|
| **1** | $t_k=2\pi k/\log\varphi$: $\|\zeta(\tfrac12+it_k)\|$ vs. zufällige Höhen; parallel $\|\zeta_F(\tfrac12+it_k)\|$ | Resonanzgitter auf der kritischen Linie — diagnostisch |
| **2** | $\theta_p=(\log p)/(\log\varphi)\bmod 1$ für $p\le B$; Histogramm / $\chi^2$-Test | Fibonacci-Sampling der Prim-Frequenzen — Uniformität = Nullbefund |
| **3** | In Fenstern $F_k\le p<F_{k+1}$: Zählung ABCE vs. CEAB (mod-$12$-4-Gramm auf aufeinanderfolgenden Primzahlen), $D_k=(N_{\mathrm{ABCE}}-N_{\mathrm{CEAB}})/(N_{\mathrm{ABCE}}+N_{\mathrm{CEAB}})$ | EABC-Fibonacci-Fenster — minimale Wigner-Zählung |

```bash
python3 eabc_zeta_fibonacci_witnesses.py --prime-bound 50000
pytest tests/test_eabc_zeta_fibonacci_witnesses.py -q
```

**Label:** gesamter §9.7 = **Experiment** (Schicht **C**); **kein Theorem**, **kein RH-Beweis**.

### 9.8 Meromorphe Normalform und Resonanztürme (klassisch)

Aus $F_n^{-s}=5^{s/2}\,\varphi^{-ns}\,(1-(-1)^n\varphi^{-2n})^{-s}$ folgt die **exakte meromorphe Normalform** durch Binomialentwicklung $(1-u)^{-s}=\sum_{m\ge0}(s)_m/m!\,u^m$ mit $u=(-1)^n\varphi^{-2n}$:

$$\boxed{\;\zeta_F(s)=5^{s/2}\sum_{m\ge0}\frac{(s)_m}{m!}\,\frac{(-1)^m\,\varphi^{-(s+2m)}}{1-(-1)^m\,\varphi^{-(s+2m)}}\;}$$

| $m$ | Term | Resonanz / Polstruktur (imaginäre Achse) |
|-----|------|------------------------------------------|
| **0** | $5^{s/2}\,\varphi^{-s}/(1-\varphi^{-s})$ | **Hauptresonator** — $s=2\pi ik/\log\varphi$ |
| **gerade $m>0$** | Paritäts-Spiegelterm | $s=-2m+2\pi ik/\log\varphi$ |
| **ungerade $m$** | Paritäts-Spiegelterm (Phasenverschiebung) | $s=-2m+(2k+1)\pi i/\log\varphi$ |

**Lesart:** Vertikale **Resonanztürme** — $m=0$ ist der goldene Hauptkamm; $m>0$ sind **parity-verschobene Spiegelresonatoren** (Binet-Korrekturen). Saubere mathematische Version von Expansion/Kontraktion/Parität.

**Label:** Box + Tabelle = **Theorem** (klassische Analysis); EABC-Kopplung = **nicht** enthalten.

**Numerik:** `eabc_zeta_fibonacci_witnesses.py` (Zeuge 4) vergleicht direkte Partialsumme $\sum F_n^{-s}$ mit meromorpher Partialsumme ($m\le m_{\max}$).

### 9.9 Goldene Fourier-Zeugen und mod-210-Fibonacci-Schalen (Experiment)

**Goldene Koordinaten** auf Primvierlingen $p$ (Start von $(p,p+2,p+6,p+8)$):

$$\theta_\varphi(p)=\frac{\log p}{\log\varphi}\bmod 1,\qquad
\chi(p)=\begin{cases}+1& p\equiv 5\pmod{12}\;\text{(ABCE)}\\ -1& p\equiv 11\pmod{12}\;\text{(CEAB)}\end{cases}$$

$$\boxed{\;C_m(N)=\sum_{p\le N}\chi(p)\,e^{2\pi i m\,\theta_\varphi(p)},\qquad
Z_m(N)=\frac{|C_m(N)|}{\sqrt{Q(N)}}\;}$$

mit $Q(N)$ = Anzahl klassifizierter Vierlinge $\le N$.

**Fibonacci-Schalen:** $F_k\le p<F_{k+1}$; je Schale mod-210-Triple $(D_{11,k},D_{101,k},D_{191,k})$ mit $D_{r,k}=N_{\mathrm{ABCE},r}-N_{\mathrm{CEAB},r}$ für $r\in\{11,101,191\}$ (Wigner-Zellen, vgl. `eabc_wigner_zellen.py`, `eabc_quadruplets.csv`).

**Zentralfrage (Experiment):** Bleibt die Signatur $(+,+,-)$ — d. h. $D_{11}>0$, $D_{101}>0$, $D_{191}<0$ — entlang Fibonacci-Schalen **stabiler** als in linearen Fenstern gleicher Breite?

| Komponente | Schicht | Label |
|------------|---------|-------|
| Meromorphe Normalform, Resonanztürme | — | **Theorem** (klassisch) |
| $\theta_\varphi$, $C_m$, $Z_m$ auf Vierlingen | **C** | **Experiment** |
| mod-210-Signatur vs. Fibonacci-Schalen | **C** | **Hypothese** / **Experiment** |
| Kopplung $\zeta_F$-Resonanz ↔ EABC | **C** | **Analogie** — **nicht** behauptet |

**Nullmodell-Hierarchie** (`collatz_eabc_zirkulationshypothese.md` §4.8.2): $\zeta_F$ und der goldene Log-Kamm sind **Stufe 0** (reguläres Referenzgitter); die Fourier-Zeugen $C_m,Z_m$ und mod-210-Schalen testen, ob **arithmetische** Vierlingsordnung dieselbe chirale Signatur trägt wie das ideale goldene Gitter — **ohne** Theorem über $\Sigma_A$ oder HL.

```bash
python3 eabc_zeta_fibonacci_witnesses.py --quadruplet-bound 500000 --m-max 8
pytest tests/test_eabc_zeta_fibonacci_witnesses.py -q
```

Implementierung: Zeuge 4–6 in `eabc_zeta_fibonacci_witnesses.py` → `eabc_zeta_fibonacci_witnesses.json`.

**Label:** gesamter §9.9 = **Experiment** (Schicht **C**).

**Experiment (N=10⁸, 2025-06):** `eabc_zeta_fibonacci_witnesses.py --quadruplet-bound 100000000 --m-max 8` auf $Q(10^8)=4768$ klassifizierten Vierlingen. Goldene Fourier-Zeugen: $Z_0=0{,}695$, $Z_1=0{,}737$, $Z_2=1{,}023$, …, $Z_8=0{,}490$ — $Z_1$ fällt gegenüber $N=5\cdot10^6$ ($Z_1\approx2{,}03$) auf $\mathcal{O}(1)$ zurück (kein persistenter Peak). **Aggregat** mod-210 über alle Vierlinge: $(D_{11},D_{101},D_{191})=(+44,+22,-19)$, Signatur **$(+,+,-)$** (bei $5\cdot10^6$ noch $(+,+,+)$). Fibonacci-Schalen: $2/28$ Schalen mit $(+,+,-)$ ($7{,}1\%$) vs. lineare Fenster gleicher Breite $1/23$ ($4{,}3\%$); per-Shell-Stabilität bleibt schwach, das Aggregat aligniert mit `eabc_quadruplets`-Skala. **Label:** Experiment (Schicht **C**) — diagnostisch, kein Theorem.

### 9.10 Kurzform

$$\boxed{\;\text{Fibonacci-Zeta = geordneter Resonator; Riemann-Zeta = arithmetisch gestörter Resonator — kontrolliertes Gegenmodell, kein „näher an Riemann“.}\;}$$

| Aussage | Schicht | Tao-Label |
|---------|---------|-----------|
| Binet, Matrix-Eigenwerte, $\zeta_F$-Definition, Resonanzgitter | — | **Theorem** / **Definition** |
| $\zeta_F$ als Stufe-0-Gegenmodell (§4.8.2) | **C** | **Hypothese** / **Modell** |
| Streifen-, kritische-Linie-, EABC-Paarung | **C** | **Analogie** |
| Drei numerische Zeugen (§9.7) | **C** | **Experiment** |
| Meromorphe Normalform (§9.8) | — | **Theorem** (klassisch) |
| Goldene Fourier-Zeugen, mod-210-Schalen (§9.9) | **C** | **Experiment** / **Hypothese** |
| RH oder EABC aus Fibonacci | — | **nicht behauptet** |

---

## 10. Kurzform (Erweiterung)

$$\boxed{\;\mathcal{S}=\mathbb{Z}\cup(\tfrac12+\mathbb{Z})\;\text{ — halbzahlige Gitterlesart für Streifen-Exponenten.}\;}$$

$$\boxed{\;\zeta(2n)=(-1)^{n+1}\,\dfrac{B_{2n}(2\pi)^{2n}}{2(2n)!},\qquad \zeta(1-2n)=-\dfrac{B_{2n}}{2n}\;\text{(klassisch).}\;}$$

$$\boxed{\;\text{Bernoulli-Subtraktion}=\text{Euler–Maclaurin-/Funktionalgleichungs-Kompensation}\Rightarrow\text{ endliche Spezialwerte.}\;}$$

$$\boxed{\;\mathbb{R}\subset\mathbb{C}\subset\mathbb{H}\subset\mathbb{O}\;\Leftrightarrow\;\text{Schicht-C-Ikone: }C_4\text{-Schatten }\mathbb{H}\text{, Oktanion-Programm ohne B-Theorem.}\;}$$

| Aussage | Schicht | Tao-Label |
|---------|---------|-----------|
| §6.2, §8.1–8.2 (klassisch) | — | **Theorem** |
| §6 Gitter, §8.4 Streifen-Quotient | **C** | **Analogie** |
| §7 Hurwitz/Oktanionen | **C** | **Analogie** / **Ikone** |
| Bernoulli-Uhr | **C** | **Definition** (Gedankenmodell) |
| §9 Fibonacci-Zeta | **C** (+ klassisch) | **Hypothese** / **Analogie** |

$$\boxed{\;\zeta_F(s)=\sum_{n\ge1}F_n^{-s}\;\Leftrightarrow\;\text{goldener Log-Kamm }n\log\varphi\text{ mit Binet-Streifen }\varphi^{-ns}+(-1)^n\varphi^{-2n}.\;}$$

---

*Verknüpft mit `collatz_eabc_zirkulationshypothese.md` §4.8.2 (Stufe 0: $\zeta_F$), §4.9; epistemischer Rahmen: `collatz_eabc_epistemik_schichten.md`; Bernoulli: `collatz_bernoulli_schalen.pdf`, `BernoulliClock.lean`; Oktanionen: `divisionsalgebren.tex`, `collatz_oktonionen_beweis.pdf`; Fibonacci-Numerik: `collatz_eabc_zeta_fibonacci_check.py`, `eabc_zeta_fibonacci_witnesses.py` (§9.7–9.9); Wigner-Zellen: `eabc_wigner_zellen.py`, `eabc_quadruplets.csv`.*
