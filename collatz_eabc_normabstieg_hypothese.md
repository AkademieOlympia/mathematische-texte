# EABC-Normabstiegs-Hypothese

**Stand:** Juni 2026 · **Kanonsiche Formulierung** der Normabstiegs-Hypothese über die
Hurwitz-Divisionsalgebren. Tao-Labels: Definition | Theorem | Conjecture | Heuristik |
Forschungsfrage | Forschungsvision | Experiment.

**Epistemische Warnung:** Kein Collatz-Beweis. Kein ERPC.

**Querverweise:**
- `collatz_eabc_euklidische_hebung.md` — algebraischer Kontext, Hurwitz-Kette (PR #54)
- `collatz_eabc_invarianzprogramm.md` — EABC-Fluktuationsfeld $\delta(x)$, $\kappa$-Klassifikation
- `collatz_generalangriff_2026.md` — strategischer Pointer
- `collatz_eabc_gauss_defekt_test.py` — Experiment zur trivialen Teilhypothese (§8)
- `collatz_eabc_gauss_faktor_eabc_test.py` — Experiment Z[i]-Faktor → EABC (§9)

---

## 1. Divisionsalgebren, Gitter und Norm

Sei $A\in\{\mathbb{R},\mathbb{C},\mathbb{H},\mathbb{O}\}$ eine **normierte Divisionsalgebra** über
$\mathbb{R}$ (Hurwitz-Kette). Jedem $A$ ist ein **Gitter** $\Lambda_A\subset A$ zugeordnet
(Maximalordnung, wo definiert):

| Algebra | Gitter $\Lambda_A$ | Norm $N:A\to\mathbb{R}_{\ge 0}$ |
|---------|-------------------|--------------------------------|
| $\mathbb{R}$ | $\Lambda_{\mathbb{R}}=\mathbb{Z}$ | $N(x)=|x|$ |
| $\mathbb{C}$ | $\Lambda_{\mathbb{C}}=\mathbb{Z}[i]$ | $N(a+bi)=a^2+b^2$ |
| $\mathbb{H}$ | $\Lambda_{\mathbb{H}}=\mathbb{H}_{\mathrm H}$ (Hurwitz) | $N(q)=\sum a_i^2$ |
| $\mathbb{O}$ | $\Lambda_{\mathbb{O}}=\mathbb{O}_{\mathrm H}$ (Hurwitz-Oktanionen) | $N$ multiplikativ |

**Theorem (Hurwitz).** Die Norm ist **multiplikativ** auf allen vier Algebren:
$N(xy)=N(x)N(y)$.

**Spezialfall $\mathbb{Z}$.** Für $A=\mathbb{R}$, $\Lambda_A=\mathbb{Z}$ ist $N$ die Betragsnorm;
dies ist der **eindimensionale Spezialfall** der universellen Schablone.

**Label:** Definition (Gitter, Norm); Multiplikativität = **Theorem**.

---

## 2. Projektion, Quotient und Hebungsschritt

Sei $\Pi_A:A\to\Lambda_A$ eine **Projektion** (nächster-Gitter-Punkt bzw. Quotientenwahl im
euklidischen Spezialfall). Für $x,y\in\Lambda_A$, $y\neq 0$, $y$ invertierbar in $A$:

\[
z := x\,y^{-1}\in A,\qquad
q := \Pi_A(z)\in\Lambda_A,\qquad
r := x - q\,y\in\Lambda_A.
\]

Der **Hebungsschritt** (euklidischer Operator) ist
\[
E_A : \Lambda_A\times(\Lambda_A\setminus\{0\}) \longrightarrow \Lambda_A\times\Lambda_A,
\qquad
E_A(x,y) := (y,\, r).
\]

Im Fall $A=\mathbb{R}$, $\Lambda_A=\mathbb{Z}$: $q=\lfloor x/y\rfloor$ (mit Restkorrektur),
$r=x-qy$ — der **klassische euklidische Algorithmus**.

**Label:** **Definition**.

---

## 3. Zulässigkeit und euklidischer Normabstieg

Ein Schritt $E_A(x,y)=(y,r)$ heißt **zulässig** (euklidisch), wenn
\[
N(r) < N(y).
\]

Die Iteration $E_A^{(k)}$ erzeugt eine **Abstiegskette der Restnormen**; in euklidischen Ringen
terminiert sie beim $\gcd$.

| Algebra | $N(r)<N(y)$ | Label |
|---------|-------------|-------|
| $\mathbb{Z}$ | **Theorem** | klassisch |
| $\mathbb{Z}[i]$ | **Theorem** | Gauß-Euklid |
| $\mathbb{H}_{\mathrm H}$ | **Theorem** | Hurwitz-Euklid |
| $\mathbb{O}_{\mathrm H}$ | **Forschungsfrage** | kein assoziativer $\gcd$ |

**Label:** Zulässigkeit = **Definition**; $N(r)<N(y)$ = **Theorem** (R,C,H) bzw. **Forschungsfrage** (O).

---

## 4. Defektoperator und iterierte Defektreduktion

### Definition 1 (Schritt-Defekt)

\[
D_A(x,y) := r = x - \Pi_A(x\,y^{-1})\,y \in \Lambda_A.
\]
Der **Defekt** eines Hebungsschritts ist der **Rest** $r$, nicht die Divisibilitätsrelation
$x\mid y$.

### Definition 2 (Gitter-Defekt)

\[
\widetilde{D}_A(x) := x - \Pi_A(x) \in A.
\]
Misst die Abweichung von $x$ vom Gitter (Peano-Projektion in $\mathbb{R}$).

### Definition 3 (Lokaler Normdefekt)

\[
\varepsilon(x) := N(x) - N(\Pi_A(x)), \qquad x\in\Lambda_A
\]
(oder $x\in A$ vor Gitterprojektion). Auf dem Gitter $\Lambda_A$ gilt $\varepsilon(x)=0$;
der interessante Wert liegt im **Defektfeld** $x-\Pi_A(x)$ und dessen Normstatistik.

**Heuristik.** Die **iterierte Defektreduktion** unter $E_A$ ist das primäre Objekt — **nicht**
Divisibilität oder Faktorzerlegung allein. In $\mathbb{O}$ bleibt Normabstieg möglich, während
klassischer $\gcd$ wegfällt.

**Label:** **Definition**; Defektreduktion statt Divisibilität = **Heuristik / Forschungsvision**.

---

## 4a. Universelle Invariante über R, C, H, O

| Algebra | Ordnung / Struktur | Was „funktioniert“ euklidisch |
|---------|-------------------|------------------------------|
| $\mathbb{R}$ | totale Ordnung | **Normabstieg** $|r|<|y|$ |
| $\mathbb{C}$ | Norm + Rotation | **Normabstieg** $N(r)<N(y)$ |
| $\mathbb{H}$ | Norm + Rotation (nicht-kommutativ) | **Normabstieg** $N(r)<N(y)$ |
| $\mathbb{O}$ | Norm + Alternativität | **Normabstieg** (offen); kein assoziativer $\gcd$ |

**Durchbruch (Heuristik):** Nicht $D_A(x,y)$ allein, sondern die Frage: **Welche universelle
Größe bleibt über R, C, H, O erhalten?** Antwort-Kandidat: **multiplikative Norm** $N$ und
**Normabstieg** unter dem Hebungsoperator $E_A$. Euklidisches $\mathbb{Z}$ ist der
1D-Spezialfall dieser Kette.

**Label:** **Heuristik / Forschungsvision**.

---

## 4b. Bernoulli-Lücke und Defektfeld (Heuristik)

Bernoulli-Summation misst Diskretisierungsdefekte:
\[
\sum_{n=a}^{b} f(n) - \int_a^b f(x)\,dx.
\]
**Heuristik:** Peano-Diskretisierung, Bernoulli-Kontinuumslücke und EABC-Orientierungen sind
drei Gesichter eines **Defektfeldes**; lokale Normdefekte $\varepsilon(x)$ sind die
arithmetische Feinauflösung auf dem Gitter $\Lambda_A$.

**Label:** **Heuristik** (kein etablierter Isomorphismus).

---

## 4c. Einheitliches Bild (Forschungsvision)

$$\text{Peano (diskret)} \;\longleftrightarrow\; \text{Bernoulli (Kontinuumslücke)}
\;\longleftrightarrow\; \text{EABC (Orientierungen)} \;\longleftrightarrow\;
\text{Normabstieg (lokale Defekte } \varepsilon).$$

**Label:** **Forschungsvision**.

---

## 5. Norm-irreduzible Elemente

$u\in\Lambda_A\setminus\{0\}$ heißt **norm-irreduzibel**, wenn aus $u=vw$ mit
$N(u)=N(v)N(w)$ stets $N(v)=1$ oder $N(w)=1$ folgt (bis auf Einheiten).

**Theorem.** In $\mathbb{Z}$ sind die norm-irreduziblen Elemente genau die **Primzahlen**
und ihre Negationen — d.\,h. Irreduzible **bis auf Einheiten** $\pm 1$.

**Label:** **Definition** (norm-irreduzibel); Spezialfall $\mathbb{Z}$ = **Theorem**.

---

## 6. EABC und minimale irreduzible Normdefekte

### Conjecture (EABC-Primdefekt)

In der **reellen Spezialisierung** ($A=\mathbb{R}$, $\Lambda_A=\mathbb{Z}$) gilt heuristisch:
\[
p\ \text{prim in }\mathbb{Z}
\quad\stackrel{?}{\Longleftrightarrow}\quad
\widetilde{D}_{\mathbb{R}}(p)\ \text{minimal irreduzibler Normdefektzustand}.
\]

In höheren Algebren entstehen **Defektklassen**, die auf die vier Hurwitz-Basisrichtungen
projizieren:
\[
(E,A,B,C) \longleftrightarrow (1,\, i,\, j,\, k)
\]
(vgl. `collatz_eabc_euklidische_hebung.md` §3, `PAPER_HURWITZ_RESONANZ.md`).

**Label:** Prim = minimaler Defekt in $\mathbb{Z}$ = **Conjecture / Forschungsvision**;
EABC-Kodierung = **Definition**.

---

## 7. Hauptkonjektur: einheitlicher Normabstiegs-Rahmen

### Conjecture (EABC-Normabstiegs-Hypothese)

Es existiert ein **einheitlicher Normabstiegs-Rahmen** über $A\in\{\mathbb{R},\mathbb{C},\mathbb{H},\mathbb{O}\}$:
Gitter $\Lambda_A$, multiplikative Norm $N$, Projektion $\Pi_A$, Schritt-Defekt $D_A(x,y)$ und
Operator $E_A$, sodass:

1. $\mathbb{Z}$ der 1D-Spezialfall ist;
2. Primzahlen in $\mathbb{Z}$ als **eindimensionale Schatten** höherdimensionaler
   norm-irreduzibler Defekte lesbar sind;
3. die native **EABC-Klassifikation** $\kappa:\mathbb{P}_{>3}\to\{E,A,B,C\}$ (mod $12$) die
   **sichtbare Projektion** dieses Rahmens auf die Peano-Achse ist.

**Label:** **Conjecture** (Gesamtrahmen); mod-$12$-EABC = **Definition** (`collatz_eabc_invarianzprogramm.md`).

---

## 8. Testbare Teilhypothese: Gauß-Zerlegung und EABC

### Theorem (Zerlegungsverhalten in $\mathbb{Z}[i]$)

Für eine rationale Primzahl $p$:

- $p\equiv 1\pmod 4$ $\Rightarrow$ $p$ **zerlegt** (split) in $\mathbb{Z}[i]$;
- $p\equiv 3\pmod 4$ $\Rightarrow$ $p$ **träge** (inert) in $\mathbb{Z}[i]$;
- $p=2$ ramifiziert (Sonderfall, außerhalb $\mathbb{P}_{>3}$).

**Label:** **Theorem** — klassische algebraische Zahlentheorie.

### Conjecture / Heuristik (Gauß–EABC-Brücke)

Die **erste konkrete Brücke** zwischen höherdimensionaler Normstruktur und EABC:

\[
p\equiv 1\pmod 4 \quad\stackrel{?}{\longleftrightarrow}\quad \kappa(p)\in\{E,A\}
\quad (4n+1\text{-Primzahlen}),
\]
\[
p\equiv 3\pmod 4 \quad\stackrel{?}{\longleftrightarrow}\quad \kappa(p)\in\{B,C\}
\quad (4n+3\text{-Primzahlen}).
\]

**Epistemische Einordnung:**

| Aussage | Status |
|---------|--------|
| Split/inert in $\mathbb{Z}[i]$ für $p\equiv 1,3\pmod 4$ | **Theorem** |
| Grobe Zuordnung split $\leftrightarrow E\cup A$, inert $\leftrightarrow B\cup C$ | Für $p>3$: **arithmetisch exakt** (mod-$12$-Restklassen); als *Bedeutungs*-Brücke **Heuristik** |
| Feinzuordnung split $\to E$ vs.\ $A$ (bzw.\ inert $\to B$ vs.\ $C$) | **nicht** durch Split/inert allein bestimmt — offen |

**Experiment:** `collatz_eabc_gauss_defekt_test.py` $\to$ `collatz_eabc_gauss_defekt.json`.

**Epistemische Warnung:** Die grobe bipartite Zuordnung split $\leftrightarrow E\cup A$ ist
**arithmetisch trivial** (mod $4$ × mod $12$). Sie beweist keine geometrische EABC-Brücke.

---

## 9. Stärkere Konjektur und Forschungspriorität Z[i]

### Conjecture (Defektfeld-Singularitäten)

$$\boxed{
\text{Primzahlen sind Singularitäten eines universellen Defektfeldes.}
}$$

Richtung der Kausalität: **Defekte $\to$ Primzahlen**, nicht umgekehrt. In $\mathbb{Z}[i]$
„sieht“ das Defektfeld die **Geometrie**; die Projektion auf $\mathbb{Z}$ liefert $4n\pm 1$;
EABC ist die **feinere Auflösung** derselben Orientierungsgeometrie auf der Peano-Achse.

**Label:** **Conjecture / Forschungsvision**.

### Forschungsfrage (Experiment, nicht-trivial)

Für $p\equiv 1\pmod 4$, $p=a^2+b^2$ (kanonisch $0<a\le b$): projiziere $(a,b)$ mod $12$ auf
EABC-Klassen ($\kappa$ via `eabc_from_lean`, nur für Reste $\in\{1,5,7,11\}$).

**Frage:** Ist die Verteilung $(a,b)\mapsto(E,A,B,C)$ **asymmetrisch**? Erscheinen
Vierlings-Biases oder Chiralitäten $\Rightarrow$ erste echte Brücke $\mathbb{Z}[i]\to$EABC.

**Experiment:** `collatz_eabc_gauss_faktor_eabc_test.py` $\to$
`collatz_eabc_gauss_faktor_eabc.json`.

**Strukturelle Voraussage:** Für $p\equiv 1\pmod 4$ ist genau eine Faktorleg gerade
$\Rightarrow$ höchstens **eine** EABC-sichtbare Leg pro Paar; beide gleichzeitig unmöglich.
Der Test prüft $\kappa(p)$ vs.\ ungerade Leg, gerade-Leg-mod-$12$-Verteilung und Shuffle-Null.

**Priorität:** $\mathbb{Z}[i]$ vor $\mathbb{O}$ — assoziatives Gitter, Gauß-Euklid etabliert,
direkter Primideal-Zugang.

---

## 10. Forschungsprogramm (vier Schritte)

1. **Formalisation** von $E_A$, $D_A(x,y)$, $\widetilde{D}_A(x)$ und $\varepsilon(x)$ für
   $A\in\{\mathbb{R},\mathbb{C},\mathbb{H}\}$ (Python-Stub `collatz_eabc_euklid_hebung.py`;
   Lean-Ziel in `collatz_eabc_core`).
2. **Gauß–EABC-Test (§8):** empirische Korrelation split/inert mit $\kappa$ — ehrlich als
   trivial einstufen; **Z[i]-Faktortest (§9)** als nicht-triviale Fortsetzung.
3. **Hurwitz-Erweiterung:** Defektklassen in $\mathbb{H}_{\mathrm H}$ und Vergleich mit
   EABC-Primideal-Resonanz (`PAPER_HURWITZ_RESONANZ.md`).
4. **Oktanionen / Projektion:** kanonische $\Pi_{\mathbb{O}}$, Defekt-Abstieg ohne assoziativen
   $\gcd$; Verknüpfung zum EABC-Fluktuationsfeld (`collatz_eabc_invarianzprogramm.md`).

**Label:** **Forschungsprogramm** — kein Theorem-Anspruch.

---

## 11. Zusammenfassung

$$\boxed{
\text{Primzahlen sind Singularitäten eines universellen Defektfeldes}
\quad\text{(Conjecture);}
$$
$$\boxed{
\text{Primzahlen sind eindimensionale Schatten irreduzibler Normdefekte.}
\quad\text{(Spezialfall $\mathbb{Z}$).}
$$

Die **EABC-Normabstiegs-Hypothese** behauptet einen einheitlichen algebraischen Rahmen
$(\Lambda_A,N,\Pi_A,E_A,D_A,\varepsilon)$ über die Hurwitz-Kette, in dem die mod-$12$-EABC-Klassen
die sichtbare 1D-Projektion höherdimensionaler Defektstruktur sind. Die Gauß–EABC-Brücke (§8)
ist **arithmetisch trivial**; der **Z[i]-Faktortest** (§9) ist der erste nicht-triviale
empirische Zugang.

---

*Epistemische Einordnung: Theoreme zu Euklidizität und Gauß-Zerlegung sind etabliert;
EABC-Normabstiegs-Hypothese und Prim-als-Defekt sind explizit offen.*
