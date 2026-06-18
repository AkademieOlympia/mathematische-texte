# Quaternionen-EABC-Maßhypothese ($\Sigma_p$, $\mu_p$, $I(\mu_p)$)

**Status:** Forschungshypothese + Experiment  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Experiment:** `collatz_eabc_hurwitz_orbit_test.py` → `collatz_eabc_hurwitz_orbit.json`  
**Σ→p-Defekt:** `collatz_eabc_shell_defekt_test.py` → `collatz_eabc_shell_defekt.json`  
**Zerlegungsregimen:** `collatz_eabc_Z_decomposition_test.py` → `collatz_eabc_Z_decomposition.json`  
**Dirichlet-Erzeuger:** `collatz_eabc_dirichlet_D.py` → `collatz_eabc_dirichlet_D.json`  
**Tao-Labels:** Definition | Theorem | Conjecture | Heuristik | Experiment

**Querverweise:** `collatz_eabc_kommutator_assoziator.md` (Kommutator/Chiralität $\chi_p$; Assoziator-Spektrum auf $\mathbb{O}$) ·
`collatz_eabc_normabstieg_hypothese.md` · `collatz_eabc_euklidische_hebung.md` ·
`collatz_eabc_plattenuebergang.md` (kanonsiche $\Sigma_n$-Plattenhypothese, PR #54) ·
`collatz_eabc_oktonion_singularitaet.md` (8D-Forschungsprogramm, PR #54) ·
`collatz_eabc_zerlegungsregimen.md` ($Z(n)$, $\Delta Z$, operative Regimen-Zählung) ·
`collatz_eabc_invarianzprogramm.md` · `collatz_eabc_gauss_spaltung_hypothese.md` ·
`collatz_eabc_eisenstein_spaltung.md` · `collatz_eabc_bernoulli_uebersetzung.md` (Branch `collatz/eabc-bernoulli-sensor`) ·
`collatz_morley_tm_numerik.py` (Morley-Parallel, getrennte Spur)

> **Hinweis:** `collatz_eabc_hurwitz_spaltung.md` enthält die **kompakte** Orbit-Theorie §1–§10;
> dieses Dokument erweitert um $\Sigma_p$, $H_p$, $K_p$ und Dirichlet-Spur.

---

## 1. Normschale $\Sigma_p$ und $r_4(p)$

**Definition.** Sei $\mathbb{H}_{\mathrm H}$ die Hurwitz-Maximalordnung der Quaternionen und
$N:\mathbb{H}_{\mathrm H}\to\mathbb{Z}_{\ge 0}$ die **Norm**
$$N(q)=a^2+b^2+c^2+e^2\quad\text{für }q=a+bi+cj+ek.$$

Die **Normschale** (norm shell) der Primzahl $p$ ist
$$\Sigma_p := \{q\in\mathbb{H}_{\mathrm H} : N(q)=p\}.$$

**Theorem (Jacobi, Vier-Quadrate).** Für Primzahl $p$ gilt
$$r_4(p) := \#\{(a,b,c,e)\in\mathbb{Z}^4 : a^2+b^2+c^2+e^2=p\} = 8(p+1).$$

**Theorem (Hurwitz-Schalengröße).** Für ungerade Primzahlen $p>2$:
$$|\Sigma_p| = 3\,r_4(p) = 24(p+1).$$
Für $p=2$: $|\Sigma_2|=r_4(2)=24$.

**Label:** $\Sigma_p$, $N$ = **Definition**; $r_4(p)=8(p+1)$ = **Theorem**; $|\Sigma_p|$ = **Theorem** + Zählung.

---

## 2. Einheiten $U_{\mathrm H}$, Orbits $O_i$ — Schale vs. Bahn

**Theorem (Hurwitz-Einheiten).**
$$U_{\mathrm H} := U(\mathbb{H}_{\mathrm H}),\qquad |U_{\mathrm H}|=24.$$

**Definition (Doppelbahn).** Für $q\in\Sigma_p$:
$$O(q) := \{u\,q\,v : u,v\in U_{\mathrm H}\}.$$

**Definition (Bahnzerlegung).**
$$\Sigma_p = \bigsqcup_i O_i.$$

| Objekt | Symbol | Bedeutung |
|--------|--------|-----------|
| Normschale | $\Sigma_p$ | alle $q$ mit $N(q)=p$ — **kein** Orbit |
| Einheitengruppe | $U_{\mathrm H}$ | $|U_{\mathrm H}|=24$ |
| Wahre Bahn | $O(q)=\{uqv\}$ | Doppelwirkung links und rechts |
| Zerlegung | $\Sigma_p=\sqcup O_i$ | disjunkte Vereinigung |

**Label:** $|U_{\mathrm H}|=24$ = **Theorem**; Zerlegung $\Sigma_p=\sqcup O_i$ = **Theorem** + **Experiment**
(`collatz_eabc_hurwitz_orbit_test.py::uh_orbit_partition`).

---

## 3. Glatt-EABC auf allen vier Koordinaten

Schreibe $q=(a,b,c,e)\in\Sigma_p$. Auf jeder Koordinate $|x|$ mit $x\neq 0$:

**Definition (glatt-EABC-Zerlegung).**
$$|x| = 2^{\alpha}\,3^{\beta}\, x',\qquad \gcd(x',6)=1,$$
$$\kappa(x')\in\{E,A,B,C\}\quad\text{via `eabc_from_lean.class_of`}.$$

Null-Bein: Label `0` (kein glatter Kern).

**Label:** glatt-EABC = **Definition** (analog Gauß §3, Eisenstein §3).

---

## 4. Volle $\Gamma$-Signatur — **nicht** orbit-invariant

**Definition (volle Signatur pro Bein).**
$$\Gamma_{\mathrm{leg}}(x) := (\alpha,\beta,\kappa(x'))\quad\text{oder }(0,0,0)\text{ bei }x=0.$$

**Definition.**
$$\Gamma(q) := \bigl(\Gamma_{\mathrm{leg}}(a),\,\Gamma_{\mathrm{leg}}(b),\,\Gamma_{\mathrm{leg}}(c),\,\Gamma_{\mathrm{leg}}(e)\bigr)
\in \bigl(\mathbb{N}_0\times\mathbb{N}_0\times\{E,A,B,C,0\}\bigr)^4.$$

Kompakt auch: $(\alpha_a,\beta_a,\kappa_a',\ldots,\alpha_e,\beta_e,\kappa_e')$.

**Experiment (verworfene Annahme).** $\Gamma(q)$ ist **keine** $U_{\mathrm H}$-Invariante:
unter $q\mapsto uqv$ permutieren die Einheiten die Achsen; $\Gamma(q)$ ändert sich typischerweise.
Im Gauß-Fall existiert ein **kanonischer** Repräsentant — in $\mathbb{H}_{\mathrm H}$ nicht.

**Label:** $\Gamma(q)$ = **Definition**; Nicht-Invarianz = **Experiment** (epistemische Korrektur).

---

## 5. Schalenmaß $M_p$, $\mu_p$

**Definition.**
$$M_p(\gamma) := \#\{q\in\Sigma_p : \Gamma(q)=\gamma\},$$
$$\mu_p(\gamma) := \frac{M_p(\gamma)}{|\Sigma_p|} = \frac{M_p(\gamma)}{\sum_\eta M_p(\eta)}.$$

Damit ist $\mu_p$ eine Wahrscheinlichkeitsverteilung auf dem diskreten $\Gamma$-Raum.
Das fundamentale EABC-Objekt in vier Dimensionen ist **$\mu_p$**, nicht ein einzelnes $\Gamma(q)$.

**Klassische Forschungskette (p→Σ):**
$$p \;\longmapsto\; \Sigma_p \;\longmapsto\; \mu_p \;\longmapsto\; I(\mu_p).$$

**Invertierte Kette (Σ→p), §11:** $n \mapsto \Sigma_n \mapsto \mu_n \mapsto I(\mu_n)$ ist fundamental;
Primzahlen **selektieren** Schalenparameter $n$, statt geometrisch zu *sein*.

**Label:** $M_p$, $\mu_p$ = **Definition**.

---

## 6. Orbitweise Maße $\mu^{(i)}$

Für jede Bahn $O_i\subset\Sigma_p$:

**Definition.**
$$M^{(i)}(\gamma) := \#\{q\in O_i : \Gamma(q)=\gamma\},\qquad
\mu^{(i)}(\gamma) := \frac{M^{(i)}(\gamma)}{|O_i|}.$$

Vergleich $\mu_p$ (Schale) vs. $\mu^{(i)}$ (einzelne Bahn) dient der Trennung von
Schalen- vs. Bahn-Statistik.

**Label:** $\mu^{(i)}$ = **Definition**; Orbit-Vergleich = **Experiment**.

---

## 7. Invarianten $H(p)$, $\chi(p)$, $K_p$ auf $\mu_p$

Einzelne $\Gamma(q)$ sind keine Invarianten. Kandidaten auf der **Verteilung** $\mu_p$:

| Symbol | Definition | Interpretation |
|--------|------------|----------------|
| **$\chi_p$** | $\displaystyle\sum_\gamma \chi(\gamma)\,\mu_p(\gamma)$ | mittlere Chiralität auf $\Sigma_p$ |
| **$H_p$** | $\displaystyle-\sum_\gamma \mu_p(\gamma)\log\mu_p(\gamma)$ | niedrig $\Rightarrow$ orientiert; hoch $\Rightarrow$ isotrop |
| **$K_p$** | $\mathrm{Cov}_{\mu_p}\bigl(\kappa(a'),\kappa(b'),\kappa(c'),\kappa(e')\bigr)$ | Bein-Kovarianzmatrix |

**Definition (Punkt-Chiralität).**
$$\chi(q) := \#\{\text{Legs}\in\{E,C\}\} - \#\{\text{Legs}\in\{A,B\}\}$$
(nur über nicht-Null-Legs).

**Morley-Parallel (getrennte Spur):** `collatz_morley_tm_numerik.py` — Sensoren $G_M$, $W_M$ auf
Dreiecken; **nur heuristische Analogie** zu $I(\mu_p)$, kein Beweisbaustein.

**Label:** $\chi_p$, $H_p$, $K_p$ = **Experiment**.

---

## 8. Hauptasymptotik: $\displaystyle\lim_{p\to\infty}\mu_p = \mu_\infty$?

**Vier Teilfragen (offen):**

1. **Existenz:** Existiert $\displaystyle\lim_{p\to\infty}\mu_p$ im Sinne schwacher Konvergenz auf dem $\Gamma$-Raum?
2. **Chiralität:** Bleibt $\chi_p$ von Null getrennt ($\chi_p\not\to 0$)?
3. **Entropie:** Konvergiert $H_p$ gegen ein universelles Niveau (isotrop vs. orientiert)?
4. **Restklassen:** Gibt es ein **universelles Spektrum** unabhängig von $p\bmod 12$?

**Zeta-Brücke (Forschungsrichtung):** Asymptotik von $M_p(\gamma)$ für große $p$ könnte mit
Koeffizienten der $\theta$-Serie / Epstein-Zeta für $\mathbb{H}_{\mathrm H}$ korrelieren —
**nicht bewiesen**.

**Label:** Limit-Verhalten = **Conjecture**; numerische Trends nach $p$-Größe = **Experiment**
(`H_trend`, `asymptotic_buckets` im JSON).

---

## 9. EABC-Maßhypothese (Conjecture)

> **Conjecture (Quaternionen-EABC-Maßhypothese).**
> Für die Hurwitz-Normschalen $\Sigma_p$ existiert eine asymptotische Verteilung
> $\mu_\infty$ auf dem glatt-EABC-$\Gamma$-Raum der vier Beine, sodass
> $$\mu_p \xrightarrow{p\to\infty} \mu_\infty$$
> schwach konvergiert und nontrivial strukturierte Invarianten $I(\mu_p)=(\chi_p,H_p,K_p,\ldots)$
> von der uniformen Null-Hypothese (Produkt der Bein-Marginalen, $\chi_p\approx 0$, maximales $H_p$)
> **messbar abweichen** — oder die Hypothese ist durch Shuffle-Null und Marginal-Tests **falsifiziert**.

**Kanalpaare:** AB-Kanal $(\kappa_a',\kappa_b')$, CE-Kanal $(\kappa_c',\kappa_e')$.
**Unabhängigkeitstest:** $\chi^2(\mu_p \,\|\, \prod_i \mathrm{marg}_i)$.
**$V_4$-Brücke:** `collatz_eabc_invarianzprogramm.md` — Paarflip-Struktur auf EABC-Klassen.

**Label:** Maßhypothese = **Conjecture**.

---

## 10. Forschungsprogramm und Falsifikation

**Zentrale Objekte (nicht $q$ einzeln):**
- $\Sigma_p$ — Normschale
- $\mu_p$ — Schalenverteilung
- $I(\mu_p)$ — Invariantenfunktional $(\chi_p, H_p, K_p, \chi^2_{\mathrm{indep}}, \ldots)$

**Forschungsprogramm:**
1. Verifiziere $r_4(p)=8(p+1)$ und $|\Sigma_p|$ für wachsende $p$.
2. Zerlege $\Sigma_p=\sqcup O_i$ unter $U_{\mathrm H}$; vergleiche $\mu_p$ vs. $\mu^{(i)}$.
3. Messe volle $\Gamma=(\alpha,\beta,\kappa)$ pro Bein; stratifiziere nach $(\alpha,\beta)$-Mustern.
4. Tracke $H_p$, $\chi_p$, $K_p$ und $\mu_p$-Marginalen in $p$-Größen-Buckets (klein/mittel/groß).
5. Vergleiche mit Gauß- und Eisenstein-Punkt-$\Gamma$ als 2D-Referenz.

**Falsifikationskriterien:**
- $\chi_p\approx 0$ ohne Restklassen-Struktur
- $H_p$ maximal (isotrop) ohne $p\bmod 12$-Bias
- $\mu_p\approx$ Produkt der vier Bein-Marginalen (Shuffle-Null nicht übertroffen)
- $\Rightarrow$ Hurwitz-Schale trägt **keine** zusätzliche EABC-Orientierung jenseits mod-$12$-Marginalen

**Null-Kontrolle (Experiment):** Zufällige Permutation der $\Gamma$-Labels auf $\Sigma_p$ (Shuffle-Test).

---

## 11. Perspektivwechsel Σ→p: Schalen fundamental, Primzahlen als Defekt-Singularitäten

### 11.1 Inversion der Forschungsrichtung

Bisher (§5): **p→Σ** — gegeben Primzahl $p$, studiere Schale $\Sigma_p$ und Maß $\mu_p$.

**Neue Lesart (Σ→p):** Die **Normschalen** $\Sigma_n$ ($n\ge 1$) sind das fundamentale geometrische Objekt in $\mathbb{H}_{\mathrm H}$.
Primzahlen sind **keine** Geometrie selbst, sondern **Selektoren** des Parameters $n$:
$$n \;\longmapsto\; \Sigma_n \;\longmapsto\; \mu_n \;\longmapsto\; I(\mu_n).$$

Die Frage „Was ist eine Primzahl?" wird ersetzt durch:
> **Welche Schalen $\Sigma_n$ tragen außergewöhnliche EABC-Invarianten?**

**Label:** Σ→p-Perspektive = **Heuristik** (Tao: Conjecture/Heuristik); kein Ersatz für §5-Theoreme.

### 11.2 Spektral-Analogie: Primzahlen wie $\mathrm{Spec}(\Delta)$, nicht wie die Geometrie

Analogie (explizit **nicht** bewiesen): Auf einer Riemannschen Mannigfaltigkeit $(M,g)$ ist die **Geometrie** $g$
fundamental; das **Laplace-Spektrum** $\mathrm{Spec}(\Delta)$ kodiert globale Information, ist aber **abgeleitet**.

Hier:
- **Geometrie:** Vereinigung aller Hurwitz-Normschalen $\displaystyle \mathbb{H}_{\mathrm H}^{\mathrm{shell}} := \bigcup_{n\ge 1} \Sigma_n$.
- **Spektrum:** Die Folge $(\mu_n)_{n\ge 1}$ und ihre Invarianten $I(\mu_n)$.
- **Primzahlen:** Erscheinen wie **eigenwertartige Singularitäten** — nicht als Grundobjekt der Geometrie,
  sondern als $n$, an denen $I(\mu_n)$ vom typischen Verhalten **zusammengesetzter** $n$ abweicht.

**Label:** Spektral-Analogie = **Heuristik**.

### 11.3 Invariantenfunktional und Defekt $D(n)$

**Definition (Schalen-Invarianten).** Für jedes $n$ mit $\Sigma_n\neq\emptyset$:
$$I(\mu_n) := \bigl(H_n,\,\chi_n,\,K_n,\,\ldots\bigr),$$
wobei $H_n=-\sum_\gamma \mu_n(\gamma)\log\mu_n(\gamma)$, $\chi_n=\sum_\gamma \chi(\gamma)\mu_n(\gamma)$,
$K_n=\mathrm{Cov}_{\mu_n}(\kappa_a',\kappa_b',\kappa_c',\kappa_e')$ wie in §7.

**Definition (Baseline).** $I_{\mathrm{avg}}(n)$ ist ein Referenzniveau, z.\,B.:
- rollierendes Mittel über Nachbarn $m\in[n-w,n+w]\setminus\{n\}$;
- Stratum-Mittel über $\omega(m)=\omega(n)$ (Anzahl verschiedener Primfaktoren);
- globales Mittel $\bar I$ über alle untersuchten $n$ (Proxy für $\mu_\infty$).

**Definition (Schalen-Defekt).**
$$D(n) := I(\mu_n) - I_{\mathrm{avg}}(n).$$

**Heuristik (Prim-Singularität):** Primzahlen sind die $n$, an denen $\|D(n)\|$ (z.\,B. euklidische Norm
der Komponenten $(D_H,D_\chi,D_K)$) **außergewöhnlich groß** ist — vergleichbar mit Ausreißern im Spektrum.

**Experiment:** `collatz_eabc_shell_defekt_test.py` vergleicht $|D(n)|$ für Prim- vs. zusammengesetztes $n$
auf kleinem $n$ (Geschwindigkeit); Ergebnis ist **explorativ**, kein Theorem.

**Label:** $D(n)$ = **Definition**; Prim-Defekt-Singularität = **Conjecture** / **Heuristik**.

### 11.4 Bernoulli-Korrektur: $\mu_n - \mu_\infty$

**Heuristik (Bernoulli / Endlich-Schalen-Korrektur).** Falls $\mu_n \xrightarrow{n\to\infty} \mu_\infty$ (§8),
dann ist $\mu_n - \mu_\infty$ eine **Diskretisierungskorrektur** — analog Endglied-Korrekturen in
Bernoulli-Summation (vgl. `collatz_eabc_normabstieg_hypothese.md` §5, Defektfeld).

Primzahlen wären dann die $n$, bei denen diese Korrektur **nicht** durch Mittelung über viele $n$ absorbiert wird.

**Label:** Bernoulli-Heuristik = **Heuristik** (nicht bewiesen).

### 11.5 Quaternionen vs. Gauß: erster Ort für Entropie und Maßtheorie

| Aspekt | Gauß $\mathbb{Z}[i]$ | Hurwitz $\mathbb{H}_{\mathrm H}$ |
|--------|----------------------|----------------------------------|
| Objekt | Punkt-Spaltung | **Schale** $\Sigma_n$ + Maß $\mu_n$ |
| Invariante | $\Gamma(z)$ pro Punkt | $I(\mu_n)=(H_n,\chi_n,K_n)$ |
| Primzahl-Rolle | split/inert (Theorem) | **Defekt-Singularität** (Conjecture) |
| Maßtheorie | marginal | **zentral** (§5–§7) |

In vier Dimensionen ist $\mu_n$ das natürliche EABC-Objekt; Entropie $H_n$ und Kovarianz $K_n$
existieren **erst hier** in voller Form — nicht in der 2D-Punkt-Statistik.

### 11.6 Hauptvermutung (Σ→p)

> **Conjecture (Σ→p-EABC-Schalen-Primhypothese).**
> Die Abbildung $n \mapsto \Sigma_n \mapsto \mu_n$ ist das **fundamentale** EABC-Objekt in $\mathbb{H}_{\mathrm H}$.
> Rational-prime $p$ sind genau die (oder dominieren unter) die $n$, bei denen die EABC-Geometrie von $\Sigma_n$
> — gemessen durch $I(\mu_n)$ und den Defekt $D(n)$ — **signifikant** vom typischen Verhalten
> zusammengesetzter Schalen abweicht.

**Falsifikation (Experiment):** Auf wachsenden $n$-Intervallen zeigt $|D(n)|$ **keinen** systematischen
Prim-Überhang gegenüber $\omega(n)$-stratifizierten Baselines → naive Singularitätslesart verworfen oder verfeinert.

**Label:** Σ→p-Hauptvermutung = **Conjecture**; numerischer Prim-vs.-Composite-Vergleich = **Experiment**.

---

## 12. EABC-Spektralgeometrische Hauptvermutung (kanonisch)

### 12.1 Perspektivinversion: $\mathcal H = \bigsqcup_n \Sigma_n$ fundamental

**Heuristik (geometrisch präzise).** Das eigentliche Objekt ist nicht die Primzahl $p$, sondern die
**disjunkte Vereinigung aller Normschalen**
$$\mathcal H^{\mathrm{shell}} := \bigsqcup_{n\ge 1} \Sigma_n,\qquad
\Sigma_n = \{q\in\mathbb{H}_{\mathrm H} : N(q)=n\}.$$
Die Zahl $n$ ist **kein Ursprung** der Schale, sondern ihr **Label** (Normniveau).
Die fundamentale Abbildung ist
$$\boxed{\; n \longmapsto \Sigma_n \longmapsto \mu_n \;}.$$
Primzahlen sind **nicht** die primären Objekte der Theorie.

**Label:** Perspektivinversion = **Heuristik** (präzisiert §11.1).

### 12.2 Spektral-Analogie: Eigenwerte vs. Geometrie

Auf einer Riemannschen Mannigfaltigkeit $(M,g)$ ist die **Metrik** $g$ fundamental;
$\mathrm{Spec}(\Delta)=\{\lambda_1,\lambda_2,\ldots\}$ ist **abgeleitet** — Projektion der Geometrie,
nicht die Geometrie selbst.

**Analogie (explizit nicht bewiesen):**
- **Geometrie:** $\mathcal H^{\mathrm{shell}}=\bigsqcup_n\Sigma_n$ und die Maße $\mu_n$.
- **Spektrum:** die Folge $(I(\mu_n))_{n\ge 1}$ für ein Invariantenfunktional $I$.
- **Primzahlen:** erscheinen wie **spektrale Anomalien** — nicht als Grundobjekt, sondern als $n$,
  an denen $I(\mu_n)$ vom typischen Verhalten **zusammengesetzter** Schalen abweicht.

**Label:** Spektral-Analogie = **Heuristik**.

### 12.3 Pipeline und Observable

**Definition (Invariantenfunktional).** Für jedes $n$ mit $\Sigma_n\neq\emptyset$ und kanonischem
EABC-Maß $\mu_n$ (§5):
$$I(\mu_n) := \bigl(H_n,\,\chi_n,\,K_n,\,\ldots\bigr),$$
mit $H_n$, $\chi_n$, $K_n$ wie in §7.

Die Folge $I(2), I(3), I(4), \ldots$ bildet das **Spektrum der EABC-Geometrie** der Normschalen.

**Label:** $I(\mu_n)$ = **Definition**.

### 12.4 Anomaliefunktion $D(n)$ und Referenz $I_{\mathrm{ref}}$

**Definition (Referenzmodell).** $I_{\mathrm{ref}}(n)$ ist ein Referenzniveau für die Observable
$I(\mu_n)$. Kandidaten (Experiment `collatz_eabc_shell_defekt_test.py`):

| Symbol | $I_{\mathrm{ref}}(n)$ | Interpretation |
|--------|----------------------|----------------|
| **rolling** | Mittel über Nachbarn $m\in[n-w,n+w]\setminus\{n\}$ | lokale Glättung |
| **cumulative** | $\displaystyle\frac{1}{n-1}\sum_{m<n} I(\mu_m)$ | Präfix-Mittel der Spektralfolge |
| **$\omega$-Stratum** | Mittel über $\{m : \omega(m)=\omega(n),\, m\neq n\}$ | arithmetische Schichtung |
| **$\tau$-Stratum** | Mittel über $\{m : \tau(m)=\tau(n),\, m\neq n\}$ | Teileranzahl-Baseline |
| **$\mu_\infty$-Proxy** | globales Mittel $\bar I$ über untersuchtes $n$-Intervall | endliches Grenzmaß |

**Definition (Anomalie / Schalen-Defekt).**
$$D(n) := I(\mu_n) - I_{\mathrm{ref}}(n).$$

**Epistemische Warnung:** Die Wahl von $I_{\mathrm{ref}}$ entscheidet, ob beobachtete Anomalien
**neue Geometrie** oder **repackagierte Arithmetik** sind (z.\,B. $\omega(p)=1$, $\tau(p)=2$ für alle
Primzahlen $p$). Das Experiment muss Baselines **ehrlich vergleichen**.

**Label:** $D(n)$, $I_{\mathrm{ref}}$-Varianten = **Definition** + **Experiment**.

### 12.5 Emergenz statt Voraussetzung

**Klassisch:** $\text{Primzahl} \Rightarrow \text{arithmetische Eigenschaft}$.

**Spektralgeometrische Lesart (Conjecture):**
$$\text{geometrische Eigenschaft von } \Sigma_n \;\Rightarrow\; \text{Primzahlstruktur emergiert}.$$

Primzahligkeit wird **nicht vorausgesetzt**, sondern soll aus der Struktur der Anomalien $D(n)$
**charakterisierbar** sein — falls die Vermutung zutrifft.

**Label:** Emergenz-Aussage = **Conjecture**.

### 12.6 Bernoulli-Brücke (Heuristik)

Bernoulli-Zahlen messen die Differenz zwischen diskreter Summe und kontinuierlichem Integral.
**Heuristik:** $\mu_n$ beschreibt die **diskrete** Schale; $\mu_\infty$ (§8–§9) die ideale
**kontinuierliche** Grenzgeometrie. Der Defekt $\mu_n - \mu_\infty$ ist ein Kandidat für
Bernoulli-artige Endlich-Korrekturen.

Explorativ im Experiment: $V_n := \|I(\mu_n) - I_\infty\|$ als Proxy; Korrelation mit $|D(n)|$
(vgl. `collatz_eabc_bernoulli_uebersetzung.md`, Branch `collatz/eabc-bernoulli-sensor`).

**Label:** Bernoulli-Brücke = **Heuristik**.

### 12.7 Boxed Hauptvermutung

> **Conjecture (EABC-Spektralgeometrische Hauptvermutung).**
> Sei $\Sigma_n = \{q\in\mathbb{H}_{\mathrm H} : N(q)=n\}$ die Hurwitz-Normschale der Norm $n$.
> Auf jeder Normschale sei ein kanonisches EABC-Maß $\mu_n$ definiert (§5).
> Zu jedem Maß existiere eine Familie von Invarianten $I(\mu_n)$.
> Die Folge $(I(\mu_n))_{n\ge 1}$ bildet das **Spektrum der EABC-Geometrie**.
>
> **Primzahlen sind nicht die fundamentalen Objekte der Theorie.**
> Vielmehr sind Primzahlen diejenigen Normniveaus $n$, auf denen die zugehörigen EABC-Invarianten
> **signifikant** vom typischen Verhalten der Normschalen zusammengesetzter Zahlen abweichen.
>
> Formal existiert eine Anomaliefunktion
> $$D(n) = I(\mu_n) - I_{\mathrm{ref}}(n),$$
> wobei $I_{\mathrm{ref}}$ ein geeignetes Referenzmodell beschreibt.
> Die Vermutung lautet, dass die Primzahlen durch die Struktur der Anomalien $D(n)$
> **charakterisiert werden können** — sofern $I_{\mathrm{ref}}$ geometrisch und nicht rein
> arithmetisch trivial gewählt ist.
>
> Die Frage „Was ist eine Primzahl?" wird ersetzt durch:
> **„Welche Normschalen tragen außergewöhnliche EABC-Invarianten?"**

**Falsifikation (Experiment):** Kein $I_{\mathrm{ref}}$ zeigt auf wachsenden $n$-Intervallen
stabile Prim-vs.-Composite-Separation in $|D(n)|$, nach Kontrolle für $\omega(n)$, $\tau(n)$ und
Schalengröße $|\Sigma_n|$.

**Label:** EABC-Spektralgeometrische Hauptvermutung = **Conjecture**;
`collatz_eabc_shell_defekt_test.py` = **Experiment**.

---

## 13. EABC-Spektralgeometrische Erzeugerhypothese ($\hat D(s)$)

### 13.1 Ontologie: Höhenfunktion statt Arithmetik zuerst

**Heuristik (Ontologie-Shift).** Die kanonische Grundmenge ist nicht die Menge der Primzahlen,
sondern die **disjunkte Vereinigung aller Normschalen**
$$\mathcal H^{\mathrm{shell}} := \bigsqcup_{n\ge 1} \Sigma_n.$$
Die **Höhenfunktion** (Norm) ist
$$N : \mathcal H^{\mathrm{shell}} \longrightarrow \mathbb{Z}_{\ge 1},\qquad q \longmapsto N(q).$$
Die arithmetische Zahl $n$ ist das **Label** des Normniveaus, nicht das ontologische Ursprungsobjekt.
Die fundamentale Abbildung bleibt
$$\boxed{\; n \longmapsto \Sigma_n \longmapsto \mu_n \;},$$
und erst daraus werden Invarianten $I(\mu_n)$ und Anomalien $D(n)$ abgeleitet (§12).

**Label:** Höhen-Ontologie $\mathcal H^{\mathrm{shell}}$, $N$ = **Heuristik** (präzisiert §12.1).

### 13.2 Hauptobjekt-Upgrade: Dirichlet-Reihe statt $D(n)$ allein

Tiefe Strukturen in der Zahlentheorie erscheinen historisch oft erst in **Erzeugerfunktionen**
(Euler-Produkte, $\zeta$, $L$-Funktionen). Entsprechend ist das entscheidende Testobjekt nicht
die punktweise Folge $D(n)$, sondern die **Dirichlet-Reihe**
$$\boxed{\;\hat D(s) := \sum_{n\ge 1} \frac{D(n)}{n^s}\;},$$
wobei $D(n)$ die skalare Anomalie aus §12.4 ist (Experiment: euklidische Norm
$\|I(\mu_n)-I_{\mathrm{ref}}(n)\|$; bevorzugte Referenzen: **rolling** oder **$\mu_\infty$-Proxy**).

**Definition.** Für reelles oder komplexes $s$ und endliches $N$:
$$\hat D_N(s) := \sum_{2\le n\le N} \frac{D(n)}{n^s}$$
ist die **Partialsumme** (explorativ; Konvergenz bei $s\le 1$ nicht vorausgesetzt).

**Label:** $\hat D(s)$, $\hat D_N(s)$ = **Definition** + **Experiment**
(`collatz_eabc_dirichlet_D.py`).

### 13.3 Bernoulli-Brücke (Conjecture / Heuristik)

Für die Riemann-Zeta-Funktion gilt bekanntlich für $m\ge 1$
$$\zeta(1-2m) = -\frac{B_{2m}}{2m},$$
wobei $B_{2m}$ die Bernoulli-Zahlen sind (z.\,B. $B_2=\tfrac16$, $B_4=-\tfrac{1}{30}$).

**Heuristik / offene Frage:** Existiert eine Skalierung oder Renormierung, sodass Werte von
$\hat D(s)$ bei negativen ganzzahligen $s$ **Bernoulli-artige** Endlich-Korrekturen tragen?
Explorativ fragt das Experiment nach $\hat D(-2)$, $\hat D(-4)$ und vergleicht mit $B_2$, $B_4$
— **ohne** Gleichheitsbehauptung (vgl. `collatz_eabc_bernoulli_uebersetzung.md`, Bernoulli-Sensor).

**Epistemische Warnung:** $\zeta(-2)=\zeta(-4)=0$ (triviale Nullstellen); die natürliche
Bernoulli-Korrespondenz liegt bei **ungeraden** negativen $s$ ($\zeta(-1)$, $\zeta(-3)$, …).
Die Wahl $s\in\{-2,-4,\ldots\}$ ist dennoch ein numerischer Stresstest für Wachstum und
Renormierungsbedarf der Partialsummen.

**Label:** Bernoulli-Analogie für $\hat D$ = **Conjecture** / **Heuristik** + **Experiment**.

### 13.4 Boxed Erzeugerhypothese

> **Conjecture (EABC-Spektralgeometrische Erzeugerhypothese).**
> Sei $\mathcal H^{\mathrm{shell}}=\bigsqcup_n\Sigma_n$ mit Höhenfunktion $N$ und kanonischen
> EABC-Maßen $\mu_n$ wie in §5–§12. Sei $D(n)=I(\mu_n)-I_{\mathrm{ref}}(n)$ die Anomaliefunktion
> mit geometrisch motivierter Referenz $I_{\mathrm{ref}}$ (rolling oder $\mu_\infty$-Proxy).
>
> Die **Dirichlet-Erzeugerfunktion**
> $$\hat D(s)=\sum_{n\ge 1}\frac{D(n)}{n^s}$$
> ist das **entscheidende Spektralobjekt** der Theorie: Primzahl- und Schalenstruktur sollen sich
> nicht primär in einzelnen $D(n)$, sondern in der **globalen analytischen Signatur** von $\hat D$
> (Pole, Residuen, Bernoulli-artige Werte bei negativen $s$, Vergleich zu $\zeta$ und $L$-Funktionen)
> manifestieren.
>
> **Falsifikation (Experiment):** Wenn $\hat D_N(s)$ für alle untersuchten $s$ und wachsendes $N$
> lediglich ein **Repackaging** bekannter arithmetischer Dirichlet-Reihen ist (z.\,B. Vielfaches
> von $\sum \tau(n)/n^s$, $\sum \omega(n)/n^s$ oder Partialsummen von $\zeta(s)$ ohne zusätzliche
> EABC-Geometrie), dann liefert die Erzeugerhypothese **keine neue** analytische Struktur.

**Label:** EABC-Spektralgeometrische Erzeugerhypothese = **Conjecture**;
`collatz_eabc_dirichlet_D.py` = **Experiment**.

### 13.5 Kritische Bewertungsleiter (epistemisch)

| Stufe | Objekt | Neuheit vs. Arithmetik |
|-------|--------|-------------------------|
| 1 | $p \mapsto \Sigma_p$ | **nicht** neu — klassische Hurwitz-Geometrie |
| 2 | $n \mapsto \Sigma_n \mapsto \mu_n$ | Perspektivwechsel (§12), aber zählbar |
| 3 | $\mu_n \mapsto I(\mu_n)$ | erste **Invariante** der EABC-Geometrie |
| 4 | $I(\mu_n) \mapsto D(n)$ | erste **Anomalie**; $I_{\mathrm{ref}}$-Wahl kritisch |
| 5 | $D(n) \mapsto \hat D(s)$ | **entscheidender Test** — analytische Globalisierung |

Die Leiter endet nicht bei $D(n)$: erst die Erzeugerfunktion vergleicht die EABC-Spektralgeometrie
mit dem Kanon analytischer Zahlentheorie ($\zeta$, Bernoulli, $L$-Funktionen).

**Label:** Bewertungsleiter = **Heuristik** (Methodologie).

---

## 14. Plattenübergang — Quaternion als Testbed vor $\mathbb{O}$

**Heuristik (Plattenfolge).** Die Normschalen $(\Sigma_n)_{n\ge 1}$ auf $\mathbb{H}_{\mathrm H}$ bilden eine
diskrete **Schichtfolge**; Primzahlen sind Normniveaus, an denen die produktinduzierte Zerlegung
$\Sigma_a\times\Sigma_b\to\Sigma_{ab}$ **unterbrochen** wird. Kanonische Formulierung:
`collatz_eabc_plattenuebergang.md`.

**Experiment (Quaternion, $n\le 200$, ohne Neuberechnung):**

| Kanal | Quelle | Befund |
|-------|--------|--------|
| $Z(n)$, $\Delta Z(n)$ | `collatz_eabc_Z_decomposition.json` | $Z_{\mathrm{EABC}}\approx Z_{\mathrm{fact}}$ ($r=0.88$); $Z(p)=0$; $\Delta Z$-Prim-Ratio 1.64 aber $\omega$-Repackaging |
| $D(n)=I(\mu_n)-I_{\mathrm{ref}}$ | `collatz_eabc_shell_defekt.json` | rolling: kein Prim-Überhang (Ratio 0.90); Top-10-$|D|$: 0 Prim |

**Fazit:** Quaternionen liefern das **implementierte Testbed** für die Plattenhypothese, bevor volles
$\Lambda_{\mathbb{O}}$ budgetiert wird. Bisher **kein** Signal, das über klassische Produktirreduzibilität
und $\omega(n)$ hinausgeht — die 8D-Frage bleibt offen (`collatz_eabc_oktonion_singularitaet.md` §3.5).

**Label:** Plattenübergang = **Conjecture**; Quaternion-Tests = **Experiment** (epistemisch hohl bzw. Prim-$D$ nicht gestützt).

---

## Vergleich Gauß / Eisenstein / Hurwitz

| Aspekt | Gauß $\mathbb{Z}[i]$ | Eisenstein $\mathbb{Z}[\omega]$ | Hurwitz $\mathbb{H}_{\mathrm H}$ |
|--------|----------------------|--------------------------------|----------------------------------|
| Normbeine | $2$ | $2$ | $4$ |
| Kanonische Wahl | ja | ja | **nein** |
| Geometrisches Objekt | Spaltung | Spaltung | **Schale** $\Sigma_p$ + **Orbits** $O_i$ |
| $\Gamma$-Raum | $\{E,A,B,C\}^2$ | $\{E,A,B,C\}^2$ | volle $(\alpha,\beta,\kappa)^4$ |
| EABC-Objekt | Punkt-$\Gamma$ | Punkt-$\Gamma$ | **$\mu_p$** |

---

## Epistemische Einordnung

| Aussage | Label |
|---------|-------|
| $r_4(p)=8(p+1)$ | **Theorem** |
| $|\Sigma_p|$, $U_{\mathrm H}$-Bahnen | **Theorem** + **Experiment** |
| $\Gamma(q)$ orbit-invariant | **verworfen** |
| $\lim_{p\to\infty}\mu_p=\mu_\infty$ | **Conjecture** (§9) |
| $H_p$-Bias nach $p\bmod 12$ | **Experiment** |
| $\chi_p$-Drift | **Experiment** |
| Σ→p-Perspektive, Spektral-Analogie | **Heuristik** |
| EABC-Spektralgeometrische Hauptvermutung (§12) | **Conjecture** |
| $D(n)=I(\mu_n)-I_{\mathrm{ref}}(n)$ | **Definition** |
| $I_{\mathrm{ref}}$-Varianten (rolling, $\omega$, $\tau$, $\mu_\infty$) | **Definition** + **Experiment** |
| Primzahlen als $|D(n)|$-Singularitäten / Emergenz | **Conjecture** |
| Bernoulli $V_n$-Korrelation mit $|D(n)|$ | **Heuristik** + **Experiment** |
| `collatz_eabc_shell_defekt_test.py` | **Experiment** |
| EABC-Spektralgeometrische Erzeugerhypothese (§13) | **Conjecture** |
| $\hat D(s)=\sum D(n)/n^s$ | **Definition** + **Experiment** |
| Bernoulli-Brücke $\hat D(-2m)$ vs.\ $B_{2m}$ | **Heuristik** + **Experiment** |
| `collatz_eabc_dirichlet_D.py` | **Experiment** |

---

*Kanonsiche Quaternionen-EABC-Maßhypothese: Objekte sind $\Sigma_n$, $\mu_n$, $I(\mu_n)$, $D(n)$, $\hat D(s)$ — nicht einzelne $q$; Primzahlen emergieren als mögliche Spektralanomalien (§12); die analytische Signatur ist $\hat D(s)$ (§13). Oktanionische Erweiterung: `collatz_eabc_oktonion_singularitaet.md`.*
