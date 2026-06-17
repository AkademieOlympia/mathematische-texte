# Oktonionische EABC-Singularitätshypothese ($\Sigma_n^{(8)}$, $\mu_n$, $D(n)$)

**Status:** Forschungshypothese + Forschungsprogramm (kein Experiment auf $\mu_n$ in 8D)  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Stub:** `collatz_eabc_oktonion_shell_stub.py` → `collatz_eabc_oktonion_shell_stub.json` (nur $r_8(n)$, $n\le 10$)  
**Produktbäume:** `collatz_eabc_product_tree_stub.py` → `collatz_eabc_product_tree_stub.json` (Catalan/Klammerung, H vs. O)  
**Tao-Labels:** Definition | Theorem | Conjecture | Heuristik | Experiment (Quaternion-Referenz)

**Querverweise:**
- `collatz_eabc_plattenuebergang.md` — **kanonsiche Plattenübergangs-Hypothese** ($\Sigma_n$-Folge, Produktirreduzibilität)
- `collatz_eabc_quaternion_mass_hypothese.md` §12–§14 (Spektralgeometrie, $\hat D(s)$, Platten-Testbed) — **Quaternion-Referenzdaten**
- `collatz_eabc_zerlegungsregimen.md` / `collatz_eabc_Z_decomposition_test.py` — $Z(n)$, $\Delta Z$, operative Regimen-Zählung
- `collatz_eabc_euklidische_hebung.md` §4, §8–§9 (Hurwitz-Kette $\mathbb{R}\subset\mathbb{C}\subset\mathbb{H}\subset\mathbb{O}$)
- `collatz_eabc_normabstieg_hypothese.md` (Norm-Defekt-Abstieg, Gauß–EABC-Brücke)
- `collatz_eabc_bernoulli_uebersetzung.md` §17 (Forschungsvision Defekt-Tetraeder, $\mathcal{D}_{\mathrm{krit}}$)
- `collatz_generalangriff_2026.md` (strategischer Pointer)
- `collatz_oktonionen_beweis.pdf` (algebraischer Hintergrund, nicht EABC-Maßtheorie)

---

## 1. Gitter $\Lambda_{\mathbb{O}}$ und Normschale $\Sigma_n^{(8)}$

**Definition (Oktanionische Normschale, $\mathbb{Z}^8$-Stub).** Auf der **8-dimensionalen** Divisionsalgebra $\mathbb{O}$
sei zunächst das **ganzzahlige Koordinatengitter**
$$\Lambda_{\mathbb{O}}^{(8)} := \mathbb{Z}^8 \subset \mathbb{O}$$
(mit einer festen Identifikation $\mathbb{O}\cong\mathbb{R}^8$ über die Standardbasis $e_1,\ldots,e_8$).
Die **Norm** ist
$$N(x) = x_1^2 + x_2^2 + \cdots + x_8^2.$$
Die **Normschale** des Niveaus $n\ge 1$ ist
$$\Sigma_n^{(8)} := \{x\in\Lambda_{\mathbb{O}}^{(8)} : N(x)=n\}.$$

**Theorem (Jacobi, Acht-Quadrate / $\theta_3^8$).** Die Anzahl
$$r_8(n) := \#\{(x_1,\ldots,x_8)\in\mathbb{Z}^8 : \sum_{i=1}^8 x_i^2 = n\}$$
ist der Koeffizient von $q^n$ in $\theta_3(q)^8$; sie wächst asymptotisch wie $n^{3}$ mit arithmetischen
Korrekturen über Teilerfunktionen (Eisensteinreihe $E_4$). Explizite geschlossene Formeln sind
tabelliert (Hardy–Wright, OEIS A000118); der Stub verifiziert $r_8(n)$ für $n\le 10$ per Zählung.

**Ehrliche Grenze (Hurwitz-Oktanionen).** Die **Hurwitz-Maximalordnung** $\mathbb{O}_{\mathrm H}$
(mit halb-ganzzahligen Koordinaten und Summenbedingung) ist das geometrisch korrekte Gitter —
aber **ohne Assoziativität** kein klassischer euklidischer Ring (`collatz_eabc_euklidische_hebung.md` §4).
Vollständige Hurwitz-Schalen $|\Sigma_n^{\mathrm{Hurwitz}}|$ und $U_{\mathbb{O}}$-Bahnen sind **nicht**
im Stub implementiert; siehe Forschungsprogramm §8.

**Label:** $\Sigma_n^{(8)}$, $r_8(n)$ = **Definition** / **Theorem**; Hurwitz-Verfeinerung = **offen**.

---

## 2. Diskrete 7-Sphäre $S^7$ — nicht die Quaternionen-3-Sphäre $S^3$

| Algebra | Dimension | Normschale als Sphäre | Einheitengruppe | EABC-Beine |
|---------|-----------|---------------------|-----------------|------------|
| $\mathbb{H}$ (Quaternionen) | $4$ | $S^3\subset\mathbb{R}^4$ (Hurwitz: $|\Sigma_p|=24(p+1)$ für ungerades $p$) | $|U_{\mathrm H}|=24$ | $4$ Koordinaten |
| $\mathbb{O}$ (Oktanionen) | $8$ | $S^7\subset\mathbb{R}^8$ | $|U_{\mathbb{O}}|=240$ (E8-Wurzelsystem) | $8$ Koordinaten |

**Heuristik (dimensionsbedingter Unterschied).** In $\mathbb{H}$ ist die Normschale eine **3-Sphäre** mit
bekannter Hopf-Faserung $S^3\to S^2$. In $\mathbb{O}$ lebt die Schale auf $S^7$; die klassische
**Hopf-Faserung** lautet $S^7\to S^4$ (Oktanionen-Projektion). Die **Faserstruktur auf Schalen**
ist hier **Heuristik**, kein etabliertes Theorem für EABC-Maße $\mu_n$.

**Label:** $S^7$ vs. $S^3$ = **Heuristik** (geometrische Analogie).

---

## 3. Frage: Was unterscheidet die Prim-Normschale?

**Klassische Lesart (verworfen als alleinige Erklärung):** Primzahl $p$ $\Rightarrow$ kleine Punktmenge
$|\Sigma_p|$ oder spezielle Gitterpunkte.

**EABC-Lesart (Organisation, nicht Punktzahl):** Das fundamentale Objekt ist nicht $|\Sigma_n|$,
sondern das **Schalenmaß**
$$M_n(\gamma) := \#\{x\in\Sigma_n^{(8)} : \Gamma(x)=\gamma\},\qquad
\mu_n(\gamma) := \frac{M_n(\gamma)}{|\Sigma_n^{(8)}|},$$
wobei $\Gamma(x)$ die glatt-EABC-Signatur auf **allen acht Koordinaten** ist (analog §3–§5 der
Quaternionen-Maßhypothese). Primzahlen sollen sich durch **Verteilungsstruktur** $\mu_p$ unterscheiden,
nicht durch bloße Kardinalität $r_8(p)$.

**Theorem (Multiplikativität).** $N(xy)=N(x)N(y)$ auf $\mathbb{O}$.

**Heuristik (Faktorisierung zusammengesetzter Schalen).** Für $n=ab$ mit $\gcd(a,b)=1$ (in geeigneten
Fällen) entsteht $\Sigma_n^{(8)}$ aus **Produktstruktur** der Norm:
$$\Sigma_a \star \Sigma_b \;\leadsto\; \text{Faltung von Maßen } \mu_a, \mu_b.$$
Prim-Schalen wären **irreduzibel** in diesem Bild — keine nichttriviale Faktorisierung $n=ab$.

**Label:** EABC-Organisation = **Heuristik**; $N(xy)=N(x)N(y)$ = **Theorem**.

---

## 3.5 Platten als Normschichtfolge ($\Sigma_n$-Index)

**Definition (Platte).** Die Folge $(\Sigma_n^{(8)})_{n\ge 1}$ auf $\Lambda_{\mathbb{O}}$ ist die **diskrete
Plattenstruktur** in $\mathbb{R}^8$; der Index $n$ ist das Normniveau (nicht „Primzahl als Punkt").
Siehe `collatz_eabc_plattenuebergang.md` für die vollständige **EABC-Plattenübergangs-Hypothese**.

**Conjecture (Plattenübergang).** $\Sigma_p^{(8)}$ ist eine **nicht zusammensetzbare Übergangsschicht**:
keine nichttriviale Produktzerlegung $\Sigma_a\times\Sigma_b\to\Sigma_p$ für $p=ab$, $a,b>1$.
Das ist die 8D-Lesart von „Prim unterbricht produktinduzierte Plattenzerlegung".

**Quaternion-Testbed ($n\le 200$):** In $\mathbb{H}_{\mathrm H}$ zeigen $Z(n)$, $\Delta Z(n)$ und $D(n)$
mit rolling-$I_{\mathrm{ref}}$ **kein** EABC-Signal über klassische Faktorisierung hinaus
(`collatz_eabc_Z_decomposition.json`, `collatz_eabc_shell_defekt.json`). Die oktonionische Hypothese
bleibt **offen**, bis $\mu_n$ in 8D existiert.

**Label:** Plattenfolge = **Definition**; Plattenübergangs-Hypothese = **Conjecture**;
Quaternion-Befund = **Experiment** (vgl. `collatz_eabc_plattenuebergang.md` §4–§6).

---

## 3.6 Produktbäume und Nicht-Assoziativität (8D-Novelty)

**Definition (Mehrstufige Produktstruktur).** Für $n=f_1\cdots f_k$, $k\ge 2$, und binären Produktbaum $T$
(Catalan-Zahl $C_{k-1}$) betrachte die induzierte Abbildung
$$\Phi_T : \Sigma_{f_1}\times\cdots\times\Sigma_{f_k} \longrightarrow \Sigma_n,\qquad (x_1,\ldots,x_k)\mapsto x_1\star_T\cdots\star_T x_k,$$
wobei $\star_T$ die Klammerung entlang $T$ bezeichnet.

**Theorem (Norm bleibt; Geometrie nicht).** $N(\Phi_T(\mathbf{x}))=N(x_1)\cdots N(x_k)=n$ für alle $T$
(Norm-Multiplikativität). Für $\mathbb{O}$ kann $\Phi_T\neq\Phi_{T'}$ bei gleicher Komposition und $T\neq T'$.

**Schlüssel-Kontrast (H vs. O).**

| Algebra | Assoziativität | Catalan-Bäume pro $n=f_1\cdots f_k$ | Effektive $Z^{\mathrm{tree}}$ |
|---------|----------------|-------------------------------------|-------------------------------|
| $\mathbb{H}$ | ja | $C_{k-1}$, kollabiert zu **1** | $\approx Z_{\mathrm{fact}}$ (binär); Klammerung **keine** neue Geometrie |
| $\mathbb{O}$ | nein | bis zu $C_{k-1}$ **verschiedene** Kanäle | potenziell $>Z_{\mathrm{fact}}$; **oktonion-spezifisch** |

**Boxed Frage (Rekonstruierbarkeit).**
> $$\boxed{\;\text{Welche } \Sigma_n^{(8)} \text{ sind aus kleineren Platten + Produktbäumen } (T,\Gamma) \text{ nicht EABC-rekonstruierbar — und tritt das erst jenseits von } \mathbb{H} \text{ auf?}\;}$$

**Stub:** `collatz_eabc_product_tree_stub.py` zählt für $\mathbb{H}$ Catalan-Summen und zeigt Assoziativitäts-Kollaps;
für $\mathbb{O}$ nur theoretische Catalan-Budgets ohne $\mu_n$-Test.

**Label:** $\Phi_T$, $Z^{\mathrm{tree}}$ = **Definition**; H-Kollaps = **Theorem**; O-Kanalvielfalt = **Heuristik** / **offen**.

---

## 4. Drei Konjekturen (Oktanionen) und Quaternion-Vergleich

Analog zu `collatz_eabc_quaternion_mass_hypothese.md` §8–§12 werden drei Richtungen unterschieden:

| # | Konjektur (Oktanionen) | Quaternion-Experiment ($n\le 200$) | Epistemischer Stand |
|---|------------------------|-------------------------------------|---------------------|
| **(1) Irreduzibilität** | Prim-Schalen $\Sigma_p^{(8)}$ ohne Faktorisierung $p=ab$ in der Maßstruktur | Trivial (arithmetisch); Maß-Faktorisierung **offen** | **Definition** / **Heuristik** |
| **(2) Maximale Isotropie** | $\mu_p \approx \mu_\infty$ (maximale Entropie $H_p$, minimale Orientierung) | $H_p$ bei Prim **nicht** maximal: `prime_mean_delta_H = +0.41` vs. Composite $-0.12$ (Abweichung von $\mu_\infty$-Proxy) — Prim **weniger** isotrop als Mittel | **(2) eher widerlegt** als Prim-Isotropie |
| **(3) Maximale Anomalie** | $|D(p)|$ groß (Spektral-Defekt $D(n)=I(\mu_n)-I_{\mathrm{ref}}(n)$) | Rolling-$I_{\mathrm{ref}}$: `mean_abs_D_prime / composite = 0.90`, `primes_larger_on_mean: false`; Top-10-Outlier: **0 Prim** (rolling) | **(3) widerlegt** für bevorzugte rolling-Baseline |

**Kanonische Bewertung (Quaternion, $n\le 200$):**

- **Konjektur (3) Anomalie $D(n)$:** Für die in §12.4 **bevorzugte** Referenz **rolling** zeigt
  `collatz_eabc_shell_defekt_test.py` **keinen** systematischen Prim-Überhang in $|D(n)|$
  (`ratio_mean_prime_over_composite ≈ 0.90$). Die naive **Prim-Defekt-Singularität** ist damit
  **falsifiziert** (explorativ, kleines $n$).
- **Konjektur (2) Isotropie:** Prim-Schalen weichen von $\mu_\infty$ in $H_n$ **stärker ab** als
  zusammengesetzte — entgegen „maximal isotrop bei Prim“. Eher **Orientierung/Struktur** als
  Uniformität.
- **Artefakt-Warnung:** `cumulative`- und `mu_infinity`-Referenzen können künstlichen Prim-Überhang
  erzeugen (`best_ratio ≈ 1.37$), weil $\omega(p)=1$, $\tau(p)=2$ arithmetisch fixiert sind — **kein**
  Beweis neuer Geometrie (`collatz_eabc_shell_defekt.json`, Feld `best_I_ref.epistemic_note`).

**Dirichlet-Erzeuger $\hat D(s)$ (Quaternion §13):** `collatz_eabc_dirichlet_D.py` liefert **keinen**
Bernoulli-Match: $\hat D_N(-2)/B_2 \sim 2\times 10^7$, Verdict „kein Bernoulli-Match behauptet“;
$\zeta$-Repackaging nicht stabil. Die Erzeugerhypothese bleibt **offen**, nicht gestützt.

**Label:** Konjekturen (1)–(3) = **Conjecture**; Quaternion-Vergleich = **Experiment** (Referenz).

---

## 5. Hopf-Faserung $S^7\to S^4$ und Schalen-Heuristik

**Heuristik (Faser auf der Schale).** Für $x\in\Sigma_n^{(8)}\subset S^7$ betrachte die Hopf-Projektion
$$\pi_{\mathrm{Hopf}} : S^7 \longrightarrow S^4$$
(induziert durch die Oktanionen-Multiplikation mit rein imaginären Einheiten).
Die **Fasern** sind $S^3$-Bündel über $S^4$. Auf einer festen Normschale könnte $\mu_n$ eine
**Faser-Statistik** tragen (Verteilung der $\pi_{\mathrm{Hopf}}$-Bilder in $S^4$) — analog zu
$U_{\mathrm H}$-Bahnen in $\mathbb{H}$, aber mit $|U_{\mathbb{O}}|=240$ statt $24$.

**Ehrliche Grenze:** Ohne kanonische Achsenwahl (§4 der Quaternionen-Maßhypothese: $\Gamma(q)$ **nicht**
orbit-invariant) ist eine $U_{\mathbb{O}}$-invariante $\Gamma$-Signatur **noch unklarer** als in $\mathbb{H}$.
Hopf-Heuristik ist **Forschungsrichtung**, kein implementierter Test.

**Label:** Hopf-Faser-Heuristik = **Heuristik**.

---

## 6. Boxed: Oktonionische EABC-Singularitätshypothese

> **Conjecture (Oktonionische EABC-Singularitätshypothese).**
> Sei $\mathcal{O}^{\mathrm{shell}} := \bigsqcup_{n\ge 1} \Sigma_n^{(8)}$ die Vereinigung aller
> Oktanionischen Normschalen (Gitter $\Lambda_{\mathbb{O}}$) mit Höhenfunktion $N$.
> Die fundamentale Abbildung ist
> $$\boxed{\; n \longmapsto \Sigma_n^{(8)} \longmapsto \mu_n \;},$$
> mit glatt-EABC-Maß $\mu_n$ auf dem $\Gamma$-Raum der **acht** Koordinatenbeine.
>
> **Primzahlen** sind die Normniveaus $n=p$, an denen **mindestens eine** der folgenden Singularitäten
> gegenüber zusammengesetztem $n$ auftritt:
> 1. **Irreduzibilität** der Schalen-Maßstruktur (keine Faktorisierung $\mu_n \approx \mu_a \star \mu_b$);
> 2. **Extremale Isotropie** $\mu_p \approx \mu_\infty$;
> 3. **Extremale Anomalie** $|D(p)|$ groß in $I(\mu_n)=(H_n,\chi_n,K_n,\ldots)$.
>
> **Quaternion-Referenz ($n\le 200$):** Richtung **(3) Anomalie** ist für rolling-$I_{\mathrm{ref}}$
> **nicht gestützt**; Richtung **(2) Isotropie** (maximal bei Prim) **ebenfalls nicht**.
> Die oktonionische Hypothese ist daher **bewusst offen** und verlangt erst $\mu_n$ in 8D — oder
> muss auf andere Invarianten (Hopf-Faser, $U_{\mathbb{O}}$-Orbits) verfeinert werden.

**Label:** Oktonionische EABC-Singularitätshypothese = **Conjecture** (mit falsifizierter Quaternion-Richtung).

---

## 7. Ehrliche Grenzen

| Grenze | Konsequenz |
|--------|------------|
| **Assoziativität verloren** | Kein klassischer $\gcd$, links/rechts-Ideale nicht äquivalent (`collatz_eabc_euklidische_hebung.md` §4) |
| **$\gcd$ / Ideale problematisch** | Faktorisierung $\Sigma_a\star\Sigma_b$ algebraisch heikel; Maß-Faltung nur **Heuristik** |
| **Volle Enumeration** | $r_8(n)$ wächst $\sim n^3$; $|\Sigma_n|$ für großes $n$ **nicht** voll enumerierbar |
| **$\Gamma$ nicht kanonisch** | Acht Beine, $240$ Einheiten — orbit-invariante Signatur **unbekannt** |
| **Kein $\mu_n$-Experiment** | Nur $r_8(n)$-Stub für $n\le 10$; keine $H_n$, $\chi_n$, $D(n)$ in 8D |

**Label:** Grenzen = **etabliertes Problem** / **epistemische Warnung**.

---

## 8. Forschungsprogramm

1. **Hurwitz-Gitter $\mathbb{O}_{\mathrm H}$:** Schalengröße und Einheitengruppe $|U_{\mathbb{O}}|=240$
   exakt zählen (Literatur: Coxeter, Conway–Smith); Stub auf $\mathbb{Z}^8$ erweitern.
2. **$r_8(n)$-Asymptotik:** Koeffizienten von $\theta_3^8$ / $E_4$ als Referenz für $|\Sigma_n|$-Wachstum
   ohne Punktenumeration (`collatz_eabc_oktonion_shell_stub.py`).
3. **Maß $\mu_n$ in 8D:** Glatt-EABC auf acht Beinen; $H_n$, $\chi_n$, $K_n$ — Rechenaufwand
   explizit budgetieren (erst $n\le 20$? Stichproben auf $\Sigma_n$?).
4. **$U_{\mathbb{O}}$-Bahnen:** Zerlegung $\Sigma_n = \sqcup O_i$ analog `collatz_eabc_hurwitz_orbit_test.py`.
5. **Hopf-Statistik:** Faser-Verteilung auf $\Sigma_n^{(8)}$ als alternatives $I(\mu_n)$.
6. **Vergleich $\mathbb{H}\to\mathbb{O}$:** Welche der drei Konjekturen überlebt den Dimensionswechsel?
   Quaternion-Daten sprechen gegen **(2)** und **(3)**; **(1)** Irreduzibilität bleibt offen.
7. **$\hat D^{(8)}(s)$:** Dirichlet-Reihe der 8D-Anomalie — erst nach existierendem $D(n)$ sinnvoll
   (vgl. §13 Quaternion; Bernoulli-Brücke dort **ohne** Match).

**Hurwitz-Kette:** `collatz_eabc_euklidische_hebung.md` §9, `collatz_eabc_normabstieg_hypothese.md` §1.

---

## 9. Vergleichstabelle Hurwitz-Kette

| Stufe | Algebra | Schale | Maß | Prim-Singularität (Experiment) |
|-------|---------|--------|-----|--------------------------------|
| 1D | $\mathbb{R}$ / $\mathbb{Z}$ | Punkt | — | Peano-Defekt (Conjecture) |
| 2D | $\mathbb{C}$ / $\mathbb{Z}[i]$ | Punkt-Spaltung | marginal | split/inert **Theorem** |
| 4D | $\mathbb{H}$ / $\mathbb{H}_{\mathrm H}$ | $\Sigma_n$, $S^3$ | $\mu_n$ | $|D(n)|$ Prim-Überhang **nicht** (rolling, $n\le 200$) |
| 8D | $\mathbb{O}$ / $\mathbb{O}_{\mathrm H}$ | $\Sigma_n^{(8)}$, $S^7$ | $\mu_n$ (**offen**) | Hypothese **offen**; Quaternion gegen (2),(3) |

---

## Epistemische Einordnung

| Aussage | Label |
|---------|-------|
| $r_8(n)$, $\theta_3^8$ | **Theorem** |
| $\Sigma_n^{(8)}$, $\mu_n$, $I(\mu_n)$, $D(n)$ | **Definition** |
| Oktonionische EABC-Singularitätshypothese (boxed) | **Conjecture** |
| Hopf $S^7\to S^4$ auf Schalen | **Heuristik** |
| Faktorisierung $\Sigma_a\star\Sigma_b$ | **Heuristik** |
| Quaternion: Prim-Anomalie $D(n)$ (rolling) | **Experiment**, **falsifiziert** (explorativ) |
| Quaternion: Prim-Isotropie $\mu_p\approx\mu_\infty$ | **Experiment**, **nicht gestützt** |
| $\hat D(s)$ Bernoulli-Match | **Experiment**, **kein Match** |
| `collatz_eabc_oktonion_shell_stub.py` | **Experiment** (nur $r_8$, kleines $n$) |

---

*Kanonsiche Notiz: In 8D ist die Frage nicht „wie viele Punkte auf $\Sigma_p$?", sondern „welche Organisation trägt $\mu_p$?" — analog §12 der Quaternionen-Maßhypothese. Quaternion-Daten ($n\le 200$) falsifizieren bisher die Richtungen maximale Isotropie und maximale Anomalie; die oktonionische Hypothese bleibt als Forschungsprogramm mit expliziten Grenzen.*
