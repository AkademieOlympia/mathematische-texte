# Oktonionische EABC-Assoziator-Spektralhypothese ($M_n(t)$, $M_n^E(t)$, $\hat D_E(s)$)

**Status:** Forschungshypothese + Sampling-Experiment (kein volles $\Sigma_n$-Zählen)  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Experiment:** `collatz_eabc_oktonion_spectrum.py` → `collatz_eabc_oktonion_spectrum.json`  
**Vorgänger:** `collatz_eabc_oktonion_associator.py` ($\mathfrak{a}_E(n)$, Mittelwert-Stub)  
**Tao-Labels:** Definition | Theorem | Conjecture | Heuristik | Experiment

**Querverweise:**
- `collatz_eabc_holonomie.md` — **kanonisch:** projektive Holonomie $\mathcal H_E=d_E\circ\Gamma$; $V_4$-„Assoziator“ deprecated
- `collatz_eabc_kommutator_assoziator.md` — Kommutator (Chiralität $\chi$, $\mathbb{H}$) vs. Assoziator ($M_n^E$, $\mathbb{O}$)
- `collatz_eabc_oktonion_singularitaet.md` §3.7–§3.8 — Assoziator, Spektrum
- `collatz_eabc_plattenuebergang.md` §2.6 — $\mathfrak{a}_E(n)$, Platten-Klammerung
- `collatz_eabc_oktonion_associator.py` — algebraischer und EABC-Assoziator
- `collatz_eabc_oktonion_spectrum.py` — Histogramme $M_n(t)$, $M_n^E(t)$, Prim-vs.-Composite-KL
- `collatz_eabc_dirichlet_D.py` — Quaternion-Analogon $\hat D(s)$ für $D(n)=I(\mu_n)-I_{\mathrm{ref}}$
- `collatz_eabc_quaternion_mass_hypothese.md` §12–§13 — Spektralgeometrie, Dirichlet-Erzeuger
- `collatz_eabc_plattenuebergang.md` — kanonsiche $\Sigma_n$-Plattenfolge
- `collatz_generalangriff_2026.md` — strategischer Pointer

**Epistemischer Schlüssel:** Nichtassoziativität allein $\not\Rightarrow$ Primzahlinformation.
Erst die **Verknüpfung** Nichtassoziativität + **Normniveau** $N(x)N(y)N(z)=n$ + **EABC-Projektion**
$\Gamma_E$ liefert ein schalenindiziertes Observable.

---

## 1. Gitter $\Lambda_{\mathbb{O}}$ und Platten $\Sigma_n$ in $\mathbb{R}^8$

**Definition (Oktanionisches Gitter, $\mathbb{Z}^8$-Stub).** Mit fester Identifikation $\mathbb{O}\cong\mathbb{R}^8$
über die Standardbasis $e_1,\ldots,e_8$ sei
$$\Lambda_{\mathbb{O}}^{(8)} := \mathbb{Z}^8 \subset \mathbb{O}.$$
Die **quadratische Norm** ist
$$N(x) = x_1^2 + \cdots + x_8^2.$$

**Definition (Platte / Normschale).** Für $n\ge 1$
$$\Sigma_n := \{x\in\Lambda_{\mathbb{O}}^{(8)} : N(x)=n\}.$$
Die Folge $(\Sigma_n)_{n\ge 1}$ ist die **diskrete Plattenstruktur** in $\mathbb{R}^8$ (Index $n$ = Normniveau).

**Theorem (Norm-Multiplikativität).** $N(xy)=N(x)N(y)$ auf $\mathbb{O}$.

**Ehrliche Grenze:** Hurwitz-Maximalordnung $\mathbb{O}_{\mathrm H}$ mit $|U_{\mathbb{O}}|=240$ ist das
geometrisch korrekte Gitter — Experimente nutzen vorerst $\mathbb{Z}^8$ (`collatz_eabc_oktonion_shell_stub.py`).

**Label:** $\Lambda_{\mathbb{O}}$, $\Sigma_n$ = **Definition**; $N(xy)=N(x)N(y)$ = **Theorem**.

---

## 2. Assoziator $[x,y,z]$ — Algebraeigenschaft, arithmetische Relevanz nur über Normniveaus

**Definition (algebraischer Assoziator).**
$$[x,y,z] := (xy)z - x(yz).$$

**Theorem.** $[x,y,z]=0$ für alle $x,y,z\in\mathbb{R},\mathbb{C},\mathbb{H}$; in $\mathbb{O}$ typischerweise
**nichtnull** (Nichtassoziativität).

**Warnung (Universalität).** Der Assoziator ist eine **Eigenschaft der Algebra $\mathbb{O}$**, nicht einer
einzelnen Normschale $\Sigma_n$. Für beliebiges Tripel $(e_1,e_2,e_4)$ gilt $\|[e_1,e_2,e_4]\|>0$
unabhängig vom Schalenindex.

**Definition (normgebundenes Tripelmenge).** Für $n\in\mathbb{N}$
$$\mathfrak{A}_n := \bigl\{(x,y,z)\in\Lambda_{\mathbb{O}}^3 : N(x)N(y)N(z)=n\bigr\}.$$
Arithmetische Relevanz entsteht erst durch die **Produktnorm-Bedingung** $N(x)N(y)N(z)=n$,
die Tripel an das Plattenniveau $n$ bindet.

**Label:** $[x,y,z]$ = **Definition**; Assoziativität in $\mathbb{H}$ = **Theorem**;
Schalenbindung über $\mathfrak{A}_n$ = **Definition**; bloße Nichtassoziativität = **keine Prim-Information**.

---

## 3. Algebraisches Assoziator-Spektrum $M_n(t)$

**Definition (Assoziator-Norm-Stufe).**
$$\alpha(x,y,z) := N\bigl([x,y,z]\bigr) = \sum_{i=1}^8 [x,y,z]_i^2 \in \mathbb{Z}_{\ge 0}.$$

**Definition (Spektrum $M_n$).**
$$M_n(t) := \#\bigl\{(x,y,z)\in\mathfrak{A}_n : \alpha(x,y,z)=t\bigr\},\qquad t\in\mathbb{Z}_{\ge 0}.$$

**Heuristik (partielle Zeta / Erzeuger).** Für $m_n(t):=M_n(t)/|\mathfrak{A}_n|$ (normalisiert) betrachte
$$S_n(s) := \sum_{t\ge 1} \frac{m_n(t)}{t^s}$$
als **arithmetisches Spektralobjekt** der Schale $n$ — **Conjecture** (noch ohne analytischen Beweis).

**Label:** $\alpha$, $M_n(t)$ = **Definition**; $S_n(s)$ = **Conjecture** (Erzeugerform).

---

## 4. EABC-Projektion $\Gamma_E$ und Distanz $d_E$ auf dem $\Gamma$-Raum

**Definition (glatt-EABC-Projektion).** Für $x\in\Lambda_{\mathbb{O}}^{(8)}$ sei
$$\Gamma_E(x) := \bigl(\alpha_1,\beta_1,\ldots,\alpha_8,\beta_8\bigr)\in\mathbb{Z}^{16},$$
wobei $(\alpha_i,\beta_i)$ die glatt-gestrippte EABC-Zerlegung der $i$-ten Koordinate ist
($|c_i|=2^{\alpha_i}3^{\beta_i}c_i'$, $\gcd(c_i',6)=1$; $(0,0)$ bei $c_i=0$; vgl. `kappa_glatt`).

**Definition (projektive EABC-Holonomie / $\Gamma_E$-Differenz).**
$$\Delta\Gamma_E(x,y,z) := \Gamma_E\bigl((xy)z\bigr) - \Gamma_E\bigl(x(yz)\bigr)\in\mathbb{Z}^{16},\qquad
\mathcal H_E\bigl((xy)z,x(yz)\bigr) = d_E\bigl(\Gamma_E((xy)z),\,\Gamma_E(x(yz))\bigr).$$
*(Auf $V_4$ ist der analoge Klammertest trivial; siehe `collatz_eabc_holonomie.md`.)*

**Definition (EABC-Holonomie-Norm-Stufe).**
$$\alpha_E(x,y,z) := N\bigl(\Delta\Gamma_E(x,y,z)\bigr) = \sum_{j=1}^{16} \Delta\Gamma_{E,j}^2.$$

**Definition (Distanz $d_E$ auf $\Gamma$-Raum).** Für $\gamma,\eta\in\mathbb{Z}^{16}$
$$d_E(\gamma,\eta) := \|\gamma-\eta\|_2 = \sqrt{\sum_{j=1}^{16}(\gamma_j-\eta_j)^2}.$$
Dann $\alpha_E(x,y,z)=d_E\bigl(\Gamma_E((xy)z),\,\Gamma_E(x(yz))\bigr)^2$.

**Label:** $\Gamma_E$, $\alpha_E$, $d_E$ = **Definition**.

---

## 5. EABC-Assoziator-Spektrum $M_n^E(t)$

**Definition.**
$$M_n^E(t) := \#\bigl\{(x,y,z)\in\mathfrak{A}_n : \alpha_E(x,y,z)=t\bigr\},\qquad t\in\mathbb{Z}_{\ge 0}.$$

**Definition (Invariantenfunktional auf Spektrum).** Mit $p_n^E(t):=M_n^E(t)/\sum_s M_n^E(s)$ (falls Support nichtleer):
$$I(M_n^E) := \bigl(H_n^E,\,|{\mathrm{supp}}_n^E|,\,\mathbb{E}_n^E[\alpha_E],\,\mathrm{Var}_n^E[\alpha_E]\bigr),$$
wobei $H_n^E:=-\sum_t p_n^E(t)\log p_n^E(t)$ die Entropie des EABC-Spektrums ist.

**Experiment** (`collatz_eabc_oktonion_spectrum.py`): Stichproben-Tripel aus $\mathfrak{A}_n$ auf $\mathbb{Z}^8$,
Histogramme für $M_n(t)$ und $M_n^E(t)$, Vergleich Prim vs. zusammengesetzt (KL-Distanz zum Composite-Referenzprofil).

**Label:** $M_n^E$, $I(M_n^E)$ = **Definition**; Sampling-Lauf = **Experiment**.

---

## 6. Boxed: Oktonionische EABC-Assoziator-Spektralhypothese

> **Conjecture (Oktonionische EABC-Assoziator-Spektralhypothese).**
> Sei $n\ge 2$ ein Normniveau (Plattenindex). Prim-Normen $n=p$ besitzen ein **charakteristisches**
> EABC-Assoziator-Spektrum $M_p^E(t)$, das sich vom **typischen** Spektrum zusammengesetzter $n$
> unterscheidet:
> $$\boxed{\;M_p^E \not\approx M_n^E\;\text{für typisches zusammengesetztes }n\;}$$
> im Sinne von $I(M_p^E)$, KL-Profil oder Support-Struktur — **nicht** im Sinne bloßer Maximierung
> von $\mathfrak{a}_E$ (Mittelwert kann bei Prim **kleiner** sein: Stabilität statt Chaos).

**Theorem (Trivialitätsgrenze).** Ohne EABC-Projektion ist $M_n(t)$ rein algebraisch und unterscheidet
Prim von Composite nur über **Faktorisierbarkeit** von $n=N(x)N(y)N(z)$ — kein neues Phänomen.

**Theorem (Definitions-Artefakt $\mathfrak{a}_E$).** Für Prim $p$ ohne nichttriviale Tripel-Zerlegung
$a,b,c\ge 2$ ist $\mathfrak{a}_E(p)=0$ per Konvention — **minimiert trivial**, nicht „charakteristisch".

**Label:** Spektralhypothese = **Conjecture**; Trivialitätsgrenzen = **Theorem** / **epistemische Warnung**.

---

## 7. Dirichlet-Erzeuger $\hat D_E(s)$ aus $D(n)=I(M_n^E)$

**Definition (Spektral-Anomalie).** Mit Referenz $I_{\mathrm{ref}}^E(n)$ auf der Composite-Schicht
(z.\,B. rolling-Mittel von $I(M_m^E)$ über $m<n$):
$$D_E^{\mathrm{spec}}(n) := I(M_n^E) - I_{\mathrm{ref}}^E(n).$$

Skalare Komponente für Dirichlet-Test (Experiment):
$$D(n) := \|D_E^{\mathrm{spec}}(n)\|_2 \quad\text{oder}\quad D(n):=H_n^E - H_{\mathrm{ref}}^E(n).$$

**Definition (Dirichlet-Erzeuger).**
$$\boxed{\;\hat D_E(s) := \sum_{n\ge 2} \frac{D(n)}{n^s}\;},$$
analog zu $\hat D(s)$ in `collatz_eabc_quaternion_mass_hypothese.md` §13 und `collatz_eabc_dirichlet_D.py`.

**Heuristik (analytischer Test).** Globale Signatur von $\hat D_E(s)$ (Bernoulli-Koeffizienten, $\zeta$-Vergleich)
ist der **entscheidende** Test — nicht einzelne $M_n^E(t)$-Stützpunkte. Quaternion-Referenz:
`collatz_eabc_dirichlet_D.py` zeigt **keinen** Bernoulli-Match für $\hat D(s)$ aus $\mu_n$-Defekt.

**Label:** $D(n)$, $\hat D_E(s)$ = **Definition**; Bernoulli-Brücke = **Experiment** (offen in 8D).

---

## 8. Falsifikationskriterien

| Kriterium | Verdict | Bedeutung |
|-----------|---------|-----------|
| **Triviale Symmetrien** | Falsifikation | $M_n^E$ hängt nur von $\omega(n)$, $\tau(n)$ oder Faktorisierbarkeit ab |
| **$\mathfrak{a}_E$-Mittel allein** | Falsifikation | Prim minimiert per Definition ($\mathfrak{a}_E(p)=0$); kein Spektrum |
| **Bloßer Assoziator in $\mathbb{O}$** | Falsifikation | $[x,y,z]\neq 0$ ohne $\mathfrak{A}_n$-Bindung — universell, nicht prim-spezifisch |
| **KL(Prim $\|\|$ Composite-Ref) $\approx 0$** für alle $p$ | Falsifikation (explorativ) | Kein unterscheidbares EABC-Spektralprofil |
| **Reproduzierbare Prim-Anomalien** in $I(M_p^E)$ | **Substanz** | KL $>0$, stabiler Support-/Entropie-Unterschied über Seeds und $n$-Fenster |

**Experiment-Verdict** (siehe `collatz_eabc_oktonion_spectrum.json`, Feld `falsification`):
wird beim Lauf von `collatz_eabc_oktonion_spectrum.py` ausgegeben.

**Label:** Falsifikationskriterien = **Definition** + **Experiment**.

---

## Epistemische Einordnung

| Aussage | Label |
|---------|-------|
| $\Lambda_{\mathbb{O}}$, $\Sigma_n$, $\mathfrak{A}_n$ | **Definition** |
| $[x,y,z]$, $N(xy)=N(x)N(y)$ | **Theorem** |
| $M_n(t)$, $M_n^E(t)$, $\alpha_E$, $d_E$ | **Definition** |
| $S_n(s)$, $\hat D_E(s)$ | **Conjecture** / **Definition** |
| Oktonionische EABC-Assoziator-Spektralhypothese (boxed) | **Conjecture** |
| Nichtassoziativität $\not\Rightarrow$ Prim | **Theorem** / **epistemische Warnung** |
| `collatz_eabc_oktonion_spectrum.py` | **Experiment** (Sampling, $\mathbb{Z}^8$) |

---

*Kanonsiche Notiz: Die Spektralhypothese verlangt die Dreifach-Verknüpfung Nichtassoziativität + Normniveau $n$ + EABC-$\Gamma_E$. Mittelwert $\mathfrak{a}_E$ allein ist zu schwach; das Histogramm $M_n^E$ und $I(M_n^E)$ sind die operative Testgröße — mit expliziten Trivialitäts- und Falsifikationskriterien (§8). Quaternion-Dirichlet-Referenz: `collatz_eabc_dirichlet_D.py`.*
