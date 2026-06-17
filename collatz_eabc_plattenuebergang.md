# EABC-Plattenübergangs-Hypothese ($\Sigma_n$-Folge, Produktirreduzibilität)

**Status:** Forschungshypothese + Experiment (Quaternion-Testbed, Oktanion-Ziel)  
**Branch:** `collatz/eabc-euklidische-hebung` (PR #54)  
**Experimente:**
- `collatz_eabc_Z_decomposition_test.py` → `collatz_eabc_Z_decomposition.json` (Commit `e9013a2`)
- `collatz_eabc_shell_defekt_test.py` → `collatz_eabc_shell_defekt.json`
- `collatz_eabc_product_tree_stub.py` → `collatz_eabc_product_tree_stub.json` (Produktbäume, H vs. O)
- `collatz_eabc_oktonion_associator.py` → `collatz_eabc_oktonion_associator.json` ($\mathfrak{a}_E(n)$, EABC-Associator)

**Tao-Labels:** Definition | Theorem | Conjecture | Heuristik | Experiment

**Querverweise:**
- `collatz_eabc_zerlegungsregimen.md` — operative $Z(n)$-, $\Delta Z$-Definitionen und Regimen-Zählung
- `collatz_eabc_quaternion_mass_hypothese.md` §14 — Quaternionen als Testbed vor vollem $\mathbb{O}$
- `collatz_eabc_oktonion_singularitaet.md` §3.5 — 8D-Plattenprogramm
- `collatz_eabc_euklidische_hebung.md` §4, §8–§9 (Hurwitz-Kette $\mathbb{R}\subset\mathbb{C}\subset\mathbb{H}\subset\mathbb{O}$)
- `collatz_generalangriff_2026.md` — strategischer Pointer

---

## 1. Platten = Normschichten $\Sigma_n$ auf $\Lambda_{\mathbb{O}}$

**Definition (Platte / Normschicht).** Sei $\Lambda_{\mathbb{O}}$ ein geeignetes **ganzzahliges oktanionisches Gitter**
(mit fester Identifikation $\mathbb{O}\cong\mathbb{R}^8$; präzise Hurwitz-Variante $\mathbb{O}_{\mathrm H}$
siehe `collatz_eabc_oktonion_singularitaet.md` §1) und $N:\mathbb{O}\to\mathbb{Z}_{\ge 0}$ die **quadratische Norm**
$$N(x)=x_1^2+\cdots+x_8^2.$$

Für jedes $n\in\mathbb{N}$ definiere die **Normschale** (Platte)
$$\Sigma_n := \{x\in\Lambda_{\mathbb{O}} : N(x)=n\}.$$

Der Index $n$ ist der **Plattenindex** (Normniveau). Die Folge
$$(\Sigma_n)_{n\ge 1}$$
bildet eine **diskrete Plattenstruktur** in $\mathbb{R}^8$: Schichten konstanter Höhe $N=n$.

**Perspektivwechsel (Organisation, nicht Punktzahl).** Die Vermutung lautet **nicht** primär
„$p$ ist prim“, sondern:
$$\Sigma_p \text{ ist eine nicht zusammensetzbare Übergangsschicht}$$
innerhalb der Folge $\Sigma_1,\Sigma_2,\Sigma_3,\ldots$

**Label:** $\Lambda_{\mathbb{O}}$, $\Sigma_n$, Plattenindex $n$ = **Definition**.

---

## 2. Multiplikativität und Produktabbildungen

**Theorem (Norm-Multiplikativität).** Auf $\mathbb{O}$ (und damit auf jedem Teilschritt der Hurwitz-Kette)
gilt
$$N(xy)=N(x)\,N(y).$$

**Definition (Produktabbildung).** Für $n=ab$ mit $a,b\in\mathbb{N}$ existiert die natürliche Abbildung
$$\Phi_{a,b}:\;\Sigma_a\times\Sigma_b \longrightarrow \Sigma_n,\qquad (x,y)\mapsto xy.$$
Dadurch besitzt $\Sigma_n$ eine **produktinduzierte Zerlegungsstruktur** aus kleineren Platten $\Sigma_a$, $\Sigma_b$.

**Heuristik (EABC-Maß-Faltung).** Das fundamentale EABC-Objekt ist nicht $|\Sigma_n|$, sondern das
Schalenmaß $\mu_n$ auf dem $\Gamma$-Raum (vgl. `collatz_eabc_quaternion_mass_hypothese.md` §5).
Schematisch für zusammengesetzte Schalen:
$$\Sigma_{ab} \;\leadsto\; \mu_{ab} \approx \mu_a \otimes \mu_b \quad\text{(Vier- bzw. Acht-Bein-Tensorbild).}$$

**Theorem (Produktirreduzibilität = klassische Primzahl).** Für Primzahl $p$ existiert keine nichttriviale
Zerlegung $p=ab$ mit $a,b>1$. Daher besitzt $\Sigma_p$ **keine nichttriviale Produktzerlegung** aus
kleineren Normschalen — das ist die geometrische Lesart von
$$p\text{ prim}\;\iff\;\Sigma_p\text{ ist produkt-irreduzibel}.$$

**Label:** $N(xy)=N(x)N(y)$ = **Theorem**; $\mu_a\otimes\mu_b$-Schema = **Heuristik**;
produkt-irreduzibel $\Leftrightarrow$ prim = **Theorem** (klassische Arithmetik, **keine neue Geometrie** allein).

---

## 2.5 Produktbäume und Klammerung (Mehrstufige Zerlegung)

**Warnung (Nicht-Assoziativität).** Die Abbildung $\Phi_{a,b}:\Sigma_a\times\Sigma_b\to\Sigma_{ab}$ ist für die **Norm**
multiplikativ korrekt ($N(xy)=N(x)N(y)$ auf $\mathbb{O}$). Für **Oktanionen** gilt im Allgemeinen jedoch
$$(xy)z \neq x(yz)$$
(`collatz_eabc_euklidische_hebung.md` §4). Norm und Geometrie der Produktabbildung entkoppeln sich:
gleiche Faktoren, **verschiedene Klammerung** — potenziell verschiedene EABC-Kanäle.

**Definition (Produktbaum).** Sei $n=f_1\cdots f_k$ mit $k\ge 2$, $f_i\ge 2$. Ein **binärer Produktbaum**
$T$ auf den Blättern $(f_1,\ldots,f_k)$ ist eine vollständige Klammerung der Multiplikation.
Es gibt $C_{k-1}$ solcher Bäume (Catalan-Zahl).

**Definition (Erweiterte Zerlegungsmenge).**
$$\mathcal{Z}_n^{\mathrm{tree}} := \bigl\{(f_1,\ldots,f_k,\,T,\,\Gamma)\;:\; n=f_1\cdots f_k,\; k\ge 2,\; T\text{ binärer Produktbaum}\bigr\},$$
wobei $\Gamma$ die EABC-Signatur des induzierten Kanals bezeichnet (Schalenmaße $\mu_{f_i}$ entlang $T$).
Setze $Z^{\mathrm{tree}}(n):=\#\mathcal{Z}_n^{\mathrm{tree}}$ (operational erst in 8D vollständig).

**Theorem (Assoziative Stufen $\mathbb{R},\mathbb{C},\mathbb{H}$).** Auf $\mathbb{H}$ (und allen Teilschritten der
Hurwitz-Kette mit Assoziativität) sind alle Klammerungen einer festen Faktorkomposition **äquivalent**:
$$(x_1\cdots x_k)\text{ unabhängig von }T.$$
Pro Komposition $n=f_1\cdots f_k$ kollabiert die Catalan-Vielfalt $C_{k-1}$ auf **eine effektive Klasse**.
Für **binäre** Faktorisierungen ($k=2$) ist $C_1=1$ ohnehin — der Produktbaum-Zusatz ändert $Z(n)$ auf $\mathbb{H}$
**nicht** gegenüber $\mathcal{Z}_{\mathrm{fact}}$.

**Heuristik (Oktanionen).** Auf $\mathbb{O}$ kann dieselbe Komposition $n=abc$ **zwei** geometrisch verschiedene
Abbildungen liefern: $(\Sigma_a\times\Sigma_b)\times\Sigma_c \to \Sigma_n$ vs.
$\Sigma_a\times(\Sigma_b\times\Sigma_c)\to\Sigma_n$. Hier ist $Z^{\mathrm{tree}}(n)$ potenziell **größer** als
$Z_{\mathrm{fact}}(n)$ — die **oktonion-spezifische Novelty** der Plattenhypothese.

**Experiment** (`collatz_eabc_product_tree_stub.py`, $n\le 50$):
- **H:** $Z_{\mathrm{EABC}}\approx Z_{\mathrm{fact}}$ bleibt unverändert; Catalan-Summen pro $n$ sind $\ge 1$,
  aber assoziativ auf **1 effektive Klasse/Komposition** reduziert.
- **O:** nur theoretische Catalan-Zählung (keine $\mu_n$-Enumeration); ab $k\ge 3$ erste echte Klammerungsvielfalt.

**Label:** Produktbaum, $\mathcal{Z}_n^{\mathrm{tree}}$ = **Definition**; H-Kollaps = **Theorem**;
O-Klammerungs-Heuristik = **Heuristik**; Stub = **Experiment** (H) / **Theorie** (O).

---

## 2.6 EABC-Associator und $\mathfrak{a}_E(n)$ (oktonion-spezifisches Observable)

**Korrektur.** Die Produktabbildung $\Phi_T:\Sigma_{n_1}\times\cdots\times\Sigma_{n_k}\to\Sigma_n$ ist in $\mathbb{O}$
**nicht** wohldefiniert ohne Angabe des Baums $T$. Zählen allein $Z(n)=|\mathcal{Z}_n|$ reicht nicht —
es zählen **Äquivalenzklassen** $(T,\Gamma)$ mit messbarer Baum-Abhängigkeit.

**Definition (Assoziator und EABC-Version).**
$$[x,y,z]=(xy)z-x(yz),\qquad
\Gamma_{\mathrm{assoc}}(x,y,z)=\Gamma((xy)z)-\Gamma(x(yz))$$
auf glatt-gestrippten EABC-Koordinaten ($16$ Exponenten $(\alpha_i,\beta_i)_{i=1}^8$).

**Definition.** $\mathfrak{a}_E(n):=$ Mittel von $\|\Gamma_{\mathrm{assoc}}\|$ über repräsentative Stichproben
zu Faktorisierungen $n=abc$ ($a,b,c\ge 2$) und Klammerung $(xy)z$ vs. $x(yz)$.

**Boxed Frage.**
> $$\boxed{\;\text{Minimieren Prim-Platten } \Sigma_p \text{ den mittleren EABC-Assoziator } \mathfrak{a}_E \text{ — oder maximieren zusammengesetzte Platten ihn?}\;}$$

**Experiment** (`collatz_eabc_oktonion_associator.py`, $n\le 50$):
- Prim: $\mathfrak{a}_E(p)$ undefiniert (kein $abc$), Konvention $0$ → **minimiert trivial**.
- Zusammengesetzt: $\mathfrak{a}_E(n)>0$ wo Tripel-Faktorisierung existiert.
- Quaternion-Teilalgebra: Assoziator $=0$; generisches $(e_1,e_2,e_4)$: beide Normen $>0$.

**Grenze:** Sampling, $\mathbb{Z}^8$-Stub, ein Klammerungspaar — keine volle Catalan-Mittelung über $\mu_n$.

**Label:** $\mathfrak{a}_E$ = **Definition** + **Experiment**; Prim-Minimum = **Theorem** (trivial);
ob Zusammengesetzte **maximieren** = **offen** (Profil variiert mit Faktorisierung).

---

## 3. Boxed: EABC-Plattenübergangs-Hypothese

> **Conjecture (EABC-Plattenübergangs-Hypothese).**
> Sei $\Lambda_{\mathbb{O}}$ mit multiplikativer Norm $N$ und Plattenfolge $(\Sigma_n)_{n\ge 1}$.
> **Primzahlen sind Normniveaus, an denen die produktinduzierte Plattenzerlegung unterbrochen wird.**
>
> Die EABC-Struktur misst, **wie sich die Orientierung** der Normschalen (Verteilung $\mu_n$ auf
> $V_4=\{E,A,B,C\}$ über alle Koordinatenbeine) **vor und nach** einem solchen Übergang verändert.
>
> Damit wird „Primzahl" zunächst eine **geometrische Aussage** (Produktirreduzibilität der Schicht).
> **Neu** wird die Hypothese erst, wenn zusätzliche Invarianten $Z(n)$, $\Delta Z(n)$ oder $D(n)$
> ein Profil zeigen, das **nicht** durch klassische Faktorisierung ($\omega(n)$, $\tau(n)$) erklärt ist.

**Label:** EABC-Plattenübergangs-Hypothese = **Conjecture** (mit expliziter Trivialitätsgrenze §2).

---

## 4. Testgrößen: $Z(n)$, $\Delta Z(n)$ (Zerlegungsmodi)

**Definition.** Sei $\mathcal{Z}_n$ die Menge der **EABC-Zerlegungsmodi** von $\Sigma_n$.
Operational (Hurwitz-Quaternionen, $n\le 200$, siehe `collatz_eabc_zerlegungsregimen.md`):

| Symbol | Bedeutung |
|--------|-----------|
| $\mathcal{Z}_{\mathrm{fact}}(n)$ | Ungeordnete Faktorisierungen $n=ab$, $a,b>1$ (klassische Baseline) |
| $\mathcal{Z}_{\mathrm{EABC}}(n)$ | Verschiedene Kanäle $(\sigma_a,\sigma_b)$ über alle Faktorisierungen |
| $Z(n)$ | Primär $Z(n):=\#\mathcal{Z}_{\mathrm{EABC}}(n)$ |
| $\Delta Z(n)$ | $Z(n)-Z(n-1)$ |

**Experiment** (`collatz_eabc_Z_decomposition.json`, `e9013a2`, $n\le 200$):

| Kennzahl | Wert | Lesart |
|----------|------|--------|
| $|\Delta Z|$-Mittel Prim / zusammengesetzt | Ratio **1.641** | Scheinbarer Prim-Sprung |
| Pearson $Z$ vs. $\omega(n)$ | $r=0.815$ | **Starkes** $\omega$-Repackaging |
| Pearson $Z_{\mathrm{fact}}$ vs. $Z_{\mathrm{EABC}}$ | $r=0.878$; 79.7 % exakte Gleichheit | **Kritische Falsifikation** |
| $Z(p)$ für Prim $p$ | **0** (keine Faktorisierung $p=ab$) | Definitions-Artefakt |
| $\mathcal{Z}_{\mathrm{regime}}$: Mittel Prim / comp | **1.00** / **0.50** | Trivial: Prim $\Rightarrow$ Regime $1$ per Konstruktion |

**Ehrliche Epistemik:** $\Delta Z$ zeigt Sprünge an Primzahlen (oft $Z(p)=0$ nach $Z(p-1)>0$),
aber das **repackagiert** überwiegend $\omega(n)$ und klassische Faktorisierung — **kein neues**
EABC-Phänomen, solange $Z_{\mathrm{EABC}}\approx Z_{\mathrm{fact}}$.

**Label:** $Z$, $\Delta Z$ = **Definition**; Quaternion-Lauf = **Experiment**; Verdict = **epistemisch hohl**
(gegenüber klassischer Arithmetik).

---

## 5. Testgröße: $D(n)=I(\Sigma_n)-I_{\mathrm{ref}}(n)$ (Spektral-Defekt)

**Definition** (`collatz_eabc_quaternion_mass_hypothese.md` §12). Mit $I(\mu_n)=(H_n,\chi_n,K_n,\ldots)$ und
Referenz $I_{\mathrm{ref}}(n)$:
$$D(n)=I(\mu_n)-I_{\mathrm{ref}}(n).$$

**Experiment** (`collatz_eabc_shell_defekt.json`, $n\le 200$):

| $I_{\mathrm{ref}}$ | $|D|$-Mittel Prim / comp | Ratio | Bewertung |
|--------------------|--------------------------|-------|-----------|
| **rolling** (bevorzugt) | 1.135 / 1.260 | **0.90** | **Kein** Prim-Überhang |
| cumulative | 1.37+ | >1 | Artefakt: $\omega(p)=1$ fixiert |
| $\omega$-Stratum | — | hoch | Repackaging bekannter Arithmetik |

Top-10-$|D|$-Ausreißer (rolling): **0/10** Prim. Verdict: *„rolling: kein Prim-Überhang"*.

**Ehrliche Epistemik:** Die scharfe Frage „zeigt $D(p)$ ein anderes Profil als $D(n)$ für zusammengesetztes $n$?"
ist für die **geometrisch motivierte** rolling-Baseline **nicht gestützt** (explorativ, kleines $n$).

**Label:** $D(n)$ = **Definition**; rolling-Experiment = **Experiment**, Prim-Anomalie **nicht gestützt**.

---

## 6. Boxed: Scharfe EABC-Forschungsfragen

> **Frage 1 (Zerlegungssprünge).**
> $$\boxed{\;\text{Sind die Sprünge }\Delta Z(n)\text{ an Primzahlen anders als an typischen zusammengesetzten Zahlen — **über** }\omega(n)\text{ hinaus?}\;}$$
> **Quaternion ($n\le 200$):** Nein im Wesentlichen; $Z_{\mathrm{EABC}}\approx Z_{\mathrm{fact}}$, $Z(p)=0$ artefaktisch.

> **Frage 2 (Spektralprofil).**
> $$\boxed{\;D(p)\text{ zeigt ein anderes statistisches Profil als }D(n)\text{ für zusammengesetzte }n\text{ — bei geometrisch gewähltem }I_{\mathrm{ref}}\text{?}\;}$$
> **Quaternion (rolling):** Nein; Ratio $\approx 0.90$; keine Prim-Dominanz in Top-Ausreißern.

> **Frage 3 (echte Geometrie).**
> $$\boxed{\;\mathcal{Z}_n\text{ oder }\mu_n\text{ besitzen Struktur, die **nicht** vollständig durch klassische Primfaktorzerlegung bestimmt ist.}\;}$$
> **Offen in 8D** (`collatz_eabc_oktonion_singularitaet.md`); Quaternion-Kanal bisher **nicht** positiv.

> **Frage 4 (Produktbäume / EABC-Rekonstruierbarkeit).**
> $$\boxed{\;\text{Welche Schalen } \Sigma_n \text{ verlieren EABC-Rekonstruierbarkeit aus kleineren Platten, wenn Produktbäume }(T,\Gamma)\text{ mitgezählt werden — und nur in }\mathbb{O}\text{?}\;}$$
> **H ($n\le 50$):** Klammerung trivial (Assoziativität); $Z_{\mathrm{EABC}}\approx Z_{\mathrm{fact}}$ unverändert.
> **O:** offen; erste Klammerungsvielfalt ab $k\ge 3$-Faktor-Kompositionen (Catalan $\ge 2$).

> **Frage 5 (EABC-Associator / $\mathfrak{a}_E$).**
> $$\boxed{\;\mathfrak{a}_E(p)\text{ vs. }\mathfrak{a}_E(n)\text{ für zusammengesetztes }n\text{: minimieren Prim-Platten den Associator?}\;}$$
> **Experiment ($n\le 50$):** Ja **trivial** ($\mathfrak{a}_E(p)=0$ per Definition); zusammengesetzte $n$ zeigen
> positive $\mathfrak{a}_E$ — ob **Maximum** bei bestimmten Faktorprofilen liegt, **offen**.

---

## 7. Quaternion-Testbed vor vollem $\mathbb{O}$

Die Hurwitz-Quaternionen ($\Lambda=\mathbb{H}_{\mathrm H}$, vier Koordinatenbeine) sind das **implementierte**
Testbed für Plattenübergänge, bevor $\mu_n$ in 8D budgetiert wird:
- gleiche Definition $\Sigma_n$, $\mu_n$, $Z(n)$, $D(n)$;
- volle Enumeration $|\Sigma_n|$ für $n\le 200$ machbar;
- Ergebnisse in §4–§5 gelten **explorativ** für die Plattenhypothese in 4D.

**Label:** Quaternion-Testbed = **Experiment**; 8D-Übertragung = **offen**.

---

## Epistemische Einordnung

| Aussage | Label |
|---------|-------|
| $\Sigma_n$ als Platte, Index $n$ | **Definition** |
| $N(xy)=N(x)N(y)$ | **Theorem** |
| $\Phi_{a,b}:\Sigma_a\times\Sigma_b\to\Sigma_n$ | **Definition** |
| $p$ prim $\Leftrightarrow$ $\Sigma_p$ produkt-irreduzibel | **Theorem** (klassisch) |
| EABC-Plattenübergangs-Hypothese (boxed §3) | **Conjecture** |
| $\mu_a\otimes\mu_b$-Faltung | **Heuristik** |
| $Z(n)$, $\Delta Z(n)$ | **Definition** + **Experiment** |
| $\Delta Z$-Prim-Sprünge ohne $\omega$-Neuheit | **Experiment**, **falsifiziert** (Quaternion) |
| $D(n)$, rolling-$I_{\mathrm{ref}}$ | **Definition** + **Experiment** |
| Prim-Profil in $|D(n)|$ (rolling) | **Experiment**, **nicht gestützt** |
| Produktbäume $T$, $Z^{\mathrm{tree}}(n)$ | **Definition** |
| H: Klammerungs-Kollaps (Assoziativität) | **Theorem** |
| O: Klammerung als geometrische Novelty | **Heuristik** (offen) |
| `collatz_eabc_product_tree_stub.py` | **Experiment** (H) / **Theorie** (O) |
| `collatz_eabc_oktonion_associator.py` | **Experiment** ($\mathfrak{a}_E$, $n\le 50$) |
| $\mathfrak{a}_E(n)$, EABC-Associator | **Definition** + **Experiment** |
| 8D-Platten $\Sigma_n^{(8)}$ | **Forschungsprogramm** |

---

*Kanonsiche Plattenhypothese: $(\Sigma_n)$ ist die diskrete Schichtfolge; Primzahlen unterbrechen die produktinduzierte Zerlegung — geometrisch trivial, EABC-neu erst wenn $Z$ oder $D$ über klassische Arithmetik hinausgehen. Quaternion-Experimente ($n\le 200$) zeigen bisher kein solches Signal.*
