# EABC-Zirkulation und Spektralgeometrie (kanonische reine Mathematik)

**Status:** Definition + Vermutung + Hypothese + Experiment  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Vermutung | Hypothese | Experiment | Analogie

**Epistemische Abgrenzung:** Dieses Dokument formuliert das EABC-Programm **ohne Physik-Metaphern**. Die fundamentale Verschiebung ist **Korrelation → Zirkulation** (diskrete Differentialgeometrie und algebraische Topologie auf einem endlichen Graphen), nicht Bell→Sagnac-Physik. Bell-Korrelationen und Sagnac-Analogien bleiben als **sekundäre** didaktische Einträge dokumentiert (`collatz_eabc_bell_holonomie.md`, `collatz_eabc_sagnac.md`).

**Querverweise:**
- `collatz_eabc_fehlerterm_hypothese.md` — **kanonische Endform:** $N_\pm$, $D_E$, Hauptvermutung, Fehlerterm-Hypothese
- `collatz_eabc_sagnac.md` — **didaktischer Einstieg** (gegenläufige Zyklen; verweist hierher als mathematischer Kern)
- `collatz_eabc_transport.md` — $G_E$, Transport $T_n$, Übergangsmatrix
- `collatz_eabc_zyklus_holonomie.md` — Hierarchie Klasse→Kante→Pfad→Zyklus→Holonomie
- `collatz_eabc_sagnac_circulation.py` — $C_E(X)$, $\omega(e)$, diskrete 1-Form $\alpha$
- `collatz_eabc_graph_laplacian.py` — $\mathrm{Spec}(L_E)$, Spektrallücke, Normschalen
- `collatz_eabc_holonomie_fehlerterm.py` — Numerik $N_\pm$, $D_E$, $\widetilde{D}_E$
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54 / PR #59

---

## 1. Korrelation auf Kanten: Paarstatistik $E(a,b)$

**Setup.** Auf dem gerichteten EABC-Transportgraphen $G_E=(V,E)$ mit $V=\{E,A,B,C\}$ betrachten wir **Paarstatistiken** entlang der Primfolge:
$$X_n := \kappa(p_n),\qquad \tau_n := (X_n, X_{n+1}).$$

**Kanten-Korrelation (Bell-Lesart, kombinatorisch).** Für Kantenlabels $a,b\in V$:
$$E(a,b) := \frac{\#\{n:\,\tau_n=(a,b),\,p_{n+1}\le X\}}{\#\{n:\,p_{n+1}\le X\}}$$
(bzw. bedingte Versionen auf gemeinsamen Trägern; Details: `collatz_eabc_bell_holonomie.md`).

**Interpretation (rein mathematisch).** $E(a,b)$ misst **bivariate Häufigkeit** von aufeinanderfolgenden EABC-Klassen — eine **Korrelationsstatistik auf gerichteten Kanten**, nicht eine physikalische Messachse.

$$\boxed{\;\text{Stufe 1: Korrelation} = \text{Paarstatistik } E(a,b) \text{ auf Kanten von } G_E.\;}$$

**Label:** $E(a,b)$ = **Definition** (kombinatorisch); physikalische Bell-Nichtlokalität = **nicht behauptet**.

---

## 2. Zyklusstatistik: $D_E(X) = \sum_\gamma \mathrm{sgn}(\gamma)$

**Geschlossene 5-Fenster** auf der Primfolge:
$$C_n^{(5)} := (X_n, X_{n+1}, X_{n+2}, X_{n+3}, X_{n+4}).$$

**Erkannte Orientierungen** (zwei Klassen in $H_1(C_4,\mathbb{Z})$):
- **$\gamma^+$:** ABCEA ($E\!\to\!A\!\to\!B\!\to\!C\!\to\!E$ als zyklisch verschobene Lesart)
- **$\gamma^-$:** CEABC

**Zählgrößen** (Prim-Obergrenze $X$):
$$N_+(X) := \#\{n:\,p_{n+4}\le X,\; C_n^{(5)}=\mathrm{ABCEA}\},$$
$$N_-(X) := \#\{n:\,p_{n+4}\le X,\; C_n^{(5)}=\mathrm{CEABC}\}.$$

**Zyklus-Observable (algebraische Topologie, nicht Bell):**
$$D_E(X) := N_+(X) - N_-(X) = \sum_{\gamma \le X} \mathrm{sgn}(\gamma),$$
wobei $\mathrm{sgn}(\gamma^+)=+1$, $\mathrm{sgn}(\gamma^-)=-1$, und die Summe nur über **erkannte** Orientierungen läuft.

$$\boxed{\;D_E(X) = \sum_{\gamma} \mathrm{sgn}(\gamma)\;\text{— Zyklusstatistik in } H_1(C_4,\mathbb{Z})\text{, nicht Kantenkorrelation.}\;}$$

**Label:** $N_\pm$, $D_E$ = **Definition**.

---

## 3. Der gerichtete $C_4$-Zyklus und Homologie

**Zustandsraum:** $V=\{E,A,B,C\}$.

**Kanonischer gerichteter 4-Zyklus:**
$$E \xrightarrow{} A \xrightarrow{} B \xrightarrow{} C \xrightarrow{} E.$$

**Zwei Orientierungen in $H_1(C_4,\mathbb{Z})$:**
| Orientierung | Wort | Vorzeichen |
|--------------|------|------------|
| $\gamma^+$ | ABCEA | $+1$ |
| $\gamma^-$ | CEABC | $-1$ |

Beide Wörter tragen dasselbe Lückenmuster $(2,4,2,4)$ mod $12$; sie unterscheiden sich durch **zyklische Verschiebung** der Orientierung auf demselben $C_4$-Gerüst.

**Kantenorientierung $\omega(e)\in\{+1,-1,0\}$** auf kanonischen Zykluskanten:
$$\omega(E\!\to\!A)=\omega(A\!\to\!B)=\omega(B\!\to\!C)=\omega(C\!\to\!E)=+1,$$
$$\omega(A\!\to\!E)=\omega(B\!\to\!A)=\omega(C\!\to\!B)=\omega(E\!\to\!C)=-1.$$

**Label:** $C_4$, $\omega$, $\gamma^\pm$ = **Definition**.

---

## 4. Diskrete 1-Form und Zirkulation $C_E(X)$

**Diskrete 1-Form** $\alpha(e)$ auf Zykluskanten:
$$\alpha(i\!\to\!j) := \frac{\omega(i,j)}{4}\quad\text{auf } E\!\to\!A,\,A\!\to\!B,\,B\!\to\!C,\,C\!\to\!E \text{ und Rückkanten}.$$

**Linienintegral** entlang geschlossenem Wort $\gamma$:
$$\oint_\gamma \alpha := \mathrm{sgn}(\gamma)\sum_{e\in\gamma}|\alpha(e)|.$$

Dann $\oint_{\mathrm{ABCEA}}\alpha = +1$, $\oint_{\mathrm{CEABC}}\alpha = -1$.

**Kumulative Zirkulation** (Summe über erkannte Zyklen bis $X$):
$$C_E(X) := \sum_{\gamma \le X} \oint_\gamma \alpha = \sum_{\gamma \le X} \mathrm{sgn}(\gamma) = D_E(X).$$

**Normalisierung:**
$$S_E(X) := \frac{C_E(X)}{N_+(X)+N_-(X)} = \frac{D_E(X)}{N_+(X)+N_-(X)}.$$

$$\boxed{\;C_E(X) = \sum_{\gamma \le X} \oint_\gamma \alpha = D_E(X).\;}$$

**Label:** $\alpha$, $C_E$ = **Definition**.

**Experiment:** `collatz_eabc_sagnac_circulation.py` — `circulation_C_E(X)`.

---

## 5. Fehlerterm: Erwartungswert, Zentrierung, Prime-Race-Struktur

Unter **keiner bevorzugten Orientierung** (Symmetrie-Hypothese mod $12$):
$$\mathbb{E}[C_E(X)] = 0.$$

**Zentrierter Fehlerterm:**
$$D_E(X) = C_E(X) - \mathbb{E}[C_E(X)] \quad\text{(bei } \mathbb{E}=0\text{: } D_E = C_E\text{).}$$

**Strukturelle Analogie** (nicht Theorem): dieselbe Fehlerterm-Architektur wie bei
- Chebyshev-Races ($\pi(x;4,3)-\pi(x;4,1)$),
- mod-$q$-Races,
- $\pi(x)$-Fehlertermen gegenüber der Hauptterm-Asymptotik.

**Hauptvermutung (Hauptterm).**
$$N_+(X) \sim N_-(X)\quad\Rightarrow\quad \lim_{X\to\infty} S_E(X) = 0.$$

**Fehlerterm-Hypothese (stärker).** $D_E$ trägt nichttrivialen Chebyshev-artigen Bias, gesteuert durch Nullstellen der Dirichlet-$L$-Funktionen modulo $12$.

**Normalisierter Fehlerterm:**
$$\widetilde{D}_E(X) := \frac{D_E(X)}{\sqrt{N_+(X)+N_-(X)}}.$$

$$\boxed{\;\mathbb{E}[C_E]=0 \text{ im Hauptterm; } D_E(X) \text{ trägt die interessante Fehlerterm-Struktur.}\;}$$

**Label:** Hauptvermutung = **Vermutung**; Fehlerterm-Hypothese = **Hypothese**; Prime-Race-Analogie = **Analogie**.

---

## 6. Spektralgeometrie: Graph-Laplace $L_E$ und $\mathrm{Spec}(L_E)$

**Gerichteter Transportgraph** aus der Primfolge: Kanten $i\to j$ gezählt, wenn $\kappa(p_n)=i$, $\kappa(p_{n+1})=j$.

**Adjazenzmatrix** $A_E$ (gewichtet nach Übergangszählungen $T_{ij}(X)$).

**Graph-Laplace** (gerichtet):
$$L_E = D_{\mathrm{out}} - A_E,$$
mit $D_{\mathrm{out}}=\mathrm{diag}(\sum_j A_{ij})$.

**Symmetrisierte Variante** (für reelle Eigenwerte):
$$L_E^{\mathrm{sym}} = D - \tfrac{1}{2}(A_E + A_E^\top),\qquad D=\mathrm{diag}(A_E\mathbf{1}).$$

**Forschungskette:**
$$\boxed{\;\text{Primzahlen} \;\to\; \text{EABC-Zyklen} \;\to\; G_E \;\to\; \mathrm{Spec}(L_E).\;}$$

**Spektralinvarianten:** Eigenwerte $\lambda_0\le\cdots\le\lambda_{|V|-1}$, Spektrallücke $\lambda_1-\lambda_0$, Fiedler-Vektor (bei Symmetrisierung).

**Label:** $L_E$, $\mathrm{Spec}(L_E)$ = **Definition**; Bezug zu Primstruktur = **Hypothese** / **Experiment**.

**Experiment:** `collatz_eabc_graph_laplacian.py` — `spectral_report(X)`.

---

## 7. Normschalen: $\Sigma_n \to G_n \to L_n \to \mathrm{Spec}(L_n)$

**Schalenfolge** (`collatz_eabc_quaternion_mass_hypothese.md`):
$$\mathcal{H} = \bigsqcup_{n\ge 1} \Sigma_n,$$
mit Normniveau $n$, lokalem Transportgraphen $G_n$ und Laplace $L_n$.

**Kette:**
$$\Sigma_n \;\longrightarrow\; G_n \;\longrightarrow\; L_n \;\longrightarrow\; \mathrm{Spec}(L_n).$$

**Primzahlen als Normschalen-Anomalien:** $n$ ist **Prim-Normniveau**, wenn das Spektrum $\mathrm{Spec}(L_n)$ vom Referenzspektrum $I_{\mathrm{ref}}(n)$ abweicht — analog zum Defekt $D(n)=I(\mu_n)-I_{\mathrm{ref}}(n)$.

**Brücke zur Zirkulation:** $C_E(X)$ bzw. $D_E(X)$ als **diskrete Zirkulation** auf der globalen Primfolge; $\mathrm{Spec}(L_E(X))$ als **spektrale** Lesart desselben Transportobjekts bei Obergrenze $X$.

**Label:** Normschalen-Kette = **Definition** / **Analogie**; Prim = Spektralanomalie = **Hypothese** (Emergenz, nicht Voraussetzung).

---

## 8. Boxed Forschungsfragen

$$\boxed{\;\text{Besitzt die durch die Primzahlfolge induzierte EABC-Zirkulation einen nichttrivialen Fehlerterm?}\;}$$

$$\boxed{\;\text{Welche Spektralinvarianten der EABC-Zirkulation unterscheiden Prim- von Nichtprim-Normschalen?}\;}$$

**Verwandte offene Fragen:**
- Verhält sich $\widetilde{D}_E(X)$ wie reines Rauschen oder zeigt stabile Vorzeichenasymmetrie (Chebyshev-artig)?
- Korreliert die Spektrallücke $\lambda_1-\lambda_0$ von $L_E(X)$ mit $|D_E(X)|$ oder $S_E(X)$?
- Welche $L$-Funktions-Nullstellen mod $12$ erklären Oszillationen in $D_E$?

---

## 9. Hierarchie (rein mathematisch)

$$\boxed{\;\text{Korrelation } E(a,b) \;\to\; \text{Zirkulation } C_E \;\to\; \text{Fehlerterm } D_E \;\to\; \mathrm{Spec}(L_E).\;}$$

| Stufe | Objekt | Symbol |
|------:|--------|--------|
| 1 | Primzahlen | $p_n$ |
| 2 | EABC-Klasse | $X_n=\kappa(p_n)$ |
| 3 | Kanten-Korrelation | $E(a,b)$, $T_{ij}$ |
| 4 | $C_4$-Orientierung | $\gamma^\pm$, $\omega(e)$ |
| 5 | 1-Form / Zirkulation | $\alpha$, $C_E(X)$ |
| 6 | Fehlerterm | $D_E(X)$, $\widetilde{D}_E(X)$ |
| 7 | Graph-Laplace | $L_E$, $\mathrm{Spec}(L_E)$ |
| 8 | Normschale | $\Sigma_n$, $G_n$, $L_n$ |

---

## 10. Python-Symbolzuordnung

| LaTeX | Python | Modul |
|-------|--------|-------|
| $C_E(X)$, $D_E(X)$ | `C_E`, `D_E` | `collatz_eabc_sagnac_circulation` |
| $\omega(e)$, $\alpha$ | `edge_omega`, `discrete_one_form` | `collatz_eabc_sagnac_circulation` |
| $T_{ij}$, $A_E$ | `transition_counts`, `adjacency_matrix` | `collatz_eabc_graph_laplacian` |
| $L_E$, $\mathrm{Spec}(L_E)$ | `laplacian_directed`, `eigenvalues` | `collatz_eabc_graph_laplacian` |
| Spektrallücke | `spectral_gap` | `collatz_eabc_graph_laplacian` |

---

## 11. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $E(a,b)$, $T_{ij}$, $N_\pm$, $D_E$, $C_E$, $\alpha$, $\omega$ | **Definition** |
| $L_E$, $\mathrm{Spec}(L_E)$, Spektrallücke | **Definition** |
| $N_+\sim N_-$ $\Rightarrow$ $S_E\to 0$ | **Vermutung** |
| $D_E$ mit $L$-Funktionen mod $12$ | **Hypothese** |
| Prime-Race / $\pi(x)$-Fehlerterm-Analogie | **Analogie** |
| Prim-Normschalen vs. Referenzspektrum | **Hypothese** |
| Bell / Sagnac | **Analogie** (sekundär, didaktisch) |
| Physikalische Rotation / Nichtlokalität | **nicht behauptet** |

---

*Kanonsiche Lesart: Das EABC-Programm ruht auf **Korrelation → Zirkulation** — Paarstatistik auf Kanten, Zyklusstatistik $D_E$ in $H_1(C_4,\mathbb{Z})$, diskrete 1-Form $\alpha$, Fehlerterm-Struktur wie bei Prime Races, und Spektralgeometrie via $\mathrm{Spec}(L_E)$. Didaktische Physik-Metaphern (Sagnac, Bell) verweisen hierher.*
