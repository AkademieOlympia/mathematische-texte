# EABC-Zirkulationshypothese

**Status:** Detailhypothese (Zählgrößen, Fehlerterm, Prime Race) — **primäre kanonische Formulierung:** `collatz_eabc_diskrete_geometrie.md`  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Vermutung | Hypothese | Experiment | Analogie

$$\boxed{\;\text{Primärdokument: } \texttt{collatz\_eabc\_diskrete\_geometrie.md} \text{ — } \Phi_E,\; E^\pm,\; \langle\omega_E,h\rangle.\;}$$

**Querverweise:**
- `collatz_eabc_diskrete_geometrie.md` — **kanonisch:** $G_E$, $E^+$, $E^-$, $\Phi_E$, EABC-Vermutung, drei Ebenen
- `collatz_eabc_holonomie_stufen.md` — drei Stufen (Analogie / echte Holonomie / Wilson) + Fall A/B/C in $N$
- `collatz_eabc_epistemik_physik.md` — **kanonische Abgrenzung:** Holonomie/Zirkulation ja; Zwillingsparadoxon/Zeitdilatation nein
- `collatz_eabc_epistemik_schichten.md` — Schichten A/B/C/R; asymptotische Chiralität methodisch in §4.1; Stufen 0–3 in §4.2
- `collatz_eabc_zirkulation_spektral.md` — Spektralgeometrie, diskrete 1-Form $\alpha$, $\mathrm{Spec}(L_E)$
- `collatz_eabc_fehlerterm_hypothese.md` — **Teilhypothese:** Fehlerterm $D_E$, $\widetilde{D}_E$ (eingebettet in §5)
- `collatz_eabc_sagnac.md` — **Intuition only:** Sagnac-Bild für $\gamma^\pm$ (kein physikalischer Kern)
- `collatz_eabc_zyklus_holonomie.md` — Hierarchie Klasse→Kante→Pfad→Zyklus→Holonomie
- `collatz_eabc_transport.md` — $G_E$, Transport $T_n$
- `collatz_eabc_bell_holonomie.md` — **sekundäre** Analogie Bell/CHSH (nicht primäre Lesart)
- `collatz_eabc_holonomie_beweisversuch.md` — analytischer Beweisversuch
- `collatz_eabc_sagnac_circulation.py` — $C_E(X)$, $\omega(e)$, $\alpha$
- `collatz_eabc_holonomie_fehlerterm.py` — $N_\pm$, $D_E$, $S_E$, $\widetilde{D}_E$
- `eabc_quadruplets_1e10.py` — Vierlings-Zählung bis $X$, $W_E$, $R_\beta$ (`Z_E` nur Diagnose-Alias), mod-$420$-Diagnostik
- `eabc_quadruplets_fit_alpha.py` — Stufen 0–3: Diagnostik (nicht Theorie); $\alpha_E$-Plateau, H₀a/H₀b–H₃
- `eabc_quadruplets_plot.py` — Vierfeld-Diagnose-Plot ($W_E$, $R_{1/2}$, $\alpha_{\mathrm{loc}}$, $R_\beta$)
- `collatz_eabc_graph_laplacian.py` — $\mathrm{Spec}(L_E)$
- `collatz_eabc_evolution_analytik.md` — **Evolutionspfad** Bell→Sagnac→$C_E$→$\mathrm{Spec}(L_E)$, Wachstumsszenarien, Dirichlet-Stub
- `collatz_eabc_D_growth.py` — Wachstumsdiagnostik $D_E(X)$, Charakter-Koeffizienten $a_\chi$
- `collatz_eabc_kritische_abbildung.md` — **Modellabbildung:** Geschwindigkeitsmodell $s_v(x)$, ABCEA/CEABC in der komplexen Ebene; **Holonomie-Sensor** $v_j=\gamma_{\mathrm{ref}}/\ell_j$
- `collatz_eabc_chirale_polarisation.md` — **Helizität** $\lambda=\pm 1$, $N_R/N_L$, $\phi_R/\phi_L$, Stufe-2-Upgrade
- `collatz_eabc_brachistochrone.md` — $T_R$, $T_L$, Birefringenz-Analogie
- `collatz_eabc_wigner_analog.md` — Wigner-Analogie: $W_E$ (4-Pfad) vs. $D_E$ (5-Zyklus), $W_{ab}$
- `collatz_eabc_uebergangsraum.md` — **Detail:** $C_4\cong S^1$, $\langle\omega_E,h\rangle$, $L_{\mathrm{mag}}$
- `collatz_eabc_signierte_massstruktur.md` — signierte Maßstruktur auf $G_E$, arithmetische Wigner-Negativität
- `collatz_eabc_hodge_eabc.py` — `Phi_E`, `flux_density_limit`, `harmonic_holonomy_component`, `magnetic_laplacian`
- `collatz_eabc_kritische_abbildung.py` — Numerik $x_{n,v}$, Schaltkreis-Trajektorien, `edge_velocities_from_gaps`, `holonomy_sensor_trajectory`
- `collatz_generalangriff_2026.md` — Gesamtarchitektur PR #54 / PR #59

---

## 0. Boxed These (verweist auf kanonisches Dokument)

$$\boxed{\;\text{EABC = Theorie des priminduzierten Flusses auf } C_4\cong S^1 \text{ — nicht Zählstatistik ABCE/CEAB.}\;}$$

Vollständige paper-ready Definitionen ($G_E$, $E^\pm$, $\Phi_E$, Vermutungen): **`collatz_eabc_diskrete_geometrie.md`**.

$$\boxed{\;\text{Shift: } E(a,b) \;\to\; \oint_\gamma \alpha \quad\text{(Paarstatistik } \to \text{ Zyklusgeometrie).}\;}$$

**Epistemische Abgrenzung:** Das EABC-Programm ist **keine** Bell-Korrelationstheorie und **kein** quantenphysikalischer Effekt. Es ist **diskrete Zirkulation** auf einem gerichteten Restklassengraphen.

---

## 1. Setup

**Zustandsraum:**
$$V = \{E, A, B, C\},\qquad E\equiv 1,\; A\equiv 5,\; B\equiv 7,\; C\equiv 11 \pmod{12}.$$

**Primfolge-Labels:**
$$X_n := \kappa(p_n),$$
wobei $(p_n)_{n\ge 1}$ die aufsteigende Folge aller Primzahlen $>3$ ist.

**Gerichteter Pfad auf $V$:** Aufeinanderfolgende Transportkanten
$$\tau_n := (X_n, X_{n+1}) = (\kappa(p_n), \kappa(p_{n+1}))$$
bilden einen **gerichteten Pfad** auf dem EABC-Transportgraphen $G_E=(V,E)$.

**Elementarer Zyklus:**
$$A \longrightarrow B \longrightarrow C \longrightarrow E \longrightarrow A.$$

**Label:** $V$, $X_n$, $\tau_n$, elementarer Zyklus = **Definition**.

---

## 2. Zwei Orientierungen

Auf demselben $C_4$-Gerüst $A\!\to\!B\!\to\!C\!\to\!E\!\to\!A$ betrachten wir zwei **Orientierungen** in $H_1(C_4,\mathbb{Z})$:

| Orientierung | Wort | Vorzeichen |
|--------------|------|------------|
| $\gamma^+$ | ABCEA | $+1$ | $R$ ($\lambda=+1$) |
| $\gamma^-$ | CEABC | $-1$ | $L$ ($\lambda=-1$) |

**Chirale Lesart** (`collatz_eabc_chirale_polarisation.md`): $C_E = \sum \omega(\gamma) = N_R - N_L$ als diskreter **Chiralitätsfluss** / Polarisationsoperator — nicht nur Umlaufzahl.

Beide Wörter tragen dasselbe Lückenmuster $(2,4,2,4)$ mod $12$; sie unterscheiden sich durch **zyklische Verschiebung** der Orientierung auf demselben Zyklus.

**Geschlossenes 5-Fenster:**
$$C_n^{(5)} := (X_n, X_{n+1}, X_{n+2}, X_{n+3}, X_{n+4}).$$

**Label:** $\gamma^\pm$, ABCEA, CEABC = **Definition**.

---

## 3. Zählgrößen und Zirkulation

Für Prim-Obergrenze $X$:

$$N_+(X) := \#\{n:\,p_{n+4}\le X,\; C_n^{(5)}=\mathrm{ABCEA}\},$$
$$N_-(X) := \#\{n:\,p_{n+4}\le X,\; C_n^{(5)}=\mathrm{CEABC}\}.$$

**Legacy-Aliase:** $N_{\mathrm{ABCEA}}:=N_+$, $N_{\mathrm{CEABC}}:=N_-$.

**Zirkulation:**
$$C_E(X) := N_+(X) - N_-(X).$$

**Normalisierte Observable:**
$$S_E(X) := \frac{N_+(X) - N_-(X)}{N_+(X) + N_-(X)} = \frac{C_E(X)}{N_+(X)+N_-(X)}.$$

**Fehlerterm (Alias):**
$$D_E(X) := N_+(X) - N_-(X) = C_E(X).$$

**Normalisierter Fehlerterm:**
$$\widetilde{D}_E(X) := \frac{D_E(X)}{\sqrt{N_+(X)+N_-(X)}}.$$

**Legacy-Aliase:** $\Delta_E(X):=D_E(X)$; $\chi_{\mathrm{Hol}}(X):=S_E(X)$; $\mathrm{Hol}_E:=\lim_{X\to\infty} S_E(X)$.

$$\boxed{\;C_E(X) = N_+(X) - N_-(X) = D_E(X).\;}$$

**Label:** $N_\pm$, $C_E$, $S_E$, $D_E$, $\widetilde{D}_E$ = **Definition**.

---

## 4. Hauptterm-Vermutung und Zentralvermutung (Flussdichte)

**Hauptterm-Vermutung.** Asymptotische Symmetrie der gegenläufigen Zyklusorientierungen:
$$N_+(X) \sim N_-(X)\qquad (X\to\infty),$$
und damit
$$S_E(X) \to 0.$$

**Zentralvermutung (arithmetische Orientierungsklasse).** Normalisierter Magnetfluss (`collatz_eabc_diskrete_geometrie.md` §2):
$$\boxed{\;\Phi_E = \lim_{X\to\infty} W_E(X) = \lim_{X\to\infty}\frac{C_E(X)}{N_+(X)+N_-(X)} \;\stackrel{?}{\neq}\; 0.\;}$$

| Grenzfall | Folgerung |
|-----------|-----------|
| $\lim S_E = 0$ | keine bevorzugte globale Orientierung (Hauptterm) |
| $\lim S_E \neq 0$ | nichttriviale **arithmetische Orientierungsklasse** |

**Numerik:** `flux_density_limit` in `collatz_eabc_hodge_eabc.py`.

**Label:** Hauptterm-Vermutung = **Vermutung**; Orientierungsklassen-Vermutung = **Vermutung** / **Forschungsfrage**.

### 4.1 Asymptotische Chiralität (methodische Lesart)

$$\boxed{\;\text{Die EABC-Holonomie-Vermutung behauptet eine asymptotisch stabile, nichtverschwindende Richtungspräferenz der priminduzierten Zyklen auf dem EABC-Kreisgraphen.}\;}$$

**Kerngröße** (Schicht **B**, Definition):

$$W_E(X) = \frac{N_+(X)-N_-(X)}{N_+(X)+N_-(X)}$$

| Status | Aussage | Lean (`HolonomyCore.lean`) |
|--------|---------|---------------------------|
| **GREEN** | $-1 \le W_E(X) \le 1$ für alle $X$ | `W_E_bounds` (bewiesen) |
| **RED** | $\displaystyle\lim_{X\to\infty} W_E(X) = \Phi_E \neq 0$ | `HasNonzeroHolonomyLimit`, `EABC_holonomy_limit_conjecture` (`sorry`) |

**Asymptotische Chiralität** heißt: $W_E(X)$ **stabilisiert sich asymptotisch gegen einen Grenzwert** $\Phi_E$ — eine **asymptotisch nichtverschwindende Orientierungspräferenz** auf $C_4$. Das ist weder „Ruhe“ der Einzelprimzahlen noch eine metaphorische „statistische Ruhe im Unendlichen“; lokal bleibt die Primfolge chaotisch.

| Lesart | Inhalt | Label |
|--------|--------|-------|
| **Definition** | $N_\pm(X)$ zählt ABCEA- bzw. CEABC-Fenster unter der Bedingung $p_{n+4}\le X$; $W_E(X)$ ist der Quotient oben. Primzahlen **laufen** nicht — es gibt nur eine **priminduzierte EABC-Zählung** auf $G_E$. | **Definition** (Schicht **B**) |
| **Vermutung (stark)** | $\exists\,\Phi\neq 0:\; W_E(X)\to\Phi$ — asymptotisch nichtverschwindende Orientierungspräferenz / asymptotische Chiralität. | **Vermutung** (Schicht **R**) |
| **Vermutung (konservativ, Hauptterm)** | $N_+(X)\sim N_-(X)\;\Rightarrow\;\Phi_E=0$ — Symmetrie der gegenläufigen Zyklusorientierungen. | **Vermutung** (Schicht **B**/**R**) |
| **Experiment** | Numerik bei endlichem $X$ (z. B. $X\approx 10^6$, $W_E\approx 0{,}12$) ist **Stichprobe**, kein Beweis des Grenzwerts. | **Experiment** |

**Methodische Präzisionen (populärwissenschaftliche Lesart vermeiden):**

1. **Kein Lauf-Bild.** Metaphern wie „Primzahlen laufen im Kreis“ gehören zur Ikone-Spur (Schicht **C**). Mathematisch fixiert ist nur $N_\pm(X)$ und der daraus gebildete Quotient $W_E(X)$ auf Schicht **B**.
2. **Stabilisierung des Quotienten.** Selbst wenn $W_E(X)\to\Phi_E\neq 0$, bleibt die Primfolge lokal chaotisch; **$W_E(X)$ stabilisiert sich asymptotisch gegen einen Grenzwert**, nicht die Einzelereignisse $C_n^{(5)}$.
3. **Keine unumkehrbare Tendenz.** Eine nichttriviale Grenzpräferenz $\Phi_E\neq 0$ ist eine **asymptotisch nichtverschwindende Orientierungspräferenz** — keine behauptete „unumkehrbare“ physikalische Tendenz.
4. **Zwei konkurrierende asymptotische Lesarten.** Die **starke EABC-Vermutung** ($\Phi_E\neq 0$) und die **konservative Nullhypothese** ($N_+\sim N_-\Rightarrow\Phi_E=0$, vgl. Hauptterm oben) sind bewusst getrennt zu halten — Experimente entscheiden nicht zwischen ihnen.

### 4.2 Fehlerterm-Skalierung und Hypothesenhierarchie (H₀a/H₀b–H₃)

$$\boxed{\;\text{Das EABC-Programm untersucht zunächst die Größenordnung der chiralen Differenz }D(X)=A(X)-C(X)\text{, und erst in einem zweiten Schritt die Existenz einer asymptotischen Orientierung }W_E(X)=D(X)/Q(X).\;}$$

**Methodische Leitlinie (Tao-Stil):** Vier **Stufen** — Definition → numerisches Signal → Vermutung → Orientierung. Keine Exponentenannahme vor Stufe 2; Holonomie (Stufe 3) **am Ende**, nicht am Anfang. Die Hypothesen H₀a–H₃ referenzieren jeweils eine Stufe; sie ersetzen sie nicht.

#### Stufe 0 — Größenordnung (keine Exponentenannahme)

**Definition** (Schicht **B**; Vierlings-/Fensterzählung):

$$A(X) := N_+(X)\quad(\mathrm{ABCEA}),\qquad C(X) := N_-(X)\quad(\mathrm{CEABC}),$$
$$D(X) := A(X)-C(X) = D_E(X),\qquad Q(X) := A(X)+C(X).$$

$$W_E(X) := \frac{D(X)}{Q(X)},\qquad R_\beta(X) := \frac{D(X)}{Q(X)^\beta}\quad(\beta\in\mathbb{R}_+).$$

$W_E(X)=D(X)/Q(X)$ ist **Definition** — kein Grenzwert, kein Exponent. Man untersucht $R_\beta(X)$ **ohne** die Existenz eines Exponenten vorauszusetzen: bleibt $R_{1/2}$ beschränkt? wächst $R_{2/3}$? Kein Zufallsmodell, keine Vermutung auf dieser Stufe.

**Kein** $Z_E=D/\sqrt{Q}$ auf der Definitionsebene — nur $D$, $Q$, $W_E$, $R_\beta$. Spezialfälle: $R_1(X)=W_E(X)$; $R_{1/2}$, $R_{2/3}$, $R_{3/4}$, $R_{0.9}$ sind **diagnostische** Normierungen (Stufe 0), keine Kerngrößen.

| Begriff | Rolle | Observable |
|---------|-------|------------|
| **$D$, $Q$, $R_\beta$** | Größenordnung | $R_\beta=D/Q^\beta$ |
| **$W_E$** | normierte Differenz (Definition) | $D/Q$ |
| **$\Phi_E$** | asymptotische Orientierung (Stufe 3) | $\lim W_E$ (falls existent) |
| **$\alpha_E$** | Skalierungsexponent (Stufe 2) | $\inf\{\beta : R_\beta\text{ beschränkt}\}$ |
| **$\alpha_{\mathrm{eff}}(X)$** | numerisches Signal (Stufe 1) | $d\log|D|/d\log Q$ |

**Label:** **Definition** (Stufe 0).

#### Heuristik: $Z_E=R_{1/2}$, $\alpha_E=\tfrac{1}{2}$ (Diskussion, nach Stufe 0)

**Nur Diskussion** — kein Theorem, keine Definition:

- $Z_E := R_{1/2}(X)=D(X)/\sqrt{Q(X)}$ — in CSV/Plot als Spalte `R_1_2` bzw. Alias $Z_E$ (**Heuristik/Diagnose**).
- $\alpha_E=\tfrac{1}{2}$ entspräche einem **naiven Zufallsmodell** ($|D|\approx\sqrt{Q}$) — **Heuristik**, keine bewiesene Referenz.
- Bei Primzahlrennen können Differenzen **deutlich größer** als $\sqrt{Q}$ werden, ohne nichttrivialen normierten Grenzwert.

#### Stufe 1 — Effektive Skalierung (numerisch)

$$\alpha_{\mathrm{eff}}(X) := \frac{\log|D(X)|}{\log Q(X)}\quad (D\neq 0,\; Q>1),$$

$$\alpha_{\mathrm{loc}}(X_i,X_{i+1}) := \frac{\Delta\log|D|}{\Delta\log Q}\quad\text{zwischen aufeinanderfolgenden Checkpoints.}$$

**Frage (Stufe 1):** Stabilisiert sich $\alpha_{\mathrm{eff}}(X)$? — kann oszillieren, driften oder extrem langsam konvergieren. **Nur numerisches Signal**, kein asymptotischer Satz.

**Wichtiger Vorbehalt:** $\alpha_{\mathrm{loc}}$ **muss nicht konvergieren** — in der analytischen Zahlentheorie kann ein Exponent über lange Bereiche stabil erscheinen und später driften oder oszillieren (vgl. explizite Formeln, $L$-Funktionen).

**Label:** **Experiment** (numerische Ableitung).

#### Stufe 2 — Asymptotischer Exponent (Vermutung, erst wenn Daten es nahelegen)

$$\boxed{\;\text{Existiert }\alpha_E\text{ mit }|D(X)|=Q(X)^{\alpha_E+o(1)}\text{?}\;}$$

**Kritischer Exponent** (mathematisches Objekt über Renormierungsobservablen):

$$\alpha_E := \inf\{\beta\in\mathbb{R}_+ : R_\beta(X)\text{ bleibt beschränkt für }X\to\infty\}.$$

Erst hier „nichttriviale asymptotische Skalierung" formulieren — **entschärft** gegenüber einer sofortigen Exponentenannahme. Numerisch: man vergleicht, bei welchem $\beta$ die Kurve $R_\beta(X)$ im beobachteten Fenster noch **stabil** vs. **wächst**; Zwischenwerte zeigen die **Skala des Bias**.

**Testarchitektur** (Stufe 0/2):

| Observable | $\beta$ | Lesart |
|------------|---------|--------|
| $R_{1/2}$ | $\tfrac{1}{2}$ | Wurzelrauschen-Heuristik (H₀a/H₁) |
| $R_{2/3}$, $R_{3/4}$, $R_{0.9}$ | Zwischenwerte | sublinearer Bias (H₂-Diagnostik) |
| $R_1=W_E$ | $1$ | Orientierung (Stufe 3; H₀b/H₃) |

**Drei Skalierungsfälle** (an $\alpha_E$ gebunden, **unter H₂**):

| Fall | $\alpha_E$ | $R_{1/2}$ | $W_E$ | Lesart |
|------|------------|-----------|-------|--------|
| 1. Rauschen | $\le \tfrac{1}{2}$ | $O(1)$ | $\to 0$ | Wurzelrauschen (H₀a) |
| 2. Sublinearer Bias | $\tfrac{1}{2} < \alpha_E < 1$ | $\to\infty$ | $\to 0$ | arithmetischer Bias (H₁–H₂, nicht H₃) |
| 3. Starke Holonomie | $= 1$ | $\to\infty$ (langsamer als Fall 2) | $\to \Phi_E \neq 0$ | Holonomie-Grenzfall (H₃) |

Unter H₂ folgen formal: $W_E(X) \sim Q(X)^{\alpha_E-1}$, $R_{1/2}(X) \sim Q(X)^{\alpha_E-1/2}$, $R_\beta(X) \sim Q(X)^{\alpha_E-\beta}$.

**Wichtig:** Auch bei $\Phi_E=0$ kann $\alpha_E>\tfrac{1}{2}$ gelten — ein **eigenständiges arithmetisches Phänomen**, nicht äquivalent zur Holonomiefrage.

**Label:** $\alpha_E$ = **Definition** (Stufe 2); H₂ = **Hypothese**; Kurvenvergleich / $\alpha_E$-Schätzung = **Experiment**.

#### Stufe 3 — Orientierung (zweiter Schritt des Programms)

$$\boxed{\;\text{Orientierung: }W_E(X)\to 0\quad\text{oder}\quad W_E(X)\to\Phi_E\neq 0\text{?}\;}$$

$$\Phi_E := \lim_{X\to\infty} W_E(X)\quad\text{(falls der Grenzwert existiert).}$$

$W_E(X)$ und $\Phi_E$ messen die **normierte Orientierung** — unabhängig von jeder Skalierungsexponentenannahme. Holonomie **am Ende**, nicht am Anfang.

- **H₀b** (keine Orientierung): $W_E(X)\to 0$ — analytische Nullhypothese.
- **H₃** (starke Holonomie): $W_E(X)\to\Phi_E\neq 0$ — Vermutung (Lean-**RED**: `HasNonzeroHolonomyLimit`).

#### Implikationskette (Referee-präzise)

$$\boxed{\;\Phi_E\neq 0 \;\Rightarrow\; D(X)\sim\Phi_E\,Q(X) \;\Rightarrow\; \alpha_E=1.\;}$$

**Nicht umgekehrt:** $\alpha_E=1$ impliziert **keinen** Grenzwert von $W_E$.

**Gegenbeispiel:** $D(X)=Q(X)\sin(\log\log Q(X))$ — $\alpha_E=1$, aber $W_E(X)$ **ohne** Grenzwert.

| Parameter | Stärke | Bedeutung |
|-----------|--------|-----------|
| **$\Phi_E$** | Orientierungsparameter (stärker) | normierter Grenzwert $W_E$ |
| **$\alpha_E$** | Skalenparameter (schwächer) | Größenordnung von $|D|$ relativ zu $Q$ |

Weitere Vorwärtsimplikationen (H₀a/H₀b–H₃):

- **H₃ $\Rightarrow$ H₂** (mit $\alpha_E=1$), aber **H₂ $\nRightarrow$ H₃** (jedes $\alpha_E\in(\tfrac{1}{2},1)$ lässt $W_E\to 0$).
- **H₂ $\Rightarrow$ H₁** auf hinreichend großen Skalen ($R_{1/2}\to\infty$), aber **H₁ $\nRightarrow$ H₂** (transienter Bias reicht nicht).
- **H₀b $\nRightarrow$ H₀a** und **H₀a $\nRightarrow$ H₀b** — analytische Orientierung und Wurzelnormierung sind **getrennte** Fragen.
- **H₂ mit $\alpha_E>\tfrac{1}{2}$ $\nRightarrow$ H₃** — auch $\Phi_E=0$ ist mit $\alpha_E>\tfrac{1}{2}$ vereinbar.

#### Boxed Kernfragen (entschärft)

1. **Stufe 0/1:** Bleibt $R_{1/2}$ beschränkt? Stabilisiert sich $\alpha_{\mathrm{eff}}$? — **Größenordnung zuerst.**
2. **Stufe 2:** Existiert $\alpha_E$ mit $|D(X)|=Q(X)^{\alpha_E+o(1)}$? — **nur wenn Daten es nahelegen.**
3. **Stufe 3:** $W_E(X)\to\Phi_E\neq 0$? — **Orientierung danach.**

#### Hypothesen (H₀a/H₀b–H₃)

Schrittweise Eskalation entlang der Stufen: **0/1** Größenordnung → **2** Skalierung → **3** Orientierung.

| Hypothese | Aussage | Stufe | Label |
|-----------|---------|-------|-------|
| **H₀a** (Wurzelrauschen) | $\alpha_E \le \tfrac{1}{2}$ — $R_{1/2}(X)$ **bleibt beschränkt** | 0/2 | **Heuristik** (naives Zufallsmodell) |
| **H₁** (mehr als Wurzelrauschen) | $\alpha_E > \tfrac{1}{2}$ empirisch — $R_{1/2}(X)\to\infty$ | 0/1 | **Experiment** |
| **H₂** (asymptotischer Bias) | $|D(X)| \asymp Q(X)^{\alpha_E}$ mit $\alpha_E > \tfrac{1}{2}$ asymptotisch | 2 | **Hypothese** (arithmetisch stark) |
| **H₀b** (keine Orientierung) | $W_E(X)\to 0$ — analytische Nullhypothese; **H₀b $\nRightarrow$ $\alpha_E\le\tfrac{1}{2}$** | 3 | **Nullhypothese** (analytisch) |
| **H₃** (starke Holonomie) | $W_E(X)\to\Phi_E\neq 0$ — impliziert $D(X)\sim\Phi_E Q(X)$, $\alpha_E=1$ | 3 | **Vermutung** (Lean-**RED**: `HasNonzeroHolonomyLimit`) |

#### Prime-Race-Analogie

Wie bei $\Delta(x)=\pi(x;a,q)-\pi(x;b,q)$ fragt man zuerst nach der **Größenordnung** der Differenz relativ zu einer Normierung — nicht nach einem Grenzwert des Quotienten. $D(X)=N_+(X)-N_-(X)$ ist ein orientierter Prime Race auf demselben EABC-Zyklus. Die asymptotische Orientierungsfrage (Stufe 3: $W_E\to 0$ oder $\Phi_E\neq 0$) ist von der Skalierungsfrage (Stufe 2: $\alpha_E$ via $R_\beta$) **methodisch getrennt**.

**Numerik und Diagnose-Plot.** `eabc_quadruplets_1e10.py` erzeugt Checkpoints mit $Q$, $D$, $W_E$, $R_\beta$ für $\beta\in\{\tfrac{1}{2},\tfrac{2}{3},\tfrac{3}{4},0{,}9,1\}$ (CSV-Spalten `R_1_2`, `R_2_3`, `R_3_4`, `R_9_10`, `R_1`; Alias `Z_E`=`R_1_2` nur Diagnose), $\alpha_{\mathrm{eff}}$ (`eabc_quadruplets.csv`). `eabc_quadruplets_fit_alpha.py` strukturiert die Ausgabe nach **Stufen 0–3**; schätzt $\alpha_E$ heuristisch aus $R_\beta$-Plateaus (**Experiment**, kein Theorem). `eabc_quadruplets_plot.py` (oder `fit_alpha --plot`) erzeugt ein Vierfeld-Diagramm (`eabc_quadruplets_diagnose.png`):

| Panel | Kurve | Stufe / Lesart |
|-------|-------|----------------|
| 1 | $R_{1/2}(X)$ ($Z_E$ Alias) | Stufe 0/1: Größenordnung (H₀a / H₁) |
| 2 | $\alpha_{\mathrm{loc}}(X)$ | Stufe 1: $\alpha_{\mathrm{eff}}$-Diagnose, kein Satz |
| 3 | $R_\beta$ für $\beta\in\{\tfrac{1}{2},\tfrac{2}{3},\tfrac{3}{4},0{,}9,1\}$ | Stufe 2: Skalierungsdiagnostik |
| 4 | $W_E(X)=R_1(X)$ | Stufe 3: H₀b / H₃ (Orientierung) |

**Label:** Stufe 0 = **Definition**; Stufe 1 = **Experiment**; Stufe 2 / H₂ = **Hypothese**; Stufe 3 / H₃ = **Vermutung** (Lean-RED); H₀a = **Heuristik**; H₀b = **Nullhypothese**.

---

## 5. EABC-Hypothese (Fehlerterm)

**EABC-Hypothese (Fehlerterm).** Der Fehlerterm
$$D_E(X) = N_+(X) - N_-(X)$$
ist **kein** reines Rauschen. Er trägt einen **nichttrivialen Chebyshev-artigen Bias** bzw. **Oszillationsstruktur**, verknüpft mit den **Fehlertermen der Primzahlverteilung modulo $12$** (bzw. den zugehörigen Dirichlet-$L$-Funktionen auf $(\mathbb{Z}/12\mathbb{Z})^\times$).

Qualitativ: wie beim klassischen Chebyshev-Bias mod $4$ kann die **absolute Differenz** $D_E(X)$ vorzeichenbehaftet und strukturiert oszillieren, während der **normierte Hauptterm** $S_E(X)$ verschwindet.

**Zentrale Frage:** Verhält sich $\widetilde{D}_E(X)$ wie reines Rauschen, oder zeigt sie stabile Vorzeichenasymmetrie / Oszillation gegenüber Isotropie- und Shuffle-Nullmodellen?

**Label:** EABC-Fehlerterm-Hypothese = **Hypothese** (stärker als Hauptvermutung allein).

**Details und Numerik:** `collatz_eabc_fehlerterm_hypothese.md` §3–5; `collatz_eabc_holonomie_fehlerterm.py`; Stufen 0–3 und H₀a/H₀b–H₃: §4.2; `eabc_quadruplets_1e10.py`, `eabc_quadruplets_fit_alpha.py`, `eabc_quadruplets_plot.py`.

---

## 6. Abgrenzung: keine Bell-Korrelation, kein Quanteneffekt

**Nicht** Bell-Korrelationstheorie, sondern **diskrete Zirkulation** auf dem gerichteten Restklassengraphen $G_E$.

**ABCEA vs. CEABC** ist **kein** Quanteneffekt — es ist ein **orientierter Zyklus-Race** (Prime Race zwischen zwei Orientierungen desselben EABC-Zyklus).

$$\boxed{\;D_E(X)\;\text{ist ein Prime Race zwischen zwei Orientierungen desselben EABC-Zyklus.}\;}$$

**Sekundäre Analogien** (didaktisch, nicht primär):
- Bell/CHSH auf $G_E$: `collatz_eabc_bell_holonomie.md`
- Sagnac-Bild ($\gamma^\pm$): `collatz_eabc_sagnac.md` — **Intuition only** (Zirkulationstheorie, **keine** SRT; Skalierungsparameter $v$: `collatz_eabc_kritische_abbildung.md` §0.1)

**Label:** Zirkulation auf $G_E$ = **Definition**; Bell/Sagnac = **Analogie** (sekundär).

---

## 7. Verknüpfungen

Die EABC-Zirkulationshypothese verbindet:

| Schicht | Objekt | Dokument |
|---------|--------|----------|
| Fehlerterm | $D_E(X)$, $\widetilde{D}_E(X)$ | hier §5; `collatz_eabc_fehlerterm_hypothese.md` |
| Dirichlet-$L$ mod $12$ | Nullstellen, Charaktere | `collatz_eabc_holonomie_beweisversuch.md` |
| Graph-Zirkulation | $C_E$, $\alpha$, $\omega(e)$ | `collatz_eabc_zirkulation_spektral.md` §4 |
| Graph-Spektrum | $\mathrm{Spec}(L_E)$ | `collatz_eabc_zirkulation_spektral.md` §6; `collatz_eabc_graph_laplacian.py` |
| Wachstum / Dirichlet | $D_E(X)$-Szenarien A–D ($X$), Fall A/B/C ($N$) | `collatz_eabc_evolution_analytik.md` §3; `collatz_eabc_holonomie_stufen.md` §4; `collatz_eabc_D_growth.py` |
| Kritische Abbildung | $s_v(x)$, $x_{n,v}$, Holonomie-Schaltkreis, **Sensor** $v_j=\gamma_{\mathrm{ref}}/\ell_j$ | `collatz_eabc_kritische_abbildung.md`; `collatz_eabc_kritische_abbildung.py` |
| Chirale Polarisation | $N_R$, $N_L$, $\phi_R$, $\phi_L$, $U_E$ | `collatz_eabc_chirale_polarisation.md`; `collatz_eabc_chirale_transport.py` |
| Brachistochrone | $T_R$, $T_L$, Birefringenz | `collatz_eabc_brachistochrone.md`; `collatz_eabc_brachistochrone.py` |

$$\boxed{\;\text{Fehlerterm } D_E \;\leftrightarrow\; L\text{-Funktionen mod } 12 \;\leftrightarrow\; \text{Zirkulation } C_E \;\leftrightarrow\; \mathrm{Spec}(L_E).\;}$$

---

## 8. Boxed Schlussfolgerungen

$$\boxed{\;S_E(X) \to 0\;\text{im Hauptterm, aber}\;D_E(X)\;\text{kann nichttrivial sein.}\;}$$

$$\boxed{\;\text{ABCEA/CEABC-Asymmetrie} \;\Rightarrow\; \text{Fehlerterm/Bias, nicht Hauptterm-Holonomie.}\;}$$

$$\boxed{\;\text{Besitzt die durch die Primzahlfolge induzierte EABC-Zirkulation einen nichttrivialen Fehlerterm?}\;}$$

$$\boxed{\;\text{Welche Spektralinvarianten der EABC-Zirkulation unterscheiden Prim- von Nichtprim-Normschalen?}\;}$$

---

## 9. Epistemische Tabelle

| Aussage | Label |
|---------|-------|
| $V$, $X_n$, $\tau_n$, elementarer Zyklus | **Definition** |
| $\gamma^\pm$, $N_\pm$, $C_E$, $S_E$, $D_E$, $\widetilde{D}_E$ | **Definition** |
| $S_E(X)\to 0$ | **Vermutung** (Hauptterm) |
| $D_E$ mit Chebyshev-Bias / $L$-Funktionen mod $12$ | **Hypothese** |
| Verhalten von $\widetilde{D}_E$ | **Experiment** |
| Bell / Sagnac | **Analogie** (sekundär, didaktisch) |
| Physikalische Rotation / Nichtlokalität | **nicht behauptet** |

---

## 10. Python-Symbolzuordnung

| LaTeX | Python | Modul |
|-------|--------|-------|
| $C_E(X)$ | `C_E` | `collatz_eabc_sagnac_circulation` |
| $D_E(X)$, $\Delta_E(X)$ | `D_E`, `Delta_E` | `collatz_eabc_holonomie_fehlerterm` |
| $S_E(X)$ | `S_E` | `collatz_eabc_holonomie_fehlerterm` |
| $\widetilde{D}_E(X)$ | `D_tilde_E` | `collatz_eabc_holonomie_fehlerterm` |
| $N_\pm(X)$ | `N_plus`, `N_minus` | `collatz_eabc_holonomie_fehlerterm` |
| $\omega(e)$, $\alpha$ | `edge_omega`, `discrete_one_form` | `collatz_eabc_sagnac_circulation` |
| $\mathrm{Spec}(L_E)$ | `eigenvalues_symmetrized` | `collatz_eabc_graph_laplacian` |
| $\Phi_E$, $W_E(X)$, Flussdichte | `Phi_E`, `flux_density_limit` | `collatz_eabc_hodge_eabc` |
| $Q$, $D$, $W_E$, $Z_E$ (Vierlinge) | `Q_total`, `diff`, `W_E`, `Z_E` | `eabc_quadruplets_1e10` |
| Skalierungsexponent $\alpha_E$ | `fit_alpha`, `estimate_alpha_E` | `eabc_quadruplets_fit_alpha` |
| $\langle\omega_E,h\rangle$ | `inner_product_omega_h` | `collatz_eabc_hodge_eabc` |
| $L_{\mathrm{mag}}$ | `magnetic_laplacian` | `collatz_eabc_hodge_eabc` |
| $N_R$, $N_L$, $\phi_R$, $\phi_L$ | `N_R`, `N_L`, `phi_R`, `phi_L` | `collatz_eabc_chirale_transport` |
| $T_R$, $T_L$ | `T_R`, `T_L` | `collatz_eabc_brachistochrone` |

---

*Kanonsiche Hypothese: EABC als Zirkulationstheorie — Shift von Paarstatistik $E(a,b)$ zu Zyklusstatistik $\oint_\gamma \alpha$. Spektralgeometrie und diskrete 1-Form: `collatz_eabc_zirkulation_spektral.md`. Sagnac nur als Intuition: `collatz_eabc_sagnac.md`.*
