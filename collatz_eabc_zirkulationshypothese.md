# EABC-Zirkulationshypothese

**Status:** Detailhypothese (Zählgrößen, Fehlerterm, Prime Race) — **primäre kanonische Formulierung:** `collatz_eabc_diskrete_geometrie.md`  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59)  
**Tao-Labels:** Definition | Vermutung | Hypothese | Experiment | Analogie

$$\boxed{\;\text{Primärdokument: } \texttt{collatz\_eabc\_diskrete\_geometrie.md} \text{ — } \Phi_E,\; E^\pm,\; \langle\omega_E,h\rangle.\;}$$

**Querverweise:**
- `collatz_eabc_diskrete_geometrie.md` — **kanonisch:** $G_E$, $E^+$, $E^-$, $\Phi_E$, EABC-Vermutung, drei Ebenen
- `collatz_eabc_holonomie_stufen.md` — drei Stufen (Analogie / echte Holonomie / Wilson) + Fall A/B/C in $N$
- `collatz_eabc_epistemik_physik.md` — **kanonische Abgrenzung:** Holonomie/Zirkulation ja; Zwillingsparadoxon/Zeitdilatation nein
- `collatz_eabc_epistemik_schichten.md` — Schichten A/B/C/R; Lakatos-Einordnung und Ebenen 0–4 in §4; state-centric Hierarchie §0 + §4.7; asymptotische Chiralität in §4.1; erster Belastungstest harter Kern in §4.3; Literaturpositionierung in §4.5
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
- `eabc_witness_54044321_verify.py` — Verifikation Zeuge $p=54\,044\,321$ (§4.4, mod-$60060$-Kanal)
- `eabc_quadruplets_fit_alpha.py` — Stufen 0–3: Diagnostik (nicht Theorie); $\alpha_E$-Plateau, H₀a/H₀b–H₃
- `eabc_quadruplets_plot.py` — Vierfeld-Diagnose-Plot ($W_E$, $R_{1/2}$, $\alpha_{\mathrm{loc}}$, $R_\beta$)
- `eabc_occupancy_tree.py` — monoidale Vierlings-Besetzungszustände $Z=(O,T,n)$, Williams-Baum (§4.6–4.7)
- `collatz_eabc_core/CollatzEabc/OccupancyMonoid.lean` — Lean: kommutatives Monoid $(Z,\oplus,Z_0)$, `blockScan_append` (§4.6–4.7)
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

## 4. Einleitung: Paradigmenwechsel und orientierter Zirkulationsfehler

$$\boxed{\;\text{Der größte Fortschritt ist nicht die Einführung neuer Größen, sondern die Entkopplung des Programms von }\Phi_E.\;}$$

Ein negatives $\Phi_E$-Ergebnis ist **Antwort auf eine von mehreren Fragen**, nicht Projekt-Zusammenbruch.

### Forschungsprogramm vs. Theorie

**Gegenwärtiger Status:** EABC ist wissenschaftstheoretisch eher **Forschungsprogramm** als bereits **Theorie**.

| Begriff | Inhalt |
|---------|--------|
| **Theorie** | erklärt einen zusammenhängenden Phänomenbereich durch einen stabilen Prinzipiensatz |
| **Forschungsprogramm** | definiert Gegenstandsbereich, Observablen, Invarianten, Hierarchie offener Fragen |

Genau das ist hier passiert — explizit in diesem Abschnitt (neben/ergänzend zur Lakatos-Einordnung unten). Das Programm legt fest, **was** untersucht wird ($G_E$, $D_E$, Skalierungsobservablen), **welche** Fragen in welcher Reihenfolge sinnvoll sind, und **welche** Aussagen den harten Kern nicht berühren dürfen.

### Konservative Paper-Einleitung

> Wir betrachten einen gerichteten Restklassengraphen $G_E$, der aus den vier nichttrivialen Restklassen modulo $12$ besteht. Aus der durch die Primzahlfolge induzierten Orientierung ergeben sich Zykluszählfunktionen $N_+(X)$ und $N_-(X)$ sowie die orientierte Zyklusdifferenz $D_E(X)=N_+(X)-N_-(X)$. Ziel der vorliegenden Untersuchung ist die Analyse der Größenordnung und Skalierung dieser Differenz. Hierzu werden normierte Observablen $R_\beta(X)$, lokale Skalierungsgrößen $\alpha_{\mathrm{eff}}(X)$ sowie die normierte Orientierung $W_E(X)$ betrachtet. Die Frage nach der Existenz eines asymptotischen Orientierungsparameters $\Phi_E=\lim W_E(X)$ wird als mögliche Endstufe einer allgemeineren Untersuchung der Skalierungsstruktur von $D_E(X)$ verstanden und **nicht** als Ausgangsannahme vorausgesetzt.

**Referee-Perspektive:** „Ich glaube noch nicht $\Phi_E\neq 0$ — aber $D_E(X)$ ist definiert und ihre Größenordnung untersuchbar." Die methodische Verschiebung geht von „Ist die Vermutung wahr?" zu „Welche Eigenschaften besitzt die Observable?"

$$\boxed{\;\text{Das EABC-Programm untersucht primär die Wachstumsordnung der orientierten Zyklusdifferenz }D_E(X)\text{, während die Holonomie }\Phi_E\text{ als mögliche Endstufe dieser Skalierungstheorie erscheint.}\;}$$

**Programm (primär):** Fehlertermtheorie des **orientierten Zirkulationsfehlers** $D_E(X)=A(X)-C(X)$ auf dem EABC-Kreisgraphen — Prime-Race-/Fehlertermstruktur zwischen den gegenläufigen Zyklusorientierungen $\gamma^+$ (ABCEA) und $\gamma^-$ (CEABC). Die Holonomiefrage ($\Phi_E\neq 0$) ist **Ebene 4** (§4.2), nicht Eingangsthese.

**Paradigmenwechsel** (Verschiebung des mathematischen Objekts, nicht nur methodische Verbesserung):

| Früher | Jetzt |
|--------|-------|
| **Grenzwerttheorie** $W_E\to\Phi_E\neq 0$ (Ergodentheorie / Dichtetheorie) | **Fehlertermtheorie** $D_E(X)=N_+(X)-N_-(X)$ (analytische Zahlentheorie) |
| Eingangsfrage: existiert nichttrivialer Grenzwert? | Eingangsfrage: welche Wachstumsordnung hat $|D_E|$ relativ zu $Q$? |

Die Verschiebung $\Phi_E \to D_E$ ändert das **mathematische Objekt** des Programms: nicht der normierte Quotient steht im Zentrum, sondern die absolute orientierte Zyklusdifferenz als primäre arithmetische Observable. Der Forschungsgegenstand ist von einer einzelnen Vermutung entkoppelt.

### Lakatos-Einordnung und Programmversionen

Im Sinne eines **Lakatos'schen Forschungsprogramms** (harter Kern, Schutzmantel, progressive vs. degenerative Verschiebungen):

| Rolle | Objekt | Ebene |
|-------|--------|-------|
| **Harter Kern** | gerichteter EABC-Kreisgraph $G_E$; orientierte Zyklusdifferenz $D_E(X)=N_+(X)-N_-(X)$; daraus induzierte Skalierungsobservablen | 0–1 |
| **Primäre Theorie** | $D_E(X)$, $Q(X)=N_+(X)+N_-(X)$ | 1 |
| **Sekundäre Schutzmantel** | $R_\beta$, $\alpha_{\mathrm{loc}}$, $\alpha_{\mathrm{eff}}$, $\alpha_E$ | 2 |
| **Orientierung** | $W_E(X)=D_E(X)/Q(X)$ | 3 |
| **Endfrage** | $\Phi_E=\lim_{X\to\infty} W_E(X)$ | 4 |

Der **harte Kern** ist nicht mehr „Holonomie", sondern: (1) der gerichtete EABC-Kreisgraph $G_E$, (2) die orientierte Zyklusdifferenz $D_E(X)$, (3) die daraus induzierte Familie von Skalierungsobservablen. Erst darauf bauen alle weiteren Fragen.

$$\boxed{\;\text{Das Programm hängt nicht mehr an einer Vermutung.}\;}$$

**Programmversionen** (Paradigmenwechsel in epistemischer Lesart):

| | **V1** (früher) | **V2** (jetzt) |
|---|-----------------|----------------|
| Eingang | $\Phi_E\neq 0$ als Leitvermutung (Grenzwerttheorie) | $D_E(X)$, $|D_E|$ relativ zu $Q$ als Leitfrage (Fehlertermtheorie) |
| Hilfsgrößen | $W_E$, Skalierung nachträglich zur Stützung von $\Phi_E$ | $R_\beta$, $\alpha_{\mathrm{loc}}$, $\alpha_E$ dienen der $D_E$-Theorie |
| $\Phi_E$ | Zielobjekt / Eingangsfrage | **Endfrage** der $D_E$-Skalierungstheorie (falls $\alpha_E=1$ und der Grenzwert existiert) |

**Stabilisierungskette** (Vorwärtsführung; $\Rightarrow$ für bewiesene Implikationen, $\leadsto$ für offene Fragen):

$$G_E \;\leadsto\; (D_E,\,Q_E) \;\leadsto\; (R_\beta,\,\alpha_{\mathrm{eff}},\,\alpha_{\mathrm{loc}}) \;\leadsto\; (\alpha_E,\,W_E) \;\leadsto\; \Phi_E.$$

Selbst $\Phi_E=0$ bleibt **wissenschaftlich interessant**: offen bleiben $\alpha_E=\tfrac{1}{2}$?, $\alpha_E>\tfrac{1}{2}$?, $\alpha_{\mathrm{loc}}$-Plateaus?, kritische $R_\beta$ — als eigenständige Prime-Race-/Fehlertermstruktur auf dem EABC-Zyklus (vgl. §4.2).

**Hauptterm-Vermutung.** Asymptotische Symmetrie der gegenläufigen Zyklusorientierungen:
$$N_+(X) \sim N_-(X)\qquad (X\to\infty),$$
und damit
$$S_E(X) \to 0.$$

**Sekundäre Frage (Ebene 4).** Normalisierter Magnetfluss / arithmetische Orientierungsklasse (`collatz_eabc_diskrete_geometrie.md` §2):
$$\Phi_E := \lim_{X\to\infty} W_E(X) = \lim_{X\to\infty}\frac{C_E(X)}{N_+(X)+N_-(X)} \;\stackrel{?}{\neq}\; 0.$$

$$\boxed{\;\text{Ob die Wachstumsordnung von }D_E\text{ bis zur linearen Skala }(|D_E|\asymp Q)\text{ reicht und damit eine nichtverschwindende Holonomie }\Phi_E\text{ erzeugt — sekundäre Frage (Ebene 4).}\;}$$

| Grenzfall | Folgerung |
|-----------|-----------|
| $\lim S_E = 0$ | keine bevorzugte globale Orientierung (Hauptterm) |
| $\lim S_E \neq 0$ | nichttriviale **arithmetische Orientierungsklasse** (H₃) |

**Numerik:** `flux_density_limit` in `collatz_eabc_hodge_eabc.py`.

**Label:** Hauptterm-Vermutung = **Vermutung**; Orientierungsfrage (H₃) = **Vermutung** / **Forschungsfrage** (Ebene 4).

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

### 4.2 Fünf-Ebenen-Architektur und Hypothesenhierarchie (H₀a/H₀b–H₃)

$$\boxed{\;\text{Nicht }\Phi_E\text{ ist der Anfang, sondern }D_E(X).\;}$$

**Vorwärtskette** (bewiesene Implikationen: $\Rightarrow$; offene Forschungsfragen: $\leadsto$):
$$G_E \;\leadsto\; D_E(X) \;\leadsto\; \alpha_E \;\leadsto\; W_E(X) \;\leadsto\; \Phi_E.$$
$$\Phi_E \neq 0 \;\Longrightarrow\; D_E(X)\sim\Phi_E\,Q(X) \;\Longrightarrow\; \alpha_E=1 \quad\text{(nur diese Richtung; keine Äquivalenz).}$$

$$\boxed{\;\begin{array}{c}\textbf{Ebene 0 — Geometrie:}\; G_E=(E,A,B,C),\;\gamma^+=\mathrm{ABCEA},\;\gamma^-=\mathrm{CEABC}\\[4pt]\textbf{Ebene 1 — Zirkulationsfehler:}\; D_E(X)=N_+(X)-N_-(X)\;\text{(primäre arithmetische Observable)}\\[4pt]\textbf{Ebene 2 — Skalierung:}\; R_\beta(X)=D_E(X)/Q(X)^\beta,\;\alpha_{\mathrm{loc}}(X)=\Delta\log|D_E|/\Delta\log Q,\;\alpha_E=\lim\alpha_{\mathrm{loc}}(X)?\\[4pt]\textbf{Ebene 3 — Orientierung:}\; W_E(X)=D_E(X)/Q(X)\\[4pt]\textbf{Ebene 4 — Holonomie:}\; \Phi_E=\lim_{X\to\infty} W_E(X)\end{array}\;}$$

**Referee-sicher:** Scheitern von $\Phi_E\neq 0$ zerstört das Programm nicht. Selbst bei $\Phi_E=0$ bleiben eigenständige Probleme: $\alpha_E=\tfrac{1}{2}$? $\alpha_E>\tfrac{1}{2}$? $\alpha_{\mathrm{loc}}$-Plateaus? kritische $R_\beta$ — als Prime-Race-/Fehlertermstruktur auf dem EABC-Zyklus.

**Methodische Leitlinie (Tao-Stil):** Fünf **Ebenen** — Geometrie → Fehlerterm → Skalierung → Orientierung → Holonomie. Keine Exponentenannahme vor Ebene 2; Holonomie-Hypothese (Ebene 4) **am Ende**, nicht am Anfang. Die Hypothesen H₀a–H₃ referenzieren jeweils eine Ebene; sie ersetzen sie nicht.

#### Ebene 0 — Geometrie

**Kreisgraph** $G_E=(V,E)$ mit $V=\{E,A,B,C\}$ (Restklassen $1,5,7,11$ mod $12$). Zwei **orientierte Zyklen** in $H_1(C_4,\mathbb{Z})$:

| Orientierung | Wort | Vorzeichen |
|--------------|------|------------|
| $\gamma^+$ | ABCEA | $+1$ |
| $\gamma^-$ | CEABC | $-1$ |

Beide tragen dasselbe Lückenmuster $(2,4,2,4)$; sie unterscheiden sich durch **zyklische Verschiebung** der Orientierung auf demselben $C_4$-Gerüst. Vollständige Definitionen: §1–2, `collatz_eabc_diskrete_geometrie.md`.

**Label:** **Definition** (Ebene 0).

#### Ebene 1 — Zirkulationsfehler (primäre Observable)

$$A(X) := N_+(X)\quad(\mathrm{ABCEA}),\qquad C(X) := N_-(X)\quad(\mathrm{CEABC}),$$
$$D_E(X) := A(X)-C(X) = N_+(X)-N_-(X),\qquad Q(X) := A(X)+C(X).$$

$D_E(X)$ ist die **primäre arithmetische Observable** des EABC-Programms — absolute orientierte Zyklusdifferenz, nicht normierter Quotient. Prime-Race-Analogie: wie $\Delta(x)=\pi(x;a,q)-\pi(x;b,q)$ fragt man zuerst nach der **Größenordnung** von $|D_E|$.

**Label:** **Definition** (Ebene 1).

#### Ebene 2 — Skalierung

$$R_\beta(X) := \frac{D_E(X)}{Q(X)^\beta}\quad(\beta\in\mathbb{R}_+).$$

Man untersucht $R_\beta(X)$ **ohne** die Existenz eines Exponenten vorauszusetzen: bleibt $R_{1/2}$ beschränkt? wächst $R_{2/3}$? Kein Zufallsmodell, keine Vermutung auf dieser Ebene.

**Lokale Skalierung** (numerisch primär):
$$\alpha_{\mathrm{loc}}(X_i) := \frac{\log|D_E(X_{i+1})| - \log|D_E(X_i)|}{\log Q(X_{i+1}) - \log Q(X_i)}\quad\text{zwischen Checkpoints }X_i, X_{i+1}.$$

**Grobe Gesamtindikation** (sekundär):
$$\alpha_{\mathrm{eff}}(X) := \frac{\log|D_E(X)|}{\log Q(X)}\quad (D_E\neq 0,\; Q>1).$$

**$\alpha_{\mathrm{loc}}$ ist für die Numerik wichtiger:** Drift, Übergangsbereiche und transiente Plateaus werden sichtbar, die $\alpha_{\mathrm{eff}}(X)$ als punktuelle Gesamtindikation überdeckt. **Wichtiger Vorbehalt:** $\alpha_{\mathrm{loc}}$ **muss nicht konvergieren** — in der analytischen Zahlentheorie kann ein Exponent über lange Bereiche stabil erscheinen und später driften oder oszillieren.

**Asymptotischer Exponent** (Vermutung, erst wenn Daten es nahelegen):
$$\boxed{\;\text{Existiert }\alpha_E\text{ mit }|D_E(X)|=Q(X)^{\alpha_E+o(1)}\text{?}\;}$$
$$\alpha_E := \inf\{\beta\in\mathbb{R}_+ : R_\beta(X)\text{ bleibt beschränkt für }X\to\infty\}.$$

Formal: $\alpha_E=\lim_{X\to\infty}\alpha_{\mathrm{loc}}(X)$ **falls** der Limes existiert — offene Frage.

| Begriff | Rolle | Ebene |
|---------|-------|-------|
| **$D_E$, $Q$** | primäre Zählgrößen | 1 |
| **$R_\beta$** | Skalierungs-Observable | 2 |
| **$\alpha_{\mathrm{loc}}$** | lokale Skalierung (Numerik primär) | 2 |
| **$\alpha_{\mathrm{eff}}$** | grobe Gesamtindikation | 2 |
| **$\alpha_E$** | kritischer Exponent | 2 |

**Kein** $Z_E=D/\sqrt{Q}$ auf der Definitionsebene — nur $D_E$, $Q$, $R_\beta$. Spezialfälle: $R_1(X)=W_E(X)$; $R_{1/2}$, $R_{2/3}$, $R_{3/4}$, $R_{0.9}$ sind **diagnostische** Normierungen (Ebene 2), keine Kerngrößen.

**Heuristik** (Diskussion, nach Ebene 2): $Z_E := R_{1/2}(X)=D_E(X)/\sqrt{Q(X)}$ — in CSV/Plot als Spalte `R_1_2` bzw. Alias $Z_E$ (**Heuristik/Diagnose**). $\alpha_E=\tfrac{1}{2}$ entspräche einem **naiven Zufallsmodell** ($|D_E|\approx\sqrt{Q}$) — **Heuristik**, keine bewiesene Referenz.

**Testarchitektur** (Ebene 2):

| Observable | $\beta$ | Lesart |
|------------|---------|--------|
| $R_{1/2}$ | $\tfrac{1}{2}$ | Wurzelrauschen-Heuristik (H₀a/H₁) |
| $R_{2/3}$, $R_{3/4}$, $R_{0.9}$ | Zwischenwerte | sublinearer Bias (H₂-Diagnostik) |
| $R_1=W_E$ | $1$ | Orientierung (Ebene 3; H₀b/H₃) |

**Drei Skalierungsfälle** (an $\alpha_E$ gebunden, **unter H₂**):

| Fall | $\alpha_E$ | $R_{1/2}$ | $W_E$ | Lesart |
|------|------------|-----------|-------|--------|
| 1. Rauschen | $\le \tfrac{1}{2}$ | $O(1)$ | $\to 0$ | Wurzelrauschen (H₀a) |
| 2. Sublinearer Bias | $\tfrac{1}{2} < \alpha_E < 1$ | $\to\infty$ | $\to 0$ | arithmetischer Bias (H₁–H₂, nicht H₃) |
| 3. Starke Holonomie | $= 1$ | $\to\infty$ (langsamer als Fall 2) | $\to \Phi_E \neq 0$ | Holonomie-Grenzfall (H₃) |

Unter H₂ folgen formal: $W_E(X) \sim Q(X)^{\alpha_E-1}$, $R_{1/2}(X) \sim Q(X)^{\alpha_E-1/2}$, $R_\beta(X) \sim Q(X)^{\alpha_E-\beta}$.

**Wichtig:** Auch bei $\Phi_E=0$ kann $\alpha_E>\tfrac{1}{2}$ gelten — ein **eigenständiges arithmetisches Phänomen**, nicht äquivalent zur Holonomiefrage.

**Label:** $\alpha_E$ = **Definition** (Ebene 2); H₂ = **Hypothese**; $\alpha_{\mathrm{loc}}$-Diagnose / $\alpha_E$-Schätzung = **Experiment**.

#### Ebene 3 — Orientierung

$$W_E(X) := \frac{D_E(X)}{Q(X)} = \frac{N_+(X)-N_-(X)}{N_+(X)+N_-(X)}.$$

$W_E(X)$ ist **Definition** — kein Grenzwert, kein Exponent. Es misst die **normierte Orientierung** — unabhängig von jeder Skalierungsexponentenannahme.

$$\boxed{\;\text{Orientierung: }W_E(X)\to 0\quad\text{oder}\quad W_E(X)\to\Phi_E\neq 0\text{?}\;}$$

- **H₀b** (keine Orientierung): $W_E(X)\to 0$ — analytische Nullhypothese.
- **H₃-Vorstufe:** $W_E(X)$ stabilisiert sich gegen einen **nichtverschwindenden** Grenzwert — noch ohne Existenzbehauptung für $\Phi_E$.

**Label:** **Definition** (Ebene 3); H₀b = **Nullhypothese**.

#### Ebene 4 — Holonomie

$$\Phi_E := \lim_{X\to\infty} W_E(X)\quad\text{(falls der Grenzwert existiert).}$$

Holonomie **am Ende**, nicht am Anfang. Entspricht Lean-**RED** `HasNonzeroHolonomyLimit`, `EABC_holonomy_limit_conjecture`.

- **H₃** (starke Holonomie): $W_E(X)\to\Phi_E\neq 0$ — Vermutung.

**Label:** **Vermutung** (Ebene 4; Lean-RED).

#### Implikationskette (Referee-präzise)

$$\boxed{\;\Phi_E\neq 0 \;\Rightarrow\; D_E(X)\sim\Phi_E\,Q(X) \;\Rightarrow\; \alpha_E=1.\;}$$

**Nicht umgekehrt:** $\alpha_E=1$ impliziert **keinen** Grenzwert von $W_E$.

**Gegenbeispiel:** $D_E(X)=Q(X)\sin(\log\log Q(X))$ — $\alpha_E=1$, aber $W_E(X)$ **ohne** Grenzwert.

| Parameter | Stärke | Bedeutung |
|-----------|--------|-----------|
| **$\Phi_E$** | Orientierungsparameter (stärker) | normierter Grenzwert $W_E$ |
| **$\alpha_E$** | Skalenparameter (schwächer) | Größenordnung von $|D_E|$ relativ zu $Q$ |

Weitere Vorwärtsimplikationen (H₀a/H₀b–H₃):

- **H₃ $\Rightarrow$ H₂** (mit $\alpha_E=1$), aber **H₂ $\nRightarrow$ H₃** (jedes $\alpha_E\in(\tfrac{1}{2},1)$ lässt $W_E\to 0$).
- **H₂ $\Rightarrow$ H₁** auf hinreichend großen Skalen ($R_{1/2}\to\infty$), aber **H₁ $\nRightarrow$ H₂** (transienter Bias reicht nicht).
- **H₀b $\nRightarrow$ H₀a** und **H₀a $\nRightarrow$ H₀b** — analytische Orientierung und Wurzelnormierung sind **getrennte** Fragen.
- **H₂ mit $\alpha_E>\tfrac{1}{2}$ $\nRightarrow$ H₃** — auch $\Phi_E=0$ ist mit $\alpha_E>\tfrac{1}{2}$ vereinbar.

#### Boxed Kernfragen (entschärft)

1. **Ebene 1/2:** Welche Wachstumsordnung hat $D_E(X)$? Bleibt $R_{1/2}$ beschränkt? Was zeigt $\alpha_{\mathrm{loc}}$? — **$D_E$ zuerst.**
2. **Ebene 2:** Existiert $\alpha_E$ mit $|D_E(X)|=Q(X)^{\alpha_E+o(1)}$? — **nur wenn Daten es nahelegen.**
3. **Ebene 3/4:** Reicht die Wachstumsordnung bis zur linearen Skala und folgt $W_E(X)\to\Phi_E\neq 0$? — **Holonomie danach.**

#### Hypothesen (H₀a/H₀b–H₃)

Schrittweise Eskalation entlang der Ebenen: **1/2** Fehlerterm/Skalierung → **3** Orientierung → **4** Holonomie.

| Hypothese | Aussage | Ebene | Label |
|-----------|---------|-------|-------|
| **H₀a** (Wurzelrauschen) | $\alpha_E \le \tfrac{1}{2}$ — $R_{1/2}(X)$ **bleibt beschränkt** | 2 | **Heuristik** (naives Zufallsmodell) |
| **H₁** (mehr als Wurzelrauschen) | $\alpha_E > \tfrac{1}{2}$ empirisch — $R_{1/2}(X)\to\infty$ | 2 | **Experiment** |
| **H₂** (asymptotischer Bias) | $|D_E(X)| \asymp Q(X)^{\alpha_E}$ mit $\alpha_E > \tfrac{1}{2}$ asymptotisch | 2 | **Hypothese** (arithmetisch stark) |
| **H₀b** (keine Orientierung) | $W_E(X)\to 0$ — analytische Nullhypothese; **H₀b $\nRightarrow$ $\alpha_E\le\tfrac{1}{2}$** | 3 | **Nullhypothese** (analytisch) |
| **H₃** (starke Holonomie) | $W_E(X)\to\Phi_E\neq 0$ — impliziert $D_E(X)\sim\Phi_E Q(X)$, $\alpha_E=1$ | 4 | **Vermutung** (Lean-**RED**: `HasNonzeroHolonomyLimit`) |

#### Prime-Race-Analogie ($\pi(x)$–$\mathrm{Li}(x)$)

Klassisch trennt man beim Primzahlgesetz **Hauptterm** und **Fehlerterm**:

| Klassisch | EABC-Analog |
|-----------|-------------|
| Hauptterm $\mathrm{Li}(x)\sim x/\log x$ | Symmetrie $N_+(X)\sim N_-(X)$, $Q(X)\to\infty$ |
| Fehlerterm $\pi(x)-\mathrm{Li}(x)$ (Größenordnung zuerst) | $D_E(X)=N_+(X)-N_-(X)$ (absolute Zyklusdifferenz) |
| Prime Race $\Delta(x)=\pi(x;a,q)-\pi(x;b,q)$ | orientierter Race $\gamma^+$ vs. $\gamma^-$ auf demselben $C_4$-Zyklus |
| Quotient $\Delta(x)/(\pi(x;a,q)+\pi(x;b,q))$ (Grenzwert sekundär) | $W_E(X)=D_E/Q$, $\Phi_E=\lim W_E$ (Ebene 3/4) |

Wie bei $\pi(x)-\mathrm{Li}(x)$ und bei $\Delta(x)=\pi(x;a,q)-\pi(x;b,q)$ fragt man zuerst nach der **Größenordnung** von $|D_E|$ relativ zu $Q$ — nicht nach einem Grenzwert des Quotienten. $D_E(X)$ ist ein orientierter Prime Race auf demselben EABC-Zyklus. Die asymptotische Orientierungsfrage (Ebene 4: $W_E\to 0$ oder $\Phi_E\neq 0$) ist von der Skalierungsfrage (Ebene 2: $\alpha_E$ via $R_\beta$) **methodisch getrennt** — analog zur Trennung von $\pi(x)-\mathrm{Li}(x)$ und mod-$q$-Races.

**Numerik und Diagnose-Plot.** `eabc_quadruplets_1e10.py` erzeugt Checkpoints mit $Q$, $D_E$, $W_E$, $R_\beta$ für $\beta\in\{\tfrac{1}{2},\tfrac{2}{3},\tfrac{3}{4},0{,}9,1\}$ (CSV-Spalten `R_1_2`, `R_2_3`, `R_3_4`, `R_9_10`, `R_1`; Alias `Z_E`=`R_1_2` nur Diagnose), $\alpha_{\mathrm{eff}}$ (`eabc_quadruplets.csv`). `eabc_quadruplets_fit_alpha.py` strukturiert die Ausgabe nach **Ebenen 0–4**; schätzt $\alpha_E$ heuristisch aus $R_\beta$-Plateaus (**Experiment**, kein Theorem). `eabc_quadruplets_plot.py` (oder `fit_alpha --plot`) erzeugt ein Vierfeld-Diagramm (`eabc_quadruplets_diagnose.png`):

| Panel | Kurve | Ebene / Lesart |
|-------|-------|----------------|
| 1 | $R_{1/2}(X)$ ($Z_E$ Alias) | Ebene 2: Skalierung (H₀a / H₁) |
| 2 | $\alpha_{\mathrm{loc}}(X)$ | Ebene 2: $\alpha_{\mathrm{loc}}$-Diagnose (primär; kein Satz) |
| 3 | $R_\beta$ für $\beta\in\{\tfrac{1}{2},\tfrac{2}{3},\tfrac{3}{4},0{,}9,1\}$ | Ebene 2: Skalierungsdiagnostik |
| 4 | $W_E(X)=R_1(X)$ | Ebene 3/4: H₀b / H₃ (Orientierung / Holonomie) |

**Label:** Ebene 0/1 = **Definition**; Ebene 2 / H₂ = **Hypothese**; $\alpha_{\mathrm{loc}}$-Diagnose = **Experiment**; Ebene 4 / H₃ = **Vermutung** (Lean-RED); H₀a = **Heuristik**; H₀b = **Nullhypothese**.

### 4.3 Numerischer Stand und erster Belastungstest (bis $10^{10}$)

$$\boxed{\;\text{Zeigt }D_E(X)\text{ überhaupt ein reproduzierbares, nichttriviales Wachstum?}\;}$$

**Scope des Tests:** nur $(G_E,\,N_+,\,N_-,\,D_E)$ und die daraus abgeleiteten Skalierungsobservablen — **kein** $\Phi_E$, **keine** Holonomie. Der Belastungstest prüft den harten Kern (Ebene 0–2), nicht die Endfrage Ebene 4.

**Methodische Leitlinie (Tao-Stil):** Alle folgenden Befunde sind **Experimente** an endlichem $X$, keine Theoreme und keine Vermutungsbestätigungen. Die **Nullhypothese H₀a** lautet $|D_E(X)|=O(\sqrt{Q(X)})$, d. h. $R_{1/2}(X)=D_E(X)/\sqrt{Q(X)}$ bleibt beschränkt.

**Pipeline** (Vollauf bis $X=10^{10}$):

```bash
python3 eabc_quadruplets_1e10.py --X 10000000000
python3 eabc_quadruplets_fit_alpha.py eabc_quadruplets.csv --plot --plot-loglog
python3 eabc_quadruplets_plot.py eabc_quadruplets.csv --plot-loglog
```

**Checkpoint-Strategie:** geometrische Serie $10^e$, $2\times 10^e$, $5\times 10^e$ für $e=6,\ldots,10$ — **13 Checkpoints** in `eabc_quadruplets.csv`.

#### Boxed Hauptaussagen (Belastungstest $10^{10}$, Experiment)

$$\boxed{\;\text{Bis }10^{10}\text{ gibt es keinen numerischen Hinweis auf Holonomie }(\Phi_E\neq 0).\;}$$

$$\boxed{\;\text{Bis }10^{10}\text{ gibt es auch keinen klaren Hinweis auf einen sublinearen Bias }(\alpha_E>\tfrac{1}{2}).\;}$$

$$\boxed{\;\text{Bis }10^{10}\text{ verhält sich }D_E(X)\text{ vollständig kompatibel mit einem Fehlerterm auf Wurzelskala }(|D_E|=O(\sqrt{Q})).\;}$$

#### Datenstand (13 Checkpoints, Experiment)

**Kompakte Referenztabelle** (Schlüsselgrößen; vollständige Spalten in `eabc_quadruplets.csv`):

| $X$ | $Q$ | $D_E$ | $W_E$ | $R_{1/2}$ | $W_{420}$ |
|-----|-----|-------|-------|-----------|-----------|
| $10^6$ | 166 | $+2$ | 0,0120 | 0,15 | 0,218 |
| $2\times 10^6$ | 295 | $+9$ | 0,0305 | 0,52 | 0,184 |
| $5\times 10^6$ | 546 | $+24$ | 0,0440 | 1,03 | 0,220 |
| $10^7$ | 899 | $+1$ | 0,0011 | 0,03 | 0,140 |
| $2\times 10^7$ | 1468 | $-26$ | $-0,0177$ | $-0,68$ | 0,196 |
| $5\times 10^7$ | 2847 | $+27$ | 0,0095 | 0,51 | 0,129 |
| $10^8$ | 4768 | $+48$ | 0,0101 | 0,70 | 0,083 |
| $2\times 10^8$ | 8097 | $+63$ | 0,0078 | 0,70 | 0,048 |
| $5\times 10^8$ | 16331 | $+137$ | 0,0084 | 1,07 | 0,043 |
| $10^9$ | 28388 | $+52$ | 0,0018 | 0,31 | 0,026 |
| $2\times 10^9$ | 49262 | $-16$ | $-0,0003$ | $-0,07$ | 0,017 |
| $5\times 10^9$ | 102950 | $-58$ | $-0,0006$ | $-0,18$ | 0,008 |
| $10^{10}$ | 180529 | $+25$ | 0,00014 | 0,06 | 0,0065 |

**$D_E$:** schwaches Wachstum der absoluten Differenz; bei $X=10^{10}$ ist $D_E=25$ bei $Q=180529$ — nichttrivial ($D_E\neq 0$), aber weit unter linearer Skala ($|D_E|\ll Q$).

**$W_E$-Verlauf** (Orientierung, Ebene 3): $|W_E(10^6)|\approx 1{,}20\times 10^{-2}$, $|W_E(10^8)|\approx 1{,}01\times 10^{-2}$, $|W_E(10^{10})|\approx 1{,}38\times 10^{-4}$ — **Vorzeichenwechsel** ($D_E$ wechselt das Vorzeichen bei $2\times 10^7$, $2\times 10^9$, $5\times 10^9$), kein Stabilisierungstrend. Das spricht für **H₀b** ($W_E\to 0$), nicht für H₃.

**$R_{1/2}$-Folge** (13 Werte, gerundet):
$$0{,}15,\;0{,}52,\;1{,}03,\;0{,}03,\;-0{,}68,\;0{,}51,\;0{,}70,\;0{,}70,\;1{,}07,\;0{,}31,\;-0{,}07,\;-0{,}18,\;0{,}06.$$
Keine monotone Drift — Schwankungen der Größenordnung $O(1)$, kompatibel mit **H₀a** (Wurzelrauschen). Die erste große Überraschung wäre **nicht** $\Phi_E\neq 0$, sondern ein **reproduzierbarer Drift** von $R_{1/2}$ weg von $O(1)$; bis $10^{10}$ nicht erkennbar.

**$\widehat{\alpha}_E\approx 0{,}3645$** (Log-log-Fit über alle Checkpoints): **trügerisch** — dominiert von $\log|D_E|$-Sprüngen bei kleinen $|D_E|$ (insbesondere $D_E=+1$ bei $10^7$). Die $\alpha_{\mathrm{loc}}$-Liste zeigt Rauschen ohne stabilisierenden Exponenten; kein empirischer Beleg für $\alpha_E>\tfrac{1}{2}$.

**$W_{420}$ (Nebenkanal):** Spannweite der sechs regulären mod-$420$-Klassen, normiert — fast **monoton fallend** ($0{,}218\to 0{,}0065$). Systematischer als $W_E$, aber **separate Diagnostik** (mod-$420$-Prime-Race, nicht ABCEA/CEABC-Zyklus); nur **Experiment**, kein Theorem.

#### Referee-Diagnose (6 Punkte, Experiment)

| # | Befund | Label |
|---|--------|-------|
| 1 | **Gesichert:** $D_E(X)=N_+(X)-N_-(X)$ ist wohldefiniert und auf endlichen Skalen berechenbar | **Definition** |
| 2 | **Numerisch:** $D_E\neq 0$ auf allen 13 Checkpoints (Vorzeichenwechsel inklusive) | **Experiment** |
| 3 | **Nicht sichtbar:** $\alpha_E>\tfrac{1}{2}$ — $R_{1/2}$ bleibt $O(1)$, kein Drift | **Experiment** (H₀a kompatibel) |
| 4 | **Nicht sichtbar:** $\Phi_E\neq 0$ — $W_E$ instabil, $|W_E|\to 0$ auf großen Skalen | **Experiment** (H₀b kompatibel) |
| 5 | **Hypothesen-Ranking:** H₀a und H₀b haben **stärkere numerische Unterstützung** als H₂ (sublinearer Bias) und H₃ (Holonomie) | **Experiment** |
| 6 | **Offen:** ob $|D_E|=O(\sqrt{Q})$ asymptotisch gilt — bis $10^{10}$ weder widerlegt noch bewiesen | **Forschungsfrage** |

**Schluss (Experiment, nicht Theorem):** Der Belastungstest bis $10^{10}$ bestätigt den harten Kern ($D_E$ definiert, nichttrivial, berechenbar) und liefert **keinen** numerischen Anlass, H₂ oder H₃ gegen H₀a/H₀b zu bevorzugen. Holonomie ($\Phi_E\neq 0$) bleibt **Ebene 4** und wird erst relevant, wenn $R_{1/2}$ oder $W_E$ ein klares asymptotisches Signal zeigen — bis $10^{10}$ nicht der Fall.

**Verknüpfung:** vgl. `collatz_eabc_epistemik_schichten.md` §0 (10$^{10}$-Referenzpunkt); Diagnose-Plot: `eabc_quadruplets_diagnose.png`.

**Label:** Belastungstest = **Experiment**; H₀a/H₀b = **Heuristik/Nullhypothese**; H₂/H₃ = **Hypothese/Vermutung** (numerisch nicht gestützt bis $10^{10}$).

### 4.4 mod-$60060$-Zeuge: Primvierling $p=54\,044\,321$ (13-verfeinerte Schale)

$$\boxed{\;p=54\,044\,321,\quad p\bmod 60060=50\,381,\quad \text{Orientierung ABCE.}\;}$$

**Label:** **Beispiel/Zeuge** — struktureller Beleg für A-Startklasse und späten mod-$60060$-Kanal; **kein** Holonomie-Theorem, **kein** Beleg für $\Phi_E\neq 0$.

#### Setup

| Größe | Wert | Lesart |
|-------|------|--------|
| Primvierling | $(54\,044\,321,\;54\,044\,323,\;54\,044\,327,\;54\,044\,329)$ | alle vier Positionen prim |
| $p\bmod 12$ | $5$ | **A-Startklasse** → 4-Fenster **ABCE** ($\gamma^+$-Vierling) |
| $p\bmod 60060$ | $50\,381$ | **13-verfeinerte Schale** $60060=2\cdot3\cdot5\cdot7\cdot11\cdot13$ |
| Orientierung | ABCE | **nicht** CEAB ($p\bmod 12=11$ wäre $\gamma^-$-Vierling) |

Die mod-$420$-Diagnostik in §4.3 und `eabc_quadruplets_1e10.py` nutzt sechs reguläre HL-Kanäle. Die **13-verfeinerte Schale** mod $60060$ zerlegt jeden mod-$420$-Kanal weiter; es gibt $378$ HL-zulässige Startrestklassen mod $60060$ (Hardy–Littlewood: alle vier Positionen $p,p{+}2,p{+}6,p{+}8$ sind $\not\equiv 0$ mod jedes $q\in\{2,3,5,7,11,13\}$; vgl. `eabc_hl_coefficient_hypotheses.py`, `admissible_mod`).

#### CERN-Motiv (didaktisch, **Analogie**)

Wie eine theoretisch zugelassene, lange leere Detektorkammer, die schließlich ein Ereignis registriert: Der Kanal $r=50\,381$ ist **arithmetisch zulässig** (HL-Sieb), aber unter den **spätesten** mod-$60060$-Kanälen, die bis $X=10^8$ erstmals einen Primvierling tragen. Erst bei $p=54\,044\,321$ erscheint dort ein Vierling — der **letzte** der zehn spätesten HL-Kanäle in dieser Reichweite.

**Späteste zehn HL-Kanäle mod $60060$ bis $10^8$** (Erstauftreten $p_1(r)$ absteigend; **Experiment**, reproduzierbar):

| Rang | Kanal $r$ | Orient. | erstes $p$ |
|------|-----------|---------|------------|
| 1 | $50\,381$ | ABCE | $54\,044\,321$ |
| 2 | $34\,211$ | CEAB | $33\,067\,211$ |
| 3 | $54\,911$ | CEAB | $31\,706\,531$ |
| 4 | $46\,301$ | ABCE | $26\,953\,181$ |
| 5 | $18\,281$ | ABCE | $26\,684\,921$ |
| 6 | $36\,011$ | CEAB | $26\,282\,231$ |
| 7 | $39\,041$ | ABCE | $21\,720\,701$ |
| 8 | $18\,371$ | CEAB | $20\,799\,131$ |
| 9 | $11\,651$ | CEAB | $20\,612\,231$ |
| 10 | $2\,411$ | CEAB | $20\,062\,451$ |

Unter diesen zehn ist **nur** Kanal $50\,381$ ein **ABCE**-Start ($p\bmod 12=5$); die übrigen neun sind **CEAB** ($p\bmod 12=11$). Der Zeuge füllt damit den spätesten **A-Start**-HL-Kanal in dieser Liste — konsistent mit §2 ($\gamma^+$ vs. $\gamma^-$), ohne eine Holonomie-Aussage.

#### Verifikation

```bash
python3 eabc_witness_54044321_verify.py
python3 eabc_witness_54044321_verify.py --scan   # optional: 10-späteste-Kanäle bis 10^8
```

Lean (`PatternCount.lean`): Kommentar-Zeuge $p=54\,044\,321$ — ABCE-Vierling, mod-$60060$-Kanal $50\,381$ (**Beispiel**).

**Abgrenzung:** Dieser Abschnitt dokumentiert **Geometrie und Restklassenstruktur** (Ebene 0), nicht $D_E$-Skalierung (Ebene 1–2) und nicht $\Phi_E$ (Ebene 4). Ein spät gefüllter HL-Kanal ist **kein** Widerspruch zu H₀a/H₀b aus §4.3.

### 4.5 Literaturpositionierung und Novelty

$$\boxed{\;\text{EABC ist keine neue Zahlentheorie ex nihilo, sondern eine neue gerichtete Observable auf einem bekannten Literaturfeld.}\;}$$

**Kernthese (Positionierung).** Das EABC-Programm führt **keine** neue modulare Zahlentheorie „von Grund auf“ ein. Es definiert vielmehr eine **gerichtete Observable** auf einem etablierten Gegenstandsbereich:

| Literaturfeld | Bezug zu EABC |
|---------------|---------------|
| **Prime Races** / Chebyshev-Bias | $D_E(X)$ als orientierter Race zwischen zwei Zählfunktionen |
| **Consecutive-prime biases** mod $q$ | Primfolge-induzierte Kanten auf $(\mathbb{Z}/12\mathbb{Z})^\times$ |
| **Hardy–Littlewood $k$-Tupel** | Vierlings-/Fünferfenster, HL-zulässige Startrestklassen (vgl. §4.4) |
| **Sieve gap cycles** / primoriale Kanäle | mod-$420$-, mod-$60060$-Verfeinerung, späte Kanalfüllung |

$$\boxed{\;\text{Neu ist nicht mod }12,\text{ sondern die gerichtete EABC-Zählung.}\;}$$

#### Definition — orientierter Prime-Race-Parameter

Auf dem reduzierten Restklassengraphen
$$(\mathbb{Z}/12\mathbb{Z})^\times = \{1,5,7,11\} \;\cong\; \{E,A,B,C\}$$
mit dem elementaren Zyklus $A\to B\to C\to E\to A$ (vgl. §1–2) ist
$$D_E(X) := N_+(X) - N_-(X)$$
ein **orientierter Prime-Race-Parameter**: $N_+(X)$ zählt geschlossene 5-Fenster vom Typ **ABCEA** ($\gamma^+$), $N_-(X)$ vom Typ **CEABC** ($\gamma^-$). Die beiden Wörter sind **chiral entgegengesetzte Orientierungen** desselben $C_4$-Gerüsts — nicht zwei verschiedene Lückenmuster, sondern zyklisch verschobene Richtungen auf demselben Kreis.

**Label:** $D_E$, ABCEA vs. CEABC als chirale Gegenrichtungen = **Definition** (Ebene 0–1).

#### Literatur — Nachbarn und Lesart

Während mod-$q$-Biases aufeinanderfolgender Primzahlen und Prime Races gut untersucht sind, wird hier eine **spezifische gerichtete Observable** auf $(\mathbb{Z}/12\mathbb{Z})^\times$ eingeführt. Die beiden chiralen Fünferfenster ABCEA und CEABC definieren eine Zirkulationszählung $D_E(X)$. Diese Observable ist als **Spezialisierung** der Literatur zu consecutive-prime biases und Prime Races lesbar; in genau dieser **graphorientierten, EABC-chiralen Form** scheint sie jedoch **nicht dokumentiert** zu sein (vgl. Prior-Art-Suche unten).

| # | Referenz | arXiv / Quelle | Relevanz für EABC |
|---|----------|----------------|-------------------|
| 1 | **Lemke-Oliver & Soundararajan** | [arXiv:1603.03720](https://arxiv.org/abs/1603.03720) | Unerwartete Biases in aufeinanderfolgenden Prim-Mustern mod $q$; HL-Heuristik für $k$-Tupel |
| 2 | **Rubinstein & Sarnak** | *Experimental Mathematics* **3** (1994), 173–197 | Klassischer Prime Race / Chebyshev-Bias unter GRH/GSH; Vorlage für $D_E$ als orientierte Differenz |
| 3 | **Goldston, Graham, Pintz, Yildirim** | [arXiv:1510.00743](https://arxiv.org/abs/1510.00743) | Eratosthenes-Sieb-Lückenzyklen, primoriale Kanäle; Verfeinerungsschalen (§4.4) |
| 4 | **Pomerance, Rubinstein & Sarnak** | [arXiv:1709.06168](https://arxiv.org/abs/1709.06168) | Set von Primzahlen mit vorgegebenen Anfangslücken; HL-/Sieve-Kontext |
| 5 | **Lamzouri** | [arXiv:1101.0836](https://arxiv.org/abs/1101.0836) | Prime-Race-Oszillationen, Vorzeichenwechsel (vgl. §4.3) |
| 6 | *(weitere Prime-Race-Literatur)* | [arXiv:2405.03540](https://arxiv.org/abs/2405.03540) | Aktuelle Race-/Bias-Forschung; Einordnung als Nachbarfeld |

**Label:** Tabelle = **Literatur** (extern); Zuordnung zu $D_E$ = **Analogie** / **Spezialisierung**, kein Zitat-Beweis.

#### Prior Art — arXiv-Suche (Explizit)

Eine gezielte Suche nach **ABCEA**, **CEABC**, **EABC-orientierter Zählung** und der Kombination „mod $12$ + chiral + consecutive primes + 5-cycle“ in arXiv und einschlägiger Primzahl-Literatur lieferte **keine direkte Vorarbeit** für genau diese graphorientierte, chirale Fünferfenster-Zählung $D_E(X)$.

| Suchgegenstand | Ergebnis |
|----------------|----------|
| ABCEA / CEABC als orientierte 5-Zyklus-Wörter | **kein** direkter Treffer in etablierter Literatur |
| $D_E(X)=N_+ - N_-$ auf $G_E$ | **kein** dokumentiertes Analogon unter diesem Namen/Setup |
| mod-$12$-Restklassengraph + chirale Zirkulation | Nachbarliteratur (L&O, Rubinstein–Sarnak), aber **nicht** in EABC-Form |

$$\boxed{\;\text{Keine direkte Prior Art für die orientierte ABCEA/CEABC-Zählung }D_E(X)\text{ gefunden (arXiv-Suche).}\;}$$

Das ist **kein** Behaupten mathematischer Neuheit im Sinne eines Theorems — nur eine **ehrliche Einordnung**: Die Observable sitzt im Nachbarfeld bekannter Biases, die **spezifische gerichtete Zählung** auf dem EABC-Zyklus ist jedoch (Stand Dokumentation) **nicht zitiert**.

**Label:** Prior-Art-Statement = **Literatur-Recherche** (offen für Widerspruch durch Referee).

#### Zeuge — numerische Fallstudie (Querverweis §4.4)

Der mod-$60060$-Zeuge $p=54\,044\,321$ ($p\bmod 60060=50\,381$, $50\,381\equiv A\pmod{12}$) illustriert die HL-/Kanal-Seite der Positionierung: **spätes Erstauftreten** eines zulässigen **13-verfeinerten ABCE-Kanals** unter den spätesten HL-Kanälen bis $10^8$. Details, Verifikation und Tabelle der zehn spätesten Kanäle: **§4.4**.

| Größe | Wert | Lesart |
|-------|------|--------|
| $p$ | $54\,044\,321$ | Primvierling-Start |
| $p\bmod 60060$ | $50\,381$ | 13-verfeinerte HL-Schale |
| $p\bmod 12$ | $5$ ($A$) | $\gamma^+$-Startklasse → ABCE-Vierling |
| Rolle | spätester ABCE-Kanal (Rang 1) | **Experiment** / **Beispiel** |

Der Zeuge stützt die **geometrische** Lesart (Ebene 0) und die HL-Nachbarschaft — **nicht** $\Phi_E\neq 0$ oder $D_E$-Skalierung über $10^{10}$ (§4.3).

**Label:** Zeuge = **Experiment** / **Beispiel** (konsistent mit §4.4).

#### Collatz — defensive Abgrenzung (Programmarchitektur)

$$\boxed{\;\text{Collatz}\;\leftrightarrow\;\text{EABC}\;\leftrightarrow\;\text{Quaternion}\;}$$

ist **derzeit programminterne Heuristik/Architektur**, **nicht** etablierte Literatur. Die Verknüpfung Collatz–EABC–Quaternion dient der Forschungsorganisation (`collatz_generalangriff_2026.md`, Schicht **C** in `collatz_eabc_epistemik_schichten.md`) — sie darf **nicht** als zitierte mathematische oder physikalische These ausgegeben werden.

| Aussage | Status | Label |
|---------|--------|-------|
| $D_E(X)$, $G_E$, ABCEA/CEABC | definiert, programmintern kanonisch | **Definition** |
| Collatz $\Rightarrow$ EABC $\Rightarrow$ Quaternion | interne Brückenheuristik | **Heuristik** / **Ikone** (Schicht **C**) |
| Collatz–EABC–Quaternion in der Literatur | **nicht** etabliert | **explizit nicht behauptet** |

**Label:** Collatz-Box = **Heuristik** (defensiv); keine Literaturbehauptung.

### 4.6 Monoid der Vierlings-Besetzungszustände (Cook–Mertz / Williams)

$$\boxed{\;\Phi_E \longrightarrow D_E \longrightarrow Z\;}$$
**Perspektivwechsel:** von einer asymptotischen Observable zu einer **algebraischen Streaming-Struktur**. EABC erscheint hier nicht primär als Holonomie- oder Fehlerterm-Theorie, sondern als **kommutative Streaming-Kompression arithmetischer Muster**; $N_\pm$, $D_E$, $Q_E$ werden **Auswertungen** eines algebraischen Endzustands (vgl. unten).

---

#### Ebene 1 — Abstraktes Monoid (unabhängig von Primzahlen)

**Definition (Zustandsmonoid, abstrakt).** Sei $Z=(O,T,n)$ mit

- $O$ in einer Menge mit assoziativer, kommutativer Vereinigung $\cup$;
- $T$ in einer Menge mit idempotenter Meet-Operation $\wedge$ (hier: punktweises $\min$ auf Triggerwerten);
- $n$ in einem kommutativen Monoid $(\mathcal N,+)$ (hier: $\mathbb N$-Zähler / `Counter` pro Kanal).

Definiere
$$Z_1 \oplus Z_2 = \bigl(O_1 \cup O_2,\; T_1 \wedge T_2,\; n_1 + n_2\bigr), \qquad Z_0 = (\varnothing,\, \top,\, 0),$$
wobei $\top$ das neutrale Element der Meet-Operation ist (in der Implementierung: leeres Trigger-Dictionary als $\top$ auf dem **Support** der später besetzten Kanäle).

$$\boxed{(Z,\oplus,Z_0)\ \text{ ist ein kommutatives Monoid.}}$$

**Beweis (Komponentenweise).** Assoziativität und Kommutativität von $\oplus$ reduzieren sich auf die Standardgesetze von $\cup$, $\wedge$ und $+$:
$$(A\cup B)\cup C = A\cup(B\cup C), \quad \wedge(\wedge(a,b),c)=\wedge(a,\wedge(b,c)), \quad (n_1+n_2)+n_3=n_1+(n_2+n_3).$$

**Folgerung (Theorem).** $\bigoplus_{i=1}^{k} Z_i$ ist **wohldefiniert**. Lineare Auswertung, binäre Tree-Reduction, parallele Auswertung und MapReduce-artige Verarbeitung liefern **exakt denselben Endzustand** — mathematisch stärker als „der Python-Code funktioniert“.

**Label:** Ebene 1 = **Definition** + **Theorem** (reine Algebra).

---

#### Ebene 2 — EABC-Interpretation

Erst auf Ebene 2 wird $Z$ zum **EABC-Zustandsraum** (HL-Kanäle mod $60060$, $|C_M|=378$):

| Komponente | EABC-Lesart |
|------------|-------------|
| $O$ | Menge **offener** (besetzter) Muster / Kanäle $O\subseteq C_M$ |
| $T(c)$ | **früheste Trigger-** / Erstbesetzungsinformation (Prim $p$ auf Kanal $c$) |
| $n(c)$ | Zählvektor ($N_+$, $N_-$, Ereignisse pro Kanal/Faser — je nach Auswertung) |

**Scan:** Primstrom $\to$ lokale Blockzustände $Z_i$ via `scan_block`; **Merge:** $Z(N)=\bigoplus_i Z_i$ via `merge_state` / `reduce_tree`.

**Zwei Architekturen:**

| Bisher (Holonomie-Hierarchie) | Neu (Streaming) |
|------------------------------|-----------------|
| Primzahlen $\to N_\pm \to D_E \to W_E$ | Primstrom $\to$ lokale $Z_i$ $\to$ Monoid-Merge $\to Z_{\mathrm{final}}$ |
| asymptotische Observable zuerst | algebraischer Endzustand zuerst |

Der Endzustand $Z_{\mathrm{final}}$ enthält genügend Information, um $N_+$, $N_-$, $D_E$, $Q_E$ **zu rekonstruieren** (Auswertungsfunktional, nicht Teil des Monoid-Beweises).

#### Williams-Pipeline

$$\text{lange Prim-/Vierlingshistorie} \longrightarrow \text{lokale Blockzustände} \longrightarrow \text{monoidale Tree-Evaluation} \longrightarrow \text{kompakter EABC-Endzustand}$$

Holonomie ($\Phi_E$, Ebene 4) ist **sekundär**; sie setzt eine stabile asymptotische Lesart von Auswertungen auf $Z$ voraus (vgl. §4.1–4.3).

**Implementierung:** `eabc_occupancy_tree.py` — `identity_state` ($Z_0$), `merge_state` ($\oplus$), `state_equal`, `scan_block`, `reduce_tree`; Tests: `tests/test_eabc_occupancy_tree.py` (Assoziativität, Kommutativität).

```bash
python3 eabc_occupancy_tree.py --N 100000000 --simple
pytest tests/test_eabc_occupancy_tree.py -q
```

**Label:** Ebene 2 + Williams-Pipeline = **Definition** / **Heuristik** (algorithmische Lesart); $\rho(N)=|O(N)|/|C_M|$ bei $N=10^8$: **Experiment** (`test_scan_occupancy_reference_1e8`).

**Lean (Ebene 1, GREEN):** `collatz_eabc_core/CollatzEabc/OccupancyMonoid.lean` — `CommMonoid (OccupancyState α β)`, `foldMerge_perm`, `blockScan_append` (Streaming-Faktorisierung $F(P_1\sqcup P_2)=F(P_1)\oplus F(P_2)$ auf Blocklisten).

---

### 4.7 Zustandszentrierte Mathematik (state-centric shift)

**Shift:** objektzentriert $\longrightarrow$ **zustandszentriert**

| Alt (objektzentriert) | Neu (zustandszentriert) |
|-----------------------|-------------------------|
| „Welche Eigenschaft haben Primzahlen?" | „Welcher Zustand genügt zur Beschreibung des Prozesses?" |

**Hierarchie (kanonisch):**
$$\text{Primstrom} \longrightarrow \text{Streaming-Monoid} \longrightarrow \text{Observablen} \longrightarrow \text{Interpretation}$$
$$P \;\longrightarrow\; Z \;\longrightarrow\; (D_E,\,Q_E,\,W_E,\,\ldots) \;\longrightarrow\; \Phi_E$$

$$\boxed{\;Z\ \text{ist zentral — nicht } D_E,\ \text{nicht } \Phi_E.\;}$$

**Harter Kern (Lakatos):** Das Monoid $(Z,\oplus,Z_0)$ ist **vollständig beweisbar** — unabhängig von Primzahlen, HL, RH, Holonomie. Das ist der algebraische **Hard Core** des Programms; alles Weitere (EABC-Lesart, $D_E$, $\Phi_E$) ist Schutzgürtel oder Interpretationsschicht.

**Satz (Streaming-Faktorisierung, paper-ready).** Es existiert ein Monoid-Homomorphismus
$$F:\ \text{PrimeStream} \longrightarrow M = (Z,\oplus,Z_0)$$
(kommutatives Monoid), sodass für jede Partition $P = P_1 \sqcup \cdots \sqcup P_k$ gilt:
$$\boxed{\;F(P) = F(P_1) \oplus \cdots \oplus F(P_k).\;}$$
**Label:** **Theorem** (Ebene 1; Lean: `blockScan_append`, `foldMerge_append` in `OccupancyMonoid.lean`).

**Observablen als Funktoren** auf dem Endzustand $Z_{\mathrm{final}}$ (Ebene 2, Struktur — nicht Teil des Monoid-Beweises):
$$f_D(Z) = D_E,\qquad f_Q(Z) = Q_E,\qquad f_W(Z) = W_E.$$

$$\boxed{\;\text{Die wichtigste Aussage ist nicht }\Phi_E\neq 0,\ \text{sondern}\quad F(P_1 \sqcup P_2) = F(P_1) \oplus F(P_2).\;}$$

**Offene Forschungsfrage (Minimalität).** Ist $Z=(O,T,n)$ **minimal**? Definiere Kompressionskomplexität
$$K_E(N) = \min\{\text{Bitgröße}(S) : S \Rightarrow (D_E,\,Q_E,\,W_E)\}.$$
**Offen:** $K_E(N) = O(\log N)$? $N^\alpha$? $O(1)$?

**Williams — Vorsicht:** Williams beweist **keine** Monoide. Die Verbindung ist **Tree-Evaluation-Philosophie**: der Zustand genügt, die volle Historie nicht. Dieselbe Architektur wie stream $\to$ state $\to$ merge $\to$ query (vgl. §4.6).

**Parallele zu Streaming-Algorithmen:**

| Phase | EABC | Streaming |
|-------|------|-----------|
| Eingabe | Primstrom $P$ | Datenstrom |
| Zustand | lokale $Z_i$, Merge zu $Z_{\mathrm{final}}$ | Sketch / Aggregat |
| Antwort | $f_D(Z)$, $f_Q(Z)$, … | Query auf Endzustand |

**Label:** §4.7 = **Definition** / **Theorem** (Monoid) + **Forschungsfrage** ($K_E$) + **Analogie** (Williams, Streaming).

---

## 5. EABC-Hypothese (Fehlerterm)

**EABC-Hypothese (Fehlerterm).** Der Fehlerterm
$$D_E(X) = N_+(X) - N_-(X)$$
ist **kein** reines Rauschen. Er trägt einen **nichttrivialen Chebyshev-artigen Bias** bzw. **Oszillationsstruktur**, verknüpft mit den **Fehlertermen der Primzahlverteilung modulo $12$** (bzw. den zugehörigen Dirichlet-$L$-Funktionen auf $(\mathbb{Z}/12\mathbb{Z})^\times$).

Qualitativ: wie beim klassischen Chebyshev-Bias mod $4$ kann die **absolute Differenz** $D_E(X)$ vorzeichenbehaftet und strukturiert oszillieren, während der **normierte Hauptterm** $S_E(X)$ verschwindet.

**Zentrale Frage:** Verhält sich $\widetilde{D}_E(X)$ wie reines Rauschen, oder zeigt sie stabile Vorzeichenasymmetrie / Oszillation gegenüber Isotropie- und Shuffle-Nullmodellen?

**Label:** EABC-Fehlerterm-Hypothese = **Hypothese** (stärker als Hauptvermutung allein).

**Details und Numerik:** `collatz_eabc_fehlerterm_hypothese.md` §3–5; `collatz_eabc_holonomie_fehlerterm.py`; Ebenen 0–4 und H₀a/H₀b–H₃: §4.2; `eabc_quadruplets_1e10.py`, `eabc_quadruplets_fit_alpha.py`, `eabc_quadruplets_plot.py`.

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
