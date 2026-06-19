# EABC Epistemik: Schichten A / B / C und rote RadiationSpace-Schicht

**Status:** Kanonisches Schichten-Framework — **kein** Physikanspruch  
**Branch:** `collatz/eabc-05-holonomie-fehlerterm` (PR #59), `collatz/eabc-h12-epistemik-chiralitaet-doku` (PR #72), `collatz/eabc-h13-quadruplets-alpha` (PR #73)  
**Tao-Labels:** Definition | Theorem | Modellabbildung | Analogie | Ikone | Forschungsfrage | Scaffold

**Vier Wahrheitstypen (Lean `FlussPhiE.lean`):**

$$\boxed{\;\textbf{Theorem} \;\neq\; \textbf{Struktur} \;\neq\; \textbf{Brücke} \;\neq\; \textbf{Ikone}\;}$$

| Schicht | Label | Lean | Beispiel |
|---------|-------|------|----------|
| **A** | Theorem | bewiesen, kein `sorry` | `E_plus_union_E_minus`, `Phi_E_eq_inner_product_discrete` |
| **B** | Struktur | Definition / `Prop` | `W_E_up_to`, `phi_E_conjecture`, `HasPhi_E` |
| **R** | Forschungsbrücke | `sorry`, asymptotisch | `phi_E_conjecture_statement`, `Phi_E_zero_of_symmetry` |
| **C** | Ikone | nur Markdown | Thomson, Sagnac, RadiationSpace-Lesart |

**Querverweise:**
- `collatz_eabc_diskrete_geometrie.md` — **kanonisch:** $G_E$, $E^\pm$, $\Phi_E$, EABC-Vermutung (Schicht **B**)
- `collatz_eabc_epistemik_physik.md` — Physik vs. Modell; Wegfunktion, Nicht-SRT (Schicht **C**-Abgrenzung)
- `collatz_generalangriff_2026.md` — Gesamtarchitektur, Brücke **L**, Hurwitz-Nebenzweig (Schicht **A**)
- `collatz_mathlib_eabc_kandidaten.md` — Mathlib-Inventar; Lean `RadiationSpace` (rote Schicht)
- `collatz_eabc_euklidische_hebung.md` — algebraischer Nebenzweig $\mathbb{R}\subset\mathbb{C}\subset\mathbb{H}\subset\mathbb{O}$ (Schicht **A**)
- `collatz_eabc_core/CollatzEabc/RadiationSpace.lean` — interpretativer Lean-Scaffold (rot)
- `collatz_eabc_core/CollatzEabc/OccupancyMonoid.lean` — Streaming-Monoid $(Z,\oplus,Z_0)$, `blockScan_append` (Schicht **A**, §4.6–4.7)
- `collatz_eabc_zirkulationshypothese.md` — Monoid §4.6; state-centric Hierarchie §4.7; Level-2 §4.8; Spin-Liquid-Analogie (methodisch, nicht ontologisch) §4.8.1

---

## 0. Forschungsprogramm vs. Theorie

**Gegenwärtiger Status:** Das EABC-Programm ist wissenschaftstheoretisch eher **Forschungsprogramm** als bereits **Theorie**.

| Begriff | Inhalt |
|---------|--------|
| **Theorie** | erklärt einen zusammenhängenden Phänomenbereich durch einen stabilen Prinzipiensatz |
| **Forschungsprogramm** | definiert Gegenstandsbereich, Observablen, Invarianten, Hierarchie offener Fragen |

EABC erfüllt die zweite Beschreibung: Es legt $G_E$, $D_E(X)$, Skalierungsobservablen ($R_\beta$, $\alpha_{\mathrm{loc}}$, $\alpha_E$) und eine **Hierarchie offener Fragen** fest — von der Größenordnung von $|D_E|$ bis zur möglichen Endfrage $\Phi_E$. Die Lakatos-Einordnung (§2, vgl. `collatz_eabc_zirkulationshypothese.md` §4) ergänzt diese Lesart; sie widerspricht ihr nicht.

**Vierter Perspektivwechsel (§4.7):** objektzentriert $\longrightarrow$ **zustandszentriert** — nicht „Welche Eigenschaft haben Primzahlen?", sondern „Welcher Zustand genügt zur Beschreibung des Prozesses?"

**Kanonische Hierarchie (state-centric):**
$$\text{Primstrom} \longrightarrow \text{Streaming-Monoid} \longrightarrow \text{Observablen} \longrightarrow \text{Interpretation}$$
$$P \;\longrightarrow\; Z \;\longrightarrow\; (D_E,\,Q_E,\,W_E,\,\ldots) \;\longrightarrow\; \Phi_E$$

$$\boxed{\;Z\ \text{ist das fundamentale Objekt; } D_E,\,Q_E,\,W_E\ \text{sind Funktoren } f_D,\,f_Q,\,f_W\ \text{auf } Z;\ \Phi_E\ \text{ist nur Interpretationsschicht (Ebene 4).}\;}$$

**Dritter Perspektivwechsel (§4.6):** $\Phi_E \longrightarrow D_E \longrightarrow Z$ — vom asymptotischen Grenzwert zur **algebraischen Streaming-Struktur** $(Z,\oplus,Z_0)$ (Ebene 1: Monoid-Theorem in `OccupancyMonoid.lean`; Ebene 2: EABC-Interpretation). $N_\pm$, $D_E$, $Q_E$ werden Auswertungen eines Endzustands, nicht Primärdefinition.

$$\boxed{\;\text{Der größte Fortschritt ist nicht die Einführung neuer Größen, sondern die Entkopplung des Programms von }\Phi_E.\;}$$

$$\boxed{\;\text{Lakatos-Hard-Core (voll beweisbar): } (Z,\oplus,Z_0)\ \text{kommutatives Monoid — unabhängig von Primzahlen, HL, RH, Holonomie.}\;}$$

**Harter Kern (empirisch/numerisch):** (1) gerichteter EABC-Kreisgraph $G_E$, (2) orientierte Zyklusdifferenz $D_E(X)=N_+(X)-N_-(X)$, (3) daraus induzierte Skalierungsobservablen — **nachgelagert** an $Z$ via $f_D,\,f_Q,\,f_W$. Erst darauf bauen Holonomie- und $\Phi_E$-Fragen (Interpretationsschicht).

**Erster numerischer Belastungstest** (vgl. `collatz_eabc_zirkulationshypothese.md` §4.3): Vollauf bis $X=10^{10}$ (13 Checkpoints) — prüft nur $(G_E,N_\pm,D_E)$, $R_{1/2}=D_E/\sqrt{Q}$ und Skalierungsobservablen; **kein** $\Phi_E$, **keine** Holonomie. **Ergebnis (Experiment):** bis $10^{10}$ kein numerischer Hinweis auf Holonomie ($\Phi_E\neq 0$) oder sublinearen Bias ($\alpha_E>\tfrac{1}{2}$); $D_E$ kompatibel mit Wurzelfehlerterm ($|D_E|=O(\sqrt{Q})$). H₀a/H₀b numerisch bevorzugt gegenüber H₂/H₃. Holonomie bleibt Endfrage (Ebene 4).

**Literaturpositionierung und Novelty** (vgl. `collatz_eabc_zirkulationshypothese.md` §4.5): EABC ist **keine** neue Zahlentheorie ex nihilo, sondern eine **gerichtete Observable** $D_E(X)=N_+-N_-$ auf dem bekannten Feld Prime Races / consecutive-prime biases / HL-$k$-Tupel / Sieb-Lückenzyklen. Neu ist nicht mod $12$, sondern die chirale ABCEA/CEABC-Zählung auf $G_E$. Die Collatz$\leftrightarrow$EABC$\leftrightarrow$Quaternion-Brücke ist **programminterne Heuristik** (Schicht **C**), nicht etablierte Literatur.

**Stabilisierungskette (klassisch, Observable-first):** ($\leadsto$ = führt zu Fragen; $\Rightarrow$ = bewiesene Implikation)
$$G_E \;\leadsto\; (D_E,\,Q_E) \;\leadsto\; (R_\beta,\,\alpha_{\mathrm{eff}},\,\alpha_{\mathrm{loc}}) \;\leadsto\; (\alpha_E,\,W_E) \;\leadsto\; \Phi_E.$$

**Stabilisierungskette (state-centric, kanonisch ab §4.7):**
$$P \;\Longrightarrow\; (Z,\oplus,Z_0) \;\Longrightarrow\; (D_E,\,Q_E,\,W_E) \;\leadsto\; \Phi_E.$$
(Monoid-Schritt = **Theorem**; Observablen = Funktoren auf $Z$; $\Phi_E$ = Interpretation.)

**Referee-Perspektive:** „Ich glaube noch nicht $\Phi_E\neq 0$ — aber $F(P_1\sqcup P_2)=F(P_1)\oplus F(P_2)$ ist bewiesen und $D_E(X)$ ist definiert." Verschiebung: von „Ist die Vermutung wahr?" zu „Welcher Zustand beschreibt den Prozess?" und „Welche Eigenschaften besitzt die Observable?"

---

## 0.1 Leitprinzip: zwei parallele Spuren

$$\boxed{\;\textbf{Beweis-Spur} \;\|\; \textbf{Ikone-Spur}\;}$$

| Spur | Ziel | Erlaubt | Verboten |
|------|------|---------|----------|
| **Beweis** | Lean-`Prop`, kombinatorische Theoreme, spektrale Identitäten | Definition, Lemma, `sorry` nur mit Label „offen“ | Physik als Theorem |
| **Ikone** | Didaktik, Modellabbildung, historische Parallelen | Analogie, Metapher, Forschungsfrage | `theorem` mit Physikbehauptung |

**Regel:** Ikone-Spur-Inhalte dürfen **nie** in die Beweis-Spur migrieren, ohne explizite mathematische Formalisierung und Schicht-Upgrade (typisch B $\leftarrow$ C).

---

## 1. Schicht A — hart (kombinatorisch / algebraisch)

**Label:** **Theorem** | **Definition** | **Conjecture** — Lean-provable, **unabhängig von Physik**

**Inhalt (vollständig):**

| Thema | Objekt / Behauptung | Artefakt |
|-------|---------------------|----------|
| $24I_3$ | Ikosaeder-/Dreieckssymmetrie, Hurwitz-Einheiten | `collatz_hurwitz_polytop_eabc.tex`, `collatz_eabc_hurwitz_orbit_test.py` |
| Prime-Defekte | Defekt $D_A(x)=x-\Pi_A(x)$, minimale irreduzible Zustände | `collatz_eabc_euklidische_hebung.md`, `collatz_eabc_normabstieg_hypothese.md` |
| $M^+ = 24I_3 + w_p v v^\top$ | Metrik-/Projektions-Erweiterung (algebraisch) | `collatz_morley_metrik_erweiterung.md` |
| Basel | $\zeta(2)$, Bernoulli-Normschalen, Fensterkonstanten | `collatz_eabc_bernoulli_uebersetzung.md` |
| Fenstermassen | Gleitfenster $C_n^{(5)}$, ABCEA/CEABC-Zählung | `collatz_eabc_holonomie_fehlerterm.py` |
| Covering | Überlagerung $C_4$, Windungszahl | `collatz_eabc_diskrete_geometrie.md` §0–2 |
| Fraktalskalen | Skalenhierarchie $\Sigma_n$, Normniveaus | `collatz_eabc_plattenuebergang.md`, `collatz_eabc_quaternion_mass_hypothese.md` |
| Holonomie / Zirkulation | $D_E=N_+-N_-$, Lücken $(2,4,2,4)$, Taubenloch | `CollatzEabc/HolonomieFehlerterm.lean`, `collatz_eabc_fehlerterm_hypothese.md` |

**Lean-Anker:** `CollatzEabc.Core`, `HolonomieFehlerterm`, `Kappa`, `ArithLanguage`, `FlussPhiE` (**Schicht A** — kombinatorischer Teil).

**Epistemische Grenze:** Alles in Schicht A ist **arithmetisch oder kombinatorisch** — keine Thomson-Schale, kein Strahlungsraum, keine Quaternionenrotation als „physikalisch“.

---

## 2. Schicht B — geometrisch (diskrete Differential- / Spektralgeometrie)

**Label:** **Definition** | **Theorem** (endlich) | **Vermutung** (asymptotisch) — **vollständig Mathematik**

**Inhalt (vollständig):**

| Thema | Objekt | Artefakt |
|-------|--------|----------|
| Holonomien | $\mathrm{Hol}_E$, $\mathcal{P}\exp\oint A$ (Programm) | `collatz_eabc_holonomie_stufen.md` |
| Kreisgraphen | $G_E=(V,E)$, $C_4\cong S^1$, $E^+/E^-$ | `collatz_eabc_diskrete_geometrie.md`, `FlussPhiE.lean` |
| Harmonische Formen | $h\in H^1(C_4,\mathbb{Z})$, $\langle\omega_E,h\rangle$ | `collatz_eabc_uebergangsraum.md`, `collatz_eabc_hodge_eabc.py` |
| Fluss | $\omega_E$, $C_E=\oint_\gamma\omega_E$, $W_E(X)$, $\alpha_E$, Ebenen 0–4, H₀a/H₀b–H₃, Literaturpositionierung §4.5, Level-2 §4.8, Spin-Liquid §4.8.1 | `collatz_eabc_zirkulationshypothese.md` §4.1–4.5, §4.8–4.8.1 |
| Orientierungsklassen | $\Phi_E=\lim W_E$, chirale 1-Form | `collatz_eabc_chirale_polarisation.md` |

**Kernfrage (Schicht B als Struktur, Schicht R als Brücke):**

$$\boxed{\;\text{Das EABC-Programm untersucht primär die Wachstumsordnung der orientierten Zyklusdifferenz }D_E(X)\text{, während die Holonomie }\Phi_E\text{ als mögliche Endstufe dieser Skalierungstheorie erscheint.}\;}$$

$$\boxed{\;\text{Nicht }\Phi_E\text{ ist der Anfang, sondern }D_E(X).\;}$$

$$\boxed{\;\text{Ob die Wachstumsordnung von }D_E\text{ bis zur linearen Skala reicht und damit eine nichtverschwindende Holonomie }\Phi_E\text{ erzeugt — sekundäre Frage (Ebene 4).}\;}$$

**Paradigmenwechsel:** Früher Grenzwerttheorie $W_E\to\Phi_E\neq 0$ (Ergodentheorie/Dichtetheorie); jetzt Fehlertermtheorie $D_E(X)=A(X)-C(X)$ (analytische Zahlentheorie). Die Verschiebung ändert das **mathematische Objekt**, nicht nur die Methode.

$$\boxed{\;\text{Das Programm hängt nicht mehr an einer Vermutung.}\;}$$

**Lakatos-Einordnung** (vgl. `collatz_eabc_zirkulationshypothese.md` §4; ergänzt §0 Forschungsprogramm vs. Theorie):

| Rolle | Objekt | Ebene |
|-------|--------|-------|
| **Harter Kern** | $G_E$; $D_E(X)=N_+-N_-$; induzierte Skalierungsobservablen | 0–1 |
| **Primäre Theorie** | $D_E$, $Q$ | 1 |
| **Sekundäre Schutzmantel** | $R_\beta$, $\alpha_{\mathrm{loc}}$ (Numerik primär), $\alpha_{\mathrm{eff}}$, $\alpha_E$ | 2 |
| **Orientierung** | $W_E$ | 3 |
| **Endfrage** | $\Phi_E$ | 4 |

**V1 → V2:** V1 startete mit $\Phi_E\neq 0$ und nutzte Hilfsgrößen zur Stützung; V2 macht $\Phi_E$ zur **Endfrage** der $D_E$-Skalierungstheorie. Stabilisierungskette: $G_E \leadsto (D_E,Q_E) \leadsto (R_\beta,\alpha_{\mathrm{eff}},\alpha_{\mathrm{loc}}) \leadsto (\alpha_E,W_E) \leadsto \Phi_E$. Selbst $\Phi_E=0$ bleibt interessant ($\alpha_E$, $R_\beta$, $\alpha_{\mathrm{loc}}$).

**Vorwärtskette** ($\leadsto$ offene Fragen; $\Rightarrow$ bewiesene Implikationen): $G_E \leadsto D_E(X) \leadsto \alpha_E \leadsto W_E(X) \leadsto \Phi_E$.

**Fünf Ebenen** (vgl. `collatz_eabc_zirkulationshypothese.md` §4.2):

| Ebene | Inhalt |
|-------|--------|
| **0** | Geometrie: $G_E$, $\gamma^+$ (ABCEA), $\gamma^-$ (CEABC) |
| **1** | Zirkulationsfehler: $D_E(X)=N_+-N_-$ (primäre Observable) |
| **2** | Skalierung: $R_\beta$, $\alpha_{\mathrm{loc}}$ (Numerik primär), $\alpha_E$ |
| **3** | Orientierung: $W_E(X)=D_E/Q$ |
| **4** | Holonomie (Endfrage): $\Phi_E=\lim W_E$ |

$$W_E(X) = \frac{N_+(X)-N_-(X)}{N_+(X)+N_-(X)},\qquad
R_\beta(X) = \frac{D_E(X)}{Q(X)^\beta},\qquad
\alpha_E := \inf\{\beta : R_\beta\text{ bleibt beschränkt}\}$$

$$\boxed{\;\Phi_E \neq 0 \;\Rightarrow\; D_E(X)\sim\Phi_E Q(X) \;\Rightarrow\; \alpha_E=1.\;}$$

**Nicht äquivalent:** $\alpha_E=1$ impliziert keinen Grenzwert von $W_E$. Gegenbeispiel: $D_E(X)=Q(X)\sin(\log\log Q(X))$. $\Phi_E$ (Orientierung) ist **stärker** als $\alpha_E$ (Skala). Scheitern von $\Phi_E\neq 0$ zerstört das Programm nicht — eigenständige Theorie von $D_E$, $R_\beta$, $\alpha_{\mathrm{loc}}$, $\alpha_E$. Selbst bei $\Phi_E=0$: Fragen zu $\alpha_E=\tfrac{1}{2}$?, $\alpha_E>\tfrac{1}{2}$?, $\alpha_{\mathrm{loc}}$-Plateaus?, kritische $R_\beta$ bleiben offen.

| Status | Aussage |
|--------|---------|
| **GREEN** | $-1 \le W_E(X) \le 1$ (bewiesen: `W_E_bounds`) |
| **RED** | $\lim_{X\to\infty} W_E(X) = \Phi_E \neq 0$ — **H₃** (Holonomie; Ebene 4; Vermutung: `HasNonzeroHolonomyLimit`, `EABC_holonomy_limit_conjecture`) |

Lean: `phi_E_conjecture` (**B**, `Prop`) vs. `phi_E_conjecture_statement` (**R**, `sorry`). Methodische Lesart — $W_E(X)$ stabilisiert sich asymptotisch gegen einen Grenzwert; keine „statistische Ruhe im Unendlichen“, keine „unumkehrbare Tendenz“: `collatz_eabc_zirkulationshypothese.md` §4.1. **Fünf Ebenen in §4.2** (0–4): (0) harter Kern $G_E$; (1) primär $D_E$, $Q$; (2) sekundär $R_\beta$, $\alpha_{\mathrm{loc}}$, $\alpha_{\mathrm{eff}}$, $\alpha_E$; (3) Orientierung $W_E$; (4) Endhypothese $\Phi_E$. **Erster Belastungstest** (§4.3): nur harter Kern, nicht Holonomie. Prime-Race-Analogie $\pi(x)$–$\mathrm{Li}(x)$: Größenordnung von $|D_E|$ vor Grenzwert von $W_E$. H₀b ($W_E\to 0$) ist analytisch und **nicht** mit H₀a identisch; $\alpha_E>\tfrac{1}{2} \nRightarrow \Phi_E\neq 0$ (vgl. §4.2). Numerik (#73): `eabc_quadruplets_1e10.py`, `eabc_quadruplets_fit_alpha.py`, `eabc_quadruplets_plot.py`.

**Abgrenzung zu Schicht C:** Sagnac, AB-Phase, magnetischer Laplace sind **Analogien** auf $G_E$ — geometrisch sauber in B, physikalisch ikonisch erst in C. Priminduzierte Zählung $N_\pm(X)$ ist **B**, nicht „laufende“ Primzahlen (**C**-Ikone).

---

## 3. Schicht C — ikonische Physik (Interpretation, keine mathematische Konsequenz)

**Label:** **Analogie** | **Modellabbildung** | **Ikone** — **nicht** als Theorem exportierbar

**Inhalt (vollständig):**

| Ikone | EABC-Lesart | Warnung |
|-------|-------------|---------|
| Thomson | Kugel-/Schalenbild für Normniveaus $\Sigma_n$ | Emergenz-Hypothese, kein Beweis |
| Elektronenschalen | $n\mapsto\mu_n$, diskrete Spektralschichten | `collatz_eabc_quaternion_mass_hypothese.md` §12 |
| Strahlungsraum | Platzhalter für künftige $\varepsilon$-Schnittstelle | → **rote Schicht** |
| Shell-Sprünge | Übergänge $\Sigma_n\to\Sigma_{n+1}$, $Z(n)$, $\Delta Z(n)$ | `collatz_eabc_plattenuebergang.md` |
| Quaternionenrotationen | $SU(2)$-Phase auf Transportbündel | `collatz_eabc_chirale_transport.py` |
| Near-Zero als „physikalisch“ | near-zero-Moden in $\mathrm{Spec}(L_{\mathrm{mag}})$ | **Experiment**, nicht Theorem |
| Spin Liquid / Korrelationsgeometrie | $\langle X\rangle\approx 0$, Struktur in $\langle XX^T\rangle$; Herbertsmithite-Muster | **Analogie** — methodisch, **nicht ontologisch** (§3.1, `collatz_eabc_zirkulationshypothese.md` §4.8.1) |

**Explizit:** Interpretationen in Schicht C sind **keine mathematischen Konsequenzen** aus A oder B. Ein numerischer Befund (z. B. near-zero-Mode) bleibt **Experiment**, bis er in B formalisiert ist.

**Verknüpfung:** `collatz_eabc_epistemik_physik.md` — Wegfunktion ohne Zeitdilatation, kein Zwillingsparadoxon, Sagnac didaktisch.

### 3.1 Vier-Ebenen-Observablen-Stack und Spin-Liquid-Analogie (methodisch, nicht ontologisch)

$$\boxed{\;\text{Spin-Liquid-Analogie = Schicht \textbf{C}, \textbf{methodisch} — \textbf{nicht ontologisch}.}\;}$$

Ergänzend zu den **fünf Ebenen** 0–4 in Schicht **B** (Geometrie → Fehlerterm → Skalierung → Orientierung → Holonomie; vgl. `collatz_eabc_zirkulationshypothese.md` §4.2) führt das Level-2-Programm einen **vierstufigen Observablen-Stack** ein:

| Level | Inhalt | Schicht / Label |
|-------|--------|-----------------|
| **0** | $(Z,\oplus)$ Streaming-Algebra | **A** — **Theorem** (`OccupancyMonoid.lean`, §4.6–4.7) |
| **1** | $D_E$, $W_E$, $\Phi_E$ — Drift / Holonomie | **B**/**R** — **Experiment** (§4.3: Level-1-Null) |
| **2** | $\Sigma_A$, $\mathrm{Spec}(\Sigma_A)$, $\Delta_F$ — Fluktuationsgeometrie | **B** — **Definition** + **Experiment** (§4.8) |
| **3** | Witness, Emergenz, Korrelationsmodelle | **C** — **Analogie** / **Hypothese** (§4.4, §4.8.1) |

**Spin-Liquid-Analogie** (vollständig: `collatz_eabc_zirkulationshypothese.md` §4.8.1): In der modernen Vielteilchenphysik kann $\langle X\rangle\approx 0$ bei stark strukturiertem $\langle XX^T\rangle$ vorliegen (Spin Liquids, QHE, topologische Isolatoren, …). EABC-Parallel: Level-1 ($W_E\to 0$, $D_E=O(\sqrt{Q})$) verschwindet, Level-2 ($\Sigma_A^{\mathrm{prime}}\neq\Sigma_A^{\mathrm{rand}}$) bleibt von Nullmodellen unterscheidbar. Der mod-$60060$-**Witness** §4.4 illustriert dasselbe epistemische Muster auf der HL-/Kanalseite.

$$\boxed{\;\text{Kein Anspruch: }\Sigma_A\text{ sind Spinons oder Primzahlen bilden ein Spin Liquid. Nur methodische Parallele: Struktur in Fluktuationsgeometrie, nicht im Mittelwert.}\;}$$

**Nullmodell-Leiter** (vollständig: `collatz_eabc_zirkulationshypothese.md` §4.8.2): Level-2 testet $\Sigma_A^{\mathrm{prime}}$ gegen eine **Hierarchie** realistischerer Ensembles — nicht gegen $\mu_A\approx 0$:

| Stufe | Null | Frage |
|-------|------|-------|
| 1 | $\Sigma_A^{\mathrm{perm}}$ | Abweichung von marginaltreuer Zufallsreihenfolge? |
| 2 | $\Sigma_A^{\mathrm{Markov}}$ | Abweichung bei erhaltenen lokalen Übergängen? |
| 3 | $\Sigma_A^{\mathrm{HL}}$ | Abweichung bei HL-konsistenten Korrelationen? |
| 4 | Theorie | Mechanismus / neue Mathematik |

$$\boxed{\;\text{Gegner ist falsches Nullmodell — nicht der verschwindende Mittelwert.}\;}$$

**Label:** gesamter Abschnitt §3.1 = **Analogie** (Schicht **C**); verknüpft mit **Experiment** (Level 1–2) und **Beispiel** (Witness).

---

## 4. Rote Schicht — interpretativer Scaffold (keine Theoreme)

**Farbe:** explizit **nicht** Teil des Beweisprogramms.

**Zweck:** Zukünftige Schnittstelle zwischen **metrischen Defekten** (Schicht A: $D_A$, Prime-Defekte) und **strahlungsartigen Strukturen** (ikonisch: Rydberg $E_n\sim 1/n^2$, $\varepsilon$-Feld).

$$\boxed{\;\text{RadiationSpace} = \text{markierter, unbewiesener Verbindungsraum — kein Theorem.}\;}$$

**Lean-Stil (kanonisch):**

```lean
/-
Interpretative scaffold. No physical meaning asserted.
Serves as future interface between metric defects and radiative structures.
-/
class RadiationSpace where
  carrier : Type
  epsilon : State → carrier
```

**Verboten in Lean:** `theorem radiation_space_exists`, `theorem rydberg_from_eabc`, jede `sorry`-freie Physikbehauptung.

**Erlaubt:** Klasse/Struktur, Kommentare, `State`-Stub, Dokumentation.

**Artefakt:** `collatz_eabc_core/CollatzEabc/RadiationSpace.lean`

**Zukünftige Schnittstelle (nur Kommentar, nicht formalisiert):**
- $\varepsilon : \mathrm{State}\to\mathrm{carrier}$ — effektives „Strahlungsfeld“ auf EABC-Zuständen
- Rydberg-Skalierung $E_n\sim 1/n^2$ als **Lesart** für Schalenindex $n$, nicht als abgeleitetes Gesetz

---

## 5. Dokument-Index mit Schicht-Tags

| Dokument | Schicht | Kurzrolle |
|----------|---------|-----------|
| `collatz_eabc_diskrete_geometrie.md` | **B** | Kanonische Geometrie $G_E$, $\Phi_E$ |
| `collatz_eabc_epistemik_physik.md` | **C** (+ Abgrenzung) | Physik vs. Modell, Nicht-SRT |
| `collatz_eabc_epistemik_schichten.md` | **Meta** | Dieses Framework |
| `collatz_generalangriff_2026.md` | **A**–**C** | Forschungsorganisation, Brücke L |
| `collatz_mathlib_eabc_kandidaten.md` | **A**–**B** | Mathlib + Lean-Module |
| `collatz_eabc_euklidische_hebung.md` | **A** | Hurwitz-Kette, Defekt/Hebung |
| `collatz_eabc_holonomie_stufen.md` | **B** | Stufen 1–3, Wachstum A/B/C in $N$ |
| `collatz_eabc_zirkulationshypothese.md` | **B** (+ **C** §4.8.1) | $N_\pm$, $C_E$, $D_E$; Level-2 $\Sigma_A$, $\Delta_F$; Spin-Liquid-Analogie (methodisch) |
| `collatz_eabc_fehlerterm_hypothese.md` | **A**–**B** | Prime Race, $D_E$ |
| `collatz_eabc_kritische_abbildung.md` | **B** (+ C-Analogie) | Holonomie-Sensor, Wegparadoxon |
| `collatz_eabc_chirale_polarisation.md` | **B**–**C** | Helizität als Modell |
| `collatz_eabc_quaternion_mass_hypothese.md` | **A**–**C** | Spektralgeometrie / Schalen-Ikone |
| `collatz_eabc_plattenuebergang.md` | **A**–**C** | Normschichtfolge |
| `collatz_eabc_brachistochrone.md` | **C** | Fermat-Modell |
| `collatz_eabc_sagnac.md` | **C** | Didaktische Intuition |
| `CollatzEabc/RadiationSpace.lean` | **rot** | Scaffold only |
| `CollatzEabc/FlussPhiE.lean` | **A**–**B**–**R** | $C_4$, $h$ bewiesen (**A**); $W_E$, $\Phi_E$ (**B**); EABC-Vermutung (**R**) |
| `CollatzEabc/HolonomyCore.lean` | **B**–**R** | **GREEN** `W_E_bounds`; **RED** `HasNonzeroHolonomyLimit` |
| `CollatzEabc/HolonomieFehlerterm.lean` | **A**–**B**–**R** | Lücken bewiesen (**A**); Zähldefinitionen (**B**); Prime `sorry` (**R**) |

---

## 6. Historische Parallelen (nur Dokumentation)

**Label:** **Ikone** / **Heuristik** — **nicht** formalisiert, **nicht** in Lean.

| Figur | Parallele zum EABC-Programm | Schicht |
|-------|----------------------------|---------|
| **Riemann** | Flächen / Überlagerungen / analytische Fortsetzung ↔ $C_4$-Überlagerung, $\Phi_E$ als globale Phase | B |
| **Dirac** | Diskrete Symmetrie + globale Phase (monopole/Spin-Lesart) ↔ chirale Transportphases $\phi_R-\phi_L$ | C |
| **Grothendieck** | Kohomologie als Schnittstelle lokaler Daten ↔ $H^1(C_4)$, $\langle\omega_E,h\rangle$, Defektfunktoren | A–B |

Diese Parallelen dienen der **Orientierung** in der Ikone-Spur. Sie begründen **keine** Theoreme und ersetzen **keine** Beweisstrategie in A/B.

---

## 7. Migrationsregeln (Schicht-Upgrade)

| Von | Nach | Bedingung |
|-----|------|-----------|
| C → B | Ikone als diskrete Definition auf $G_E$ | Explizite mathematische Objekte + Lean-`Prop` |
| B → A | Asymptotik → endliche kombinatorische Identität | Beweis ohne Spektral-Limes |
| rot → B | `RadiationSpace` mit beweisbarer Abbildung | Kein $\varepsilon$ ohne mathematische Domäne |
| beliebig → Theorem | — | **Verboten** ohne Peer-Review der Schicht |

---

## 8. Boxed Zusammenfassung

$$\boxed{\;\text{A = hart (24I}_3\text{, Defekte, Fenster, Holonomie kombinatorisch).}\;}$$

$$\boxed{\;\text{B = Struktur (Kreisgraph, Fluss, Harmonik, Grenzwert-`Prop`).}\;}$$

$$\boxed{\;\text{R = Forschungsbrücke (Prim-Enumeration, asymptotik, } \Phi_E \neq 0 \text{ — `sorry` erlaubt).}\;}$$

$$\boxed{\;\text{C = ikonische Physik (Thomson, Schalen, Quaternionen) — Interpretation, kein Export.}\;}$$

$$\boxed{\;\text{Rot = RadiationSpace-Scaffold — Schnittstelle für } \varepsilon \text{ und Rydberg-Lesart, kein Theorem.}\;}$$

$$\boxed{\;\textbf{Beweis} \;\|\; \textbf{Ikone}\;\text{ — parallele Spuren, strikt getrennt.}\;}$$

---

*Schichten-Framework PR #59 + euklidische-hebung. Physikmetaphern bleiben in C/rot; Lean beweist nur A/B (und kombinatorische Teile von B).*
