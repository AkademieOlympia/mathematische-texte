# EABC Φ_E — analytischer Lean-Beweisversuch

**Status:** Skeleton (PR #63; Branch `collatz/eabc-05-holonomie-fehlerterm`)  
**Lean (Fehlerterm):** `collatz_eabc_core/CollatzEabc/HolonomieFehlerterm.lean`  
**Lean (Pattern-Zählung):** `collatz_eabc_core/CollatzEabc/PatternCount.lean`  
**Lean (Hodge-Layer):** `collatz_eabc_core/CollatzEabc/FlussPhiE.lean`  
**Lean (minimaler Kern):** `collatz_eabc_core/CollatzEabc/HolonomyCore.lean` (`EABC`-Namespace)  
**Epistemik:** `collatz_eabc_epistemik_schichten.md` — **A / B / R / C**; Lakatos §4 (#72)  
**Kanonische Geometrie:** `collatz_eabc_diskrete_geometrie.md`  
**Numerik:** `eabc_quadruplets_1e10.py`, `eabc_quadruplets_fit_alpha.py` (#73); `collatz_eabc_hodge_eabc.py` (`Phi_E`, `flux_density_limit`, `inner_product_omega_h`)

---

## Epistemische Schichten in `FlussPhiE.lean`

| Schicht | Label | Lean-Rolle | Sorry? |
|---------|-------|------------|--------|
| **A** | Theorem | Kombinatorisch bewiesen | nein |
| **B** | Struktur | Definitionen, `Prop`-Skelette | nein |
| **R** | Forschungsbrücke | Asymptotik, Prim-Enumeration, EABC-Vermutung | ja |
| **C** | Ikone | Nur Markdown (Physik-Analogien) | — |

**Leitsatz:** Theorem ≠ Struktur ≠ Brücke ≠ Ikone

---

## Programm

**Status:** EABC ist wissenschaftstheoretisch eher **Forschungsprogramm** als bereits **Theorie** — es definiert Gegenstandsbereich, Observablen, Invarianten und eine Hierarchie offener Fragen (vgl. `collatz_eabc_epistemik_schichten.md` §0, `collatz_eabc_zirkulationshypothese.md` §4).

EABC = priminduzierter Fluss auf $C_4 \cong S^1$ mit

- $G_E = (V,E)$, $V=\{E,A,B,C\}$ — **harter Kern** (Ebene 0)
- $E^+ = \{EA, AB, BC, CE\}$, $E^- = \{EC, CB, BA, AE\}$
- orientierte Zyklen $\gamma^+ = \mathrm{ABCEA}$, $\gamma^- = \mathrm{CEABC}$
- $D_E(X) = N_+(X) - N_-(X)$ (Ebene 1: primäre Observable; **harter Kern**)
- Skalierungsobservablen $R_\beta$, $\alpha_{\mathrm{loc}}$, $\alpha_E$ (Ebene 2; aus $D_E$ induziert)
- $W_E(X) = D_E(X)/Q(X)$ (Ebene 3)
- $\Phi_E = \lim_{X\to\infty} W_E(X)$ (Ebene 4: **Endfrage**, nicht Ausgangsannahme)

$$\boxed{\;\text{Der größte Fortschritt ist nicht die Einführung neuer Größen, sondern die Entkopplung des Programms von }\Phi_E.\;}$$

$$\boxed{\;\text{Das EABC-Programm untersucht primär die Wachstumsordnung der orientierten Zyklusdifferenz }D_E(X)\text{, während die Holonomie }\Phi_E\text{ als mögliche Endstufe dieser Skalierungstheorie erscheint.}\;}$$

$$\boxed{\;\text{Nicht }\Phi_E\text{ ist der Anfang, sondern }D_E(X).\;}$$

**Lakatos-Einordnung** (`collatz_eabc_zirkulationshypothese.md` §4): harter Kern $G_E$, $D_E$, induzierte Skalierungsobservablen (Ebene 0–1); primär $D_E$, $Q$ (Ebene 1); sekundär $R_\beta$, $\alpha_{\mathrm{loc}}$ (Numerik primär), $\alpha_{\mathrm{eff}}$, $\alpha_E$ (Ebene 2); Orientierung $W_E$ (Ebene 3); **Endfrage** $\Phi_E$ (Ebene 4).

**V1 → V2:** V1 startete mit $\Phi_E\neq 0$ und Hilfsgrößen zur Stützung; V2 macht $\Phi_E$ zur **Endfrage** der $D_E$-Skalierungstheorie. Stabilisierungskette: $G_E \leadsto (D_E,Q_E) \leadsto (R_\beta,\alpha_{\mathrm{eff}},\alpha_{\mathrm{loc}}) \leadsto (\alpha_E,W_E) \leadsto \Phi_E$.

**Paradigmenwechsel:** Früher Grenzwerttheorie $W_E\to\Phi_E\neq 0$; jetzt Fehlertermtheorie $D_E(X)=A(X)-C(X)$ (analytische Zahlentheorie). Referee-Perspektive: $D_E(X)$ ist definiert und ihre Größenordnung untersuchbar — unabhängig davon, ob $\Phi_E\neq 0$ je bestätigt wird.

**Vorwärtskette** ($\leadsto$ offene Fragen; $\Rightarrow$ bewiesene Implikationen): $G_E \leadsto D_E(X) \leadsto \alpha_E \leadsto W_E(X) \leadsto \Phi_E$.

**Endfrage (Schicht R):** $\Phi_E \neq 0$ bzw. $\langle\omega_E, h\rangle \neq 0$ für kanonisches harmonisches $h$. Dies ist **H₃** (starke Holonomie: $W_E(X)\to\Phi_E\neq 0$) — **Ebene 4** des Programms, nach Geometrie (Ebene 0), Zirkulationsfehler $D_E$ (Ebene 1), Skalierung $\alpha_{\mathrm{loc}}/\alpha_E$ (Ebene 2) und Orientierung $W_E$ (Ebene 3). Entspricht Lean-**RED** `HasNonzeroHolonomyLimit`. Vorgelagert (`collatz_eabc_zirkulationshypothese.md` §4.2): zuerst $D_E$ und $R_\beta$, dann $\alpha_{\mathrm{loc}}$, dann $\alpha_E$, dann Orientierung (H₀b vs. H₃). $\alpha_E>\tfrac{1}{2} \nRightarrow \Phi_E\neq 0$. Scheitern von H₃ zerstört das Programm nicht; selbst bei $\Phi_E=0$ bleiben Fragen zu $\alpha_E$, $\alpha_{\mathrm{loc}}$-Plateaus und kritischen $R_\beta$ offen.

$$\boxed{\;\Phi_E \neq 0 \;\Rightarrow\; D_E(X)\sim\Phi_E Q(X) \;\Rightarrow\; \alpha_E = 1.\;}$$

**Nicht äquivalent:** $\alpha_E=1$ allein garantiert weder $\Phi_E\neq 0$ noch $\lim W_E$. Gegenbeispiel: $D_E(X)=Q(X)\sin(\log\log Q(X))$ — $\alpha_E=1$, aber $W_E$ ohne Grenzwert. $\Phi_E$ (Orientierungsparameter) ist **stärker** als $\alpha_E$ (Skalenparameter).

---

## Lean-Modul `HolonomyCore.lean` (sauberer Split)

**GREEN/RED im Modul:** GREEN = `Node` … `W_E_bounds` (bewiesen); RED = `HasNonzeroHolonomyLimit`, `EABC_holonomy_limit_conjecture` (`sorry`) — formal **H₃** (Ebene 4: Holonomie-Grenzfall, $\alpha_E=1$).

Minimale Architektur im Namespace `EABC` — ohne Prim-Enumeration, ohne Hodge-Layer:

| Objekt | Lean-Name | Status |
|--------|-----------|--------|
| Kreisgraph-Knoten | `Node`, `next`, `prev` | **Definition** |
| Zykluszählung | `CycleCounts`, `circulation`, `size` | **Definition** |
| $W_E$ auf endlicher Stichprobe | `phiApprox`, `W_E` | **Definition** |
| $-1 \le W_E \le 1$ | `phiApprox_bounds`, `W_E_bounds` | **Theorem** (bewiesen) |
| Fluss bis Schranke $X$ | `EABCFlow`, `C_E`, `S_E` | **Definition** |
| $\Phi_E \neq 0$ (asymptotisch) | `HasNonzeroHolonomyLimit`, `EABC_holonomy_limit_conjecture` | **Vermutung** (`sorry`; **H₃**, Ebene 4, $\alpha_E=1$) |

**Grenze Beweis vs. Vermutung:** `W_E(X)=\frac{N_+(X)-N_-(X)}{N_+(X)+N_-(X)}$ liegt formal in $[-1,1]$;
die asymptotische Aussage $\lim_{X\to\infty} W_E(X)=\Phi_E\neq 0$ steht in `EABC_holonomy_limit_conjecture` (`Tendsto` in $\mathbb{R}$), nicht als exakte rationale Endkonstante.
Die zu starke Variante (eventuell konstantes rationales $\Phi$) ist nur als auskommentiertes Scaffold erhalten.

`FlussPhiE.lean` ergänzt denselben Fluss um C₄-Kanten, harmonisches $h$, Prim-Brücken (`HolonomieFehlerterm`).  
**Zentrale Observablen** $D_E$, $Q_E$, $W_E$, $R_\beta$, $\tilde D_E$ liegen in `HolonomieFehlerterm`; `FlussPhiE` aliasiert $C_E$, $S_E$ und die Grenzwert-`Prop`s.

```bash
cd collatz_eabc_core
lake build CollatzEabc.HolonomyCore
lake build CollatzEabc.PatternCount
lake build CollatzEabc.HolonomieFehlerterm
lake build CollatzEabc.FlussPhiE
```

---

## Lean-Modul `HolonomieFehlerterm.lean` (Fehlerterm-Schicht)

**Fokus:** Ebene 0–2 ($D_E$, $Q_E$, $R_\beta$, $\tilde D_E$) — nicht Holonomie-Beweis (Ebene 4).

### Schicht A — Theorem (bewiesen, endliche Folge / Fenster)

| Aussage | Lean-Name | Status |
|---------|-----------|--------|
| $D_E = N_+ - N_-$ | `D_E`, `D_E_eq_diff` | **Theorem** |
| $Q_E = N_+ + N_-$ | `Q_E`, `Q_E_eq_sum` | **Theorem** |
| $-1 \le \chi_{\mathrm{Hol}} \le 1$ | `chi_Hol_bounds` | **Theorem** |
| $N_+=N_- \Rightarrow \chi_{\mathrm{Hol}}=0$ | `chi_Hol_zero_of_balance` | **Theorem** |
| Lückenmuster $(2,4,2,4)$ | `abcea_gap_pattern`, `ceabc_gap_pattern` | **Theorem** |
| Taubenloch auf ABCE-Fenstern | `pigeonhole_three_bits`, `bell_triple_sum_ge_one` | **Theorem** |
| Testfall $D_E=0$ | `test_manual_D_E_zero` | **Theorem** |
| $N_+(X)=N_-(X) \Rightarrow W_E(X)=0$ | `W_E_up_to_zero_of_balance` | **Theorem** |
| $-1 \le W_E(X) \le 1$ | `W_E_up_to_bounds` | **Theorem** |
| $N_\pm(X)$ Gleitfenster auf κ-Primfolge | `N_plus_sliding_up_to`, `N_minus_sliding_up_to` | **Theorem** |
| $D_E(X)$, $Q_E(X)$ Primvierlinge | `D_E_up_to_eq_quadruplet`, `Q_E_up_to_eq_quadruplet` | **Theorem** |
| Referenz X=1000 (Gleitfenster) | `example` (`native_decide`, 4/4) | **Theorem** |
| Referenz X=10^6 (Vierlinge) | `#eval` 84 / 82 (`PatternCount`) | **computabel** |

### Schicht B — Struktur (Definitionen)

| Objekt | Lean-Name | Status |
|--------|-----------|--------|
| Gleitfenster-Zählung | `N_plus`, `N_minus`, `countSlidingWord` | **Struktur** |
| $D_E$, $Q_E$, $\chi_{\mathrm{Hol}}$ auf Listen | `D_E`, `Q_E`, `chi_Hol` | **Struktur** |
| $R_\beta$, $\tilde D_E$ auf Listen | `R_beta`, `D_tilde_E` | **Struktur** (`noncomputable`) |
| κ-Primfolge bis $X$ | `primeEabcClassesUpTo` (`PatternCount`) | **Struktur** (computabel) |
| Primvierling-Zählung bis $X$ | `N_plus_up_to`, `N_minus_up_to`, `D_E_up_to`, `Q_E_up_to` | **Struktur** (computabel) |
| Gleitfenster κ-Folge bis $X$ | `N_plus_sliding_up_to`, `N_minus_sliding_up_to` | **Struktur** (computabel) |
| $W_E(X)$, $R_\beta(X)$, $\tilde D_E(X)$ | `W_E_up_to`, `R_beta_up_to`, `D_tilde_E_up_to` | **Struktur** |
| Bell/CHSH-Skeleton | `EabcWindowObservables`, `LocalRealismOnGE`, `chshSum` | **Struktur** |

### Schicht R — Forschungsbrücken (`sorry`)

| Brücke | Lean-Name | Status |
|--------|-----------|--------|
| $\mathrm{Hol}_E = 0$ (asymptotisch) | `Hol_E_zero` | **Brücke** (`sorry`) |
| $R_{1/2}(X) = O(1)$, Primdichte | (noch nicht formalisiert) | **Brücke** (offen) |

```bash
cd collatz_eabc_core
lake build CollatzEabc.PatternCount
lake build CollatzEabc.HolonomieFehlerterm   # 1 sorry (Hol_E_zero)
```

---

## Lean-Modul `FlussPhiE.lean`

### Schicht B — Struktur

| Objekt | Lean-Name | Status |
|--------|-----------|--------|
| Gerichtete Kanten $E^\pm$ | `C4DirectedEdge`, `E_plus`, `E_minus` | **Struktur** |
| Kantenquelle/-ziel | `edgeSrc`, `edgeTgt` | **Struktur** |
| Vorwärtszyklus (Liste) | `forwardCycleVertices` | **Struktur** |
| Kanonisches $h$ | `h_canonical`, `C4HarmonicForm` | **Struktur** |
| $C_E$, $S_E$ bis $X$ | `C_E_up_to`, `S_E_up_to` | **Struktur** (Alias zu `D_E_up_to`, `Q_E_up_to`) |
| $W_E$, $R_\beta$, $\tilde D_E$ bis $X$ | `W_E_up_to`, `R_beta_up_to`, `D_tilde_E_up_to` | **Struktur** (in `HolonomieFehlerterm`) |
| $\Phi_E$ als Grenzwert | `HasPhi_E` | **Struktur** |
| $\Phi_E = 0$ (Nullhypothese) | `Phi_E_eq_zero` | **Struktur** |
| $\alpha_E = 1$ (Skalenparameter) | `HasAlpha_E_one` | **Struktur** |
| $\Phi_E \neq 0$ (Vermutung) | `phi_E_conjecture` | **Struktur** (`Prop`) |
| Diskrete Paarung | `innerProductOmegaH`, `circulationOmega`, `OmegaE` | **Struktur** |

### Schicht A — Theorem (bewiesen)

| Aussage | Lean-Name | Status |
|---------|-----------|--------|
| $E^+ \sqcup E^-$ | `E_plus_union_E_minus` | **Theorem** |
| Zyklus geschlossen | `forward_cycle_closed` | **Theorem** |
| Harmonisches $h$ existiert | `harmonic_form_exists` | **Theorem** |
| $h \notin \mathrm{im}\,\delta$ | `h_canonical_not_coboundary` | **Theorem** |
| $N_+=N_-$ auf Liste $\Rightarrow$ $\chi_{\mathrm{Hol}}=0$ | `chi_Hol_zero_of_balance` (in `HolonomieFehlerterm`) | **Theorem** |
| $N_+(X)=N_-(X) \Rightarrow W_E(X)=0$ | `W_E_up_to_zero_of_balance` (in `HolonomieFehlerterm`) | **Theorem** |
| $-1 \le W_E(X) \le 1$ | `W_E_up_to_bounds` (in `HolonomieFehlerterm`) | **Theorem** |
| $\langle\omega,h\rangle = C_E$ (diskret) | `Phi_E_eq_inner_product_discrete` | **Theorem** |

### Schicht R — Forschungsbrücken (`sorry`)

| Brücke | Lean-Name | Status |
|--------|-----------|--------|
| Asymptotische Symmetrie $\Rightarrow \Phi_E=0$ | `Phi_E_zero_of_symmetry` | **Brücke** (`sorry`) |
| Prim-$\omega_E$ $\leftrightarrow$ asymptotische Paarung | `Phi_E_eq_inner_product` | **Brücke** (`sorry`) |
| HL-Symmetrie $\Rightarrow \mathrm{Hol}_E=0 \Rightarrow \Phi_E=0$ | `hol_E_zero_of_HL` | **Brücke** (`sorry`) |
| $\Phi_E\neq 0 \Rightarrow \alpha_E=1$ | `phi_E_ne_zero_implies_alpha_E_one` | **Brücke** (`sorry`) |
| **EABC-Vermutung** $\Phi_E \neq 0$ | `phi_E_conjecture_statement` | **Brücke** (`sorry`; **H₃**, Ebene 4) |

---

## Bewiesen vs. offen

### Schicht A — bewiesen (kombinatorisch)

1. Partition $E^+ \sqcup E^- =$ alle acht gerichteten Kanten (`E_plus_union_E_minus`).
2. Geschlossener Vorwärtszyklus (`forward_cycle_closed`).
3. Existenz und Nicht-Korand-Eigenschaft des harmonischen Generators $h$ (`harmonic_form_exists`, `h_canonical_not_coboundary`).
4. Diskrete Identität $\langle\omega_E, h\rangle = \sum_{E^+}\omega - \sum_{E^-}\omega$ (`Phi_E_eq_inner_product_discrete`).
5. Endliche Symmetrie $N_+=N_- \Rightarrow W_E=0$ (`chi_Hol_zero_of_balance`, `W_E_up_to_zero_of_balance`; Testfall `test_manual_D_E_zero`).
6. Schranken $-1 \le W_E(X) \le 1$ für alle $X$ (`W_E_up_to_bounds`) — rein kombinatorisch, unabhängig von Prim-Asymptotik.
7. $D_E$, $Q_E$, $R_\beta$, $\tilde D_E$ formal definiert (Listen + Prim-Skeleton).

### Schicht R — `sorry` (analytisch / Prim-Enumeration)

1. `N_plus_up_to`, `N_minus_up_to`, `Hol_E_zero` — Primvierling-Zählung **computabel** (`PatternCount`); Asymptotik `Hol_E_zero` offen.
2. `Phi_E_zero_of_symmetry` — Filter-Grenzwert bei asymptotisch gleichen $N_\pm$.
3. `Phi_E_eq_inner_product` — Prim-induzierte $\omega_E$ und asymptotische Paarung.
4. `hol_E_zero_of_HL` — Hardy–Littlewood-artige Symmetriehypothese $\Rightarrow \Phi_E=0$.
5. `phi_E_ne_zero_implies_alpha_E_one` — $\Phi_E \neq 0 \Rightarrow |D_E| \sim c \cdot Q$ (eine Richtung, keine Äquivalenz).
6. `phi_E_conjecture_statement` — **EABC-Vermutung** $\Phi_E \neq 0$ (**H₃**, Ebene 4).
7. `EABC_holonomy_limit_conjecture` in `HolonomyCore` — minimaler Kern.

---

## Ehrliche Einschätzung: analytisch erreichbar (Fokus $10^{10}$, Fehlerterm)

**Was Lean heute leisten kann (Schicht A):** Alle kombinatorischen Identitäten auf endlichen Folgen und Fenstern — Schranken, Symmetrie $\Rightarrow$ Null, Lückenmuster, Taubenloch, diskrete Hodge-Paarung. Das ist vollständig und unabhängig von Primzahltheorie.

**Was numerisch bis $10^{10}$ sinnvoll ist (Ebene 1–2):** $D_E(X)$, $Q_E(X)$, $R_\beta(X)$, $\tilde D_E(X)$, Fits von $\alpha_{\mathrm{loc}}$, $\alpha_E$ aus `eabc_quadruplets_fit_alpha.py`. Hier liegt der empirische Kern — Wachstumsordnung und Skalierungsexponenten, nicht der Holonomie-Grenzwert.

**Was analytisch offen und schwer ist (Schicht R):**
- Explizite $\kappa$-Folge $\to$ Gleitfenster-Zählung in Lean (braucht projektspezifische Prim-Transport-Formalisation, nicht nur `PrimeCounting`).
- Asymptotik $D_E(X) = o(Q_E(X))$ oder $D_E(X) \sim c \cdot Q_E(X)^\alpha$ — entspricht HL-artigen Primverteilungsaussagen über mod-12-Klassen; derzeit **Hypothese**, nicht Theorem.
- $\Phi_E \neq 0$ (**H₃**) ist die schwächste empirische Stütze bei großen $X$ (Fluktuationen um 0); ein Lean-Beweis wäre ein Collatz-äquivalentes Analyseproblem.

**Programm-Fokus $10^{10}$:** Primär $D_E$-Skalierung und $R_\beta$-Plateaus dokumentieren; $\Phi_E$ als Endfrage (Ebene 4) zurückstellen. Lean-Roadmap: `N_plus_up_to`/`N_minus_up_to` als **computable** `#eval`-fähige Primvierling-Zählung (`PatternCount`, Abgleich Python: $10^6 \Rightarrow 84/82$); Asymptotik bleibt `sorry` mit klarer Epistemik.

---

## Build

```bash
cd collatz_eabc_core
lake build                              # gesamtes Paket grün
lake build CollatzEabc.HolonomyCore     # 1 sorry (H₃)
lake build CollatzEabc.HolonomieFehlerterm  # 1 sorry (Hol_E_zero)
lake build CollatzEabc.FlussPhiE        # 5 sorry (Schicht R)
```

Sorries in Schicht R sind erlaubt; Ziel ist ein kompilierender Beweisrahmen mit klarer A/B/R-Trennung.

---

## Analytischer Roadmap (Stand 18. Juni 2026)

**Numerischer Kontext ($10^{10}$):** $R_{1/2}=O(1)$, $W_E\to 0$, kein Hinweis auf $\Phi_E\neq 0$ (H₃).
Der analytische Fokus liegt daher auf **Fehlerterm-Struktur** (Ebene 1–2), nicht auf `HasNonzeroHolonomyLimit`.

### Was „analytisch lösen" hier bedeutet

| Lesart | Inhalt | Lean-Status |
|--------|--------|-------------|
| **Kombinatorisch** | $G_E$, $E^\pm$, $h$, Gleitfenster, $-1\le W_E\le 1$, Symmetrie $\Rightarrow W_E=0$ | **Schicht A — bewiesen** |
| **Strukturell** | $D_E$, $Q_E$, $R_\beta$, $D̃_E$, $\Phi_E$, $\alpha_E$ als Definitionen/`Prop` | **Schicht B — formalisiert** |
| **Prim-Asymptotik** | $N_\pm(X)$ aus κ-Folge, HL-Symmetrie, $\Phi_E=0$ | **Schicht R — `sorry`** |
| **Holonomie H₃** | $\lim W_E=\Phi_E\neq 0$ | **Schicht R — `sorry`** (Ebene 4, Endfrage) |

**Ehrliche Einschätzung:** „Analytisch lösen" im Collatz-/EABC-Kontext heißt primär:
(1) die **arithmetische Zählung** $N_\pm(X)$ über Primzahlen mod 12 formalisieren,
(2) **Fehlerterm-Wachstum** $|D_E|$ vs. $Q^\beta$ mit zahlentheoretischen Werkzeugen angehen,
(3) erst danach Grenzwerte $W_E\to\Phi_E$. Schritt (3) ist numerisch derzeit **nicht** gestützt.

### Ziel-Tabelle

| Ziel | Lean-Status | Nächster Schritt |
|------|-------------|------------------|
| Kreisgraph $G_E$, $E^\pm$, harmonisches $h$ | **GREEN** (`FlussPhiE`) | — |
| $D_E$, $Q_E$, $\chi_{\mathrm{Hol}}$ auf Listen | **GREEN** (`HolonomieFehlerterm`) | — |
| $-1\le\chi_{\mathrm{Hol}}\le 1$ | **GREEN** `chi_Hol_bounds` | — |
| $N_+=N_-\Rightarrow\chi_{\mathrm{Hol}}=0$ | **GREEN** `chi_Hol_zero_of_balance` | — |
| $W_E$, $Q_E$, $R_\beta$, $D̃_E$ bis Primgrenze $X$ | **B** (`HolonomieFehlerterm` + `PatternCount`) | — |
| $-1\le W_E(X)\le 1$ | **GREEN** `W_E_up_to_bounds` | — |
| $N_+(X)=N_-(X)\Rightarrow W_E(X)=0$ | **GREEN** `W_E_up_to_zero_of_balance` | — |
| $\Phi_E=0$ (Nullhypothese) | **B** `Phi_E_eq_zero` | `Phi_E_zero_of_symmetry` (Filter) |
| $\Phi_E\neq 0$ (H₃) | **R** `phi_E_conjecture_statement` | nicht priorisieren bei $10^{10}$-Daten |
| $\Phi_E\neq 0\Rightarrow\alpha_E=1$ | **R** `phi_E_ne_zero_implies_alpha_E_one` | aus `Tendsto` + $|D_E|/Q\to|\Phi_E|$ |
| $N_+(X)$, $N_-(X)$ Primvierling-Zählung | **B** `PatternCount` (`#eval` 84/82 bei $10^6$) | Gleitfenster-Variante `N_\pm_sliding_up_to` |
| $\mathrm{Hol}_E=0$ | **R** `Hol_E_zero` | HL-artige Prim-Symmetrie |
| Lückenmuster $(2,4,2,4)$, Taubenloch | **GREEN** | — |
| `HasNonzeroHolonomyLimit` (minimaler Kern) | **R** `HolonomyCore` | getrennt halten von Fehlerterm |

### Implementiert in dieser Session (vs. Roadmap)

**Lean (`PatternCount.lean`, neu):**
- Primvierling-Zählung `N_plus_up_to`, `N_minus_up_to`, `D_E_up_to`, `Q_E_up_to` (computabel)
- `#eval` bei $10^6$: N₊=84, N₋=82 (Abgleich `eabc_quadruplets_1e10.py`)

**Lean (`HolonomieFehlerterm.lean`):**
- `Q_E`, `R_beta`, `D_tilde_E` (Listen); `W_E_up_to`, `R_beta_up_to`, `D_tilde_E_up_to` (Primgrenze)
- Gleitfenster-Variante `N_plus_sliding_up_to`, `N_minus_sliding_up_to`
- **Theorem:** `chi_Hol_bounds`, `chi_Hol_zero_of_balance`, `W_E_up_to_bounds`, `W_E_up_to_zero_of_balance`

**Lean (`FlussPhiE.lean`):**
- **Struktur:** `HasAlpha_E_one` ($|D_E|/Q\to c\neq 0$)
- **Brücke (R):** `phi_E_ne_zero_implies_alpha_E_one` — explizit `sorry`, keine Holonomie-Behauptung

**Nicht implementiert (bewusst):**
- Beweis von H₃ / `HasNonzeroHolonomyLimit`
- Analytische Abschätzung $R_{1/2}=O(1)$ (`Hol_E_zero` bleibt `sorry`)

### Mathlib-Brücken (priorisiert)

| Brücke | Mathlib-Kandidat | Ziel |
|--------|------------------|------|
| $\pi(x)$, Primzählung | `Mathlib.NumberTheory.PrimeCounting` | effizientere Enumeration (optional) |
| Dirichlet-Charaktere mod 12 | `Mathlib.NumberTheory.DirichletCharacter.*` | HL-Symmetrie / $D_E$-Fehlerterm |
| Filter / `Tendsto` | `Mathlib.Order.Filter.AtTopBot.Tendsto` | `Phi_E_zero_of_symmetry` |
| Potenzen / $Q^\beta$ | `Mathlib.Analysis.SpecialFunctions.Pow.Real` | `R_beta_up_to` (genutzt) |
| Nat.density | — | **fehlt** in Mathlib; nur `Finset`-Zählungen |

### Empfohlene Reihenfolge (analytisch)

1. **`Hol_E_zero`** — als `Prop` (`Tendsto (fun X => W_E_up_to X) atTop (nhds 0)`).
2. **`Phi_E_zero_of_symmetry`** — aus `W_E_up_to_zero_of_balance` + `Filter.Eventually` + `Tendsto_const`.
3. **Fehlerterm-Schranken** — $|D_E|\ll Q$ oder $R_{1/2}=O(1)$ als `Prop`-Ziel (Schicht R).
4. **`phi_E_ne_zero_implies_alpha_E_one`** — nur wenn H₃ relevant wird.
5. **H₃ / `HasNonzeroHolonomyLimit`** — Ebene 4, nach geklärter $D_E$-Skalierung.

---

## Nächste Schritte

1. `Hol_E_zero` als explizite `Prop` (`Tendsto (fun X => W_E_up_to X) atTop (nhds 0)`).
2. `Phi_E_zero_of_symmetry` aus `W_E_up_to_zero_of_balance` + Filter-`Eventually`.
3. Numerische Brücke Lean↔Python für $D_E$, $Q_E$, $R_{1/2}$ bis $10^{10}$ (CSV-Checkpoint-Vergleich).
4. Optional: effizientere Primvierling-Enumeration (Sieb, wie Python).
5. Optional: `Mathlib.Combinatorics.SimpleGraph` für $G_E$; magnetischer Laplace $L_{\mathrm{mag}}$.
