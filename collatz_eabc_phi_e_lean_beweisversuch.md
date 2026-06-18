# EABC Φ_E — analytischer Lean-Beweisversuch

**Status:** Skeleton (PR #63, `collatz/eabc-h03-diskrete-geometrie-fluss`; Basis PR #59)  
**Lean (Hodge-Layer):** `collatz_eabc_core/CollatzEabc/FlussPhiE.lean`  
**Lean (minimaler Kern):** `collatz_eabc_core/CollatzEabc/HolonomyCore.lean` (`EABC`-Namespace)  
**Epistemik:** `collatz_eabc_epistemik_schichten.md` — **A / B / R / C**  
**Kanonische Geometrie:** `collatz_eabc_diskrete_geometrie.md`  
**Numerik:** `collatz_eabc_hodge_eabc.py` (`Phi_E`, `flux_density_limit`, `inner_product_omega_h`)

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

EABC = priminduzierter Fluss auf $C_4 \cong S^1$ mit

- $G_E = (V,E)$, $V=\{E,A,B,C\}$
- $E^+ = \{EA, AB, BC, CE\}$, $E^- = \{EC, CB, BA, AE\}$
- $C_E(X) = N_+(X) - N_-(X)$, $S_E(X) = N_+(X) + N_-(X)$
- $W_E(X) = C_E(X)/S_E(X)$, $\Phi_E = \lim_{X\to\infty} W_E(X)$

**Vermutung (Schicht R):** $\Phi_E \neq 0$ bzw. $\langle\omega_E, h\rangle \neq 0$ für kanonisches harmonisches $h$. Dies ist **H₃** (starke Holonomie: $R_1(X)=W_E(X)\to\Phi_E\neq 0$, äquivalent $D\sim\Phi_E Q$) und entspricht der Lean-**RED**-Vermutung `HasNonzeroHolonomyLimit`. Vorgelagert (Ebenen 0–3, `collatz_eabc_zirkulationshypothese.md` §4.2): H₀a ($R_{1/2}=Z_E$ beschränkt), H₀b ($W_E\to 0$), H₁ ($R_{1/2}\to\infty$, empirisch), H₂ ($|D|\asymp Q^\alpha$, $\alpha>\tfrac{1}{2}$, asymptotisch).

---

## Lean-Modul `HolonomyCore.lean` (sauberer Split)

**GREEN/RED im Modul:** GREEN = `Node` … `W_E_bounds` (bewiesen); RED = `HasNonzeroHolonomyLimit`, `EABC_holonomy_limit_conjecture` (`sorry`) — formal **H₃** (Holonomie).

Minimale Architektur im Namespace `EABC` — ohne Prim-Enumeration, ohne Hodge-Layer:

| Objekt | Lean-Name | Status |
|--------|-----------|--------|
| Kreisgraph-Knoten | `Node`, `next`, `prev` | **Definition** |
| Zykluszählung | `CycleCounts`, `circulation`, `size` | **Definition** |
| $W_E$ auf endlicher Stichprobe | `phiApprox`, `W_E` | **Definition** |
| $-1 \le W_E \le 1$ | `phiApprox_bounds`, `W_E_bounds` | **Theorem** (bewiesen) |
| Fluss bis Schranke $X$ | `EABCFlow`, `C_E`, `S_E` | **Definition** |
| $\Phi_E \neq 0$ (asymptotisch) | `HasNonzeroHolonomyLimit`, `EABC_holonomy_limit_conjecture` | **Vermutung** (`sorry`; **H₃**) |

**Grenze Beweis vs. Vermutung:** `W_E(X)=\frac{N_+(X)-N_-(X)}{N_+(X)+N_-(X)}$ liegt formal in $[-1,1]$;
die asymptotische Aussage $\lim_{X\to\infty} W_E(X)=\Phi_E\neq 0$ steht in `EABC_holonomy_limit_conjecture` (`Tendsto` in $\mathbb{R}$), nicht als exakte rationale Endkonstante.
Die zu starke Variante (eventuell konstantes rationales $\Phi$) ist nur als auskommentiertes Scaffold erhalten.

`FlussPhiE.lean` ergänzt denselben Fluss um C₄-Kanten, harmonisches $h$, Prim-Brücken (`HolonomieFehlerterm`).

```bash
cd collatz_eabc_core
lake build CollatzEabc.HolonomyCore
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
| $C_E$, $S_E$, $W_E$ bis $X$ | `C_E_up_to`, `S_E_up_to`, `W_E_up_to` | **Struktur** (Primteil `sorry` in `HolonomieFehlerterm`) |
| $\Phi_E$ als Grenzwert | `HasPhi_E` | **Struktur** |
| $\Phi_E = 0$ (Hypothese) | `Phi_E_eq_zero` | **Struktur** |
| $\Phi_E \neq 0$ (Vermutung) | `phi_E_conjecture` | **Struktur** (`Prop`) |
| Diskrete Paarung | `innerProductOmegaH`, `circulationOmega`, `OmegaE` | **Struktur** |

### Schicht A — Theorem (bewiesen)

| Aussage | Lean-Name | Status |
|---------|-----------|--------|
| $E^+ \sqcup E^-$ | `E_plus_union_E_minus` | **Theorem** |
| Zyklus geschlossen | `forward_cycle_closed` | **Theorem** |
| Harmonisches $h$ existiert | `harmonic_form_exists` | **Theorem** |
| $h \notin \mathrm{im}\,\delta$ | `h_canonical_not_coboundary` | **Theorem** |
| $N_+=N_-$ auf Liste $\Rightarrow$ $\chi_{\mathrm{Hol}}=0$ | `chi_Hol_zero_of_balance` | **Theorem** |
| $N_+(X)=N_-(X) \Rightarrow W_E(X)=0$ | `W_E_up_to_zero_of_balance` | **Theorem** |
| $\langle\omega,h\rangle = C_E$ (diskret) | `Phi_E_eq_inner_product_discrete` | **Theorem** |

### Schicht R — Forschungsbrücken (`sorry`)

| Brücke | Lean-Name | Status |
|--------|-----------|--------|
| Asymptotische Symmetrie $\Rightarrow \Phi_E=0$ | `Phi_E_zero_of_symmetry` | **Brücke** (`sorry`) |
| Prim-$\omega_E$ $\leftrightarrow$ asymptotische Paarung | `Phi_E_eq_inner_product` | **Brücke** (`sorry`) |
| HL-Symmetrie $\Rightarrow \mathrm{Hol}_E=0 \Rightarrow \Phi_E=0$ | `hol_E_zero_of_HL` | **Brücke** (`sorry`) |
| **EABC-Vermutung** $\Phi_E \neq 0$ | `phi_E_conjecture_statement` | **Brücke** (`sorry`; **H₃**) |

### Abhängigkeit `HolonomieFehlerterm.lean` (Schicht R)

| Objekt | Lean-Name | Status |
|--------|-----------|--------|
| $N_+(X)$, $N_-(X)$ auf κ-Folge | `N_plus_up_to`, `N_minus_up_to` | **Brücke** (`sorry`) |
| $\mathrm{Hol}_E = 0$ | `Hol_E_zero` | **Brücke** (`sorry`) |

---

## Bewiesen vs. offen

### Schicht A — bewiesen (kombinatorisch)

1. Partition $E^+ \sqcup E^- =$ alle acht gerichteten Kanten (`E_plus_union_E_minus`).
2. Geschlossener Vorwärtszyklus (`forward_cycle_closed`).
3. Existenz und Nicht-Korand-Eigenschaft des harmonischen Generators $h$ (`harmonic_form_exists`, `h_canonical_not_coboundary`).
4. Diskrete Identität $\langle\omega_E, h\rangle = \sum_{E^+}\omega - \sum_{E^-}\omega$ (`Phi_E_eq_inner_product_discrete`).
5. Endliche Symmetrie $N_+=N_- \Rightarrow W_E=0$ (`chi_Hol_zero_of_balance`, `W_E_up_to_zero_of_balance`; Testfall in `HolonomieFehlerterm.test_manual_D_E_zero`).

### Schicht R — `sorry` (analytisch / Prim-Enumeration)

1. `N_plus_up_to`, `N_minus_up_to`, `Hol_E_zero` — in `HolonomieFehlerterm.lean`.
2. `Phi_E_zero_of_symmetry` — Filter-Grenzwert bei asymptotisch gleichen $N_\pm$.
3. `Phi_E_eq_inner_product` — Prim-induzierte $\omega_E$ und asymptotische Paarung.
4. `hol_E_zero_of_HL` — Hardy–Littlewood-artige Symmetriehypothese $\Rightarrow \Phi_E=0$.
5. `phi_E_conjecture_statement` — **EABC-Vermutung** $\Phi_E \neq 0$ (**H₃**).

---

## Build

```bash
cd collatz_eabc_core
lake build CollatzEabc.HolonomyCore   # minimaler Kern (1 sorry)
lake build CollatzEabc.FlussPhiE      # Hodge-Layer (4 sorry)
```

Sorries in Schicht R sind erlaubt; Ziel ist ein kompilierender Beweisrahmen mit klarer A/B/R-Trennung.

---

## Nächste Schritte

1. `N_plus_up_to` / `N_minus_up_to` aus κ-Folge + `PrimeCounting` definieren (Schicht R).
2. `Phi_E_zero_of_symmetry` aus `W_E_up_to_zero_of_balance` + Filter-`Eventually` (Schicht R).
3. Prim-$\omega_E$ aus Gleitfenster-Zählung; `Phi_E_eq_inner_product` als Grenzwertbrücke (Schicht R).
4. Optional: `Mathlib.Combinatorics.SimpleGraph` für $G_E$; magnetischer Laplace $L_{\mathrm{mag}}$ (vgl. `collatz_mathlib_eabc_kandidaten.md`).
