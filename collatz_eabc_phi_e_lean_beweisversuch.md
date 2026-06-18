# EABC Φ_E — analytischer Lean-Beweisversuch

**Status:** Skeleton (PR #59, `collatz/eabc-05-holonomie-fehlerterm`)  
**Lean:** `collatz_eabc_core/CollatzEabc/FlussPhiE.lean`  
**Kanonische Geometrie:** `collatz_eabc_diskrete_geometrie.md`  
**Numerik:** `collatz_eabc_hodge_eabc.py` (`Phi_E`, `flux_density_limit`, `inner_product_omega_h`)

---

## Programm

EABC = priminduzierter Fluss auf $C_4 \cong S^1$ mit

- $G_E = (V,E)$, $V=\{E,A,B,C\}$
- $E^+ = \{EA, AB, BC, CE\}$, $E^- = \{EC, CB, BA, AE\}$
- $C_E(X) = N_+(X) - N_-(X)$, $S_E(X) = N_+(X) + N_-(X)$
- $W_E(X) = C_E(X)/S_E(X)$, $\Phi_E = \lim_{X\to\infty} W_E(X)$

**Vermutung:** $\Phi_E \neq 0$ bzw. $\langle\omega_E, h\rangle \neq 0$ für kanonisches harmonisches $h$.

---

## Lean-Modul `FlussPhiE.lean`

| Objekt | Lean-Name | Status |
|--------|-----------|--------|
| Gerichtete Kanten $E^\pm$ | `C4DirectedEdge`, `E_plus`, `E_minus` | **Definition** |
| Vorwärtszyklus | `forwardCycleVertices`, `forward_cycle_closed` | **Theorem** |
| Kanonisches $h$ | `h_canonical`, `harmonic_form_exists` | **Theorem** |
| $h \notin \mathrm{im}\,\delta$ | `h_canonical_not_coboundary` | **Theorem** |
| $C_E$, $S_E$, $W_E$ bis $X$ | `C_E_up_to`, `S_E_up_to`, `W_E_up_to` | **Definition** (Primteil `sorry` via `HolonomieFehlerterm`) |
| $\Phi_E = 0$ | `Phi_E_eq_zero`, `HasPhi_E` | **Definition** |
| $\Phi_E \neq 0$ | `phi_E_conjecture` | **Vermutung** (`phi_E_conjecture_statement` = `sorry`) |
| $\langle\omega,h\rangle = C_E$ (diskret) | `innerProductOmegaH`, `Phi_E_eq_inner_product_discrete` | **Theorem** |
| $N_+=N_-$ auf Liste $\Rightarrow$ $\chi_{\mathrm{Hol}}=0$ | `chi_Hol_zero_of_balance` | **Theorem** |
| $N_+(X)=N_-(X) \Rightarrow W_E(X)=0$ | `W_E_up_to_zero_of_balance` | **Theorem** |
| asymptotische Symmetrie $\Rightarrow \Phi_E=0$ | `Phi_E_zero_of_symmetry` | **Skeleton** (`sorry`) |
| Prim-Brücke $\Phi_E \leftrightarrow \langle\omega,h\rangle$ | `Phi_E_eq_inner_product` | **Skeleton** (`sorry`) |
| HL-Symmetrie $\Rightarrow \mathrm{Hol}_E=0$ | `hol_E_zero_of_HL` | **Skeleton** (`sorry`) |

---

## Bewiesen vs. offen

### Bewiesen (kombinatorisch)

1. Partition $E^+ \sqcup E^- =$ alle acht gerichteten Kanten (`E_plus_union_E_minus`).
2. Existenz und Nicht-Korand-Eigenschaft des harmonischen Generators $h$ (`harmonic_form_exists`, `h_canonical_not_coboundary`).
3. Diskrete Identität $\langle\omega_E, h\rangle = \sum_{E^+}\omega - \sum_{E^-}\omega$ (`Phi_E_eq_inner_product_discrete`).
4. Endliche Symmetrie $N_+=N_- \Rightarrow W_E=0$ auf Listen (`chi_Hol_zero_of_balance`; Testfall in `HolonomieFehlerterm.test_manual_D_E_zero`).

### `sorry` (analytisch / Prim-Enumeration)

1. `N_plus_up_to`, `N_minus_up_to`, `Hol_E_zero` — in `HolonomieFehlerterm.lean`.
2. `Phi_E_zero_of_symmetry` — Filter-Grenzwert bei asymptotisch gleichen $N_\pm$.
3. `Phi_E_eq_inner_product` — Prim-induzierte $\omega_E$ und asymptotische Paarung.
4. `hol_E_zero_of_HL` — Hardy–Littlewood-artige Symmetriehypothese $\Rightarrow \Phi_E=0$.
5. `phi_E_conjecture_statement` — EABC-Vermutung $\Phi_E \neq 0$.

---

## Build

```bash
cd collatz_eabc_core
lake build CollatzEabc.FlussPhiE
```

Sorries sind erlaubt; Ziel ist ein kompilierender Beweisrahmen.

---

## Nächste Schritte

1. `N_plus_up_to` / `N_minus_up_to` aus κ-Folge + `PrimeCounting` definieren.
2. `Phi_E_zero_of_symmetry` aus `W_E_up_to_zero_of_balance` + Filter-`Eventually`.
3. Prim-$\omega_E$ aus Gleitfenster-Zählung; `Phi_E_eq_inner_product` als Grenzwertbrücke.
4. Optional: `Mathlib.Combinatorics.SimpleGraph` für $G_E$; magnetischer Laplace $L_{\mathrm{mag}}$ (vgl. `collatz_mathlib_eabc_kandidaten.md`).
