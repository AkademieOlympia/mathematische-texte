# Mathlib-Inventar: EABC-Ansatz und Riemann-Kugel

**Stand:** 16. Juni 2026  
**Mathlib:** v4.29.0 (Lean `leanprover/lean4:v4.29.0`), lokal via  
`ptolemaeus-lean/Ptolemaeus/.lake/packages/mathlib` (Rev. `8a178386ff`)  
**Methode:** Repo-Grep + semantische Suche + `rg` über lokales Mathlib-Paket  
**Hinweis:** Kein Collatz-Beweis — reines Werkzeug-Inventar.

---

## Repo-Kontext (Kurzüberblick)

| Thema | Wichtige Projektdateien |
|-------|-------------------------|
| mod-12 EABC, Transfer-Operator | `collatz_schlussartikel_arxiv.tex`, `collatz_eabc_spektrum.tex`, `collatz_mixing_test.py` |
| Lean (ℤ₂, Uniformität) | `collatz_z2_attraktor.lean`, `collatz_uniformity.lean`, `collatz_density_appendix.lean` |
| Lean Holonomie / Bell / $D_E$ (PR #54) | `collatz_eabc_core/CollatzEabc/HolonomieFehlerterm.lean` — Lücken $(2,4,2,4)$ bewiesen; Prime/CHSH `sorry` |
| Lean $\Phi_E$, harmonisches $h$, $W_E$ (PR #59) | `collatz_eabc_core/CollatzEabc/FlussPhiE.lean` — $C_4$, $E^\pm$, $h$ bewiesen; Prim-Grenzwert `sorry` |
| Bell / $D_E$ Theorie + Numerik | `collatz_eabc_bell_holonomie.md`, `collatz_eabc_fehlerterm_hypothese.md`, `collatz_eabc_bell_inequality_test.py` |
| Hurwitz-Polytop / 24 Einheiten | `collatz_hurwitz_polytop_eabc.tex`, `Rechtsorbit.py`, `Quantencomputer.py` |
| Riemann-Kugel / SL(2) / Möbius | `Edinburg.py`, `collatz_sinh_hyperbel.tex`, `Arithmetik zu Topologie.tex` |
| Mathlib-Update (Vorgänger) | `collatz_lean_mathlib_update.tex`, `collatz_offene_punkte.md` |

Numerische Kernwerte (heuristisch, nicht formalisiert): $|\lambda_2|\approx 0{,}35$, $\Lambda\approx -0{,}830$, Eigenwerte $\pm 2\cos(\pi/8)$, $\pm 2\cos(3\pi/8)$ am 4×4-Operator.

---

## §1 Bereits genutzt im Projekt

| Mathlib-Modul / Theorem | Relevanz | Priorität |
|-------------------------|----------|-----------|
| `Mathlib.NumberTheory.Padics.PadicVal` (`padicValNat.mul/pow/div_of_dvd`) | 2-adische Valuation, LTE-Worst-Familie $2^{k+1}3^r-1$ | — (aktiv) |
| `Mathlib.NumberTheory.Padics.PadicIntegers` | $\mathbb{Z}_2$, `dist2`, `ExceptionSet` | — (aktiv) |
| `Mathlib.Topology.MetricSpace.Ultra.Basic` | Ultrametrik von $\|\cdot\|_2$ | — (aktiv) |
| `Mathlib.Topology.Basic` / `Closure` | Hülle der Ausnahmemenge $E$ | — (aktiv) |
| `Mathlib.Data.Nat.Parity` / `Nat.ModEq` | Odd-Klassen, $3^j\equiv 3\pmod 8$ | — (aktiv) |
| `Mathlib.Topology.Algebra.InfiniteSum.Basic` | $(1-p)^n\to 0$, Mischschranken | — (aktiv) |
| `Mathlib.Data.Nat.Factorization.Basic` (`Nat.card_multiples`) | C-Ketten-Dichte $N/2^k$ | — (aktiv) |
| `Mathlib.Analysis.SpecificLimits.Basic/Normed` | Geometrische Reihe, `tail_series_formula` | — (aktiv) |

**Lean-Dateien mit diesen Imports:** `collatz_uniformity.lean`, `collatz_density_appendix.lean`, `collatz_z2_attraktor.lean`.

---

## §2 Neue Kandidaten für EABC (mod-12, Operator, Mixing)

| # | Mathlib-Name | 1-Satz-Relevanz | Priorität |
|---|--------------|-----------------|-----------|
| 1 | `Mathlib.LinearAlgebra.Matrix.Irreducible.Defs` — `Matrix.IsIrreducible`, `Matrix.IsPrimitive`, `Matrix.toQuiver` | Graph-theoretische Irreduzibilität nichtnegativer Matrizen; direkte Sprache für mod-12-RPF-Operator (4×4, endlich). | **hoch** |
| 2 | `Mathlib.LinearAlgebra.Matrix.Gershgorin` — `eigenvalue_mem_ball` | Obere Schranken für Eigenwerte der 4×4-Transfermatrix; Hilfsmittel für $|\lambda_2|<1$ ohne vollen Perron-Frobenius. | **hoch** |
| 3 | `Mathlib.LinearAlgebra.Matrix.Charpoly.Eigs` — `mem_spectrum_iff_isRoot_charpoly` | Verbindet Spektrum der rationalen mod-12-Matrix mit Nullstellen des Charpoly (exakte Rechnung in $\mathbb{Q}$). | **hoch** |
| 4 | `Mathlib.LinearAlgebra.Eigenspace.Basic` / `Matrix` | Eigenräume, `HasEigenvalue`; für die vier Eigenwerte $\pm 2\cos(\pi/8)$, $\pm 2\cos(3\pi/8)$ nach numerischer Bestätigung. | **mittel** |
| 5 | `Mathlib.Probability.Kernel.IonescuTulcea.Traj` — `traj`, `trajMeasure` | Markov-Ketten-Infrastruktur: Collatz als Übergangskern auf $\mathbb{Z}/2^m\mathbb{Z}$ oder EABC-Zuständen. | **hoch** |
| 6 | `Mathlib.Probability.Kernel.Irreducible` — `Kernel.IsIrreducible` | Irreduzibilität stochastischer Kerne; probabilistisches Pendant zu `Matrix.IsIrreducible`. | **mittel** |
| 7 | `Mathlib.Probability.Distributions.Geometric` — `geometricPMFRealSum` | Geometrische Verteilung für Wartezeiten bis C-Kette / Ausnahme; passt zu $(1-p)^n\to 0$. | **mittel** |
| 8 | `Mathlib.Data.Nat.Multiplicity` — `emultiplicity_pow/mul` | Volles LTE-Gerüst für $p=2$, $a=3$; nächster Schritt nach `padicValNat`. | **hoch** |
| 9 | `Mathlib.Data.ZMod.Basic` | Formale mod-$n$-Arithmetik für EABC-Klassen und Übergänge auf `ZMod 12`. | **hoch** |
| 10 | `Mathlib.Combinatorics.SimpleGraph.Basic` / `Finite` | Endliche Graphen für mod-12-Übergänge; Brücke zu `Matrix.toQuiver`. | **mittel** |
| 11 | `Mathlib.Dynamics.Ergodic.Ergodic` — `Ergodic`, `PreErgodic` | Definition maßerhaltender ergodischer Systeme; Voraussetzung für Drift-Argumente. | **mittel** |
| 12 | `Mathlib.Analysis.InnerProductSpace.MeanErgodic` — `tendsto_birkhoffAverage_orthogonalProjection` | Von-Neumann-Mittelwertsatz ($L^2$); einziger rigoroser Ergodensatz, nur indirekt für Collatz. | **niedrig** |
| 13 | `Mathlib.Dynamics.BirkhoffSum.Basic` | Birkhoff-Summen $S_n f = \sum_{k<n} f\circ T^k$; algebraisches Gerüst für Drift, nicht punktweise Konvergenz. | **mittel** |
| 14 | `Mathlib.NumberTheory.LSeries.HurwitzZetaValues` — `riemannZeta_two_mul_nat` | $\zeta(2n)$ via Bernoulli; Verbindung zu Bernoulli-Normschalen (nicht Collatz-Kern). | **niedrig** |
| 15 | `Mathlib.Analysis.SpecialFunctions.Stirling` — `log_stirlingSeq_diff_le` | Robbins-Schranke für asymptotische Schalenanalyse. | **niedrig** |
| 16 | `Mathlib.Analysis.SpecificLimits.Normed` — `tsum_choose_mul_geometric_of_norm_lt_one` | Verallgemeinerte Tail-Reihe für Mischschranken mit beliebigem $r<1$. | **mittel** |
| 17 | `Mathlib.Probability.StrongLaw` | Starkes Gesetz der großen Zahlen; nur für i.i.d.-Modelle, nicht für deterministische Collatz-Bahnen. | **niedrig** |

### Empfohlene erste Lean-Schritte (EABC-Operator)

1. `ZMod 12` + explizite 4×4-Matrix `M : Matrix (Fin 4) (Fin 4) ℚ` aus EABC-Übergängen.
2. `Matrix.IsIrreducible M` via `toQuiver` + Pfadexistenz (reine `Finset`-Kombinatorik).
3. `Gershgorin` / `Charpoly.Eigs` für $|\lambda_2|<1$.
4. Optional: `Probability.Kernel` für stochastische Lesart.

### Holonomie-Fehlerterm / Bell (PR #54, Phase 1–2)

| Ziel | Mathlib-Kandidat | Status im Projekt |
|------|------------------|-------------------|
| Primzählung bis $X$ | `Mathlib.NumberTheory.PrimeCounting` | **offen** — `N_plus_up_to` in `HolonomieFehlerterm.lean` = `sorry` |
| Dirichlet-Charaktere mod $12$ | `Mathlib.NumberTheory.DirichletCharacter.*` | **Hypothese** in `collatz_eabc_fehlerterm_hypothese.md`, nicht formalisiert |
| Bell/CHSH / Chebyshev-Bias | — | **nicht in Mathlib**; kombinatorisches Taubenloch in `HolonomieFehlerterm.lean` (**Theorem**) |
| $\Phi_E$, $W_E$, harmonisches $h$ auf $C_4$ | `Mathlib.Combinatorics.SimpleGraph.Basic` (optional) | **teilweise** — `FlussPhiE.lean`: $h$, $\langle\omega,h\rangle$ diskret **bewiesen**; $\Phi_E$-Grenzwert **`sorry`** |

Siehe `collatz_eabc_phi_e_lean_beweisversuch.md`, `collatz_eabc_diskrete_geometrie.md` §2 und `collatz_eabc_bell_holonomie.md` §12.

---

## §3 Riemann-Kugel / Hurwitz / Quaternionen

| # | Mathlib-Name | 1-Satz-Relevanz | Priorität |
|---|--------------|-----------------|-----------|
| 1 | `Mathlib.Topology.Compactification.OnePoint.ProjectiveLine` — `OnePoint.equivProjectivization` | $\mathrm{OnePoint}\,K \simeq \mathbb{P}^1(K)$ via $GL(2)$-Wirkung; algebraische Riemann-Kugel als Projektivlinie. | **hoch** |
| 2 | `Mathlib.Topology.Compactification.OnePoint.Sphere` — `onePointEquivSphereOfFinrankEq` | Homöomorphie Einpunktkompaktifikation $\leftrightarrow$ Sphäre (Stereographie). | **hoch** |
| 3 | `Mathlib.Geometry.Manifold.Instances.Sphere` — `stereographic` | Stereographische Projektion, analytische Sphären-Mannigfaltigkeit; „Bamberg-Kugel = S²“. | **hoch** |
| 4 | `Mathlib.Analysis.Complex.UpperHalfPlane.MoebiusAction` | $GL(2,\mathbb{R})$-Wirkung als Möbius-Transformationen; Collatz-Schritte als $SL(2)$ in `collatz_sinh_hyperbel.tex`. | **hoch** |
| 5 | `Mathlib.LinearAlgebra.Matrix.SpecialLinearGroup` / `ProjectiveSpecialLinearGroup` | $SL(2)$, $PSL(2)$; Kanal für hyperbolische Collatz-Geometrie. | **mittel** |
| 6 | `Mathlib.LinearAlgebra.Projectivization.Action` | $GL(V)$-Wirkung auf $\mathbb{P}(V)$; formale Basis für Projektionszeuge-Diagramm. | **mittel** |
| 7 | `Mathlib.Algebra.Quaternion` — `Quaternion`, `normSq` | Quaternionen-Algebra $\mathbb{H}$ über beliebigem Ring; Rechenkern für Hurwitz-Einheiten-Skripte. | **hoch** |
| 8 | `Mathlib.Analysis.Quaternion` / `Normed.Algebra.QuaternionExponential` | Analytische Struktur auf $\mathbb{H}$; für Normschalen und Exponentialabbildung. | **mittel** |
| 9 | `Mathlib.GroupTheory.SpecificGroups.Quaternion` | Endliche Quaternionengruppe $Q_8$; **nicht** dasselbe wie Hurwitz-Gitter — Namenskonflikt beachten. | **niedrig** |
| 10 | `Mathlib.NumberTheory.ModularForms.Basic` / `SlashActions` | Modulformen, $SL(2,\mathbb{Z})$-Wirkung; Selberg-/Eisenstein-Kontext in Hurwitz-Polytop-Text. | **mittel** |
| 11 | `Mathlib.NumberTheory.LSeries.HurwitzZeta` / `HurwitzZetaValues` | Hurwitz-Zeta $\zeta(s,a)$; **analytische** Hurwitz-Funktion, **nicht** Hurwitz-Quaternionen-Gitter. | **niedrig** (Namensfalle) |
| 12 | `Mathlib.NumberTheory.Padics.Complex` | 2-adische Komplexe; mögliche Brücke $\mathbb{Z}_2 \leftrightarrow$ komplexe Dynamik. | **niedrig** |

### Wichtige Unterscheidung: zwei „Hurwitz“-Bedeutungen

| Begriff im Projekt | Mathlib-Entsprechung | Status |
|--------------------|----------------------|--------|
| Hurwitz-**Zeta** / Bernoulli-Schalen | `LSeries.HurwitzZeta*` | vorhanden |
| Hurwitz-**Quaternionen** / 24 Einheiten / $H_4$ | — | **fehlt** (s. §4) |
| Satz von Hurwitz (normierte Divisionsalgebren) | — | **fehlt** |
| Oktonionen $\mathbb{O}$ | — | **fehlt** |

### Repo-Verbindung Riemann-Kugel ↔ EABC

- `Edinburg.py`: explizit „Bamberg-Kugel = S² = Riemann-Kugel“.
- `collatz_hurwitz_polytop_eabc.tex`: Projektionsdiagramm Gauß/Eisenstein/Hurwitz ↔ EABC-Kanäle.
- Mathlib deckt **Sphäre, Stereographie, $\mathbb{P}^1$, Möbius, Quaternionen-Algebra** ab.
- Die **24 Hurwitz-Einheiten** und **Divisionsalgebren-Satz** müssen projektspezifisch formalisiert werden.

---

## §4 Fehlende Lücken (bestätigt nicht in Mathlib)

| Thema | Gesucht | Ergebnis |
|-------|---------|----------|
| Vollständiger **Perron-Frobenius-Satz** | `Perron`, Spektralradius primärer Eigenwert | Nur `IsIrreducible`/`IsPrimitive`-**Definitionen**; kein Hauptsatz |
| **Kingman**-Subadditivität | `Kingman`, `subadditive ergodic` | Nicht vorhanden |
| **Punktweiser Birkhoff-Satz** | `Birkhoff` + `ae` vs. punktweise | Nur Von-Neumann ($L^2$) und Birkhoff-**Summen** |
| **Collatz-Ergodizität** | — | Nicht vorhanden |
| **von-Staudt-Clausen** | `vonStaudt` | Nicht in v4.29 (`rg` leer) |
| Asymptotische **Nat.density** | `Nat.density` | Nicht vorhanden; nur `Finset`-Zählungen |
| **Hurwitz-Ganzzahlen** / Quaternionen-Gitter | `HurwitzInteger` | Nicht vorhanden |
| **Oktonionen** / Cayley-Dickson | `Octonion` | Nicht vorhanden |
| **Satz von Hurwitz** (Divisionsalgebren $\mathbb{R},\mathbb{C},\mathbb{H},\mathbb{O}$) | — | Nicht vorhanden |
| Dedizierte **RiemannSphere**-Definition | `RiemannSphere` | Nicht als eigenständiger Typ; äquivalent via `OnePoint` + `Sphere` |
| **AdjacencyMatrix** für `SimpleGraph` | `adjacencyMatrix` | Kein dediziertes Modul gefunden |
| mod-12-**EABC-Klassifikation** | — | Projektspezifisch |
| Transfer-Operator-Eigenwerte $\pm 2\cos(\pi/8)$ | — | Projektspezifisch (numerisch in Python) |

---

## §5 Priorisierte Import-Roadmap für Lean

### Phase A — mod-12-Operator (höchste Priorität)

```lean
import Mathlib.Data.ZMod.Basic
import Mathlib.Data.Matrix.Basic
import Mathlib.LinearAlgebra.Matrix.Irreducible.Defs
import Mathlib.LinearAlgebra.Matrix.Gershgorin
import Mathlib.LinearAlgebra.Matrix.Charpoly.Eigs
```

**Ziel:** `theorem mod12_matrix_irreducible : Matrix.IsIrreducible M` und `theorem lambda2_lt_one : |λ₂| < 1` (nach Konstruktion von `M`).

### Phase B — LTE und Padics (parallel zu Phase A)

```lean
import Mathlib.Data.Nat.Multiplicity
import Mathlib.NumberTheory.Padics.PadicVal.Basic  -- bereits genutzt
```

**Ziel:** Vollständiges LTE-Lemma `padicVal_two_three_pow_sub_one`.

### Phase C — Markov / Mischung

```lean
import Mathlib.Probability.Distributions.Geometric
import Mathlib.Probability.Kernel.IonescuTulcea.Traj
import Mathlib.Probability.Kernel.Irreducible
```

**Ziel:** Stochastische Lesart der C-Ketten-Wartezeiten; optional Markov-Kern auf `Fin 4`.

### Phase D — ℤ₂-Attraktor (Fortsetzung bestehender Arbeit)

```lean
-- bereits in collatz_z2_attraktor.lean:
import Mathlib.NumberTheory.Padics.PadicIntegers
import Mathlib.Topology.MetricSpace.Ultra.Basic
```

**Ziel:** Stufe E (`dist₂(U^k(n),E) → 0`) — Mathlib endet vor punktweiser Uniformität.

### Phase E — Riemann-Kugel / Projektion (geometrischer Zweig)

```lean
import Mathlib.Topology.Compactification.OnePoint.ProjectiveLine
import Mathlib.Geometry.Manifold.Instances.Sphere
import Mathlib.Analysis.Complex.UpperHalfPlane.MoebiusAction
import Mathlib.Algebra.Quaternion
import Mathlib.LinearAlgebra.Projectivization.Action
```

**Ziel:** Formale Brücke Stereographie $\leftrightarrow$ $\mathbb{P}^1$; Quaternionen-Multiplikation für Hurwitz-Einheiten (24 Einheiten **neu definieren**).

### Phase F — Ergodentheorie (nur indirekt)

```lean
import Mathlib.Dynamics.Ergodic.Ergodic
import Mathlib.Analysis.InnerProductSpace.MeanErgodic
```

**Ziel:** Hilfsmittel für Hilbertraum-Modelle; **nicht** für punktweisen Collatz-Beweis ausreichend.

---

## Top-10 Kandidaten (Kurzliste)

| Rang | Mathlib-Modul | Warum |
|------|---------------|-------|
| 1 | `LinearAlgebra.Matrix.Irreducible.Defs` | Direktes Werkzeug für mod-12-RPF |
| 2 | `Data.ZMod.Basic` | Formale mod-12-EABC-Klassen |
| 3 | `LinearAlgebra.Matrix.Gershgorin` | Spektralschranke $|\lambda_2|<1$ |
| 4 | `Data.Nat.Multiplicity` | Volles LTE $p=2$ |
| 5 | `Probability.Kernel.IonescuTulcea.Traj` | Markov-Ketten-Maße |
| 6 | `LinearAlgebra.Matrix.Charpoly.Eigs` | Exakte Eigenwerte der 4×4-Matrix |
| 7 | `Topology.Compactification.OnePoint.ProjectiveLine` | Riemann-Kugel = $\mathbb{P}^1$ |
| 8 | `Geometry.Manifold.Instances.Sphere` | Stereographische Projektion |
| 9 | `Algebra.Quaternion` | Rechenkern Hurwitz-Einheiten |
| 10 | `Analysis.Complex.UpperHalfPlane.MoebiusAction` | SL(2)/Möbius-Collatz-Geometrie |

---

## Referenzen

- [Mathlib4 Docs](https://leanprover-community.github.io/mathlib4_docs/)
- Projekt: `collatz_lean_mathlib_update.tex`, `collatz_offene_punkte.md`
- Build: `cd ptolemaeus-lean/Ptolemaeus && lake env lean ../../collatz_z2_attraktor.lean`
