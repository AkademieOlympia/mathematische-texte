/-
  CollatzEabc.FlussPhiE — Φ_E, W_E, harmonische 1-Form auf C₄ ≅ S¹.

  **Nicht sorry-frei.** Epistemische Schichten (vgl. `collatz_eabc_epistemik_schichten.md`):

  | Schicht | Label            | Inhalt in diesem Modul                          |
  |---------|------------------|-------------------------------------------------|
  | **A**   | Theorem          | Kombinatorisch bewiesen (kein `sorry`)          |
  | **B**   | Struktur         | Definitionen, `Prop`-Skelette ohne Beweisanspruch |
  | **R**   | Forschungsbrücke | Asymptotik, Prim-Enumeration, EABC-Vermutung (`sorry`) |
  | **C**   | Ikone            | Nur Markdown — nicht formalisiert               |

  **Theorem ≠ Struktur ≠ Brücke ≠ Ikone**

  Referenz:
    collatz_eabc_diskrete_geometrie.md
    collatz_eabc_phi_e_lean_beweisversuch.md
    collatz_eabc_hodge_eabc.py
-/

import Mathlib.Data.Finset.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Order.Filter.AtTopBot.Defs
import Mathlib.Order.Filter.AtTopBot.Tendsto
import Mathlib.Topology.Basic
import CollatzEabc.HolonomieFehlerterm

namespace CollatzEabc

open Filter List Topology

/-!
## Schicht B — Struktur: G_E = C₄, orientierte Kanten E⁺ / E⁻

Vertices V = {E,A,B,C} = `EabcLetter` (E=0, A=1, B=2, C=3).
Forward cycle: E → A → B → C → E.
-/

/-- Gerichtete Kante auf dem Kreisgraphen G_E. -/
inductive C4DirectedEdge | EA | AB | BC | CE | EC | CB | BA | AE
  deriving DecidableEq, Repr, Fintype

/-- Vorwärtsorientierte Kantenmenge E⁺ = {EA, AB, BC, CE}. -/
def E_plus : Finset C4DirectedEdge :=
  {C4DirectedEdge.EA, C4DirectedEdge.AB, C4DirectedEdge.BC, C4DirectedEdge.CE}

/-- Rückwärtsorientierte Kantenmenge E⁻ = {EC, CB, BA, AE}. -/
def E_minus : Finset C4DirectedEdge :=
  {C4DirectedEdge.EC, C4DirectedEdge.CB, C4DirectedEdge.BA, C4DirectedEdge.AE}

/-- Quelle einer gerichteten Kante. -/
def edgeSrc : C4DirectedEdge → EabcLetter
  | .EA => ⟨0, by decide⟩
  | .AB => ⟨1, by decide⟩
  | .BC => ⟨2, by decide⟩
  | .CE => ⟨3, by decide⟩
  | .EC => ⟨3, by decide⟩
  | .CB => ⟨2, by decide⟩
  | .BA => ⟨1, by decide⟩
  | .AE => ⟨0, by decide⟩

/-- Ziel einer gerichteten Kante. -/
def edgeTgt : C4DirectedEdge → EabcLetter
  | .EA => ⟨1, by decide⟩
  | .AB => ⟨2, by decide⟩
  | .BC => ⟨3, by decide⟩
  | .CE => ⟨0, by decide⟩
  | .EC => ⟨0, by decide⟩
  | .CB => ⟨1, by decide⟩
  | .BA => ⟨2, by decide⟩
  | .AE => ⟨3, by decide⟩

/-- Elementarer Vorwärtszyklus als Knotenfolge E→A→B→C→E. -/
def forwardCycleVertices : List EabcLetter :=
  [⟨0, by decide⟩, ⟨1, by decide⟩, ⟨2, by decide⟩, ⟨3, by decide⟩, ⟨0, by decide⟩]

/-- Diskrete 1-Kozykel (Koeffizienten in ℤ). -/
abbrev C4Cochain := C4DirectedEdge → ℤ

/-- Diskreter Korand: (δf)(e) = f(tgt e) − f(src e) entlang einer Kante (Vorzeichenkonvention). -/
def coboundary (f : EabcLetter → ℤ) (e : C4DirectedEdge) : ℤ :=
  f (edgeTgt e) - f (edgeSrc e)

/-- Kanonischer harmonischer Generator: +1 auf E⁺, −1 auf E⁻. -/
def h_canonical : C4Cochain
  | .EA | .AB | .BC | .CE => 1
  | .EC | .CB | .BA | .AE => -1

/-- Kanonische harmonische 1-Form h auf C₄ (packagiert). -/
structure C4HarmonicForm where
  coeff : C4Cochain
  is_canonical : coeff = h_canonical

/-!
## Schicht A — Theorem: kombinatorische Kernidentitäten (bewiesen)
-/

/-- **[Schicht A]** E⁺ ⊔ E⁻ partitionieren die acht gerichteten Kanten. -/
theorem E_plus_union_E_minus :
    E_plus ∪ E_minus = Finset.univ ∧ Disjoint E_plus E_minus := by
  constructor <;> decide

/-- **[Schicht A]** Vorwärtszyklus E→A→B→C→E ist geschlossen. -/
theorem forward_cycle_closed :
    forwardCycleVertices.getLast? = forwardCycleVertices.head? := rfl

/-- **[Schicht A]** Kanonische harmonische 1-Form auf C₄ existiert (H¹(C₄) ≅ ℤ). -/
theorem harmonic_form_exists : Nonempty C4HarmonicForm :=
  ⟨⟨h_canonical, rfl⟩⟩

/-- **[Schicht A]** h ist kein Korand einer 0-Kette (nichttrivial in H¹). -/
theorem h_canonical_not_coboundary :
    ¬ ∃ f : EabcLetter → ℤ, ∀ e, h_canonical e = coboundary f e := by
  intro ⟨f, hf⟩
  have hA : f 1 = f 0 + 1 := by
    have := hf .EA; simp [h_canonical, coboundary, edgeSrc, edgeTgt] at this; linarith
  have hB : f 2 = f 1 + 1 := by
    have := hf .AB; simp [h_canonical, coboundary, edgeSrc, edgeTgt] at this; linarith
  have hC : f 3 = f 2 + 1 := by
    have := hf .BC; simp [h_canonical, coboundary, edgeSrc, edgeTgt] at this; linarith
  have hE : f 0 = f 3 + 1 := by
    have := hf .CE; simp [h_canonical, coboundary, edgeSrc, edgeTgt] at this; linarith
  linarith

/-!
## Schicht B — Struktur: C_E, S_E, W_E, Φ_E (Prim-Skeleton via HolonomieFehlerterm)
-/

/-- C_E(X) = N₊(X) − N₋(X) auf der κ-Folge bis Primgrenze X. -/
def C_E_up_to (X : ℕ) : ℤ :=
  D_E_up_to X

/-- S_E(X) = N₊(X) + N₋(X) = Q_E(X). -/
def S_E_up_to (X : ℕ) : ℕ :=
  Q_E_up_to X

-- W_E, R_beta, D_tilde_E bis X: siehe `HolonomieFehlerterm`.

/-- Φ_E existiert als Grenzwert von W_E(X) für X → ∞. -/
def HasPhi_E (φ : ℝ) : Prop :=
  Tendsto (fun X : ℕ => (W_E_up_to X : ℝ)) atTop (nhds φ)

/-- **[Schicht B]** Nullhypothese: asymptotische Symmetrie / Φ_E = 0. -/
def Phi_E_eq_zero : Prop :=
  HasPhi_E 0

/-- **[Schicht B]** α_E = 1: |D_E(X)| wächst asymptotisch linear in Q_E(X). -/
def HasAlpha_E_one : Prop :=
  ∃ C : ℝ, C ≠ 0 ∧
    Tendsto (fun X : ℕ => |((D_E_up_to X : ℤ) : ℝ)| / (Q_E_up_to X : ℝ)) atTop (nhds C)

/-- **[Schicht B]** EABC-Vermutung: Φ_E ≠ 0 — stabile arithmetische Orientierungsklasse. -/
def phi_E_conjecture : Prop :=
  ∃ φ : ℝ, HasPhi_E φ ∧ φ ≠ 0

/-!
## Schicht A — Theorem: endliche Symmetrie ⇒ W_E = 0
-/

/-- **[Schicht A]** Auf endlicher Folge: N₊ = N₋ impliziert χ_Hol = 0. -/
theorem W_E_list_zero_of_balance (classes : List EabcLetter)
    (h : N_plus classes = N_minus classes) : chi_Hol classes = 0 :=
  chi_Hol_zero_of_balance classes h

/-!
## Schicht B — Struktur: diskrete 1-Form ω_E und Paarung ⟨ω_E, h⟩
-/

/-- Gewichte einer diskreten 1-Form ω_E auf Kanten (endliche Stichprobe). -/
abbrev OmegaE := C4DirectedEdge → ℤ

/-- ⟨ω_E, h⟩ = Σ_{e∈E⁺} ω(e) − Σ_{e∈E⁻} ω(e) für kanonisches h. -/
def innerProductOmegaH (ω : OmegaE) : ℤ :=
  (E_plus.sum fun e => ω e * h_canonical e) +
    (E_minus.sum fun e => ω e * h_canonical e)

/-- Zirkulation C_E aus Kantengewichten (diskret). -/
def circulationOmega (ω : OmegaE) : ℤ :=
  (E_plus.sum ω) - (E_minus.sum ω)

lemma h_coeff_E_plus (e : C4DirectedEdge) (he : e ∈ E_plus) : h_canonical e = 1 := by
  rcases e <;> simp [E_plus, h_canonical] at he ⊢

lemma h_coeff_E_minus (e : C4DirectedEdge) (he : e ∈ E_minus) : h_canonical e = -1 := by
  rcases e <;> simp [E_minus, h_canonical] at he ⊢

/-!
## Schicht A — Theorem: diskrete Paarung ⟨ω_E, h⟩ = C_E
-/

/-- **[Schicht A]** Diskrete Paarung ⟨ω,h⟩ = Zirkulation C_E bei ω auf E⁺∪E⁻. -/
theorem Phi_E_eq_inner_product_discrete (ω : OmegaE) :
    innerProductOmegaH ω = circulationOmega ω := by
  unfold innerProductOmegaH circulationOmega
  have hsum : E_plus.sum (fun e => ω e * h_canonical e) = E_plus.sum ω := by
    refine Finset.sum_congr rfl fun e he => ?_
    simp [h_coeff_E_plus e he, mul_one]
  have hsum' : E_minus.sum (fun e => ω e * h_canonical e) = -E_minus.sum ω := by
    rw [← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun e he => ?_
    simp [h_coeff_E_minus e he]
  rw [hsum, hsum']
  simp [sub_eq_add_neg]

/-!
## Schicht R — Forschungsbrücken (`sorry`): Asymptotik, Prim-Enumeration, EABC-Vermutung

Diese Sektion verbindet Schicht-A-Endlichkeit mit Schicht-B-Grenzwerten.
Prim-Zählung (`N_plus_up_to`, `N_minus_up_to`, `Hol_E_zero`) in `HolonomieFehlerterm`.
-/

/-- **[Schicht R]** Asymptotische Symmetrie N₊∼N₋ ⇒ Φ_E=0 (Primteil `sorry`). -/
theorem Phi_E_zero_of_symmetry
    (h : ∀ᶠ X in atTop, N_plus_up_to X = N_minus_up_to X) : Phi_E_eq_zero := by
  unfold Phi_E_eq_zero HasPhi_E
  sorry

/-- **[Schicht R]** Brücke Folge → Prim-Obergrenze: asymptotische Paarung Φ_E ↔ ⟨ω_E,h⟩ (offen). -/
theorem Phi_E_eq_inner_product :
    ∀ φ : ℝ, HasPhi_E φ →
      ∃ ω : ℕ → OmegaE,
        sorry := by
  intro φ hφ
  sorry

/-- **[Schicht R]** HL-Symmetrie (N₊ ≈ N₋) folgt Hol_E = Φ_E = 0. -/
theorem hol_E_zero_of_HL (hHL : Hol_E_zero) : Phi_E_eq_zero := by
  sorry

/-- **[Schicht R]** Φ_E ≠ 0 ⇒ |D_E| ∼ c·Q asymptotisch ⇒ α_E = 1 (keine Äquivalenz). -/
theorem phi_E_ne_zero_implies_alpha_E_one
    {φ : ℝ} (hφ : HasPhi_E φ) (hne : φ ≠ 0) : HasAlpha_E_one := by
  sorry

/-- **[Schicht R]** EABC-Vermutung Φ_E ≠ 0 — stabile Orientierungsklasse (unbewiesen). -/
theorem phi_E_conjecture_statement : phi_E_conjecture := by
  sorry

end CollatzEabc
