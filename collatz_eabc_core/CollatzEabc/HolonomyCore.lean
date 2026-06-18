/-
EABC Holonomy Core

Status:
- Kreisgraph / Orientierung / Zirkulation: formal scaffold (abstrakter `CycleCounts`-Kern)
- bewiesen: `-1 ≤ W_E(X) ≤ 1` auf endlichen Stichproben (`phiApprox_bounds`, `W_E_bounds`)
- Primgrenze X: dieselbe Schranke in `CollatzEabc.PatternCount.W_E_bounds`
- Vermutung: `lim_{X→∞} W_E(X) = Φ_E ≠ 0` (asymptotisch, nicht exakte Konstante)
- keine physikalische Aussage

Hierarchie: D_E, Q_E (`PatternCount`) → W_E → Φ_E (hier: Schicht R).

GREEN LAYER:
  Node, next, prev, CycleCounts, circulation, size, phiApprox,
  phiApprox_bounds, W_E_bounds.
RED LAYER:
  HasNonzeroHolonomyLimit and EABC_holonomy_limit_conjecture.
  This is an explicit conjectural bridge, not a theorem.
-/

import Mathlib
import CollatzEabc.PatternCount
import Mathlib.Order.Filter.AtTopBot.Tendsto
import Mathlib.Topology.Basic

namespace EABC

inductive Node where
  | E | A | B | C
deriving DecidableEq, Repr

open Node
open Filter Topology

def next : Node → Node
  | E => A
  | A => B
  | B => C
  | C => E

def prev : Node → Node
  | E => C
  | A => E
  | B => A
  | C => B

structure CycleCounts where
  Nplus : Nat
  Nminus : Nat

def circulation (c : CycleCounts) : Int :=
  (c.Nplus : Int) - (c.Nminus : Int)

def size (c : CycleCounts) : Nat :=
  c.Nplus + c.Nminus

def phiApprox (c : CycleCounts) : Rat :=
  if _h : size c = 0 then 0
  else ((circulation c : Int) : Rat) / ((size c : Nat) : Rat)

theorem phiApprox_bounds (c : CycleCounts) :
    -1 ≤ phiApprox c ∧ phiApprox c ≤ 1 := by
  unfold phiApprox circulation size
  by_cases h : c.Nplus + c.Nminus = 0
  · rw [dif_pos h]
    norm_num
  ·
    have hpos_nat : 0 < c.Nplus + c.Nminus := Nat.pos_of_ne_zero h
    have hpos_rat : (0 : Rat) < ((c.Nplus + c.Nminus : Nat) : Rat) := by
      exact_mod_cast hpos_nat
    rw [dif_neg h]
    constructor
    ·
      have hineq_int :
          -((c.Nplus + c.Nminus : Nat) : Int)
            ≤ (c.Nplus : Int) - (c.Nminus : Int) := by
        omega
      have hineq_rat :
          (-(((c.Nplus + c.Nminus : Nat) : Int) : Rat))
            ≤ (((c.Nplus : Int) - (c.Nminus : Int) : Int) : Rat) := by
        exact_mod_cast hineq_int
      have hdiv := div_le_div_of_nonneg_right hineq_rat (le_of_lt hpos_rat)
      have hneg :
          (-(((c.Nplus + c.Nminus : Nat) : Int) : Rat)) / ((c.Nplus + c.Nminus : Nat) : Rat) = -1 := by
        rw [show (-(((c.Nplus + c.Nminus : Nat) : Int) : Rat))
              = -((c.Nplus + c.Nminus : Nat) : Rat) from by push_cast; rfl]
        field_simp [ne_of_gt hpos_rat]
      rw [← hneg]
      exact hdiv
    ·
      have hineq_int :
          (c.Nplus : Int) - (c.Nminus : Int)
            ≤ ((c.Nplus + c.Nminus : Nat) : Int) := by
        omega
      have hineq_rat :
          (((c.Nplus : Int) - (c.Nminus : Int) : Int) : Rat)
            ≤ (((c.Nplus + c.Nminus : Nat) : Int) : Rat) := by
        exact_mod_cast hineq_int
      have hdiv := div_le_div_of_nonneg_right hineq_rat (le_of_lt hpos_rat)
      have hone :
          (((c.Nplus + c.Nminus : Nat) : Int) : Rat) / ((c.Nplus + c.Nminus : Nat) : Rat) = 1 := by
        rw [show (((c.Nplus + c.Nminus : Nat) : Int) : Rat)
              = ((c.Nplus + c.Nminus : Nat) : Rat) from by push_cast; rfl]
        field_simp [ne_of_gt hpos_rat]
      rw [← hone]
      exact hdiv

structure EABCFlow where
  counts : Nat → CycleCounts

def C_E (F : EABCFlow) (X : Nat) : Int :=
  circulation (F.counts X)

def S_E (F : EABCFlow) (X : Nat) : Nat :=
  size (F.counts X)

def W_E (F : EABCFlow) (X : Nat) : Rat :=
  phiApprox (F.counts X)

theorem W_E_bounds (F : EABCFlow) (X : Nat) :
    -1 ≤ W_E F X ∧ W_E F X ≤ 1 := by
  unfold W_E
  exact phiApprox_bounds (F.counts X)

/-- Primgrenze: `PatternCount.W_E_bounds` (Gleitfenster-Träger auf κ-Primfolge). -/
abbrev W_E_prime_bounds := CollatzEabc.W_E_bounds

/-- **Schicht R:** `W_E(X)` konvergiert gegen ein von Null verschiedenes `Φ ∈ ℝ`. -/
def HasNonzeroHolonomyLimit (F : EABCFlow) : Prop :=
  ∃ Φ : ℝ, Φ ≠ 0 ∧
    Tendsto (fun X => ((W_E F X : Rat) : ℝ)) atTop (𝓝 Φ)

theorem EABC_holonomy_limit_conjecture (F : EABCFlow) :
    HasNonzeroHolonomyLimit F := by
  sorry

-- Too-strong scaffold (exact rational constant eventually):
-- def HasNonzeroHolonomy (F : EABCFlow) : Prop :=
--   ∃ Φ : Rat, Φ ≠ 0 ∧ ∀ᶠ X in atTop, W_E F X = Φ
-- theorem EABC_holonomy_conjecture (F : EABCFlow) :
--     HasNonzeroHolonomy F := by sorry

end EABC
