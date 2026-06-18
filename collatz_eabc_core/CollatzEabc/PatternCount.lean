/-
  CollatzEabc.PatternCount — Kernmodul Ebene A+B: diskrete Zählgrößen D_E, Q_E, N_±.

  Hierarchie:

    Residue → windows5, isABCEA, isCEABC
           → N_plus_up_to / N_minus_up_to  (`kappaPrimeStreamUpTo`, Black Box)
           → D_E / Q_E  (**Kern**)
           → W_E
           → Φ_E  (Conjecture in `FlussPhiE` / `HolonomyCore`)

  D_E, Q_E vor W_E, α_E, Φ_E.

  **Primärträger:** Gleitfenster `windows5` auf `kappaPrimeStreamUpTo X`.
  **Sekundärträger:** Primvierlinge — `N_plus_quadruplet_up_to` / `N_minus_quadruplet_up_to`.
-/

import Mathlib.Data.Nat.Prime.Defs
import Mathlib.Data.List.Basic
import Mathlib.Data.Int.Basic
import Mathlib.Data.Rat.Defs
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import CollatzEabc.Kappa

namespace CollatzEabc

open List

/-!
### Ebene A — reine Kombinatorik (KEINE Primzahltheorie)

`Residue` ≅ `EabcLetter` (Fin 4: E=0, A=1, B=2, C=3).

| Zyklus | Wort  | Orientierung      |
|--------|-------|-------------------|
| γ⁺     | ABCEA | A → B → C → E → A |
| γ⁻     | CEABC | C → E → A → B → C |
-/

abbrev Residue := EabcLetter

def wordABCEA : List Residue := [⟨1, by decide⟩, ⟨2, by decide⟩, ⟨3, by decide⟩,
  ⟨0, by decide⟩, ⟨1, by decide⟩]

def wordCEABC : List Residue := [⟨3, by decide⟩, ⟨0, by decide⟩, ⟨1, by decide⟩,
  ⟨2, by decide⟩, ⟨3, by decide⟩]

def isABCEA (xs : List Residue) : Bool :=
  decide (xs = wordABCEA)

def isCEABC (xs : List Residue) : Bool :=
  decide (xs = wordCEABC)

/-- Primäre Zähl-API: alle 5er-Gleitfenster. -/
def windows5 (stream : List Residue) : List (List Residue) :=
  let k := 5
  if _h : stream.length < k then []
  else
    (List.range (stream.length + 1 - k)).map fun i =>
      List.take k (stream.drop i)

/-!
### Ebene B — Zählgrößen; `kappaPrimeStreamUpTo` = Black Box (computabel)
-/

def isPrime (n : ℕ) : Bool :=
  decide (Nat.Prime n)

def kappaPrimeStreamUpTo (X : ℕ) : List Residue :=
  (List.range (X + 1)).filterMap fun p =>
    if _h : 3 < p ∧ Nat.Prime p then classOfLetter p else none

/-- **[Schicht R]** mod-12-Korrektheit der κ-Primfolge (offen). -/
def kappaPrimeStream_mod12_correct (_X : ℕ) : Prop :=
  sorry

abbrev primeEabcClassesUpTo (X : ℕ) : List EabcLetter :=
  kappaPrimeStreamUpTo X

def primeEabcCountUpTo (X : ℕ) : ℕ :=
  (kappaPrimeStreamUpTo X).length

def N_plus_up_to (X : ℕ) : ℕ :=
  ((windows5 (kappaPrimeStreamUpTo X)).filter isABCEA).length

def N_minus_up_to (X : ℕ) : ℕ :=
  ((windows5 (kappaPrimeStreamUpTo X)).filter isCEABC).length

def D_E (X : ℕ) : ℤ :=
  (N_plus_up_to X : ℤ) - N_minus_up_to X

def Q_E (X : ℕ) : ℕ :=
  N_plus_up_to X + N_minus_up_to X

def W_E (X : ℕ) : ℚ :=
  if h : Q_E X = 0 then 0
  else (D_E X : ℚ) / Q_E X

def W_E_of_pos (X : ℕ) (_h : 0 < Q_E X) : ℚ :=
  (D_E X : ℚ) / Q_E X

noncomputable def R_half (X : ℕ) : ℝ :=
  let q := Q_E X
  if _h : q = 0 then 0
  else (D_E X : ℝ) / Real.sqrt q

abbrev D_E_up_to := D_E
abbrev Q_E_up_to := Q_E
abbrev W_E_up_to := W_E
noncomputable abbrev D_tilde_E_up_to := R_half

/-!
### Schicht A — erste formale Sätze (GREEN)
-/

theorem Q_E_eq (X : ℕ) : Q_E X = N_plus_up_to X + N_minus_up_to X := rfl

abbrev Q_E_eq_sum := Q_E_eq

theorem D_E_eq_diff (X : ℕ) : D_E X = (N_plus_up_to X : ℤ) - N_minus_up_to X := rfl

theorem abs_D_E_le_Q_E (X : ℕ) : |D_E X| ≤ Q_E X := by
  unfold D_E Q_E
  have h₁ :
      -(N_plus_up_to X + N_minus_up_to X : ℤ) ≤
        (N_plus_up_to X : ℤ) - N_minus_up_to X := by omega
  have h₂ :
      (N_plus_up_to X : ℤ) - N_minus_up_to X ≤
        (N_plus_up_to X + N_minus_up_to X : ℤ) := by omega
  exact abs_le.mpr ⟨h₁, h₂⟩

abbrev D_E_abs_le_Q_E := abs_D_E_le_Q_E

theorem W_E_of_pos_eq (X : ℕ) (h : 0 < Q_E X) : W_E_of_pos X h = W_E X := by
  unfold W_E_of_pos W_E Q_E D_E
  simp [ne_of_gt h]

theorem W_E_bounds (X : ℕ) : -1 ≤ W_E X ∧ W_E X ≤ 1 := by
  unfold W_E Q_E D_E
  by_cases h : Q_E X = 0
  · rw [dif_pos h]
    norm_num
  · have hpos_nat : 0 < Q_E X := Nat.pos_of_ne_zero h
    have hpos_rat : (0 : ℚ) < (Q_E X : ℚ) := by exact_mod_cast hpos_nat
    rw [dif_neg h]
    constructor
    · have hineq_int :
          -(Q_E X : Int) ≤ (N_plus_up_to X : Int) - (N_minus_up_to X : Int) := by
        unfold Q_E; omega
      have hineq_rat :
          (-(Q_E X : Int) : ℚ) ≤ ((N_plus_up_to X : Int) - (N_minus_up_to X : Int) : ℚ) := by
        exact_mod_cast hineq_int
      have hdiv := div_le_div_of_nonneg_right hineq_rat (le_of_lt hpos_rat)
      have hneg : (-(Q_E X : Int) : ℚ) / Q_E X = -1 := by
        rw [show (-(Q_E X : Int) : ℚ) = -(Q_E X : ℚ) from by push_cast; rfl]
        field_simp [ne_of_gt hpos_rat]
      rw [← hneg]
      exact hdiv
    · have hineq_int :
          (N_plus_up_to X : Int) - (N_minus_up_to X : Int) ≤ (Q_E X : Int) := by
        unfold Q_E; omega
      have hineq_rat :
          (((N_plus_up_to X : Int) - (N_minus_up_to X : Int) : Int) : ℚ) ≤ (Q_E X : ℚ) := by
        exact_mod_cast hineq_int
      have hdiv := div_le_div_of_nonneg_right hineq_rat (le_of_lt hpos_rat)
      have hone : ((Q_E X : Int) : ℚ) / Q_E X = 1 := by
        rw [show ((Q_E X : Int) : ℚ) = (Q_E X : ℚ) from by push_cast; rfl]
        field_simp [ne_of_gt hpos_rat]
      rw [← hone]
      exact hdiv

theorem W_E_bounds_of_pos (X : ℕ) (h : 0 < Q_E X) : |W_E_of_pos X h| ≤ 1 := by
  rw [W_E_of_pos_eq X h]
  rcases W_E_bounds X with ⟨hlo, hhi⟩
  exact abs_le.mpr ⟨hlo, hhi⟩

abbrev W_E_abs_le_one_of_Q_pos := W_E_bounds_of_pos

theorem W_E_zero_of_balance (X : ℕ) (h : N_plus_up_to X = N_minus_up_to X) : W_E X = 0 := by
  unfold W_E Q_E D_E
  simp [h]

/-!
### Sekundärträger: Primvierlinge (vgl. `eabc_quadruplets_1e10.py`)
-/

def isPrimeQuadruplet (p : ℕ) : Bool :=
  isPrime p && isPrime (p + 2) && isPrime (p + 6) && isPrime (p + 8)

def isQuadrupletPlus (p : ℕ) : Bool :=
  isPrimeQuadruplet p && p % 12 = 5

def isQuadrupletMinus (p : ℕ) : Bool :=
  isPrimeQuadruplet p && p % 12 = 11

def N_plus_quadruplet_up_to (X : ℕ) : ℕ :=
  ((List.range (X + 1)).filter isQuadrupletPlus).length

def N_minus_quadruplet_up_to (X : ℕ) : ℕ :=
  ((List.range (X + 1)).filter isQuadrupletMinus).length

def D_E_quadruplet_up_to (X : ℕ) : ℤ :=
  (N_plus_quadruplet_up_to X : ℤ) - N_minus_quadruplet_up_to X

def Q_E_quadruplet_up_to (X : ℕ) : ℕ :=
  N_plus_quadruplet_up_to X + N_minus_quadruplet_up_to X

example : N_plus_up_to 1000 = 4 := by native_decide

example : N_minus_up_to 1000 = 4 := by native_decide

example : D_E 1000 = 0 := by native_decide

example : N_plus_quadruplet_up_to 1000 = 3 := by native_decide

example : N_minus_quadruplet_up_to 1000 = 2 := by native_decide

/-!
Zeuge §4.4 (`collatz_eabc_zirkulationshypothese.md`):
  p = 54044321, Primvierling (54044321, 54044323, 54044327, 54044329),
  p % 12 = 5 → ABCE (`isQuadrupletPlus`), p % 60060 = 50381.
  Verifikation: `eabc_witness_54044321_verify.py`.
-/

-- Gleitfenster (primär), X = 10⁶
#eval N_plus_up_to 1000000

#eval N_minus_up_to 1000000

#eval D_E 1000000

-- Primvierlinge (sekundär, Python 84/82), X = 10⁶
#eval N_plus_quadruplet_up_to 1000000

#eval N_minus_quadruplet_up_to 1000000

end CollatzEabc
