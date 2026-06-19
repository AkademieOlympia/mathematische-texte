/-
  CollatzEabc.V4Characters — V₄ ≅ C₂×C₂ Charaktere (§4.8.23).

  Schicht **A**: algebraisch feste Charakterwerte auf {E,A,B,C}.
  Kein Quaternion-/Oktanion-Bild; keine physikalische Ladungsontologie.

  Querverweise: `Mod12Matrix.lean`, `collatz_eabc_zirkulationshypothese.md` §4.8.23.
-/

import Mathlib

namespace CollatzEabc

/-- EABC-Labels als Klein-Vierer-Gruppe V₄. -/
inductive EabcLabel | E | A | B | C
  deriving DecidableEq, Repr

open EabcLabel

/-- Nichttrivialer Charakter χ_A: E↦1, A↦-1, B↦1, C↦-1. -/
def chi_A : EabcLabel → Int
  | E => 1
  | A => -1
  | B => 1
  | C => -1

/-- Nichttrivialer Charakter χ_B: E↦1, A↦1, B↦-1, C↦-1. -/
def chi_B : EabcLabel → Int
  | E => 1
  | A => 1
  | B => -1
  | C => -1

/-- Nichttrivialer Charakter χ_C: E↦1, A↦-1, B↦-1, C↦1. -/
def chi_C : EabcLabel → Int
  | E => 1
  | A => -1
  | B => -1
  | C => 1

/-- Trivialer Charakter (Identität). -/
def chi_trivial : EabcLabel → Int
  | _ => 1

private def mulLabel : EabcLabel → EabcLabel → EabcLabel
  | E, x => x
  | A, E => A
  | A, A => E
  | A, B => C
  | A, C => B
  | B, E => B
  | B, A => C
  | B, B => E
  | B, C => A
  | C, E => C
  | C, A => B
  | C, B => A
  | C, C => E

theorem mulLabel_assoc (a b c : EabcLabel) :
    mulLabel (mulLabel a b) c = mulLabel a (mulLabel b c) := by
  cases a <;> cases b <;> cases c <;> rfl

theorem mulLabel_one_left (a : EabcLabel) : mulLabel E a = a := by
  cases a <;> rfl

theorem mulLabel_one_right (a : EabcLabel) : mulLabel a E = a := by
  cases a <;> rfl

theorem chi_A_mul (x y : EabcLabel) :
    chi_A (mulLabel x y) = chi_A x * chi_A y := by
  cases x <;> cases y <;> simp [chi_A, mulLabel]

theorem chi_B_mul (x y : EabcLabel) :
    chi_B (mulLabel x y) = chi_B x * chi_B y := by
  cases x <;> cases y <;> simp [chi_B, mulLabel]

theorem chi_C_mul (x y : EabcLabel) :
    chi_C (mulLabel x y) = chi_C x * chi_C y := by
  cases x <;> cases y <;> simp [chi_C, mulLabel]

theorem chi_values_sq (χ : EabcLabel → Int) (h : ∀ x, χ x * χ x = 1) (x : EabcLabel) :
    χ x * χ x = 1 := h x

theorem chi_A_sq (x : EabcLabel) : chi_A x * chi_A x = 1 := by
  cases x <;> norm_num [chi_A]

theorem chi_B_sq (x : EabcLabel) : chi_B x * chi_B x = 1 := by
  cases x <;> norm_num [chi_B]

theorem chi_C_sq (x : EabcLabel) : chi_C x * chi_C x = 1 := by
  cases x <;> norm_num [chi_C]

end CollatzEabc
