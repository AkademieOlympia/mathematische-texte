/-
  CollatzEabc.Mod12Matrix — mod-12 EABC-Übergangsmatrix und Irreduzibilität.

  EABC-Klassen: E ≡ 1, A ≡ 5, B ≡ 7, C ≡ 11 (mod 12).
  Referenz: collatz_mixing_test.py, collatz_schlussartikel_arxiv.tex (§ RPF).
-/

import Mathlib.Data.ZMod.Basic
import Mathlib.LinearAlgebra.Matrix.Irreducible.Defs
import Mathlib.Data.Rat.Defs
import Mathlib.Tactic

namespace CollatzEabc

abbrev EabcIndex := Fin 4

def eabcResidue : EabcIndex → ZMod 12
  | 0 => 1
  | 1 => 5
  | 2 => 7
  | 3 => 11

private def eabcRowEntry : EabcIndex → ℚ
  | 0 => (3 : ℚ) / 8
  | 1 => (1 : ℚ) / 8
  | 2 => (3 : ℚ) / 8
  | 3 => (1 : ℚ) / 8

/-- Vereinfachte 4×4-Übergangsmatrix (uniforme Zeilen, Schlussartikel). -/
def eabcMod12TransitionSimplified : Matrix EabcIndex EabcIndex ℚ :=
  Matrix.of fun _ j => eabcRowEntry j

private lemma eabcRowEntry_nonneg (j : EabcIndex) : 0 ≤ eabcRowEntry j := by
  fin_cases j <;> norm_num [eabcRowEntry]

private lemma eabcRowEntry_pos (j : EabcIndex) : 0 < eabcRowEntry j := by
  fin_cases j <;> norm_num [eabcRowEntry]

theorem eabcMod12TransitionSimplified_nonneg (i j : EabcIndex) :
    0 ≤ eabcMod12TransitionSimplified i j := by
  simp [eabcMod12TransitionSimplified, eabcRowEntry_nonneg]

theorem eabcMod12TransitionSimplified_pos (i j : EabcIndex) :
    0 < eabcMod12TransitionSimplified i j := by
  simp [eabcMod12TransitionSimplified, eabcRowEntry_pos]

/-- Die vereinfachte mod-12-Markov-Matrix ist irreduzibel (Kante i → j für alle Paare). -/
theorem eabc_mod12_matrix_irreducible :
    Matrix.IsIrreducible eabcMod12TransitionSimplified where
  nonneg := eabcMod12TransitionSimplified_nonneg
  connected i j := by
    letI : Quiver EabcIndex := Matrix.toQuiver eabcMod12TransitionSimplified
    refine ⟨Quiver.Path.nil.cons ⟨eabcMod12TransitionSimplified_pos i j⟩, ?_⟩
    simp

private def empiricalRowEntry : EabcIndex → EabcIndex → ℚ
  | 0, 0 => (1 : ℚ) / 3
  | 0, 1 => (1 : ℚ) / 6
  | 0, 2 => (1 : ℚ) / 3
  | 0, 3 => (1 : ℚ) / 6
  | 1, 0 => (1 : ℚ) / 3
  | 1, 1 => (1 : ℚ) / 6
  | 1, 2 => (1 : ℚ) / 3
  | 1, 3 => (1 : ℚ) / 6
  | 2, 0 => 0
  | 2, 1 => (1 : ℚ) / 2
  | 2, 2 => 0
  | 2, 3 => (1 : ℚ) / 2
  | 3, 0 => 0
  | 3, 1 => (1 : ℚ) / 2
  | 3, 2 => 0
  | 3, 3 => (1 : ℚ) / 2

/-- Empirische Matrix aus collatz_mixing_test.py (n ≤ 5·10⁵). -/
def eabcMod12TransitionEmpirical : Matrix EabcIndex EabcIndex ℚ :=
  Matrix.of empiricalRowEntry

theorem eabcMod12TransitionEmpirical_nonneg (i j : EabcIndex) :
    0 ≤ eabcMod12TransitionEmpirical i j := by
  fin_cases i <;> fin_cases j <;> simp [eabcMod12TransitionEmpirical, empiricalRowEntry] <;> norm_num

end CollatzEabc
