/-
  CollatzEabc.BernoulliClock — Bernoulli-Tripel und Radialgewichte (Gedankenmodell).

  Kein Collatz-Anspruch: nur Indexarithmetik und Mathlib-`bernoulli`.
  Referenz: collatz_kepler_gedankenexperiment.tex, Abschnitt Bernoulli-Uhr.
-/

import Mathlib.NumberTheory.Bernoulli
import CollatzEabc.PrefProjection

namespace CollatzEabc

/-- Index einer Bernoulli-Uhr-Zelle (m ≥ 1). -/
structure BernoulliCell where
  m : ℕ
  hm : 1 ≤ m

/-- Radialgewicht r_m = 2^{-m}. -/
def radialAt (m : ℕ) : ℚ :=
  (1 : ℚ) / (2 ^ m)

/-- Phasenindex k_m = m mod 4 (Dreiertakt auf dem EABC-Ring). -/
def dreierPhase (m : ℕ) : Fin 4 :=
  ⟨m % 4, Nat.mod_lt _ (by decide)⟩

/-- Indizes (2m-2, 2m, 2m+2) des Bernoulli-Triplets. -/
def bernoulliTripletIndices (c : BernoulliCell) : ℕ × ℕ × ℕ :=
  let m := c.m
  (2 * m - 2, 2 * m, 2 * m + 2)

/-- Bernoulli-Tripel (B_{2m-2}, B_{2m}, B_{2m+2}). -/
def bernoulliTriplet (c : BernoulliCell) : ℚ × ℚ × ℚ :=
  let (i, j, k) := bernoulliTripletIndices c
  (bernoulli i, bernoulli j, bernoulli k)

/-- Urfall m = 1: (B_0, B_2, B_4) = (1, 1/6, -1/30). -/
theorem bernoulliTriplet_one :
    bernoulliTriplet ⟨1, by decide⟩ = (1, 6⁻¹, -(30 : ℚ)⁻¹) := by
  unfold bernoulliTriplet bernoulliTripletIndices
  simp only [bernoulli_zero, bernoulli_two]
  have h4 : bernoulli 4 = -(30 : ℚ)⁻¹ := by
    rw [bernoulli, bernoulli'_four]
    norm_num
  simp [h4]

/-- r_m stimmt mit ρ(w) für |w| = m überein. -/
theorem radialAt_eq_rho_length (m : ℕ) (w : List EabcLetter) (hw : w.length = m) :
    radialAt m = rho w := by
  simp [radialAt, rho, hw]

end CollatzEabc
