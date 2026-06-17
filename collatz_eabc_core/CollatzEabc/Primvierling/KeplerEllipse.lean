/-
  CollatzEabc.Primvierling.KeplerEllipse — Kepler-Ellipse als Ebene-B′/C-Repräsentation.

  Die Ellipse kodiert die in `Mod12` bewiesene mod-12-Chiralität als geometrische Uhr;
  sie erzeugt keine Primzahlen und beweist keine Existenz von Primvierlingen.

  Kette: Lean-Bijektion → chirale EABC-Phasen → Kepler-Ellipse → diskrete Uhr (t ∈ {0,1,2,3}).
-/

import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Data.Complex.Basic
import Mathlib.Data.Real.Sqrt
import Mathlib.Tactic
import CollatzEabc.Primvierling.Basic
import CollatzEabc.Primvierling.Chirality
import CollatzEabc.Primvierling.Mod12

namespace CollatzEabc
namespace Primvierling

/-- Halbachse der Primvierlingsellipse (Normalform). -/
def a_PV : ℝ := 4

/-- Nebenhalbachse der Primvierlingsellipse (Normalform). -/
def b_PV : ℝ := 2

/-- Exzentrizität e_PV = √(1 − b²/a²) = √3/2 (≠ ρ_PV = 3/2). -/
noncomputable def e_PV : ℝ :=
  Real.sqrt (1 - b_PV ^ 2 / a_PV ^ 2)

/-- Phasenwinkel der EABC-Familien auf der Kepler-Ellipse. -/
noncomputable def flavorAngle : EClass → ℝ
  | EClass.E => 0
  | EClass.A => Real.pi / 2
  | EClass.B => Real.pi
  | EClass.C => 3 * Real.pi / 2

/-- Diskreter Phasentick: θ(t) = (π/2)·t für t ∈ {0,1,2,3}. -/
noncomputable def phaseTick (t : Fin 4) : ℝ :=
  (Real.pi / 2) * t

/-- EABC-Flavor zum Phasentick t (geometrische Uhr E → A → B → C). -/
def tickFlavor (t : Fin 4) : EClass :=
  match t with
  | ⟨0, _⟩ => EClass.E
  | ⟨1, _⟩ => EClass.A
  | ⟨2, _⟩ => EClass.B
  | ⟨3, _⟩ => EClass.C

/-- Chiralitätsphasenverschiebung: CEAB = ABCE um π. -/
noncomputable def chiralityPhaseShift : ℝ :=
  Real.pi

/-- Reelle Koordinaten (Re z, Im z) der Kepler-Ellipse. -/
noncomputable def ellipsePointR (M θ : ℝ) : ℝ × ℝ :=
  (M + a_PV * Real.cos θ, b_PV * Real.sin θ)

/-- Kepler-Ellipse E_PV(θ) = M + a·cos θ + i·b·sin θ in ℂ. -/
noncomputable def ellipsePoint (M θ : ℝ) : ℂ :=
  Complex.ofReal (ellipsePointR M θ).1 + Complex.I * Complex.ofReal (ellipsePointR M θ).2

theorem ellipsePointR_eq (M θ : ℝ) :
    ellipsePointR M θ = ⟨(ellipsePoint M θ).re, (ellipsePoint M θ).im⟩ := by
  dsimp [ellipsePointR, ellipsePoint]
  simp [Complex.add_re, Complex.add_im, Complex.mul_re, Complex.mul_im,
    Complex.I_re, Complex.I_im, Complex.ofReal_re, Complex.ofReal_im]

/-- Chirale Parametrisierung: ABCE = E₊(θ), CEAB = E₊(θ + π). -/
noncomputable def ellipsePointChiral (M : ℝ) (ch : Chirality) (θ : ℝ) : ℂ :=
  match ch with
  | Chirality.ABCE => ellipsePoint M θ
  | Chirality.CEAB => ellipsePoint M (θ + chiralityPhaseShift)

/-- Markierter Primvierling: Ellipse bei Chiralität aus `chiralityWord`. -/
noncomputable def ellipsePointMarked (p : ℕ) (θ : ℝ) : Option ℂ :=
  chiralityWord p |>.map fun ch => ellipsePointChiral (center p : ℝ) ch θ

theorem a_PV_eq_four : a_PV = 4 := rfl

theorem b_PV_eq_two : b_PV = 2 := rfl

theorem e_PV_eq_sqrt_three_halves :
    e_PV = Real.sqrt 3 / 2 := by
  unfold e_PV a_PV b_PV
  norm_num

theorem e_PV_lt_one : e_PV < 1 := by
  rw [e_PV_eq_sqrt_three_halves]
  have h1 : Real.sqrt 3 < Real.sqrt 4 :=
    Real.sqrt_lt_sqrt (by norm_num : (0 : ℝ) ≤ 3) (by norm_num : (3 : ℝ) < 4)
  have h2 : Real.sqrt 4 = 2 := by norm_num
  linarith

theorem flavorAngle_E : flavorAngle EClass.E = 0 := rfl

theorem flavorAngle_A : flavorAngle EClass.A = Real.pi / 2 := rfl

theorem flavorAngle_B : flavorAngle EClass.B = Real.pi := rfl

theorem flavorAngle_C : flavorAngle EClass.C = 3 * Real.pi / 2 := rfl

theorem chiralityPhaseShift_eq_pi : chiralityPhaseShift = Real.pi := rfl

theorem phaseTick_zero : phaseTick 0 = 0 := by
  simp [phaseTick]

theorem phaseTick_one : phaseTick 1 = Real.pi / 2 := by
  simp [phaseTick]

theorem phaseTick_two : phaseTick 2 = Real.pi := by
  simp [phaseTick]

theorem phaseTick_three : phaseTick 3 = 3 * Real.pi / 2 := by
  simp [phaseTick]
  ring

theorem flavorAngle_eq_phaseTick (t : Fin 4) :
    flavorAngle (tickFlavor t) = phaseTick t := by
  fin_cases t <;> simp [flavorAngle, phaseTick, tickFlavor] <;> ring

theorem phaseTick_flavor_clock (t : Fin 4) :
    flavorAngle (tickFlavor t) = (Real.pi / 2) * t :=
  flavorAngle_eq_phaseTick t

theorem ellipsePointChiral_CEAB_eq (M θ : ℝ) :
    ellipsePointChiral M Chirality.CEAB θ =
      ellipsePointChiral M Chirality.ABCE (θ + chiralityPhaseShift) := by
  simp [ellipsePointChiral, chiralityPhaseShift_eq_pi]

private theorem ellipsePoint_neg_shift (M θ : ℝ) :
    ellipsePoint M (θ + Real.pi) =
      Complex.ofReal (2 * M) - ellipsePoint M θ := by
  dsimp [ellipsePoint, ellipsePointR]
  rw [Real.cos_add_pi, Real.sin_add_pi]
  simp [Complex.ofReal_add, Complex.ofReal_sub, Complex.ofReal_mul, Complex.ofReal_neg]
  ring

/-- CEAB ist die π-Gegenrotation von ABCE um den Anker M. -/
theorem ellipsePoint_chirality_pi_shift (M θ : ℝ) :
    ellipsePointChiral M Chirality.CEAB θ =
      Complex.ofReal (2 * M) - ellipsePointChiral M Chirality.ABCE θ := by
  simpa [ellipsePointChiral, chiralityPhaseShift_eq_pi] using ellipsePoint_neg_shift M θ

theorem ellipsePointChiral_ABCE_zero (M : ℝ) :
    ellipsePointChiral M Chirality.ABCE 0 = Complex.ofReal (M + a_PV) := by
  simp [ellipsePointChiral, ellipsePoint, ellipsePointR, a_PV_eq_four, Real.cos_zero, Real.sin_zero]

theorem ellipsePointChiral_CEAB_zero (M : ℝ) :
    ellipsePointChiral M Chirality.CEAB 0 = Complex.ofReal (M - a_PV) := by
  rw [ellipsePoint_chirality_pi_shift, ellipsePointChiral_ABCE_zero]
  simp [a_PV_eq_four]
  ring_nf

theorem ellipsePointChiral_zero_differ_by_twice_a (M : ℝ) :
    ellipsePointChiral M Chirality.ABCE 0 - ellipsePointChiral M Chirality.CEAB 0 =
      Complex.ofReal (2 * a_PV) := by
  rw [ellipsePointChiral_ABCE_zero, ellipsePointChiral_CEAB_zero]
  simp [a_PV_eq_four]
  ring_nf

/-- ABCE-Ellipse bei markiertem Start p ≡ 5 (mod 12). -/
theorem ellipsePointMarked_ABCE (p : ℕ) (θ : ℝ) (h : chiralityWord p = some Chirality.ABCE) :
    ellipsePointMarked p θ = some (ellipsePointChiral (center p : ℝ) Chirality.ABCE θ) := by
  simp [ellipsePointMarked, h]

/-- CEAB-Ellipse bei markiertem Start p ≡ 11 (mod 12). -/
theorem ellipsePointMarked_CEAB (p : ℕ) (θ : ℝ) (h : chiralityWord p = some Chirality.CEAB) :
    ellipsePointMarked p θ = some (ellipsePointChiral (center p : ℝ) Chirality.CEAB θ) := by
  simp [ellipsePointMarked, h]

/-- Chiralitätswort bestimmt die Ellipsenparametrisierung (π-Verschiebung für CEAB). -/
theorem chiralityWord_determines_phase_shift (p : ℕ) (θ : ℝ) :
    chiralityWord p = some Chirality.CEAB →
      ellipsePointMarked p θ =
        some (ellipsePointChiral (center p : ℝ) Chirality.ABCE (θ + chiralityPhaseShift)) := by
  intro hCE
  rw [ellipsePointMarked_CEAB p θ hCE, ellipsePointChiral_CEAB_eq]

/-- Haupttheorem (Ebene B′): Kepler-Ellipse ist phasentreue geometrische Repräsentation
    der mod-12-Chiralität aus `Mod12.chiralityWord_bijection`. -/
theorem kepler_ellipse_chirality_representation (p : ℕ) (θ : ℝ) :
    (chiralityWord p = some Chirality.ABCE →
      ellipsePointMarked p θ = some (ellipsePointChiral (center p : ℝ) Chirality.ABCE θ)) ∧
    (chiralityWord p = some Chirality.CEAB →
      ellipsePointMarked p θ =
        some (ellipsePointChiral (center p : ℝ) Chirality.CEAB θ) ∧
      ellipsePointChiral (center p : ℝ) Chirality.CEAB θ =
        ellipsePointChiral (center p : ℝ) Chirality.ABCE (θ + chiralityPhaseShift)) := by
  refine ⟨?_, ?_⟩
  · intro hAB
    simpa using ellipsePointMarked_ABCE p θ hAB
  · intro hCE
    refine ⟨?_, ?_⟩
    · simpa using ellipsePointMarked_CEAB p θ hCE
    · exact ellipsePointChiral_CEAB_eq (center p : ℝ) θ

/-- Für Primvierlinge: Bijektion Startrestklasse ↔ Ellipsenchiralität (via `chiralityWord`). -/
theorem kepler_ellipse_prime_quadruplet (p : ℕ) (hq : IsPrimeQuadruplet p) (h : 3 < p) :
    chiralityWord p = some Chirality.ABCE ∨ chiralityWord p = some Chirality.CEAB := by
  exact chiralityWord_of_prime_quadruplet p hq h

theorem chiralityWord_mod5_ellipse_ABCE (p : ℕ) (θ : ℝ) (h : p % 12 = 5) :
    ellipsePointMarked p θ = some (ellipsePointChiral (center p : ℝ) Chirality.ABCE θ) := by
  simp only [ellipsePointMarked, chiralityWord_eq_ABCE_of_mod5 p h, Option.map_some]

theorem chiralityWord_mod11_ellipse_CEAB (p : ℕ) (θ : ℝ) (h : p % 12 = 11) :
    ellipsePointMarked p θ = some (ellipsePointChiral (center p : ℝ) Chirality.CEAB θ) := by
  simp only [ellipsePointMarked, chiralityWord_eq_CEAB_of_mod11 p h, Option.map_some]

end Primvierling
end CollatzEabc
