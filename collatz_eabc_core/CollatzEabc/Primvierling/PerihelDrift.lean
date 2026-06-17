/-
  CollatzEabc.Primvierling.PerihelDrift — arithmetische Perihelverschiebung (Ebene B′/C).

  Mathematische Analogie: die ideale EABC-Uhr schließt nach vier Ticks exakt (Π = 0);
  uniforme oder chirale Phasenkorrekturen modellieren eine „Perihelverschiebung“ der
  Kepler-Ellipsenuhr — ohne physikalischen oder GR-Anspruch.

  Referenz: `Projektionszeuge_Primvierling.tex`, Abschnitt „Arithmetische Perihelverschiebung“.
-/

import Mathlib.Analysis.Complex.Exponential
import Mathlib.Analysis.Complex.Trigonometric
import CollatzEabc.Primvierling.KeplerEllipse

namespace CollatzEabc
namespace Primvierling

/-- Ideale EABC-Uhr: θ_m = (π/2)·m. -/
noncomputable def idealPhase (m : ℕ) : ℝ :=
  (Real.pi / 2) * m

/-- Perihelverschiebung Π(θ_m, θ_{m+4}) = Θ_{m+4} − Θ_m − 2π. -/
noncomputable def perihelShift (θ_m θ_m4 : ℝ) : ℝ :=
  θ_m4 - θ_m - 2 * Real.pi

/-- Phasenkorrektur δ_m über den diskreten Index m. -/
structure PhaseDrift where
  δ : ℕ → ℝ

/-- Korrigierte Uhr: Θ_m = (π/2)·m + δ_m. -/
noncomputable def driftPhase (d : PhaseDrift) (m : ℕ) : ℝ :=
  idealPhase m + d.δ m

/-- Π(m) = δ_{m+4} − δ_m für korrigierte Phasen. -/
noncomputable def perihelShiftDrift (d : PhaseDrift) (m : ℕ) : ℝ :=
  d.δ (m + 4) - d.δ m

/-- Uniforme Drift ε: θ_m = m·(π/2 + ε). -/
structure UniformDrift where
  ε : ℝ

noncomputable def uniformDriftPhase (u : UniformDrift) (m : ℕ) : ℝ :=
  m * (Real.pi / 2 + u.ε)

/-- Chirale Korrekturen an den vier EABC-Ticks. -/
structure ChiralPhaseDrift where
  δ_E : ℝ
  δ_A : ℝ
  δ_B : ℝ
  δ_C : ℝ

/-- Summe der chiralen Korrekturen; Schließung genau wenn Δ_peri = 0. -/
noncomputable def chiralDriftSum (d : ChiralPhaseDrift) : ℝ :=
  d.δ_E + d.δ_A + d.δ_B + d.δ_C

def chiralDriftClosed (d : ChiralPhaseDrift) : Prop :=
  chiralDriftSum d = 0

noncomputable def chiralDriftAt (d : ChiralPhaseDrift) (t : Fin 4) : ℝ :=
  match tickFlavor t with
  | EClass.E => d.δ_E
  | EClass.A => d.δ_A
  | EClass.B => d.δ_B
  | EClass.C => d.δ_C

/-- Ein Uhrschritt um π/2 + ε als komplexe Rotation. -/
noncomputable def phaseRotationOp (θ : ℝ) : ℂ → ℂ :=
  fun z => Complex.exp (θ * Complex.I) * z

noncomputable def clockStepAngle (ε : ℝ) : ℝ :=
  Real.pi / 2 + ε

/-- Vier aufeinanderfolgende Uhrschritte: T^4 = R_{2π + 4ε}. -/
noncomputable def fourClockRotation (ε : ℝ) : ℂ → ℂ :=
  phaseRotationOp (4 * clockStepAngle ε)

/-- Kandidaten für arithmetische Drift (nur dokumentiert, nicht bewiesen). -/
structure DriftCandidate where
  /-- ε(p) ~ 1 / log p -/
  invLog : ℕ → ℝ
  /-- ε(p) ~ 1 / p -/
  invPrime : ℕ → ℝ
  /-- ε(p) ~ α · ΔW(p) -/
  deltaWeight : ℝ → ℝ
  /-- ε(p) ~ |arg Φ_pref(w) − (π/2)|w|| -/
  prefPhaseGap : ℕ → ℝ

/-- Standard-Kandidaten (heuristisch, keine Primzahltheorem-Behauptung). -/
noncomputable def defaultDriftCandidate : DriftCandidate where
  invLog := fun p => 1 / Real.log (p : ℝ)
  invPrime := fun p => 1 / (p : ℝ)
  deltaWeight := fun _ => 0
  prefPhaseGap := fun _ => 0

theorem idealPhase_eq_phaseTick (t : Fin 4) :
    idealPhase t = phaseTick t := by
  simp [idealPhase, phaseTick]

def idealPhaseDrift : PhaseDrift :=
  ⟨fun _ => 0⟩

theorem perihelShift_ideal_zero (m : ℕ) :
    perihelShift (idealPhase m) (idealPhase (m + 4)) = 0 := by
  unfold perihelShift idealPhase
  push_cast
  ring

theorem perihelShiftDrift_zero (d : PhaseDrift) (m : ℕ)
    (h : d.δ (m + 4) = d.δ m) :
    perihelShiftDrift d m = 0 := by
  unfold perihelShiftDrift
  rw [h, sub_self]

theorem perihelShiftDrift_ideal (m : ℕ) :
    perihelShiftDrift idealPhaseDrift m = 0 :=
  perihelShiftDrift_zero idealPhaseDrift m rfl

theorem perihelShift_driftPhase (d : PhaseDrift) (m : ℕ) :
    perihelShift (driftPhase d m) (driftPhase d (m + 4)) = perihelShiftDrift d m := by
  unfold perihelShift driftPhase idealPhase perihelShiftDrift
  push_cast
  ring

theorem uniform_perihel_shift (u : UniformDrift) (m : ℕ) :
    perihelShift (uniformDriftPhase u m) (uniformDriftPhase u (m + 4)) = 4 * u.ε := by
  unfold perihelShift uniformDriftPhase
  push_cast
  ring

theorem uniform_perihel_shift_ne_zero (u : UniformDrift) (m : ℕ) (hε : u.ε ≠ 0) :
    perihelShift (uniformDriftPhase u m) (uniformDriftPhase u (m + 4)) ≠ 0 := by
  rw [uniform_perihel_shift]
  exact mul_ne_zero (by norm_num) hε

theorem ideal_four_step_angle : 4 * clockStepAngle 0 = 2 * Real.pi := by
  unfold clockStepAngle
  ring

theorem phaseRotation_two_pi (z : ℂ) :
    phaseRotationOp (2 * Real.pi) z = z := by
  unfold phaseRotationOp
  rw [Complex.exp_ofReal_mul_I, Real.cos_two_pi, Real.sin_two_pi]
  simp

theorem fourClock_ideal_eq_id (z : ℂ) :
    fourClockRotation 0 z = z := by
  unfold fourClockRotation
  rw [ideal_four_step_angle, phaseRotation_two_pi]

theorem fourClock_drift_angle (ε : ℝ) :
    4 * clockStepAngle ε = 2 * Real.pi + 4 * ε := by
  unfold clockStepAngle
  ring

theorem fourClock_drift_angle_ne_identity (ε : ℝ) (hε : ε ≠ 0) :
    4 * clockStepAngle ε ≠ 2 * Real.pi := by
  unfold clockStepAngle
  intro h
  have h4 : 4 * ε = 0 := by linarith
  exact hε (by linarith)

theorem chiralDriftAt_tickFlavor (d : ChiralPhaseDrift) (t : Fin 4) :
    chiralDriftAt d t =
      match tickFlavor t with
      | EClass.E => d.δ_E
      | EClass.A => d.δ_A
      | EClass.B => d.δ_B
      | EClass.C => d.δ_C := rfl

theorem chiralDriftClosed_iff_sum_zero (d : ChiralPhaseDrift) :
    chiralDriftClosed d ↔ chiralDriftSum d = 0 :=
  Iff.rfl

end Primvierling
end CollatzEabc
