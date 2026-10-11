import Mathlib.NumberTheory.ModularForms.LevelOne.DimensionFormula
import Riemann.E8ThetaSlash

open Complex UpperHalfPlane
open scoped MatrixGroups

noncomputable section

namespace E8ArithmeticBridge

open ModularForm EisensteinSeries

/-!
`Θ_{E8} = E₄`.

`e8Theta_modularForm` und `E₄` liegen in `ModularForm 𝒮ℒ 4`.
Mathlib beweist, dass dieser Raum Rang 1 hat und dass `E₄` an `i∞` gegen `1` geht.
Der imaginäre Strahl `iy` liegt in diesem Spitzenfilter, und dort geht `Θ_{E8}`
nach `e8Theta_limit_at_cusp` ebenfalls gegen `1`. Der Skalar ist damit `1`.
-/

lemma e8Theta_modularForm_apply (z : ℍ) : e8Theta_modularForm z = e8Theta z := rfl

/-- Die Theta-Reihe dieses Gitters ist die normierte Eisenstein-Reihe vom Gewicht 4. -/
theorem e8_theta_eq_E4 (τ : ℍ) : e8Theta τ = E₄ τ := by
  have hne : E₄ ≠ 0 := E_ne_zero (by norm_num) ⟨2, rfl⟩
  have hspan : ∀ w : ModularForm 𝒮ℒ 4, ∃ c : ℂ, c • E₄ = w :=
    (finrank_eq_one_iff_of_nonzero' E₄ hne).mp
      (Module.rank_eq_one_iff_finrank_eq_one.mp levelOne_weight_four_rank_one)
  obtain ⟨c, hc⟩ := hspan e8Theta_modularForm
  let ι : {y : ℝ // 0 < y} → ℍ := fun y =>
    ⟨((y : ℝ) : ℂ) * Complex.I, by
      rw [mul_im, ofReal_re, ofReal_im, Complex.I_re, Complex.I_im]
      simpa using y.2⟩
  have hι : Filter.Tendsto ι Filter.atTop atImInfty := by
    intro s hs
    rw [atImInfty_mem] at hs
    obtain ⟨A, hA⟩ := hs
    refine Filter.eventually_atTop.mpr ⟨⟨max A 1, by positivity⟩, fun y hy => ?_⟩
    apply hA
    have him : (ι y).im = (y : ℝ) := by
      dsimp [ι]
      rw [mul_im, ofReal_re, ofReal_im, Complex.I_re, Complex.I_im]
      ring
    rw [him]
    exact le_trans (le_max_left _ _) hy
  have hθ : Filter.Tendsto (fun y => e8Theta (ι y)) Filter.atTop (nhds 1) :=
    e8Theta_limit_at_cusp
  have hE : Filter.Tendsto (fun y => (E₄ : ℍ → ℂ) (ι y)) Filter.atTop (nhds 1) :=
    (tendsto_E_atImInfty (k := 4)).comp hι
  have hray (y : {y : ℝ // 0 < y}) :
      e8Theta (ι y) = c * (E₄ : ℍ → ℂ) (ι y) := by
    rw [← e8Theta_modularForm_apply]
    have hcoe := congr_fun (congrArg DFunLike.coe hc) (ι y)
    simpa [FunLike.coe_smul, smul_eq_mul] using hcoe.symm
  have hmul : Filter.Tendsto (fun y => c * (E₄ : ℍ → ℂ) (ι y))
      Filter.atTop (nhds (c * 1)) :=
    tendsto_const_nhds.mul hE
  have hc1 : c = 1 := by
    have huniq := tendsto_nhds_unique (hθ.congr hray) hmul
    simpa using huniq.symm
  have hform : e8Theta_modularForm = E₄ := by
    rw [← hc, hc1, one_smul]
  rw [← e8Theta_modularForm_apply]
  simpa using congr_fun (congrArg DFunLike.coe hform) τ

end E8ArithmeticBridge
