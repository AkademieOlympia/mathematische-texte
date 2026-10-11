import Riemann.E8ThetaCusp

open Complex UpperHalfPlane
open scoped Manifold MatrixGroups ModularForm Real

noncomputable section

namespace E8ArithmeticBridge

open E8PQC ModularForm ModularGroup

/-!
Slash-Invarianz von `Θ_{E8}` unter `SL(2, ℤ)`.

Mathlib erzeugt `𝒮ℒ` durch `S` und `T`. `slash_action_generators_SL2Z` transportiert
die beiden Identitäten `e8Theta_S` und `e8_theta_T` auf jede Matrix.
Die Betragsschranke für `Im τ ≥ 1` ist dieselbe Gauß-Majorante wie an der Spitze:
`|exp(π i τ ‖v‖²)| = exp(-π (Im τ) ‖v‖²)`.
-/

lemma e8Theta_S_slash : e8Theta ∣[(4 : ℤ)] S = e8Theta := by
  ext τ
  rw [SL_slash_apply, modular_S_smul, denom_S, e8Theta_S]
  rw [show ((τ : ℂ) ^ (4 : ℤ) * e8Theta τ) * (τ : ℂ) ^ (-4 : ℤ) =
      e8Theta τ * ((τ : ℂ) ^ (4 : ℤ) * (τ : ℂ) ^ (-4 : ℤ)) by ring]
  rw [← zpow_add₀ τ.ne_zero, show (4 : ℤ) + -4 = 0 by norm_num, zpow_zero, mul_one]

lemma e8Theta_T_slash : e8Theta ∣[(4 : ℤ)] T = e8Theta := by
  ext τ
  rw [SL_slash_apply, modular_T_smul]
  have hpt : (1 : ℝ) +ᵥ τ =
      ⟨(τ : ℂ) + 1, by simpa [add_im, one_im] using τ.im_pos⟩ := by
    ext
    simp [coe_vadd, add_comm]
  rw [hpt, e8_theta_T, denom_apply, T]
  simp

/-- Die Slash-Wirkung vom Gewicht 4 ist auf ganz `SL(2, ℤ)` die Identität. -/
theorem e8Theta_slash_action_eq (γ : SL(2, ℤ)) : e8Theta ∣[(4 : ℤ)] γ = e8Theta :=
  SlashInvariantForm.slash_action_generators_SL2Z e8Theta_S_slash e8Theta_T_slash γ

lemma e8Theta_summand_norm (v : E8LatticeSet) (τ : ℍ) :
    ‖cexp (π * Complex.I * (τ : ℂ) * (e8NormSq (v : Fin 8 → ℝ) : ℂ))‖ =
      Real.exp (-(Real.pi * τ.im) * e8NormSq (v : Fin 8 → ℝ)) := by
  rw [Complex.norm_exp]
  congr 1
  let N : ℂ := (e8NormSq (v : Fin 8 → ℝ) : ℂ)
  have hw : (Complex.I * (τ : ℂ)).re = -τ.im := by
    simp [mul_re, Complex.I_re, Complex.I_im, coe_re, coe_im]
  have hwim : (Complex.I * (τ : ℂ)).im = τ.re := by
    simp [mul_im, Complex.I_re, Complex.I_im, coe_re, coe_im]
  have hπ : (↑Real.pi * (Complex.I * (τ : ℂ))).re = -(Real.pi * τ.im) := by
    rw [mul_re, ofReal_re, ofReal_im, hw, hwim]
    ring
  have hπim : (↑Real.pi * (Complex.I * (τ : ℂ))).im = Real.pi * τ.re := by
    rw [mul_im, ofReal_re, ofReal_im, hw, hwim]
    ring
  rw [show π * Complex.I * (τ : ℂ) * N = (↑Real.pi * (Complex.I * (τ : ℂ))) * N by ring]
  rw [mul_re, hπ, hπim]
  simp [N, ofReal_re, ofReal_im]

lemma e8Theta_norm_sub_one_le (y0 : ℝ) (hy0 : 0 < y0) (τ : ℍ) (hτ : y0 ≤ τ.im) :
    ‖e8Theta τ - 1‖ ≤
      Real.exp (-(2 * Real.pi) * (τ.im - y0)) *
        ∑' v : {x : E8LatticeSet // x ∉ ({e8LatticeCoord 0} : Finset E8LatticeSet)},
          Real.exp (-(Real.pi * y0) * e8NormSq (v : Fin 8 → ℝ)) := by
  let s : Finset E8LatticeSet := {e8LatticeCoord 0}
  let f : E8LatticeSet → ℂ := fun v =>
    cexp (π * Complex.I * (τ : ℂ) * (e8NormSq (v : Fin 8 → ℝ) : ℂ))
  have him : 0 < τ.im := lt_of_lt_of_le hy0 hτ
  have hsummable : Summable f := by
    refine Summable.of_norm_bounded (e8_lattice_gaussian_summable him) ?_
    intro v
    exact le_of_eq (e8Theta_summand_norm v τ)
  have hzero : f (e8LatticeCoord 0) = 1 := by
    simp [f, e8LatticeCoord_zero, e8NormSq, Complex.exp_zero]
  have hsplit := hsummable.sum_add_tsum_subtype_compl s
  rw [Finset.sum_singleton, hzero] at hsplit
  have hΘ : e8Theta τ = ∑' v, f v := by simp [e8Theta, f]
  rw [hΘ, ← hsplit, add_sub_cancel_left]
  refine (norm_tsum_le_tsum_norm ?_).trans ?_
  · exact (hsummable.subtype (· ∉ s)).norm
  · have htail : Summable fun v : {x // x ∉ s} =>
        Real.exp (-(Real.pi * τ.im) * e8NormSq (v : Fin 8 → ℝ)) :=
      (e8_lattice_gaussian_summable him).subtype (· ∉ s)
    have htail0 : Summable fun v : {x // x ∉ s} =>
        Real.exp (-(Real.pi * y0) * e8NormSq (v : Fin 8 → ℝ)) :=
      (e8_lattice_gaussian_summable hy0).subtype (· ∉ s)
    have hdom : Summable fun v : {x // x ∉ s} =>
        Real.exp (-(2 * Real.pi) * (τ.im - y0)) *
          Real.exp (-(Real.pi * y0) * e8NormSq (v : Fin 8 → ℝ)) :=
      htail0.mul_left _
    have hle := Summable.tsum_le_tsum
      (fun v => e8_cusp_tail_le y0 τ.im hy0 hτ v) htail hdom
    rw [tsum_mul_left] at hle
    have hnorm : (∑' v : {x // x ∉ s}, ‖f (v : E8LatticeSet)‖) =
        ∑' v : {x // x ∉ s},
          Real.exp (-(Real.pi * τ.im) * e8NormSq (v : Fin 8 → ℝ)) := by
      refine tsum_congr fun v => ?_
      simpa [f] using e8Theta_summand_norm (v : E8LatticeSet) τ
    rw [hnorm]
    simpa [s] using hle

lemma e8Theta_isBoundedAtImInfty : IsBoundedAtImInfty e8Theta := by
  rw [isBoundedAtImInfty_iff]
  let S : ℝ :=
    ∑' v : {x : E8LatticeSet // x ∉ ({e8LatticeCoord 0} : Finset E8LatticeSet)},
      Real.exp (-(Real.pi * 1) * e8NormSq (v : Fin 8 → ℝ))
  have hS : 0 ≤ S := tsum_nonneg fun _ => Real.exp_nonneg _
  refine ⟨1 + S, 1, fun τ hτ => ?_⟩
  have hle := e8Theta_norm_sub_one_le 1 one_pos τ hτ
  have hexp : Real.exp (-(2 * Real.pi) * (τ.im - 1)) ≤ 1 :=
    Real.exp_le_one_iff.mpr <| by nlinarith [Real.pi_pos, hτ]
  calc
    ‖e8Theta τ‖ = ‖e8Theta τ - 1 + 1‖ := by ring_nf
    _ ≤ ‖e8Theta τ - 1‖ + ‖(1 : ℂ)‖ := norm_add_le _ _
    _ = ‖e8Theta τ - 1‖ + 1 := by simp
    _ ≤ Real.exp (-(2 * Real.pi) * (τ.im - 1)) * S + 1 := by
      simpa [S] using add_le_add_right hle 1
    _ ≤ 1 * S + 1 := by gcongr
    _ = 1 + S := by ring

lemma e8Theta_mdifferentiable : MDiff e8Theta := by
  rw [mdifferentiable_iff]
  refine e8Theta_differentiableOn.congr fun z hz => ?_
  rw [Function.comp_apply, ofComplex_apply_of_im_pos hz, e8Theta]
  simp [e8Theta_term]

/-- `Θ_{E8}` als Modulform vom Gewicht 4 und Stufe 1. -/
def e8Theta_modularForm : ModularForm 𝒮ℒ 4 where
  toFun := e8Theta
  slash_action_eq' γ hγ := by
    obtain ⟨γ, rfl⟩ := hγ
    exact e8Theta_slash_action_eq γ
  holo' := e8Theta_mdifferentiable
  bdd_at_cusps' hc := by
    obtain ⟨γ, rfl⟩ := isCusp_SL2Z_iff'.mp hc
    rw [OnePoint.isBoundedAt_iff_exists_SL2Z hc]
    refine ⟨γ, rfl, ?_⟩
    rw [e8Theta_slash_action_eq]
    exact e8Theta_isBoundedAtImInfty

end E8ArithmeticBridge
