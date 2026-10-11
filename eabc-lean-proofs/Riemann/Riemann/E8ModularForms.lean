import Mathlib

open Complex

open scoped Real FourierTransform RealInnerProductSpace UpperHalfPlane

noncomputable section

namespace E8ModularForms

/-!
G3A-M.1. Die achtdimensionale Gauß-Transformation für `Re a > 0` ist
`fourier_gaussian_innerProductSpace` mit `b = π a`. Der Exponent `finrank / 2`
ist die natürliche Zahl `4`. Die Quadratwurzel bleibt in der eindimensionalen
Mathlib-Rechnung.

`poisson_summation_e8` und die S-Transformation auf ganz `ℍ` stehen hier nicht.
-/

/-- `ℱ[x ↦ exp(-π a ‖x‖²)](ξ) = a^{-4} exp(-π ‖ξ‖² / a)` für `0 < Re a`. -/
theorem complexGaussian_fourier_8 (a : ℂ) (ha : 0 < a.re) (ξ : EuclideanSpace ℝ (Fin 8)) :
    𝓕 (fun x : EuclideanSpace ℝ (Fin 8) => cexp (-((π : ℂ) * a) * ‖x‖ ^ 2)) ξ =
      a ^ (-4 : ℤ) * cexp (-π * ‖ξ‖ ^ 2 / a) := by
  have hb : 0 < ((π : ℂ) * a).re := by
    rw [mul_re, ofReal_re, ofReal_im, zero_mul, sub_zero]
    exact mul_pos Real.pi_pos ha
  have ha0 : a ≠ 0 := by
    rintro rfl
    simp at ha
  rw [fourier_gaussian_innerProductSpace (b := (π : ℂ) * a) hb, finrank_euclideanSpace_fin]
  have hexp : ((8 : ℕ) / 2 : ℂ) = 4 := by norm_num
  have hdiv : (π : ℂ) / ((π : ℂ) * a) = a⁻¹ := by
    field_simp [ha0]
  have harg : -π ^ 2 * ‖ξ‖ ^ 2 / ((π : ℂ) * a) = -π * ‖ξ‖ ^ 2 / a := by
    field_simp [ha0]
  rw [hexp, hdiv, harg]
  congr 1
  rw [Complex.cpow_ofNat, inv_pow, ← Complex.cpow_ofNat, ← Complex.cpow_neg]
  rw [show (-(4 : ℂ)) = ((-4 : ℤ) : ℂ) by norm_num, Complex.cpow_intCast]

/-- `ℱ[x ↦ exp(π i τ ‖x‖²)](ξ) = τ^{-4} exp(-π i ‖ξ‖² / τ)` für `τ ∈ ℍ`. -/
theorem thetaGaussian_fourier_8 (τ : ℍ) (ξ : EuclideanSpace ℝ (Fin 8)) :
    𝓕 (fun x : EuclideanSpace ℝ (Fin 8) => cexp (π * I * (τ : ℂ) * ‖x‖ ^ 2)) ξ =
      (τ : ℂ) ^ (-4 : ℤ) * cexp (-π * I * ‖ξ‖ ^ 2 / (τ : ℂ)) := by
  have ha : 0 < (-I * (τ : ℂ)).re := by
    rw [neg_mul, neg_re, mul_re, I_re, I_im, zero_mul, zero_sub, neg_neg, one_mul]
    exact τ.im_pos
  have h := complexGaussian_fourier_8 (-I * (τ : ℂ)) ha ξ
  have hfun : (fun x : EuclideanSpace ℝ (Fin 8) => cexp (π * I * (τ : ℂ) * ‖x‖ ^ 2)) =
      fun x => cexp (-((π : ℂ) * (-I * (τ : ℂ))) * ‖x‖ ^ 2) := by
    ext x
    congr 1
    ring
  rw [hfun, h]
  have hτ : (τ : ℂ) ≠ 0 := τ.ne_zero
  have hI : (-I) ≠ 0 := by simp
  have hI4 : (-I) ^ (4 : ℕ) = 1 := by
    have h2 : (-I) ^ (2 : ℕ) = -1 := by rw [pow_two, neg_mul_neg, I_mul_I]
    rw [show (4 : ℕ) = 2 * 2 by norm_num, pow_mul, h2]
    norm_num
  congr 1
  · rw [← Complex.cpow_intCast, ← Complex.cpow_intCast]
    rw [show ((-4 : ℤ) : ℂ) = -(4 : ℂ) by norm_num]
    rw [Complex.cpow_neg, Complex.cpow_neg, Complex.cpow_ofNat, Complex.cpow_ofNat, mul_pow,
      hI4, one_mul]
  · congr 1
    field_simp [hτ, hI]
    rw [I_sq]
    ring

end E8ModularForms
