import Riemann.E8ArithmeticBridge

open Complex
open scoped Real

noncomputable section

namespace E8ArithmeticBridge

open E8PQC

/-!
Holomorphie von `Θ_{E8}` auf `{z | 0 < Im z}`.

Die Ableitung eines Summanden ist `π i ‖v‖² exp(π i z ‖v‖²)`.
Auf `{Im z > ε}` majorisiert `π ‖v‖² exp(-π ε ‖v‖²)` diese Ableitung.
`x exp(-a x)` bleibt beschränkt, also reicht die Gauß-Summe zum Parameter `ε / 2`.
Die Reihe ist dann auf jeder solchen Halbebene differenzierbar durch
`hasDerivAt_tsum_of_isPreconnected`.
-/

lemma e8NormSq_nonneg (x : Fin 8 → ℝ) : 0 ≤ e8NormSq x :=
  Finset.sum_nonneg fun _ _ => sq_nonneg _

lemma e8Gaussian_eq_exp (t : ℝ) (x : EuclideanSpace ℝ (Fin 8)) :
    e8Gaussian t x = (Real.exp (-(Real.pi * t) * ‖x‖ ^ 2) : ℂ) := by
  unfold e8Gaussian
  rw [Complex.ofReal_exp]
  congr 1
  push_cast
  ring

/-- Die reelle Gauß-Summe über `E₈` konvergiert für `t > 0`. -/
lemma e8_lattice_gaussian_summable {t : ℝ} (ht : 0 < t) :
    Summable fun v : E8LatticeSet =>
      Real.exp (-(Real.pi * t) * e8NormSq (v : Fin 8 → ℝ)) := by
  have hR : RapidDecay (e8Gaussian t) :=
    rapidDecay_of_ofReal fun k => e8Gaussian_decay ht k
  have hF : RapidDecay (e8Gaussian t ∘ e8CoordCLE.toContinuousLinearMap) :=
    rapidDecay_comp_leftInverse e8CoordCLE.toContinuousLinearMap
      e8CoordCLE.symm.toContinuousLinearMap e8CoordCLE_symm_comp
      e8CoordCLE.toHomeomorph.isClosedEmbedding.tendsto_cocompact hR
  have hFeq : e8Gaussian t ∘ e8CoordCLE.toContinuousLinearMap =
      e8Gaussian t ∘ e8CoordMapOnEuclidean := by
    ext x
    simp [e8CoordCLE_apply]
  have hdecay : RapidDecay (e8Gaussian t ∘ e8CoordMapOnEuclidean) := by
    rwa [hFeq] at hF
  have hsumC := rapidDecay_summable hdecay
  have hpoint (n : Fin 8 → ℤ) :
      (e8Gaussian t ∘ e8CoordMapOnEuclidean) (zdPoint n) =
        (Real.exp (-(Real.pi * t) *
          e8NormSq (e8LatticeCoord n : Fin 8 → ℝ)) : ℂ) := by
    rw [Function.comp_apply, e8CoordMapOnEuclidean_zdPoint, e8Gaussian_eq_exp,
      e8NormSq_eq_norm_sq]
    rfl
  simp_rw [hpoint] at hsumC
  have hreal : Summable fun n : Fin 8 → ℤ =>
      Real.exp (-(Real.pi * t) * e8NormSq (e8LatticeCoord n : Fin 8 → ℝ)) :=
    (summable_ofReal
      (f := fun n : Fin 8 → ℤ =>
        Real.exp (-(Real.pi * t) * e8NormSq (e8LatticeCoord n : Fin 8 → ℝ)))).mp hsumC
  have hbij : Function.Bijective e8LatticeCoord :=
    ⟨e8LatticeCoord_injective, e8LatticeCoord_surjective⟩
  let e : (Fin 8 → ℤ) ≃ E8LatticeSet := Equiv.ofBijective e8LatticeCoord hbij
  have hcomp :
      (fun v : E8LatticeSet => Real.exp (-(Real.pi * t) * e8NormSq (v : Fin 8 → ℝ))) ∘ e =
        fun n => Real.exp (-(Real.pi * t) * e8NormSq (e8LatticeCoord n : Fin 8 → ℝ)) := by
    ext n
    simp [e]
  rw [← hcomp] at hreal
  exact (e.summable_iff).mp hreal

/-- `π r exp(-π ε r) ≤ (2 / (ε e)) exp(-π (ε/2) r)` für `0 ≤ r`. -/
lemma e8_poly_gaussian_bound (ε r : ℝ) (hε : 0 < ε) (_hr : 0 ≤ r) :
    Real.pi * r * Real.exp (-(Real.pi * ε) * r) ≤
      (2 / (ε * Real.exp 1)) * Real.exp (-(Real.pi * (ε / 2)) * r) := by
  set α' : ℝ := Real.pi * (ε / 2)
  have hα' : 0 < α' := mul_pos Real.pi_pos (half_pos hε)
  have hy : α' * r * Real.exp (-(α' * r)) ≤ Real.exp (-1) :=
    Real.mul_exp_neg_le_exp_neg_one (α' * r)
  have hr1 : r * Real.exp (-(α' * r)) ≤ Real.exp (-1) / α' := by
    have hdiv := (div_le_div_iff_of_pos_right hα').mpr hy
    simpa [mul_assoc, mul_div_cancel_left₀ _ hα'.ne'] using hdiv
  have htwo : Real.exp (-(Real.pi * ε) * r) =
      Real.exp (-(α' * r)) * Real.exp (-(α' * r)) := by
    rw [← Real.exp_add]
    congr 1
    simp only [α']
    ring
  have hconst : Real.pi * (Real.exp (-1) / α') = 2 / (ε * Real.exp 1) := by
    simp only [α', Real.exp_neg]
    field_simp
  calc
    Real.pi * r * Real.exp (-(Real.pi * ε) * r)
        = Real.pi * (r * Real.exp (-(α' * r))) * Real.exp (-(α' * r)) := by
          rw [htwo]
          ring
    _ ≤ Real.pi * (Real.exp (-1) / α') * Real.exp (-(α' * r)) := by
          gcongr
    _ = (2 / (ε * Real.exp 1)) * Real.exp (-(Real.pi * (ε / 2)) * r) := by
          rw [hconst]
          congr 1
          simp only [α']
          ring_nf

/-- Summand `exp(π i z ‖v‖²)`. -/
def e8Theta_term (v : E8LatticeSet) (z : ℂ) : ℂ :=
  cexp (π * I * z * (e8NormSq (v : Fin 8 → ℝ) : ℂ))

/-- Ableitung des Summanden: `π i ‖v‖² exp(π i z ‖v‖²)`. -/
lemma e8Theta_term_hasDerivAt (v : E8LatticeSet) (z : ℂ) :
    HasDerivAt (e8Theta_term v)
      (e8Theta_term v z * (π * I * (e8NormSq (v : Fin 8 → ℝ) : ℂ))) z := by
  unfold e8Theta_term
  convert (hasDerivAt_const_mul (π * I * (e8NormSq (v : Fin 8 → ℝ) : ℂ))).cexp using 2
  · ring_nf
  · ring_nf

lemma e8Theta_term_norm (v : E8LatticeSet) (z : ℂ) :
    ‖e8Theta_term v z‖ =
      Real.exp (-(Real.pi * z.im) * e8NormSq (v : Fin 8 → ℝ)) := by
  unfold e8Theta_term
  rw [Complex.norm_exp]
  have hIre : (I * z).re = -z.im := by
    simp [mul_re, I_re, I_im]
  have hIim : (I * z).im = z.re := by
    simp [mul_im, I_re, I_im]
  rw [show π * I * z * (e8NormSq (v : Fin 8 → ℝ) : ℂ) =
      ((Real.pi * e8NormSq (v : Fin 8 → ℝ) : ℝ) : ℂ) * (I * z) by
    have hcast : ((Real.pi * e8NormSq (v : Fin 8 → ℝ) : ℝ) : ℂ) =
        π * (e8NormSq (v : Fin 8 → ℝ) : ℂ) := by push_cast; rfl
    rw [hcast]
    ring]
  simp [mul_re, ofReal_re, ofReal_im, hIre, hIim]
  ring

lemma e8Theta_term_deriv_norm (v : E8LatticeSet) (z : ℂ) :
    ‖e8Theta_term v z * (π * I * (e8NormSq (v : Fin 8 → ℝ) : ℂ))‖ =
      Real.pi * e8NormSq (v : Fin 8 → ℝ) *
        Real.exp (-(Real.pi * z.im) * e8NormSq (v : Fin 8 → ℝ)) := by
  rw [norm_mul, e8Theta_term_norm]
  have hr : 0 ≤ e8NormSq (v : Fin 8 → ℝ) := e8NormSq_nonneg _
  have hcoeff : ‖π * I * (e8NormSq (v : Fin 8 → ℝ) : ℂ)‖ =
      Real.pi * e8NormSq (v : Fin 8 → ℝ) := by
    rw [norm_mul, norm_mul, Complex.norm_I, Complex.norm_real, Complex.norm_real,
      Real.norm_eq_abs, Real.norm_eq_abs, abs_of_nonneg Real.pi_nonneg, abs_of_nonneg hr,
      mul_one]
  rw [hcoeff]
  ring

/-- `π ‖v‖² exp(-π ε ‖v‖²)` ist über `E₈` summierbar. -/
lemma e8Theta_deriv_majorant_summable (ε : ℝ) (hε : 0 < ε) :
    Summable fun v : E8LatticeSet =>
      Real.pi * e8NormSq (v : Fin 8 → ℝ) *
        Real.exp (-(Real.pi * ε) * e8NormSq (v : Fin 8 → ℝ)) := by
  have hgauss := e8_lattice_gaussian_summable (half_pos hε)
  refine Summable.of_nonneg_of_le (fun v => ?_) (fun v => ?_)
      ((hgauss.mul_left (2 / (ε * Real.exp 1))))
  · have hr : 0 ≤ e8NormSq (v : Fin 8 → ℝ) := e8NormSq_nonneg _
    positivity
  · simpa [mul_assoc] using
      e8_poly_gaussian_bound ε (e8NormSq (v : Fin 8 → ℝ)) hε (e8NormSq_nonneg _)

lemma e8Theta_hasDerivAt_of_im_gt {ε : ℝ} (hε : 0 < ε) {z : ℂ} (hz : ε < z.im) :
    HasDerivAt (fun w : ℂ => ∑' v : E8LatticeSet, e8Theta_term v w)
      (∑' v : E8LatticeSet,
        e8Theta_term v z * (π * I * (e8NormSq (v : Fin 8 → ℝ) : ℂ))) z := by
  let t : Set ℂ := {w | ε < w.im}
  have ht : IsOpen t := isOpen_lt continuous_const Complex.continuous_im
  have hpre : IsPreconnected t := (convex_halfSpace_im_gt ε).isPreconnected
  let y0 : ℂ := (ε + 1) * I
  have hy0im : y0.im = ε + 1 := by
    simp [y0, add_im, mul_im, ofReal_re, ofReal_im, I_re, I_im, one_re, one_im]
  have hy0 : y0 ∈ t := by
    change ε < y0.im
    rw [hy0im]
    linarith
  have hy0pos : 0 < y0.im := by
    rw [hy0im]
    linarith
  have hg0 : Summable fun v : E8LatticeSet => e8Theta_term v y0 := by
    refine Summable.of_norm_bounded
      (e8_lattice_gaussian_summable hy0pos) ?_
    intro v
    rw [e8Theta_term_norm]
  have hu := e8Theta_deriv_majorant_summable ε hε
  let g' : E8LatticeSet → ℂ → ℂ :=
    fun v w => e8Theta_term v w * (π * I * (e8NormSq (v : Fin 8 → ℝ) : ℂ))
  refine hasDerivAt_tsum_of_isPreconnected (t := t) (g' := g') hu ht hpre ?_ ?_ hy0 hg0 hz
  · intro v w hw
    exact e8Theta_term_hasDerivAt v w
  · intro v w hw
    rw [e8Theta_term_deriv_norm]
    have hr : 0 ≤ e8NormSq (v : Fin 8 → ℝ) := e8NormSq_nonneg _
    have hrate : Real.pi * ε * e8NormSq (v : Fin 8 → ℝ) ≤
        Real.pi * w.im * e8NormSq (v : Fin 8 → ℝ) := by
      gcongr
      exact hw.le
    have hexp : Real.exp (-(Real.pi * w.im) * e8NormSq (v : Fin 8 → ℝ)) ≤
        Real.exp (-(Real.pi * ε) * e8NormSq (v : Fin 8 → ℝ)) := by
      refine Real.exp_le_exp.mpr ?_
      linarith
    have hπ : 0 ≤ Real.pi * e8NormSq (v : Fin 8 → ℝ) := by positivity
    exact mul_le_mul_of_nonneg_left hexp hπ

/-- `Θ_{E8}`, auf die offene obere Halbebene fortgesetzt, ist holomorph. -/
theorem e8Theta_differentiableOn :
    DifferentiableOn ℂ
      (fun z : ℂ => ∑' v : E8LatticeSet, e8Theta_term v z) {z | 0 < z.im} := by
  intro z hz
  have hε : 0 < z.im / 2 := half_pos hz
  have hzε : z.im / 2 < z.im := by linarith
  exact (e8Theta_hasDerivAt_of_im_gt hε hzε).differentiableAt.differentiableWithinAt

end E8ArithmeticBridge
