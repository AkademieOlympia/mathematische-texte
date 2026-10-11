import Riemann.E8ThetaHolomorphy

open Complex
open scoped Real

noncomputable section

namespace E8ArithmeticBridge

open E8PQC

/-!
Spitzenwert `Θ(iy) → 1` für `y → +∞`.

Für `v ≠ 0` ist `‖v‖²` eine positive gerade ganze Zahl, also mindestens `2`.
Der Rest nach Abspalten von `v = 0` wird durch `exp(-2 π (y - y₀))` mal die
Gauß-Summe zum festen `y₀` begrenzt. Darstellungszahlen kommen nicht vor.
-/

lemma e8LatticeCoord_zero : (e8LatticeCoord 0 : Fin 8 → ℝ) = 0 := by
  simp [e8LatticeCoord]

lemma e8ToEuclidean_eq_zero_iff (x : Fin 8 → ℝ) : e8ToEuclidean x = 0 ↔ x = 0 := by
  unfold e8ToEuclidean
  constructor
  · intro h
    simpa [WithLp.ofLp_toLp] using congrArg WithLp.ofLp h
  · rintro rfl
    simp

lemma e8NormSq_eq_zero_iff (x : Fin 8 → ℝ) : e8NormSq x = 0 ↔ x = 0 := by
  rw [e8NormSq_eq_norm_sq, sq_eq_zero_iff, norm_eq_zero, e8ToEuclidean_eq_zero_iff]

/-- Für `v ≠ 0` ist `‖v‖²` eine positive gerade ganze Zahl, also mindestens `2`. -/
lemma e8_nonzero_normSq_ge_two (v : E8LatticeSet) (hv : (v : Fin 8 → ℝ) ≠ 0) :
    2 ≤ e8NormSq (v : Fin 8 → ℝ) := by
  obtain ⟨k, hk⟩ := e8_even_norm v
  have hne : e8NormSq (v : Fin 8 → ℝ) ≠ 0 := (e8NormSq_eq_zero_iff _).not.mpr hv
  have hpos : 0 < e8NormSq (v : Fin 8 → ℝ) :=
    lt_of_le_of_ne (e8NormSq_nonneg _) hne.symm
  rw [hk] at hpos ⊢
  have hkpos : 0 < k := Int.cast_pos.mp <| by
    have : 0 < (2 : ℝ) * (k : ℝ) := by simpa [Int.cast_mul] using hpos
    exact (mul_pos_iff_of_pos_left two_pos).mp this
  have hk1 : (1 : ℤ) ≤ k := (Int.add_one_le_iff).mpr hkpos
  calc
    (2 : ℝ) = (2 : ℝ) * (1 : ℝ) := by ring
    _ ≤ (2 : ℝ) * (k : ℝ) :=
      mul_le_mul_of_nonneg_left (show (1 : ℝ) ≤ (k : ℝ) by exact_mod_cast hk1) two_pos.le
    _ = ((2 * k : ℤ) : ℝ) := by simp [Int.cast_mul]

lemma e8Theta_I_mul (y : ℝ) (hy : 0 < y) :
    e8Theta ⟨(y : ℂ) * I, by
      rw [mul_im, ofReal_re, ofReal_im, I_re, I_im]
      simpa using hy⟩ =
      ∑' v : E8LatticeSet,
        (Real.exp (-(Real.pi * y) * e8NormSq (v : Fin 8 → ℝ)) : ℂ) := by
  simp only [e8Theta]
  refine tsum_congr fun v => ?_
  rw [Complex.ofReal_exp]
  congr 1
  push_cast
  calc
    ↑π * I * (↑y * I) * ↑(e8NormSq (v : Fin 8 → ℝ))
        = ↑π * ↑y * (I * I) * ↑(e8NormSq (v : Fin 8 → ℝ)) := by ring
    _ = ↑π * ↑y * (-1) * ↑(e8NormSq (v : Fin 8 → ℝ)) := by rw [I_mul_I]
    _ = -(↑π * ↑y) * ↑(e8NormSq (v : Fin 8 → ℝ)) := by ring

lemma e8_cusp_tail_le (y0 y : ℝ) (_hy0 : 0 < y0) (hy : y0 ≤ y)
    (v : {x : E8LatticeSet // x ∉ ({e8LatticeCoord 0} : Finset E8LatticeSet)}) :
    Real.exp (-(Real.pi * y) * e8NormSq (v : Fin 8 → ℝ)) ≤
      Real.exp (-(2 * Real.pi) * (y - y0)) *
        Real.exp (-(Real.pi * y0) * e8NormSq (v : Fin 8 → ℝ)) := by
  have hv : (v : Fin 8 → ℝ) ≠ 0 := by
    intro h
    have : (v : E8LatticeSet) = e8LatticeCoord 0 := by
      apply Subtype.ext
      simpa [e8LatticeCoord_zero] using h
    exact v.property (by simpa using this)
  have hN : 2 ≤ e8NormSq (v : Fin 8 → ℝ) := e8_nonzero_normSq_ge_two _ hv
  have hsplit : Real.exp (-(Real.pi * y) * e8NormSq (v : Fin 8 → ℝ)) =
      Real.exp (-(Real.pi * (y - y0)) * e8NormSq (v : Fin 8 → ℝ)) *
        Real.exp (-(Real.pi * y0) * e8NormSq (v : Fin 8 → ℝ)) := by
    rw [← Real.exp_add]
    congr 1
    ring
  rw [hsplit]
  refine mul_le_mul_of_nonneg_right ?_ (Real.exp_nonneg _)
  refine Real.exp_le_exp.mpr ?_
  have hdelta : 0 ≤ y - y0 := sub_nonneg.mpr hy
  have hcoef : -(Real.pi * (y - y0)) ≤ 0 := by nlinarith [Real.pi_pos, hdelta]
  have hmul := mul_le_mul_of_nonpos_left hN hcoef
  linarith

lemma e8Theta_cusp_norm_le (y0 y : ℝ) (hy0 : 0 < y0) (hy : y0 ≤ y) :
    ‖e8Theta ⟨(y : ℂ) * I, by
        rw [mul_im, ofReal_re, ofReal_im, I_re, I_im]
        linarith⟩ - 1‖ ≤
      Real.exp (-(2 * Real.pi) * (y - y0)) *
        ∑' v : {x : E8LatticeSet // x ∉ ({e8LatticeCoord 0} : Finset E8LatticeSet)},
          Real.exp (-(Real.pi * y0) * e8NormSq (v : Fin 8 → ℝ)) := by
  have hypos : 0 < y := lt_of_lt_of_le hy0 hy
  have hfy := e8_lattice_gaussian_summable hypos
  have hfy0 := e8_lattice_gaussian_summable hy0
  let s : Finset E8LatticeSet := {e8LatticeCoord 0}
  have hzero :
      Real.exp (-(Real.pi * y) * e8NormSq (e8LatticeCoord 0 : Fin 8 → ℝ)) = 1 := by
    simp [e8LatticeCoord_zero, e8NormSq, Real.exp_zero]
  have hsplit := hfy.sum_add_tsum_subtype_compl s
  rw [Finset.sum_singleton, hzero] at hsplit
  have hΘ := e8Theta_I_mul y hypos
  rw [hΘ, ← Complex.ofReal_tsum
    (fun v : E8LatticeSet => Real.exp (-(Real.pi * y) * e8NormSq (v : Fin 8 → ℝ)))]
  have htail :
      (∑' v : E8LatticeSet, Real.exp (-(Real.pi * y) * e8NormSq (v : Fin 8 → ℝ))) =
        1 + ∑' v : {x // x ∉ s},
          Real.exp (-(Real.pi * y) * e8NormSq (v : Fin 8 → ℝ)) := by
    exact hsplit.symm
  rw [htail, Complex.ofReal_add, Complex.ofReal_one, add_sub_cancel_left]
  rw [Complex.norm_real, Real.norm_eq_abs,
    abs_of_nonneg (tsum_nonneg fun v => Real.exp_nonneg _)]
  have htail_y : Summable fun v : {x // x ∉ s} =>
      Real.exp (-(Real.pi * y) * e8NormSq (v : Fin 8 → ℝ)) :=
    hfy.subtype (· ∉ s)
  have htail_y0 : Summable fun v : {x // x ∉ s} =>
      Real.exp (-(Real.pi * y0) * e8NormSq (v : Fin 8 → ℝ)) :=
    hfy0.subtype (· ∉ s)
  have hdom : Summable fun v : {x // x ∉ s} =>
      Real.exp (-(2 * Real.pi) * (y - y0)) *
        Real.exp (-(Real.pi * y0) * e8NormSq (v : Fin 8 → ℝ)) :=
    htail_y0.mul_left _
  have hle := Summable.tsum_le_tsum (fun v => e8_cusp_tail_le y0 y hy0 hy v) htail_y hdom
  rw [tsum_mul_left] at hle
  simpa [s] using hle

/-- `Θ(iy) → 1` für `y → +∞`. -/
theorem e8Theta_limit_at_cusp :
    Filter.Tendsto
      (fun y : {y : ℝ // 0 < y} =>
        e8Theta ⟨((y : ℝ) : ℂ) * I, by
          rw [mul_im, ofReal_re, ofReal_im, I_re, I_im]
          simpa using y.2⟩)
      Filter.atTop (nhds 1) := by
  let θ : {y : ℝ // 0 < y} → ℂ :=
    fun y => e8Theta ⟨((y : ℝ) : ℂ) * I, by
      rw [mul_im, ofReal_re, ofReal_im, I_re, I_im]
      simpa using y.2⟩
  let S : ℝ :=
    ∑' v : {x : E8LatticeSet // x ∉ ({e8LatticeCoord 0} : Finset E8LatticeSet)},
      Real.exp (-(Real.pi * 1) * e8NormSq (v : Fin 8 → ℝ))
  have hupper : Filter.Tendsto
      (fun y : {y : ℝ // 0 < y} =>
        S * Real.exp (-(2 * Real.pi) * ((y : ℝ) - 1)))
      Filter.atTop (nhds 0) := by
    have hlin : Filter.Tendsto (fun t : ℝ => (2 * Real.pi) * (t - 1))
        Filter.atTop Filter.atTop := by
      have hsub : Filter.Tendsto (fun t : ℝ => t - 1) Filter.atTop Filter.atTop := by
        rw [Filter.tendsto_atTop]
        intro b
        refine Filter.eventually_atTop.mpr ⟨b + 1, fun t ht => by linarith⟩
      exact hsub.const_mul_atTop (mul_pos two_pos Real.pi_pos)
    have hexp : Filter.Tendsto
        (fun t : ℝ => Real.exp (-((2 * Real.pi) * (t - 1))))
        Filter.atTop (nhds 0) :=
      (Real.tendsto_exp_neg_atTop_nhds_zero).comp hlin
    have hcoe : Filter.Tendsto (fun y : {y : ℝ // 0 < y} => (y : ℝ))
        Filter.atTop Filter.atTop := by
      rw [Filter.tendsto_atTop]
      intro b
      refine Filter.eventually_atTop.mpr ⟨⟨max b 1, by positivity⟩, ?_⟩
      intro y hy
      exact le_trans (le_max_left _ _) hy
    have hmul := (tendsto_const_nhds (x := S)).mul (hexp.comp hcoe)
    rw [mul_zero] at hmul
    refine hmul.congr fun y => ?_
    simp [Function.comp_apply]
  have hnorm : Filter.Tendsto (fun y => ‖θ y - 1‖) Filter.atTop (nhds 0) := by
    refine tendsto_of_tendsto_of_tendsto_of_le_of_le'
      tendsto_const_nhds hupper
      (Filter.Eventually.of_forall fun _ => norm_nonneg _) ?_
    filter_upwards [Filter.eventually_ge_atTop ⟨(1 : ℝ), one_pos⟩] with y hy
    have hy1 : 1 ≤ (y : ℝ) := hy
    simpa [θ, S, mul_comm] using e8Theta_cusp_norm_le 1 (y : ℝ) one_pos hy1
  have hsub : Filter.Tendsto (fun y => θ y - 1) Filter.atTop (nhds 0) :=
    tendsto_zero_iff_norm_tendsto_zero.mpr hnorm
  have hadd := hsub.add (tendsto_const_nhds (x := (1 : ℂ)))
  simp only [zero_add] at hadd
  refine hadd.congr fun y => ?_
  simp [θ, sub_add_cancel] 

end E8ArithmeticBridge
