import Riemann.E8ThetaIdentification

open Complex UpperHalfPlane
open scoped Real MatrixGroups ArithmeticFunction.sigma

noncomputable section

namespace E8ArithmeticBridge

open E8PQC ModularForm EisensteinSeries ArithmeticFunction

/-!
Der letzte Baustein von Gate 3A.

`‖v‖²` ist für `v ∈ E₈` eine gerade ganze Zahl. Die Schale `‖v‖² = 2n` ist
endlich, weil die Koordinatenabbildung `e8CoordCLE` ein linearer Isomorphismus
ist: beschränkte Norm heißt beschränkte ganzzahlige Koordinaten. Die
Theta-Reihe gruppiert sich deshalb faserweise zu `∑ r_{E8}(2n) q^n`, und die
Eindeutigkeit der `q`-Entwicklung liefert die Koeffizienten. Der Vergleich mit
`E₄` ergibt `r_{E8}(2n) = 240 σ₃(n)`.
-/

lemma e8_half_exists (v : E8LatticeSet) :
    ∃ n : ℕ, e8NormSq (v : Fin 8 → ℝ) = ((2 * n : ℕ) : ℝ) := by
  obtain ⟨k, hk⟩ := e8_even_norm v
  have h0 : 0 ≤ e8NormSq (v : Fin 8 → ℝ) := e8NormSq_nonneg _
  rw [hk, Int.cast_mul, Int.cast_ofNat] at h0
  have hk0 : (0 : ℝ) ≤ (k : ℝ) :=
    (mul_nonneg_iff_of_pos_left (by norm_num : (0 : ℝ) < 2)).mp h0
  have hkint : 0 ≤ k := by exact_mod_cast hk0
  refine ⟨k.toNat, ?_⟩
  rw [hk, ← Int.toNat_of_nonneg hkint]
  norm_cast

lemma e8_half_unique (v : E8LatticeSet) {n m : ℕ}
    (hn : e8NormSq (v : Fin 8 → ℝ) = ((2 * n : ℕ) : ℝ))
    (hm : e8NormSq (v : Fin 8 → ℝ) = ((2 * m : ℕ) : ℝ)) : n = m := by
  have hnat : 2 * n = 2 * m := Nat.cast_injective (hn.symm.trans hm)
  omega

noncomputable def e8HalfIndex (v : E8LatticeSet) : ℕ :=
  Classical.choose (e8_half_exists v)

lemma e8HalfIndex_spec (v : E8LatticeSet) :
    e8NormSq (v : Fin 8 → ℝ) = ((2 * e8HalfIndex v : ℕ) : ℝ) :=
  Classical.choose_spec (e8_half_exists v)

lemma e8HalfIndex_eq_iff (v : E8LatticeSet) (n : ℕ) :
    e8HalfIndex v = n ↔ e8NormSq (v : Fin 8 → ℝ) = ((2 * n : ℕ) : ℝ) := by
  constructor
  · intro h
    rw [← h]
    exact e8HalfIndex_spec v
  · intro h
    exact e8_half_unique v (e8HalfIndex_spec v) h

/-- Ganzzahlige Schranke für die Koordinaten eines Vektors der Norm `m`. -/
def e8ShellBound (m : ℕ) : ℤ :=
  Int.ceil (‖e8CoordCLE.symm.toContinuousLinearMap‖ * Real.sqrt (m : ℝ))

lemma e8_shell_coord_mem (m : ℕ) (n : Fin 8 → ℤ)
    (hn : e8NormSq (e8LatticeCoord n : Fin 8 → ℝ) = (m : ℝ)) (i : Fin 8) :
    n i ∈ Finset.Icc (-e8ShellBound m) (e8ShellBound m) := by
  set C : ℝ := ‖e8CoordCLE.symm.toContinuousLinearMap‖
  set B : ℝ := C * Real.sqrt (m : ℝ)
  have hsqrt : ‖e8CoordCLE (zdPoint n)‖ = Real.sqrt (m : ℝ) := by
    have hx : e8CoordCLE (zdPoint n) =
        e8ToEuclidean (e8LatticeCoord n : Fin 8 → ℝ) := by
      rw [e8CoordCLE_apply, e8CoordMapOnEuclidean_zdPoint]
      simp [e8LatticeCoord]
    rw [hx, ← Real.sqrt_sq (norm_nonneg _), ← e8NormSq_eq_norm_sq, hn]
  have hop : ‖zdPoint n‖ ≤ C * ‖e8CoordCLE (zdPoint n)‖ := by
    calc ‖zdPoint n‖
        = ‖e8CoordCLE.symm (e8CoordCLE (zdPoint n))‖ := by
            rw [e8CoordCLE.symm_apply_apply]
      _ = ‖e8CoordCLE.symm.toContinuousLinearMap (e8CoordCLE (zdPoint n))‖ := by
            rfl
      _ ≤ ‖e8CoordCLE.symm.toContinuousLinearMap‖ * ‖e8CoordCLE (zdPoint n)‖ :=
            e8CoordCLE.symm.toContinuousLinearMap.le_opNorm _
      _ = C * ‖e8CoordCLE (zdPoint n)‖ := rfl
  have habs : |(n i : ℝ)| ≤ ‖zdPoint n‖ := by
    have hpi : ‖zdPoint n‖ ^ 2 = ∑ j, (n j : ℝ) ^ 2 := by
      unfold zdPoint
      rw [PiLp.norm_sq_eq_of_L2]
      simp [Real.norm_eq_abs, sq_abs]
    calc |(n i : ℝ)|
        = Real.sqrt ((n i : ℝ) ^ 2) := by rw [Real.sqrt_sq_eq_abs]
      _ ≤ Real.sqrt (‖zdPoint n‖ ^ 2) := by
            apply Real.sqrt_le_sqrt
            rw [hpi]
            exact Finset.single_le_sum (f := fun j => (n j : ℝ) ^ 2)
              (fun j _ => sq_nonneg (n j : ℝ)) (Finset.mem_univ i)
      _ = ‖zdPoint n‖ := Real.sqrt_sq (norm_nonneg _)
  have hN : |(n i : ℝ)| ≤ (e8ShellBound m : ℝ) := by
    calc |(n i : ℝ)|
        ≤ ‖zdPoint n‖ := habs
      _ ≤ C * ‖e8CoordCLE (zdPoint n)‖ := hop
      _ = B := by rw [hsqrt]
      _ ≤ (e8ShellBound m : ℝ) := by
            simpa [e8ShellBound, C, B] using Int.le_ceil B
  rw [Finset.mem_Icc]
  exact ⟨by exact_mod_cast (abs_le.mp hN).1, by exact_mod_cast (abs_le.mp hN).2⟩

/-- Die Normschale `{v ∈ E₈ | ‖v‖² = m}` ist endlich. -/
lemma e8_shell_finite (m : ℕ) :
    Finite {v : E8LatticeSet // e8NormSq (v : Fin 8 → ℝ) = (m : ℝ)} := by
  have hbij : Function.Bijective e8LatticeCoord :=
    ⟨e8LatticeCoord_injective, e8LatticeCoord_surjective⟩
  let e := Equiv.ofBijective e8LatticeCoord hbij
  have happly (x : Fin 8 → ℤ) : e x = e8LatticeCoord x :=
    Equiv.ofBijective_apply e8LatticeCoord hbij x
  let Box := Fin 8 → Finset.Icc (-e8ShellBound m) (e8ShellBound m)
  let toBox : {v : E8LatticeSet // e8NormSq (v : Fin 8 → ℝ) = (m : ℝ)} → Box :=
    fun v i => ⟨(e.symm v.1) i, by
      apply e8_shell_coord_mem m (e.symm v.1) _ i
      rw [← happly, e.apply_symm_apply]
      exact v.2⟩
  have hinj : Function.Injective toBox := by
    intro v w h
    apply Subtype.ext
    have hcoord : e.symm v.1 = e.symm w.1 := by
      ext i
      exact congrArg Subtype.val (congr_fun h i)
    rw [← e.apply_symm_apply v.1, ← e.apply_symm_apply w.1, hcoord]
  exact Finite.of_injective toBox hinj

lemma e8_fiber_finite (n : ℕ) :
    Finite {v : E8LatticeSet // e8HalfIndex v = n} := by
  have : Finite {v : E8LatticeSet //
      e8NormSq (v : Fin 8 → ℝ) = ((2 * n : ℕ) : ℝ)} :=
    e8_shell_finite (2 * n)
  exact Finite.of_equiv _
    (Equiv.subtypeEquivRight fun v => (e8HalfIndex_eq_iff v n).symm)

lemma e8_representationNumber_eq_fiber (n : ℕ) :
    e8_representationNumber (2 * n) =
      Nat.card {v : E8LatticeSet // e8HalfIndex v = n} :=
  Nat.card_congr (Equiv.subtypeEquivRight fun v => (e8HalfIndex_eq_iff v n).symm)

lemma e8Theta_summable (τ : ℍ) :
    Summable fun v : E8LatticeSet =>
      cexp (π * Complex.I * (τ : ℂ) * (e8NormSq (v : Fin 8 → ℝ) : ℂ)) :=
  Summable.of_norm_bounded (e8_lattice_gaussian_summable τ.im_pos) fun v =>
    le_of_eq (e8Theta_summand_norm v τ)

lemma e8Theta_term_eq_qPow (τ : ℍ) (v : E8LatticeSet) :
    cexp (π * Complex.I * (τ : ℂ) * (e8NormSq (v : Fin 8 → ℝ) : ℂ)) =
      Function.Periodic.qParam (1 : ℝ) (τ : ℂ) ^ e8HalfIndex v := by
  have h1 : ((1 : ℝ) : ℂ) = 1 := by norm_cast
  rw [Function.Periodic.qParam, h1, div_one, ← Complex.exp_nat_mul]
  congr 1
  rw [e8HalfIndex_spec]
  push_cast
  ring

lemma e8Theta_fiber_tsum (τ : ℍ) (n : ℕ) :
    ∑' v : {v : E8LatticeSet // e8HalfIndex v = n},
        cexp (π * Complex.I * (τ : ℂ) * (e8NormSq (v.1 : Fin 8 → ℝ) : ℂ)) =
      (e8_representationNumber (2 * n) : ℂ) •
        Function.Periodic.qParam (1 : ℝ) (τ : ℂ) ^ n := by
  have : Finite {v : E8LatticeSet // e8HalfIndex v = n} := e8_fiber_finite n
  let _ : Fintype {v : E8LatticeSet // e8HalfIndex v = n} := Fintype.ofFinite _
  rw [tsum_fintype]
  have hterm : ∀ v : {v : E8LatticeSet // e8HalfIndex v = n},
      cexp (π * Complex.I * (τ : ℂ) * (e8NormSq (v.1 : Fin 8 → ℝ) : ℂ)) =
        Function.Periodic.qParam (1 : ℝ) (τ : ℂ) ^ n := by
    intro v
    rw [e8Theta_term_eq_qPow, v.2]
  simp only [hterm, Finset.sum_const, Finset.card_univ]
  rw [← Nat.card_eq_fintype_card, e8_representationNumber_eq_fiber, nsmul_eq_mul,
    smul_eq_mul]

lemma e8Theta_hasSum_shells (τ : ℍ) :
    HasSum (fun n : ℕ => (e8_representationNumber (2 * n) : ℂ) •
        Function.Periodic.qParam (1 : ℝ) (τ : ℂ) ^ n) (e8Theta τ) := by
  let term : E8LatticeSet → ℂ := fun v =>
    cexp (π * Complex.I * (τ : ℂ) * (e8NormSq (v : Fin 8 → ℝ) : ℂ))
  have hs : HasSum term (e8Theta τ) := by
    simpa [e8Theta, term] using (e8Theta_summable τ).hasSum
  exact (hs.tsum_fiberwise e8HalfIndex).congr_fun fun n => (e8Theta_fiber_tsum τ n).symm

/-- `Θ(τ) = ∑_n r_{E8}(2n) q(τ)^n`. -/
theorem e8Theta_eq_sum_shells (τ : ℍ) :
    e8Theta τ = ∑' n : ℕ, (e8_representationNumber (2 * n) : ℂ) *
        Function.Periodic.qParam (1 : ℝ) (τ : ℂ) ^ n := by
  rw [← (e8Theta_hasSum_shells τ).tsum_eq]
  simp [smul_eq_mul]

/-- Koeffizient der `q`-Entwicklung: `a_n(Θ_{E8}) = r_{E8}(2n)`. -/
theorem e8Theta_qExpansion_coeff (n : ℕ) :
    (qExpansion 1 (e8Theta_modularForm : ℍ → ℂ)).coeff n =
      (e8_representationNumber (2 * n) : ℂ) := by
  have hf : ∀ τ : ℍ, HasSum
      (fun m : ℕ => (e8_representationNumber (2 * m) : ℂ) •
        Function.Periodic.qParam (1 : ℝ) (τ : ℂ) ^ m)
      (e8Theta_modularForm τ) := by
    intro τ
    rw [e8Theta_modularForm_apply]
    exact e8Theta_hasSum_shells τ
  exact (ModularFormClass.qExpansion_coeff_unique one_pos one_mem_strictPeriods_SL hf n).symm

/-- `r_{E8}(2n) = 240 σ₃(n)` für `n ≥ 1`. -/
theorem e8_representation_numbers (n : ℕ+) :
    e8_representationNumber (2 * (n : ℕ)) = 240 * σ 3 (n : ℕ) := by
  have hfun : (e8Theta_modularForm : ℍ → ℂ) = (E₄ : ℍ → ℂ) := by
    ext τ
    rw [e8Theta_modularForm_apply, e8_theta_eq_E4]
  have hcoeff :
      (qExpansion 1 (e8Theta_modularForm : ℍ → ℂ)).coeff (n : ℕ) =
        (qExpansion 1 (E₄ : ℍ → ℂ)).coeff (n : ℕ) :=
    congrArg (fun f : ℍ → ℂ => (qExpansion (1 : ℝ) f).coeff (n : ℕ)) hfun
  rw [e8Theta_qExpansion_coeff, e4_qExpansion_coeff_eq_sigma] at hcoeff
  exact_mod_cast hcoeff

end E8ArithmeticBridge
