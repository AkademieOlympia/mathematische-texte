import Riemann.E8ThetaQExpansion

open Complex
open scoped ArithmeticFunction.sigma ArithmeticFunction.zeta LSeries.notation

noncomputable section

namespace E8ArithmeticBridge

open ArithmeticFunction LSeries

/-!
Gate 3B.

* D.2 vergleicht die Koeffizienten für jedes `s : ℂ`.
* D.1 überträgt die absolute Konvergenz der Faltung `σ₃ = ζ * n³`
  für `Re s > 4` auf die Darstellungsreihe.
* Die Mellin-Identität `e8_mellin_product_identity` ist die Produktzerlegung.
-/

/-- Koeffizientensubstitution: `∑ r_{E8}(2n) n^{-s} = 240 ∑ σ₃(n) n^{-s}`. -/
theorem e8RepresentationDirichlet_eq_sigma (s : ℂ) :
    (∑' n : ℕ+, (e8_representationNumber (2 * (n : ℕ)) : ℂ) / (n : ℂ) ^ s) =
      (240 : ℂ) * ∑' n : ℕ+, ((σ 3 (n : ℕ) : ℕ) : ℂ) / (n : ℂ) ^ s := by
  rw [← tsum_mul_left]
  refine tsum_congr fun n => ?_
  rw [e8_representation_numbers n, Nat.cast_mul, Nat.cast_ofNat, mul_div_assoc]

/-- Absolute Konvergenz von `∑ σ₃(n) n^{-s}` für `Re s > 4`. -/
lemma sigma_three_dirichlet_summable {s : ℂ} (hs : 4 < s.re) :
    Summable (fun n : ℕ+ => ((σ 3 (n : ℕ) : ℕ) : ℂ) / (n : ℂ) ^ s) := by
  have hs1 : 1 < s.re := by linarith
  have hs3 : 1 < (s - 3).re := by
    rw [Complex.sub_re]
    norm_num
    linarith
  let σc : ArithmeticFunction ℂ := σ 3
  let ζc : ArithmeticFunction ℂ := ζ
  let pk : ArithmeticFunction ℂ := pow 3
  have hconv : σc = ζc * pk := by
    ext n
    rw [show σc = ((σ 3 : ArithmeticFunction ℕ) : ArithmeticFunction ℂ) from rfl,
      show ζc = ((ζ : ArithmeticFunction ℕ) : ArithmeticFunction ℂ) from rfl,
      show pk = ((pow 3 : ArithmeticFunction ℕ) : ArithmeticFunction ℂ) from rfl,
      ← zeta_mul_pow_eq_sigma]
    simp [natCoe_apply, mul_apply]
  have hterm_pow : ∀ n, term (↗pk) s n = term (1 : ℕ → ℂ) (s - 3) n := by
    intro n
    rcases eq_or_ne n 0 with rfl | hn
    · simp
    · rw [term_of_ne_zero hn, term_of_ne_zero hn]
      have hnC : (n : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hn
      simp only [pk, pow_apply, hn, and_false, ite_false, Pi.one_apply, natCoe_apply,
        Nat.cast_pow]
      have hpow3 : (n : ℂ) ^ (3 : ℕ) = (n : ℂ) ^ (3 : ℂ) :=
        (Complex.cpow_natCast (n : ℂ) 3).symm
      rw [hpow3, (Complex.cpow_sub (x := (n : ℂ)) (y := (3 : ℂ)) (z := s) hnC).symm,
        show (3 : ℂ) - s = -(s - 3) by ring, Complex.cpow_neg, inv_eq_one_div]
  have hsum_pow : LSeriesSummable (↗pk) s := by
    rw [LSeriesSummable, funext hterm_pow]
    exact LSeriesSummable_one_iff.mpr hs3
  have hsum_zeta : LSeriesSummable (↗ζc) s :=
    LSeriesSummable_zeta_iff.mpr hs1
  have hsum_sigma : LSeriesSummable (↗σc) s := by
    rw [hconv]
    exact ArithmeticFunction.LSeriesSummable_mul hsum_zeta hsum_pow
  have hterm : ∀ n : ℕ, term (↗σc) s n = ((σ 3 n : ℕ) : ℂ) / (n : ℂ) ^ s := by
    intro n
    rcases eq_or_ne n 0 with rfl | hn
    · simp [σc]
    · rw [term_of_ne_zero hn]
      simp [σc, natCoe_apply]
  have hnat : Summable (fun n : ℕ => ((σ 3 n : ℕ) : ℂ) / (n : ℂ) ^ s) := by
    rw [LSeriesSummable] at hsum_sigma
    exact (summable_congr hterm).mp hsum_sigma
  refine ((summable_pnat_iff_summable_nat
      (f := fun n : ℕ => ((σ 3 n : ℕ) : ℂ) / (n : ℂ) ^ s)).mpr hnat).congr
    fun n => ?_
  simp

/-- G3B-D.1: absolute Konvergenz von `∑ r_{E8}(2n) n^{-s}` für `Re s > 4`. -/
theorem e8RepresentationDirichlet_summable {s : ℂ} (hs : 4 < s.re) :
    Summable (fun n : ℕ+ =>
      (e8_representationNumber (2 * (n : ℕ)) : ℂ) / (n : ℂ) ^ s) := by
  refine ((sigma_three_dirichlet_summable hs).mul_left (240 : ℂ)).congr fun n => ?_
  rw [e8_representation_numbers n, Nat.cast_mul, Nat.cast_ofNat, mul_div_assoc]

/-- G3B-MAIN: `∑ r_{E8}(2n) n^{-s} = 240 ζ(s) ζ(s-3)` für `Re s > 4`. -/
theorem e8RepresentationDirichlet_eq_zeta {s : ℂ} (hs : 4 < s.re) :
    (∑' n : ℕ+, (e8_representationNumber (2 * (n : ℕ)) : ℂ) / (n : ℂ) ^ s) =
      240 * riemannZeta s * riemannZeta (s - 3) := by
  rw [e8RepresentationDirichlet_eq_sigma, e8_mellin_product_identity s hs]
  ring

end E8ArithmeticBridge
