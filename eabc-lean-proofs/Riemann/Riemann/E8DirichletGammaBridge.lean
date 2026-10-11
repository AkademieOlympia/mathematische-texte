import Mathlib.Analysis.SpecialFunctions.Gamma.Deligne
import Mathlib.Analysis.SpecialFunctions.Gamma.Deriv
import Riemann.E8DirichletFunctionalEquation

open Complex Filter
open scoped Topology

noncomputable section

namespace E8ArithmeticBridge

/-!
Gate 3D.

`G(s) = Γ_ℝ(s) Γ_ℝ(s-3)` verbindet `Λ_{E_8}` mit `D_{E_8}`.
Die Brücke `Λ_{E_8}(s) = G(s) D_{E_8}(s)` braucht `G(s) ≠ 0`:
`riemannZeta_def_of_ne_zero` dividiert durch `Γ_ℝ`, und Mathlibs
totale Gamma-Funktion hat an ihren Polen den Wert `0`.
Dort ist das Produkt `Γ_ℝ(s) ζ(s)` der Wert `0`, nicht `Λ(s)`.

Dieselbe Bedingung an `s` und an `4 - s` trägt die multiplikative
Funktionalgleichung. Ohne sie ist die Gleichung falsch, siehe
`e8GammaFactor_mul_dirichlet_not_symmetric`.
-/

/-- \(G(s) = \Gamma_{\mathbb{R}}(s)\Gamma_{\mathbb{R}}(s-3)\). -/
def e8GammaFactor (s : ℂ) : ℂ :=
  Gammaℝ s * Gammaℝ (s - 3)

lemma e8GammaFactor_eq_zero_iff (s : ℂ) :
    e8GammaFactor s = 0 ↔
      (∃ n : ℕ, s = -(2 * n)) ∨ ∃ n : ℕ, s - 3 = -(2 * n) := by
  rw [e8GammaFactor, mul_eq_zero, Gammaℝ_eq_zero_iff, Gammaℝ_eq_zero_iff]

lemma not_gammaR_zero_and_shift (s : ℂ) :
    ¬ (Gammaℝ s = 0 ∧ Gammaℝ (s - 3) = 0) := by
  intro ⟨hs, hs3⟩
  rw [Gammaℝ_eq_zero_iff] at hs hs3
  obtain ⟨n, rfl⟩ := hs
  obtain ⟨m, hm⟩ := hs3
  have hsum : (2 : ℂ) * (m : ℂ) = (2 : ℂ) * n + 3 := by
    refine neg_injective ?_
    rw [neg_add]
    exact hm.symm
  have hmulC : (2 : ℂ) * ((m : ℂ) - n) = 3 := by
    rw [mul_sub, hsum]
    ring
  have hmul : (2 : ℤ) * ((m : ℤ) - n) = 3 := by
    exact_mod_cast hmulC
  exact (by decide : ¬ Even (3 : ℤ)) ⟨(m : ℤ) - n, by rw [← two_mul]; exact hmul.symm⟩

lemma continuousAt_Gammaℝ {s : ℂ} (hs : Gammaℝ s ≠ 0) : ContinuousAt Gammaℝ s := by
  have hs' : ∀ n : ℕ, s ≠ -(2 * n) := fun n hn =>
    hs (Gammaℝ_eq_zero_iff.mpr ⟨n, hn⟩)
  unfold Gammaℝ
  refine ContinuousAt.mul ?_ ?_
  · have hexp : ContinuousAt (fun z : ℂ => -z / 2) s :=
      (continuous_id.neg.div_const (2 : ℂ)).continuousAt
    exact ContinuousAt.comp (f := fun z : ℂ => -z / 2) (x := s)
      (continuousAt_const_cpow (b := -s / 2) (ofReal_ne_zero.mpr Real.pi_ne_zero)) hexp
  · have hdiv : ContinuousAt (fun z : ℂ => z / 2) s :=
      (continuous_id.div_const (2 : ℂ)).continuousAt
    refine ContinuousAt.comp (f := fun z : ℂ => z / 2) (x := s)
      (continuousAt_Gamma (s / 2) ?_) hdiv
    intro m hm
    apply hs' m
    have hsdiv : s = 2 * (s / 2) := by field_simp
    rw [hsdiv, hm]
    ring

lemma not_continuousAt_Gammaℝ {s : ℂ} (hs : Gammaℝ s = 0) : ¬ ContinuousAt Gammaℝ s := by
  rw [Gammaℝ_eq_zero_iff] at hs
  obtain ⟨n, rfl⟩ := hs
  intro hcont
  let s0 : ℂ := -(2 * (n : ℂ))
  have hπ : ContinuousAt (fun z : ℂ => (Real.pi : ℂ) ^ (-z / 2)) s0 := by
    have hexp : ContinuousAt (fun z : ℂ => -z / 2) s0 :=
      (continuous_id.neg.div_const (2 : ℂ)).continuousAt
    exact ContinuousAt.comp (f := fun z : ℂ => -z / 2) (x := s0)
      (continuousAt_const_cpow (b := -s0 / 2) (ofReal_ne_zero.mpr Real.pi_ne_zero)) hexp
  have hπne : (Real.pi : ℂ) ^ (-s0 / 2) ≠ 0 :=
    cpow_ne_zero_iff.mpr (Or.inl (ofReal_ne_zero.mpr Real.pi_ne_zero))
  have hquot : ContinuousAt (fun z : ℂ => Gammaℝ z / ((Real.pi : ℂ) ^ (-z / 2))) s0 :=
    hcont.div hπ hπne
  have hhalf : ContinuousAt (fun z : ℂ => Gamma (z / 2)) s0 := by
    refine hquot.congr (EventuallyEq.of_eq ?_)
    ext z
    rw [Gammaℝ_def]
    exact mul_div_cancel_left₀ _
      (cpow_ne_zero_iff.mpr (Or.inl (ofReal_ne_zero.mpr Real.pi_ne_zero)))
  have hΓ : ContinuousAt Gamma (-(n : ℂ)) := by
    have hmul : ContinuousAt (fun t : ℂ => (2 : ℂ) * t) (-(n : ℂ)) :=
      (continuous_const.mul continuous_id).continuousAt
    have hhalf' : ContinuousAt (fun z : ℂ => Gamma (z / 2)) (2 * -(n : ℂ)) := by
      simpa [s0, neg_mul, mul_comm] using hhalf
    refine (ContinuousAt.comp (f := fun t : ℂ => (2 : ℂ) * t) (x := -(n : ℂ)) hhalf' hmul).congr
      (EventuallyEq.of_eq ?_)
    ext t
    simp only [Function.comp_apply]
    field_simp
  exact not_continuousAt_Gamma_neg_nat n hΓ

/-- Jede Nullstelle der totalen Funktion `G` ist eine Unstetigkeitsstelle. -/
theorem not_continuousAt_e8GammaFactor {s : ℂ} (hs : e8GammaFactor s = 0) :
    ¬ ContinuousAt e8GammaFactor s := by
  rw [e8GammaFactor, mul_eq_zero] at hs
  rcases hs with hs | hs
  · have hother : Gammaℝ (s - 3) ≠ 0 := by
      intro hother
      exact not_gammaR_zero_and_shift s ⟨hs, hother⟩
    intro hcont
    have hshift : ContinuousAt (fun z : ℂ => Gammaℝ (z - 3)) s := by
      have hsub : ContinuousAt (fun z : ℂ => z - 3) s :=
        (continuous_sub_right (3 : ℂ)).continuousAt
      exact ContinuousAt.comp (f := fun z : ℂ => z - 3) (x := s)
        (continuousAt_Gammaℝ (s := s - 3) hother) hsub
    have hΓ : ContinuousAt Gammaℝ s := by
      refine (hcont.div hshift hother).congr ?_
      filter_upwards [(hshift.eventually_ne hother)] with z hz
      simp only [e8GammaFactor, Pi.div_apply]
      exact mul_div_cancel_right₀ _ hz
    exact not_continuousAt_Gammaℝ hs hΓ
  · have hother : Gammaℝ s ≠ 0 := by
      intro hother
      exact not_gammaR_zero_and_shift s ⟨hother, hs⟩
    intro hcont
    have hshift : ContinuousAt (fun z : ℂ => Gammaℝ (z - 3)) s := by
      refine (hcont.div (continuousAt_Gammaℝ hother) hother).congr ?_
      filter_upwards [(continuousAt_Gammaℝ hother).eventually_ne hother] with z hz
      simp only [e8GammaFactor, Pi.div_apply]
      exact mul_div_cancel_left₀ _ hz
    have hcomp : ContinuousAt Gammaℝ (s - 3) := by
      have hadd : ContinuousAt (fun z : ℂ => z + 3) (s - 3) :=
        (continuous_id.add continuous_const).continuousAt
      have hshift' : ContinuousAt (fun z : ℂ => Gammaℝ (z - 3)) ((s - 3) + 3) := by
        simpa [sub_add_cancel] using hshift
      refine (ContinuousAt.comp (f := fun z : ℂ => z + 3) (x := s - 3) hshift' hadd).congr ?_
      filter_upwards with z
      simp [Function.comp_apply]
    exact not_continuousAt_Gammaℝ hs hcomp

/-- G3D-G.1: \(\Lambda_{E_8}(s) = G(s)\, D_{E_8}(s)\), sobald \(G(s) \neq 0\). -/
theorem e8CompletedDirichlet_eq_gamma_mul {s : ℂ} (hG : e8GammaFactor s ≠ 0) :
    e8CompletedDirichlet s = e8GammaFactor s * e8MeromorphicDirichlet s := by
  have hGs : Gammaℝ s ≠ 0 := fun h => hG <| by simp [e8GammaFactor, h]
  have hG3 : Gammaℝ (s - 3) ≠ 0 := fun h => hG <| by simp [e8GammaFactor, h]
  have hs0 : s ≠ 0 := by
    intro rfl
    exact hGs (Gammaℝ_eq_zero_iff.mpr ⟨0, by simp⟩)
  have hs3 : s ≠ 3 := by
    intro hs3
    apply hG3
    rw [hs3]
    exact Gammaℝ_eq_zero_iff.mpr ⟨0, by simp⟩
  have hζ : Gammaℝ s * riemannZeta s = completedRiemannZeta s := by
    rw [riemannZeta_def_of_ne_zero hs0]
    exact mul_div_cancel₀ _ hGs
  have hζ3 : Gammaℝ (s - 3) * riemannZeta (s - 3) = completedRiemannZeta (s - 3) := by
    rw [riemannZeta_def_of_ne_zero (sub_ne_zero.mpr hs3)]
    exact mul_div_cancel₀ _ hG3
  unfold e8CompletedDirichlet e8GammaFactor e8MeromorphicDirichlet
  rw [← hζ, ← hζ3]
  ring

/-- G3D-G.2: Quotient, ebenfalls unter \(G(s) \neq 0\). -/
theorem e8MeromorphicDirichlet_eq_completed_div {s : ℂ} (hG : e8GammaFactor s ≠ 0) :
    e8MeromorphicDirichlet s = e8CompletedDirichlet s / e8GammaFactor s := by
  refine (eq_div_iff hG).mpr ?_
  rw [mul_comm]
  exact (e8CompletedDirichlet_eq_gamma_mul hG).symm

/-- G3D-FE.U.1: multiplikative Funktionalgleichung, ohne Division. -/
theorem e8MeromorphicDirichlet_functional_equation_mul {s : ℂ}
    (hG : e8GammaFactor s ≠ 0) (hG4 : e8GammaFactor (4 - s) ≠ 0) :
    e8GammaFactor s * e8MeromorphicDirichlet s =
      e8GammaFactor (4 - s) * e8MeromorphicDirichlet (4 - s) := by
  rw [← e8CompletedDirichlet_eq_gamma_mul hG, ← e8CompletedDirichlet_eq_gamma_mul hG4]
  exact (e8CompletedDirichlet_functional_equation s).symm

/-- G3D-FE.U.2: \(D(4-s) = \frac{G(s)}{G(4-s)} D(s)\). -/
theorem e8MeromorphicDirichlet_functional_equation_div {s : ℂ}
    (hG : e8GammaFactor s ≠ 0) (hG4 : e8GammaFactor (4 - s) ≠ 0) :
    e8MeromorphicDirichlet (4 - s) =
      (e8GammaFactor s / e8GammaFactor (4 - s)) * e8MeromorphicDirichlet s := by
  refine mul_left_cancel₀ hG4 ?_
  rw [← e8MeromorphicDirichlet_functional_equation_mul hG hG4]
  field_simp [hG4]

/-- Ohne die Nichtverschwindensbedingung ist die multiplikative Symmetrie falsch. -/
theorem e8GammaFactor_mul_dirichlet_not_symmetric :
    ¬ ∀ s : ℂ, e8GammaFactor s * e8MeromorphicDirichlet s =
      e8GammaFactor (4 - s) * e8MeromorphicDirichlet (4 - s) := by
  intro h
  have h0 : e8GammaFactor (-2) * e8MeromorphicDirichlet (-2) = 0 := by
    have hΓ : Gammaℝ (-2) = 0 := Gammaℝ_eq_zero_iff.mpr ⟨1, by norm_num⟩
    simp [e8GammaFactor, hΓ]
  have h6 : e8GammaFactor (4 - (-2)) * e8MeromorphicDirichlet (4 - (-2)) ≠ 0 := by
    have h62 : (4 : ℂ) - (-2) = 6 := by norm_num
    rw [h62]
    have hG : e8GammaFactor 6 ≠ 0 := by
      refine mul_ne_zero (Gammaℝ_ne_zero_of_re_pos (by norm_num))
        (Gammaℝ_ne_zero_of_re_pos (by norm_num))
    have hD : e8MeromorphicDirichlet 6 ≠ 0 := by
      refine mul_ne_zero (mul_ne_zero (by norm_num)
        (riemannZeta_ne_zero_of_one_lt_re (by norm_num)))
        (riemannZeta_ne_zero_of_one_lt_re (by norm_num))
    exact mul_ne_zero hG hD
  have hsym := h (-2)
  rw [h0] at hsym
  exact h6 hsym.symm

end E8ArithmeticBridge
