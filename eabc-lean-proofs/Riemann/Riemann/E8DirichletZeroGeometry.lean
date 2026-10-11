import Mathlib.NumberTheory.LSeries.Nonvanishing
import Mathlib.NumberTheory.LSeries.RiemannZeta
import Riemann.E8DirichletGammaBridge

open Complex

noncomputable section

namespace E8ArithmeticBridge

/-!
Gate 3E.

`e8MeromorphicDirichlet` ist die totalisierte Funktion
`240 * riemannZeta s * riemannZeta (s - 3)`. Ihre Nullstellen
sind deshalb genau die Nullstellen der beiden Zeta-Faktoren.
Pole und hebbare Singularitäten bleiben in Gate 3C: der Wert
bei `s = 1` ist `0`, weil `ζ(-2) = 0`, während der punktierte
Grenzwert dort `240 * deriv riemannZeta (-2)` ist.

`E8NontrivialZero` spricht nur über die offenen Streifen
`0 < Re < 1` und `3 < Re < 4`. Daraus folgt `Re s ≠ 2`, ohne
die Riemannsche Vermutung. Die Äquivalenz mit Mathlibs
`RiemannHypothesis` ist eine bedingte Übersetzung, kein Beweis
der Vermutung.
-/

/-- Algebraische Nullstellen der totalisierten Funktion. -/
theorem e8MeromorphicDirichlet_eq_zero_iff (s : ℂ) :
    e8MeromorphicDirichlet s = 0 ↔
      riemannZeta s = 0 ∨ riemannZeta (s - 3) = 0 := by
  unfold e8MeromorphicDirichlet
  constructor
  · intro h
    rcases mul_eq_zero.mp h with hprod | hs3
    · rcases mul_eq_zero.mp hprod with h240 | hs
      · norm_num at h240
      · exact Or.inl hs
    · exact Or.inr hs3
  · rintro (hs | hs3)
    · simp [hs]
    · simp [hs3]

/-- Bei `s = 1` ist der totalisierte Wert `0`, wegen `ζ(-2) = 0`. -/
theorem e8MeromorphicDirichlet_one : e8MeromorphicDirichlet 1 = 0 := by
  rw [e8MeromorphicDirichlet_eq_zero_iff]
  right
  have h13 : (1 - 3 : ℂ) = -2 := by norm_num
  simpa [h13] using riemannZeta_neg_two_mul_nat_add_one 0

theorem e8_zero_embedding_left {ρ : ℂ} (h : riemannZeta ρ = 0) :
    e8MeromorphicDirichlet ρ = 0 := by
  rw [e8MeromorphicDirichlet_eq_zero_iff]
  exact Or.inl h

theorem e8_zero_embedding_right {ρ : ℂ} (h : riemannZeta ρ = 0) :
    e8MeromorphicDirichlet (3 + ρ) = 0 := by
  rw [e8MeromorphicDirichlet_eq_zero_iff]
  refine Or.inr ?_
  rw [show (3 + ρ - 3 : ℂ) = ρ by ring]
  exact h

lemma e8_reflection_identity (ρ : ℂ) : (4 : ℂ) - ρ = 3 + (1 - ρ) := by
  ring

/-- Nullstelle von `ζ` im Streifen `0 < Re < 1`, oder von `ζ(· - 3)` dort. -/
def E8NontrivialZero (s : ℂ) : Prop :=
  (riemannZeta s = 0 ∧ 0 < s.re ∧ s.re < 1) ∨
    (riemannZeta (s - 3) = 0 ∧ 0 < (s - 3).re ∧ (s - 3).re < 1)

lemma re_sub_three (s : ℂ) : (s - 3).re = s.re - 3 := by
  rw [Complex.sub_re]
  norm_num

/-- Nichttriviale E8-Nullstellen liegen nicht auf der Symmetrieachse `Re s = 2`. -/
theorem e8_nontrivial_zero_ne_two {s : ℂ} (hs : E8NontrivialZero s) : s.re ≠ 2 := by
  rcases hs with ⟨_, hlo, hhi⟩ | ⟨_, hlo, hhi⟩
  · intro h
    linarith
  · rw [re_sub_three] at hlo hhi
    intro h
    linarith

lemma not_trivial_zero_of_re_pos {s : ℂ} (hlo : 0 < s.re) :
    ¬ ∃ n : ℕ, s = -2 * (n + 1) := by
  intro ⟨n, hn⟩
  have hre : s.re = -2 * ((n : ℝ) + 1) := by
    rw [hn]
    simp only [mul_re, neg_re, neg_im, add_re, add_im, natCast_re, natCast_im, re_ofNat,
      im_ofNat, one_re, one_im]
    ring
  have hpos : (0 : ℝ) < (n : ℝ) + 1 := by positivity
  have hneg : s.re < 0 := by
    rw [hre]
    nlinarith
  linarith

lemma re_eq_half_of_strip (h_rh : RiemannHypothesis) {s : ℂ}
    (hz : riemannZeta s = 0) (hlo : 0 < s.re) (hhi : s.re < 1) : s.re = 1 / 2 := by
  have hs1 : s ≠ 1 := by
    intro h
    rw [h, Complex.one_re] at hhi
    linarith
  exact h_rh s hz (not_trivial_zero_of_re_pos hlo) hs1

/-- Mathlibs Riemannsche Vermutung ist die Aussage, dass jede
nichttriviale E8-Nullstelle auf `Re = 1/2` oder `Re = 7/2` liegt. -/
theorem riemann_hypothesis_iff_e8_zero_line :
    RiemannHypothesis ↔
      ∀ s, E8NontrivialZero s → s.re = 1 / 2 ∨ s.re = 7 / 2 := by
  constructor
  · intro h_rh s hs
    rcases hs with ⟨hz, hlo, hhi⟩ | ⟨hz, hlo, hhi⟩
    · exact Or.inl (re_eq_half_of_strip h_rh hz hlo hhi)
    · right
      have hhalf : (s - 3).re = 1 / 2 := re_eq_half_of_strip h_rh hz hlo hhi
      rw [re_sub_three] at hhalf
      linarith
  · intro hE8 s hz htriv hs1
    by_cases hpos : 0 < s.re
    · by_cases hlt : s.re < 1
      · rcases hE8 s (Or.inl ⟨hz, hpos, hlt⟩) with hhalf | hseven
        · exact hhalf
        · linarith
      · have hge : 1 ≤ s.re := le_of_not_gt hlt
        exact absurd hz (riemannZeta_ne_zero_of_one_le_re hge)
    · have hle : s.re ≤ 0 := le_of_not_gt hpos
      have hΓ : Gammaℝ s ≠ 0 := by
        intro hΓ
        rw [Gammaℝ_eq_zero_iff] at hΓ
        obtain ⟨n, rfl⟩ := hΓ
        cases n with
        | zero =>
          simp [riemannZeta_zero] at hz
        | succ k =>
          apply htriv
          refine ⟨k, ?_⟩
          simp only [Nat.cast_add, Nat.cast_one]
          ring
      have hs0 : s ≠ 0 := by
        intro h0
        exact hΓ (by simpa [h0] using (Gammaℝ_eq_zero_iff.mpr ⟨0, by simp⟩))
      have hΛ : completedRiemannZeta s = 0 := by
        rw [← mul_div_cancel₀ (completedRiemannZeta s) hΓ,
          ← riemannZeta_def_of_ne_zero hs0, hz, mul_zero]
      have h1s : 1 - s ≠ 0 := by
        intro h
        apply hs1
        linear_combination -h
      have hΓ1 : Gammaℝ (1 - s) ≠ 0 :=
        Gammaℝ_ne_zero_of_re_pos <| by
          rw [Complex.sub_re, Complex.one_re]
          linarith
      have hζ1 : riemannZeta (1 - s) = 0 := by
        rw [riemannZeta_def_of_ne_zero h1s, completedRiemannZeta_one_sub, hΛ, zero_div]
      have hre1 : 1 ≤ (1 - s).re := by
        rw [Complex.sub_re, Complex.one_re]
        linarith
      exact absurd hζ1 (riemannZeta_ne_zero_of_one_le_re hre1)

end E8ArithmeticBridge
