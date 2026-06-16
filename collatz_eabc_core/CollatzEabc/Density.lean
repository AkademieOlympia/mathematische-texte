/-
  CollatzEabc.Density — Dichte-Lemmata für C-Ketten in der Collatz-Dynamik.

  Quelle: `collatz_density_appendix.lean` (Repo-Root).
  Alle Theoreme ohne `sorry`. An Mathlib 4.29 angepasst.
-/

import Mathlib.NumberTheory.Padics.PadicVal.Basic
import Mathlib.Algebra.Ring.Parity
import Mathlib.Topology.Algebra.InfiniteSum.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Data.Nat.Log
import Mathlib.Analysis.SpecificLimits.Normed
import Mathlib.Algebra.Order.Ring.Basic

theorem padicVal_ge_iff_dvd (m k : ℕ) (hm : m ≠ 0) :
    padicValNat 2 m ≥ k ↔ 2 ^ k ∣ m := by
  rw [ge_iff_le, ← padicValNat_dvd_iff_le hm]

theorem padicVal_succ_ge_iff_dvd (n k : ℕ) :
    padicValNat 2 (n + 1) ≥ k ↔ 2 ^ k ∣ (n + 1) := by
  apply padicVal_ge_iff_dvd
  exact Nat.succ_ne_zero n

private theorem two_pow_succ_mul_sub_one (k m : ℕ) :
    2 ^ (k + 1) * m + (2 ^ (k + 1) - 1) + 1 = 2 ^ (k + 1) * (m + 1) := by
  have hle : 1 ≤ 2 ^ (k + 1) := Nat.one_le_pow (k + 1) 2 (by decide)
  rw [add_assoc, Nat.sub_add_cancel hle, mul_add, mul_one]

private theorem two_pow_succ_mul_odd_form (k m : ℕ) :
    2 ^ (k + 1) * m + (2 ^ (k + 1) - 1) = 2 * (2 ^ k * m + (2 ^ k - 1)) + 1 := by
  induction k with
  | zero => simp
  | succ k ih =>
    rw [pow_succ, pow_succ]
    ring_nf at ih ⊢
    omega

private theorem two_pow_succ_mul_pred_add_sub_one (k j : ℕ) (hj : 0 < j) :
    2 ^ (k + 1) * j - 1 = 2 ^ (k + 1) * j.pred + (2 ^ (k + 1) - 1) := by
  have hjs : j.pred + 1 = j := Nat.succ_pred (Nat.ne_of_gt hj)
  have hp : 1 ≤ 2 ^ (k + 1) := Nat.one_le_pow (k + 1) 2 (by decide)
  conv_lhs => rw [← hjs, mul_add, mul_one]
  rw [Nat.add_sub_assoc hp]

theorem C_chain_start_class (k : ℕ) :
    ∀ n : ℕ, Odd n ∧ padicValNat 2 (n + 1) ≥ k + 1 ↔
      ∃ m : ℕ, n = 2 ^ (k + 1) * m + (2 ^ (k + 1) - 1) := by
  intro n
  constructor
  · intro ⟨_, hval⟩
    rw [padicVal_succ_ge_iff_dvd] at hval
    obtain ⟨j, hj⟩ := hval
    have hjpos : 0 < j := by
      by_contra hle
      push Not at hle
      have hj0 : j = 0 := by omega
      rw [hj0, mul_zero] at hj
      omega
    refine ⟨j.pred, ?_⟩
    have hn : n = 2 ^ (k + 1) * j - 1 := by omega
    rw [hn, two_pow_succ_mul_pred_add_sub_one k j hjpos]
  · intro ⟨m, hm⟩
    constructor
    · rw [hm, two_pow_succ_mul_odd_form k m]
      simpa using odd_two_mul_add_one (2 ^ k * m + (2 ^ k - 1))
    · rw [padicVal_succ_ge_iff_dvd, hm]
      refine ⟨m + 1, ?_⟩
      exact two_pow_succ_mul_sub_one k m

theorem density_C_chains_finite (N k : ℕ) :
    (Finset.filter (fun n => padicValNat 2 (n + 1) ≥ k + 1)
      (Finset.filter Odd (Finset.range (2 * N + 1)))).card
    = N / 2 ^ k := by
  have hset_eq : Finset.filter (fun n => padicValNat 2 (n + 1) ≥ k + 1)
      (Finset.filter Odd (Finset.range (2 * N + 1))) =
      (Finset.range (2 * N + 1)).filter (fun e => 2 ^ (k + 1) ∣ e + 1) := by
    ext n; simp only [Finset.mem_filter, Finset.mem_range]
    constructor
    · intro h
      exact ⟨h.1.1, (padicVal_succ_ge_iff_dvd n (k + 1)).mp h.2⟩
    · intro h
      have hodd : Odd n := by
        rcases Nat.even_or_odd n with heven | hod
        · rcases heven with ⟨t, rfl⟩
          have hdvd' : 2 ^ (k + 1) ∣ 2 * t + 1 := by simpa [two_mul] using h.2
          have h2 : 2 ∣ 2 * t + 1 :=
            dvd_trans (dvd_pow_self 2 (Nat.succ_ne_zero k)) hdvd'
          norm_num at h2
        · exact hod
      exact And.intro (And.intro h.1 hodd) ((padicVal_succ_ge_iff_dvd n (k + 1)).mpr h.2)
  rw [hset_eq]
  have h_count := Nat.card_multiples (2 * N + 1) (2 ^ (k + 1))
  rw [h_count, show (2 : ℕ) ^ (k + 1) = 2 * 2 ^ k from by ring, ← Nat.div_div_eq_div_mul]
  have h_div : (2 * N + 1) / 2 = N := by omega
  rw [h_div]

theorem density_C_chains_bound (k : ℕ) :
    ∀ N : ℕ, 0 < N →
    (Finset.filter (fun n => padicValNat 2 (n + 1) ≥ k + 1)
      (Finset.filter Odd (Finset.range (2 * N + 1)))).card
    ≤ N / 2 ^ k + 1 := by
  intro N _; rw [density_C_chains_finite]; exact Nat.le_add_right _ _

theorem collatz_norm_identity (n : ℤ) : (3 * n + 1) ^ 2 = 9 * n ^ 2 + 6 * n + 1 := by ring

theorem exception_probability_decay (p : ℝ) (_hp0 : 0 < p) (hp1 : p < 1) (n : ℕ) :
    (1 - p) ^ n ≤ (1 - p) ^ 0 := by
  rw [pow_zero]
  exact pow_le_one₀ (by linarith : 0 ≤ 1 - p) (by linarith : 1 - p ≤ 1)

theorem exception_probability_tendsto_zero (p : ℝ) (hp0 : 0 < p) (hp1 : p < 1) :
    Filter.Tendsto (fun n : ℕ => (1 - p) ^ n) Filter.atTop (nhds 0) := by
  have h : norm (1 - p) < 1 := by rw [Real.norm_of_nonneg (by linarith)]; linarith
  exact tendsto_pow_atTop_nhds_zero_of_norm_lt_one h
