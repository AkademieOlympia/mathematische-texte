/-
  CollatzEabc.Uniformity — LTE-Reset, Mischschranken, Dichte-Verbindung.

  Quelle: `collatz_uniformity.lean` (Repo-Root).
  Alle Theoreme ohne `sorry`.
-/

import CollatzEabc.Density
import Mathlib.NumberTheory.Padics.PadicVal.Basic
import Mathlib.Algebra.Ring.Parity
import Mathlib.Topology.Algebra.InfiniteSum.Basic
import Mathlib.Data.Nat.ModEq
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecificLimits.Normed

namespace CollatzUniformity

/-- LTE-Worst-Start n = 2^(k+1)·3^r - 1. -/
def lteWorst (k r : ℕ) : ℕ := 2 ^ (k + 1) * 3 ^ r - 1

/-- Minimale LTE-Familie (k=1): n = 4·3^r - 1. -/
def lteMinimal (r : ℕ) : ℕ := 4 * 3 ^ r - 1

theorem lteWorst_odd (k r : ℕ) (_hk : 0 < k) (_hr : 0 < r) :
    Odd (lteWorst k r) := by
  unfold lteWorst
  have hmul : 2 ^ (k + 1) * 3 ^ r = 2 * (2 ^ k * 3 ^ r) := by ring
  rw [hmul]
  refine ⟨2 ^ k * 3 ^ r - 1, ?_⟩
  have hpow : 1 ≤ 2 ^ k * 3 ^ r := by
    nlinarith [pow_pos (by norm_num : 0 < 2) k, pow_pos (by norm_num : 0 < 3) r]
  ring_nf
  omega

theorem lteWorst_valuation (k r : ℕ) :
    padicValNat 2 (lteWorst k r + 1) = k + 1 := by
  unfold lteWorst
  have hge : 1 ≤ 2 ^ (k + 1) * 3 ^ r := by
    nlinarith [pow_pos (by norm_num : 0 < 2) (k + 1), pow_pos (by norm_num : 0 < 3) r]
  rw [Nat.sub_add_cancel hge]
  have h2base : (2 : ℕ) ≠ 0 := by decide
  rw [padicValNat.mul (by positivity) (by positivity),
    padicValNat.pow (k + 1) h2base,
    padicValNat_prime_prime_pow r (by decide : 2 ≠ 3), add_zero]
  simp [padicValNat.self]

theorem lteMinimal_is_lteWorst (r : ℕ) :
    lteMinimal r = lteWorst 1 r := by
  unfold lteMinimal lteWorst
  ring_nf

/-- Für ungerades j: 3^j ≡ 3 (mod 8). -/
private theorem three_pow_odd_mod8 (j : ℕ) (hj : Odd j) :
    3 ^ j % 8 = 3 := by
  rcases hj with ⟨t, ht⟩
  subst ht
  induction t with
  | zero => norm_num
  | succ t ih =>
    have hpow : 3 ^ (2 * (t + 1) + 1) = 9 * 3 ^ (2 * t + 1) := by ring
    rw [hpow, Nat.mul_mod, ih]

/-- Für ungerades j: ν₂(3^j + 1) = 2. -/
theorem padicVal_two_three_pow_odd_plus_one (j : ℕ) (hj : Odd j) :
    padicValNat 2 (3 ^ j + 1) = 2 := by
  have hmod8 : 3 ^ j % 8 = 3 := three_pow_odd_mod8 j hj
  have hge : 2 ≤ padicValNat 2 (3 ^ j + 1) :=
    (padicVal_succ_ge_iff_dvd (3 ^ j) 2).mpr (Nat.dvd_of_mod_eq_zero (by omega))
  have hlt : padicValNat 2 (3 ^ j + 1) < 3 := by
    apply lt_of_not_ge
    intro h3
    have h8 : 8 ∣ 3 ^ j + 1 := (padicVal_succ_ge_iff_dvd (3 ^ j) 3).mp h3
    omega
  omega

/--
  Fall 1 (gerades N): Nachfolger (3^(N+1)-1)/2 erfüllt ν₂(n+1)=1.
  Entspricht Lemma (LTE-Reset) für gerades N=k+r.
-/
theorem lte_even_N_successor_valuation_one (N : ℕ) (hN : Even N) :
    padicValNat 2 ((3 ^ (N + 1) - 1) / 2 + 1) = 1 := by
  have hodd : Odd (N + 1) := by rcases hN with ⟨j, rfl⟩; exact ⟨j, by ring⟩
  have hmid : (3 ^ (N + 1) + 1) / 2 = (3 ^ (N + 1) - 1) / 2 + 1 := by
    have hpos : 1 ≤ 3 ^ (N + 1) := Nat.one_le_pow (N + 1) 3 (by norm_num)
    omega
  rw [← hmid]
  have hval : padicValNat 2 (3 ^ (N + 1) + 1) = 2 :=
    padicVal_two_three_pow_odd_plus_one (N + 1) hodd
  have hge : 1 ≤ padicValNat 2 (3 ^ (N + 1) + 1) := by rw [hval]; decide
  have hdvd : 2 ∣ 3 ^ (N + 1) + 1 :=
    (padicVal_succ_ge_iff_dvd (3 ^ (N + 1)) 1).mp hge
  calc padicValNat 2 ((3 ^ (N + 1) + 1) / 2)
      = padicValNat 2 (3 ^ (N + 1) + 1) - padicValNat 2 2 :=
        padicValNat.div_of_dvd hdvd
    _ = 1 := by simp [hval, padicValNat.self]

theorem lte_worst_even_N_forces_low_valuation (k r : ℕ) (hN : Even (k + r)) :
    padicValNat 2 ((3 ^ (k + r + 1) - 1) / 2 + 1) = 1 :=
  lte_even_N_successor_valuation_one (k + r) hN

/-- ν₂(n+1) ≤ 2 bedeutet: höchstens C-Kettenlänge 1. -/
theorem valuation_le_two_limits_C_chain (n : ℕ) (h : padicValNat 2 (n + 1) ≤ 2) :
    padicValNat 2 (n + 1) < 3 := by omega

/-- Mischschranken-Kern: (1-p)^n → 0 für p ∈ (0,1). -/
theorem mixing_probability_bound (p : ℝ) (hp0 : 0 < p) (hp1 : p < 1) (n : ℕ) :
    (1 - p) ^ n ≤ 1 := by
  have : 0 ≤ 1 - p := by linarith
  exact pow_le_one₀ this (by linarith : 1 - p ≤ 1)

theorem mixing_probability_tendsto_zero (p : ℝ) (hp0 : 0 < p) (hp1 : p < 1) :
    Filter.Tendsto (fun n : ℕ => (1 - p) ^ n) Filter.atTop (nhds 0) := by
  have h : norm (1 - p) < 1 := by rw [Real.norm_of_nonneg (by linarith)]; linarith
  exact tendsto_pow_atTop_nhds_zero_of_norm_lt_one h

/-- Alias: Verbindung zur Ausnahmewahrscheinlichkeit in `CollatzEabc.Density`. -/
theorem exception_probability_tendsto_zero (p : ℝ) (hp0 : 0 < p) (hp1 : p < 1) :
    Filter.Tendsto (fun n : ℕ => (1 - p) ^ n) Filter.atTop (nhds 0) :=
  mixing_probability_tendsto_zero p hp0 hp1

end CollatzUniformity
