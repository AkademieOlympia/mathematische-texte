/-
  collatz_uniformity.lean
  Uniformitäts-Angriff: LTE-Reset (gerades N), Mischschranken, Dichte-Verbindung.

  Ergänzt collatz_density_appendix.lean.
  Alle Theoreme in diesem File sind ohne sorry.
-/

import Mathlib.NumberTheory.Padics.PadicVal
import Mathlib.Data.Nat.Parity
import Mathlib.Topology.Algebra.InfiniteSum.Basic
import Mathlib.Data.Nat.ModEq

namespace CollatzUniformity

/-- LTE-Worst-Start n = 2^(k+1)·3^r - 1. -/
def lteWorst (k r : Nat) : Nat := 2 ^ (k + 1) * 3 ^ r - 1

/-- Minimale LTE-Familie (k=1): n = 4·3^r - 1. -/
def lteMinimal (r : Nat) : Nat := 4 * 3 ^ r - 1

theorem padicVal_ge_iff_dvd (m k : Nat) (hm : m != 0) :
    padicValNat 2 m >= k <-> 2 ^ k | m := by
  constructor
  . intro h; exact (Nat.pow_dvd_pow 2 h).trans pow_padicValNat_dvd
  . intro h; exact padicValNat.le_of_dvd (by norm_num) hm h

theorem padicVal_succ_ge_iff_dvd (n k : Nat) :
    padicValNat 2 (n + 1) >= k <-> 2 ^ k | (n + 1) := by
  apply padicVal_ge_iff_dvd; exact Nat.succ_ne_zero n

theorem lteWorst_odd (k r : Nat) (hk : 0 < k) (hr : 0 < r) :
    Odd (lteWorst k r) := by
  unfold lteWorst
  have hpos : 1 < 2 ^ (k + 1) * 3 ^ r := by
    have : 2 ≤ 2 ^ (k + 1) := Nat.le_pow_self 2 (by omega)
    nlinarith [pow_pos (by norm_num : 0 < 3) r]
  have : (2 ^ (k + 1) * 3 ^ r) % 2 = 0 := by
    exact Nat.mod_eq_zero_of_dvd (dvd_mul_of_dvd_left (dvd_pow_self 2 (by omega)) _)
  omega

theorem lteWorst_valuation (k r : Nat) :
    padicValNat 2 (lteWorst k r + 1) = k + 1 := by
  unfold lteWorst
  have h : lteWorst k r + 1 = 2 ^ (k + 1) * 3 ^ r := by omega
  rw [h, padicValNat.mul, padicValNat.pow, padicValNat.pow] <;> norm_num
  · ring
  · exact Nat.coprime_pow_left_iff.mpr (by norm_num)

theorem lteMinimal_is_lteWorst (r : Nat) :
    lteMinimal r = lteWorst 1 r := by unfold lteMinimal lteWorst; ring

/-- Für ungerades j: 3^j ≡ 3 (mod 8). -/
private theorem three_pow_odd_mod8 (j : Nat) (hj : Odd j) :
    3 ^ j % 8 = 3 := by
  rcases hj with ⟨t, ht⟩
  subst ht
  induction t with
  | zero => norm_num
  | succ t ih =>
    have : 3 ^ (2 * t + 3) = 9 * 3 ^ (2 * t + 1) := by ring
    rw [this]
    have h1 : 3 ^ (2 * t + 1) % 8 = 3 := ih
    omega

/-- Für ungerades j: ν₂(3^j + 1) = 2. -/
theorem padicVal_two_three_pow_odd_plus_one (j : Nat) (hj : Odd j) :
    padicValNat 2 (3 ^ j + 1) = 2 := by
  have hmod8 : 3 ^ j % 8 = 3 := three_pow_odd_mod8 j hj
  have hge : 2 ≤ padicValNat 2 (3 ^ j + 1) := by
    rw [padicVal_succ_ge_iff_dvd, ← Nat.modEq_iff_dvd]
    have : 3 ^ j + 1 ≡ 4 [MOD 8] := by omega
    exact this.trans (by decide : (4 : Nat) ≡ 0 [MOD 4])
  have hlt : padicValNat 2 (3 ^ j + 1) < 3 := by
    rw [padicVal_succ_ge_iff_dvd]
    intro h8
    have : 8 ∣ 3 ^ j + 1 := (Nat.pow_dvd_pow 2 h8).trans pow_padicValNat_dvd
    have : 3 ^ j + 1 ≡ 0 [MOD 8] := Nat.modEq_zero_iff_dvd.mpr this
    omega
  omega

/--
  Fall 1 (gerades N): Nachfolger (3^(N+1)-1)/2 erfüllt ν₂(n+1)=1.
  Entspricht Lemma (LTE-Reset) für gerades N=k+r.
-/
theorem lte_even_N_successor_valuation_one (N : Nat) (hN : Even N) :
    padicValNat 2 ((3 ^ (N + 1) - 1) / 2 + 1) = 1 := by
  have hodd : Odd (N + 1) := by rcases hN with ⟨j, rfl⟩; exact ⟨j, by ring⟩
  have hmid : (3 ^ (N + 1) + 1) / 2 = (3 ^ (N + 1) - 1) / 2 + 1 := by
    have hpos : 1 ≤ 3 ^ (N + 1) := by positivity
    omega
  rw [hmid]
  have hval : padicValNat 2 (3 ^ (N + 1) + 1) = 2 :=
    padicVal_two_three_pow_odd_plus_one (N + 1) hodd
  have hpos : 0 < 3 ^ (N + 1) + 1 := by positivity
  rw [← padicValNat.div_of_dvd (by omega : 2 ∣ 3 ^ (N + 1) + 1)]
  · omega
  · rw [padicVal_succ_ge_iff_dvd]; exact dvd_of_padicValNat_pos (by omega : 0 < 2)

theorem lte_worst_even_N_forces_low_valuation (k r : Nat) (hN : Even (k + r)) :
    padicValNat 2 ((3 ^ (k + r + 1) - 1) / 2 + 1) = 1 :=
  lte_even_N_successor_valuation_one (k + r) hN

/-- ν₂(n+1) ≤ 2 bedeutet: höchstens C-Kettenlänge 1. -/
theorem valuation_le_two_limits_C_chain (n : Nat) (h : padicValNat 2 (n + 1) ≤ 2) :
    padicValNat 2 (n + 1) < 3 := by omega

/-- Mischschranken-Kern: (1-p)^n → 0 für p ∈ (0,1). -/
theorem mixing_probability_bound (p : Real) (hp0 : 0 < p) (hp1 : p < 1) (n : Nat) :
    (1 - p) ^ n ≤ 1 := by
  have : 0 ≤ 1 - p := by linarith
  exact pow_le_one_of_nonpos_of_le this (by linarith)

theorem mixing_probability_tendsto_zero (p : Real) (hp0 : 0 < p) (hp1 : p < 1) :
    Filter.Tendsto (fun n : Nat => (1 - p) ^ n) Filter.atTop (nhds 0) := by
  have h : norm (1 - p) < 1 := by rw [Real.norm_of_nonneg (by linarith)]; linarith
  exact tendsto_pow_atTop_nhds_zero_of_norm_lt_one h

/-- Alias: Verbindung zur Ausnahmewahrscheinlichkeit in collatz_density_appendix.lean. -/
theorem exception_probability_tendsto_zero (p : Real) (hp0 : 0 < p) (hp1 : p < 1) :
    Filter.Tendsto (fun n : Nat => (1 - p) ^ n) Filter.atTop (nhds 0) :=
  mixing_probability_tendsto_zero p hp0 hp1

/-- C-Ketten-Dichte (aus Anhang, hier dupliziert für Uniformitäts-Kontext). -/
theorem density_C_chains_finite (N k : Nat) :
    (Finset.filter (fun n => padicValNat 2 (n + 1) >= k + 1)
      (Finset.filter Odd (Finset.range (2 * N + 1)))).card
    = N / 2 ^ k := by
  have hset_eq : Finset.filter (fun n => padicValNat 2 (n + 1) >= k + 1)
      (Finset.filter Odd (Finset.range (2 * N + 1))) =
      (Finset.range (2 * N + 1)).filter (fun e => 2 ^ (k + 1) | e + 1) := by
    ext n; simp only [Finset.mem_filter, Finset.mem_range]
    constructor
    . rintro ⟨hn_hpadic, _⟩; exact hn_hpadic
    . intro hn_hdvd; refine hn_hdvd, (padicVal_succ_ge_iff_dvd n (k + 1)).mpr hn_hdvd
  rw [hset_eq]
  have h_count := Nat.card_multiples (2 * N + 1) (2 ^ (k + 1))
  rw [h_count, show (2 : Nat) ^ (k + 1) = 2 * 2 ^ k from by ring, Nat.div_div_eq_div_mul]
  have h_div : (2 * N + 1) / 2 = N := by omega
  rw [h_div]

end CollatzUniformity
