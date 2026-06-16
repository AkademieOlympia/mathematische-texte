/-
  collatz_density.lean (Anhang-Auszug, ASCII)
  Dichte-Lemmata fuer C-Ketten in der Collatz-Dynamik.
-/

import Mathlib.NumberTheory.Padics.PadicVal
import Mathlib.Data.Nat.Parity
import Mathlib.Topology.Algebra.InfiniteSum.Basic

theorem padicVal_ge_iff_dvd (m k : Nat) (hm : m != 0) :
    padicValNat 2 m >= k <-> 2 ^ k | m := by
  constructor
  . intro h
    exact (Nat.pow_dvd_pow 2 h).trans pow_padicValNat_dvd
  . intro h
    exact padicValNat.le_of_dvd (by norm_num) hm h

theorem padicVal_succ_ge_iff_dvd (n k : Nat) :
    padicValNat 2 (n + 1) >= k <-> 2 ^ k | (n + 1) := by
  apply padicVal_ge_iff_dvd
  exact Nat.succ_ne_zero n

theorem C_chain_start_class (k : Nat) :
    forall n : Nat, Odd n /\ padicValNat 2 (n + 1) >= k + 1 <->
    exists m : Nat, n = 2 ^ (k + 1) * m + (2 ^ (k + 1) - 1) := by
  intro n
  constructor
  . intro hodd_hval
    rw [padicVal_succ_ge_iff_dvd] at hodd_hval
    obtain m, hm := hodd_hval
    exact m, by omega
  . intro m_hm
    constructor
    . rw [m_hm]; simp [Nat.Odd]; omega
    . rw [padicVal_succ_ge_iff_dvd, m_hm]; ring_nf; exact m, by ring

theorem density_C_chains_finite (N k : Nat) :
    (Finset.filter (fun n => padicValNat 2 (n + 1) >= k + 1)
      (Finset.filter Odd (Finset.range (2 * N + 1)))).card
    = N / 2 ^ k := by
  have hset_eq : Finset.filter (fun n => padicValNat 2 (n + 1) >= k + 1)
      (Finset.filter Odd (Finset.range (2 * N + 1))) =
      (Finset.range (2 * N + 1)).filter (fun e => 2 ^ (k + 1) | e + 1) := by
    ext n; simp only [Finset.mem_filter, Finset.mem_range]
    constructor
    . rintro hn_hpadic; exact hn_hpadic
    . intro hn_hdvd
      refine hn_hdvd, (padicVal_succ_ge_iff_dvd n (k + 1)).mpr hn_hdvd
  rw [hset_eq]
  have h_count := Nat.card_multiples (2 * N + 1) (2 ^ (k + 1))
  rw [h_count]
  rw [show (2 : Nat) ^ (k + 1) = 2 * 2 ^ k from by ring, Nat.div_div_eq_div_mul]
  have h_div : (2 * N + 1) / 2 = N := by omega
  rw [h_div]

theorem density_C_chains_bound (k : Nat) :
    forall N : Nat, 0 < N ->
    (Finset.filter (fun n => padicValNat 2 (n + 1) >= k + 1)
      (Finset.filter Odd (Finset.range (2 * N + 1)))).card
    <= N / 2 ^ k + 1 := by
  intro N _; rw [density_C_chains_finite]; exact Nat.le_add_right _ _

private theorem geom_series_deriv (x : Real) (hx : |x| < 1) :
    tsum j : Nat, (j : Real) * x ^ j = x / (1 - x) ^ 2 := by
  have h := hasSum_pow_mul_geometric_of_abs_lt_one 1 hx
  simp only [pow_one] at h; exact h.tsum_eq

theorem tail_series_formula (K : Nat) :
    tsum k : Nat, (K + k + 1 : Real) * (1 / 2) ^ (K + k + 1)
    = (K + 2 : Real) * (1 / 2) ^ K := by
  have h_split : forall k : Nat, (K + k + 1 : Real) * (1 / 2) ^ (K + k + 1) =
      (K + 1 : Real) * (1 / 2) ^ (K + 1) * (1 / 2) ^ k +
      (k : Real) * (1 / 2) ^ (K + 1) * (1 / 2) ^ k := by
    intro k; push_cast; ring
  simp_rw [h_split]; rw [tsum_add]
  . rw [tsum_mul_left, tsum_geometric_two]
    have h_deriv : tsum k : Nat, (k : Real) * (1 / 2 : Real) ^ k = 2 := by
      have hgd := geom_series_deriv (1 / 2 : Real) (by norm_num)
      have : (1 / 2 : Real) / (1 - 1 / 2) ^ 2 = 2 := by norm_num
      linarith
    have h_second : tsum k : Nat, (k : Real) * (1 / 2 : Real) ^ (K + 1) * (1 / 2) ^ k =
        (1 / 2 : Real) ^ (K + 1) * 2 := by
      simp_rw [fun k => by ring]; rw [tsum_mul_left, h_deriv]
    rw [h_second]; push_cast; ring
  . apply Summable.const_smul; exact summable_geometric_two
  . apply Summable.congr (summable_pow_mul_geometric_of_norm_lt_one 1 (by norm_num))
    intro k; simp [pow_one]

theorem collatz_norm_identity (n : Int) : (3 * n + 1) ^ 2 = 9 * n ^ 2 + 6 * n + 1 := by ring

theorem collatz_norm_before_halving (n : Int) :
    (3 * n + 1) ^ 2 > n ^ 2 <-> n ^ 2 + 3 * n > 0 := by
  constructor <;> intro h <;> nlinarith [sq_nonneg n]

theorem exception_probability_decay (p : Real) (hp0 : 0 < p) (hp1 : p < 1) (n : Nat) :
    (1 - p) ^ n <= (1 - p) ^ 0 := by
  apply pow_le_one_of_nonpos_of_le <;> linarith

theorem exception_probability_tendsto_zero (p : Real) (hp0 : 0 < p) (hp1 : p < 1) :
    Filter.Tendsto (fun n : Nat => (1 - p) ^ n) Filter.atTop (nhds 0) := by
  have h : norm (1 - p) < 1 := by rw [Real.norm_of_nonneg (by linarith)]; linarith
  exact tendsto_pow_atTop_nhds_zero_of_norm_lt_one h
