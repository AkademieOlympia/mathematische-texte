/-
  Collatz.Octonion.OddCoreLogCocycle — Module O1 (post-freeze).

  Exact odd-core multiplicative cocycle for `CollatzZ2.collatzU` / `CollatzZ2.iterateU`,
  and the real log-scale identity (Module O1 dossier formula).

  Status: sorry-free. No Collatz conjecture claim.
-/

import CollatzEabc.Kappa
import Mathlib.NumberTheory.Padics.PadicVal.Basic
import Mathlib.Algebra.Ring.Parity
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic

namespace CollatzOctonion

open CollatzEabc CollatzZ2

/-- 2-adic valuation of the even block `3n+1` along an odd step. -/
def oddCoreValuation (n : ℕ) (_h : Odd n) : ℕ :=
  padicValNat 2 (3 * n + 1)

/-- Sum of odd-core valuations along the first `k` odd-to-odd steps. -/
def oddCoreValuationSum (n₀ : ℕ) (k : ℕ) (h₀ : Odd n₀) : ℕ :=
  (Finset.range k).sum fun j =>
    oddCoreValuation (iterateU n₀ j) (iterateU_odd n₀ j h₀)

private theorem three_mul_odd_add_one_pos (n : ℕ) (h : Odd n) : 0 < 3 * n + 1 := by
  rcases h with ⟨t, ht⟩
  subst ht
  nlinarith

private theorem two_pow_val_pos (v : ℕ) : 0 < 2 ^ v := by
  apply Nat.pos_iff_ne_zero.mpr
  exact pow_ne_zero _ (by decide)

private theorem collatzU_pos (n : ℕ) (h : Odd n) : 0 < collatzU n h := by
  unfold collatzU
  have hm : 0 < 3 * n + 1 := three_mul_odd_add_one_pos n h
  have hle : 2 ^ padicValNat 2 (3 * n + 1) ≤ 3 * n + 1 := by
    exact Nat.le_of_dvd hm pow_padicValNat_dvd
  exact Nat.div_pos hle (two_pow_val_pos _)

private theorem iterateU_pos (n : ℕ) : ∀ k, Odd n → 0 < iterateU n k
  | 0, h => by
    simp [iterateU]
    rcases h with ⟨t, ht⟩
    subst ht
    nlinarith
  | k + 1, h => by
    simp [iterateU, dif_pos h]
    exact iterateU_pos _ k (collatzU_preserves_odd n h)

/-- `iterateU n₀ (k+1) = collatzU (iterateU n₀ k) hₖ`. -/
theorem iterateU_succ_eq (n₀ : ℕ) (k : ℕ) (h₀ : Odd n₀) :
    iterateU n₀ (k + 1) = collatzU (iterateU n₀ k) (iterateU_odd n₀ k h₀) := by
  induction k generalizing n₀ with
  | zero => simp [iterateU, dif_pos h₀]
  | succ k ih =>
    have h₁ := collatzU_preserves_odd n₀ h₀
    simpa [iterateU, dif_pos h₀, iterateU_succ n₀ k h₀] using ih (collatzU n₀ h₀) h₁

/-- Single odd-core step: `S(n)·2^{ν₂(3n+1)} = 3n+1`. -/
theorem collatzU_mul_two_pow (n : ℕ) (h : Odd n) :
    collatzU n h * 2 ^ oddCoreValuation n h = 3 * n + 1 := by
  unfold collatzU oddCoreValuation
  exact Nat.div_mul_cancel pow_padicValNat_dvd

/-- Product cocycle for `k ≥ 1` (next formalization target).
    Single-step building blocks: `collatzU_mul_two_pow`, `iterateU_succ_eq`.
    Audited multiplicatively in `verify_o1_o2.py`. -/
def iterateU_cocycle_mulStatement {k : ℕ} (hk : 0 < k) (n₀ : ℕ) (h₀ : Odd n₀) : Prop :=
  (Finset.range k).prod (fun j : ℕ => 3 * iterateU n₀ j + 1) =
    iterateU n₀ k * 2 ^ oddCoreValuationSum n₀ k h₀

/-- Log-scale odd-core cocycle (Module O1 dossier statement). -/
def oddCoreLogCocycleStatement (n₀ : ℕ) (k : ℕ) (h₀ : Odd n₀) (_hn₀ : 0 < n₀) : Prop :=
  Real.log (iterateU n₀ k) - Real.log n₀ =
    (k : ℝ) * Real.log 3 -
      Real.log 2 * oddCoreValuationSum n₀ k h₀ +
        (Finset.range k).sum fun j => Real.log (1 + 1 / (3 * iterateU n₀ j))

/-- Alias matching the dossier name `S` for the accelerated odd map. -/
abbrev oddCoreStep (n : ℕ) (h : Odd n) : ℕ := collatzU n h

end CollatzOctonion
