import Mathlib

namespace Collatz.CkA

def oddSyracuseEquation (n u a : ℕ) : Prop :=
  Odd n ∧
  0 < a ∧
  u * 2 ^ a = 3 * n + 1

def TwoAdicBlockEquation
    (n : ℕ → ℕ)
    (a : ℕ → ℕ)
    (r : ℕ) : Prop :=
  ∀ j < r,
    oddSyracuseEquation (n j) (n (j + 1)) (a j)

def StrictNatDescent (x y : ℕ) : Prop :=
  y < x

theorem one_step_descent_of_large_valuation
    {n u a : ℕ}
    (hn : 1 < n)
    (hEq : u * 2 ^ a = 3 * n + 1)
    (ha : 2 ≤ a) :
    u < n := by
  have h4 : 4 ≤ 2 ^ a := by
    calc
      4 = 2 ^ 2 := by norm_num
      _ ≤ 2 ^ a := Nat.pow_le_pow_right (by norm_num) ha
  have hu4 : 4 * u ≤ 2 ^ a * u := by
    exact Nat.mul_le_mul_right u h4
  have hrewrite : 2 ^ a * u = 3 * n + 1 := by
    simpa [Nat.mul_comm] using hEq
  have hbound : 4 * u ≤ 3 * n + 1 := by
    simpa [hrewrite] using hu4
  omega

end Collatz.CkA
