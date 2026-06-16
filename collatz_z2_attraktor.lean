/-
  collatz_z2_attraktor.lean
  Gerüst: Ausnahmemenge E ⊂ Z₂ und dist₂(n,E) für den Collatz-Uniformitätsansatz.

  Nächste Mathlib-Nutzung: padicValNat (bereits in collatz_uniformity.lean),
  später PadicInt.norm für |x-y|₂ auf Z₂.
-/

import Mathlib.NumberTheory.Padics.PadicInt
import Mathlib.NumberTheory.Padics.PadicVal
import Mathlib.Data.Nat.Parity

namespace CollatzZ2

/-- Odd-to-odd Collatz-Schritt auf ungeraden n (natürliche Implementierung). -/
def U_odd (n : Nat) (h : Odd n) : Nat :=
  let m := 3 * n + 1
  m / 2 ^ (padicValNat 2 m)

/-- 2-adische Distanz auf ℕ via ν₂(a-b): dist₂(a,b) = 2^{-ν₂(a-b)} (als ℚ). -/
noncomputable def dist2Nat (a b : Nat) : Rat :=
  if h : a = b then 0
  else (2 : Rat) ^ (-(padicValNat 2 (a - b) : Int))

/-- Endliche Ausnahme-Approximation E_K (Platzhalter: noch keine Collatz-Bahn-Analyse). -/
def ExceptionSetApprox (K : Nat) : Finset Nat :=
  Finset.empty

/-- Austritt aus E_K nach endlich vielen Schritten (Uniformitäts-Vermutung). -/
def exitsExceptionApprox (n K : Nat) : Prop :=
  -- TODO: Iteration von U_odd formalisieren (Function.iterate).
  True

/-- Triviale 2-adische Nachbarschaft von 1: n ≡ 1 (mod 2) und n in A_triv — noch zu formalisieren. -/
def inTrivialAttraktorApprox (n : Nat) : Prop :=
  n % 2 = 1 ∧ n = 1

/-- Widerlegte Form: dist₂(n,1) ≥ c·2^{-log n} — hier nur als negiertes Ziel skizziert. -/
theorem dist_to_one_not_the_right_attraktor :
    ∃ n : Nat, dist2Nat n 1 = 1 / 2 := by
  -- LTE-Familie: n = 4·3^r - 1 liefert ν₂(n-1) = 2 für r ≥ 1.
  sorry

/-- PadicInt-Einbettung: nächster Schritt nach dist2Nat. -/
noncomputable def natToPadic (n : Nat) : PadicInt 2 :=
  (n : PadicInt 2)

-- TODO: dist₂(x,E) = inf_{e∈E} ‖x-e‖₂ mit PadicInt.norm (Mathlib-Lücke im Gerüst)

end CollatzZ2
