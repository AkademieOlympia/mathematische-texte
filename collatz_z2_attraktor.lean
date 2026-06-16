/-
  collatz_z2_attraktor.lean
  Ausnahmemenge E ⊂ ℤ₂ und dist₂(n,E) für den Collatz-Uniformitätsansatz.

  PadicInt-Grundstruktur (Mathlib 4.29): ‖x-y‖₂, Einbettung ℕ ↪ ℤ₂, sInf dist₂(·,E).
  Valuations-Lemmas für ℕ: vgl. collatz_uniformity.lean (lteWorst, padicValNat).
-/

import Mathlib.NumberTheory.Padics.PadicIntegers
import Mathlib.NumberTheory.Padics.PadicVal.Basic
import Mathlib.Algebra.Ring.Parity
import Mathlib.Data.Real.Archimedean
import Mathlib.Topology.MetricSpace.Ultra.Basic

namespace CollatzZ2

open scoped BigOperators

/-- 2-adische ganze Zahlen ℤ₂ (p = 2). -/
abbrev Z2 := PadicInt 2

/-- Ungerade natürliche Zahl als 2-adische Einheit in ℤ₂. -/
noncomputable def oddZ2Embed (n : ℕ) (_h : Odd n) : Z2 :=
  (n : Z2)

/-- 2-adische Norm |x-y|₂ für x,y ∈ ℤ₂. -/
noncomputable def dist2 (x y : Z2) : ℝ :=
  ‖x - y‖

/-- Odd-to-odd Collatz-Schritt auf ungeraden n (natürliche Implementierung). -/
def collatzU (n : ℕ) (_h : Odd n) : ℕ :=
  let m := 3 * n + 1
  m / 2 ^ (padicValNat 2 m)

/-- Alias für Abwärtskompatibilität mit dem ersten Gerüst. -/
abbrev U_odd := collatzU

/-- Natürliche Zahl mit endlicher 2-adischer Entwicklung (jedes n ∈ ℕ liegt in ℤ₂). -/
def isNatural (x : Z2) : Prop :=
  ∃ n : ℕ, x = (n : Z2)

/-- E_K: Punkte mit Komplexität > K — Platzhalter bis Collatz-Bahn-Analyse formalisiert ist. -/
def highComplexity (_K : ℕ) (_x : Z2) : Prop :=
  -- TODO: echte Komplexitätsfunktion auf ℤ₂; derzeit nur Typ-Skelett.
  False

/-- Endliche Ausnahme-Approximation E_K (Platzhalter: noch keine Collatz-Bahn-Analyse). -/
def ExceptionSetApprox (_K : ℕ) : Finset ℕ :=
  Finset.empty

/-- dist₂(n, E) für n ∈ ℕ: Infimum der 2-adischen Abstände zu E ⊂ ℤ₂. -/
noncomputable def distToExceptionSet (n : ℕ) (E : Set Z2) : ℝ :=
  sInf ((fun e => dist2 (n : Z2) e) '' E)

/-- 2-adische Distanz auf ℕ via ν₂(a-b) (rationales Pendant zu dist2). -/
noncomputable def dist2Nat (a b : ℕ) : ℚ :=
  if _h : a = b then 0
  else (2 : ℚ) ^ (-(padicValNat 2 (a - b) : ℤ))

/-- PadicInt-Einbettung natürlicher Zahlen. -/
noncomputable def natToPadic (n : ℕ) : Z2 :=
  (n : Z2)

/-! ### 2-adische Metrik: Grund-Eigenschaften (PadicInt.instMetricSpace) -/

theorem dist2_nonneg (x y : Z2) : 0 ≤ dist2 x y :=
  norm_nonneg _

theorem dist2_comm (x y : Z2) : dist2 x y = dist2 y x := by
  unfold dist2
  rw [← dist_eq_norm, ← dist_eq_norm]
  exact norm_sub_rev x y

theorem dist2_eq_zero_iff {x y : Z2} : dist2 x y = 0 ↔ x = y := by
  have hdist : dist x y = dist2 x y := by
    unfold dist2
    rw [dist_eq_norm, norm_sub_rev]
  simpa [hdist] using dist_eq_zero (x := x) (y := y)

theorem dist2_self (x : Z2) : dist2 x x = 0 := by
  simpa using dist2_eq_zero_iff.mpr rfl

/-- Ultrametrische Ungleichung (‖x-z‖ ≤ max ‖x-y‖ ‖y-z‖). -/
theorem dist2_ultrametric (x y z : Z2) : dist2 x z ≤ max (dist2 x y) (dist2 y z) := by
  unfold dist2
  simpa [dist_eq_norm] using dist_triangle_max x y z

/-! ### ℕ ↪ ℤ₂: Norm von Differenzen und padicValNat -/

private lemma valuation_natCast (n : ℕ) : (n : Z2).valuation = padicValNat 2 n := by
  unfold PadicInt.valuation
  simp [Padic.valuation_natCast]

private theorem nat_cast_sub (a b : ℕ) (hab : b ≤ a) :
    (a : Z2) - (b : Z2) = ((a - b : ℕ) : Z2) := by
  rw [← Nat.cast_sub hab]

/-- Für a ≠ b: ‖(a:ℤ₂) - (b:ℤ₂)‖ = 2^{-ν₂(a-b)} (a ≥ b). -/
theorem dist2_nat_cast_of_le {a b : ℕ} (hab : b ≤ a) (_hne : a ≠ b) :
    dist2 (a : Z2) (b : Z2) = (2 : ℝ) ^ (-(padicValNat 2 (a - b) : ℤ)) := by
  unfold dist2
  rw [nat_cast_sub a b hab]
  have hnz : (a - b : ℕ) ≠ 0 := by omega
  have hx : ((a - b : ℕ) : Z2) ≠ 0 := Nat.cast_ne_zero.mpr hnz
  rw [PadicInt.norm_eq_zpow_neg_valuation hx, valuation_natCast]
  rfl

theorem dist2_nat_one (n : ℕ) (hn : 1 < n) :
    dist2 (n : Z2) (1 : Z2) = (2 : ℝ) ^ (-(padicValNat 2 (n - 1) : ℤ)) := by
  have hab : 1 ≤ n := by omega
  have hne : n ≠ 1 := by omega
  exact dist2_nat_cast_of_le hab hne

/-! ### LTE-Familie und Abstand zu 1 -/

/-- LTE-Minimalstart n = 4·3^r - 1 (vgl. CollatzUniformity.lteMinimal). -/
def lteMinimal (r : ℕ) : ℕ := 4 * 3 ^ r - 1

theorem lteMinimal_odd (r : ℕ) (_hr : 0 < r) : Odd (lteMinimal r) := by
  unfold lteMinimal
  refine ⟨2 * 3 ^ r - 1, ?_⟩
  have hpow : 1 ≤ 3 ^ r := Nat.one_le_pow r 3 (by norm_num)
  ring_nf
  omega

private lemma odd_two_mul_sub_one (r : ℕ) : Odd (2 * 3 ^ r - 1) := by
  refine ⟨3 ^ r - 1, ?_⟩
  have hpow : 1 ≤ 3 ^ r := Nat.one_le_pow r 3 (by norm_num)
  ring_nf
  omega

theorem lteMinimal_dist_to_one (r : ℕ) (_hr : 0 < r) :
    dist2 (lteMinimal r : Z2) (1 : Z2) = (1 / 2 : ℝ) := by
  have hn : 1 < lteMinimal r := by
    unfold lteMinimal
    have : 4 ≤ 4 * 3 ^ r := by nlinarith [pow_pos (by norm_num : 0 < 3) r]
    omega
  rw [dist2_nat_one (lteMinimal r) hn]
  have hval : padicValNat 2 (lteMinimal r - 1) = 1 := by
    unfold lteMinimal
    have h : 4 * 3 ^ r - 1 - 1 = 2 * (2 * 3 ^ r - 1) := by omega
    rw [h]
    have hodd := odd_two_mul_sub_one r
    have hnz : 2 * 3 ^ r - 1 ≠ 0 := by
      have : 1 < 2 * 3 ^ r := by nlinarith [pow_pos (by norm_num : 0 < 3) r]
      omega
    have hv : padicValNat 2 (2 * 3 ^ r - 1) = 0 := by
      rw [padicValNat.eq_zero_iff]
      refine Or.inr (Or.inr ?_)
      rcases hodd with ⟨k, hk⟩
      omega
    rw [padicValNat.mul (by norm_num) hnz, padicValNat_self, hv]
  rw [hval]
  norm_num

theorem lteMinimal_dist2Nat_to_one (r : ℕ) (_hr : 0 < r) :
    dist2Nat (lteMinimal r) 1 = (1 / 2 : ℚ) := by
  have hn : 1 < lteMinimal r := by
    unfold lteMinimal
    have : 4 ≤ 4 * 3 ^ r := by nlinarith [pow_pos (by norm_num : 0 < 3) r]
    omega
  unfold dist2Nat
  have hne : lteMinimal r ≠ 1 := by
    unfold lteMinimal
    have : 4 ≤ 4 * 3 ^ r := by nlinarith [pow_pos (by norm_num : 0 < 3) r]
    omega
  simp [hne]
  have hval : padicValNat 2 (lteMinimal r - 1) = 1 := by
    unfold lteMinimal
    have h : 4 * 3 ^ r - 1 - 1 = 2 * (2 * 3 ^ r - 1) := by omega
    rw [h]
    have hodd := odd_two_mul_sub_one r
    have hnz : 2 * 3 ^ r - 1 ≠ 0 := by
      have : 1 < 2 * 3 ^ r := by nlinarith [pow_pos (by norm_num : 0 < 3) r]
      omega
    have hv : padicValNat 2 (2 * 3 ^ r - 1) = 0 := by
      rw [padicValNat.eq_zero_iff]
      refine Or.inr (Or.inr ?_)
      rcases hodd with ⟨k, hk⟩
      omega
    rw [padicValNat.mul (by norm_num) hnz, padicValNat_self, hv]
  rw [hval]
  norm_num

/-- dist₂(n,1) ist NICHT die Uniformitätsform (LTE-Familie liefert konstant 1/2). -/
theorem dist_to_one_not_uniform_bound :
    ∃ n : ℕ, dist2 (n : Z2) (1 : Z2) = (1 / 2 : ℝ) := by
  refine ⟨lteMinimal 1, ?_⟩
  exact lteMinimal_dist_to_one 1 (by norm_num)

/-- Alias für das erste Gerüst (Rat statt ℝ). -/
theorem dist_to_one_not_the_right_attraktor :
    ∃ n : ℕ, dist2Nat n 1 = 1 / 2 := by
  refine ⟨lteMinimal 1, ?_⟩
  exact lteMinimal_dist2Nat_to_one 1 (by norm_num)

/-! ### isNatural und distToExceptionSet -/

theorem isNatural_nat (n : ℕ) : isNatural (n : Z2) :=
  ⟨n, rfl⟩

theorem distToExceptionSet_empty (n : ℕ) :
    distToExceptionSet n (∅ : Set Z2) = 0 := by
  unfold distToExceptionSet
  simp [Real.sInf_empty]

/-- Jeder Abstand zu E liegt unterhalb des Infimums über alle e ∈ E. -/
theorem distToExceptionSet_le_dist (n : ℕ) (E : Set Z2) (e : Z2) (he : e ∈ E) :
    distToExceptionSet n E ≤ dist2 (n : Z2) e := by
  unfold distToExceptionSet
  refine csInf_le ?_ ⟨e, he, rfl⟩
  exact ⟨0, fun r hr => by
    rcases hr with ⟨_, _, rfl⟩
    exact dist2_nonneg _ _⟩

/-- Austritt aus E_K nach endlich vielen Schritten (Uniformitäts-Vermutung) — noch offen. -/
def exitsExceptionApprox (_n _K : ℕ) : Prop :=
  -- TODO: Iteration von collatzU formalisieren (Function.iterate).
  True

/-- Triviale 2-adische Nachbarschaft von 1: noch zu formalisieren als ‖x-1‖ < 1. -/
def inTrivialAttraktorApprox (n : ℕ) : Prop :=
  n % 2 = 1 ∧ n = 1

end CollatzZ2
