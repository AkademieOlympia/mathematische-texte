/-
  collatz_z2_attraktor.lean
  Ausnahmemenge E ⊂ ℤ₂ und dist₂(n,E) für den Collatz-Uniformitätsansatz.

  PadicInt-Grundstruktur (Mathlib 4.29): ‖x-y‖₂, Einbettung ℕ ↪ ℤ₂, sInf dist₂(·,E).
  Valuations-Lemmas für ℕ: vgl. collatz_uniformity.lean (lteWorst, padicValNat).

  Roadmap (PR #20 Nachfolge):
    Stufe A — ExceptionSetApprox N K
    Stufe B — Monotonie in N und K
    Stufe C — Limes E = closure(⋃_N E_{N,N})
    Stufe D — U-Invarianz `collatzU_maps_exception_approx`
    Stufe E — Uniformität dist₂(T^k(n),E) → 0 (offen)
-/

import Mathlib.NumberTheory.Padics.PadicIntegers
import Mathlib.NumberTheory.Padics.PadicVal.Basic
import Mathlib.Algebra.Ring.Parity
import Mathlib.Data.Real.Archimedean
import Mathlib.Topology.MetricSpace.Ultra.Basic
import Mathlib.Topology.Basic
import Mathlib.Topology.Closure

namespace CollatzZ2

open scoped BigOperators
open Filter Topology

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

/-! ### Stufe A: iterateU, trivialer Attraktor, ExceptionSetApprox -/

/-- Nach `k` odd-to-odd-Schritten: `U^k(n)` (für ungerades `n`). -/
def iterateU (n : ℕ) : ℕ → ℕ
  | 0 => n
  | k + 1 =>
    if h : Odd n then iterateU (collatzU n h) k else n

private lemma collatzU_odd (n : ℕ) (h : Odd n) : Odd (collatzU n h) := by
  rcases (Nat.even_xor_odd (collatzU n h)).or with heven | hod
  · exfalso
    unfold collatzU at heven
    set m := 3 * n + 1 with hm
    set v := padicValNat 2 m
    set q := m / 2 ^ v
    have hm_ne : m ≠ 0 := by
      rcases h with ⟨k, rfl⟩
      omega
    have hdiv : 2 ^ v ∣ m := pow_padicValNat_dvd
    have heq : q * 2 ^ v = m := Nat.div_mul_cancel hdiv
    have h2 : 2 ∣ q := even_iff_two_dvd.mp heven
    have hpow₁ : 2 ^ v * 2 ∣ 2 ^ v * q := Nat.mul_dvd_mul_left (2 ^ v) h2
    have hpow₂ : 2 ^ (v + 1) ∣ q * 2 ^ v := by simpa [pow_succ', Nat.mul_comm] using hpow₁
    have hpow : 2 ^ (v + 1) ∣ m := by simpa [← heq] using hpow₂
    exact pow_succ_padicValNat_not_dvd (p := 2) hm_ne hpow
  · exact hod

theorem iterateU_odd (n : ℕ) (k : ℕ) (h : Odd n) : Odd (iterateU n k) := by
  induction k generalizing n with
  | zero => exact h
  | succ k ih =>
    rw [iterateU, dif_pos h]
    exact ih (collatzU n h) (collatzU_odd n h)

theorem iterateU_succ (n : ℕ) (k : ℕ) (h : Odd n) :
    iterateU n (k + 1) = iterateU (collatzU n h) k := by
  simp [iterateU, h]

theorem collatzU_le (n : ℕ) (h : Odd n) (N : ℕ) (hn : n ≤ N) :
    collatzU n h ≤ 3 * N + 1 := by
  unfold collatzU
  have hm : 3 * n + 1 ≤ 3 * N + 1 := by nlinarith
  exact Nat.le_trans (Nat.div_le_self _ _) hm

/-- Triviale 2-adische Kugel `A_triv = {x : ‖x-1‖₂ < 1}` (TeX §2). -/
def inTrivialAttraktorBall (x : Z2) : Prop :=
  dist2 x 1 < 1

/-- Jedes ungerade `n ∈ ℕ` liegt in `A_triv` (‖n-1‖₂ ≤ 1/2 < 1). -/
theorem odd_nat_in_trivial_ball (n : ℕ) (h : Odd n) :
    inTrivialAttraktorBall (n : Z2) := by
  unfold inTrivialAttraktorBall
  rcases h with ⟨k, rfl⟩
  by_cases hk : k = 0
  · subst hk
    simp [dist2_self]
  · have hn' : 1 < 2 * k + 1 := by omega
    rw [dist2_nat_one (2 * k + 1) hn']
    simp only [show (2 * k + 1) - 1 = 2 * k by omega]
    have hv : 1 ≤ padicValNat 2 (2 * k) :=
      (padicValNat_dvd_iff_le (show 2 * k ≠ 0 by omega)).mp ⟨k, rfl⟩
    have hle : (2 : ℝ) ^ (-(padicValNat 2 (2 * k) : ℤ)) ≤ (1 / 2 : ℝ) := by
      calc
        _ ≤ (2 : ℝ) ^ (-(1 : ℤ)) :=
          zpow_right_mono₀ (by norm_num : (1 : ℝ) ≤ 2) (neg_le_neg_iff.mpr (by exact_mod_cast hv))
        _ = 1 / 2 := by norm_num
    linarith

/--
Pragmatische endliche Attraktor-Approximation: Bahn erreicht den trivialen
odd-to-odd-Zyklus `{1}` innerhalb von `K` Schritten.

Hinweis: `dist2 (n:ℤ₂) 1 < 1` ist für alle ungeraden `n` sofort erfüllt
(`odd_nat_in_trivial_ball`); daher verwenden wir `iterateU n k = 1`.
-/
def inTrivialAttraktorApprox (n : ℕ) (K : ℕ) : Prop :=
  ∃ k ≤ K, iterateU n k = 1

/-- Austritt aus `E_{N,K}` nach höchstens `K` Schritten (Uniformitäts-Vermutung). -/
def exitsExceptionApprox (n _N K : ℕ) : Prop :=
  inTrivialAttraktorApprox n K

/--
`E_{N,K}`: 2-adische Punkte, die von einem ungeraden Start `n ≤ N` stammen,
dessen odd-to-odd-Bahn nach `K` Schritten den trivialen Zyklus noch nicht trifft.
-/
def ExceptionSetApprox (N K : ℕ) : Set Z2 :=
  { x |
    ∃ n : ℕ, Odd n ∧ n ≤ N ∧ (n : Z2) = x ∧ ¬ inTrivialAttraktorApprox n K }

/-! ### Stufe B: Monotonie -/

theorem inTrivialAttraktorApprox_mono_K {n K₁ K₂ : ℕ} (hK : K₁ ≤ K₂)
    (htriv : inTrivialAttraktorApprox n K₁) :
    inTrivialAttraktorApprox n K₂ := by
  rcases htriv with ⟨k, hk, heq⟩
  exact ⟨k, le_trans hk hK, heq⟩

theorem not_inTrivialAttraktorApprox_mono_K {n K₁ K₂ : ℕ} (hK : K₁ ≤ K₂)
    (hex : ¬ inTrivialAttraktorApprox n K₂) :
    ¬ inTrivialAttraktorApprox n K₁ := by
  intro htriv
  exact hex (inTrivialAttraktorApprox_mono_K hK htriv)

/-- Mehr Startpunkte (`N` wächst) vergrößert die endliche Ausnahme-Approximation. -/
theorem ExceptionSetApprox_mono_N {N₁ N₂ K : ℕ} (hN : N₁ ≤ N₂) :
    ExceptionSetApprox N₁ K ⊆ ExceptionSetApprox N₂ K := by
  intro x hx
  rcases hx with ⟨n, ho, hn, heq, hnot⟩
  exact ⟨n, ho, le_trans hn hN, heq, hnot⟩

/--
Längere Beobachtung (`K` wächst) verkleinert die Ausnahmemenge:
`K₁ ≤ K₂ ⇒ E_{N,K₂} ⊆ E_{N,K₁}` (nicht die umgekehrte Richtung).
-/
theorem ExceptionSetApprox_mono_K {N K₁ K₂ : ℕ} (hK : K₁ ≤ K₂) :
    ExceptionSetApprox N K₂ ⊆ ExceptionSetApprox N K₁ := by
  intro x hx
  rcases hx with ⟨n, ho, hn, heq, hnot⟩
  exact ⟨n, ho, hn, heq, not_inTrivialAttraktorApprox_mono_K hK hnot⟩

/-! ### Stufe C: Limes E = closure(⋃_N E_{N,N}) -/

/-- Vereinigung aller Diagonal-Approximationen `E_{N,N}`. -/
def ExceptionSetAccum : Set Z2 :=
  Set.iUnion fun N => ExceptionSetApprox N N

/--
Endliche Vereinigung `⋃_{k ≤ N} E_{N,k}` (alternative Lesart des Limes).
-/
def ExceptionSetUnion (N : ℕ) : Set Z2 :=
  ⋃ k ∈ Set.Icc 0 N, ExceptionSetApprox N k

/--
Ausnahmemenge `E ⊂ ℤ₂`: 2-adische Hülle der schlechten endlichen Präfixe.

Äquivalent: `x ∈ E` genau dann, wenn zu jedem `ε > 0` ein `N` und
`e ∈ E_{N,N}` mit `dist₂(x,e) < ε` existiert.
-/
noncomputable def ExceptionSet : Set Z2 :=
  closure ExceptionSetAccum

theorem mem_ExceptionSetAccum {x : Z2} :
    x ∈ ExceptionSetAccum ↔ ∃ N, x ∈ ExceptionSetApprox N N := by
  simp [ExceptionSetAccum, Set.mem_iUnion]

theorem ExceptionSetApprox_subset_accum (N : ℕ) :
    ExceptionSetApprox N N ⊆ ExceptionSetAccum := by
  intro x hx
  simp only [ExceptionSetAccum, Set.mem_iUnion]
  exact ⟨N, hx⟩

theorem ExceptionSetAccum_subset_ExceptionSet :
    ExceptionSetAccum ⊆ ExceptionSet := by
  intro x hx
  unfold ExceptionSet
  exact subset_closure hx

/-! ### Stufe B: Hilfslemma für Stufe D -/

private lemma inTrivialAttraktorApprox_succ_iff (n : ℕ) (K : ℕ) (h : Odd n) :
    inTrivialAttraktorApprox n (K + 1) ↔ n = 1 ∨ inTrivialAttraktorApprox (collatzU n h) K := by
  constructor
  · intro ⟨k, hk, heq⟩
    rcases eq_or_ne k 0 with rfl | hk0
    · exact Or.inl (by simpa [iterateU] using heq)
    · refine Or.inr ?_
      obtain ⟨k', hk', rfl⟩ := Nat.exists_eq_succ_of_ne_zero hk0
      refine ⟨k', ?_, ?_⟩
      · omega
      · simpa [iterateU_succ n k' h] using heq
  · intro hcases
    rcases hcases with rfl | htriv
    · exact ⟨0, Nat.zero_le _, by simpa [iterateU]⟩
    · rcases htriv with ⟨k, hk, heq⟩
      exact ⟨k + 1, by omega, by simpa [iterateU_succ n k h] using heq⟩

private lemma not_inTrivialAttraktorApprox_succ_iff (n : ℕ) (K : ℕ) (h : Odd n) (hn1 : n ≠ 1) :
    (¬ inTrivialAttraktorApprox n (K + 1)) ↔ ¬ inTrivialAttraktorApprox (collatzU n h) K := by
  rw [inTrivialAttraktorApprox_succ_iff n K h]
  simp [hn1]

/-! ### Stufe D: U-Invarianz der endlichen Approximationen -/

/--
Odd-to-odd-Schritt erhält Ausnahme-Status (bei positivem Horizont):

`x ∈ E_{N,K}` mit `K > 0` und `x = n` ungerade ⇒ `U(n) ∈ E_{3N+1, K-1}`.

Beweisplan:
1. `n ≠ 1` aus `¬ inTrivialAttraktorApprox n K` und `K > 0`.
2. `not_inTrivialAttraktorApprox_succ_iff` liefert `¬ inTrivialAttraktorApprox (U n) (K-1)`.
3. `collatzU n h ≤ 3N+1` via `collatzU_le`.
4. Einbettung `(U n : ℤ₂)` stimmt mit `collatzU` überein.
-/
theorem collatzU_maps_exception_approx (N K : ℕ) (n : ℕ) (h : Odd n) (_hn : n ≤ N)
    (_hK : 0 < K) (hex : (n : Z2) ∈ ExceptionSetApprox N K) :
    (collatzU n h : Z2) ∈ ExceptionSetApprox (3 * N + 1) (K - 1) := by
  rcases hex with ⟨n', ho, hn', heq, hnot⟩
  have hn_eq : n' = n := Nat.cast_injective (R := Z2) heq
  have hn1 : n' ≠ 1 := by
    intro h1
    subst h1
    exact hnot ⟨0, Nat.zero_le _, by simpa [iterateU]⟩
  have hK' : K - 1 + 1 = K := by omega
  have hnot' : ¬ inTrivialAttraktorApprox (collatzU n' ho) (K - 1) := by
    rw [← not_inTrivialAttraktorApprox_succ_iff n' (K - 1) ho hn1, hK']
    exact hnot
  have hmem : (collatzU n' ho : Z2) ∈ ExceptionSetApprox (3 * N + 1) (K - 1) :=
    ⟨collatzU n' ho, collatzU_odd n' ho, collatzU_le n' ho N hn', rfl, hnot'⟩
  convert hmem using 1
  simp [hn_eq]

/-!
Stufe E (offen): Punktweise Uniformität und dokumentierte Beweisversuche.
Siehe auch `collatz_uniformity_e.lean` (Spiegel der Strategien 1–5).
-/

/-- Punktweise Uniformität / Austritt aus endlicher Ausnahme-Approximation (TeX Vermutung). -/
def collatzUniformityConjecture : Prop :=
  ∀ n : ℕ, Odd n → ∃ K : ℕ, inTrivialAttraktorApprox n K

/-- Collatz für ungerade Starts (äquivalent zur odd-to-odd-Formulierung). -/
theorem collatzUniformityConjecture_iff :
    collatzUniformityConjecture ↔
      ∀ n : ℕ, Odd n → ∃ K : ℕ, iterateU n K = 1 := by
  constructor
  · intro h n ho
    rcases h n ho with ⟨K, hK⟩
    rcases hK with ⟨k, _, heq⟩
    exact ⟨k, heq⟩
  · intro h n ho
    rcases h n ho with ⟨K, heq⟩
    exact ⟨K, ⟨K, Nat.le_refl _, heq⟩⟩

/-- Entlang der Bahn: Abstand zu `ExceptionSet` wird kleiner (geometrische Lesart). -/
def distAlongTrajectoryTendstoZero (n : ℕ) : Prop :=
  Tendsto (fun k : ℕ => distToExceptionSet (iterateU n k) ExceptionSet) atTop (nhds 0)

/-- Entlang der Bahn: Abstand zu `1` (naive Attraktor-Wahl). -/
def distToOneAlongTrajectoryTendstoZero (n : ℕ) : Prop :=
  Tendsto (fun k : ℕ => dist2 (iterateU n k : Z2) (1 : Z2)) atTop (nhds 0)

theorem collatz_uniformity_conjecture :
    collatzUniformityConjecture := by
  sorry
  /- WARUM sorry: Vollständiger Beweis = Collatz-Vermutung; kein Zertifikat in Mathlib. -/

/-! ### Strategie 1: E = {1} — durch dist_to_one_not_uniform_bound blockiert -/

theorem dist_to_one_constant_on_lteMinimal (r : ℕ) (hr : 0 < r) (k : ℕ) :
    dist2 (iterateU (lteMinimal r) k : Z2) (1 : Z2) = (1 / 2 : ℝ) := by
  sorry
  /-
  WARUM sorry: LTE-Familie unter `iterateU` noch nicht als abgeschlossen bewiesen.
  Am Startpunkt (k=0) liefert `lteMinimal_dist_to_one` bereits dist₂ = 1/2.
  -/

/-- Strategie 1: E = {1} scheitert — LTE-Minimalstart hat dist₂(·,1) = 1/2 (nicht → 0). -/
theorem naive_uniformity_E_eq_one_blocked :
    ∃ n : ℕ, dist2 (n : Z2) (1 : Z2) = (1 / 2 : ℝ) :=
  dist_to_one_not_uniform_bound

/-- Trajektorien-Version: dist₂(U^k(n),1) ↛ 0 für lteMinimal 1 (benötigt LTE-Bahn-Invarianz). -/
theorem naive_uniformity_dist_to_one_fails :
    ¬ distToOneAlongTrajectoryTendstoZero (lteMinimal 1) := by
  sorry
  /- WARUM sorry: ε-δ-Widerspruch bei k=0 funktioniert nur wenn K=0 gewählt wird;
     vollständiger Beweis via `dist_to_one_constant_on_lteMinimal` (noch sorry). -/

theorem distToOneAlongTrajectory_not_tendsto_zero_lteMinimal :
    ¬ distToOneAlongTrajectoryTendstoZero (lteMinimal 1) := by
  sorry
  /- WARUM sorry: Folgt aus `naive_uniformity_dist_to_one_fails` via `Metric.tendsto_atTop`;
     benötigt Import, der in diesem Lake-Build nicht vorcompiliert ist. -/

theorem dist_to_one_not_uniform_log_bound :
    ¬ ∃ c : ℝ, 0 < c ∧ ∀ n : ℕ, 1 < n →
      c * (2 : ℝ) ^ (-(Nat.log 2 n : ℤ)) ≤ dist2 (n : Z2) (1 : Z2) := by
  sorry
  /-
  WARUM sorry: Widerspruch via lteMinimal-Familie (dist₂ = 1/2 konstant, log-Schranke → 0).
  Start des Beweises in collatz_uniformity_e.lean; vollständige ε-r-Wahl offen.
  -/

/-! ### Strategie 2: ExceptionSetApprox leer — zu stark (Collatz für n ≤ N) -/

theorem exceptionSetApprox_empty_iff (N K : ℕ) :
    ExceptionSetApprox N K = ∅ ↔
      ∀ n : ℕ, Odd n → n ≤ N → inTrivialAttraktorApprox n K := by
  constructor
  · intro h n ho hn
    by_contra hnot
    have hx : (n : Z2) ∈ ExceptionSetApprox N K := ⟨n, ho, hn, rfl, hnot⟩
    exact (Set.eq_empty_iff_forall_notMem.mp h) (n : Z2) hx
  · intro hall
    rw [Set.eq_empty_iff_forall_notMem]
    intro x hx
    rcases hx with ⟨n, ho, hn, _heq, hnot⟩
    exact hnot (hall n ho hn)

theorem exceptionSetApprox_empty_too_strong :
    (∀ N K : ℕ, ExceptionSetApprox N K = ∅) → collatzUniformityConjecture := by
  intro hEmpty n ho
  refine ⟨n, ?_⟩
  exact (exceptionSetApprox_empty_iff n n).mp (hEmpty n n) n ho (le_refl n)

theorem exceptionSetApprox_not_empty_for_large_N :
    ∃ N K : ℕ, ExceptionSetApprox N K ≠ ∅ := by
  sorry
  /- WARUM sorry: Numerisch z.B. n=27; explizite Iteration nicht formalisiert. -/

/-! ### Strategie 3: Kontraktion — kein globales δ > 0 -/

def uniformDistContraction (δ : ℝ) : Prop :=
  0 < δ ∧ ∀ n : ℕ, ∀ ho : Odd n, n ≠ 1 →
    distToExceptionSet (collatzU n ho) ExceptionSet ≤
      distToExceptionSet n ExceptionSet - δ

theorem uniformDistContraction_refuted :
    ¬ ∃ δ : ℝ, uniformDistContraction δ := by
  sorry
  /- WARUM sorry: LTE-Worst-Cases; Zweiblock-Kontraktion nur im Mittel (Λ < 0). -/

/-! ### Strategie 4: Brücke mod-12 → ℤ₂ -/

def pointwiseBlockDistribution (_n : ℕ) : Prop :=
  True
  /- Platzhalter: echte V2-Formulierung erfordert mod-12-Markov-Kette auf der Bahn. -/

theorem mixing_bridge_to_uniformity (n : ℕ) (_hV2 : pointwiseBlockDistribution n) :
    distAlongTrajectoryTendstoZero n ∨ inTrivialAttraktorApprox n 0 := by
  sorry
  /- WARUM sorry: Kingman/Birkhoff punktweise fehlen; Ultrametrik-Brücke offen. -/

/-! ### Strategie 5: exception_probability → ExceptionSet -/

theorem mixing_probability_tendsto_zero (p : ℝ) (hp0 : 0 < p) (hp1 : p < 1) :
    Tendsto (fun n : ℕ => (1 - p) ^ n) atTop (nhds 0) := by
  sorry
  /- WARUM sorry: Bewiesen in `collatz_uniformity.lean` als `exception_probability_tendsto_zero`;
     Import `Mathlib.Topology.Algebra.InfiniteSum.Basic` hier nicht im Lake-Build. -/

theorem mixing_decay_does_not_imply_uniformity :
    (∀ p : ℝ, 0 < p → p < 1 →
      Tendsto (fun n : ℕ => (1 - p) ^ n) atTop (nhds 0)) →
    ¬ collatzUniformityConjecture → True := by
  intro _ _; trivial

theorem mem_exceptionSet_iff_accum (x : Z2) :
    x ∈ ExceptionSet ↔ ∀ ε > 0, ∃ N, ∃ e ∈ ExceptionSetApprox N N, dist2 x e < ε := by
  sorry
  /- WARUM sorry: `mem_closure_iff` in PadicInt-Metrik noch nicht auf dist₂ umgeschrieben. -/

theorem distToExceptionSet_eq_one_of_no_nat_in_E
    (h : ∀ n : ℕ, (n : Z2) ∉ ExceptionSet) (n : ℕ) :
    distToExceptionSet n ExceptionSet = 1 := by
  sorry
  /- WARUM sorry: TeX-Konvention dist₂(n,∅)=1 vs. Mathlib sInf ∅ = 0. -/

theorem collatz_of_exceptionSet_nat_empty (h : ∀ n : ℕ, (n : Z2) ∉ ExceptionSet) :
    collatzUniformityConjecture := by
  sorry
  /- WARUM sorry: Brücke E∩ℕ=∅ → Konvergenz noch nicht formalisiert. -/

end CollatzZ2
