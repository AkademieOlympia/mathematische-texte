/-
  CollatzEabc.EABCCompression — Lean-Theorem-Backlog Stufe 0–4.

  Stufe 0 (Schicht **A**): 𝒞_E periodisch, Range = 4, nicht injektiv.
  Stufe 1 (Schicht **A**/**B**): H(EABC_N) ≤ 2, I_E(N) ≥ log₂ N − 2.
  Stufe 2 (Schicht **B**): D_E, Q_E, R_N, C_F — Definitionen (keine Sätze).
  Stufe 3 (Schicht **B**/**R**): structureStable — offener Satz (`Prop`, kein Beweis).
  Stufe 4: Physik (ℏ, σ), SU(2), L_i — nicht formalisiert.

  §4.8.8 / §4.8.16 / §4.8.17 in `collatz_eabc_zirkulationshypothese.md`.
-/

import CollatzEabc.Kappa
import CollatzEabc.Mod12Matrix
import Mathlib.Data.Finset.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Data.Int.CardIntervalMod
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Log.NegMulLog
import Mathlib.Algebra.BigOperators.Field
import Mathlib.Analysis.Convex.Jensen
import Mathlib.Tactic

namespace CollatzEabc

-- Interpretation: compression, not physics

/-!
## Diagnose-Box 1 — Kompression formal (Stufe 0–1)
-/

/-!
### Stufe 0 — Reine Algebra
-/

abbrev C_E (n : ℕ) : Option EabcLetter := classOfLetter n
abbrev eabcClass (n : ℕ) : Option EabcLetter := classOfLetter n

def eabcNatRep : EabcLetter → ℕ
  | 0 => 1 | 1 => 5 | 2 => 7 | 3 => 11

def eabcResidueNat (n : ℕ) : ℕ :=
  match n % 12 with
  | 1 => 1 | 5 => 5 | 7 => 7 | 11 => 11 | _ => 0

theorem C_E_on_residue (c : EabcLetter) :
    C_E (eabcNatRep c) = some c := by
  fin_cases c <;> simp [C_E, classOfLetter, eabcNatRep]

theorem mod12_period (n : ℕ) : C_E (n + 12) = C_E n := by
  simp [C_E, classOfLetter, Nat.add_mod]

theorem eabcClass_period (n : ℕ) : eabcClass (n + 12) = eabcClass n :=
  mod12_period n

theorem C_E_not_injective :
    ∃ a b : ℕ, a ≠ b ∧ C_E a = C_E b ∧ C_E a ≠ none := by
  refine ⟨1, 13, by decide, ?_, ?_⟩
  · simp [C_E, classOfLetter]
  · simp [C_E, classOfLetter]

theorem not_injective_C_E : ¬ Function.Injective (fun n : ℕ => C_E n) := by
  intro hinj
  rcases C_E_not_injective with ⟨a, b, hab, heq, _⟩
  exact hab (hinj (by simpa using heq))

theorem eabcResidue_roundtrip (c : EabcLetter) :
    eabcNatRep c % 12 = match c with
      | 0 => 1 | 1 => 5 | 2 => 7 | 3 => 11 := by
  fin_cases c <;> simp [eabcNatRep]

theorem C_E_fiber_infinite (c : EabcLetter) :
    {n : ℕ | C_E n = some c}.Infinite := by
  let r := eabcNatRep c
  have hr : C_E r = some c := C_E_on_residue c
  refine Set.infinite_of_injective_forall_mem (f := fun k => r + 12 * k) ?_ ?_
  · intro i j hij
    have : 12 * i = 12 * j := by simpa [Nat.add_mul, Nat.mul_add] using hij
    exact Nat.mul_left_cancel (by decide : (0 : ℕ) < 12) this
  · intro k
    suffices C_E (r + 12 * k) = some c from Set.mem_setOf.mpr this
    induction k with
    | zero => exact hr
    | succ k ih =>
      rw [show r + 12 * (k + 1) = (r + 12 * k) + 12 by ring]
      exact (mod12_period (r + 12 * k)).trans ih

theorem eabcClass_surjective (c : EabcLetter) :
    ∃ n, eabcClass n = some c :=
  ⟨eabcNatRep c, C_E_on_residue c⟩

theorem C_E_range_card :
    Set.ncard {c : EabcLetter | ∃ n, C_E n = some c} = 4 := by
  have hsub : {c : EabcLetter | ∃ n, C_E n = some c} = Set.univ := by
    ext c; simp [eabcClass_surjective c]
  simp [hsub, Set.ncard_univ]

theorem eabcClass_range_card : (Finset.univ : Finset EabcLetter).card = 4 := by decide

/-!
### Stufe 1 — Kompression / Entropie
-/

def eabcClassCount (S : Finset ℕ) (c : EabcLetter) : ℕ :=
  (S.filter (fun n => eabcClass n = some c)).card

noncomputable def log2 (x : ℝ) : ℝ := Real.log x / Real.log 2

noncomputable def shannonEntropyEABC (S : Finset ℕ) : ℝ :=
  let total := ∑ c : EabcLetter, eabcClassCount S c
  if total = 0 then 0 else
    -∑ c : EabcLetter,
      let count := eabcClassCount S c
      if count = 0 then (0 : ℝ) else (count / total : ℝ) * log2 (count / total)

def eabcInitialWindow (N : ℕ) : Finset ℕ :=
  (Finset.range N).map ⟨fun i => i + 1, by intro a b h; exact Nat.add_one_inj.mp h⟩

def eabcPeriodWindow (k : ℕ) : Finset ℕ := Finset.range (12 * k)

private lemma classOfLetter_eq_some_iff (n : ℕ) (c : EabcLetter) :
    classOfLetter n = some c ↔ n % 12 = eabcNatRep c := by
  fin_cases c <;> simp [classOfLetter, eabcNatRep] <;> split <;> simp

private lemma eabcClassCount_period (k : ℕ) (c : EabcLetter) :
    eabcClassCount (eabcPeriodWindow k) c = k := by
  set r := eabcNatRep c with hr
  have hr_lt : r < 12 := by fin_cases c <;> simp [hr, eabcNatRep]
  rw [eabcClassCount, eabcPeriodWindow]
  simp only [eabcClass, classOfLetter_eq_some_iff, hr]
  rw [← Nat.count_eq_card_filter_range, Nat.count_modEq_card]
  simp [hr, hr_lt, Nat.mod_eq_of_lt hr_lt, Nat.mul_mod, Nat.zero_mod]

private lemma log2_four : log2 (4 : ℝ) = 2 := by
  unfold log2
  have hpos : 0 < Real.log 2 := Real.log_pos (by norm_num : (1 : ℝ) < 2)
  have hlog4 : Real.log (4 : ℝ) = Real.log 2 * 2 := by
    rw [show (4 : ℝ) = 2 * 2 by norm_num, Real.log_mul (by norm_num) (by norm_num), two_mul]
  field_simp [hlog4, ne_of_gt hpos]

private lemma log2_inv_four : log2 (1 / 4 : ℝ) = -2 := by
  unfold log2
  rw [show (1 / 4 : ℝ) = (4 : ℝ)⁻¹ by norm_num, Real.log_inv (by norm_num : (4 : ℝ) ≠ 0), log2_four]
  field_simp

theorem shannonEntropyEABC_eq_two (k : ℕ) (hk : 0 < k) :
    shannonEntropyEABC (eabcPeriodWindow k) = 2 := by
  have hcount : ∀ c, eabcClassCount (eabcPeriodWindow k) c = k := fun c => eabcClassCount_period k c
  have htotal : ∑ c, eabcClassCount (eabcPeriodWindow k) c = 4 * k := by
    simp [Finset.sum_const, Finset.card_fin, hcount]
  have htotal_ne : (∑ c, eabcClassCount (eabcPeriodWindow k) c) ≠ 0 := by omega
  have hp : ∀ c, (eabcClassCount (eabcPeriodWindow k) c : ℝ) / (4 * k) = 1 / 4 := by
    intro c; field_simp [hcount, hk]
  have hbody :
      -∑ c, (eabcClassCount (eabcPeriodWindow k) c : ℝ) / (4 * k) *
        log2 ((eabcClassCount (eabcPeriodWindow k) c : ℝ) / (4 * k)) = 2 := by
    calc
      _ = -∑ c, (1 / 4 : ℝ) * log2 (1 / 4) := by simp [hp, Finset.sum_const, Finset.card_fin]
      _ = -4 * ((1 / 4 : ℝ) * log2 (1 / 4)) := by simp [Finset.sum_const, Finset.card_fin]
      _ = 2 := by rw [log2_inv_four]; norm_num
  simpa [shannonEntropyEABC, dif_neg htotal_ne, htotal] using hbody

private lemma entropy_prob_le_log_four
    (p : EabcLetter → ℝ) (hp_nonneg : ∀ c, 0 ≤ p c) (hp_sum : ∑ c, p c = 1) :
    -∑ c, p c * Real.log (p c) ≤ Real.log 4 := by
  have hp_mem : ∀ c, p c ∈ Set.Ici (0 : ℝ) := fun c => Set.mem_Ici.mpr (hp_nonneg c)
  have hJ :=
    Real.concaveOn_negMulLog.le_map_sum (t := Finset.univ) (w := fun _ => (4 : ℝ)⁻¹) (p := p)
      (by intro c _; positivity) (by simp [Finset.sum_const, Finset.card_fin]) (fun c _ => hp_mem c)
  have hleft :
      (4 : ℝ)⁻¹ * ∑ c, Real.negMulLog (p c) ≤ Real.negMulLog (1 / 4 : ℝ) := by
    simpa [Finset.sum_const, Finset.card_fin, hp_sum, smul_eq_mul, one_div] using hJ
  have hsum : ∑ c, Real.negMulLog (p c) = -∑ c, p c * Real.log (p c) := by
    simp [Real.negMulLog, Finset.sum_neg_distrib, mul_comm]
  have hright : Real.negMulLog (1 / 4 : ℝ) = (1 / 4 : ℝ) * Real.log 4 := by
    simp [Real.negMulLog, Real.log_div (by norm_num) (by norm_num : (4 : ℝ) ≠ 0), Real.log_one]
  rw [hsum] at hleft; rw [hright] at hleft; linarith

theorem shannonEntropyEABC_le_two (S : Finset ℕ) : shannonEntropyEABC S ≤ 2 := by
  by_cases h0 : (∑ c, eabcClassCount S c) = 0
  · simp [shannonEntropyEABC, h0]
  · have htotal_ne : (∑ c, eabcClassCount S c) ≠ 0 := h0
    set total := ∑ c, eabcClassCount S c with htotal_def
    have hne : (total : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr htotal_ne
    have hprob : ∑ c, (eabcClassCount S c : ℝ) / total = 1 := by
      have hcast : (∑ c, (eabcClassCount S c : ℝ)) = total := by simp [htotal_def]
      rw [Finset.sum_div, hcast, div_self hne]
    have hle_nat :
        -∑ c, (eabcClassCount S c : ℝ) / total * Real.log ((eabcClassCount S c : ℝ) / total) ≤
          Real.log 4 := by
      refine entropy_prob_le_log_four (fun c => (eabcClassCount S c : ℝ) / total) ?_ hprob
      intro c; exact div_nonneg (Nat.cast_nonneg _) (Nat.cast_nonneg _)
    have hlog2pos : 0 < Real.log 2 := Real.log_pos (by norm_num : (1 : ℝ) < 2)
    have hterm (c : EabcLetter) :
        (if eabcClassCount S c = 0 then (0 : ℝ) else
          (eabcClassCount S c : ℝ) / total *
            (Real.log ((eabcClassCount S c : ℝ) / total) / Real.log 2)) =
          (eabcClassCount S c : ℝ) / total *
            Real.log ((eabcClassCount S c : ℝ) / total) / Real.log 2 := by
      by_cases hc : eabcClassCount S c = 0 <;> simp [hc, zero_div, mul_div_assoc]
    have hsum :
        ∑ c, if eabcClassCount S c = 0 then (0 : ℝ) else
          (eabcClassCount S c : ℝ) / total *
            (Real.log ((eabcClassCount S c : ℝ) / total) / Real.log 2) =
          ∑ c, (eabcClassCount S c : ℝ) / total *
            Real.log ((eabcClassCount S c : ℝ) / total) / Real.log 2 := by
      refine Finset.sum_congr rfl fun c _ => ?_
      exact hterm c
    have hbound :
        -∑ c, (eabcClassCount S c : ℝ) / total *
            Real.log ((eabcClassCount S c : ℝ) / total) / Real.log 2 ≤ 2 := by
      have := (div_le_div_iff₀ hlog2pos hlog2pos).mpr hle_nat
      rwa [log2_four]
    have hval :
        shannonEntropyEABC S =
          -∑ c, (eabcClassCount S c : ℝ) / total *
            Real.log ((eabcClassCount S c : ℝ) / total) / Real.log 2 := by
      simp [shannonEntropyEABC, dif_neg htotal_ne, log2, htotal_def, hsum]
    rw [hval]
    exact hbound

noncomputable def uniformEntropy (N : ℕ) : ℝ :=
  if N = 0 then 0 else log2 N

theorem uniformEntropy_eq_log2 (N : ℕ) (hN : 0 < N) :
    uniformEntropy N = log2 N := by
  simp [uniformEntropy, Nat.ne_of_gt hN]

noncomputable def I_E (N : ℕ) : ℝ :=
  uniformEntropy N - shannonEntropyEABC (eabcInitialWindow N)

theorem I_E_ge_log2_sub_two (N : ℕ) (hN : 0 < N) :
    I_E N ≥ log2 N - 2 := by
  unfold I_E
  rw [uniformEntropy_eq_log2 N hN]
  linarith [shannonEntropyEABC_le_two (eabcInitialWindow N)]

/-!
## Diagnose-Box 2 — Strukturstabilität offen (Stufe 2–3)
-/

structure EABCMeasures where
  D_E : ℤ
  Q_E : ℤ
  R_N : ℤ
  C_F : ℤ

def eabcMeasureDist (M₁ M₂ : EABCMeasures) : ℝ :=
  |((M₁.D_E - M₂.D_E) : ℤ)| + |((M₁.Q_E - M₂.Q_E) : ℤ)| +
    |((M₁.R_N - M₂.R_N) : ℤ)| + |((M₁.C_F - M₂.C_F) : ℤ)|

def structureStable (_C : ℕ → Option EabcLetter) (M : EABCMeasures → ℝ) (ε : ℝ) : Prop :=
  ∀ M_raw M_comp : EABCMeasures,
    eabcMeasureDist M_raw M_comp ≤ ε → |M M_raw - M M_comp| ≤ ε

structure StructureStability where
  preservesObservableStructure : Prop

/-- Offener Satz (Stufe 3): $\mathcal{C}_E$ strukturstabil bezüglich $M$ — kein Beweis, nur `Prop`. -/
def C_E_structureStable_conjecture : Prop :=
  structureStable C_E (fun m => (m.D_E : ℝ)) 0

def structureStable_open : Prop := C_E_structureStable_conjecture

/-!
## Diagnose-Box 3 — Physik Hypothese (Stufe 4)

-- Stufe 4: not formalized
-- (ℏ, σ), SU(2), L_i — fehlende Konstruktion (Schicht C/R).
-/

theorem C_E_extreme_compression :
    (¬ Function.Injective (fun n : ℕ => C_E n)) ∧
    (∀ c : EabcLetter, {n : ℕ | C_E n = some c}.Infinite) ∧
    (∀ S : Finset ℕ, shannonEntropyEABC S ≤ 2) := by
  refine ⟨not_injective_C_E, ?_, ?_⟩
  · intro c; exact C_E_fiber_infinite c
  · intro S; exact shannonEntropyEABC_le_two S

end CollatzEabc
