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
      (eabcClassCount S c : ℝ) / total *
        Real.log ((eabcClassCount S c : ℝ) / total) / Real.log 2

def eabcInitialWindow (N : ℕ) : Finset ℕ :=
  (Finset.range N).map ⟨fun i => i + 1, by intro a b h; exact Nat.add_one_inj.mp h⟩

def eabcPeriodWindow (k : ℕ) : Finset ℕ := Finset.range (12 * k)

private lemma log2_four : log2 (4 : ℝ) = 2 := by
  unfold log2
  calc Real.log 4 / Real.log 2
      = Real.log (2 ^ (2 : ℕ)) / Real.log 2 := by rw [show (4 : ℝ) = 2 ^ 2 by norm_num]
    _ = (2 * Real.log 2) / Real.log 2 := by rw [Real.log_pow]; ring_nf
    _ = 2 := by field_simp

private lemma entropy_prob_le_log_four
    (p : EabcLetter → ℝ) (hp_nonneg : ∀ c, 0 ≤ p c) (hp_sum : ∑ c, p c = 1) :
    -∑ c, p c * Real.log (p c) ≤ Real.log 4 := by
  have hp_mem : ∀ c, p c ∈ Set.Ici (0 : ℝ) := fun c => Set.mem_Ici.mpr (hp_nonneg c)
  have hJ :=
    ConcaveOn.le_map_sum (𝕜 := ℝ) (β := ℝ) Real.concaveOn_negMulLog
      (t := Finset.univ) (w := fun _ => (4 : ℝ)⁻¹) (p := p)
      (h₀ := by intro c _; positivity)
      (h₁ := by simp [Finset.sum_const, Finset.card_fin])
      (hmem := fun c _ => hp_mem c)
  have hJ' : ∑ c, (4 : ℝ)⁻¹ * Real.negMulLog (p c) ≤
      Real.negMulLog (∑ c, (4 : ℝ)⁻¹ * p c) := by simpa [smul_eq_mul] using hJ
  have hsumw : ∑ c, (4 : ℝ)⁻¹ * p c = (4 : ℝ)⁻¹ * ∑ c, p c := by
    rw [← Finset.mul_sum]
  have hcent : (4 : ℝ)⁻¹ * ∑ c, p c = (1 / 4 : ℝ) := by rw [hp_sum]; norm_num
  have hleft : (4 : ℝ)⁻¹ * ∑ c, Real.negMulLog (p c) ≤ Real.negMulLog (1 / 4 : ℝ) := by
    rw [← Finset.mul_sum, hsumw, hcent] at hJ'
    exact hJ'
  have hsum : ∑ c, Real.negMulLog (p c) = -∑ c, p c * Real.log (p c) := by
    simp [Real.negMulLog, Finset.sum_neg_distrib, mul_comm]
  have hright : Real.negMulLog (1 / 4 : ℝ) = (4 : ℝ)⁻¹ * Real.log 4 := by
    unfold Real.negMulLog
    rw [show (1 / 4 : ℝ) = (4 : ℝ)⁻¹ by norm_num, Real.log_inv (4 : ℝ)]
    ring
  rw [hsum, hright] at hleft
  exact le_of_mul_le_mul_left hleft (by norm_num : 0 < (4 : ℝ)⁻¹)

theorem shannonEntropyEABC_le_two (S : Finset ℕ) : shannonEntropyEABC S ≤ 2 := by
  by_cases h0 : (∑ c, eabcClassCount S c) = 0
  · simp [shannonEntropyEABC, h0]
  · have htotal_ne : (∑ c, eabcClassCount S c) ≠ 0 := h0
    set total := ∑ c, eabcClassCount S c with htotal_def
    have hne : (total : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr htotal_ne
    have hprob : ∑ c, (eabcClassCount S c : ℝ) / total = 1 := by
      have hcast : (∑ c, (eabcClassCount S c : ℝ)) = total := by simp [htotal_def]
      rw [← Finset.sum_div, hcast, div_self hne]
    have hle_nat :
        -∑ c, (eabcClassCount S c : ℝ) / total * Real.log ((eabcClassCount S c : ℝ) / total) ≤
          Real.log 4 := by
      refine entropy_prob_le_log_four (fun c => (eabcClassCount S c : ℝ) / total) ?_ hprob
      intro c; exact div_nonneg (Nat.cast_nonneg _) (Nat.cast_nonneg _)
    have hlog2pos : 0 < Real.log 2 := Real.log_pos (by norm_num : (1 : ℝ) < 2)
    have hfactor :
        ∑ c, (eabcClassCount S c : ℝ) / total *
          Real.log ((eabcClassCount S c : ℝ) / total) / Real.log 2 =
        (∑ c, (eabcClassCount S c : ℝ) / total *
          Real.log ((eabcClassCount S c : ℝ) / total)) / Real.log 2 := by
      simpa using (Finset.sum_div (s := (Finset.univ : Finset EabcLetter))
        (fun c => (eabcClassCount S c : ℝ) / total * Real.log ((eabcClassCount S c : ℝ) / total))
        (Real.log 2)).symm
    have hcore :
        -∑ c, (eabcClassCount S c : ℝ) / total *
            Real.log ((eabcClassCount S c : ℝ) / total) / Real.log 2 ≤ log2 (4 : ℝ) := by
      rw [hfactor, ← neg_div]
      simpa [log2] using div_le_div_of_nonneg_right hle_nat hlog2pos.le
    have hbound : -∑ c, (eabcClassCount S c : ℝ) / total *
        Real.log ((eabcClassCount S c : ℝ) / total) / Real.log 2 ≤ 2 := by
      simpa [log2_four] using hcore
    have hshannon :
        shannonEntropyEABC S = -∑ c, (eabcClassCount S c : ℝ) / total *
          Real.log ((eabcClassCount S c : ℝ) / total) / Real.log 2 := by
      unfold shannonEntropyEABC
      by_cases ht : total = 0
      · exfalso; exact htotal_ne ht
      · rw [show (∑ c, eabcClassCount S c) = total from htotal_def.symm]
        simp [ht]
    rw [hshannon]
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
