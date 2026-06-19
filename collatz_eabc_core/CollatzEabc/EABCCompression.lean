/-
  CollatzEabc.EABCCompression — mod-12 EABC-Kompressionsoperator 𝒞_E.

  Schicht **A** (Kompression): periodische, nicht injektive Abbildung ℕ → {E,A,B,C}
  (Reste 1,5,7,11 modulo 12). Verknüpft mit §4.8.8 / §4.8.16 in
  `collatz_eabc_zirkulationshypothese.md`.

  Schichtung: Kompression → Observablen → numerische Tests → Interpretation
  (vgl. `collatz_eabc_epistemik_schichten.md`).

  Referenz: `Kappa.classOfLetter`, `Mod12Matrix.eabcResidue`.
-/

import CollatzEabc.Kappa
import CollatzEabc.Mod12Matrix
import Mathlib.Data.Finset.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic

namespace CollatzEabc

-- Interpretation: compression, not physics

/-!
### Box: periodische nicht-injektive Projektion

$$\boxed{\;\mathcal{C}_E(n)=n\bmod 12\in\{1,5,7,11\}\cong\{E,A,B,C\}\;:\;
\text{12-periodisch, nicht injektiv, Bildgröße }4.\;}$$
-/

/-- EABC-Kompressionsoperator 𝒞_E : ℕ → Option EabcLetter (nur Reste 1,5,7,11). -/
abbrev C_E (n : ℕ) : Option EabcLetter :=
  classOfLetter n

/-- Komprimierte Klasse als `Fin 4` (E=0, A=1, B=2, C=3); identisch mit `classOfLetter`. -/
abbrev eabcClass (n : ℕ) : Option EabcLetter :=
  classOfLetter n

/-- Mod-12-Rest in `{1,5,7,11}` als natürliche Zahl; sonst `0` (kein EABC-Rest). -/
def eabcResidueNat (n : ℕ) : ℕ :=
  match n % 12 with
  | 1 => 1
  | 5 => 5
  | 7 => 7
  | 11 => 11
  | _ => 0

private lemma eabcResidue_val (c : EabcLetter) :
    (eabcResidue c).val = match c with
      | 0 => 1 | 1 => 5 | 2 => 7 | 3 => 11 := by
  fin_cases c <;> simp [eabcResidue]

/-- 𝒞_E ist auf allen Restklassen mod 12 wohldefiniert (12-Periodizität). -/
theorem mod12_period (n : ℕ) :
    C_E (n + 12) = C_E n := by
  simp [C_E, classOfLetter, Nat.add_mod]

theorem eabcClass_period (n : ℕ) :
    eabcClass (n + 12) = eabcClass n :=
  mod12_period n

/-- 𝒞_E ist nicht injektiv: unendlich viele Zahlen teilen dieselbe EABC-Klasse. -/
theorem C_E_not_injective :
    ∃ a b : ℕ, a ≠ b ∧ C_E a = C_E b ∧ C_E a ≠ none := by
  refine ⟨1, 13, by decide, ?_, ?_⟩
  · simp [C_E, classOfLetter]
  · simp [C_E, classOfLetter]

/-- Alias: explizite Negation der Injektivität auf dem EABC-Bild. -/
theorem not_injective_C_E : ¬ Function.Injective (fun n : ℕ => C_E n) := by
  intro hinj
  rcases C_E_not_injective with ⟨a, b, hab, heq, _⟩
  exact hab (hinj (by simpa using heq))

/-- Jede EABC-Klasse hat unendlich viele Urbilder (extrem verlustreiche Kompression). -/
theorem C_E_fiber_infinite (c : EabcLetter) :
    {n : ℕ | C_E n = some c}.Infinite := by
  let r : ℕ := (eabcResidue c).val
  have hr : C_E r = some c := C_E_on_residue c
  refine Set.infinite_of_injective_forall_mem (f := fun k => r + 12 * k) ?_ ?_
  · intro i j hij
    have : 12 * i = 12 * j := by simpa [Nat.add_mul, Nat.mul_add] using hij
    simpa using Nat.mul_left_cancel (by decide : (0 : ℕ) < 12) this
  · intro k
    simpa [Set.mem_setOf_eq, mod12_period] using hr

/-- Restklassen-Bijektion zwischen `EabcLetter` und den vier EABC-Resten mod 12. -/
theorem eabcResidue_roundtrip (c : EabcLetter) :
    (eabcResidue c).val % 12 = match c with
      | 0 => 1 | 1 => 5 | 2 => 7 | 3 => 11 := by
  fin_cases c <;> simp [eabcResidue]

theorem C_E_on_residue (c : EabcLetter) :
    C_E (eabcResidue c |>.val) = some c := by
  fin_cases c <;> simp [C_E, classOfLetter, eabcResidue, ZMod.val_natCast]

/-- Surjektivität auf die vier EABC-Klassen: `#Range(𝒞_E|_{E,A,B,C}) = 4`. -/
theorem eabcClass_surjective (c : EabcLetter) :
    ∃ n, eabcClass n = some c :=
  ⟨(eabcResidue c).val, C_E_on_residue c⟩

/-- Kardinalität des EABC-Bildes ist exakt vier. -/
theorem C_E_range_card :
    Set.ncard {c : EabcLetter | ∃ n, C_E n = some c} = 4 := by
  have hsub : {c : EabcLetter | ∃ n, C_E n = some c} = Set.univ := by
    ext c
    simp only [Set.mem_setOf_eq, Set.mem_univ, iff_true]
    exact eabcClass_surjective c
  simp [hsub, Set.ncard_univ]

theorem eabcClass_range_card :
    (Finset.univ : Finset EabcLetter).card = 4 := by
  decide

/-!
### Shannon-Entropie endlicher Fenster (Schicht A)

Empirische Verteilung der vier EABC-Klassen in einem endlichen Fenster `S`.
Maximum bei Gleichverteilung: $H = \log_2 4 = 2$ Bit.
-/

/-- Häufigkeit der EABC-Klasse `c` im endlichen Fenster `S`. -/
def eabcClassCount (S : Finset ℕ) (c : EabcLetter) : ℕ :=
  (S.filter (fun n => eabcClass n = some c)).card

private def log2 (x : ℝ) : ℝ :=
  Real.log x / Real.log 2

/-- Shannon-Entropie (Bit) der vier EABC-Klassen im Fenster `S` (nur EABC-Treffer). -/
noncomputable def shannonEntropyEABC (S : Finset ℕ) : ℝ :=
  let total := ∑ c : EabcLetter, eabcClassCount S c
  if h0 : total = 0 then
    0
  else
    -∑ c : EabcLetter,
      let count := eabcClassCount S c
      if hc : count = 0 then
        0
      else
        let p : ℝ := count / total
        p * log2 p

/-- Gleichverteiltes Fenster: je `k` Treffer pro EABC-Klasse (`4k` Zahlen). -/
def eabcUniformWindow (k : ℕ) : Finset ℕ :=
  (Finset.univ : Finset EabcLetter).biUnion fun c =>
    (Finset.range k).map
      ⟨fun i => (eabcResidue c).val + 12 * i, by
        intro i j hij
        have : 12 * i = 12 * j := by
          have := congrArg (fun t => t - (eabcResidue c).val) hij
          simpa [Nat.add_sub_cancel_left] using this
        exact Nat.mul_left_cancel (by decide : (0 : ℕ) < 12) this⟩

private lemma eabcClassCount_uniform (k : ℕ) (c : EabcLetter) :
    eabcClassCount (eabcUniformWindow k) c = k := by
  simp [eabcClassCount, eabcUniformWindow, C_E_on_residue, Finset.card_map,
    Finset.card_range, Finset.mem_biUnion, Finset.mem_univ, true_and]

private lemma log2_four : log2 (4 : ℝ) = 2 := by
  have h2 : (0 : ℝ) < 2 := by norm_num
  have h4 : (4 : ℝ) = 2 ^ 2 := by norm_num
  simp [log2, h4, Real.log_pow, Real.log_two]

private lemma log2_inv_four : log2 (1 / 4 : ℝ) = -2 := by
  have hpos : (0 : ℝ) < 4 := by norm_num
  have hne : (4 : ℝ) ≠ 0 := by norm_num
  rw [log2, one_div, Real.log_inv hne, log2_four]
  field_simp

/-- Gleichverteilung auf vier Klassen: Entropie exakt 2 Bit. -/
theorem shannonEntropyEABC_eq_two (k : ℕ) (hk : 0 < k) :
    shannonEntropyEABC (eabcUniformWindow k) = 2 := by
  have hcount : ∀ c : EabcLetter, eabcClassCount (eabcUniformWindow k) c = k :=
    fun c => eabcClassCount_uniform k c
  have htotal : ∑ c : EabcLetter, eabcClassCount (eabcUniformWindow k) c = 4 * k := by
    simp [Finset.sum_const, Finset.card_fin, hcount]
  have htotal_ne : (∑ c : EabcLetter, eabcClassCount (eabcUniformWindow k) c) ≠ 0 := by
    omega
  simp only [shannonEntropyEABC, dif_neg htotal_ne, htotal]
  have hp : ∀ c : EabcLetter,
      (eabcClassCount (eabcUniformWindow k) c : ℝ) / (4 * k) = 1 / 4 := by
    intro c
    simp [hcount, hk]
  simp [Finset.sum_const, Finset.card_fin, hp, log2_inv_four]

/-- Shannon-Entropie über höchstens vier Klassen: $H(\mathrm{EABC}_N) \le 2$ Bit. -/
theorem shannonEntropyEABC_le_two (S : Finset ℕ) :
    shannonEntropyEABC S ≤ 2 := by
  by_cases h0 : (∑ c : EabcLetter, eabcClassCount S c) = 0
  · simp [shannonEntropyEABC, h0]
  · simp only [shannonEntropyEABC, dif_neg h0]
    have htotal_pos : 0 < ∑ c : EabcLetter, eabcClassCount S c := by omega
    have hbound :
        ∀ c : EabcLetter,
          -((eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c' *
              log2 ((eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c')) ≤
            (eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c' * 2 := by
      intro c
      by_cases hc : eabcClassCount S c = 0
      · simp [hc]
      · have hp_pos : 0 < (eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c' := by
          apply div_pos <;> exact_mod_cast Nat.pos_of_ne_zero hc
        have hp_le : (eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c' ≤ 1 := by
          rw [div_le_one (by exact_mod_cast htotal_pos)]
          exact_mod_cast Finset.single_le_sum (fun _ _ => Nat.zero_le _) (Finset.mem_univ c)
        have hlog_le :
            log2 ((eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c') ≤ 2 := by
          rw [log2_four]
          apply log2_le_log2_of_le hp_pos
          linarith
        have hterm_nonpos :
            (eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c' *
              log2 ((eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c') ≤ 0 := by
          have hlog_nonpos :
              log2 ((eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c') ≤ 0 := by
            rw [log2, div_nonpos_iff]
            right
            constructor
            · exact Real.log_nonpos (le_of_lt hp_pos) hp_le
            · exact Real.log_pos (by norm_num : (1 : ℝ) < 2)
          nlinarith
        linarith
    calc
      -∑ c : EabcLetter,
          (eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c' *
            log2 ((eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c')
          = ∑ c : EabcLetter,
              -((eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c' *
                log2 ((eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c')) := by
            simp [Finset.sum_neg_distrib]
      _ ≤ ∑ c : EabcLetter,
            (eabcClassCount S c : ℝ) / ∑ c' : EabcLetter, eabcClassCount S c' * 2 := by
            exact Finset.sum_le_sum fun c _ => hbound c
      _ = 2 := by
            simp [Finset.sum_mul, div_mul_eq_mul_div, mul_comm]
            field_simp

private lemma log2_le_log2_of_le {a b : ℝ} (ha : 0 < a) (hab : a ≤ b) :
    log2 a ≤ log2 b := by
  unfold log2
  exact div_le_div_of_le_left (Real.log_le_log ha hab) (Real.log_pos (by norm_num : (1 : ℝ) < 2))
    (Real.log_pos (by norm_num : (1 : ℝ) < 2))

/-!
### Strukturstabilität (Schicht R — offen)

Qualitative Erhaltung von Observablen unter 𝒞_E; numerisch via `eabc_compression_robustness.py`.
-/

/-- Daten einer strukturstabilen Kompression bzgl. Prim-Observablen (OFFEN). -/
structure StructureStability where
  /-- Observablen-Metrik bleibt qualitativ erhalten nach Anwendung von 𝒞_E. -/
  preservesObservableStructure : Prop

/-- **Offene Vermutung:** EABC-Kompression ist strukturstabil bzgl. Primstruktur (§4.8.8). -/
def structureStable : Prop :=
  Nonempty StructureStability

/-- Dokumentierter offener Satz — Nachweis via Robustheitstest §4.8.9. -/
theorem structureStable_open : structureStable := by
  sorry

end CollatzEabc
