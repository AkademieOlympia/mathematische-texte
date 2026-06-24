/-
  CollatzEabc.HolonomieFehlerterm — Phase 1–2: D_E, Lückenmuster, Bell/CHSH-Skeleton.

  **Nicht sorry-frei** (Prime-Enumeration, asymptotische Hol_E=0, vollständiger CHSH-Beweis offen).
  Kombinatorische Teile (Lückenmuster, Taubenloch auf Fenstern) sind bewiesen.

  Referenz:
    collatz_eabc_fehlerterm_hypothese.md
    collatz_eabc_bell_holonomie.md
    collatz_eabc_holonomie_fehlerterm.py
    collatz_eabc_bell_inequality_test.py
-/

import Mathlib.Data.Fin.Basic
import Mathlib.Data.List.Basic
import Mathlib.Data.Int.Basic
import Mathlib.Data.Rat.Defs
import Mathlib.Topology.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic
import CollatzEabc.Mod12Matrix
import CollatzEabc.PatternCount
import CollatzEabc.PrefProjection

namespace CollatzEabc

open List

/-!
### EABC-Wörter und mod-12-Lücken (bewiesen)

`EabcLetter` = Fin 4 mit E=0, A=1, B=2, C=3 (vgl. `PrefProjection`, `Mod12Matrix`).
-/

/-- Restklasse mod 12 der EABC-Klasse: E↦1, A↦5, B↦7, C↦11. -/
def eabcResidueNat : EabcLetter → ℕ
  | 0 => 1
  | 1 => 5
  | 2 => 7
  | 3 => 11

/-- Kanonisches Lückenmuster (2,4,2,4) entlang geschlossener 5-Zyklen. -/
def canonicalGapPattern : List ℕ := [2, 4, 2, 4]

/-- Offenes 4-Fenster ABCE (Bell-Taubenloch-Träger). -/
def wordABCE : List EabcLetter := [⟨0, by decide⟩, ⟨1, by decide⟩, ⟨2, by decide⟩,
  ⟨3, by decide⟩]

/-- Lücke (r₂ − r₁) mod 12. -/
def gapMod12 (r₁ r₂ : ℕ) : ℕ := (r₂ + 12 - r₁) % 12

/-- Lückenmuster entlang einer Restklassenfolge. -/
def gapPattern (residues : List ℕ) : List ℕ :=
  match residues with
  | [] => []
  | _ :: rs =>
    (residues.zip rs).map fun (a, b) => gapMod12 a b

/-- Restklassenfolge eines EABC-Worts. -/
def residuesOf (w : List EabcLetter) : List ℕ :=
  w.map eabcResidueNat

/-- Zählt Vorkommen von `patternWord` als Gleitfenster in `classes`. -/
def countSlidingWord (classes : List EabcLetter) (patternWord : List EabcLetter) : ℕ :=
  let k := patternWord.length
  if _h : classes.length < k then 0
  else
    let windows :=
      (List.range (classes.length + 1 - k)).filter fun i =>
        List.take k (classes.drop i) = patternWord
    windows.length

/-- N₊(W) = #{ABCEA-Fenster in endlicher Klassenfolge W}. -/
def N_plus (classes : List EabcLetter) : ℕ :=
  countSlidingWord classes wordABCEA

/-- N₋(W) = #{CEABC-Fenster in endlicher Klassenfolge W}. -/
def N_minus (classes : List EabcLetter) : ℕ :=
  countSlidingWord classes wordCEABC

/-- Fehlerterm D_E(W) = N₊ − N₋ auf endlicher Folge. -/
def D_E_on (classes : List EabcLetter) : ℤ :=
  (N_plus classes : ℤ) - N_minus classes

/-- Paar (D_E, N₊+N₋) für normalisierte Auswertung in Python (`D̃_E = D_E/√total`). -/
def D_E_pair (classes : List EabcLetter) : ℤ × ℕ :=
  (D_E_on classes, N_plus classes + N_minus classes)

/-- Q_E(W) = N₊ + N₋ (Gesamtzählung orientierter 5-Zyklen). -/
def Q_E_on (classes : List EabcLetter) : ℕ :=
  N_plus classes + N_minus classes

/-- χ_Hol(W) = D_E / Q_E = W_E auf endlicher Folge. -/
def chi_Hol (classes : List EabcLetter) : ℚ :=
  let total := Q_E_on classes
  if _h : total = 0 then 0
  else (D_E_on classes : ℚ) / total

/-- R_β(W) = D_E / Q_E^β auf endlicher Folge (β ∈ ℝ; β = 1/2 ⇒ D̃_E). -/
noncomputable def R_beta (classes : List EabcLetter) (β : ℝ) : ℝ :=
  let q := Q_E_on classes
  if _h : q = 0 then 0
  else (D_E_on classes : ℝ) / (q : ℝ) ^ β

/-- D̃_E(W) = D_E / √Q_E — Belastungstest-Observable R_{1/2}. -/
noncomputable def D_tilde_E (classes : List EabcLetter) : ℝ :=
  R_beta classes (1 / 2)

/-!
### Bewiesene kombinatorische Lemmata
-/

/-- **[Schicht A]** −1 ≤ χ_Hol ≤ 1 auf jeder endlichen Folge. -/
theorem chi_Hol_bounds (classes : List EabcLetter) :
    -1 ≤ chi_Hol classes ∧ chi_Hol classes ≤ 1 := by
  unfold chi_Hol Q_E_on D_E_on
  by_cases h : N_plus classes + N_minus classes = 0
  · rw [dif_pos h]
    norm_num
  ·
    have hpos_nat : 0 < N_plus classes + N_minus classes := Nat.pos_of_ne_zero h
    have hpos_rat : (0 : ℚ) < ((N_plus classes + N_minus classes : ℕ) : ℚ) := by
      exact_mod_cast hpos_nat
    rw [dif_neg h]
    constructor
    ·
      have hineq_int :
          -((N_plus classes + N_minus classes : ℕ) : Int)
            ≤ (N_plus classes : Int) - (N_minus classes : Int) := by omega
      have hineq_rat :
          (-(((N_plus classes + N_minus classes : ℕ) : Int) : ℚ))
            ≤ (((N_plus classes : Int) - (N_minus classes : Int) : Int) : ℚ) := by
        exact_mod_cast hineq_int
      have hdiv := div_le_div_of_nonneg_right hineq_rat (le_of_lt hpos_rat)
      have hneg :
          (-(((N_plus classes + N_minus classes : ℕ) : Int) : ℚ)) /
              ((N_plus classes + N_minus classes : ℕ) : ℚ) = -1 := by
        rw [show (-(((N_plus classes + N_minus classes : ℕ) : Int) : ℚ))
              = -((N_plus classes + N_minus classes : ℕ) : ℚ) from by push_cast; rfl]
        field_simp [ne_of_gt hpos_rat]
      rw [← hneg]
      exact hdiv
    ·
      have hineq_int :
          (N_plus classes : Int) - (N_minus classes : Int)
            ≤ ((N_plus classes + N_minus classes : ℕ) : Int) := by omega
      have hineq_rat :
          (((N_plus classes : Int) - (N_minus classes : Int) : Int) : ℚ)
            ≤ (((N_plus classes + N_minus classes : ℕ) : Int) : ℚ) := by
        exact_mod_cast hineq_int
      have hdiv := div_le_div_of_nonneg_right hineq_rat (le_of_lt hpos_rat)
      have hone :
          (((N_plus classes + N_minus classes : ℕ) : Int) : ℚ) /
              ((N_plus classes + N_minus classes : ℕ) : ℚ) = 1 := by
        rw [show (((N_plus classes + N_minus classes : ℕ) : Int) : ℚ)
              = ((N_plus classes + N_minus classes : ℕ) : ℚ) from by push_cast; rfl]
        field_simp [ne_of_gt hpos_rat]
      rw [← hone]
      exact hdiv

/-- **[Schicht A]** N₊ = N₋ ⇒ χ_Hol = 0. -/
theorem chi_Hol_zero_of_balance (classes : List EabcLetter)
    (h : N_plus classes = N_minus classes) : chi_Hol classes = 0 := by
  unfold chi_Hol Q_E_on D_E_on
  simp [h]

/-- **[Schicht A]** Q_E_on = N₊ + N₋ (definitorisch). -/
theorem Q_E_on_eq_sum (classes : List EabcLetter) :
    Q_E_on classes = N_plus classes + N_minus classes := rfl

/-- **[Schicht A]** D_E_on = N₊ − N₋ (definitorisch). -/
theorem D_E_on_eq_diff (classes : List EabcLetter) :
    D_E_on classes = (N_plus classes : ℤ) - N_minus classes := rfl

theorem abcea_residues : residuesOf wordABCEA = [5, 7, 11, 1, 5] := rfl

theorem ceabc_residues : residuesOf wordCEABC = [11, 1, 5, 7, 11] := rfl

theorem abcea_gap_pattern : gapPattern (residuesOf wordABCEA) = canonicalGapPattern := rfl

theorem ceabc_gap_pattern : gapPattern (residuesOf wordCEABC) = canonicalGapPattern := rfl

theorem holonomy_words_share_gap_pattern :
    gapPattern (residuesOf wordABCEA) = gapPattern (residuesOf wordCEABC) := by
  rw [abcea_gap_pattern, ceabc_gap_pattern]

theorem abcea_ceabc_only_start_differs :
    (residuesOf wordABCEA).head? = some 5 ∧
    (residuesOf wordCEABC).head? = some 11 := by
  constructor <;> rfl

/-!
### Taubenloch auf ABCE-Fenstern (Theorem, Fensterweise)
-/

/-- Drei Bits auf einem ABCE-Fenster: mindestens ein gleiches Paar (Taubenloch). -/
theorem pigeonhole_three_bits (a b c : Bool) :
    (a = b) ∨ (a = c) ∨ (b = c) := by
  rcases a <;> rcases b <;> rcases c <;> simp

/-- Taubenloch-Summe ≥ 1 für drei {0,1}-Observablen auf demselben Fenster. -/
theorem bell_triple_sum_ge_one (oE oA oC : Bool) :
    (if oE = oA then 1 else 0) + (if oE = oC then 1 else 0) + (if oA = oC then 1 else 0) ≥ 1 := by
  rcases pigeonhole_three_bits oE oA oC with h | h | h <;> simp [h] <;> omega

/-!
### Testfall (computabel, ohne Primzahlen)
-/

/-- Manuelle Folge: genau ein ABCEA (Index 0) und ein CEABC (Index 5). -/
def testClassesManual : List EabcLetter :=
  [⟨1, by decide⟩, ⟨2, by decide⟩, ⟨3, by decide⟩, ⟨0, by decide⟩, ⟨1, by decide⟩,
   ⟨3, by decide⟩, ⟨0, by decide⟩, ⟨1, by decide⟩, ⟨2, by decide⟩, ⟨3, by decide⟩]

theorem test_manual_N_plus : N_plus testClassesManual = 1 := by
  unfold N_plus countSlidingWord testClassesManual wordABCEA
  decide

theorem test_manual_N_minus : N_minus testClassesManual = 1 := by
  unfold N_minus countSlidingWord testClassesManual wordCEABC
  decide

theorem test_manual_counts :
    N_plus testClassesManual = 1 ∧ N_minus testClassesManual = 1 :=
  ⟨test_manual_N_plus, test_manual_N_minus⟩

theorem test_manual_D_E_zero : D_E_on testClassesManual = 0 := by
  simp [D_E_on, test_manual_N_plus, test_manual_N_minus]

/-!
### Prime-basierte Zählung (Schicht B — aus `PatternCount`)

**Primärträger:** `windows5` auf `kappaPrimeStreamUpTo X` — `N_plus_up_to`, `N_minus_up_to`,
`D_E`, `Q_E`, `W_E`, `R_half` (vgl. `collatz_eabc_holonomie_fehlerterm.py`).

**Sekundärträger:** Primvierlinge — `N_plus_quadruplet_up_to` (vgl. `eabc_quadruplets_1e10.py`).

Asymptotik `Hol_E_zero`, Primdichte und R_{1/2} = O(1) bleiben Schicht R (`sorry`).
-/

/-- Alias: Gleitfenster-Zählung = `PatternCount.N_plus_up_to`. -/
abbrev N_plus_sliding_up_to := N_plus_up_to

/-- Alias: Gleitfenster-Zählung = `PatternCount.N_minus_up_to`. -/
abbrev N_minus_sliding_up_to := N_minus_up_to

/-- R_β(X) = D_E(X) / Q_E(X)^β (Skalierungsobservable, Ebene 2). -/
noncomputable def R_beta_up_to (β : ℝ) (X : ℕ) : ℝ :=
  let q := Q_E X
  if _h : q = 0 then 0
  else (D_E X : ℝ) / (q : ℝ) ^ β

/-- Legacy-Alias zu `PatternCount.W_E_zero_of_balance`. -/
abbrev W_E_up_to_zero_of_balance := W_E_zero_of_balance

/-- Legacy-Alias zu `PatternCount.W_E_bounds`. -/
abbrev W_E_up_to_bounds := W_E_bounds

/-- Hauptvermutung Hol_E = 0 (asymptotisch; Schicht R — Experiment, kein Theorem). -/
def Hol_E_zero : Prop :=
  sorry

/-!
### Referenzwerte (computabel)

Gleitfenster (`holonomy_counts`): X=1000 ⇒ N₊=N₋=4, D_E=0.
Primvierlinge (`eabc_quadruplets_1e10.py`): X=1000 ⇒ N₊^quad=3, N₋^quad=2.
-/

example : N_plus_up_to 1000 = 4 := by native_decide

example : N_minus_up_to 1000 = 4 := by native_decide

example : D_E 1000 = 0 := by native_decide

#eval N_plus_up_to 1000000

#eval N_minus_up_to 1000000

#eval N_plus_quadruplet_up_to 1000000

#eval N_minus_quadruplet_up_to 1000000

/-!
### Bell / CHSH auf G_E (Skeleton)

**Hypothese (Lokalität auf G_E):** vier Observablen pro Fenster faktorisieren über
einen versteckten Zustand λₙ = Gleitfenster Pₙ^(4) auf der Prim-Transportkette.

**Analogie:** klassische CHSH-Schranke |S| ≤ 2; Verletzung = nicht-faktorisierbare
Holonomie-Reste (vgl. D_E ≠ 0, |S_EABC| > 2 in `collatz_eabc_bell_holonomie.md` §12).
-/

/-- Vier {0,1}-Observablen auf demselben Fenster-Index. -/
structure EabcWindowObservables where
  sigmaE : Bool
  sigmaC : Bool
  oPfad : Bool
  oHol : Bool

/-- Lokale Realität auf G_E: alle vier Bits sind Funktionen eines λₙ. -/
def LocalRealismOnGE (obs : ℕ → EabcWindowObservables) : Prop :=
  ∃ hidden : ℕ → Bool,
    ∀ n, obs n =
      ⟨hidden n, hidden n, hidden n, hidden n⟩

/-- Korrelation E(α,β) = 2·P[gleich] − 1 auf endlicher Stichprobe (Definition). -/
def correlation (xs ys : List Bool) : ℚ :=
  if _h : xs.length = 0 then 0
  else
    let agree := (xs.zip ys).filter (fun p => p.1 = p.2) |>.length
    (2 * agree : ℚ) / xs.length - 1

/-- CHSH-Summe S = E(a,b) − E(a,b') + E(a',b) + E(a',b'). -/
def chshSum (e_ab e_abp e_apb e_apbp : ℚ) : ℚ :=
  e_ab - e_abp + e_apb + e_apbp

/-- CHSH-LHV-Schranke |S| ≤ 2 unter strikter Lokalität (vollständiger Beweis offen). -/
theorem chsh_lhv_bound_skel
    (_obs : ℕ → EabcWindowObservables)
    (_hLR : LocalRealismOnGE _obs) : True :=
  trivial

/- TODO: aus `LocalRealismOnGE` die Standard-CHSH-Ungleichung ableiten;
   vgl. `collatz_eabc_bell_holonomie.md` §7.5. -/

/-- Brücke D_E ↔ CHSH (Hypothese, dokumentiert in `collatz_eabc_bell_holonomie.md` §12). -/
def de_chsh_analogy_note : String :=
  "D_E bias = non-factorizable holonomy; chi_Hol bounded; S_EABC vs 2 ~ D_tilde_E"

end CollatzEabc
