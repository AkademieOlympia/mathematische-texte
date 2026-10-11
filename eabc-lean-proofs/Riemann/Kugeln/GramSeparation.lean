import Kugeln.Kissing

/-!
BK-I2.7D in drei Schritten. Keiner ist hier bewiesen.

1. Konstruktion: eine Menge \(S\) und eine Abbildung \(\Pi : S \to \mathbb{R}^3\).
2. Diskretheit und Kontakt: das Bild ist diskret und seine Minimalschale hat zwölf Punkte.
3. fcc-Identifikation: die Gram-Matrix des Bildes stimmt bis auf orthogonale
   Transformation und Umnummerierung mit der fcc-Kontaktschale überein.

Schritt 2 und Schritt 3 setzen Schritt 1 voraus. \(\Pi\) wird hier nicht definiert.

Bewiesen ist ein Ausschluss der beiden Standardwinkel und das
Abnahmezertifikat der fcc-Rohschale. Der Zentralkosinus benachbarter Ecken
des regulären Ikosaeders in den Koordinaten \((0,\pm 1,\pm\varphi)\) ist
\(\sqrt{5}/5\). Er ist von \(1/2\) verschieden. Die ganzzahligen
Skalarprodukte eines festen Rohvektors mit den übrigen elf Vektoren der
Schale \((\pm 1,\pm 1,0)\) haben die Vielfachheiten eins, vier, zwei und
vier für die Werte \(-2,-1,0,1\). \(\sqrt{5}/3\) ist der Betrag des Kosinus
des Diederwinkels, nicht der Zentralkosinus.

Ein orthogonaler Projektor vom Rang 3 ist eine Dreiebene in
\(\mathrm{Gr}(3,8)\). Diese Datei zeichnet keine davon aus und definiert
\(\Pi\) nicht.
-/

noncomputable section

open Finset Real
open scoped RealInnerProductSpace InnerProductSpace

namespace Kugeln

lemma boolSign_mul_dichotomy (s t : Bool) :
    boolSign s * boolSign t = 1 ∨ boolSign s * boolSign t = -1 := by
  cases s <;> cases t <;> simp [boolSign]

lemma rawPoint_inner_mem (i i' : Fin 3) (s t s' t' : Bool)
    (hne : rawPoint i s t ≠ rawPoint i' s' t') :
    ∑ j, rawPoint i s t j * rawPoint i' s' t' j = -2 ∨
      ∑ j, rawPoint i s t j * rawPoint i' s' t' j = -1 ∨
        ∑ j, rawPoint i s t j * rawPoint i' s' t' j = 0 ∨
          ∑ j, rawPoint i s t j * rawPoint i' s' t' j = 1 := by
  by_cases hi : i = i'
  · subst hi
    have hsign : s ≠ s' ∨ t ≠ t' := by
      contrapose! hne
      simp [hne.1, hne.2]
    rw [sum_shift (fun j => rawPoint i s t j * rawPoint i s' t' j) i]
    simp only [rawPoint_self, rawPoint_add_one, rawPoint_add_two, mul_zero, zero_add]
    have hs : boolSign s * boolSign s' = if s = s' then (1 : ℤ) else -1 := by
      cases s <;> cases s' <;> simp [boolSign]
    have ht : boolSign t * boolSign t' = if t = t' then (1 : ℤ) else -1 := by
      cases t <;> cases t' <;> simp [boolSign]
    rw [hs, ht]
    cases s <;> cases t <;> cases s' <;> cases t' <;> simp_all
  · have hsplit : i' = i + 1 ∨ i' = i + 2 := by
      by_cases h1 : i' = i + 1
      · exact Or.inl h1
      · exact Or.inr (fin3_third (Ne.symm hi) h1)
    rw [sum_shift (fun j => rawPoint i s t j * rawPoint i' s' t' j) i]
    cases hsplit with
    | inl h =>
      subst h
      have ha : rawPoint (i + 1) s' t' i = boolSign t' := by
        simpa [fin3_succ_add_two] using rawPoint_add_two (i + 1) s' t'
      have hb : rawPoint (i + 1) s' t' (i + 1) = 0 := rawPoint_self (i + 1) s' t'
      have hc : rawPoint (i + 1) s' t' (i + 2) = boolSign s' := by
        have hstep : (i + 1) + 1 = i + 2 := by
          rw [add_assoc]
          exact congrArg (fun k : Fin 3 => i + k) (by decide : (1 + 1 : Fin 3) = 2)
        simpa [hstep] using rawPoint_add_one (i + 1) s' t'
      rw [rawPoint_self, rawPoint_add_one, rawPoint_add_two, ha, hb, hc]
      rcases boolSign_mul_dichotomy t s' with h1 | h1
      · simp [h1]
      · simp [h1]
    | inr h =>
      subst h
      have ha : rawPoint (i + 2) s' t' i = boolSign s' := by
        have hidx : (i + 2) + 1 = i := by
          rw [add_assoc]
          have h3 : ((2 + 1 : Fin 3) = 0) := by decide
          rw [h3, add_zero]
        simpa [hidx] using rawPoint_add_one (i + 2) s' t'
      have hb : rawPoint (i + 2) s' t' (i + 1) = boolSign t' := by
        simpa [fin3_two_add_two] using rawPoint_add_two (i + 2) s' t'
      have hc : rawPoint (i + 2) s' t' (i + 2) = 0 := rawPoint_self (i + 2) s' t'
      rw [rawPoint_self, rawPoint_add_one, rawPoint_add_two, ha, hb, hc]
      rcases boolSign_mul_dichotomy s t' with h1 | h1
      · simp [h1]
      · simp [h1]

private def pairInner (p q : Fin 3 × Bool × Bool) : ℤ :=
  ∑ j, rawPoint p.1 p.2.1 p.2.2 j * rawPoint q.1 q.2.1 q.2.2 j

theorem rawIndex_inner_multiplicities (p : Fin 3 × Bool × Bool) :
    (Finset.univ.filter fun q => pairInner p q = -2).card = 1 ∧
      (Finset.univ.filter fun q => pairInner p q = -1).card = 4 ∧
        (Finset.univ.filter fun q => pairInner p q = 0).card = 2 ∧
          (Finset.univ.filter fun q => pairInner p q = 1).card = 4 := by
  obtain ⟨i, s, t⟩ := p
  fin_cases i <;> fin_cases s <;> fin_cases t <;> decide

private lemma fccRaw_filter_card (p : Fin 3 × Bool × Bool) (z : ℤ) :
    (fccRaw.filter fun w => ∑ j, rawPoint p.1 p.2.1 p.2.2 j * w j = z).card =
      (Finset.univ.filter fun q => pairInner p q = z).card := by
  let f : Fin 3 × Bool × Bool → Fin 3 → ℤ := fun q => rawPoint q.1 q.2.1 q.2.2
  have himage :
      (Finset.univ.filter fun q => pairInner p q = z).image f =
        fccRaw.filter fun w => ∑ j, rawPoint p.1 p.2.1 p.2.2 j * w j = z := by
    ext w
    constructor
    · intro hw
      obtain ⟨q, hq, rfl⟩ := mem_image.mp hw
      have hqz : pairInner p q = z := (mem_filter.mp hq).2
      refine mem_filter.mpr ⟨?_, hqz⟩
      exact mem_image.mpr ⟨q, mem_univ q, rfl⟩
    · intro hw
      obtain ⟨hwraw, hwz⟩ := mem_filter.mp hw
      obtain ⟨q, -, rfl⟩ := mem_image.mp hwraw
      exact mem_image.mpr ⟨q, mem_filter.mpr ⟨mem_univ q, hwz⟩, rfl⟩
  rw [← himage, card_image_of_injective _ rawPoint_injective]

/-- Abnahmezertifikat der unnormierten fcc-Schale. Keine Abbildung \(\Pi\). -/
theorem fccRaw_inner_multiplicities {v : Fin 3 → ℤ} (hv : v ∈ fccRaw) :
    (fccRaw.filter fun w => ∑ j, v j * w j = -2).card = 1 ∧
      (fccRaw.filter fun w => ∑ j, v j * w j = -1).card = 4 ∧
        (fccRaw.filter fun w => ∑ j, v j * w j = 0).card = 2 ∧
          (fccRaw.filter fun w => ∑ j, v j * w j = 1).card = 4 := by
  obtain ⟨p, -, rfl⟩ := mem_image.mp hv
  rw [fccRaw_filter_card p (-2), fccRaw_filter_card p (-1), fccRaw_filter_card p 0,
    fccRaw_filter_card p 1]
  exact rawIndex_inner_multiplicities p

private lemma int_div_two_left_inj {z n : ℤ} : (z : ℝ) / 2 = (n : ℝ) / 2 ↔ z = n := by
  rw [div_eq_div_iff (by norm_num) (by norm_num)]
  norm_cast
  omega

/-- Dieselben Vielfachheiten für die Kosinus \(-1,-1/2,0,1/2\). -/
theorem fccRaw_cosine_multiplicities {v : Fin 3 → ℤ} (hv : v ∈ fccRaw) :
    (fccRaw.filter fun w => ((∑ j, v j * w j : ℤ) : ℝ) / 2 = -1).card = 1 ∧
      (fccRaw.filter fun w => ((∑ j, v j * w j : ℤ) : ℝ) / 2 = -1 / 2).card = 4 ∧
        (fccRaw.filter fun w => ((∑ j, v j * w j : ℤ) : ℝ) / 2 = 0).card = 2 ∧
          (fccRaw.filter fun w => ((∑ j, v j * w j : ℤ) : ℝ) / 2 = 1 / 2).card = 4 := by
  have h := fccRaw_inner_multiplicities hv
  have card_eq (n : ℤ) (c : ℝ) (hc : (n : ℝ) / 2 = c) :
      (fccRaw.filter fun w => ((∑ j, v j * w j : ℤ) : ℝ) / 2 = c).card =
        (fccRaw.filter fun w => ∑ j, v j * w j = n).card := by
    refine congrArg Finset.card (filter_congr fun w _ => ?_)
    rw [← hc]
    exact int_div_two_left_inj
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [card_eq (-2) (-1) (by norm_num)]
    exact h.1
  · rw [card_eq (-1) (-1 / 2) (by norm_num)]
    exact h.2.1
  · rw [card_eq 0 0 (by norm_num)]
    exact h.2.2.1
  · rw [card_eq 1 (1 / 2) (by norm_num)]
    exact h.2.2.2

/-- Benachbarte Standardecken \((0,1,\varphi)\) und \((1,\varphi,0)\). -/
def goldenRatio : ℝ := (1 + Real.sqrt 5) / 2

lemma goldenRatio_sq : goldenRatio ^ 2 = goldenRatio + 1 := by
  unfold goldenRatio
  have hsq : (Real.sqrt 5) ^ 2 = 5 := Real.sq_sqrt (by norm_num)
  field_simp
  nlinarith [hsq]

lemma goldenRatio_add_two : goldenRatio + 2 = (5 + Real.sqrt 5) / 2 := by
  unfold goldenRatio
  ring

def icoNeighborA : Space3 := WithLp.toLp 2 ![0, 1, goldenRatio]

def icoNeighborB : Space3 := WithLp.toLp 2 ![1, goldenRatio, 0]

lemma icoNeighbor_cosine : ⟪icoNeighborA, icoNeighborB⟫ /
    (‖icoNeighborA‖ * ‖icoNeighborB‖) = Real.sqrt 5 / 5 := by
  have hsq : (Real.sqrt 5) ^ 2 = 5 := Real.sq_sqrt (by norm_num)
  have hnorm : ‖icoNeighborA‖ ^ 2 = goldenRatio + 2 := by
    rw [EuclideanSpace.norm_eq, Real.sq_sqrt (Finset.sum_nonneg fun _ _ => sq_nonneg _)]
    simp [icoNeighborA, PiLp.toLp_apply, Fin.sum_univ_succ, goldenRatio_sq]
    ring
  have hnormB : ‖icoNeighborB‖ ^ 2 = goldenRatio + 2 := by
    rw [EuclideanSpace.norm_eq, Real.sq_sqrt (Finset.sum_nonneg fun _ _ => sq_nonneg _)]
    simp [icoNeighborB, PiLp.toLp_apply, Fin.sum_univ_succ, goldenRatio_sq]
    ring
  have hinner : ⟪icoNeighborA, icoNeighborB⟫ = goldenRatio := by
    simp [icoNeighborA, icoNeighborB, PiLp.inner_apply, RCLike.inner_apply, Fin.sum_univ_succ]
  have hpos : 0 < goldenRatio + 2 := by
    unfold goldenRatio
    positivity
  have hnormA : ‖icoNeighborA‖ = Real.sqrt (goldenRatio + 2) := by
    rw [← Real.sqrt_sq (norm_nonneg _), hnorm]
  have hnormEq : ‖icoNeighborB‖ = Real.sqrt (goldenRatio + 2) := by
    rw [← Real.sqrt_sq (norm_nonneg _), hnormB]
  rw [hinner, hnormA, hnormEq, Real.mul_self_sqrt hpos.le, goldenRatio_add_two]
  unfold goldenRatio
  field_simp
  calc
    (1 + Real.sqrt 5) * 5 = 5 + Real.sqrt 5 * 5 := by rw [add_mul, one_mul]
    _ = 5 + 5 * Real.sqrt 5 := by rw [mul_comm (Real.sqrt 5) (5 : ℝ)]
    _ = 5 * Real.sqrt 5 + 5 := by rw [add_comm]
    _ = Real.sqrt 5 * 5 + 5 := by rw [mul_comm (5 : ℝ) (Real.sqrt 5)]
    _ = Real.sqrt 5 * 5 + (Real.sqrt 5) ^ 2 := by
      conv_lhs => arg 2; rw [← hsq]
    _ = Real.sqrt 5 * (5 + Real.sqrt 5) := by rw [mul_add, pow_two]

theorem sqrt_five_div_five_ne_fcc (z : ℤ) (hz : z = -2 ∨ z = -1 ∨ z = 0 ∨ z = 1) :
    Real.sqrt 5 / 5 ≠ (z : ℝ) / 2 := by
  intro h
  have hsq : (Real.sqrt 5) ^ 2 = 5 := Real.sq_sqrt (by norm_num)
  field_simp at h
  norm_cast at h
  have hpow := congrArg (fun t : ℝ => t ^ 2) h
  rw [mul_pow, hsq] at hpow
  rcases hz with rfl | rfl | rfl | rfl <;> norm_num at hpow

/-- Die Nachbarschaft des Ikosaeders kommt unter den fcc-Kosinus nicht vor. -/
theorem icosahedral_neighbor_not_fcc_cosine :
    Real.sqrt 5 / 5 ≠ 1 / 2 := by
  simpa using sqrt_five_div_five_ne_fcc 1 (by simp)

end Kugeln
