import Mathlib.Analysis.InnerProductSpace.PiL2
import Mathlib.Geometry.Euclidean.Angle.Unoriented.Basic
import Mathlib.Tactic.FinCases

/-!
Spezifikation der Kusszahl in Dimension 3.

Eine Kusskonfiguration ist eine endliche Menge von Zentren im Abstand 2
vom Ursprung, deren paarweise Abstände mindestens 2 sind. Das ist die
Berührung einer zentralen Einheitskugel durch nicht überlappende
Einheitskugeln. Auf der Einheitskugel ist dieselbe Bedingung der
Winkelabstand mindestens \(\pi/3\).

Bewiesen wird diese Übersetzung und die untere Schranke 12 durch die
zwölf Vektoren vom Typ \((\pm 1,\pm 1,0)\), skaliert auf die Länge 2.
Die obere Schranke ist der Satz von Schütte und van der Waerden,
Math. Ann. 125 (1953), 325–334. Sie bleibt `sorry`.

Die Eulersche Zählung \(F = 2V - 4\) liefert für \(V = 13\) zweiundzwanzig
Dreiecke. Der Exzess des gleichseitigen sphärischen Dreiecks der
Seitenlänge \(\pi/3\) ist \(3\arccos(1/3) - \pi\), und zweiundzwanzig
solche Exzesse bleiben unter \(4\pi\). Diese Flächenbilanz schließt 13
nicht aus. Der Disput von 1694 war der zwischen Newton und David Gregory.

Dieses Modul importiert keine Datei aus `Riemann.*`. Es ist nicht
BK-I2.7D.
-/

noncomputable section

open InnerProductGeometry Finset Real
open scoped RealInnerProductSpace InnerProductSpace

namespace Kugeln

abbrev Space3 := EuclideanSpace ℝ (Fin 3)

/-- Zentren von Einheitskugeln, die die Einheitskugel im Ursprung berühren
und einander nicht durchdringen. -/
def IsKissingConfiguration (S : Finset Space3) : Prop :=
  (∀ x ∈ S, dist x (0 : Space3) = 2) ∧
  (∀ x ∈ S, ∀ y ∈ S, x ≠ y → dist x y ≥ 2)

theorem half_norm_of_kissing {S : Finset Space3} (hS : IsKissingConfiguration S) {x : Space3}
    (hx : x ∈ S) : ‖(1 / 2 : ℝ) • x‖ = 1 := by
  have hx2 : ‖x‖ = 2 := by simpa [dist_zero_right] using hS.1 x hx
  have hhalf : (0 : ℝ) ≤ 1 / 2 := by positivity
  rw [norm_smul, Real.norm_eq_abs, abs_of_nonneg hhalf, hx2]
  norm_num

theorem half_inner_of_kissing {S : Finset Space3} (hS : IsKissingConfiguration S) {x y : Space3}
    (hx : x ∈ S) (hy : y ∈ S) (hne : x ≠ y) :
    ⟪(1 / 2 : ℝ) • x, (1 / 2 : ℝ) • y⟫ ≤ 1 / 2 := by
  have hx2 : ‖x‖ = 2 := by simpa [dist_zero_right] using hS.1 x hx
  have hy2 : ‖y‖ = 2 := by simpa [dist_zero_right] using hS.1 y hy
  have hdist : 2 ≤ ‖x - y‖ := by simpa [dist_eq_norm] using hS.2 x hx y hy hne
  have hsq : 4 ≤ ‖x - y‖ ^ 2 := by nlinarith
  have hpol : ‖x - y‖ ^ 2 = ‖x‖ ^ 2 - 2 * ⟪x, y⟫_ℝ + ‖y‖ ^ 2 := norm_sub_sq_real x y
  rw [hx2, hy2] at hpol
  have hxy : ⟪x, y⟫_ℝ ≤ 2 := by linarith
  have hinner : ⟪x, y⟫ ≤ 2 := by simpa using hxy
  calc
    ⟪(1 / 2 : ℝ) • x, (1 / 2 : ℝ) • y⟫ = (1 / 2) * (1 / 2) * ⟪x, y⟫ := by
      rw [inner_smul_left, inner_smul_right, mul_assoc]
      simp
    _ ≤ (1 / 2) * (1 / 2) * 2 := by gcongr
    _ = 1 / 2 := by norm_num

/-- Auf der Einheitskugel ist die Kussbedingung der Winkel mindestens \(\pi/3\). -/
theorem angle_of_kissing {S : Finset Space3} (hS : IsKissingConfiguration S) {x y : Space3}
    (hx : x ∈ S) (hy : y ∈ S) (hne : x ≠ y) :
    π / 3 ≤ angle ((1 / 2 : ℝ) • x) ((1 / 2 : ℝ) • y) := by
  have hnx := half_norm_of_kissing hS hx
  have hny := half_norm_of_kissing hS hy
  have hinn := half_inner_of_kissing hS hx hy hne
  unfold angle
  rw [hnx, hny]
  have harg : ⟪(1 / 2 : ℝ) • x, (1 / 2 : ℝ) • y⟫ / (1 * 1) ≤ 1 / 2 := by
    simpa using hinn
  have hpi : Real.arccos (1 / 2) = π / 3 := by
    rw [← cos_pi_div_three]
    exact Real.arccos_cos (by positivity) (by linarith [pi_pos])
  have hmono := Real.arccos_le_arccos harg
  linarith

/-- Obere Schranke von Schütte und van der Waerden. Offen in dieser Bibliothek. -/
theorem schuette_van_der_waerden_3d (S : Finset Space3) (h : IsKissingConfiguration S) :
    S.card ≤ 12 := by
  sorry

/-! Die zwölf Zentren des Kuboktaeders. -/

def boolSign (b : Bool) : ℤ :=
  if b then 1 else -1

lemma boolSign_sq (b : Bool) : boolSign b ^ 2 = 1 := by
  cases b <;> simp [boolSign]

lemma boolSign_injective : Function.Injective boolSign := by
  intro a b h
  cases a <;> cases b <;> simp [boolSign] at h <;> rfl

/-- Eine Koordinate 0, die beiden anderen \(\pm 1\). Der Nullindex bestimmt die Achse. -/
def rawPoint (i : Fin 3) (s t : Bool) : Fin 3 → ℤ :=
  fun j => if j = i then 0 else if j = i + 1 then boolSign s else boolSign t

lemma fin3_add_ne_self {i k : Fin 3} (hk : k ≠ 0) : i + k ≠ i := by
  intro h
  apply hk
  have : i + k = i + 0 := by simpa using h
  exact add_left_cancel this

lemma fin3_succ_add_two (i : Fin 3) : i + 1 + 2 = i := by
  have h3 : ((1 + 2 : Fin 3) = 0) := by decide
  rw [add_assoc, h3, add_zero]

lemma fin3_two_add_two (i : Fin 3) : i + 2 + 2 = i + 1 := by
  have h4 : ((2 + 2 : Fin 3) = 1) := by decide
  rw [add_assoc, h4]

lemma fin3_third {i j : Fin 3} (h0 : j ≠ i) (h1 : j ≠ i + 1) : j = i + 2 := by
  fin_cases i <;> fin_cases j <;> simp_all

lemma rawPoint_self (i : Fin 3) (s t : Bool) : rawPoint i s t i = 0 := by
  simp [rawPoint]

lemma rawPoint_add_one (i : Fin 3) (s t : Bool) : rawPoint i s t (i + 1) = boolSign s := by
  have hne : (i + 1 : Fin 3) ≠ i := fin3_add_ne_self (by decide : (1 : Fin 3) ≠ 0)
  simp [rawPoint, hne]

lemma rawPoint_add_two (i : Fin 3) (s t : Bool) : rawPoint i s t (i + 2) = boolSign t := by
  have h1 : (i + 2 : Fin 3) ≠ i := fin3_add_ne_self (by decide : (2 : Fin 3) ≠ 0)
  have h2 : (i + 2 : Fin 3) ≠ i + 1 := by
    intro h
    have : (2 : Fin 3) = 1 := add_left_cancel h
    exact absurd this (by decide : (2 : Fin 3) ≠ 1)
  simp [rawPoint, h1, h2]

lemma rawPoint_eq_zero_iff (i j : Fin 3) (s t : Bool) : rawPoint i s t j = 0 ↔ j = i := by
  constructor
  · intro h
    by_contra hj
    by_cases h1 : j = i + 1
    · rw [h1, rawPoint_add_one] at h
      cases s <;> simp [boolSign] at h
    · rw [fin3_third hj h1, rawPoint_add_two] at h
      cases t <;> simp [boolSign] at h
  · intro hj
    simp [rawPoint, hj]

@[simp]
lemma sum_fin3 {α : Type*} [AddCommMonoid α] (f : Fin 3 → α) : ∑ i, f i = f 0 + f 1 + f 2 := by
  simp [Fin.sum_univ_succ, add_assoc]

lemma sum_shift {α : Type*} [AddCommMonoid α] (f : Fin 3 → α) (i : Fin 3) :
    ∑ j, f j = f i + f (i + 1) + f (i + 2) := by
  rw [← Equiv.sum_comp (Equiv.addLeft i) f, sum_fin3]
  simp [add_zero]

lemma rawPoint_sq_sum (i : Fin 3) (s t : Bool) : ∑ j, rawPoint i s t j ^ 2 = 2 := by
  rw [sum_shift (fun j => rawPoint i s t j ^ 2) i, rawPoint_self, rawPoint_add_one,
    rawPoint_add_two, boolSign_sq, boolSign_sq]
  ring

lemma rawPoint_injective :
    Function.Injective (fun p : Fin 3 × Bool × Bool => rawPoint p.1 p.2.1 p.2.2) := by
  rintro ⟨i, s, t⟩ ⟨i', s', t'⟩ h
  have hi : i = i' := by
    have hz : rawPoint i' s' t' i = 0 := by
      simpa [rawPoint_self, eq_comm] using congr_fun h i
    exact (rawPoint_eq_zero_iff i' i s' t').1 hz
  subst hi
  have hs : s = s' := boolSign_injective <| by
    have := congr_fun h (i + 1)
    simpa [rawPoint_add_one] using this
  have ht : t = t' := boolSign_injective <| by
    have := congr_fun h (i + 2)
    simpa [rawPoint_add_two] using this
  simp [hs, ht]

def fccRaw : Finset (Fin 3 → ℤ) :=
  Finset.univ.image (fun p : Fin 3 × Bool × Bool => rawPoint p.1 p.2.1 p.2.2)

lemma fccRaw_card : fccRaw.card = 12 := by
  rw [fccRaw, card_image_of_injective _ rawPoint_injective, card_univ]
  simp [Fintype.card_prod, Fintype.card_fin, Fintype.card_bool]

lemma mem_fccRaw_sq {v : Fin 3 → ℤ} (hv : v ∈ fccRaw) : ∑ j, v j ^ 2 = 2 := by
  obtain ⟨p, -, rfl⟩ := mem_image.mp hv
  exact rawPoint_sq_sum p.1 p.2.1 p.2.2

lemma int_norm_sq_polar (v w : Fin 3 → ℤ) :
    ∑ i, (v i - w i) ^ 2 = ∑ i, v i ^ 2 + ∑ i, w i ^ 2 - 2 * ∑ i, v i * w i := by
  calc
    ∑ i, (v i - w i) ^ 2 = ∑ i, (v i ^ 2 + w i ^ 2 - 2 * (v i * w i)) := by
      refine sum_congr rfl ?_
      intro i _
      ring
    _ = ∑ i, (v i ^ 2 + w i ^ 2) - ∑ i, (2 * (v i * w i)) := sum_sub_distrib _ _
    _ = ∑ i, v i ^ 2 + ∑ i, w i ^ 2 - 2 * ∑ i, v i * w i := by
      rw [sum_add_distrib, ← mul_sum]

lemma rawPoint_inner_le (i i' : Fin 3) (s t s' t' : Bool)
    (hne : rawPoint i s t ≠ rawPoint i' s' t') :
    ∑ j, rawPoint i s t j * rawPoint i' s' t' j ≤ 1 := by
  by_cases hi : i = i'
  · subst hi
    have hsign : s ≠ s' ∨ t ≠ t' := by
      contrapose! hne
      simp [hne.1, hne.2]
    rw [sum_shift (fun j => rawPoint i s t j * rawPoint i s' t' j) i]
    simp only [rawPoint_self, rawPoint_add_one, rawPoint_add_two, mul_zero, zero_add]
    cases s <;> cases t <;> cases s' <;> cases t' <;> simp_all [boolSign]
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
      cases t <;> cases s' <;> simp [boolSign]
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
      cases s <;> cases t' <;> simp [boolSign]

lemma mem_fccRaw_inner_le {v w : Fin 3 → ℤ} (hv : v ∈ fccRaw) (hw : w ∈ fccRaw) (hne : v ≠ w) :
    ∑ j, v j * w j ≤ 1 := by
  obtain ⟨p, -, rfl⟩ := mem_image.mp hv
  obtain ⟨q, -, rfl⟩ := mem_image.mp hw
  exact rawPoint_inner_le p.1 q.1 p.2.1 p.2.2 q.2.1 q.2.2 hne

def embed (v : Fin 3 → ℤ) : Space3 :=
  WithLp.toLp 2 (fun i => Real.sqrt 2 * (v i : ℝ))

lemma embed_apply (v : Fin 3 → ℤ) (i : Fin 3) : embed v i = Real.sqrt 2 * (v i : ℝ) := by
  simp [embed, PiLp.toLp_apply]

lemma embed_injective : Function.Injective embed := by
  intro v w h
  have hfun := WithLp.toLp_injective (p := 2) h
  funext i
  have hi := congr_fun hfun i
  have h2 : Real.sqrt 2 ≠ 0 := Real.sqrt_ne_zero'.mpr (by norm_num)
  exact Int.cast_injective (mul_left_cancel₀ h2 hi)

lemma embed_coord_sq (v : Fin 3 → ℤ) (i : Fin 3) :
    ‖embed v i‖ ^ 2 = 2 * (v i : ℝ) ^ 2 := by
  rw [embed_apply, Real.norm_eq_abs, abs_mul, abs_of_nonneg (Real.sqrt_nonneg _), mul_pow,
    Real.sq_sqrt (by norm_num), sq_abs]

lemma embed_norm_sq {v : Fin 3 → ℤ} (hv : ∑ j, v j ^ 2 = 2) : ‖embed v‖ ^ 2 = 4 := by
  have hnonneg : 0 ≤ ∑ i, ‖embed v i‖ ^ 2 := sum_nonneg fun i _ => sq_nonneg _
  rw [EuclideanSpace.norm_eq, Real.sq_sqrt hnonneg]
  simp_rw [embed_coord_sq, ← mul_sum]
  have hsum : ∑ i, (v i : ℝ) ^ 2 = 2 := by
    have hcast := congrArg (fun n : ℤ => (n : ℝ)) hv
    simpa [Int.cast_sum, Int.cast_pow] using hcast
  rw [hsum]
  norm_num

lemma embed_sub_norm_sq (v w : Fin 3 → ℤ) :
    ‖embed v - embed w‖ ^ 2 = 2 * ∑ i, ((v i : ℝ) - (w i : ℝ)) ^ 2 := by
  have hsub : embed v - embed w = WithLp.toLp 2
      (fun i => Real.sqrt 2 * ((v i : ℝ) - (w i : ℝ))) := by
    rw [embed, embed, ← WithLp.toLp_sub]
    congr 1
    funext i
    rw [Pi.sub_apply]
    ring
  have hnonneg : 0 ≤ ∑ i, ‖(embed v - embed w) i‖ ^ 2 := sum_nonneg fun i _ => sq_nonneg _
  rw [EuclideanSpace.norm_eq, Real.sq_sqrt hnonneg]
  have hcoord : ∀ i, ‖(embed v - embed w) i‖ ^ 2 =
      2 * ((v i : ℝ) - (w i : ℝ)) ^ 2 := by
    intro i
    have happ : (embed v - embed w) i = Real.sqrt 2 * ((v i : ℝ) - (w i : ℝ)) := by
      rw [hsub]
    rw [happ, Real.norm_eq_abs, abs_mul, abs_of_nonneg (Real.sqrt_nonneg _), mul_pow,
      Real.sq_sqrt (by norm_num), sq_abs]
  simp_rw [hcoord, ← mul_sum]

lemma embed_dist_sq {v w : Fin 3 → ℤ} (hv : ∑ j, v j ^ 2 = 2) (hw : ∑ j, w j ^ 2 = 2)
    (hin : ∑ j, v j * w j ≤ 1) : 4 ≤ ‖embed v - embed w‖ ^ 2 := by
  have hzz : (2 : ℤ) ≤ ∑ i, (v i - w i) ^ 2 := by
    rw [int_norm_sq_polar, hv, hw]
    linarith
  have hreal : (2 : ℝ) ≤ ∑ i, ((v i - w i : ℤ) : ℝ) ^ 2 := by
    exact_mod_cast hzz
  have hcast : ∑ i, ((v i : ℝ) - (w i : ℝ)) ^ 2 = ∑ i, ((v i - w i : ℤ) : ℝ) ^ 2 := by
    simp [Int.cast_sub]
  rw [embed_sub_norm_sq, hcast]
  linarith

def fccCenters : Finset Space3 :=
  fccRaw.image embed

lemma fccCenters_card : fccCenters.card = 12 := by
  rw [fccCenters, card_image_of_injective _ embed_injective, fccRaw_card]

theorem fccCenters_kissing : IsKissingConfiguration fccCenters := by
  constructor
  · intro x hx
    obtain ⟨v, hv, rfl⟩ := mem_image.mp hx
    rw [dist_zero_right]
    have hsq : ‖embed v‖ ^ 2 = 4 := embed_norm_sq (mem_fccRaw_sq hv)
    have h22 : (2 : ℝ) ^ 2 = 4 := by norm_num
    exact (sq_eq_sq₀ (norm_nonneg _) (by norm_num)).mp (by rwa [h22])
  · intro x hx y hy hne
    obtain ⟨v, hv, rfl⟩ := mem_image.mp hx
    obtain ⟨w, hw, rfl⟩ := mem_image.mp hy
    have hvw : v ≠ w := by
      intro h
      apply hne
      simp [h]
    rw [dist_eq_norm]
    have hsq : 4 ≤ ‖embed v - embed w‖ ^ 2 :=
      embed_dist_sq (mem_fccRaw_sq hv) (mem_fccRaw_sq hw) (mem_fccRaw_inner_le hv hw hvw)
    have h22 : (2 : ℝ) ^ 2 = 4 := by norm_num
    have habs := sq_le_sq.mp (show (2 : ℝ) ^ 2 ≤ ‖embed v - embed w‖ ^ 2 by rwa [h22])
    simpa [abs_of_nonneg (norm_nonneg _), abs_of_nonneg (by norm_num : (0 : ℝ) ≤ 2)] using habs

/-- Zwölf Einheitskugeln können eine zentrale Einheitskugel berühren. -/
theorem exists_kissingConfiguration_card_twelve :
    ∃ S : Finset Space3, IsKissingConfiguration S ∧ S.card = 12 :=
  ⟨fccCenters, fccCenters_kissing, fccCenters_card⟩

end Kugeln
