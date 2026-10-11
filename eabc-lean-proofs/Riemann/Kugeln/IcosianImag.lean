import Kugeln.GramSeparation
import Mathlib.Data.Finset.Card
import Mathlib.Data.List.Dedup

/-!
Negativzertifikat für den Kandidaten `Π_Im`.

Die Quelle ist die klassische binäre Ikosaedergruppe `2I` mit 120
Elementen in `ℍ`. Eine Koordinate `a + b √5` wird als Paar `(4a, 4b)`
gespeichert, also als Zähler über dem Nenner 4. Die Abbildung

`Π_Im (a, b, c, d) = (b, c, d)`

verwirft den Realteil. Das Bild hat 75 Punkte. Die kleinste positive
Schale hat zwölf Punkte und die Kosinuszeile eines Ikosaeders,

`(-1) × 1`, `(-√5 / 5) × 5`, `(√5 / 5) × 5`.

Sie hat nicht den fcc-Fingerabdruck. Dies ist kein Satz über BK-I2.7D
und definiert `Π` für die Bamberg-Quelle nicht. Gate 3 bleibt unberührt.
-/

namespace Kugeln

/-- Zählerpaar `(a, b)` für die Koordinate `(a + b √5) / 4`. -/
abbrev Coord := ℤ × ℤ

/-- Quaternion in dieser Skalierung. -/
structure Quat where
  c0 : Coord
  c1 : Coord
  c2 : Coord
  c3 : Coord
deriving DecidableEq, Repr

/-- Imaginärteil, das Bild von `Π_Im`. -/
structure Imag where
  c1 : Coord
  c2 : Coord
  c3 : Coord
deriving DecidableEq, Repr

def cadd (x y : Coord) : Coord := (x.1 + y.1, x.2 + y.2)

def cmul (x y : Coord) : Coord :=
  (x.1 * y.1 + 5 * x.2 * y.2, x.1 * y.2 + x.2 * y.1)

def csq (x : Coord) : Coord := cmul x x

/-- Summe der Quadrate der vier Zähler. Für eine Einheitsquaternion ist das `(16, 0)`. -/
def qnorm (q : Quat) : Coord :=
  cadd (cadd (csq q.c0) (csq q.c1)) (cadd (csq q.c2) (csq q.c3))

def imagOf (q : Quat) : Imag := ⟨q.c1, q.c2, q.c3⟩

/-- Summe der Quadrate der drei Imaginärkoordinaten, noch ohne den Nenner 16. -/
def imagNorm (v : Imag) : Coord :=
  cadd (cadd (csq v.c1) (csq v.c2)) (csq v.c3)

def imagDot (v w : Imag) : Coord :=
  cadd (cadd (cmul v.c1 w.c1) (cmul v.c2 w.c2)) (cmul v.c3 w.c3)

def coordGet (q : Quat) : Fin 4 → Coord
  | 0 => q.c0
  | 1 => q.c1
  | 2 => q.c2
  | 3 => q.c3

def reorder (p : Fin 4 × Fin 4 × Fin 4 × Fin 4) (q : Quat) : Quat :=
  ⟨coordGet q p.1, coordGet q p.2.1, coordGet q p.2.2.1, coordGet q p.2.2.2⟩

/-- Die zwölf geraden Permutationen, als Ursprungsindex je Zielposition. -/
def evenIdx : List (Fin 4 × Fin 4 × Fin 4 × Fin 4) :=
  [(0, 1, 2, 3), (0, 2, 3, 1), (0, 3, 1, 2), (1, 0, 3, 2), (1, 2, 0, 3), (1, 3, 2, 0),
    (2, 0, 1, 3), (2, 1, 3, 0), (2, 3, 0, 1), (3, 0, 2, 1), (3, 1, 0, 2), (3, 2, 1, 0)]

def bools : List Bool := [false, true]

def negCoord (x : Coord) : Coord := (-x.1, -x.2)

def signCoord (s : Bool) (x : Coord) : Coord := if s then x else negCoord x

def axes : List Quat :=
  [⟨(4, 0), (0, 0), (0, 0), (0, 0)⟩, ⟨(-4, 0), (0, 0), (0, 0), (0, 0)⟩,
    ⟨(0, 0), (4, 0), (0, 0), (0, 0)⟩, ⟨(0, 0), (-4, 0), (0, 0), (0, 0)⟩,
    ⟨(0, 0), (0, 0), (4, 0), (0, 0)⟩, ⟨(0, 0), (0, 0), (-4, 0), (0, 0)⟩,
    ⟨(0, 0), (0, 0), (0, 0), (4, 0)⟩, ⟨(0, 0), (0, 0), (0, 0), (-4, 0)⟩]

def halves : List Quat :=
  bools.flatMap fun a =>
    bools.flatMap fun b =>
      bools.flatMap fun c =>
        bools.flatMap fun d =>
          [⟨signCoord a (2, 0), signCoord b (2, 0), signCoord c (2, 0), signCoord d (2, 0)⟩]

def tailRaw (s1 s2 s3 : Bool) : Quat :=
  ⟨(0, 0), signCoord s1 (2, 0), signCoord s2 (-1, 1), signCoord s3 (1, 1)⟩

def permsSrc (i : Fin 4) : List (Fin 4 × Fin 4 × Fin 4 × Fin 4) :=
  evenIdx.filter fun p => p.1 = i

def tailSrc (i : Fin 4) : List Quat :=
  bools.flatMap fun s1 =>
    bools.flatMap fun s2 =>
      bools.flatMap fun s3 =>
        (permsSrc i).map fun p => reorder p (tailRaw s1 s2 s3)

/-- Gerade Permutationen von `(0, ±1/2, ±φ⁻¹/2, ±φ/2)`. -/
def icosianTail : List Quat := tailSrc 0 ++ tailSrc 1 ++ tailSrc 2 ++ tailSrc 3

/-- Die 120 Elemente von `2I`. -/
def binaryIcosahedralList : List Quat := axes ++ halves ++ icosianTail

def binaryIcosahedral : Finset Quat := binaryIcosahedralList.toFinset

theorem binaryIcosahedralList_length : binaryIcosahedralList.length = 120 := by decide

theorem binaryIcosahedralList_nodup : binaryIcosahedralList.Nodup := by decide

theorem binaryIcosahedral_card : binaryIcosahedral.card = 120 := by
  rw [binaryIcosahedral, List.toFinset_card_of_nodup binaryIcosahedralList_nodup,
    binaryIcosahedralList_length]

theorem binaryIcosahedral_norm : ∀ q ∈ binaryIcosahedralList, qnorm q = (16, 0) := by
  intro q hq
  rw [binaryIcosahedralList, List.mem_append, List.mem_append] at hq
  rcases hq with (hq | hq) | hq
  · exact (by decide : ∀ q ∈ axes, qnorm q = (16, 0)) q hq
  · exact (by decide : ∀ q ∈ halves, qnorm q = (16, 0)) q hq
  · exact (by decide : ∀ q ∈ icosianTail, qnorm q = (16, 0)) q hq

/-- `Π_Im`. -/
def imagProjection : Quat → Imag := imagOf

/-- Imaginärteile mit festem Vorzeichen der Koordinate, die im Realteil landet. -/
def minimalImagList : List Imag :=
  bools.flatMap fun s1 =>
    bools.flatMap fun s2 =>
      (permsSrc 3).map fun p => imagOf (reorder p (tailRaw s1 s2 true))

def largeImagList : List Imag :=
  bools.flatMap fun s1 =>
    bools.flatMap fun s3 =>
      (permsSrc 2).map fun p => imagOf (reorder p (tailRaw s1 true s3))

def midTailImagList : List Imag :=
  bools.flatMap fun s2 =>
    bools.flatMap fun s3 =>
      (permsSrc 1).map fun p => imagOf (reorder p (tailRaw true s2 s3))

def outerTailImagList : List Imag :=
  bools.flatMap fun s1 =>
    bools.flatMap fun s2 =>
      bools.flatMap fun s3 =>
        (permsSrc 0).map fun p => imagOf (reorder p (tailRaw s1 s2 s3))

def halfImagList : List Imag :=
  bools.flatMap fun b =>
    bools.flatMap fun c =>
      bools.flatMap fun d =>
        [⟨signCoord b (2, 0), signCoord c (2, 0), signCoord d (2, 0)⟩]

def axisImagList : List Imag :=
  [⟨(4, 0), (0, 0), (0, 0)⟩, ⟨(-4, 0), (0, 0), (0, 0)⟩, ⟨(0, 0), (4, 0), (0, 0)⟩,
    ⟨(0, 0), (-4, 0), (0, 0)⟩, ⟨(0, 0), (0, 0), (4, 0)⟩, ⟨(0, 0), (0, 0), (-4, 0)⟩]

def originImag : Imag := ⟨(0, 0), (0, 0), (0, 0)⟩

def minimalShell : Finset Imag := minimalImagList.toFinset

def largeShell : Finset Imag := largeImagList.toFinset

def midShell : Finset Imag := midTailImagList.toFinset ∪ halfImagList.toFinset

def outerShell : Finset Imag := outerTailImagList.toFinset ∪ axisImagList.toFinset

def imagImage : Finset Imag :=
  {originImag} ∪ minimalShell ∪ midShell ∪ largeShell ∪ outerShell

theorem minimalImagList_nodup : minimalImagList.Nodup := by decide

theorem largeImagList_nodup : largeImagList.Nodup := by decide

theorem midTailImagList_nodup : midTailImagList.Nodup := by decide

theorem halfImagList_nodup : halfImagList.Nodup := by decide

theorem outerTailImagList_nodup : outerTailImagList.Nodup := by decide

theorem axisImagList_nodup : axisImagList.Nodup := by decide

theorem minimalShell_card : minimalShell.card = 12 := by
  rw [minimalShell, List.toFinset_card_of_nodup minimalImagList_nodup]
  decide

theorem largeShell_card : largeShell.card = 12 := by
  rw [largeShell, List.toFinset_card_of_nodup largeImagList_nodup]
  decide

theorem midShell_card : midShell.card = 20 := by
  have hdisj : Disjoint midTailImagList.toFinset halfImagList.toFinset := by decide
  rw [midShell, Finset.card_union_of_disjoint hdisj,
    List.toFinset_card_of_nodup midTailImagList_nodup,
    List.toFinset_card_of_nodup halfImagList_nodup]
  decide

theorem outerShell_card : outerShell.card = 30 := by
  have hdisj : Disjoint outerTailImagList.toFinset axisImagList.toFinset := by decide
  rw [outerShell, Finset.card_union_of_disjoint hdisj,
    List.toFinset_card_of_nodup outerTailImagList_nodup,
    List.toFinset_card_of_nodup axisImagList_nodup]
  decide

theorem shell_norms :
    (∀ v ∈ minimalImagList, imagNorm v = (10, -2)) ∧
      (∀ v ∈ largeImagList, imagNorm v = (10, 2)) ∧
        (∀ v ∈ midTailImagList, imagNorm v = (12, 0)) ∧
          (∀ v ∈ halfImagList, imagNorm v = (12, 0)) ∧
            (∀ v ∈ outerTailImagList, imagNorm v = (16, 0)) ∧
              (∀ v ∈ axisImagList, imagNorm v = (16, 0)) ∧
                imagNorm originImag = (0, 0) := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · decide
  · decide
  · decide
  · decide
  · decide
  · decide
  · decide

theorem imagImage_card : imagImage.card = 75 := by
  have h1 : Disjoint ({originImag} : Finset Imag) minimalShell := by decide
  have h2 : Disjoint ({originImag} ∪ minimalShell) midShell := by decide
  have h3 : Disjoint ({originImag} ∪ minimalShell ∪ midShell) largeShell := by decide
  have h4 : Disjoint ({originImag} ∪ minimalShell ∪ midShell ∪ largeShell) outerShell := by decide
  rw [imagImage, Finset.card_union_of_disjoint h4, Finset.card_union_of_disjoint h3,
    Finset.card_union_of_disjoint h2, Finset.card_union_of_disjoint h1, Finset.card_singleton,
    minimalShell_card, midShell_card, largeShell_card, outerShell_card]

theorem axes_image :
    (axes.map imagProjection).toFinset = {originImag} ∪ axisImagList.toFinset := by
  decide

theorem halves_image :
    (halves.map imagProjection).toFinset = halfImagList.toFinset := by
  decide

theorem tailSrc0_image :
    ((tailSrc 0).map imagProjection).toFinset = outerTailImagList.toFinset := by
  decide

theorem tailSrc1_image :
    ((tailSrc 1).map imagProjection).toFinset = midTailImagList.toFinset := by
  decide

theorem tailSrc2_image :
    ((tailSrc 2).map imagProjection).toFinset = largeShell := by
  decide

theorem tailSrc3_image :
    ((tailSrc 3).map imagProjection).toFinset = minimalShell := by
  decide

theorem projected_image_eq :
    (binaryIcosahedralList.map imagProjection).toFinset = imagImage := by
  simp only [binaryIcosahedralList, icosianTail, List.map_append, List.toFinset_append, axes_image,
    halves_image, tailSrc0_image, tailSrc1_image, tailSrc2_image, tailSrc3_image, imagImage,
    midShell, outerShell, largeShell, minimalShell]
  ext v
  simp only [Finset.mem_union, Finset.mem_singleton]
  tauto

theorem imagMinimalShell_fingerprint :
    ∀ v ∈ minimalShell,
      (minimalShell.filter fun w => imagDot v w = (-10, 2)).card = 1 ∧
        (minimalShell.filter fun w => imagDot v w = (2, -2)).card = 5 ∧
          (minimalShell.filter fun w => imagDot v w = (-2, 2)).card = 5 := by
  decide

noncomputable section

open Real

private lemma sqrt_five_sq : (sqrt 5) ^ 2 = 5 := sq_sqrt (by norm_num)

theorem imag_antipode_cosine :
    ((-10 : ℝ) + 2 * sqrt 5) / (10 - 2 * sqrt 5) = -1 := by
  have hsq := sqrt_five_sq
  have hden : (10 - 2 * sqrt 5 : ℝ) ≠ 0 := by
    nlinarith [sqrt_nonneg 5, hsq]
  field_simp [hden]
  nlinarith [hsq]

theorem imag_neighbor_cosine_pos :
    ((-2 : ℝ) + 2 * sqrt 5) / (10 - 2 * sqrt 5) = sqrt 5 / 5 := by
  have hsq := sqrt_five_sq
  have hden : (10 - 2 * sqrt 5 : ℝ) ≠ 0 := by
    nlinarith [sqrt_nonneg 5, hsq]
  field_simp [hden]
  nlinarith [hsq]

theorem imag_neighbor_cosine_neg :
    ((2 : ℝ) - 2 * sqrt 5) / (10 - 2 * sqrt 5) = -(sqrt 5 / 5) := by
  have hsq := sqrt_five_sq
  have hden : (10 - 2 * sqrt 5 : ℝ) ≠ 0 := by
    nlinarith [sqrt_nonneg 5, hsq]
  field_simp [hden]
  nlinarith [hsq]

theorem one_div_sqrt_five : 1 / sqrt 5 = sqrt 5 / 5 := by
  have hsq := sqrt_five_sq
  field_simp
  nlinarith [hsq]

theorem imagNorm_real_values :
    ((10 : ℝ) + (-2) * sqrt 5) / 16 = (5 - sqrt 5) / 8 ∧
      ((12 : ℝ) / 16 = 3 / 4) ∧
        ((10 : ℝ) + 2 * sqrt 5) / 16 = (5 + sqrt 5) / 8 ∧
          (16 : ℝ) / 16 = 1 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> ring

theorem imagMinimal_norm_lt_others :
    (5 - sqrt 5) / 8 < 3 / 4 ∧ (3 / 4 < (5 + sqrt 5) / 8) ∧ (5 + sqrt 5) / 8 < 1 := by
  have hsq := sqrt_five_sq
  refine ⟨?_, ?_, ?_⟩ <;> nlinarith [sqrt_nonneg 5, hsq]

/-- Die Zwölfer-Minimalschale ist ikosaedrisch, nicht fcc. -/
theorem imagMinimalShell_not_fcc (v : Imag) (hv : v ∈ minimalShell) :
    (minimalShell.filter fun w => imagDot v w = (-2, 2)).card = 5 ∧
      sqrt 5 / 5 ≠ 1 / 2 := by
  exact ⟨(imagMinimalShell_fingerprint v hv).2.2, icosahedral_neighbor_not_fcc_cosine⟩

end

end Kugeln
