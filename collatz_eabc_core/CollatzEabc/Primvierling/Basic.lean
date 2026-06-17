/-
  CollatzEabc.Primvierling.Basic — arithmetische Kernstruktur Q(p)=(p,p+2,p+6,p+8).

  Konsistent mit Root-`EABC.lean` (`IsPrimeQuadruplet`, `Q`, `center`, `gapSignature`)
  und `witness.py` (`quadruplet_center`, `centered_prime_quadruplet`, `QUADRUPLET_OFFSETS`).
-/

import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Tactic

namespace CollatzEabc
namespace Primvierling

/-- Ein Primzahlvierling im engeren Sinn: p, p+2, p+6, p+8 sind prim. -/
def IsPrimeQuadruplet (p : ℕ) : Prop :=
  Nat.Prime p ∧ Nat.Prime (p + 2) ∧ Nat.Prime (p + 6) ∧ Nat.Prime (p + 8)

/-- Das Vierlingsobjekt Q(p). -/
def Q (p : ℕ) : List ℕ :=
  [p, p + 2, p + 6, p + 8]

/-- Symmetrieanker M = p + 4 (vgl. `EABC.center`, `quadruplet_center`). -/
def center (p : ℕ) : ℕ :=
  p + 4

/-- Derselbe Anker in `ℤ`. -/
def centerZ (p : ℕ) : ℤ :=
  center p

/-- Arithmetischer Schwerpunkt; identisch mit `center`. -/
def quadrupletCentroid (p : ℕ) : ℕ :=
  center p

/-- Die interne Gap-Signatur des Primzahlvierlings. -/
def gapSignature : List ℕ :=
  [2, 4, 2]

/-- Vier Glieder als explizites Tupel (vgl. `quadruplet_members`). -/
def quadrupletMembers (p : ℕ) : ℕ × ℕ × ℕ × ℕ :=
  (p, p + 2, p + 6, p + 8)

/-- i-tes Glied von Q(p). -/
def memberAt (p : ℕ) : Fin 4 → ℕ
  | 0 => p
  | 1 => p + 2
  | 2 => p + 6
  | 3 => p + 8

/-- Kanonische zentrierte Offsets (-4, -2, 2, 4). -/
def centeredOffsets : ℤ × ℤ × ℤ × ℤ :=
  (-4, -2, 2, 4)

/-- Zentrierte Offsets als Funktion auf `Fin 4`. -/
def centeredOffsetAt : Fin 4 → ℤ
  | 0 => -4
  | 1 => -2
  | 2 => 2
  | 3 => 4

/-- Zentriertes i-tes Glied: Q(p)_i − M. -/
def centeredQuadrupletAt (p : ℕ) (i : Fin 4) : ℤ :=
  (memberAt p i : ℤ) - centerZ p

/-- Zentrierte Normalform Q(p) − M (vgl. `centered_prime_quadruplet`). -/
def centeredQuadruplet (p : ℕ) : ℤ × ℤ × ℤ × ℤ :=
  (centeredQuadrupletAt p 0, centeredQuadrupletAt p 1,
   centeredQuadrupletAt p 2, centeredQuadrupletAt p 3)

theorem quadrupletCentroid_eq_center (p : ℕ) :
    quadrupletCentroid p = center p := rfl

private theorem quadruplet_sum (p : ℕ) : p + (p + 2) + (p + 6) + (p + 8) = 4 * (p + 4) := by
  ring

/-- Arithmetischer Schwerpunkt von Q(p) ist M = p + 4. -/
theorem quadrupletCentroid_eq (p : ℕ) :
    (p + (p + 2) + (p + 6) + (p + 8)) / 4 = quadrupletCentroid p := by
  rw [quadruplet_sum, quadrupletCentroid, center]
  exact Nat.mul_div_cancel_left (p + 4) (by decide : 0 < 4)

theorem center_eq_quadrupletCentroid (p : ℕ) :
    center p = quadrupletCentroid p := rfl

theorem centeredOffsets_eq_offsetAt :
    centeredOffsets =
      (centeredOffsetAt 0, centeredOffsetAt 1, centeredOffsetAt 2, centeredOffsetAt 3) := by
  simp [centeredOffsets, centeredOffsetAt]

/-- Jedes Glied von Q(p) ist M plus dem kanonischen Offset. -/
theorem members_from_center (p : ℕ) (i : Fin 4) :
    (memberAt p i : ℤ) = centerZ p + centeredOffsetAt i := by
  fin_cases i <;> simp [memberAt, centerZ, center, centeredOffsetAt] <;> omega

theorem centeredQuadrupletAt_eq_offset (p : ℕ) (i : Fin 4) :
    centeredQuadrupletAt p i = centeredOffsetAt i := by
  unfold centeredQuadrupletAt
  have h := members_from_center p i
  linarith

/-- Zentrierte Normalform ist unabhängig vom Startpunkt p. -/
theorem centered_normal_form (p : ℕ) :
    centeredQuadruplet p = centeredOffsets := by
  unfold centeredQuadruplet centeredOffsets
  refine Prod.ext ?_ (Prod.ext ?_ (Prod.ext ?_ ?_)) <;>
    simp [centeredQuadrupletAt_eq_offset, centeredOffsetAt]

theorem centeredOffsetAt_values :
    centeredOffsetAt 0 = -4 ∧
    centeredOffsetAt 1 = -2 ∧
    centeredOffsetAt 2 = 2 ∧
    centeredOffsetAt 3 = 4 := by
  simp [centeredOffsetAt]

theorem Q_eq_memberAt (p : ℕ) (i : Fin 4) :
    (Q p).get i = memberAt p i := by
  fin_cases i <;> simp [Q, memberAt]

theorem quadrupletMembers_eq (p : ℕ) :
    quadrupletMembers p = (p, p + 2, p + 6, p + 8) := rfl

theorem gapSignature_values :
    gapSignature = [2, 4, 2] := rfl

/-- Aufeinanderfolgende Abstände in Q(p) sind 2, 4, 2. -/
theorem gapSignature_from_members (p : ℕ) :
    memberAt p 1 - memberAt p 0 = 2 ∧
    memberAt p 2 - memberAt p 1 = 4 ∧
    memberAt p 3 - memberAt p 2 = 2 := by
  simp [memberAt]

end Primvierling
end CollatzEabc
