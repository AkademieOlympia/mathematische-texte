import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Data.ZMod.Basic
import Mathlib.Tactic

namespace EABC

/-- Die vier EABC-Familien als Restklassen modulo 12. -/
inductive EClass where
  | E
  | A
  | B
  | C
  deriving DecidableEq, Repr

/-- Zuordnung der EABC-Klassen zu ihren Restklassen modulo 12.

E ≡ 1, A ≡ 5, B ≡ 7, C ≡ 11 mod 12.
-/
def residue : EClass → Nat
  | EClass.E => 1
  | EClass.A => 5
  | EClass.B => 7
  | EClass.C => 11

/-- Klassifikation einer natürlichen Zahl nach den vier primfähigen Restklassen modulo 12.
Zahlen in anderen Restklassen erhalten `none`. -/
def classOf (n : Nat) : Option EClass :=
  match n % 12 with
  | 1  => some EClass.E
  | 5  => some EClass.A
  | 7  => some EClass.B
  | 11 => some EClass.C
  | _  => none

/-- Ein Primzahlvierling im engeren Sinn: p, p+2, p+6, p+8 sind prim. -/
def IsPrimeQuadruplet (p : Nat) : Prop :=
  Nat.Prime p ∧ Nat.Prime (p + 2) ∧ Nat.Prime (p + 6) ∧ Nat.Prime (p + 8)

/-- Das Vierlingsobjekt Q(p). -/
def Q (p : Nat) : List Nat :=
  [p, p + 2, p + 6, p + 8]

/-- Mittelpunkt des Vierlings. -/
def center (p : Nat) : Nat :=
  p + 4

/-- Die interne Gap-Signatur des Primzahlvierlings. -/
def gapSignature : List Nat :=
  [2, 4, 2]

/-- Chiralität: Start bei A ergibt ABCE, Start bei C ergibt CEAB. -/
inductive Chirality where
  | ABCE
  | CEAB
  deriving DecidableEq, Repr

/-- Chirale Ordnung als Liste von EABC-Klassen. -/
def chiralityOrder : Chirality → List EClass
  | Chirality.ABCE => [EClass.A, EClass.B, EClass.C, EClass.E]
  | Chirality.CEAB => [EClass.C, EClass.E, EClass.A, EClass.B]

/-- EABC-Rotation T: E → A → B → C → E.

Diese Version rotiert entlang der zyklischen Reihenfolge E-A-B-C.
Für die Vierlingsfolge ABCE erscheint sie als zyklische Verschiebung der Liste.
-/
def T : EClass → EClass
  | EClass.E => EClass.A
  | EClass.A => EClass.B
  | EClass.B => EClass.C
  | EClass.C => EClass.E

/-- Vierfache Anwendung des Operators T. -/
def T4 (x : EClass) : EClass :=
  T (T (T (T x)))

/-- Der Operator T hat Ordnung 4. -/
theorem T4_has_order_4 (x : EClass) : T4 x = x := by
  cases x <;> rfl

end EABC
