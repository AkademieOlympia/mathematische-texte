/-
  CollatzEabc.Primvierling.Chirality — EABC-Flavors und Chiralitätswörter.

  Spiegelt Root-`EABC.lean` (`EClass`, `classOf`, `Chirality`, `chiralityOrder`) und
  `witness.py` / `eabc_from_lean.py` für das sorry-freie CollatzEabc-Paket.
-/

import Mathlib.Tactic
import CollatzEabc.Primvierling.Basic

namespace CollatzEabc
namespace Primvierling

/-- Die vier EABC-Familien (Reste 1, 5, 7, 11 mod 12). -/
inductive EClass where
  | E
  | A
  | B
  | C
  deriving DecidableEq, Repr

/-- Kanonische Restklasse jeder EABC-Familie. -/
def residue : EClass → ℕ
  | EClass.E => 1
  | EClass.A => 5
  | EClass.B => 7
  | EClass.C => 11

/-- Klassifikation nach den vier primfähigen Restklassen mod 12. -/
def classOf (n : ℕ) : Option EClass :=
  match n % 12 with
  | 1 => some EClass.E
  | 5 => some EClass.A
  | 7 => some EClass.B
  | 11 => some EClass.C
  | _ => none

/-- Chiralität: markierter Start legt ABCE bzw. CEAB fest. -/
inductive Chirality where
  | ABCE
  | CEAB
  deriving DecidableEq, Repr

/-- Chirale EABC-Reihenfolge auf dem Flavor-Ring. -/
def chiralityOrder : Chirality → List EClass
  | Chirality.ABCE => [EClass.A, EClass.B, EClass.C, EClass.E]
  | Chirality.CEAB => [EClass.C, EClass.E, EClass.A, EClass.B]

/-- Chiralitätswort aus markiertem Startpunkt (vgl. `witness.chirality_word`). -/
def chiralityWord (p : ℕ) : Option Chirality :=
  match p % 12 with
  | 5 => some Chirality.ABCE
  | 11 => some Chirality.CEAB
  | _ => none

/-- EABC-Rotation T: E → A → B → C → E. -/
def T : EClass → EClass
  | EClass.E => EClass.A
  | EClass.A => EClass.B
  | EClass.B => EClass.C
  | EClass.C => EClass.E

theorem classOf_residue (c : EClass) : classOf (residue c) = some c := by
  cases c <;> simp [classOf, residue]

theorem residue_lt_twelve (c : EClass) : residue c < 12 := by
  cases c <;> decide

theorem residue_mod12 (c : EClass) : residue c % 12 = residue c :=
  Nat.mod_eq_of_lt (residue_lt_twelve c)

/-- Restklassenfolge eines Chiralitätsworts. -/
def chiralityResidues (ch : Chirality) : List ℕ :=
  (chiralityOrder ch).map residue

theorem chiralityResidues_ABCE :
    chiralityResidues Chirality.ABCE = [5, 7, 11, 1] := rfl

theorem chiralityResidues_CEAB :
    chiralityResidues Chirality.CEAB = [11, 1, 5, 7] := rfl

/-- EABC-Klassen der vier Glieder von Q(p). -/
def quadrupletFlavors (p : ℕ) : List (Option EClass) :=
  (Q p).map classOf

/-- Mod-12-Reste der vier Glieder von Q(p). -/
def quadrupletResidues (p : ℕ) : List ℕ :=
  (Q p).map (· % 12)

theorem quadrupletResidues_eq (p : ℕ) :
    quadrupletResidues p = [p % 12, (p + 2) % 12, (p + 6) % 12, (p + 8) % 12] := by
  simp [quadrupletResidues, Q]

end Primvierling
end CollatzEabc
