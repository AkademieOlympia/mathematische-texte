/-
  CollatzEabc.PatternCount — berechenbare Prim- und κ-Pattern-Zählung (Schicht B).

  Zwei Träger (vgl. `eabc_quadruplets_1e10.py`, `collatz_eabc_holonomie_fehlerterm.py`):

  1. Gleitfenster auf der EABC-Primfolge bis X (p > 3, Restklasse in {1,5,7,11}):
     `N_plus_up_to` / `N_minus_up_to` in `HolonomieFehlerterm` = `N_plus`/`N_minus` von `primeEabcClassesUpTo X`.

  2. Primvierlinge Q(p) = (p, p+2, p+6, p+8): Start p ≡ 5 (mod 12) ⇒ ABCEA-Träger,
     p ≡ 11 (mod 12) ⇒ CEABC-Träger — `N_plus_quadruplet_up_to`, `N_minus_quadruplet_up_to`.

  Schicht R (Asymptotik, Hol_E = 0, R_{1/2} = O(1)) bleibt in `HolonomieFehlerterm` als `sorry`.
-/

import Mathlib.Data.Nat.Prime.Defs
import Mathlib.Data.List.Basic
import CollatzEabc.Kappa

namespace CollatzEabc

open List

/-!
### Prim-Obergrenze und κ-Folge (Schicht B, computabel)
-/

/-- Entscheidbare Primzahlprüfung (trial division via Mathlib-`decide`). -/
def isPrime (n : ℕ) : Bool :=
  decide (Nat.Prime n)

/-- Primzahlen p mit 3 < p ≤ X und definierter EABC-Klasse (κ-Eintrag). -/
def primeEabcClassesUpTo (X : ℕ) : List EabcLetter :=
  (List.range (X + 1)).filterMap fun p =>
    if _ : 3 < p ∧ Nat.Prime p then classOfLetter p else none

/-- Anzahl EABC-Primklassen bis X. -/
def primeEabcCountUpTo (X : ℕ) : ℕ :=
  (primeEabcClassesUpTo X).length

/-!
### Primvierlinge Q(p) = (p, p+2, p+6, p+8) — Schicht B
-/

/-- Entscheidbares Primzahlvierling Q(p) (vgl. `eabc_from_lean.is_prime_quadruplet`). -/
def isPrimeQuadruplet (p : ℕ) : Bool :=
  isPrime p && isPrime (p + 2) && isPrime (p + 6) && isPrime (p + 8)

/-- N₊^quad(X): Vierlinge mit Start p ≡ 5 (mod 12) — ABCEA-Träger. -/
def isQuadrupletPlus (p : ℕ) : Bool :=
  isPrimeQuadruplet p && p % 12 = 5

/-- N₋^quad(X): Vierlinge mit Start p ≡ 11 (mod 12) — CEABC-Träger. -/
def isQuadrupletMinus (p : ℕ) : Bool :=
  isPrimeQuadruplet p && p % 12 = 11

/-- Zählt ABCEA-Träger-Vierlinge mit p ≤ X. -/
def N_plus_quadruplet_up_to (X : ℕ) : ℕ :=
  ((List.range (X + 1)).filter isQuadrupletPlus).length

/-- Zählt CEABC-Träger-Vierlinge mit p ≤ X. -/
def N_minus_quadruplet_up_to (X : ℕ) : ℕ :=
  ((List.range (X + 1)).filter isQuadrupletMinus).length

/-- D_E^quad(X) = N₊^quad − N₋^quad. -/
def D_E_quadruplet_up_to (X : ℕ) : ℤ :=
  (N_plus_quadruplet_up_to X : ℤ) - N_minus_quadruplet_up_to X

/-- Q_E^quad(X) = N₊^quad + N₋^quad. -/
def Q_E_quadruplet_up_to (X : ℕ) : ℕ :=
  N_plus_quadruplet_up_to X + N_minus_quadruplet_up_to X

/-!
### Kleine Referenzwerte (computabel, `native_decide`)
-/

example : N_plus_quadruplet_up_to 1000 = 3 := by native_decide

example : N_minus_quadruplet_up_to 1000 = 2 := by native_decide

example : D_E_quadruplet_up_to 1000 = 1 := by native_decide

end CollatzEabc
