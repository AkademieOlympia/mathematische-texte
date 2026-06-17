/-
  CollatzEabc.PrefProjection — wohldefinierte Vorabbildung Φ_pref auf EABC-Wörtern.

  Kein Collatz-Anspruch: nur kombinatorische Wortabbildung
  (Phasenindex, Radial ρ = 2^{-|w|}, Zeit t = |w|).
  Referenz: collatz_kepler_gedankenexperiment.tex, Abschnitt Φ_pref.
-/

import Mathlib.Data.Fin.Basic
import Mathlib.Data.List.Basic
import Mathlib.Data.Rat.Defs

namespace CollatzEabc

/-- EABC-Buchstaben als `Fin 4` (E=0, A=1, B=2, C=3); vgl. `EabcIndex` in `Mod12Matrix`. -/
abbrev EabcLetter := Fin 4

/-- Phasenindex ι_E: E↦0, A↦1, B↦2, C↦3. -/
def iotaE (c : EabcLetter) : ℕ :=
  c.val

/-- Kumulative Phase Θ(w) = Σ ι_E(w_j). -/
def theta (w : List EabcLetter) : ℕ :=
  (w.map iotaE).sum

/-- Phase modulo 4. -/
def theta4 (w : List EabcLetter) : Fin 4 :=
  ⟨theta w % 4, Nat.mod_lt _ (by decide)⟩

/-- Einheitsphase als Gitterpaar in ℤ×ℤ (Lean-kompatibel, ohne ℂ). -/
def phaseUnit (k : Fin 4) : Int × Int :=
  match k with
  | 0 => (1, 0)
  | 1 => (0, 1)
  | 2 => (-1, 0)
  | 3 => (0, -1)

/-- Einheitsphase des Worts: phaseUnit ∘ Θ₄. -/
def uPhase (w : List EabcLetter) : Int × Int :=
  phaseUnit (theta4 w)

/-- Radialkomponente ρ(w) = 2^{-|w|} als rationale Zahl. -/
def rho (w : List EabcLetter) : ℚ :=
  (1 : ℚ) / (2 ^ w.length)

/-- Zeitkoordinate t(w) = |w|. -/
def tLen (w : List EabcLetter) : ℕ :=
  w.length

/-- Φ_pref auf Wortebene: (uPhase, tLen, rho) — komplexe Spur z₀+ρ·u bleibt im TeX. -/
structure PrefImage where
  phase : Int × Int
  time : ℕ
  radial : ℚ
  deriving Repr

def phiPref (w : List EabcLetter) : PrefImage where
  phase := uPhase w
  time := tLen w
  radial := rho w

end CollatzEabc
