/-
  CollatzEabc.ArithLanguage — formale EABC-Grammatik L und Realisierbarkeit (Stufe 2).

  **Kein Collatz-Beweis.** BB-Verbot (sorry-frei); vollständige C-Ketten/EA-DFA
  in Python (`collatz_l_arith_test.py`); Lean-Erweiterung folgt.

  Referenz: `collatz_equivalenz_e_infty.tex`, `collatz_generalangriff_2026.md`.
-/

import CollatzEabc.Kappa

namespace CollatzEabc

/-- Enthält das Wort aufeinanderfolgende B-Buchstaben (Index 2)? -/
def hasBB (w : List EabcLetter) : Bool :=
  (w.zip w.tail).any fun (a, b) => a == ⟨2, by decide⟩ && b == ⟨2, by decide⟩

/-- Notwendige lokale Grammatikbedingung: `BB ∉ L` (vgl. Schlussartikel Prop. B→B). -/
def isGrammarValid (w : List EabcLetter) : Prop :=
  ¬ hasBB w

/--
  `w` ist arithmetisch realisierbar (`L_arith`), wenn ein ungerades `n` existiert,
  dessen volles κ-Präfix `w` ist (`Kappa.kappaPrefixWord`).
-/
def RealizableWord (w : List EabcLetter) : Prop :=
  ∃ n : ℕ, ∃ h : Odd n, ∃ hdef : kappaPrefixDefined n w.length h,
    kappaPrefixWord n w.length h hdef = w

end CollatzEabc
