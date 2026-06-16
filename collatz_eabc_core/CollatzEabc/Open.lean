/-
  CollatzEabc.Open — Stufe E: dokumentierte offene Uniformitätsvermutung.

  Kein `sorry`, kein Beweis — nur Prop-Definitionen und logische Äquivalenzen.
  Siehe `collatz_offene_punkte.md` im Repo-Root für Kontext.
-/

import CollatzEabc.Z2Attraktor

namespace CollatzZ2

/-- Uniformitätsvermutung (OFFEN — nicht bewiesen, entspricht Collatz-Engpass). -/
def collatzUniformityConjecture : Prop :=
  ∀ n : ℕ, Odd n → ∃ K : ℕ, iterateU n K = 1

/-- Äquivalente Formulierung über endliche Attraktor-Approximation. -/
def collatzUniformityConjectureApprox : Prop :=
  ∀ n : ℕ, Odd n → ∃ K : ℕ, inTrivialAttraktorApprox n K

/-- Collatz für ungerade Starts (odd-to-odd-Formulierung). -/
theorem collatzUniformityConjecture_iff :
    collatzUniformityConjecture ↔ collatzUniformityConjectureApprox := by
  constructor
  · intro h n ho
    rcases h n ho with ⟨K, heq⟩
    exact ⟨K, ⟨K, Nat.le_refl _, heq⟩⟩
  · intro h n ho
    rcases h n ho with ⟨K, hK⟩
    rcases hK with ⟨k, _, heq⟩
    exact ⟨k, heq⟩

theorem collatzUniformityConjectureApprox_iff_iterate :
    collatzUniformityConjectureApprox ↔ collatzUniformityConjecture :=
  collatzUniformityConjecture_iff.symm

end CollatzZ2
